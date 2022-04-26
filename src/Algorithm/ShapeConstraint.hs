{-# LANGUAGE FlexibleInstances #-}
{-|
Module      : Algorithm.ShapeConstraint
Description : evaluation of shape-constraints for symbolic regression
Copyright   : (c) Fabricio Olivetti de Franca, 2022
License     : GPL-3
Maintainer  : fabricio.olivetti@gmail.com
Stability   : experimental
Portability : POSIX

This package provides support functions to evaluate different 
shape-constraints for symbolic regression using modal arithmetic.

See https://direct.mit.edu/evco/article/30/1/75/99840 for more details.

@
@article{kronberger2022shape,
  title={Shape-Constrained Symbolic Regression—Improving Extrapolation with Prior Knowledge},
  author={Kronberger, Gabriel and de Fran{\c{c}}a, Fabricio Olivetti and Burlacu, Bogdan and Haider, Christian and Kommenda, Michael},
  journal={Evolutionary computation},
  volume={30},
  number={1},
  pages={75--98},
  year={2022},
  publisher={MIT Press One Rogers Street, Cambridge, MA 02142-1209, USA journals-info~…}
}
@

-}
module Algorithm.ShapeConstraint
    ( getViolationFun
    , Evaluator(..)
    , Shape(..)
    )
      where

import Numeric.ModalInterval
import Numeric.ModalInterval.Algorithms
import Data.SRTree

import Data.Map.Strict (Map(..))
import qualified Data.Map.Strict as M
import qualified Data.Vector as V
import Data.Bool (bool)
import Data.Maybe (fromMaybe)
import qualified Numeric.LinearAlgebra.Data as LA

-- | Data type describing the type of constraint:
--
-- * Range of the function
-- * Range of a partial derivative
-- * Monotonicity (non-increasing and non-decreasing)
-- * Partial monotonicity
-- * Inflection, Convex or Concave
--
-- The monotonicity constraints require the index of the
-- variable, the ranges require a tuple of the minimum and
-- maximum values, the partial monotonicity also requires
-- the range of the domain that has this property,
-- inflection, convexity and concavity requires the first
-- and second index of the partial derivatives.
data Shape   = Range (Double, Double)                     -- ^ f(x) \in [a, b]
             | DiffRng Int (Double, Double)               -- ^ d f(x) / dx \in [a, b]
             | NonIncreasing Int                          -- ^ d f(x) / dx \in [-inf, 0]
             | NonDecreasing Int                          -- ^ d f(x) / dx \in [0, inf]
             | PartialNonIncreasing Int (Double, Double)  -- ^ d f(x) / dx \in [-inf, 0], a <= x <= b
             | PartialNonDecreasing Int (Double, Double)  -- ^ d f(x) / dx \in [0, inf], a <= x <= b
             | Inflection Int Int                         -- ^ d^2 f(x)/dx^2 == 0
             | Convex Int Int                             -- ^ d^2 f(x)/dx^2 \in (0, Infinity)
             | Concave Int Int                            -- ^ d^2 f(x)/dx^2 \in (-Infinitiy,0)
                    deriving (Show, Read)

-- inner (a, b) no sharper than that
-- outer (a, b) no wider than that

-- | The `Evaluator` indicates the algorithm used to calculate the constraint:
--
-- * InnerInterval can result in false positives (i.e., it says it is feasible but it is not)
-- * OuterInterval can result in false negatives (i.e., it says it is infeasible but it is not)
-- * Sampling can result in false positives and usually requires a lot of samples
-- * Hybrid - not yet implemented
-- * Bisection can improve the results of OuterInterval
data Evaluator = InnerInterval | OuterInterval | Sampling Int | Hybrid | Bisection Int
                   deriving (Show, Read)

-- | Function that calculates how much of the constraints a symbolic tree vioaltes
type ConstraintFun = SRTree Int Double -> Double

-- | A list of domains of the variables
type Domains       = [(Double, Double)]

-- | A list of domains using Modal Arithmetic
type KDomains      = Map (Kaucher Double)

-- | Calculates the partial derivatives of a tree when needed.
ofShape :: (Floating a, Eq a, OptIntPow a) => SRTree Int a -> Shape -> SRTree Int a
t `ofShape` Range _ = t
t `ofShape` DiffRng ix _ = simplify $ deriveBy ix t
t `ofShape` NonIncreasing ix = simplify $ deriveBy ix t
t `ofShape` NonDecreasing ix = simplify $ deriveBy ix t
t `ofShape` PartialNonIncreasing ix _ = simplify $ deriveBy ix t
t `ofShape` PartialNonDecreasing ix _ = simplify $ deriveBy ix t
t `ofShape` Inflection ix iy = simplify $ deriveBy iy $ simplify $ deriveBy ix t
t `ofShape` Convex ix iy = simplify $ deriveBy iy $ simplify $ deriveBy ix t
t `ofShape` Concave ix iy = simplify $ deriveBy iy $ simplify $ deriveBy ix t
{-# INLINE ofShape #-}

toTuple :: Kaucher Double -> (Double, Double)
toTuple k = (fromMaybe (-1/0) $ inf k, fromMaybe (1/0) $ sup k)
{-# INLINE toTuple #-}

eps :: Double
eps = 1e-6
{-# INLINE eps #-}

calcPenalty :: Bool -> Double -> Double
calcPenalty b pnlty = bool 0 pnlty b
{-# INLINE calcPenalty #-}

inRng (a, b) (x, y) = calcPenalty (x < a - eps) (a - x) + calcPenalty (y > b + eps) (y - b)
{-# INLINE inRng #-}
isPositive (x, y)   = calcPenalty (x <= -eps) (abs x) + calcPenalty (y <= -eps) (abs y)
{-# INLINE isPositive #-}
isNegative (x, y)   = calcPenalty (x >= eps) (abs x) + calcPenalty (y >= eps) (abs y)
{-# INLINE isNegative #-}
isZero (x, y)       = calcPenalty (x <= -eps) (abs x) + calcPenalty (y >= eps) y
{-# INLINE isZero #-}

-- | Evaluates the constraint
evalConstraints :: Shape -> (Double, Double) -> Double
evalConstraints (Range ab) xy                 = inRng ab xy
evalConstraints (DiffRng _ ab) xy             = inRng ab xy
evalConstraints (NonIncreasing _) xy          = isNegative xy
evalConstraints (NonDecreasing _) xy          = isPositive xy
evalConstraints (PartialNonIncreasing _ _) xy = isNegative xy
evalConstraints (PartialNonDecreasing _ _) xy = isPositive xy
evalConstraints (Inflection _ _) xy           = isZero xy
evalConstraints (Convex _ _) xy               = isPositive xy
evalConstraints (Concave _ _) xy              = isNegative xy
{-# INLINE evalConstraints #-}

-- | Convert the domains for the partial monotonicity constraint
convertDomains :: Shape -> Map Int (Kaucher Double) -> Map Int (Kaucher Double)
convertDomains (PartialNonIncreasing ix (a, b)) ds = M.insert ix (a <.< b) ds
convertDomains (PartialNonDecreasing ix (a, b)) ds = M.insert ix (a <.< b) ds 
convertDomains _ ds = ds
{-# INLINE convertDomains #-}

-- | Sum a list of constraint violations
sumCnstr :: (SRTree Int (Kaucher Double) -> Map Int (Kaucher Double) -> Kaucher Double) -> Map Int (Kaucher Double) -> [SRTree Int (Kaucher Double)] -> [Shape] -> Double -> Double
sumCnstr f ds [] []         acc = acc
sumCnstr f ds (x:xs) (z:zs) acc = sumCnstr f ds xs zs (acc + cnstr)
   where 
     cnstr = evalConstraints z $ toTuple $ f x $ convertDomains z ds
{-# INLINE sumCnstr #-}

toMap :: Domains -> Map Int (Kaucher Double)
toMap    = M.fromList . zip [0..] . map (uncurry (<.<))
{-# INLINE toMap #-}

-- | Returns a constraint function given an `Evaluator` a list of `Shape` and the
-- `Domains` of the input variables.
getViolationFun :: Evaluator -> [Shape] -> Domains -> ConstraintFun
getViolationFun InnerInterval shapes domains t = sumCnstr innerApprox ds ts shapes 0.0 
  where
    t'       = fmap singleton t
    ts       = map (t' `ofShape`) shapes
    ds       = toMap domains

getViolationFun OuterInterval shapes domains t = sumCnstr outerApprox ds ts shapes 0.0
  where
    t'       = fmap singleton t
    ts       = map (t' `ofShape`) shapes
    ds       = toMap domains

getViolationFun (Sampling nSamples) shapes domains t = go shapes ts samples 0.0
  where
    ds        = toMap domains    
    samples   = map (M.fromList . zip [0..] . makeSamples nSamples . (`convertDomains` ds)) shapes
    getSize s = LA.size $ s M.! 0
    ts        = map (t `ofShape`) shapes
        
    evalSamples :: Maybe (LA.Vector Double) -> (Double, Double)
    evalSamples mxs = (fromMaybe (-1/0) mmin, fromMaybe (1/0) mmax)
      where
        mmax = LA.maxElement <$> mxs
        mmin = LA.minElement <$> mxs
    
    go [] _ _               acc = acc 
    go (x:xs) (y:ys) (z:zs) acc = go xs ys zs (acc + cnstr)
      where 
        cnstr = evalConstraints x rng
        y'    = fmap (LA.fromList . replicate (getSize z)) y
        rng   = evalSamples $ evalTreeWithMap y' z

getViolationFun e shapes domains t = error $ show e <> " evaluator not implemented"

instance OptIntPow (LA.Vector Double) where 
    (^.) = (^^)
    {-# INLINE (^.) #-}

-- | Creates the samples of the `Sampling` approach
makeSamples :: Int -> Map Int (Kaucher Double) -> [LA.Vector Double]
makeSamples n domains = LA.toColumns $ LA.fromLists $ map (map midpoint) $ head $ dropWhile (\x -> length x < n) samples
  where
    step x  = let (a,b) = bisect x in [a,b]
    f       = map (concatMap step)
    samples = map sequence $ iterate f $ map (pure.snd) $ M.toAscList domains
{-# INLINE makeSamples #-}
