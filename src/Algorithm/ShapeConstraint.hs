{-# LANGUAGE FlexibleInstances #-}
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

data Shape   = Range (Double, Double)                     -- f(x) \in [a, b]
             | DiffRng Int (Double, Double)               -- d f(x) / dx \in [a, b]
             | NonIncreasing Int                          -- d f(x) / dx \in [-inf, 0]
             | NonDecreasing Int                          -- d f(x) / dx \in [0, inf]
             | PartialNonIncreasing Int (Double, Double)  -- d f(x) / dx \in [-inf, 0], a <= x <= b
             | PartialNonDecreasing Int (Double, Double)  -- d f(x) / dx \in [0, inf], a <= x <= b
             | Inflection Int Int                         -- d^2 f(x)/dx^2 == 0
             | Convex Int Int                             -- d^2 f(x)/dx^2 \in (0, Infinity)
             | Concave Int Int                            -- d^2 f(x)/dx^2 \in (-Infinitiy,0)
                    deriving (Show, Read)

-- inner (a, b) no sharper than that
-- outer (a, b) no wider than that
data Evaluator = InnerInterval | OuterInterval | Sampling Int | Hybrid | Bisection Int
                   deriving (Show, Read)

type ConstraintFun = SRTree Int Double -> Double
type Domains       = [(Double, Double)]
type KDomains      = Map (Kaucher Double)

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

-- convertDomains :: Shape -> (Double, Double) -> Int -> Kaucher Double
convertDomains (PartialNonIncreasing ix (a, b)) ds = M.insert ix (a <.< b) ds
convertDomains (PartialNonDecreasing ix (a, b)) ds = M.insert ix (a <.< b) ds 
convertDomains _ ds = ds
{-# INLINE convertDomains #-}

go f ds [] []         acc = acc
go f ds (x:xs) (z:zs) acc = go f ds xs zs (acc + cnstr)
   where 
     cnstr = evalConstraints z $ toTuple $ f x $ convertDomains z ds
{-# INLINE go #-}
toMap    = M.fromList . zip [0..] . map (uncurry (<.<))
{-# INLINE toMap #-}
getViolationFun :: Evaluator -> [Shape] -> Domains -> ConstraintFun
getViolationFun InnerInterval shapes domains t = go innerApprox ds ts shapes 0.0 
  where
    t'       = fmap singleton t
    ts       = map (t' `ofShape`) shapes
    ds       = toMap domains

getViolationFun OuterInterval shapes domains t = go outerApprox ds ts shapes 0.0
  where
    t'       = fmap singleton t
    ts       = map (t' `ofShape`) shapes
    ds       = toMap domains
{-    
getViolationFun (Sampling nSamples) shapes domains t = f
  where
    ds       = toMap domains    
    samples  = map (makeSamples nSamples . (`convertDomains` ds)) shapes

    evalSamples :: Maybe (LA.Vector Double) -> (Double, Double)
    evalSamples mxs = (fromMaybe (-1/0) mmin, fromMaybe (1/0) mmax)
      where
        mmax = LA.maxElement <$> mxs
        mmin = LA.minElement <$> mxs

    f   = let getSize s = LA.size $ s M.! 0
              ts        = zipWith (\shape s -> LA.fromList . replicate (getSize s) <$> (t `ofShape` shape)) shapes samples
              rngs      = zipWith (\t s -> evalSamples $ evalTreeWithMap t s) ts samples
              cnstrs    = zipWith evalConstraints shapes rngs
          in  sum cnstrs
-}
getViolationFun e shapes domains t = error $ show e <> " evaluator not implemented"

instance OptIntPow (LA.Vector Double) where 
    (^.) = (^^)
    {-# INLINE (^.) #-}
{-
makeSamples :: Int -> Map Int (Kaucher Double) -> [LA.Vector Double]
makeSamples n domains = LA.toColumns $ LA.fromLists $ map (map midpoint) $ head $ dropWhile (\x -> length x < n) samples
  where
    step x  = let (a,b) = bisect x in [a,b]
    f       = map (concatMap step)
    samples = map sequence $ iterate f $ map pure domains
{-# INLINE makeSamples #-}
-}
