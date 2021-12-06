module Lib
--    ( getViolationFun
--    , Evaluator(..)
--    , Shape(..)
--    ) 
      where

import Numeric.ModalInterval
import Numeric.ModalInterval.Algorithms
import Data.SRTree

import Data.Map.Strict (Map(..))
import qualified Data.Map.Strict as M
import qualified Data.Vector as V
import Data.Bool (bool)
import Data.Maybe (fromMaybe)

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

ofShape :: SRTree Int Double -> Shape -> SRTree Int Double 
t `ofShape` Range _ = t
t `ofShape` DiffRng ix _ = simplify $ deriveBy ix t
t `ofShape` NonIncreasing ix = simplify $ deriveBy ix t
t `ofShape` NonDecreasing ix = simplify $ deriveBy ix t
t `ofShape` PartialNonIncreasing ix _ = simplify $ deriveBy ix t
t `ofShape` PartialNonDecreasing ix _ = simplify $ deriveBy ix t
t `ofShape` Inflection ix iy = simplify $ deriveBy iy $ simplify $ deriveBy ix t
t `ofShape` Convex ix iy = simplify $ deriveBy iy $ simplify $ deriveBy ix t
t `ofShape` Concave ix iy = simplify $ deriveBy iy $ simplify $ deriveBy ix t

toTuple :: Kaucher Double -> (Double, Double)
toTuple k = (fromMaybe (-1/0) $ inf k, fromMaybe (1/0) $ sup k)
{-# INLINE toTuple #-}

eps :: Double
eps = 1e-16
{-# INLINE eps #-}

calcPenalty :: Bool -> Double -> Double
calcPenalty b pnlty = bool 0 pnlty b 
{-# INLINE calcPenalty #-}

inRng (a, b) (x, y) = calcPenalty (x < a - eps) (a - x) + calcPenalty (y > b + eps) (y - b)
isPositive (x, y)   = calcPenalty (x <= eps) (abs x) + calcPenalty (y <= eps) (abs y)
isNegative (x, y)   = calcPenalty (x >= -eps) (abs x) + calcPenalty (y >= -eps) (abs y)
isZero (x, y)       = calcPenalty (x <= -eps) (abs x) + calcPenalty (y >= eps) y 

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

convertDomains :: Shape -> (Double, Double) -> Int -> Kaucher Double 
convertDomains (PartialNonIncreasing ix (a, b)) (x, y) iy = if ix == iy then a <.< b else x <.< y
convertDomains (PartialNonDecreasing ix (a, b)) (x, y) iy = if ix == iy then a <.< b else x <.< y 
convertDomains _ (a, b) _ = a <.< b

getViolationFun :: Evaluator -> [Shape] -> Domains -> ConstraintFun
getViolationFun InnerInterval shapes domains = f
  where
    toMap    = M.fromList . zip [0..]
    domains' = map (toMap . (\s -> zipWith (convertDomains s) domains [0..])) shapes
    
    f t = let ts       = map (fmap singleton . (t `ofShape`)) shapes              
              rngs     = map toTuple $ zipWith innerApprox ts domains' 
              cnstrs   = zipWith evalConstraints shapes rngs
          in  sum cnstrs
          
getViolationFun OuterInterval shapes domains = f
  where
    toMap    = M.fromList . zip [0..]
    domains' = map (toMap . (\s -> zipWith (convertDomains s) domains [0..])) shapes
      
    f t = let ts       = map (fmap singleton . (t `ofShape`)) shapes
              rngs     = map toTuple $ zipWith outerApprox ts domains'
              cnstrs   = zipWith evalConstraints shapes rngs
          in  sum cnstrs
          
getViolationFun (Sampling nSamples) shapes domains = f
  where
    toMap    = M.fromList . zip [0..]
    samples  = map (map toMap . makeSamples nSamples . (\s -> zipWith (convertDomains s) domains [0..])) shapes
        
    evalSamples mxs = (fromMaybe (-1/0) mmin, fromMaybe (1/0) mmax)
      where
        mmax = maximum <$> mxs
        mmin = minimum <$> mxs 
              
    f t = let ts       = map (t `ofShape`) shapes
              rngs     = zipWith (\t s -> evalSamples $ traverse (evalTreeWithMap t) s) ts samples
              cnstrs   = zipWith evalConstraints shapes rngs
          in  sum cnstrs

getViolationFun e shapes domains = error $ show e <> " evaluator not implemented"

makeSamples :: Int -> [Kaucher Double] -> [[Double]]
makeSamples n domains = map (map midpoint) $ head $ dropWhile (\x -> length x < n) samples
  where
    step x  = let (a,b) = bisect x in [a,b]
    f       = map (concatMap step)
    samples = map sequence $ iterate f $ map pure domains
