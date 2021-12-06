module Main where

import Lib
import Data.SRTree 

x = Var 0 
y = Var 1 

domains = [(-5, 5), (-5, 5)]

shapes = [ Range (0, 2) 
  , PartialNonDecreasing 0 (0, 5)
  , PartialNonDecreasing 1 (0, 5)
  , PartialNonIncreasing 0 (-5, 0)
  , PartialNonIncreasing 1 (-5, 0)
  ]

expr1 = (1 / (1 + x^.(-4))) + (1 / (1 + y^.(-4)))
expr2 = (x^.4 + y^.4 + 2*x^.4 * y^.4) / (1 + x^.4 + y^.4 + 2*x^.4 * y^.4)

viol1 = getViolationFun InnerInterval shapes domains
viol2 = getViolationFun OuterInterval shapes domains
viol3 = getViolationFun (Sampling 1000) shapes domains 

main :: IO ()
main = do
    print $ viol1 expr1
    print $ viol2 expr1
    print $ viol3 expr1
    print $ viol1 expr2
    print $ viol2 expr2
    print $ viol3 expr2
