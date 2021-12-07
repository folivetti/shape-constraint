module Main where

import Lib
import Data.SRTree 

x = Var 0 
y = Var 1 

sigma = Var 0 
theta = Var 1 

domains = [(-5, 5), (-5, 5)]
domains2 = [(1, 3), (1, 3)] 

shapes = [ Range (0, 2) 
  , PartialNonDecreasing 0 (0, 5)
  , PartialNonDecreasing 1 (0, 5)
  , PartialNonIncreasing 0 (-5, 0)
  , PartialNonIncreasing 1 (-5, 0)
  ]

shapes2 = [ NonIncreasing 1 ]

expr1 = (1 / (1 + x^.(-4))) + (1 / (1 + y^.(-4)))
expr2 = (x^.4 + y^.4 + 2*x^.4 * y^.4) / (1 + x^.4 + y^.4 + 2*x^.4 * y^.4)
expr3 = exp((-(theta/sigma)^.2)/2)/(sqrt(2*pi)*sigma)

viol1 = getViolationFun InnerInterval 
viol2 = getViolationFun OuterInterval 
viol3 = getViolationFun (Sampling 1000) 

main :: IO ()
main = do
    print $ viol1 shapes domains expr1
    print $ viol2 shapes domains expr1
    print $ viol3 shapes domains expr1

    print $ viol1 shapes domains expr2
    print $ viol2 shapes domains expr2
    print $ viol3 shapes domains expr2

    print $ viol1 shapes2 domains2 expr3
    print $ viol2 shapes2 domains2 expr3
    print $ viol3 shapes2 domains2 expr3
