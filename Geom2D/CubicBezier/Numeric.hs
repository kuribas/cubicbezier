-- | Some general numeric functions for the bezier library
module Geom2D.CubicBezier.Numeric where
import Numeric.GSL.Root

-- | find the real roots of a quadratic equation
quadraticRoot :: Double -> Double -> Double -> [Double]
quadraticRoot a b c = result where
  d = b*b - 4*a*c
  q = - (b + signum b * sqrt d) / 2
  x1 = q/a
  x2 = c/q
  result | d < 0     = []
         | d == 0    = [x1]
         | otherwise = [x1, x2]

-- | solve a linear equation with two variables and two systems:
-- @
-- a x + b y + c = 0
-- d x + e y + f = 0
-- @
solveLinear2x2 a b c d e f =
  case det of 0 -> Nothing
              _ -> Just ((c * e - b * f) / det, (a * f - c * d)  / det)
  where det = d * b - a * e

-- | find the root of a function numerically

-- use deckers method because it guarantees to converge
-- and has good convergence
findRoot f xacc xl xu = fst $ uniRoot Brent xacc 50 f xl xu
