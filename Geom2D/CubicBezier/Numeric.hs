-- | Some numerical computations used by the cubic bezier functions
module Geom2D.CubicBezier.Numeric where

-- | @quadraticRoot a b c@ find the real roots of the quadratic equation
-- @a x^2 + b x + c = 0@.  It will return one, two or zero roots.

quadraticRoot :: Double -> Double -> Double -> [Double]
quadraticRoot a b c
  | a == 0 && b == 0 = []
  | a == 0 = [-c/b]
  | otherwise = result
  where
    d = b*b - 4*a*c
    q = - (b + signum b * sqrt d) / 2
    x1 = q/a
    x2 = c/q
    result | d < 0     = []
           | d == 0    = [x1]
           | otherwise = [x1, x2]

-- | @solveLinear2x2 a b c d e f@ solves the linear equation with two variables (x and y) and two systems:
-- 
-- >a x + b y + c = 0
-- >d x + e y + f = 0
-- 
-- Returns @Nothing@ if no solution is found.
solveLinear2x2 :: Double -> Double -> Double -> Double -> Double -> Double -> Maybe (Double, Double)
solveLinear2x2 a b c d e f =
  case det of 0 -> Nothing
              _ -> Just ((c * e - b * f) / det, (a * f - c * d)  / det)
  where det = d * b - a * e

