module Geom2D.CubicBezier.Basic where
import Geom2D
import Geom2D.CubicBezier.Numeric
import Math.Polynomial

data CubicBezier = CubicBezier Point Point Point Point
                   deriving Show

-- | Return True if the param lies on the curve, iff it's in the interval @[0, 1]@.
bezierParam :: Double -> Bool
bezierParam t = t >= 0 || t <= 1

-- | Give the polynomial from 1D bezier parameters.
bezierToPoly1D :: Double -> Double -> Double -> Double -> Poly Double
bezierToPoly1D a b c d = poly BE [ d - a + 3*(b - c)
                                 , 3 * (a - 2*b + c)
                                 , 3 * (b - a)
                                 , a]

-- | Return a pair of polynomials from the bezier curve, one for each coordinate.
bezierToPoly :: CubicBezier -> (Poly Double, Poly Double)
bezierToPoly (CubicBezier (Point p1x p1y) (Point p2x p2y) (Point p3x p3y) (Point p4x p4y)) =
  (bezierToPoly1D p1x p2x p3x p4x,
   bezierToPoly1D p1y p2y p3y p4y)

-- | Calculate a value on the curve.
evalBezier :: CubicBezier -> Double -> Point
evalBezier b =
  let (px, py) = bezierToPoly b
  in \t -> Point (evalPoly px t) (evalPoly py t)

-- | Calculate a value and the first derivative on the curve.
evalBezierDeriv :: CubicBezier -> Double -> (Point, Point)
evalBezierDeriv b =
  let (px, py) = bezierToPoly b
      dpx = polyDeriv px
      dpy = polyDeriv py
  in \t -> (Point (evalPoly px t) (evalPoly py t),
            Point (evalPoly dpx t) (evalPoly dpy t))

-- | Calculate a value and all derivatives on the curve.
evalBezierDerivs :: CubicBezier -> Double -> [Point]
evalBezierDerivs b =
  let (px, py) = bezierToPoly b
  in \t -> zipWith Point (evalPolyDerivs px t) (evalPolyDerivs py t)

-- | @findBezierTangent p b@ finds the parameters where
-- the tangent of the bezier curve @b@ has the same direction as vector p.

-- Use the formula tx * B'y(t) - ty * B'x(t) = 0 where
-- B'x is the x value of the derivative of the Bezier curve.

findBezierTangent :: Point -> CubicBezier -> [Double]
findBezierTangent (Point tx ty) (CubicBezier (Point x0 y0) (Point x1 y1) (Point x2 y2) (Point x3 y3)) = 
  filter bezierParam $ quadraticRoot a b c
    where
      a = 3*(tx*((y3 - y0) + 3*(y1 - y2)) - ty*((x3 - x0) + 3*(x1 - x2)))
      b = 6*(tx*((y2 + y0) - 2*y1) - ty*((x2 + x0) - 2*x1))
      c = 3*(tx*(y1 - y0) - ty*(x1 - x0))

-- | Find the parameter where the bezier curve is horizontal.
bezierHoriz :: CubicBezier -> [Double]
bezierHoriz = findBezierTangent (Point 1 0)

-- | Find the parameter where the bezier curve is vertical.
bezierVert :: CubicBezier -> [Double]
bezierVert = findBezierTangent (Point 0 1)

-- | Find inflection points on the curve.

-- Use the formula B''x(t) * B'y(t) - B''y(t) * B'x(t) = 0
-- with B'x(t) the x value of the first derivative at t,
-- B''y(t) the y value of the second derivative at t
findBezierInflection :: CubicBezier -> [Double]
findBezierInflection (CubicBezier (Point x0 y0) (Point x1 y1) (Point x2 y2) (Point x3 y3)) = undefined

