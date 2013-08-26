module Geom2D.CubicBezier.Basic where
import Geom2D
import Geom2D.CubicBezier.Numeric
import Math.BernsteinPoly

data CubicBezier = CubicBezier {
  bezierC0 :: Point,
  bezierC1 :: Point,
  bezierC2 :: Point,
  bezierC3 :: Point} deriving Show

newtype Path = Path [PathSegment]
    
data PathSegment = LineSegment Line
                 | BezierSegment CubicBezier

instance AffineTransform CubicBezier where
  transform t (CubicBezier c0 c1 c2 c3) =
    CubicBezier (transform t c0) (transform t c1) (transform t c2) (transform t c3)

-- | Return True if the param lies on the curve, iff it's in the interval @[0, 1]@.
bezierParam :: Double -> Bool
bezierParam t = t >= 0 && t <= 1

-- | Convert a tolerance from the codomain to the domain of the bezier curve.
-- Should be good enough, but may not hold for high very tolerance values.

-- The magnification of error from the domain to the codomain of the
-- curve approaches the length of the tangent for small errors.  We
-- can use the maximum of the convex hull of the derivative, and double it to
-- have some margin for larger values.
bezierParamTolerance (CubicBezier p1 p2 p3 p4) eps = eps / maxDist
  where 
    maxDist = 6 * maximum [vectorDistance p1 p2,
                           vectorDistance p2 p3,
                           vectorDistance p3 p4]

-- | Reorient to the curve B(1-t).
reorient :: CubicBezier -> CubicBezier
reorient (CubicBezier p0 p1 p2 p3) = CubicBezier p3 p2 p1 p0 

-- | Give the bernstein polynomial for each coordinate.
bezierToBernstein (CubicBezier a b c d) = (listToBernstein $ map pointX coeffs,
                                           listToBernstein $ map pointY coeffs)
  where coeffs = [a, b, c, d]

-- | Calculate a value on the curve.
evalBezier :: CubicBezier -> Double -> Point
evalBezier b t = Point (bernsteinEval x t) (bernsteinEval y t)
  where (x, y) = bezierToBernstein b

-- | Calculate a value and the first derivative on the curve.
evalBezierDeriv b =
  let (px, py) = bezierToBernstein b
      px' = bernsteinDeriv px
      py' = bernsteinDeriv py
  in \t -> (Point (bernsteinEval px t) (bernsteinEval py t),
            Point (bernsteinEval px' t) (bernsteinEval py' t))

-- | Calculate a value and all derivatives on the curve.
evalBezierDerivs :: CubicBezier -> Double -> [Point]
evalBezierDerivs b t = zipWith Point (bernsteinEvalDerivs px t)
                       (bernsteinEvalDerivs py t)
  where (px, py) = bezierToBernstein b

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
findBezierInflection (CubicBezier (Point x0 y0) (Point x1 y1) (Point x2 y2) (Point x3 y3)) =
  filter bezierParam $ quadraticRoot a b c
    where
      ax = x1 - x0
      bx = x3 - x1 - ax
      cx = x3 - x2 - ax - 2*bx
      ay = y1 - y0
      by = y2 - y1 - ay
      cy = y3 - y2 - ay - 2*by
      a = bx*cy - by*cx
      b = ax*cy - ay*cx
      c = ax*by - ay*bx

-- | Find the cusps of a bezier.

-- find a cusp.  We look for points where the tangent is both horizontal
-- and vertical, which is only true for the zero vector.
findBezierCusp b = filter vertical $ bezierHoriz b
  where vertical = (== 0) . pointY . snd . evalBezierDeriv b

-- | Split a bezier curve into two curves.
splitBezier :: CubicBezier -> Double -> (CubicBezier, CubicBezier)
splitBezier (CubicBezier a b c d) t =
  let ab = interpolateVector a b t
      bc = interpolateVector b c t
      cd = interpolateVector c d t
      abbc = interpolateVector ab bc t
      bccd = interpolateVector bc cd t
      mid = interpolateVector abbc bccd t
  in (CubicBezier a ab abbc mid, CubicBezier mid bccd cd d)

-- | Return the subsegment between the two parameters.
bezierSubsegment :: CubicBezier -> Double -> Double -> CubicBezier
bezierSubsegment b t1 t2 
  | t1 > t2   = bezierSubsegment b t2 t1
  | otherwise = snd $ flip splitBezier (t1/t2) $
                fst $ splitBezier b t2

-- | Split a bezier curve into a list of beziers
-- The parameters should be in ascending order or
-- the result is unpredictable.
splitBezierN :: CubicBezier -> [Double] -> [CubicBezier]
splitBezierN b [] = [b]
splitBezierN b [t] = [a, b] where
  (a, b) = splitBezier b t
splitBezierN b (t:u:rest) =
  bezierSubsegment b 0 t :
  bezierSubsegment b t u :
  tail (splitBezierN b $ u:rest)

