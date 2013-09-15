module Geom2D.CubicBezier.Curvature
       (curvature, radiusOfCurvature, curvatureExtrema, findRadius)
where
import Geom2D
import Geom2D.CubicBezier.Basic
import Geom2D.CubicBezier.Intersection
import Math.BernsteinPoly

-- | Curvature of the Bezier curve.  A negative curvature means the curve
-- curves to the right.
curvature :: CubicBezier -> Double -> Double
curvature b t
  | t == 0 = curvature' b
  | t == 1 = negate $ curvature' $ reorient b
  | t < 0.5 = curvature' $ snd $ splitBezier b t
  | otherwise = negate $ curvature' $ reorient $ fst $ splitBezier b t

curvature' :: CubicBezier -> Double
curvature' (CubicBezier c0 c1 c2 _c3) = 2/3 * b/a^(3::Int)
  where 
    a = vectorDistance c1 c0
    b = (c1^-^c0) `vectorCross` (c2^-^c1)

-- | Radius of curvature of the Bezier curve.  This
-- is the reciprocal of the curvature.
radiusOfCurvature :: CubicBezier -> Double -> Double
radiusOfCurvature b t = 1 / curvature b t

extrema :: CubicBezier -> BernsteinPoly
extrema bez = 
  let (x, y) = bezierToBernstein bez
      x' = bernsteinDeriv x
      y' = bernsteinDeriv y
      x'' = bernsteinDeriv x'
      y'' = bernsteinDeriv y'
      x''' = bernsteinDeriv x''
      y''' = bernsteinDeriv y''
  in (y'~*y' ~+ x'~*x') ~* (x'~*y''' ~- y'~*x''') ~-
     3 *~ (x'~*y'' ~- y'~*x'') ~* (y'~*y'' ~+ x'~*x'')

-- | Find extrema of the curvature, but not inflection points.
curvatureExtrema :: CubicBezier -> Double -> [Double]
curvatureExtrema b eps
  | colinear b eps = []
  | otherwise = bezierFindRoot (extrema b) 0 1 $
                bezierParamTolerance b eps

radiusSquareEq :: CubicBezier -> Double -> BernsteinPoly
radiusSquareEq bez d =
  let (x, y) = bezierToBernstein bez
      x' = bernsteinDeriv x
      y' = bernsteinDeriv y
      x'' = bernsteinDeriv x'
      y'' = bernsteinDeriv y'
      a =  x'~*x' ~+  y'~*y'
      b =  x'~*y'' ~-  x''~*y'
  in (a~*a~*a) ~- (d*d) *~ b~*b

-- | Find points on the curve that have a certain radius of curvature.
-- Values to the left and to the right of the curve are returned.
findRadius :: CubicBezier  -- ^ the curve
           -> Double       -- ^ distance
           -> Double       -- ^ tolerance
           -> [Double]     -- ^ result
findRadius b d eps
  | colinear b eps = []  -- either empty or a huge list of t's
  | otherwise = bezierFindRoot (radiusSquareEq b d) 0 1 $
                bezierParamTolerance b eps
