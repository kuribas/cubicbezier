module Geom2D.CubicBezier.Curvature
       (curvature, radiusOfCurvature, curvatureExtrema, findRadius)
where
import Geom2D
import Geom2D.CubicBezier.Basic
import Geom2D.CubicBezier.Intersection
import Data.Maybe
import Data.List

-- | Curvature of the Bezier curve.
curvature b t
  | t == 0 = curvature' b
  | t == 1 = curvature' $ reorient b
  | t < 0.5 = curvature' $ snd $ splitBezier b t
  | otherwise = curvature' $ reorient $ fst $ splitBezier b t

curvature' b@(CubicBezier c0 c1 c2 c3) = 2/3 * b/a^3
  where 
    a = vectorDistance c1 c0
    b = (c1^-^c0) `vectorCross` (c2^-^c1)

-- | Radius of curvature of the Bezier curve.  This
-- is the reciprocal of the curvature.
radiusOfCurvature :: CubicBezier -> Double -> Double
radiusOfCurvature b t = 1 / curvature b t

extrema (CubicBezier p0 p1 p2 p3) =
  let bez = [p0, p1, p2, p3]
      x' = genericDeriv $ map pointX bez
      y' = genericDeriv $ map pointY bez
      x'' = genericDeriv x'
      y'' = genericDeriv y'
      x''' = genericDeriv x''
      y''' = genericDeriv y''
  in -- (y'^2 + x'^2) * (x'*y''' - y'*x''') -
     -- 3 * (x'*y'' - y'*x'') * (y'*y'' + x'*x'')
   (y'~*y' ~+ x'~*x') ~* (x'~*y''' ~- y'~*x''') ~-
   3 *~ (x'~*y'' ~- y'~*x'') ~* (y'~*y'' ~+ x'~*x'')

-- | Find extrema of the curvature, but not inflection points.
curvatureExtrema :: CubicBezier -> Double -> [Double]
curvatureExtrema b eps = bezierFindRoot (extrema b) 0 1 $
                         bezierParamTolerance b eps

radiusSquareEq (CubicBezier p0 p1 p2 p3) d =
  let bez = [p0, p1, p2, p3]
      x' = genericDeriv $ map pointX bez
      y' = genericDeriv $ map pointY bez
      x'' = genericDeriv x'
      y'' = genericDeriv y'
      a =  x'~*x' ~+  y'~*y'
      b =  x'~*y'' ~-  x''~*y'
  in (a~*a~*a) ~- (d*d) *~ b~*b

-- | Find points on the curve that have a certain radius of curvature.
findRadius :: CubicBezier -> Double -> Double -> [Double]
findRadius b d eps = bezierFindRoot (radiusSquareEq b d) 0 1 $
                     bezierParamTolerance b eps
