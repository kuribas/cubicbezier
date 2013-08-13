module Geom2D.CubicBezier.Curvature
where
import Geom2D
import Geom2D.CubicBezier.Basic
import Geom2D.CubicBezier.Numeric
import Data.Maybe

-- | Curvature of the Bezier curve.
curvature :: CubicBezier -> Double -> Double
curvature b =
  let bd = evalBezierDerivs b
  in \t -> let
    (_: Point x' y': Point x'' y'': _) = bd t
    in (x' * y'' - y' * x'') / (x'^2 + y'^2)**1.5

-- | Radius of curvature of the Bezier curve.  This
-- is the reciprocal of the curvature.
radiusOfCurvature :: CubicBezier -> Double -> Double
radiusOfCurvature b =
  let bd = evalBezierDerivs b
  in \t -> let
    (_: Point x' y': Point x'' y'': _) = bd t
    in (x'^2 + y'^2)**1.5 / (x' * y'' - y' * x'')

-- | Find a local maximum or minimum of the curvature
-- on the bezier.  It may not return all solutions
-- when there are inflection points or multiple local
-- extrema.

-- calculate the extrama of the curvature by
-- finding the root of the derivative of
-- the square of the curvature.
-- Inflection points ((x'*y'' - y'*x'') == 0)
-- are not handled.
curvatureExtremum b =
  let bd = evalBezierDerivs b
      curvDeriv t = let
        (_: Point x' y': Point x'' y'': Point x''' y''': _) = bd t
        in (y'^2 + x'^2) * (x'*y''' - y'*x''') - 3 * (x'*y'' - y'*x'') * (y'*y'' + x'*x'')
  in if signum (curvDeriv 0) == signum (curvDeriv 1)
     then Nothing -- we cannot use the root finder
     else Just $ findRoot curvDeriv 1e-8 0 1

if' a b c = if a then b else c

-- | Find points on the curve that have a certain radius of curvature.
-- The curve shouldn't have inflection points.  This function may not
-- find all points when there are more than two such points.
findRadius b d = let
  bd = evalBezierDerivs b
  radiusSquare t =
    let (_: Point x' y': Point x'' y'': _ ) = bd t
    in (x'^2 + y'^2)^3 - (d*(x'*y'' - y'*x''))^2

  radiusBetween t1 t2 =
    let r1 = radiusOfCurvature b t1
        r2 = radiusOfCurvature b t2
    in if' (r1 == d) [r1] $
       if' (r2 == d) [r2] $
       if' ((r1 * r2 >= 0) &&
            (if r1 < r2 then r1 < d && d < r2 else r2 < d && d < r1))
            [findRoot radiusSquare 1e-8 t1 t2]
       []
  in case curvatureExtremum b of
    Nothing -> radiusBetween 0 1
    Just ex -> radiusBetween 0 ex ++ radiusBetween ex 1
