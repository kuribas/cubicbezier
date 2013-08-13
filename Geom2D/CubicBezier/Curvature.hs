module Geom2D.CubicBezier.Curvature
where
import Geom2D
import Geom2D.CubicBezier.Basic
import Geom2D.CubicBezier.Numeric
import Data.Maybe

-- | Compare a value with the radius of curvature.  This
-- is similar to @compare value radius@, but works with 
-- infinite radius.
compareRadius :: CubicBezier -- ^ The bezier curve
              -> Double      -- ^ The distance to the curve
              -> Double      -- ^ The parameter to the bezier
              -> Ordering    -- ^ The result
compareRadius b dist =
  let bd = evalBezierDerivs b
  in \t -> let
    (_: Point x' y': Point x'' y'': _) = bd t
    in compare ((x'^2 + y'^2)^3)
       ((dist * (x' * y'' - y' * x''))^2)

-- | Curvature of the Bezier curve.
curvature :: CubicBezier -> Double -> Double
curvature b =
  let bd = evalBezierDerivs b
  in \t -> let
    (_: Point x' y': Point x'' y'': _) = bd t
    in (x' * y'' - y' * x'') / (x'^2 + y'^2)**1.5

-- | Radius of curvature of the Bezier curve.  This
-- is the reciprocal of the curvature.
radiusCurvature :: CubicBezier -> Double -> Double
radiusCurvature b =
  let bd = evalBezierDerivs b
  in \t -> let
    (_: Point x' y': Point x'' y'': _) = bd t
    in (x'^2 + y'^2)**1.5 / (x' * y'' - y' * x'')

-- | Find a local maximum or minimum of the curvature
-- on the bezier.  It doesn't return inflection points, and may
-- not find a solution when there are multiple local
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
     else Just $ findRoot curvDeriv 1e-7 0 1

-- | Find points on the curve that have a certain radius of curvature.
-- This function may not find all points when there are more than two
-- such points, or the curve has inflection points.
findRadius b d = let
  bd = evalBezierDerivs b
  radiusSquare t =
    let (_: Point x' y': Point x'' y'': _ ) = bd t
    in (x'^2 + y'^2)^3 - (d*(x'*y'' - y'*x''))^2

  radiusBetween t1 t2 =
    let r1 = radiusCurvature b t1
        r2 = radiusCurvature b t2
        between | r1 < r2   = r1 <= d && d <= r2
                | otherwise = r2 <= d && d <= r1
    in if between
       then [findRoot radiusSquare 1e-7 t1 t2]
       else []

  in case curvatureExtremum b of
    Nothing -> radiusBetween 0 1
    Just ex -> radiusBetween 0 ex ++ radiusBetween ex 1
