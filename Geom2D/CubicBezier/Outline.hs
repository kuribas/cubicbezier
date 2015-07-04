-- | Offsetting bezier curves and stroking curves.

module Geom2D.CubicBezier.Outline
       (bezierOffset)
       where
import Geom2D
import Geom2D.CubicBezier.Basic
import Geom2D.CubicBezier.Approximate

offsetPoint :: (Floating a) =>  a -> Point a -> Point a -> Point a
offsetPoint dist start tangent =
  start ^+^ (rotate90L $* dist *^ normVector tangent)

bezierOffsetPoint :: CubicBezier Double -> Double -> Double -> (DPoint, DPoint)
bezierOffsetPoint cb dist t = (offsetPoint dist p p', p')
  where (p, p') = evalBezierDeriv cb t

-- | Calculate an offset path from the bezier curve to within
-- tolerance.  If the distance is positive offset to the left,
-- otherwise to the right. A smaller tolerance may require more bezier
-- curves in the path to approximate the offset curve
bezierOffset :: CubicBezier Double -- ^ The curve
             -> Double      -- ^ Offset distance.
             -> Maybe Int   -- ^ maximum subcurves
             -> Double      -- ^ Tolerance.
             -> [CubicBezier Double]        -- ^ The offset curve
bezierOffset cb dist (Just m) tol =
  approximatePathMax m (bezierOffsetPoint cb dist) 15 tol 0 1

bezierOffset cb dist Nothing tol =
  approximatePath (bezierOffsetPoint cb dist) 15 tol 0 1
