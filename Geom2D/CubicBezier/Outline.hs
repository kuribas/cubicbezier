-- | Offsetting bezier curves and stroking curves.

module Geom2D.CubicBezier.Outline
       (bezierOffset, bezierOffsetMax)
       where
import Geom2D
import Geom2D.CubicBezier.Basic
import Geom2D.CubicBezier.Approximate
import Geom2D.CubicBezier.Curvature

offsetPoint :: Double -> Point -> Point -> Point
offsetPoint dist start tangent =
  start ^+^ (rotate90L $* dist *^ normVector tangent)

bezierOffsetPoint :: CubicBezier -> Double -> Double -> (Point, Point)
bezierOffsetPoint cb dist t = (offsetPoint dist p p', p')
  where (p, p') = evalBezierDeriv cb t

-- Approximate the bezier curve offset by dist.  A positive value
-- means to the left, a negative to the right.
offsetSegment :: Double -> Double -> CubicBezier -> [CubicBezier]
offsetSegment dist tol cb =
  approximatePath (bezierOffsetPoint cb dist) 15 tol 0 1

offsetSegmentMax :: Int -> Double -> Double -> CubicBezier -> [CubicBezier]
offsetSegmentMax m dist tol cb =
  approximatePathMax m (bezierOffsetPoint cb dist) 15 tol 0 1

-- | Calculate an offset path from the bezier curve to within
-- tolerance.  If the distance is positive offset to the left,
-- otherwise to the right. A smaller tolerance may require more bezier
-- curves in the path to approximate the offset curve
bezierOffset :: CubicBezier -- ^ The curve
             -> Double      -- ^ Offset distance.
             -> Double      -- ^ Tolerance.
             -> [CubicBezier]        -- ^ The offset curve
bezierOffset cb dist tol =
  --Path $ map BezierSegment $
  concatMap (offsetSegment dist tol) $
  splitBezierN cb $
  findRadius cb dist tol

-- | Like bezierOffset, but limit the number of subpaths for each
-- smooth subsegment.  The number should not be smaller than one.
bezierOffsetMax :: Int -> CubicBezier -> Double -> Double -> [CubicBezier]
bezierOffsetMax n cb dist tol =
  -- Path $ map BezierSegment $
  concatMap (offsetSegmentMax n dist tol) $
  splitBezierN cb $
  findRadius cb dist tol
