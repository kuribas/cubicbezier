module Geom2D.CubicBezier.Outline where
import Geom2D
import Geom2D.CubicBezier.Basic
import Geom2D.CubicBezier.Approximate

offsetPoint :: Double -> Point -> Point -> Point
offsetPoint dist start tangent =
  start ^+^ (rotateVector90Left $ dist *^ normVector tangent)

bezierOffsetPoint :: CubicBezier -> Double -> Double -> Point
bezierOffsetPoint cb dist t =
  uncurry (offsetPoint dist) $
  evalBezierDeriv cb t

-- Approximate the bezier curve offset by dist.  A positive value
-- means to the left, a negative to the right.
approximateOffset :: CubicBezier -> Double -> Double -> (CubicBezier, Double, Double)
approximateOffset cb@(CubicBezier p1 p2 p3 p4) dist eps =
  approximateCurveWithParams offsetCb points ts eps
  where tan1 = p2 ^-^ p1
        tan2 = p4 ^-^ p3
        offsetCb = CubicBezier
                   (offsetPoint dist p1 tan1)
                   (offsetPoint dist p2 tan1)
                   (offsetPoint dist p3 tan2)
                   (offsetPoint dist p4 tan2)
        points = map (bezierOffsetPoint cb dist) ts
        ts = [0.1,0.2..0.9]
