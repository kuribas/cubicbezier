module CubicBezier.Outline where
import CubicBezier.Geom2D
import CubicBezier.Bezier
import Data.List

data Side = LeftSide | RightSide

offsetPoint dist Left start tangent =
  addVector start
  (rotateVector90Left $ scaleVector dist $ unitVector tangent)

offsetPoint dist Right start tangent =
  addVector start
  (rotateVector90Right $ scaleVector dist $ unitVector tangent)

bezierOffsetPoint cb side dist t =
  uncurry (offsetPoint dist side) $
  evalBezierDeriv cb t

bezierOffset cb@(CubicBezier p1 p2 p3 p4) side dist =
  approximateCurveWithTs offsetCb points ts
  where offsetCb = CubicBezier
                   (offsetPoint dist side p1 (subtractVector p2 p1))
                   (offsetPoint dist side p2 (subtractVector p2 p1))
                   (offsetPoint dist side p3 (subtractVector p4 p3))
                   (offsetPoint dist side p4 (subtractVector p4 p3))
        points = map (bezierOffsetPoint cb side dist) ts
        ts = [0.1,0.2..0.9]

