module Geom2D.CubicBezier.Outline where
import Geom2D
import Geom2D.CubicBezier.Basic
import Geom2D.CubicBezier.Approximate

data Side = LeftSide | RightSide

offsetPoint dist LeftSide start tangent =
  addVector start
  (rotateVector90Left $ scaleVector dist $ unitVector tangent)

offsetPoint dist RightSide start tangent =
  addVector start
  (rotateVector90Right $ scaleVector dist $ unitVector tangent)

bezierOffsetPoint cb side dist t =
  uncurry (offsetPoint dist side) $
  evalBezierDeriv cb t

bezierOffset cb@(CubicBezier p1 p2 p3 p4) side dist =
  approximateCurveWithParams offsetCb points ts
  where offsetCb = CubicBezier
                   (offsetPoint dist side p1 (subVector p2 p1))
                   (offsetPoint dist side p2 (subVector p2 p1))
                   (offsetPoint dist side p3 (subVector p4 p3))
                   (offsetPoint dist side p4 (subVector p4 p3))
        points = map (bezierOffsetPoint cb side dist) ts
        ts = [0.1,0.2..0.9]
