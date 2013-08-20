module Geom2D.CubicBezier.Outline where
import Geom2D
import Geom2D.CubicBezier.Basic
import Geom2D.CubicBezier.Approximate
import Data.Map as M hiding (map) 
import Data.Function
import Data.List

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
approximateOffset cb@(CubicBezier p1 p2 p3 p4) dist tol =
  approximateCurveWithParams offsetCb points ts tol
  where tan1     = p2 ^-^ p1
        tan2     = p4 ^-^ p3
        offsetCb = CubicBezier
                   (offsetPoint dist p1 tan1)
                   (offsetPoint dist p2 tan1)
                   (offsetPoint dist p3 tan2)
                   (offsetPoint dist p4 tan2)
        points   = map (bezierOffsetPoint cb dist) ts
        ts = [i/16 | i <- [1..15]]

offsetSegment cb dist tol = offsetTol cb dist tol

-- subdivide the original curve and approximate the offset until
-- the maximum error is below tolerance
offsetTol cb dist tol
  | err <= tol = [cb_out]
  | otherwise     = offsetTol cb_l dist tol ++
                    offsetTol cb_r dist tol 
  where
    (cb_out, t, err) = approximateOffset cb dist tol
    (cb_l, cb_r) = splitBezier cb t

-- Keep a map from maxError to (t_min, t_err, curve, outline) for
-- each subsegment to keep track of the segment with the maximum
-- error.  This ensures a n log(n) execution time, rather than n^2
-- when a list is used.

offsetMax dist tol n segments
  | n <= 1 = error "minimum segments to offset is 1"
  | (n == 1) || (err < tol) = map fourth $
                              sortBy (compare `on` first) $
                              M.elems segments

    -- split the maximum curve in two and add the two segments to the map
  | otherwise = offsetMax dist tol (n-1) $
                M.insert err_l (t_min, t_err_l, cb_l, outline_l) $
                M.insert err_r (t_err, t_err_r, cb_r, outline_r) $
                newSegments
  where
    first (a,_,_,_)  = a
    fourth (_,_,_,a) = a
    ((err, (t_min, t_err, curve, outline)), newSegments) = M.deleteFindMax segments
    (cb_l, cb_r) = splitBezier curve t_err
    (outline_l, t_err_l, err_l)  = approximateOffset cb_l dist tol
    (outline_r, t_err_r, err_r)  = approximateOffset cb_r dist tol
    
offsetSegmentMax n cb dist tol =
  offsetMax dist tol n segments
  where segments              = M.singleton err (0, t_err, cb, outline)
        (outline, t_err, err) = approximateOffset cb dist tol