-- curve intersection using bezier clipping.
module Geom2D.CubicBezier.Intersection
       (bezierIntersection)
       where

import Geom2D
import Geom2D.CubicBezier.Basic
import Data.List
import Data.Function
import Data.Maybe
import Control.Monad

-- project points with equally spaced x values onto a vertical line to
-- find the outermost point.  Return the y values
project :: [Point] -> [Double]
project l = zipWith (flip (/)) [1..] $
            map (subtract (pointY $ head l) . pointY) $
            tail l
      
-- find the outermost point.  It doesn't look at the x values.
-- outerBy should be either maximumBy or minimumBy
findOuter :: Ord a => (((a, b) -> (a, b) -> Ordering) ->
                       [(Double, [Point])] -> (c, [Point]))
             -> [Point] -> [Point]
findOuter _ [] = []
findOuter _ [p] = [p]
findOuter outerBy l@(p:pts) = p : findOuter outerBy next
  where next = snd $ outerBy (compare `on` fst) $
               zip (project l) $ tails pts

-- take the y values and turn it in into a convex hull with upper en
-- lower points separated.
makeHull :: [Double] -> ([Point], [Point])
makeHull ds =
  let points = zipWith Point [i / 3 | i <- [0..3]] ds
  in (findOuter maximumBy points,
      findOuter minimumBy points)

-- test if the chords cross the fat line
-- use continuation passing style
testBelow :: Double -> [Point] -> Maybe Double -> Maybe Double
testBelow dmin [] _ = Nothing
testBelow dmin [_] _ = Nothing
testBelow dmin (p:q:rest) cont
  | pointY p >= dmin = cont
  | pointY p > pointY q = Nothing
  | pointY q < dmin = testBelow dmin (q:rest) cont
  | otherwise = Just $ intersectPt dmin p q

testBetween :: Double -> Point -> Maybe Double -> Maybe Double
testBetween dmax (Point x y) cont
  | (y <= dmax) = Just x
  | otherwise = cont

-- test if the chords cross the line y=dmax somewhere
testAbove :: Double -> [Point] -> Maybe Double
testAbove dmax [] = Nothing
testAbove dmax [_] = Nothing
testAbove dmax (p:q:rest)
  | pointY p < pointY q = Nothing
  | pointY q > dmax = testAbove dmax (q:rest)
  | otherwise = Just $ intersectPt dmax p q
                
-- find the x value where the line through the two points
-- intersect the line y=d
intersectPt d (Point x1 y1) (Point x2 y2) =
  x1 + (d  - y1) * (x2 - x1) / (y2 - y1)

-- make a hull and test over which interval the
-- curve is garuanteed to lie inside the fat line
chopHull dmin dmax ds = do
  let (upper, lower) = makeHull ds
  left_t <- testBelow dmin upper $
            testBetween dmax (head upper) $
            testAbove dmax lower
  right_t <- testBelow dmin (reverse upper) $
             testBetween dmax (last upper) $
             testAbove dmax (reverse lower)
  Just (left_t, right_t)

bezierClip p@(CubicBezier p0 p1 p2 p3) q@(CubicBezier q0 q1 q2 q3)
  tmin tmax umin umax prevClip eps reverse

  -- no intersection
  | chop_interval == Nothing = []

  -- not enough reduction, so split the curve in case we have
  -- multiple intersections
  | prevClip > 0.8 && newClip > 0.8 =
    if new_tmax - new_tmin > umax - umin -- split the longest segment
    then let
      (p1, p2) = splitBezier newP 0.5
      half_t = new_tmin + (new_tmax - new_tmin) / 2
      in bezierClip q p1 umin umax new_tmin half_t newClip eps (not reverse) ++
         bezierClip q p2 umin umax half_t new_tmax newClip eps (not reverse)
    else let
      (q1, q2) = splitBezier q 0.5
      half_t = umin + (umax - umin) / 2
      in bezierClip q1 newP umin half_t new_tmin new_tmax newClip eps (not reverse) ++
         bezierClip q2 newP half_t umax new_tmin new_tmax newClip eps (not reverse)

  -- within tolerance      
  | max (umax - umin) (new_tmax - new_tmin) < eps =
    if reverse
    then [ (umin + (umax-umin)/2,
            new_tmin + (new_tmax-new_tmin)/2) ]
    else [ (new_tmin + (new_tmax-new_tmin)/2,
            umin + (umax-umin)/2) ]

  -- iterate with the curves reversed.
  | otherwise =
      bezierClip q newP umin umax new_tmin new_tmax newClip eps (not reverse)

  where
    d = lineDistance (Line q0 q3)
    d1 = d q1
    d2 = d q2
    (dmin, dmax) | d1*d2 > 0 = (3/4 * minimum [0, d1, d2],
                                3/4 * maximum [0, d1, d2])
                 | otherwise = (4/9 * minimum [0, d1, d2],
                                4/9 * maximum [0, d1, d2])
    chop_interval = chopHull dmin dmax $
                    map d [p0, p1, p2, p3]
    Just (chop_tmin, chop_tmax) = chop_interval
    newP = bezierSubsegment p chop_tmin chop_tmax
    newClip = chop_tmax - chop_tmin
    new_tmin = tmax * chop_tmin + tmin * (1 - chop_tmin)
    new_tmax = tmax * chop_tmax + tmin * (1 - chop_tmax)

-- | Find the intersections between two Bezier curves within given
-- tolerance, using the Bezier Clip algorithm. Returns the parameters
-- for both curves.

bezierIntersection :: CubicBezier -> CubicBezier -> Double -> [(Double, Double)]
bezierIntersection p q eps = bezierClip p q 0 1 0 1 0 eps' False
  where
    eps' = min (bezierParamTolerance p eps) (bezierParamTolerance q eps)
