module Geom2D.CubicBezier.Intersection (bezierIntersection)
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
-- The type is more scary than the function.  outerBy should
-- be maximumBy or minimumBy
findOuter :: Ord a => (((a, b) -> (a, b) -> Ordering) ->
                       [(Double, [Point])] -> (c, [Point]))
             -> [Point] -> [Point]
findOuter outerBy [] = []
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

-- test if the chords cross the line y=dmin somewhere
testBelow dmin [] = Nothing
testBelow dmin [_] = Nothing
testBelow dmin (p:q:rest)
  | pointY p >= dmin = Nothing
  | pointY p > pointY q = Nothing
  | pointY q > dmin = testBelow dmin (q:rest)
  | otherwise = Just $ intersectPt dmin p q

-- test if the chords cross the line y=dmax somewhere
testAbove dmax [] = Nothing
testAbove dmax [_] = Nothing
testAbove dmax (p:q:rest)
  | pointY p <= dmax = Nothing
  | pointY p < pointY q = Nothing
  | pointY q < dmax = testBelow dmax (q:rest)
  | otherwise = Just $ intersectPt dmax p q

-- find the x value where the line through the two points
-- intersect the line y=d
intersectPt d (Point x1 y1) (Point x2 y2) =
  x1 + (d  - y1) * (x2 - x1) / (y2 - y1)

chopHull dmin dmax ds = (left_t, right_t)
  where
    (upper, lower) = makeHull ds
    left_t  = fromMaybe 0 $
              testBelow dmin lower `mplus`
              testAbove dmax upper
    right_t = fromMaybe 1 $
              (testBelow dmin $ reverse lower) `mplus`
              (testAbove dmax $ reverse upper)

fatLineDistance (Point x1 y1) (Point x2 y2) = \(Point x y) -> a*x + b*y + c
  where a = a' / d
        b = b' / d
        c = -(b' / a') / d
        a' = y1 - y2
        b' = x2 - x1
        d = sqrt(a'^2 + b'^2)

bezierClip p@(CubicBezier p0 p1 p2 p3) q@(CubicBezier q0 q1 q2 q3)
  tmin tmax umin umax prevClip eps reverse

  -- not enough reduction, so split the curve in case we have
  -- multiple intersections
  | prevClip > 0.8 && newClip > 0.8 =
    if new_tmax - new_tmin > umax - umin -- split the longest segment
    then let
      (p1, p2) = splitBezier newP 0.5
      half_t = (new_tmax-new_tmin) / 2
      in bezierClip p1 q new_tmin half_t umin umax newClip eps reverse ++
         bezierClip p2 q half_t tmax umin umax newClip eps reverse
    else let
      (q1, q2) = splitBezier q 0.5
      half_t = (umax - umin) / 2
      in bezierClip newP q1 new_tmin new_tmax umin half_t newClip eps reverse ++
         bezierClip newP q2 new_tmin new_tmax half_t umax newClip eps reverse

  -- within tolerance      
  | max (umax - umin) (new_tmax - new_tmin) < eps =
    if reverse
    then [(umax - umin)/2, (new_tmax/new_tmin)/2]
    else [(new_tmax - new_tmin)/2 , (umax - umin)/2]

  -- iterate with the curves reversed.
  | otherwise =
      bezierClip q p umin umax tmin tmax newClip eps (not reverse)

  where
    d = fatLineDistance q0 q3
    d1 = d q1
    d2 = d q2
    (dmin, dmax) | d1*d2 > 0 = (3/4 * minimum [0, d1, d2],
                                3/4 * maximum [0, d1, d2])
                 | otherwise = (4/9 * minimum [0, d1, d2],
                                4/9 * maximum [0, d1, d2])
    (chop_tmin, chop_tmax) = chopHull dmin dmax $ map d [p0, p1, p2, p3]
    newP = fst $ flip splitBezier chop_tmin $
           snd $ splitBezier p chop_tmax
    newClip = chop_tmax - chop_tmin
    new_tmin = chop_tmin * (tmax - tmin) + tmin
    new_tmax = chop_tmax * (tmax - tmin) + tmin

-- | Find the intersections between two Bezier curves to within tolerance eps.
bezierIntersection p q eps = bezierClip p q 0 1 0 1 1 eps False