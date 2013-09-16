{-# LANGUAGE BangPatterns #-}
-- | Intersection routines using Bezier Clipping.  Provides also functions for finding the roots of onedimensional bezier curves.  This can be used as a general polynomial root solver by converting from the power basis to the bernstein basis.
module Geom2D.CubicBezier.Intersection
       (bezierIntersection, bezierLineIntersections, bezierFindRoot, closest)
       where
import Geom2D
import Geom2D.CubicBezier.Basic
import Math.BernsteinPoly
import Data.Maybe


-- find the convex hull by comparing the angles of the vectors with
-- the cross product and backtracking if necessary.
findOuter' :: Bool -> Point -> Point -> [Point] -> Either [Point] [Point]
findOuter' !upper !dir !p1 l@(p2:rest)
  -- backtrack if the direction is outward
  | if upper
    then dir `vectorCross` (p2^-^p1) > 0 -- left turn
    else dir `vectorCross` (p2^-^p1) < 0 = Left l
  -- succeed
  | otherwise = case findOuter' upper (p2^-^p1) p2 rest of
    Left m -> findOuter' upper dir p1 m
    Right m -> Right (p1:m)

findOuter' _ _ p1 p = Right (p1:p)

-- find the outermost point.  It doesn't look at the x values.
findOuter :: Bool -> [Point] -> [Point]
findOuter upper (p1:p2:rest) =
  case findOuter' upper (p2^-^p1) p2 rest of
    Right l -> p1:l
    Left l -> findOuter upper (p1:l)
findOuter _ l = l    

-- take the y values and turn it in into a convex hull with upper en
-- lower points separated.
makeHull :: [Double] -> ([Point], [Point])
makeHull ds =
  let n      = fromIntegral $ length ds - 1
      points = zipWith Point [i/n | i <- [0..n]] ds
  in (findOuter True points,
      findOuter False points)

-- test if the chords cross the fat line
-- return the continuation if above the line
testBelow :: Double -> [Point] -> Maybe Double -> Maybe Double
testBelow _    [] _ = Nothing
testBelow _    [_] _ = Nothing
testBelow !dmin (p:q:rest) cont
  | pointY p >= dmin = cont
  | pointY p > pointY q = Nothing
  | pointY q < dmin = testBelow dmin (q:rest) cont
  | otherwise = Just $ intersectPt dmin p q

testBetween :: Double -> Point -> Maybe Double -> Maybe Double
testBetween !dmax (Point !x !y) cont
  | y <= dmax = Just x
  | otherwise = cont

-- test if the chords cross the line y=dmax somewhere
testAbove :: Double -> [Point] -> Maybe Double
testAbove _    [] = Nothing
testAbove _    [_] = Nothing
testAbove dmax (p:q:rest)
  | pointY p < pointY q = Nothing
  | pointY q > dmax = testAbove dmax (q:rest)
  | otherwise = Just $ intersectPt dmax p q

-- find the x value where the line through the two points
-- intersect the line y=d
intersectPt :: Double -> Point -> Point -> Double
intersectPt d (Point x1 y1) (Point x2 y2) =
  x1 + (d  - y1) * (x2 - x1) / (y2 - y1)

-- make a hull and test over which interval the
-- curve is garuanteed to lie inside the fat line
chopHull :: Double -> Double -> [Double] -> Maybe (Double, Double)
chopHull !dmin !dmax ds = do
  let (upper, lower) = makeHull ds
  left_t <- testBelow dmin upper $
            testBetween dmax (head upper) $
            testAbove dmax lower
  right_t <- testBelow dmin (reverse upper) $
             testBetween dmax (last upper) $
             testAbove dmax (reverse lower)
  Just (left_t, right_t)

bezierClip :: CubicBezier -> CubicBezier -> Double -> Double
           -> Double -> Double -> Double -> Double -> Bool
           -> [(Double, Double)]
bezierClip p@(CubicBezier !p0 !p1 !p2 !p3) q@(CubicBezier !q0 !q1 !q2 !q3)
  tmin tmax umin umax prevClip eps revCurves

  -- no intersection
  | isNothing chop_interval = []

  -- not enough reduction, so split the curve in case we have
  -- multiple intersections
  | prevClip > 0.8 && newClip > 0.8 =
    if new_tmax - new_tmin > umax - umin -- split the longest segment
    then let
      (pl, pr) = splitBezier newP 0.5
      half_t = new_tmin + (new_tmax - new_tmin) / 2
      in bezierClip q pl umin umax new_tmin half_t newClip eps (not revCurves) ++
         bezierClip q pr umin umax half_t new_tmax newClip eps (not revCurves)
    else let
      (ql, qr) = splitBezier q 0.5
      half_t = umin + (umax - umin) / 2
      in bezierClip ql newP umin half_t new_tmin new_tmax newClip eps (not revCurves) ++
         bezierClip qr newP half_t umax new_tmin new_tmax newClip eps (not revCurves)

  -- within tolerance      
  | max (umax - umin) (new_tmax - new_tmin) < eps =
    if revCurves
    then [ (umin + (umax-umin)/2,
            new_tmin + (new_tmax-new_tmin)/2) ]
    else [ (new_tmin + (new_tmax-new_tmin)/2,
            umin + (umax-umin)/2) ]

  -- iterate with the curves reversed.
  | otherwise =
      bezierClip q newP umin umax new_tmin new_tmax newClip eps (not revCurves)

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

------------------------ Line intersection -------------------------------------
-- Clipping a line uses a simplified version of the Bezier Clip algorithm,
-- and uses the (thin) line itself instead of the fat line.

-- | Find the zero of a 1D bezier curve of any degree.  Note that this
-- can be used as a bernstein polynomial root solver by converting from
-- the power basis to the bernstein basis.
bezierFindRoot :: BernsteinPoly -- ^ the bernstein coefficients of the polynomial
               -> Double  -- ^ The lower bound of the interval 
               -> Double  -- ^ The upper bound of the interval
               -> Double  -- ^ The accuracy
               -> [Double] -- ^ The roots found
bezierFindRoot p tmin tmax eps
  -- no intersection
  | chop_interval == Nothing = []

  -- not enough reduction, so split the curve in case we have
  -- multiple intersections
  | clip > 0.8 =
    let (p1, p2) = bernsteinSplit newP 0.5
        half_t = new_tmin + (new_tmax - new_tmin) / 2
    in bezierFindRoot p1 new_tmin half_t eps ++
       bezierFindRoot p2 half_t new_tmax eps

  -- within tolerance
  | new_tmax - new_tmin < eps =
      [new_tmin + (new_tmax-new_tmin)/2]

      -- iterate
  | otherwise =
        bezierFindRoot newP new_tmin new_tmax eps

  where
    chop_interval = chopHull 0 0 (bernsteinCoeffs p)
    Just (chop_tmin, chop_tmax) = chop_interval
    newP = bernsteinSubsegment p chop_tmin chop_tmax
    clip = chop_tmax - chop_tmin
    new_tmin = tmax * chop_tmin + tmin * (1 - chop_tmin)
    new_tmax = tmax * chop_tmax + tmin * (1 - chop_tmax)

-- | Find the intersections of the curve with a line.

-- Apply a transformation to the bezier that maps the line onto the
-- X-axis.  Then we only need to test the Y-values for a zero.
bezierLineIntersections :: CubicBezier -> Line -> Double -> [Double]
bezierLineIntersections b (Line p q) eps =
  bezierFindRoot (listToBernstein $ map pointY [p0, p1, p2, p3]) 0 1 $
  bezierParamTolerance b eps
  where (CubicBezier p0 p1 p2 p3) = 
          fromJust (inverse $ translate p $* rotateVec (q ^-^ p)) $* b

-- | Find the closest value(s) on the bezier to the given point, within tolerance.
closest :: CubicBezier -> Point -> Double -> [Double]
closest cb (Point px py) eps = bezierFindRoot poly 0 1 eps
  where
    (bx, by) = bezierToBernstein cb
    bx' = bernsteinDeriv bx
    by' = bernsteinDeriv by
    poly = (bx ~- listToBernstein [px, px, px, px]) ~* bx' ~+
           (by ~- listToBernstein [py, py, py, py]) ~* by'

