-- | Intersection routines using Bezier Clipping.  Provides also functions for finding the roots of onedimensional bezier curves.  This can be used as a general polynomial root solver by converting from the power basis to the bernstein basis.
module Geom2D.CubicBezier.Intersection
       (bezierIntersection, bezierLineIntersections, genericSubsegment,
        (~*), (*~), (~+), (~-), degreeElevate, genericSplit, genericEval,
        genericDeriv, bezierFindRoot)
       where
import Geom2D
import Geom2D.CubicBezier.Basic
import Data.List
import Data.Function
import Data.Maybe
import Control.Monad

infixl 7 ~*, *~
infixl 6 ~+, ~-

-- find the convex hull by comparing the angles of the vectors with
-- the cross product and backtracking if necessary.
findOuter' upper dir p1 l@(p2:rest)
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


-- | Return the subsegment between the two parameters for a 1D bezier
-- curve of any degree.
genericSubsegment :: [Double] -> Double -> Double -> [Double]
genericSubsegment b t1 t2 
  | t1 > t2   = genericSubsegment b t2 t1
  | otherwise = snd $ flip genericSplit (t1/t2) $
                fst $ genericSplit b t2

-- multiply two bezier curves
-- control point i from the product of beziers P * Q
-- is sum (P_j * Q_k) where j + k = i+1

-- | Multiply two 1D bezier curves of any degree.  The final degree
-- will be the sum of either degrees.  This operation takes O((n+m)^2)
-- with n and m the degree of the beziers.

(~*) :: [Double] -> [Double] -> [Double]
a ~* b = zipWith (flip (/)) (binCoeff (la + lb)) $
                 init $ map sum $
                 zipWith (zipWith (*)) (repeat a') (down b') ++
                 zipWith (zipWith (*)) (tail $ tails a') (repeat $ reverse b')
  where down l = tail $ scanl (flip (:)) [] l -- [[1], [2, 1], [3, 2, 1], ...
        a' = zipWith (*) a (binCoeff la)
        b' = zipWith (*) b (binCoeff lb)
        la = length a - 1
        lb = length b - 1

degreeElevate' b _ 0 = b
degreeElevate' l d times =
  degreeElevate' (head l:inner l 1) (d+1) (times-1)
  where
    inner [a] _ = [a]
    inner (a:b:rest) i =
      (i*a/fromIntegral d + b*(1 - i/fromIntegral d))
      : inner (b:rest) (i+1)

-- find the binomial coefficients of degree n.
binCoeff :: Int -> [Double]
binCoeff n = map fromIntegral $
             scanl (\x m -> x * (n-m+1) `quot` m) 1 [1..n]

-- | Degree elevate a 1D bezier curve of any degree.
degreeElevate :: (Eq a, Fractional a1, Num a) => [a1] -> a -> [a1]
degreeElevate l times = degreeElevate' l (length l) times

-- | Evaluate the 1D bezier curve.
genericEval b t = last $ fst $ genericSplit b t

-- | Find the derivative of a 1D bezier curve.
genericDeriv b = map (*n) $ zipWith (-) (tail b) b
  where n = fromIntegral $ length b - 1

-- | Split a 1D bezier curve of any degree.
genericSplit b t = (map head controls, reverse $ map last controls)
  where
    interp a b = (1-t)*a + t*b
    terp [_] = []
    terp l = let ctrs = zipWith interp l (tail l)
             in ctrs : terp ctrs
    controls = b:terp b

-- | Sum two 1D bezier curves of any degree.  The final degree will be
-- the maximum of either degrees.
(~+) :: [Double] -> [Double] -> [Double]
a ~+ b
  | la < lb = zipWith (+) (degreeElevate a (lb-la)) b
  | la > lb = zipWith (+) (degreeElevate b (la-lb)) a
  | otherwise = zipWith (+) a b
  where la = length a
        lb = length b

-- | Difference two 1D bezier curves of any degree.  The final degree will be
-- the maximum of either degrees.
(~-) :: [Double] -> [Double] -> [Double]
a ~- b
  | la < lb = zipWith (-) (degreeElevate a (lb-la)) b
  | la > lb = zipWith (-) a (degreeElevate b (la-lb))
  | otherwise = zipWith (-) a b
  where la = length a
        lb = length b

-- | Scale a 1D bezier curve of any degree by a constant.
(*~) :: Double -> [Double] -> [Double]
(*~) a = map (*a)

------------------------ Line intersection -------------------------------------
-- Clipping a line uses a simplified version of the Bezier Clip algorithm,
-- and uses the (thin) line itself instead of the fat line.

-- | Find the zero of a 1D bezier curve of any degree.  Note that this
-- can be used as a generic polynomial root solver by converting from
-- the power basis to the bernstein basis.
bezierFindRoot p tmin tmax eps
  -- no intersection
  | chop_interval == Nothing = []

  -- not enough reduction, so split the curve in case we have
  -- multiple intersections
  | clip > 0.8 =
    let (p1, p2) = genericSplit newP 0.5
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
    chop_interval = chopHull 0 0 p
    Just (chop_tmin, chop_tmax) = chop_interval
    newP = genericSubsegment p chop_tmin chop_tmax
    clip = chop_tmax - chop_tmin
    new_tmin = tmax * chop_tmin + tmin * (1 - chop_tmin)
    new_tmax = tmax * chop_tmax + tmin * (1 - chop_tmax)

-- | Find the intersections of the curve with a line.

-- Apply a transformation to the bezier that maps the line onto the
-- X-axis.  Then we only need to test the Y-values for a zero.
bezierLineIntersections :: CubicBezier -> Line -> Double -> [Double]
bezierLineIntersections b (Line p q) eps =
  bezierFindRoot (map pointY [p0, p1, p2, p3]) 0 1 $
  bezierParamTolerance b eps
  where (CubicBezier p0 p1 p2 p3) = 
          (fromJust $ inverse $ translate p $* rotateVec (q ^-^ p)) $* b