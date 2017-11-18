{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
module Geom2D.CubicBezier.Stroke
--       (penCircle, pathToPen, penStrokeOpen, penStrokeClosed, Pen,
--        bezierOffset)
       where
import Geom2D
import Geom2D.CubicBezier
import Data.Monoid
import qualified Data.Vector as V

data Pen = PenEllipse (Transform Double) (Transform  Double) (Transform Double)
         | PenPath (FindDir PenSegment)

data PenSegment = PenCorner !DPoint
                | PenCurve !(CubicBezier Double)

data FindDir b = FindDir DPoint (V.Vector (DPoint, b)) (V.Vector (DPoint, b))


clockwise :: DPoint -> DPoint -> Bool
clockwise v1 v2 =
  vc > 0 ||
  (vc == 0 && signum (pointX v1) == signum (pointX v2) &&
   signum (pointY v1) == signum (pointY v2))
  where vc = vectorCross v2 v1

square :: Double -> Double
square x = x*x

-- angle between two vectors is small
closeDirs :: DPoint -> DPoint -> Bool
closeDirs v1 v2 =
  (1 - (abs dp * dp / (vectorMagSquare v1*vectorMagSquare v2))) < 1e-5
  where dp = v1 ^.^ v2

findDir :: FindDir b -> DPoint -> b
findDir (FindDir dir1 lft rt) dir
  | clockwise dir1 dir =
      findDirPart lft dir 0 (V.length lft-1)
  | otherwise =
      findDirPart rt dir 0 (V.length rt-1)

-- binary search on the directions
findDirPart :: (V.Vector (DPoint, b)) -> DPoint -> Int -> Int -> b
findDirPart part dir min_ max_
  | min_ == max_ = snd $ part V.! min_
  | clockwise midDir dir =
      findDirPart part dir min_ (mid-1)
  | otherwise =
      findDirPart part dir mid max_
  where
    midDir = fst $ part V.! mid
    mid = min_ + (max_ - min_ + 1) `quot` 2

-- | A circular pen with unit radius.
penCircle :: Pen
penCircle = PenEllipse idTrans rotate90L rotate90R

-- | Create a pen from a path.  For predictable results the path
-- should be convex.
pathToPen :: ClosedPath Double -> Pen
pathToPen (ClosedPath []) =
  PenPath $ FindDir (Point 1.0 0.0)
  (V.fromList [(Point 1.0 0, PenCorner $ Point 0 0)])
  (V.fromList [(Point (-1.0) 0, PenCorner $ Point 0 0)])
pathToPen (ClosedPath nodes) =
  PenPath $ splitPartitions $ pathPartitions $ nodes ++ [head nodes]

pathPartitions :: [(DPoint, PathJoin Double)] -> [(DPoint, PenSegment)]
pathPartitions ((p, JoinLine):tl@((q, _):_)) =
  (q ^-^ p, PenCorner q) :
  pathPartitions tl

pathPartitions ((p1, JoinCurve p2 p3):tl@((p4, _):_)) =
  (p2 ^-^ p1, PenCurve (CubicBezier p1 p2 p3 p4)) :
  (p4 ^-^ p3, PenCorner p4) : 
  pathPartitions tl

pathPartitions _ = []

splitPartitions :: [(DPoint, b)] -> FindDir b
splitPartitions [] = error "splitPartitions: empty"
splitPartitions ps@((dir1,_):_) = FindDir dir1 (V.fromList lft) (V.fromList rt)
  where
    (lft, rt_) = span (clockwise dir1 . fst) ps
    rt = last lft:rt_

noTranslate :: Transform Double -> Transform Double
noTranslate (Transform a b _ c d _) =
  Transform a b 0 c d 0

instance (Floating a, Eq a) => AffineTransform (Pen a) a where
  {-# SPECIALIZE transform :: Transform Double -> Pen Double -> Pen Double #-}
  transform t (PenEllipse trans _ _) =
    let t2@(Transform a b c d e f) = transform t trans
    in case inverse $ noTranslate t2 of
      Nothing -> pathToPen $
        ClosedPath [
        (Point c f ^+^ p, JoinLine),
        (Point c f ^-^ p, JoinLine)]
        where
          p | a /= 0 && b /= 0 =
                sqrt(1 + a*a/(b*b)) *^ Point a d
            | d /= 0 && e /= 0 =
                sqrt(1 + d*d/(e*e)) *^ Point a d
            | a /= 0 = Point (a+d) 0
            | b /= 0 = Point 0 (b+e)
              -- singular point: create tiny pen instead of an error
            | otherwise = Point 1e-5 1e-5
      Just inv ->
        PenEllipse t2 (transform rotate90L inv) (transform rotate90R inv)

  transform t (PenPath segments) =
    PenPath $ map (transformSegment t) segments

transformSegment :: Num b => Transform b -> PenSegment b -> PenSegment b
transformSegment t (p, PenCorner q) =
  ((transform t (q^+^p) ^-^ q'), q')
  where q' = transform t q

transformSegment t (PenCurve p c) =
  ((transform t (cubicC0 c^+^p) ^-^ cubicC0 c'), PenCurve  c')
  where c' = transform t c

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
             -> Bool        -- ^ Calculate the curve faster but with
                            -- more subcurves
             -> [CubicBezier Double]        -- ^ The offset curve
bezierOffset cb dist (Just m) tol faster =
  approximatePathMax m (bezierOffsetPoint cb dist) 15 tol 0 1 faster

bezierOffset cb dist Nothing tol faster =
  approximatePath (bezierOffsetPoint cb dist) 15 tol 0 1 faster

penOffset :: Pen Double -> Point Double -> Point Double
penOffset (PenEllipse trans leftInv _) dir =
  transform trans $ normVector $ leftInv $* dir

penOffset (PenPath segments) dir =
    pathOffsetPoint (cycle segments) dir

penOffsetFun :: Pen Double -> (Double -> (DPoint, DPoint)) -> Double -> (Point Double, Point Double)
penOffsetFun pen f t =
  (px ^+^ penOffset pen px', px')
  where
    (px, px') = f t

firstPoint :: PenSegment a -> Point a
firstPoint (PenCorner _ p) = p
firstPoint (PenCurve _ c) = cubicC0 c

pathOffsetPoint :: [PenSegment Double] -> Point Double -> Point Double
pathOffsetPoint (PenCorner c p:b:rest) dir
  | vectorCross dir c > 0 = pathOffsetPoint (b:rest) dir
  | vectorCross dir (firstPoint b ^-^ p) > 0 = p
  | otherwise = pathOffsetPoint (b:rest) dir
  
pathOffsetPoint (PenCurve c curve@(CubicBezier p1 p2 p3 p4):b:rest) dir
  | vectorCross dir c > 0 = pathOffsetPoint (b:rest) dir
  | vectorCross dir (p2 ^-^ p1) > 0 = p1
  | vectorCross dir (p3 ^-^ p4) > 0 =
    case findBezierTangent dir curve of
      (t:_) -> evalBezier curve t
      [] -> p4
  | vectorCross dir (firstPoint b ^-^ p4) > 0 = p4
  | otherwise = pathOffsetPoint (b:rest) dir

pathOffsetPoint _ _ = error "unexpected end of list"

segDirs :: [(DPoint, PathJoin Double)] -> Point Double -> [(DPoint, DPoint)]
segDirs [] _ = []
segDirs [(p, JoinLine)] q = [(dp, dp)]
  where dp = q ^-^ p
segDirs [(p1, JoinCurve p2 p3 )] p4 = [(p2 ^-^ p1, p4 ^-^ p3)]
segDirs ((p, JoinLine):r@((q, _):_)) s = (dp, dp): segDirs r s
  where dp = q ^-^ p
segDirs ((p1, JoinCurve p2 p3 ):r@((p4,_):_)) q = (p2 ^-^ p1, p4 ^-^ p3):segDirs r q

penStrokeOpen :: Int -> Double -> Bool -> Pen Double -> OpenPath Double -> [ClosedPath Double]
penStrokeOpen samples tol fast pen (OpenPath segments p)  =
  union [closeOpenPath path] NonZero tol
  where
    dirs = segDirs segments (fst $ head segments)
    fdirs = map fst (tail dirs)
    fd = fst $ head dirs
    ld = snd $ last dirs
    ldirs = map snd dirs 
    pts = map fst (tail segments) ++ [p]
    leftJoins = zipWith (penJoinLeft pen) ldirs fdirs
    leftStrokes = zipWith (strokeLeft samples tol fast pen) segments pts
    rightJoins = zipWith (penJoinRight pen) ldirs fdirs
    rightStrokes = zipWith (strokeRight samples tol fast pen) segments pts
    path =
      mconcat $
      penJoinLeft pen (turnAround fd) fd :
      interleave leftStrokes leftJoins ++
      penJoinLeft pen ld (turnAround ld) :
      reverse (interleave rightStrokes rightJoins)

interleave :: [a] -> [a] -> [a]
interleave [] xs = xs
interleave xs [] = xs
interleave (x:xs) (y:ys) = x:y:interleave xs ys 

--penStrokeClosed :: ClosedPath Double -> Pen Double -> Double -> [ClosedPath Double]
penStrokeClosed :: Int -> Double -> Bool -> Pen Double -> ClosedPath Double
                -> [ClosedPath Double]
penStrokeClosed _ _ _ _ (ClosedPath [])  = [ClosedPath []]
penStrokeClosed samples tol fast pen (ClosedPath segments) =
  union [closeOpenPath leftPath, closeOpenPath rightPath] NonZero tol
  where
    dirs = segDirs segments (fst $ head segments)
    fdirs = map fst (tail dirs) ++ [fst (head dirs)]
    ldirs = map snd dirs
    pts = map fst (tail segments) ++ [fst (head segments)]
    leftJoins = zipWith (penJoinLeft pen) ldirs fdirs
    leftStrokes = zipWith (strokeLeft samples tol fast pen) segments pts
    rightJoins = zipWith (penJoinRight pen) ldirs fdirs
    rightStrokes = zipWith (strokeRight samples tol fast pen) segments pts
    leftPath =
      mconcat $ interleave leftStrokes leftJoins
    rightPath =
      mconcat $ reverse $ interleave rightStrokes rightJoins

strokeLeft :: Int -> Double -> Bool -> Pen Double -> (DPoint, PathJoin Double) -> DPoint -> OpenPath Double
strokeLeft _ _ _ pen (p, JoinLine) q =
  OpenPath [(p ^+^ offset, JoinLine)] (q ^+^ offset)
  where offset = penOffset pen (q ^-^ p)

strokeLeft samples tol fast pen (p1, JoinCurve p2 p3) p4 =
  curvesToOpen $ approximatePath
  (penOffsetFun pen (evalBezierDeriv (CubicBezier p1 p2 p3 p4)))
  samples tol 0 1 fast

strokeRight :: Int -> Double -> Bool -> Pen Double -> (DPoint, PathJoin Double) -> DPoint -> OpenPath Double
strokeRight _ _ _ pen (p, JoinLine) q =
  OpenPath [(q ^+^ offset, JoinLine)] (p ^+^ offset)
  where offset = penOffset pen (p ^-^ q)

strokeRight samples tol fast pen (p1, JoinCurve p2 p3) p4 =
  curvesToOpen $ approximatePath
  (penOffsetFun pen (evalBezierDeriv (CubicBezier p4 p3 p2 p1)))
  samples tol 0 1 fast

penJoinLeft :: Pen Double -> DPoint -> DPoint -> OpenPath Double
penJoinLeft = penJoin

penJoinRight :: Pen Double -> DPoint -> DPoint -> OpenPath Double
penJoinRight pen from to = penJoin pen (turnAround to) (turnAround from)

ellipticArc :: Transform Double -> Transform Double
            -> Point Double -> Point Double -> CubicBezier Double
ellipticArc trans leftInv from to =
  trans $* bezierArc
  (vectorAngle $ leftInv $* from)
  (vectorAngle $ leftInv $* to)

segmentsToPath :: (Eq a) => [PenSegment a] -> OpenPath a
segmentsToPath [PenCorner _ q] =
  OpenPath [] q
segmentsToPath [PenCurve _ (CubicBezier p1 p2 p3 p4)] =
  OpenPath [(p1, JoinCurve p2 p3)] p4
  
segmentsToPath (PenCorner _ p:r) =
  consOpenPath p JoinLine (segmentsToPath r)

segmentsToPath (PenCurve _ (CubicBezier p1 p2 p3 p4):r) =
  consOpenPath p1 (JoinCurve p2 p3) $
  case r of
    (PenCurve _ (CubicBezier q1 _ _ _):_)
      | p4 /= q1  -> consOpenPath p4 JoinLine $ segmentsToPath r
    _ -> segmentsToPath r

segmentsToPath [] = emptyOpenPath  

emptyOpenPath :: OpenPath a
emptyOpenPath = OpenPath [] (error "empty path")
  
penJoin :: Pen Double -> Point Double
        -> Point Double -> OpenPath Double
penJoin pen@(PenEllipse trans leftInv _) from to
  | dir == 0 = emptyOpenPath
  | dir > 0 &&
    sameQuadrant from to =
    curvesToOpen [ellipticArc trans leftInv from to]
  | otherwise =
      curvesToOpen [ellipticArc trans leftInv from next] <>
      penJoin pen next to
      where next = nextVector from
            dir = vectorCross from to

penJoin (PenPath segments) from to =
  segmentsToPath $
  nextSegments (firstSegment (cycle segments) from) to

firstSegment :: [PenSegment Double] -> Point Double -> [PenSegment Double]
firstSegment segments@(PenCorner c _:q:rest) from
  | vectorCross from c > 0 =
    firstSegment (q:rest) from
  | otherwise = segments

firstSegment segments@(PenCurve c curve@(CubicBezier p1 p2 p3 p4):q:rest) from
  | vectorCross from c > 0 = firstSegment (q:rest) from
  | vectorCross from (p2 ^-^ p1) > 0 = segments
  | vectorCross from (p4 ^-^ p3) > 0 =
      case findBezierTangent from curve of
        (t:_) -> PenCurve from (snd (splitBezier curve t)):q:rest
        _ -> q:rest
  | vectorCross from (firstPoint q ^-^ p4) > 0 =
      PenCorner (firstPoint q ^-^ p4) p4:q:rest
  | otherwise = firstSegment (q:rest) from

firstSegment _ _ = error "firstsegment: finite list"  

nextSegments :: [PenSegment Double] -> Point Double -> [PenSegment Double]
nextSegments (PenCorner c p:q:rest) to
  | vectorCross to c > 0 =
      PenCorner c p: nextSegments (q:rest) to
  | otherwise = []

nextSegments (pc@(PenCurve c curve@(CubicBezier p1 p2 p3 p4)):q:rest) to  
  | vectorCross to c > 0 = pc: nextSegments (q:rest) to
  | vectorCross to (p2 ^-^ p1) > 0 = []
  | vectorCross to (p4 ^-^ p3) > 0 =
      case findBezierTangent to curve of
        (t:_) -> [PenCurve c (fst (splitBezier curve t))]
        _ -> []
  | vectorCross to (firstPoint q ^-^ p4) > 0 =
      [PenCorner (firstPoint q ^-^ p4) p4]
  | otherwise = pc:firstSegment (q:rest) to

nextSegments _ _ = error "nextSegments: finite list"

sameQuadrant :: (Num a, Eq a) => Point a -> Point a -> Bool
sameQuadrant v w =
  signum (pointX v) /= -signum (pointX w) &&
  signum (pointY v) /= -signum (pointY w)

nextVector :: (Ord a1, Num a1, Num a) => Point a1 -> Point a
nextVector v
  | pointX v >= 0 &&
    pointY v > 0 = Point 1 0
  | pointX v > 0 &&
    pointY v <= 0 = Point 0 (-1)
  | pointX v <= 0 &&
    pointY v < 0 = Point (-1) 0
  | otherwise = Point 0 1
