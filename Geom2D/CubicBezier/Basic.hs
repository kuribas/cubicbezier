{-# LANGUAGE BangPatterns, FlexibleInstances, MultiParamTypeClasses, DeriveFunctor, ViewPatterns #-}
module Geom2D.CubicBezier.Basic
       (CubicBezier (..), QuadBezier (..), AnyBezier (..), GenericBezier(..),
        PathJoin (..), ClosedPath(..), OpenPath (..), AffineTransform (..), anyToCubic, anyToQuad,
        openPathCurves, closedPathCurves, curvesToOpen, curvesToClosed,
        bezierParam, bezierParamTolerance, reorient, bezierToBernstein,
        evalBezierDerivs, evalBezier, evalBezierDeriv, findBezierTangent, quadToCubic,
        bezierHoriz, bezierVert, findBezierInflection, findBezierCusp,
        bezierArc, arcLength, arcLengthParam, splitBezier, bezierSubsegment, splitBezierN,
        colinear)
       where
import Geom2D
import Geom2D.CubicBezier.Numeric
import Math.BernsteinPoly
import Numeric.Integration.TanhSinh
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector.Unboxed.Mutable as MV

-- | A cubic bezier curve.
data CubicBezier a = CubicBezier {
  cubicC0 :: !(Point a),
  cubicC1 :: !(Point a),
  cubicC2 :: !(Point a),
  cubicC3 :: !(Point a)}
  deriving (Eq, Show, Functor)

-- | A quadratic bezier curve.
data QuadBezier a = QuadBezier {
  quadC0 :: !(Point a),
  quadC1 :: !(Point a),
  quadC2 :: !(Point a)}
  deriving (Eq, Show, Functor)

-- Use a tuple, because it has 0(1) unzip when using unboxed vectors.
-- | A bezier curve of any degree.
data AnyBezier a = AnyBezier (V.Vector (a, a))
                   
class GenericBezier b where
  degree :: (V.Unbox a) => b a -> Int
  toVector :: (V.Unbox a) => b a -> V.Vector (a, a)
  unsafeFromVector :: (V.Unbox a) => V.Vector (a, a) -> b a

instance GenericBezier CubicBezier where
  degree _ = 3
  toVector (CubicBezier (Point ax ay) (Point bx by)
            (Point cx cy) (Point dx dy)) =
    V.create $ do
      v <- MV.new 4
      MV.write v 0 (ax, ay)
      MV.write v 1 (bx, by)
      MV.write v 2 (cx, cy)
      MV.write v 3 (dx, dy)
      return v
  unsafeFromVector v = CubicBezier
                       (uncurry Point $ v `V.unsafeIndex` 0)
                       (uncurry Point $ v `V.unsafeIndex` 1)
                       (uncurry Point $ v `V.unsafeIndex` 2)
                       (uncurry Point $ v `V.unsafeIndex` 3)

instance GenericBezier QuadBezier where
  degree _ = 2
  toVector (QuadBezier (Point ax ay) (Point bx by)
            (Point cx cy)) =
    V.create $ do
      v <- MV.new 3
      MV.write v 0 (ax, ay)
      MV.write v 1 (bx, by)
      MV.write v 2 (cx, cy)
      return v
  unsafeFromVector v = QuadBezier
                       (uncurry Point $ v `V.unsafeIndex` 0)
                       (uncurry Point $ v `V.unsafeIndex` 1)
                       (uncurry Point $ v `V.unsafeIndex` 2)

instance GenericBezier AnyBezier where
  degree (AnyBezier b) = V.length b
  toVector (AnyBezier v) = v
  unsafeFromVector = AnyBezier

data PathJoin a = JoinLine |
                  JoinCurve (Point a) (Point a)
              deriving (Show, Functor)
data OpenPath a = OpenPath [(Point a, PathJoin a)] (Point a) 
                  deriving (Show, Functor)
data ClosedPath a = ClosedPath [(Point a, PathJoin a)]
                  deriving (Show, Functor)

instance (Num a) => AffineTransform (CubicBezier a) a where
  {-# SPECIALIZE transform :: Transform Double -> CubicBezier Double -> CubicBezier Double #-}
  transform t (CubicBezier c0 c1 c2 c3) =
    CubicBezier (transform t c0) (transform t c1) (transform t c2) (transform t c3)

-- | Return the open path as a list of curves.
openPathCurves :: Fractional a => OpenPath a -> [CubicBezier a]
openPathCurves (OpenPath curves p) = go curves p
  where
    go [] _ = []
    go [(p0, jn)] p = [makeCB p0 jn p]
    go ((p0, jn):rest@((p1,_):_)) p =
      makeCB p0 jn p1 : go rest p
    makeCB p0 (JoinLine) p1 =
      CubicBezier p0 (interpolateVector p0 p1 (1/3))
      (interpolateVector p0 p1 (2/3)) p1
    makeCB p0 (JoinCurve p1 p2) p3 =
      CubicBezier p0 p1 p2 p3

-- | Return the closed path as a list of curves
closedPathCurves :: Fractional a => ClosedPath a -> [CubicBezier a]
closedPathCurves (ClosedPath []) = []
closedPathCurves (ClosedPath (cs@((p1, _):_))) =
  openPathCurves (OpenPath cs p1)

-- | Make an open path from a list of curves.  The last control point
-- of each curve except the last is ignored.
curvesToOpen :: [CubicBezier a] -> OpenPath a
curvesToOpen [] = OpenPath [] undefined
curvesToOpen [CubicBezier p0 p1 p2 p3] =
  OpenPath [(p0, JoinCurve p1 p2)] p3
curvesToOpen (CubicBezier p0 p1 p2 _:cs) =
  OpenPath ((p0, JoinCurve p1 p2):rest) lastP
  where
    OpenPath rest lastP = curvesToOpen cs

-- | Make an open path from a list of curves.  The last control point
-- of each curve is ignored.
curvesToClosed :: [CubicBezier a] -> ClosedPath a
curvesToClosed cs = ClosedPath cs2
  where
    OpenPath cs2 _ = curvesToOpen cs


-- | safely convert from `AnyBezier' to `CubicBezier`
anyToCubic :: (V.Unbox a) => AnyBezier a -> Maybe (CubicBezier a)
anyToCubic b@(AnyBezier v)
  | degree b == 3 = Just $ unsafeFromVector v
  | otherwise = Nothing

-- | safely convert from `AnyBezier' to `QuadBezier`
anyToQuad :: (V.Unbox a) => AnyBezier a -> Maybe (QuadBezier a)
anyToQuad b@(AnyBezier v)
  | degree b == 2 = Just $ unsafeFromVector v
  | otherwise = Nothing

evalBezierDerivsCubic :: Fractional a =>
                         CubicBezier a -> a -> [Point a]
evalBezierDerivsCubic (CubicBezier a b c d) t =
  [p, p', p'', p''', Point 0 0]
  where
    u = 1-t
    t2 = t*t
    t3 = t2*t
    da = 3*^(b^-^a)
    db = 3*^(c^-^b)
    dc = 3*^(d^-^c)
    p = u*^(u*^(u*^a ^+^ 3*t*^b) ^+^ 3*t2*^c) ^+^ t3*^d
    p' = u*^(u*^da ^+^ 2*t*^db) ^+^ t2*^dc
    p'' = 2*u*^(db^-^da) ^+^ 2*t*^(dc^-^db)
    p''' = 2*^(dc^-^2*^db^+^da)
{-# SPECIALIZE evalBezierDerivsCubic :: CubicBezier Double -> Double -> [DPoint] #-}    

evalBezierDerivsQuad :: Fractional a =>
                        QuadBezier a -> a -> [Point a]
evalBezierDerivsQuad (QuadBezier a b c) t = [p, p', p'', Point 0 0]
  where
    u = 1-t
    t2 = t*t
    p = u*^(u*^a ^+^ 2*t*^b) ^+^ t2*^c
    p' = 2*^(u*^(b^-^a) ^+^ t*^(c^-^b))
    p'' = 2*^(c^-^ 2*^b ^+^ a)
{-# SPECIALIZE evalBezierDerivsQuad :: QuadBezier Double -> Double -> [DPoint] #-}        

-- | Evaluate the bezier and all its derivatives using the modified horner algorithm.
evalBezierDerivs :: (GenericBezier b, V.Unbox a, Fractional a) =>
                    b a -> a -> [Point a]
evalBezierDerivs b t =
  zipWith Point (bernsteinEvalDerivs (BernsteinPoly x) t)
  (bernsteinEvalDerivs (BernsteinPoly y) t)
  where (x, y) = V.unzip $ toVector b
{-# SPECIALIZE evalBezierDerivs :: AnyBezier Double -> Double -> [DPoint] #-}
{-# NOINLINE [2] evalBezierDerivs #-}
{-# RULES "evalBezierDerivs/cubic" evalBezierDerivs = evalBezierDerivsCubic #-}
{-# RULES "evalBezierDerivs/quad"  evalBezierDerivs = evalBezierDerivsQuad #-}

-- | Return True if the param lies on the curve, iff it's in the interval @[0, 1]@.
bezierParam :: (Ord a, Num a) => a -> Bool
bezierParam t = t >= 0 && t <= 1

-- | Convert a tolerance from the codomain to the domain of the bezier
-- curve, by dividing by the maximum velocity on the curve.  The
-- estimate is conservative, but holds for any value on the curve.
bezierParamTolerance :: (GenericBezier b) => b Double -> Double -> Double
bezierParamTolerance (toVector -> v) eps = eps / maxVel
  where 
    maxVel = 3 * V.maximum (V.zipWith vectorDistance (V.map (uncurry Point) v)
                            (V.map (uncurry Point) $ V.tail v))

-- | Reorient to the curve B(1-t).
reorient :: (GenericBezier b, V.Unbox a) => b a -> b a
reorient = unsafeFromVector . V.reverse . toVector
{-# SPECIALIZE reorient :: (V.Unbox a) => AnyBezier a -> AnyBezier a #-}
{-# NOINLINE [2] reorient #-}

reorientCubic :: CubicBezier a -> CubicBezier a
reorientCubic (CubicBezier a b c d) = CubicBezier d c b a

reorientQuad :: QuadBezier a -> QuadBezier a
reorientQuad (QuadBezier a b c) = QuadBezier c b a
{-# RULES "reorient/cubic" reorient = reorientCubic #-}
{-# RULES "reorient/quad"  reorient = reorientQuad #-}

-- | Give the bernstein polynomial for each coordinate.
bezierToBernstein :: (GenericBezier b, MV.Unbox a) =>
                     b a -> (BernsteinPoly a, BernsteinPoly a)
bezierToBernstein b = (BernsteinPoly x, BernsteinPoly y)
  where (x, y) = V.unzip $ toVector b

-- | Calculate a value on the bezier curve.
evalBezier :: (GenericBezier b, MV.Unbox a, Fractional a) =>
              b a -> a -> Point a
evalBezier bc t = head $ evalBezierDerivs bc t
{-# SPECIALIZE evalBezier :: AnyBezier Double -> Double -> DPoint #-}
{-# NOINLINE [2] evalBezier #-}

evalBezierCubic :: Fractional a =>
                   CubicBezier a -> a -> Point a
evalBezierCubic (CubicBezier a b c d) t =
  u*^(u*^(u*^a ^+^ 3*t*^b) ^+^ 3*t2*^c) ^+^ t3*^d
  where
    u = 1-t
    t2 = t*t
    t3 = t2*t
{-# SPECIALIZE evalBezierCubic :: CubicBezier Double -> Double -> DPoint #-}
    
evalBezierQuad :: Fractional a =>
                  QuadBezier a -> a -> Point a
evalBezierQuad (QuadBezier a b c) t = 
  u*^(u*^a ^+^ 2*t*^b) ^+^ t2*^c
  where
    u = 1-t
    t2 = t*t
{-# SPECIALIZE evalBezierQuad :: QuadBezier Double -> Double -> DPoint #-}

{-# RULES "evalBezier/cubic" evalBezier = evalBezierCubic #-}
{-# RULES "evalBezier/quad"  evalBezier = evalBezierQuad #-}

-- | Calculate a value and the first derivative on the curve.
evalBezierDeriv :: (V.Unbox a, Fractional a) =>
                   GenericBezier b => b a -> a -> (Point a, Point a)
evalBezierDeriv bc t = (b,b')
  where
    (b:b':_) = evalBezierDerivs bc t

-- | @findBezierTangent p b@ finds the parameters where
-- the tangent of the bezier curve @b@ has the same direction as vector p.

-- Use the formula tx * B'y(t) - ty * B'x(t) = 0 where
-- B'x is the x value of the derivative of the Bezier curve.
findBezierTangent :: DPoint -> CubicBezier Double -> [Double]
findBezierTangent (Point tx ty) (CubicBezier (Point x0 y0) (Point x1 y1) (Point x2 y2) (Point x3 y3)) = 
  filter bezierParam $ quadraticRoot a b c
    where
      a = tx*((y3 - y0) + 3*(y1 - y2)) - ty*((x3 - x0) + 3*(x1 - x2))
      b = 2*(tx*((y2 + y0) - 2*y1) - ty*((x2 + x0) - 2*x1))
      c = tx*(y1 - y0) - ty*(x1 - x0)

-- | Find the parameter where the bezier curve is horizontal.
bezierHoriz :: CubicBezier Double -> [Double]
bezierHoriz = findBezierTangent (Point 1 0)

-- | Find the parameter where the bezier curve is vertical.
bezierVert :: CubicBezier Double -> [Double]
bezierVert = findBezierTangent (Point 0 1)

-- | Find inflection points on the curve.
-- Use the formula B_x''(t) * B_y'(t) - B_y''(t) * B_x'(t) = 0 with
-- B_x'(t) the x value of the first derivative at t, B_y''(t) the y
-- value of the second derivative at t
findBezierInflection :: CubicBezier Double -> [Double]
findBezierInflection (CubicBezier (Point x0 y0) (Point x1 y1) (Point x2 y2) (Point x3 y3)) =
  filter bezierParam $ quadraticRoot a b c
    where
      ax = x1 - x0
      bx = x3 - x1 - ax
      cx = x3 - x2 - ax - 2*bx
      ay = y1 - y0
      by = y2 - y1 - ay
      cy = y3 - y2 - ay - 2*by
      a = bx*cy - by*cx
      b = ax*cy - ay*cx
      c = ax*by - ay*bx

-- | Find the cusps of a bezier.

-- find a cusp.  We look for points where the tangent is both horizontal
-- and vertical, which is only true for the zero vector.
findBezierCusp :: CubicBezier Double -> [Double]
findBezierCusp b = filter vertical $ bezierHoriz b
  where vertical = (== 0) . pointY . snd . evalBezierDeriv b

-- | @bezierArc startAngle endAngle@ approximates an arc on the unit circle with
-- a single cubic béziér curve.  Maximum deviation is <0.03% for arcs
-- 90° degrees or less.
bezierArc :: Double -> Double -> CubicBezier Double
bezierArc start end = CubicBezier p0 p1 p2 p3
  where
    p0 = dirVector start
    p3 = dirVector end
    p1 = p0 ^+^ k *^ (rotate90L $* p0)
    p2 = p1 ^+^ k *^ (rotate90R $* p3)
    k = 4/3 * tan((end-start)/4)

-- | @arcLength c t tol finds the arclength of the bezier c at t, within given tolerance tol.

arcLength :: CubicBezier Double -> Double -> Double -> Double
arcLength b@(CubicBezier c0 c1 c2 c3) t eps =
  if eps / maximum [vectorDistance c0 c1,
                    vectorDistance c1 c2,
                    vectorDistance c2 c3] > 1e-10
  then (signum t *) $ fst $
       arcLengthEstimate (fst $ splitBezier b t) eps
  else arcLengthQuad b t eps

arcLengthQuad :: CubicBezier Double -> Double -> Double -> Double
arcLengthQuad b t eps = result $ absolute eps $
                        trap distDeriv 0 t
  where distDeriv t' = vectorMag $ snd $ evalD t'
        evalD = evalBezierDeriv b

outline :: CubicBezier Double -> Double
outline (CubicBezier c0 c1 c2 c3) =
  vectorDistance c0 c1 +
  vectorDistance c1 c2 +
  vectorDistance c2 c3

arcLengthEstimate :: CubicBezier Double -> Double -> (Double, (Double, Double))
arcLengthEstimate b eps = (arclen, (estimate, ol))
  where
    estimate = (4*(olL+olR) - ol) / 3
    (bl, br) = splitBezier b 0.5
    ol = outline b
    (arcL, (estL, olL)) = arcLengthEstimate bl eps
    (arcR, (estR, olR)) = arcLengthEstimate br eps
    arclen | abs(estL + estR - estimate) < eps = estL + estR
           | otherwise = arcL + arcR

-- | arcLengthParam c len tol finds the parameter where the curve c has the arclength len,
-- within tolerance tol.
arcLengthParam :: CubicBezier Double -> Double -> Double -> Double
arcLengthParam b len eps =
  arcLengthP b len ol (len/ol) 1 eps
  where ol = outline b

-- Use the Newton rootfinding method.  Start with large tolerance
-- values, and decrease tolerance as we go closer to the root.
arcLengthP :: CubicBezier Double -> Double -> Double ->
              Double -> Double -> Double -> Double
arcLengthP !b !len !tot !t !dt !eps
  | abs diff < eps = t - newDt
  | otherwise = arcLengthP b len tot (t - newDt) newDt eps
  where diff = arcLength b t (max (abs (dt*tot/50)) (eps/2)) - len
        newDt = diff / vectorMag (snd $ evalBezierDeriv b t)

-- | Convert a quadratic bezier to a cubic bezier.
quadToCubic :: (Fractional a) =>
               QuadBezier a -> CubicBezier a
quadToCubic (QuadBezier a b c) =
  CubicBezier a (1/3*^(a ^+^ 2*^b)) (1/3*^(2*^b ^+^ c)) c

-- | Split a bezier curve into two curves.
splitBezier :: (V.Unbox a, Fractional a) =>
               GenericBezier b => b a -> a -> (b a, b a)
splitBezier b t =
  (unsafeFromVector $ V.zip (bernsteinCoeffs x1) (bernsteinCoeffs y1),
   unsafeFromVector $ V.zip (bernsteinCoeffs x2) (bernsteinCoeffs y2))
  where
    (x, y) = bezierToBernstein b
    (x1, x2) = bernsteinSplit x t
    (y1, y2) = bernsteinSplit y t
{-# NOINLINE [2] splitBezier #-}
{-# SPECIALIZE splitBezier :: AnyBezier Double -> Double -> (AnyBezier Double, AnyBezier Double) #-}

-- | Split a bezier curve into two curves.
splitBezierCubic :: (Fractional a) =>  CubicBezier a -> a -> (CubicBezier a, CubicBezier a)
splitBezierCubic (CubicBezier a b c d) t =
  let ab = interpolateVector a b t
      bc = interpolateVector b c t
      cd = interpolateVector c d t
      abbc = interpolateVector ab bc t
      bccd = interpolateVector bc cd t
      mid = interpolateVector abbc bccd t
  in (CubicBezier a ab abbc mid, CubicBezier mid bccd cd d)
{-# SPECIALIZE splitBezierCubic :: CubicBezier Double -> Double -> (CubicBezier Double, CubicBezier Double) #-}     

-- | Split a bezier curve into two curves.
splitBezierQuad :: (Fractional a) =>  QuadBezier a -> a -> (QuadBezier a, QuadBezier a)
splitBezierQuad (QuadBezier a b c) t =
  let ab = interpolateVector a b t
      bc = interpolateVector b c t
      mid = interpolateVector ab bc t
  in (QuadBezier a ab mid, QuadBezier mid bc c)
{-# SPECIALIZE splitBezierQuad :: QuadBezier Double -> Double -> (QuadBezier Double, QuadBezier Double) #-}
{-# RULES "splitBezier/cubic" splitBezier = splitBezierCubic #-}
{-# RULES "splitBezier/quad"  splitBezier = splitBezierQuad #-}


-- | Return the subsegment between the two parameters.
bezierSubsegment :: (Ord a, V.Unbox a, Fractional a) => GenericBezier b =>
                    b a -> a -> a -> b a
bezierSubsegment b t1 t2 
  | t1 > t2   = bezierSubsegment b t2 t1
  | t2 == 0   = fst $ splitBezier b t1
  | otherwise = snd $ flip splitBezier (t1/t2) $
                fst $ splitBezier b t2
{-# SPECIALIZE bezierSubsegment :: CubicBezier Double -> Double -> Double -> CubicBezier Double #-}
{-# SPECIALIZE bezierSubsegment :: QuadBezier Double -> Double -> Double -> QuadBezier Double #-}

-- | Split a bezier curve into a list of beziers
-- The parameters should be in ascending order or
-- the result is unpredictable.
splitBezierN :: (Ord a, V.Unbox a, Fractional a) =>
                GenericBezier b => b a -> [a] -> [b a]
splitBezierN c [] = [c]
splitBezierN c [t] = [a, b] where
  (a, b) = splitBezier c t
splitBezierN c (t:u:rest) =
  bezierSubsegment c 0 t :
  bezierSubsegment c t u :
  tail (splitBezierN c $ u:rest)
{-# SPECIALIZE splitBezierN :: CubicBezier Double -> [Double] -> [CubicBezier Double] #-}
{-# SPECIALIZE splitBezierN :: QuadBezier Double -> [Double] -> [QuadBezier Double] #-}

-- | Return False if some points fall outside a line with a thickness of the given tolerance.

-- fat line calculation taken from the bezier-clipping algorithm (Sederberg)
colinear :: CubicBezier Double -> Double -> Bool
colinear (CubicBezier !a !b !c !d) eps = dmax - dmin < eps
  where ld = lineDistance (Line a d)
        d1 = ld b
        d2 = ld c
        (dmin, dmax) | d1*d2 > 0 = (3/4 * minimum [0, d1, d2],
                                    3/4 * maximum [0, d1, d2])
                     | otherwise = (4/9 * minimum [0, d1, d2],
                                    4/9 * maximum [0, d1, d2])
