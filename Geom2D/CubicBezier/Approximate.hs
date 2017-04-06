{-# LANGUAGE BangPatterns, MultiWayIf #-}
module Geom2D.CubicBezier.Approximate
       (approximatePath, approximateQuadPath, approximatePathMax, approximateQuadPathMax,
        approximateCubic)
       where
import Geom2D
import Geom2D.CubicBezier.Basic
import Geom2D.CubicBezier.Numeric
import Data.Maybe
import Data.List
import qualified Data.Vector.Unboxed as V
import qualified Data.Map as M
import Data.Function

interpolate :: (Num a) => a -> a -> a -> a
interpolate a b x = (1-x)*a + x*b

-- | Approximate a function with piecewise cubic bezier splines using
-- a least-squares fit, within the given tolerance.  Each subcurve is
-- approximated by using a finite number of samples.  It is recommended
-- to avoid changes in direction by subdividing the original function
-- at points of inflection.
approximatePath :: (V.Unbox a, Ord a, Floating a) =>
                   (a -> (Point a, Point a)) -- ^ The function to approximate and it's derivative
                -> Int
                   -- ^ The number of discrete samples taken to
                   -- approximate each subcurve.  More samples are
                   -- more precise but take more time to calculate.
                   -- For good precision 16 is a good candidate.
                -> a                         -- ^ The tolerance
                -> a                         -- ^ The lower parameter of the function      
                -> a                         -- ^ The upper parameter of the function
                -> Bool
                -- ^ Calculate the result faster, but with more
                -- subcurves.  Runs typically 10 times faster, but
                -- generates 50% more subcurves.   Useful for interactive use.
                -> [CubicBezier a]
approximatePath f n tol tmin tmax fast
  | err < tol = [curve]
  | otherwise = approximatePath' f n tol tmin tmax fast
  where
    (curve, err) = approx1cubic n f tmin tmax (if fast then 0 else 5)
{-# SPECIALIZE approximatePath :: (Double -> (DPoint, DPoint)) -> Int -> Double
                               -> Double -> Double -> Bool -> [CubicBezier Double]  #-}

-- | Approximate a function with piecewise quadratic bezier splines
-- using a least-squares fit, within the given tolerance.  It is
-- recommended to avoid changes in direction by subdividing the
-- original function at points of inflection.
approximateQuadPath :: (Show a, V.Unbox a, Ord a, Floating a) =>
                   (a -> (Point a, Point a)) -- ^ The function to approximate and it's derivative
                -> a                         -- ^ The tolerance
                -> a                         -- ^ The lower parameter of the function      
                -> a                         -- ^ The upper parameter of the function
                -> Bool                      
                -- ^ Calculate the result faster, but with more
                -- subcurves.
                -> [QuadBezier a]
approximateQuadPath f tol tmin tmax fast
  | err < tol = [curve]
  | otherwise = approximateQuad' f tol tmin tmax fast
  where
    curve = approx1quad f tmin tmax
    err = maxDist f curve tmin tmax (if fast then 0 else 5)
{-# SPECIALIZE approximateQuadPath :: (Double -> (DPoint, DPoint)) -> Double ->
    Double -> Double -> Bool -> [QuadBezier Double] #-}
      
-- find the distance between the function at t and the quadratic bezier.
-- calculate the value and derivative at t, and improve the closeness of t.
quadDist :: (V.Unbox a, Ord a, Floating a) =>
            (a -> (Point a, Point a)) -> QuadBezier a -> a
         -> a -> Int -> a -> a
quadDist f qb tmin tmax maxiter t =
  quadDist' f qb tmin tmax t (fst (f $ interpolate tmin tmax t)) maxiter (evalBezierDeriv qb t)

quadDist' :: (V.Unbox a, Ord a, Floating a) =>
             (a -> (Point a, Point a))
          -> QuadBezier a -> a -> a -> a -> Point a -> Int -> (Point a, Point a) -> a
quadDist' f qb tmin tmax t p maxiter (b, b')
  | maxiter <= 1 || abs (err2-err1) <= err1 * (1/8) = err2
  | otherwise = quadDist' f qb tmin tmax (t+ndist) p (maxiter-1) (b2, b2')
  where dp = p ^-^ b
        err1 = vectorMag dp
        -- distance from p to the normal at b(t) / velocity at t
        ndist = (dp ^.^ b') / (b' ^.^ b')
        (b2, b2') = evalBezierDeriv qb (t+ndist)
        err2 = vectorDistance p b2
        
-- find maximum distance using golden section search
maxDist :: (V.Unbox a, Ord a, Floating a) =>
           (a -> (Point a, Point a)) ->
           QuadBezier a -> a -> a -> Int -> a
maxDist f qb tmin tmax maxiter =
  quadDist f qb tmin tmax maxiter $
  goldSearch (quadDist f qb tmin tmax maxiter) 4

approxquad :: (Ord a, Floating a) =>
              Point a -> Point a -> Point a -> Point a -> QuadBezier a
approxquad p0 p0' p1' p1
  | abs (pointY q') < abs (pointX q'*1e-3) = 
    QuadBezier p0 (interpolateVector p0 p1 0.5) p1
  | otherwise = QuadBezier p0 (p1^+^p1'^*t) p1
  where
    q = rotateVec (flipVector p0') $* p1^-^p0
    q' = rotateVec (flipVector p0') $* p1'
    t = - pointY q / pointY q'

approx1quad :: (Ord a, Floating a) =>
               (a -> (Point a, Point a)) -> a -> a -> QuadBezier a
approx1quad f tmin tmax =
  approxquad p0 p0' p1' p1
  where (p0, p0') = f tmin
        (p1, p1') = f tmax

splitQuad :: (Show a, V.Unbox a, Ord a, Floating a) =>
                a -> a -> (a -> (Point a, Point a))
                -> a -> a -> Int -> (a, a, QuadBezier a, a, QuadBezier a)
splitQuad node offset f tmin tmax maxiter =
  splitQuad' node offset f tmin tmax maxiter maxiter
  
splitQuad' :: (Show a, V.Unbox a, Ord a, Floating a) =>
                a -> a -> (a -> (Point a, Point a))
                -> a -> a -> Int -> Int -> (a, a, QuadBezier a, a, QuadBezier a)
splitQuad' node offset f tmin tmax maxiter maxiter2
  | maxiter < 1 || (err0 < 2*err1 && err0 > err1/2) =
      (tmid, err0, curve0, err1, curve1)
  | otherwise =
    splitQuad' (if err0 < err1 then node+offset else node-offset)
    (offset/2) f tmin tmax (maxiter-1) maxiter2
  where
    tmid = interpolate tmin tmax node
    curve0 = approx1quad f tmin tmid 
    err0 = maxDist f curve0 tmin tmid maxiter2
    curve1 = approx1quad f tmid tmax 
    err1 = maxDist f curve1 tmid tmax maxiter2

approximateQuad' :: (Show a, V.Unbox a, Ord a, Floating a) =>
                    (a -> (Point a, Point a)) -> 
                    a -> a -> a -> Bool ->
                    [QuadBezier a]
approximateQuad' f tol tmin tmax fast =
  (if err0 <= tol
   then [curve0]
   else approximateQuad' f tol tmin tmid fast) ++
  (if err1 <= tol
   then [curve1]
   else approximateQuad' f tol tmid tmax fast)
  where
    (tmid, err0, curve0, err1, curve1) =
     splitQuad 0.5 0.25 f tmin tmax (if fast then 0 else 5)

approximatePath' :: (V.Unbox a, Ord a, Floating a) =>
                    (a -> (Point a, Point a)) -> Int ->
                    a -> a -> a -> Bool ->
                    [CubicBezier a]
approximatePath' f n tol tmin tmax fast =
  (if err0 <= tol
   then [curve0]
   else approximatePath' f n tol tmin tmid fast) ++
  (if err1 <= tol
   then [curve1]
   else approximatePath' f n tol tmid tmax fast)
  where
    (tmid, err0, curve0, err1, curve1) =
      splitCubic n 0.5 0.25 f tmin tmax (if fast then 0 else 5)
--{-# SPECIALIZE approximatePath' :: (Double -> (Point Double, Point Double)) -> Int -> Double -> Double -> Double -> [CubicBezier Double]  #-}      

-- | Like approximatePath, but limit the number of subcurves.
approximatePathMax :: (V.Unbox a, Floating a, Ord a) =>
                      Int                        -- ^ The maximum number of subcurves
                   -> (a -> (Point a, Point a))    -- ^ The function to approximate and it's derivative
                   -> Int
                   -- ^ The number of discrete samples taken to
                   -- approximate each subcurve.  More samples are
                   -- more precise but take more time to calculate.
                   -- For good precision 16 is a good candidate.
                   -> a                          -- ^ The tolerance
                   -> a                          -- ^ The lower parameter of the function      
                   -> a                          -- ^ The upper parameter of the function
                   -> Bool
                   -- ^ Calculate the result faster, but with more
                   -- subcurves.  Runs faster (typically 10 times),
                   -- but generates more subcurves (about 50%).
                   -- Useful for interactive use.
                   -> [CubicBezier a]
approximatePathMax m f samples tol tmin tmax fast =
  approxMax f tol m fast (splitCubic samples) segments
  where segments = M.singleton err (FunctionSegment tmin tmax outline)
        (p0, p0') = f tmin
        (p1, p1') = f tmax
        ts = V.map (\i -> fromIntegral i/(fromIntegral samples+1) `asTypeOf` tmin) $
             V.enumFromN (1::Int) samples
        points = V.map (fst . f . interpolate tmin tmax) ts
        curveCb = CubicBezier p0 (p0^+^p0') (p1^-^p1') p1
        (outline, err) =
          approximateCubic curveCb points (Just ts) (if fast then 0 else 5)
{-# SPECIALIZE approximatePathMax ::
    Int -> (Double -> (Point Double, Point Double)) -> Int                      
    -> Double -> Double -> Double -> Bool -> [CubicBezier Double] #-}
data FunctionSegment b a = FunctionSegment {
  fsTmin :: !a,  -- the least t param of the segment in the original curve
  _fsTmax :: !a,  -- the max t param of the segment in the original curve
  fsCurve :: b -- the curve segment
  }

-- | Like approximateQuadPath, but limit the number of subcurves.
approximateQuadPathMax :: (V.Unbox a, Show a, Floating a, Ord a) =>
                      Int                        -- ^ The maximum number of subcurves
                   -> (a -> (Point a, Point a))    -- ^ The function to approximate and it's derivative
                   -> a                          -- ^ The tolerance
                   -> a                          -- ^ The lower parameter of the function      
                   -> a                          -- ^ The upper parameter of the function
                   -> Bool
                   -- ^ Calculate the result faster, but with more
                   -- subcurves.  Runs faster, but generates more
                   -- subcurves.  Useful for interactive use.
                   -> [QuadBezier a]
approximateQuadPathMax m f tol tmin tmax fast =
  approxMax f tol m fast splitQuad segments
  where segments = M.singleton err (FunctionSegment tmin tmax curveQd)
        (p0, p0') = f tmin
        (p1, p1') = f tmax
        curveQd = approxquad p0 p0' p1' p1
        err = maxDist f curveQd tmin tmax (if fast then 0 else 5)
{-# SPECIALIZE approximateQuadPathMax ::
    Int -> (Double -> (Point Double, Point Double))
    -> Double -> Double -> Double -> Bool -> [QuadBezier Double] #-}
-- Keep a map from maxError to FunctionSegment for each subsegment to keep
-- track of the segment with the maximum error.  This ensures a n
-- log(n) execution time, rather than n^2 when a list is used.
approxMax :: (V.Unbox a, Ord a, Floating a) =>
             (a -> (Point a, Point a)) -> a -> Int -> Bool
          -> (a -> a -> (a -> (Point a, Point a)) -> a -> a -> Int -> (a, a, b, a, b))
          -> M.Map a (FunctionSegment b a)
          -> [b]
approxMax f tol maxCurves fast splitBez segments
  | (maxCurves <= 1) || (err < tol) =
    map fsCurve $ sortBy (compare `on` fsTmin) $
    map snd $ M.toList segments
  | otherwise = approxMax f tol (maxCurves-1) fast splitBez $
                M.insert err_l (FunctionSegment t_min t_mid curve_l) $
                M.insert err_r (FunctionSegment t_mid t_max curve_r)
                newSegments
  where
    ((err, FunctionSegment t_min t_max _), newSegments) =
      M.deleteFindMax segments
    (t_mid, err_l, curve_l, err_r, curve_r) =
      splitBez 0.5 0.25 f t_min t_max (if fast then 0 else 5)
      
splitCubic :: (V.Unbox a, Ord a, Floating a) =>
                Int -> a -> a -> (a -> (Point a, Point a))
                -> a -> a -> Int -> (a, a, CubicBezier a, a, CubicBezier a)
splitCubic n node offset f tmin tmax maxiter
  | maxiter < 1 || (err0 < 2*err1 && err0 > err1/2) =
      (tmid, err0, curve0, err1, curve1)
  | otherwise = 
      splitCubic n (if err0 < err1 then node+offset else node-offset)
      (offset/2) f tmin tmax (maxiter-1)
  where
    tmid = interpolate tmin tmax node
    (curve0, err0) = approx1cubic n f tmin tmid maxiter
    (curve1, err1) = approx1cubic n f tmid tmax maxiter
    
approx1cubic :: (V.Unbox a, Ord a, Floating a) =>
           Int -> (a -> (Point a, Point a)) -> a -> a ->
           Int -> (CubicBezier a, a)
approx1cubic n f t0 t1 maxiter =
  approximateCubic curveCb points (Just ts) maxiter
  where (p0, p0') = f t0
        (p1, p1') = f t1
        ts = V.map (\i -> fromIntegral i/(fromIntegral n+1))
             (V.enumFromN 1 n :: V.Vector Int)
        points = V.map (fst . f . interpolate t0 t1) ts
        curveCb = CubicBezier p0 (p0^+^p0') (p1^+^p1') p1
{-# SPECIALIZE approx1cubic ::  Int -> (Double -> (Point Double, Point Double)) -> Double -> Double -> Int -> (CubicBezier Double, Double) #-}

-- | @approximateCubic b pts maxiter@ finds the least squares fit of a bezier
-- curve to the points @pts@.  The resulting bezier has the same first
-- and last control point as the curve @b@, and have tangents colinear with @b@.
approximateCubic :: (V.Unbox a, Ord a, Floating a) =>
                    CubicBezier a         -- ^ Curve
                    -> V.Vector (Point a) -- ^ Points
                    -> Maybe (V.Vector a) -- ^ Params.  Approximate if Nothing
                    -> Int                -- ^ Maximum iterations
                    -> (CubicBezier a, a) -- ^ result curve and maximum error
approximateCubic curve pts mbTs maxiter =
  let ts = fromMaybe (approximateParams (cubicC0 curve) (cubicC3 curve) pts) mbTs
      curve2 = fromMaybe curve $ lsqDist curve pts ts
      (bt, bt') = V.unzip $ V.map (evalBezierDeriv curve2) ts
      err = V.maximum $ V.zipWith vectorDistance pts bt
      (c, _, _, err2, _) =
        fromMaybe (curve2, ts, undefined, err, undefined) $
        approximateCubic' curve2 pts ts maxiter err bt bt'
  in (c, err2)
{-# SPECIALIZE approximateCubic :: CubicBezier Double -> V.Vector (Point Double)
  -> Maybe (V.Vector Double) -> Int -> (CubicBezier Double, Double) #-}

-- find (a, b) which minimises ∑ᵢ(a*aᵢ + b*bᵢ + epsᵢ)²
leastSquares :: (V.Unbox a, Fractional a, Eq a) =>
                V.Vector a -> V.Vector a -> V.Vector a -> Maybe (a, a)
leastSquares as bs epses = solveLinear2x2 a b c d e f
  where
    square x = x*x
    a = V.sum $ V.map square as
    b = V.sum $ V.zipWith (*) as bs
    c = V.sum $ V.zipWith (*) as epses
    d = b
    e = V.sum $ V.map square bs
    f = V.sum $ V.zipWith (*) bs epses
{-# SPECIALIZE leastSquares ::V.Vector Double -> V.Vector Double -> V.Vector Double -> Maybe (Double, Double) #-}

-- find the least squares between the points pᵢ and B(tᵢ) for
-- bezier curve B, where pts contains the points pᵢ and ts
-- the values of tᵢ .
-- The tangent at the beginning and end is maintained.
-- Since the start and end point remains the same,
-- we need to find the new value of p2' = p1 + α₁ * (p2 - p1)
-- and p₃' = p4 + α2 * (p3 - p4)
-- minimizing (∑|B(tᵢ) - pᵢ|²) gives a linear equation
-- with two unknown values (α₁ and α₂)
lsqDist :: (V.Unbox a, Fractional a, Eq a) =>
           CubicBezier a
           -> V.Vector (Point a) -> V.Vector a -> Maybe (CubicBezier a)
lsqDist (CubicBezier (Point !p1x !p1y) (Point !p2x !p2y) (Point !p3x !p3y) (Point !p4x !p4y)) pts ts = let
  calcParams t (Point px py) = let
    t2 = t * t; t3 = t2 * t
    ax = 3 * (p2x - p1x) * (t3 - 2 * t2 + t)
    ay = 3 * (p2y - p1y) * (t3 - 2 * t2 + t)
    bx = 3 * (p3x - p4x) * (t2 - t3)
    by = 3 * (p3y - p4y) * (t2 - t3)
    cx = (p4x - p1x) * (3 * t2 - 2 * t3) + p1x - px
    cy = (p4y - p1y) * (3 * t2 - 2 * t3) + p1y - py
    in (ax * ax + ay * ay,
        ax * bx + ay * by,
        ax * cx + ay * cy,
        bx * ax + by * ay,
        bx * bx + by * by,
        bx * cx + by * cy)
  add6 (!a,!b,!c,!d,!e,!f) (!a',!b',!c',!d',!e',!f') =
    (a+a',b+b',c+c',d+d',e+e',f+f')
  ( as, bs, cs, ds, es, fs ) = V.foldl1' add6 $ V.zipWith calcParams ts pts
  in do (alpha1, alpha2) <- solveLinear2x2 as bs cs ds es fs 
        let cp1 = Point (alpha1 * (p2x - p1x) + p1x) (alpha1 * (p2y - p1y) + p1y)
            cp2 = Point (alpha2 * (p3x - p4x) + p4x) (alpha2 * (p3y - p4y) + p4y)
        Just $ CubicBezier (Point p1x p1y) cp1 cp2 (Point p4x p4y)
{-# SPECIALIZE lsqDist :: CubicBezier Double
           -> V.Vector (Point Double) -> V.Vector Double -> Maybe (CubicBezier Double) #-}

-- calculate the least Squares bezier curve by choosing approximate values
-- of t, and iterating again with an improved estimate of t, by taking the
-- the values of t for which the points are closest to the curve
approximateCubic' :: (V.Unbox a, Ord a, Floating a) =>
                     CubicBezier a
                  -> V.Vector (Point a) -> V.Vector a
                  -> Int -> a -> V.Vector (Point a)
                  -> V.Vector (Point a)
                  -> Maybe (CubicBezier a, V.Vector a, V.Vector a, a, V.Vector (Point a))
approximateCubic' (CubicBezier p1 p2 p3 p4) pts ts maxiter err bt bt' = do
  let dir1 = V.map (($* (p2^-^p1)) . rotateVec . flipVector) bt'
      dir2 = V.map (($* (p3^-^p4)) . rotateVec . flipVector) bt'
      ps = V.zipWith3 (\b b' p ->
                        rotateVec (flipVector b') $*
                        (p^-^b)) bt bt' pts
      errs = V.map (negate.pointY) ps
      as = V.zipWith (\d t -> 3*pointY d*(1-t)*(1-t)*t)
           dir1 ts
      bs = V.zipWith (\d t -> 3*pointY d*(1-t)*t*t)
           dir2 ts
  (a,b) <- leastSquares as bs errs
  let newTs = V.zipWith5 (\t p d1 d2 b' ->
                           max 0 $ min 1 $
                           t + (pointX p - 3*(1-t)*t*(a*pointX d1*(1-t) +
                                                      b*pointX d2*t)) /
                           vectorMag b')
              ts ps dir1 dir2 bt'
      newCurve = CubicBezier p1 (p2 ^+^ a*^(p2^-^p1)) (p3 ^+^ b*^(p3^-^p4)) p4
      (bt2,bt2') = V.unzip $ V.map (evalBezierDeriv newCurve) newTs
      err2 = V.zipWith vectorDistance pts bt2
      maxErr = V.maximum err2
      -- alternative method for finding the t values:
      -- newTs = V.zipWith (-) ts (V.zipWith (calcDeltaT newCurve) pts ts)
  if maxiter < 1 || abs(err - maxErr) <= err/8
    then return (newCurve, newTs, err2, maxErr, bt2)
    else approximateCubic' newCurve pts newTs (maxiter-1) maxErr bt2 bt2'
{-# SPECIALIZE approximateCubic' ::
  CubicBezier Double -> V.Vector (Point Double) -> V.Vector Double
  -> Int -> Double -> V.Vector (Point Double)
  -> V.Vector (Point Double)
  -> Maybe (CubicBezier Double, V.Vector Double, V.Vector Double, Double, V.Vector (Point Double)) #-}


-- approximate t by calculating the distances between all points
-- and dividing by the total sum
approximateParams :: (V.Unbox a, Floating a) =>
                     Point a -> Point a -> V.Vector (Point a) -> V.Vector a
approximateParams start end pts 
  | V.null pts = V.empty
  | otherwise =
      let dists = V.generate (V.length pts)
                  (\i -> if i == 0
                         then vectorDistance start (V.unsafeIndex pts 0)
                         else vectorDistance (V.unsafeIndex pts (i-1)) (V.unsafeIndex pts i))
          total = V.sum dists + vectorDistance (V.last pts) end
      in V.map (/ total) $ V.scanl1 (+) dists
{-# SPECIALIZE approximateParams ::
   Point Double -> Point Double -> V.Vector (Point Double) -> V.Vector Double #-}

-- Alternative method for finding the next t values, using
-- Newton-Rafphson.  There is no noticable difference in speed or
-- efficiency.

-- calcDeltaT :: CubicBezier -> Point -> Double -> Double
-- calcDeltaT curve (Point !ptx !pty) t = let
--   (Point bezx bezy, Point dbezx dbezy, Point ddbezx ddbezy, _) = evalBezierDerivs curve t
--   in ((bezx - ptx) * dbezx + (bezy - pty) * dbezy) /
--      (dbezx * dbezx + dbezy * dbezy + (bezx - ptx) * ddbezx + (bezy - pty) * ddbezy)

