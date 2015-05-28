{-# LANGUAGE BangPatterns #-}
module Geom2D.CubicBezier.Approximate (
  approximatePath, approximatePathMax, approximateCurve, approximateCurveWithParams)
       where
import Geom2D
import Geom2D.CubicBezier.Numeric
import Geom2D.CubicBezier.Basic
import Data.Function
import Data.List
import Data.Maybe
import qualified Data.Map as M

interpolate :: Double -> Double -> Double -> Double
interpolate a b x = (1-x)*a + x*b

-- | Approximate a function with piecewise cubic bezier splines using
-- a least-squares fit, within the given tolerance.  Each subcurve is
-- approximated by using a finite number of samples.  It is recommended
-- to avoid changes in direction by subdividing the original function
-- at points of inflection.

approximatePath :: (Double -> (DPoint, DPoint)) -- ^ The function to approximate and it's derivative
                -> Double                     -- ^ The number of discrete samples taken to approximate each subcurve
                -> Double                     -- ^ The tolerance
                -> Double                     -- ^ The lower parameter of the function      
                -> Double                     -- ^ The upper parameter of the function
                -> [CubicBezier Double]
approximatePath f n tol tmin tmax
  | err <= tol = [cb_out]
  | otherwise = approximatePath f n tol tmin terr ++
                approximatePath f n tol terr tmax
  where
    (cb_out, terr', err) = approximateCurveWithParams curveCb
                           points ts tol
    terr = interpolate tmin tmax terr'
    ts        = [i/(n+1) | i <- [1..n]]
    points    = map (fst . f . interpolate tmin tmax) ts
    (t0, t0') = f tmin
    (t1, t1') = f tmax
    curveCb = CubicBezier t0 (t0^+^t0') (t1^-^t1') t1


-- | Like approximatePath, but limit the number of subcurves.
approximatePathMax :: Int                        -- ^ The maximum number of subcurves
                   -> (Double -> (DPoint, DPoint)) -- ^ The function to approximate and it's derivative
                   -> Double                     -- ^ The number of discrete samples taken to approximate each subcurve
                   -> Double                     -- ^ The tolerance
                   -> Double                     -- ^ The lower parameter of the function      
                   -> Double                     -- ^ The upper parameter of the function
                   -> [CubicBezier Double]
approximatePathMax m f n tol tmin tmax =
  approxMax f tol m ts segments
  where segments              = M.singleton err (FunctionSegment tmin tmax t_err outline)
        (p0, p0') = f tmin
        (p1, p1') = f tmax
        ts = [i/(n+1) | i <- [1..n]]
        points = map (fst . f . interpolate tmin tmax) ts
        curveCb = CubicBezier p0 (p0^+^p0') (p1^-^p1') p1
        (outline, t_err', err) = approximateCurveWithParams curveCb
                                 points ts tol
        t_err = interpolate tmin tmax t_err'

data FunctionSegment = FunctionSegment {
  fs_t_min :: {-# UNPACK #-} !Double,  -- the least t param of the segment in the original curve
  _fs_t_max :: {-# UNPACK #-} !Double,  -- the max t param of the segment in the original curve
  _fs_t_err :: {-# UNPACK #-} !Double,  -- the param where the error is maximal
  fs_curve :: CubicBezier Double -- the curve segment
  }

-- Keep a map from maxError to FunctionSegment for each subsegment to keep
-- track of the segment with the maximum error.  This ensures a n
-- log(n) execution time, rather than n^2 when a list is used.
approxMax :: (Double -> (DPoint, DPoint)) -> Double -> Int
          -> [Double] -> M.Map Double FunctionSegment -> [CubicBezier Double]
approxMax f tol n ts segments
  | n < 1 = error "Minimum number of segments is one."
  | (n == 1) || (err < tol) =
    map fs_curve $ sortBy (compare `on` fs_t_min) $ map snd $ M.toList segments
  | otherwise = approxMax f tol (n-1) ts $
                M.insert err_l (FunctionSegment t_min t_err t_err_l curve_l) $
                M.insert err_r (FunctionSegment t_err t_max t_err_r curve_r)
                newSegments
  where
    ((err, FunctionSegment t_min t_max t_err _), newSegments) = M.deleteFindMax segments
    (fmin, fmin') = f t_min
    (fmid, fmid') = f t_err
    (fmax, fmax') = f t_max
    fcurve_l = CubicBezier fmin (fmin^+^fmin') (fmid^-^fmid') fmid
    fcurve_r = CubicBezier fmid (fmid^+^fmid') (fmax^-^fmax') fmax
    pointsl = map (fst . f . interpolate t_min t_err) ts
    pointsr = map (fst . f . interpolate t_err t_max) ts
    t_err_l = interpolate t_min t_err t_err_l'
    t_err_r = interpolate t_err t_max t_err_r'
    (curve_l, t_err_l', err_l)  = approximateCurveWithParams fcurve_l pointsl ts tol
    (curve_r, t_err_r', err_r)  = approximateCurveWithParams fcurve_r pointsr ts tol

-- | @approximateCurve b pts eps@ finds the least squares fit of a bezier
-- curve to the points @pts@.  The resulting bezier has the same first
-- and last control point as the curve @b@, and have tangents colinear with @b@.
-- return the curve, the parameter with maximum error, and maximum error.
-- Calculate to withing eps tolerance.

approximateCurve :: CubicBezier Double -> [DPoint] -> Double -> (CubicBezier Double, Double, Double)
approximateCurve curve@(CubicBezier p1 _ _ p4) pts eps =
  approximateCurveWithParams curve pts (approximateParams curve p1 p4 pts) eps

-- | Like approximateCurve, but also takes an initial guess of the
-- parameters closest to the points.  This might be faster if a good
-- guess can be made.

approximateCurveWithParams :: CubicBezier Double -> [DPoint] -> [Double] -> Double -> (CubicBezier Double, Double, Double)
approximateCurveWithParams curve pts ts eps =
  let (c, newTs) = fromMaybe (curve, ts) $
                   approximateCurve' curve pts ts 40 (bezierParamTolerance curve eps) 1
      curvePts   = map (evalBezier c) newTs
      distances  = zipWith vectorDistance pts curvePts
      (t, maxError) = maximumBy (compare `on` snd) (zip ts distances)
  in (c, t, maxError)

data LSParams = LSParams {-# UNPACK #-} !Double
                {-# UNPACK #-} !Double
                {-# UNPACK #-} !Double
                {-# UNPACK #-} !Double
                {-# UNPACK #-} !Double
                {-# UNPACK #-} !Double

addParams :: LSParams -> LSParams -> LSParams
addParams (LSParams a b c d e f) (LSParams a' b' c' d' e' f') =
  LSParams (a+a') (b+b') (c+c') (d+d') (e+e') (f+f')

-- find the least squares between the points p_i and B(t_i) for
-- bezier curve B, where pts contains the points p_i and ts
-- the values of t_i .
-- The tangent at the beginning and end is maintained.
-- Since the start and end point remains the same,
-- we need to find the new value of p2' = p1 + alpha1 * (p2 - p1)
-- and p3' = p4 + alpha2 * (p3 - p4)
-- minimizing (sum |B(t_i) - p_i|^2) gives a linear equation
-- with two unknown values (alpha1 and alpha2), which can be
-- solved easily
leastSquares :: CubicBezier Double -> [DPoint] -> [Double] -> Maybe (CubicBezier Double)
leastSquares (CubicBezier (Point !p1x !p1y) (Point !p2x !p2y) (Point !p3x !p3y) (Point !p4x !p4y)) pts ts = let
  calcParams t (Point px py) = let
    t2 = t * t; t3 = t2 * t
    ax = 3 * (p2x - p1x) * (t3 - 2 * t2 + t)
    ay = 3 * (p2y - p1y) * (t3 - 2 * t2 + t)
    bx = 3 * (p3x - p4x) * (t2 - t3)
    by = 3 * (p3y - p4y) * (t2 - t3)
    cx = (p4x - p1x) * (3 * t2 - 2 * t3) + p1x - px
    cy = (p4y - p1y) * (3 * t2 - 2 * t3) + p1y - py
    in LSParams
       (ax * ax + ay * ay)
       (ax * bx + ay * by)
       (ax * cx + ay * cy)
       (bx * ax + by * ay)
       (bx * bx + by * by)
       (bx * cx + by * cy)
  LSParams !a !b !c !d !e !f = foldl1' addParams (zipWith calcParams ts pts)
  in do (alpha1, alpha2) <- solveLinear2x2 a b c d e f
        let cp1 = Point (alpha1 * (p2x - p1x) + p1x) (alpha1 * (p2y - p1y) + p1y)
            cp2 = Point (alpha2 * (p3x - p4x) + p4x) (alpha2 * (p3y - p4y) + p4y)
        Just $ CubicBezier (Point p1x p1y) cp1 cp2 (Point p4x p4y)

offBez = CubicBezier {bezierC0 = Point (-0.9984603532054123) 0.6656402354702748, bezierC1 = Point 5.001539646794588 9.665640235470274, bezierC2 = Point (-0.15147186257614287) 7.848528137423857, bezierC3 = Point 5.848528137423857 1.8485281374238571}
offset = [(Point (-0.38561064422602453) 1.534174843357413,Point 5.459999999999999 7.23),(Point 0.20981740728624898 2.257123236825599,Point 5.040000000000001 5.520000000000001),(Point 0.8250761640484685 2.8465346207778452,Point 4.739999999999999 3.8699999999999974),(Point 1.5113436854000506 3.297312629199899,Point 4.56 2.2799999999999985),(Point 2.3027212152335723 3.5586727085985723,Point 4.499999999999999 0.7499999999999973),(Point 3.139155088167614 3.561315558394883,Point 4.560000000000001 (-0.720000000000002)),(Point 3.9078614426037217 3.327564900442083,Point 4.74 (-2.13)),(Point 4.58582886049917 2.9394762807229355,Point 5.039999999999997 (-3.4800000000000018)),(Point 5.217502221265897 2.4427069450968126,Point 5.459999999999999 (-4.770000000000001))]
Just lq = leastSquares offBez (map fst offset) ([0.1,0.2..0.9])

-- calculate the least Squares bezier curve by choosing approximate values
-- of t, and iterating again with an improved estimate of t, by taking the
-- the values of t for which the points are closest to the curve

approximateCurve' :: CubicBezier Double -> [DPoint] -> [Double] -> Int -> Double -> Double -> Maybe (CubicBezier Double, [Double])
approximateCurve' curve pts ts maxiter eps prevDeltaT = do
  newCurve <- leastSquares curve pts ts
  let deltaTs = zipWith (calcDeltaT newCurve) pts ts
      ts' = map (max 0 . min 1) $ zipWith (-) ts deltaTs
  newerCurve <- leastSquares curve pts ts'
  let deltaTs' = zipWith (calcDeltaT newerCurve) pts ts'
      newTs = interpolateTs ts ts' deltaTs deltaTs'
      thisDeltaT = maximum $ map abs $ zipWith (-) newTs ts
  if maxiter < 1 ||
     -- Because convergence may be slow initially, make sure it is converging:
     (prevDeltaT < eps/2  && thisDeltaT < prevDeltaT / 2)
    then do c <- leastSquares curve pts newTs
            return (c, newTs)
    else approximateCurve' curve pts newTs (maxiter - 1) eps thisDeltaT

-- improve convergence by making a better estimate for t
-- it is based on the observation that the ratio  
-- r = dt_2 / dt_1, with dt_2 = t_2 - t_1 and dt_1 = t_1 - t_0
-- for successive approximations of t changes little.
-- The infinite sum (dt_1 + dt_1 * r + dt_1 * r^2 + dt_1 * r^3 ...)
-- can easily be calculated by dt_1 * (1 / (1 - r))
-- which becomes dt_1^2 / (dt_1 - dt_2)
-- Only do this if it appears to converge for all values of t
-- If the value of t changes too much keep the old value.
-- This improves the convergence by a factor of about 10
interpolateTs :: [Double] -> [Double] -> [Double] -> [Double] -> [Double]
interpolateTs ts ts' deltaTs deltaTs' =
  map (max 0 . min 1) (
    if and $ zipWith (\dT dT' -> dT * dT' > 0 && dT' / dT < 1) deltaTs deltaTs'
    then zipWith3 (\t dT dT' -> let
                      newDt = (dT * dT / (dT - dT'))
                      in t - (if abs newDt > 0.2 then dT' else newDt)) ts deltaTs deltaTs'
    else zipWith (-) ts' deltaTs')

-- approximate t by calculating the distances between all points
-- and dividing by the total sum
approximateParams :: CubicBezier Double -> DPoint -> DPoint -> [DPoint] -> [Double]
approximateParams cb start end pts = let
  segments = start : (pts ++ [end])
  dists = zipWith vectorDistance segments (tail segments)
  total = sum dists
  improve p t = t - calcDeltaT cb p t
  in zipWith improve pts $ map (/ total) $ scanl1 (+) dists

-- find a value of t where B(t) is closer between the bezier curve and
-- the point (ptx, pty), by solving f' = 0 where
-- f(t) = (X(t) - x)^2 + (Y(t) - y)^2, the square of the distance between the bezier and the point
-- the reduction of t is one iteration of Newton Raphson:  f'(t)/f''(t)
-- using more iterations doesn't appear to give an improvement
-- See Curve Fitting with Piecewise Parametric Cubics by Stone & Plass
calcDeltaT :: CubicBezier Double -> DPoint -> Double -> Double
calcDeltaT curve (Point !ptx !pty) t = let
  (Point bezx bezy, Point dbezx dbezy, Point ddbezx ddbezy, _) = evalBezierDerivs curve t
  in ((bezx - ptx) * dbezx + (bezy - pty) * dbezy) /
     (dbezx * dbezx + dbezy * dbezy + (bezx - ptx) * ddbezx + (bezy - pty) * ddbezy)
