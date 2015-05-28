{-# LANGUAGE BangPatterns, FlexibleInstances, MultiParamTypeClasses, DeriveFunctor #-}
module Geom2D.CubicBezier.Basic
       (CubicBezier (..), PathJoin (..), Path (..), AffineTransform (..), 
        bezierParam, bezierParamTolerance, reorient, bezierToBernstein,
        evalBezier, evalBezierDeriv, evalBezierDerivs, findBezierTangent,
        bezierHoriz, bezierVert, findBezierInflection, findBezierCusp,
        arcLength, arcLengthParam, splitBezier, bezierSubsegment, splitBezierN,
        colinear)
       where
import Geom2D
import Geom2D.CubicBezier.Numeric
import Math.BernsteinPoly
import Numeric.Integration.TanhSinh

data CubicBezier a = CubicBezier {
  bezierC0 :: Point a,
  bezierC1 :: Point a,
  bezierC2 :: Point a,
  bezierC3 :: Point a} deriving (Show, Functor)

data PathJoin a = JoinLine (Point a) | JoinCurve (Point a) (Point a)
              deriving (Show, Functor)
data Path a = OpenPath [(Point a, PathJoin a)] (Point a)
            | ClosedPath [(Point a, PathJoin a)]
            deriving (Show, Functor)

instance (Num a) => AffineTransform (CubicBezier a) a where
  {-# SPECIALIZE transform :: Transform Double -> CubicBezier Double -> CubicBezier Double #-}
  transform t (CubicBezier c0 c1 c2 c3) =
    CubicBezier (transform t c0) (transform t c1) (transform t c2) (transform t c3)

-- | Return True if the param lies on the curve, iff it's in the interval @[0, 1]@.
bezierParam :: (Ord a, Num a) => a -> Bool
bezierParam t = t >= 0 && t <= 1
{-# SPECIALIZE bezierParam :: Double -> Bool #-}

-- | Convert a tolerance from the codomain to the domain of the bezier
-- curve, by dividing by the maximum velocity on the curve.  The
-- estimate is conservative, but holds for any value on the curve.
bezierParamTolerance :: CubicBezier Double -> Double -> Double
bezierParamTolerance (CubicBezier !p1 !p2 !p3 !p4) eps = eps / maxVel
  where 
    maxVel = 3 * max (vectorDistance p1 p2)
             (max (vectorDistance p2 p3)
              (vectorDistance p3 p4))

-- | Reorient to the curve B(1-t).
reorient :: CubicBezier a -> CubicBezier a
reorient (CubicBezier p0 p1 p2 p3) = CubicBezier p3 p2 p1 p0
{-# INLINE reorient #-}

-- | Give the bernstein polynomial for each coordinate.
bezierToBernstein :: CubicBezier Double -> (BernsteinPoly, BernsteinPoly)
bezierToBernstein (CubicBezier a b c d) = (listToBernstein $ map pointX coeffs,
                                           listToBernstein $ map pointY coeffs)
  where coeffs = [a, b, c, d]

-- | Calculate a value on the curve.
evalBezier :: Num a => CubicBezier a -> a -> Point a
evalBezier b = fst . evalBezierDeriv b
{-# SPECIALIZE evalBezier :: CubicBezier Double -> Double -> DPoint #-}

-- | Calculate a value and the first derivative on the curve.
evalBezierDeriv :: Num a => CubicBezier a -> a -> (Point a, Point a)
evalBezierDeriv cb t = (b,b')
  where
    (b,b',_,_) = evalBezierDerivs cb t
{-# SPECIALIZE evalBezierDeriv :: CubicBezier Double -> Double -> (DPoint, DPoint) #-}    
    
-- | Calculate a value and all derivatives on the curve.
evalBezierDerivs :: Num a => CubicBezier a -> a -> (Point a, Point a, Point a, Point a)
evalBezierDerivs (CubicBezier !a !b !c !d) t =
  (interp abbc bccd,
   3*^(bccd ^-^ abbc),
   6*^(cd ^-^ 2*^bc ^+^ ab),
   6*^(d ^+^ 3*^(b ^-^ c) ^-^ a))
  where
    mt = 1-t
    interp !v !w = mt*^v ^+^ t*^w
    ab = interp a b
    bc = interp b c
    cd = interp c d
    abbc = interp ab bc
    bccd = interp bc cd
{-# SPECIALIZE evalBezierDerivs :: CubicBezier Double -> Double -> (DPoint, DPoint, DPoint, DPoint) #-}
           
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

-- Use the formula B''x(t) * B'y(t) - B''y(t) * B'x(t) = 0
-- with B'x(t) the x value of the first derivative at t,
-- B''y(t) the y value of the second derivative at t
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

-- | Split a bezier curve into two curves.
splitBezier :: Num a => CubicBezier a -> a -> (CubicBezier a, CubicBezier a)
splitBezier (CubicBezier a b c d) t =
  let ab = interpolateVector a b t
      bc = interpolateVector b c t
      cd = interpolateVector c d t
      abbc = interpolateVector ab bc t
      bccd = interpolateVector bc cd t
      mid = interpolateVector abbc bccd t
  in (CubicBezier a ab abbc mid, CubicBezier mid bccd cd d)
{-# SPECIALIZE splitBezier :: CubicBezier Double -> Double -> (CubicBezier Double, CubicBezier Double) #-}     

-- | Return the subsegment between the two parameters.
bezierSubsegment :: (Fractional a, Ord a) => CubicBezier a -> a -> a -> CubicBezier a
bezierSubsegment b t1 t2 
  | t1 > t2   = bezierSubsegment b t2 t1
  | otherwise = snd $ flip splitBezier (t1/t2) $
                fst $ splitBezier b t2
{-# SPECIALIZE bezierSubsegment :: CubicBezier Double -> Double -> Double -> CubicBezier Double #-}

-- | Split a bezier curve into a list of beziers
-- The parameters should be in ascending order or
-- the result is unpredictable.
splitBezierN :: (Ord a, Fractional a) => CubicBezier a -> [a] -> [CubicBezier a]
splitBezierN c [] = [c]
splitBezierN c [t] = [a, b] where
  (a, b) = splitBezier c t
splitBezierN c (t:u:rest) =
  bezierSubsegment c 0 t :
  bezierSubsegment c t u :
  tail (splitBezierN c $ u:rest)
{-# SPECIALIZE splitBezierN :: CubicBezier Double -> [Double] -> [CubicBezier Double] #-}

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
