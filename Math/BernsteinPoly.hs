{-# LANGUAGE BangPatterns, ViewPatterns #-}
-- | Algebra on polynomials in the Bernstein form.  It is based on the
-- paper /Algebraic manipulation in the Bernstein form made simple via
-- convolutions/ by J. Sanchez-Reyes.  It's an efficient
-- implementation using the scaled basis, and using ghc rewrite rules
-- to eliminate intermediate polynomials.

module Math.BernsteinPoly
       (BernsteinPoly(..), bernsteinSubsegment, listToBernstein, zeroPoly,
        (~*), (*~), (~+), (~-), degreeElevate, bernsteinSplit, bernsteinEval,
        bernsteinEvalDeriv, binCoeff, convolve, bernsteinEvalDerivs, bernsteinDeriv)
       where
import Data.Vector.Unboxed as V
import qualified Data.Vector as B

data BernsteinPoly a = BernsteinPoly {
  bernsteinCoeffs :: V.Vector a}
                   deriving Show

data ScaledPoly a = ScaledPoly {
  scaledCoeffs :: V.Vector a }
                deriving Show
infixl 7 ~*, *~
infixl 6 ~+, ~-

{-# RULES "toScaled/fromScaled" forall a. toScaled (fromScaled a) = a;
  "fromScaled/toScaled" forall a. fromScaled (toScaled a) = a; #-}

toScaled :: (Unbox a, Num a) => BernsteinPoly a -> ScaledPoly a
toScaled (BernsteinPoly v) =
  ScaledPoly $
  V.zipWith (*) v $ binCoeff $ V.length v - 1
{-# NOINLINE[2] toScaled #-}

fromScaled :: (Unbox a, Fractional a) => ScaledPoly a -> BernsteinPoly a
fromScaled (ScaledPoly v) =
    BernsteinPoly $
    V.zipWith (/) v $ binCoeff $ V.length v - 1
{-# NOINLINE[2] fromScaled #-}

-- | Create a bernstein polynomail from a list of coëfficients.
listToBernstein :: (Unbox a, Num a) => [a] -> BernsteinPoly a
listToBernstein [] = zeroPoly
listToBernstein l = BernsteinPoly $ V.fromList l
{-# INLINE listToBernstein #-}

-- | The constant zero.
zeroPoly :: (Num a, Unbox a) => BernsteinPoly a
zeroPoly = BernsteinPoly $ V.fromList [0]
{-# SPECIALIZE zeroPoly :: BernsteinPoly Double #-}

-- | Return the subsegment between the two parameters.
bernsteinSubsegment :: (Unbox a, Ord a, Fractional a) =>
                       BernsteinPoly a -> a -> a -> BernsteinPoly a
bernsteinSubsegment b t1 t2 
  | t1 > t2   = bernsteinSubsegment b t2 t1
  | otherwise = snd $ flip bernsteinSplit (t1/t2) $
                fst $ bernsteinSplit b t2
{-# INLINE bernsteinSubsegment #-}                

-- | Calculate the convolution of two vectors.
convolve :: (Unbox a, Num a) => Vector a -> Vector a -> Vector a
convolve a b =
  V.map (\i -> V.sum $
               V.zipWith (*) a $
               V.reverse $
               V.unsafeTake i b)
  (V.enumFromN 1 $ V.length b)
  V.++ V.map (\i -> V.sum $
                    V.zipWith (*)
                    (V.unsafeDrop i a)
                    (V.reverse b))
  (V.enumFromN 1 $ V.length a-1)
{-# SPECIALIZE convolve :: Vector Double -> Vector Double -> Vector Double #-}

-- | Multiply two bernstein polynomials using convolution.  The final
-- degree will be the sum of either degrees.  This operation takes
-- O((n+m)^2) with n and m the degree of the beziers.  Note that
-- convolution can be done in O(n log n) using the FFT, which may be
-- faster for large polynomials.

(~*) :: (Unbox a, Fractional a) =>
        BernsteinPoly a -> BernsteinPoly a -> BernsteinPoly a
(toScaled -> a) ~* (toScaled -> b) =
  fromScaled $ mulScaled a b
{-# INLINE (~*) #-}  

mulScaled :: (Unbox a, Num a) => ScaledPoly a -> ScaledPoly a -> ScaledPoly a
mulScaled (ScaledPoly a) (ScaledPoly b) =
  ScaledPoly $ convolve a b
{-# INLINE mulScaled #-}    

-- | Give the binomial coefficients of degree n.
binCoeff :: (Num a, Unbox a) => Int -> V.Vector a
binCoeff n = V.map fromIntegral $
             V.scanl (\x m -> x * (n-m+1) `quot` m)
             1 (V.enumFromN 1 n)
{-# INLINE binCoeff #-}             

-- | Degree elevate a bernstein polynomial a number of times.
degreeElevateScaled :: (Unbox a, Num a)
                       => ScaledPoly a -> Int -> ScaledPoly a
degreeElevateScaled b@(ScaledPoly p) times
  | times <= 0 = b
  | otherwise = ScaledPoly $ convolve (binCoeff times) p
{-# SPECIALIZE degreeElevateScaled :: ScaledPoly Double ->
    Int -> ScaledPoly Double #-}                

degreeElevate :: (Unbox a, Fractional a)
                 => BernsteinPoly a -> Int -> BernsteinPoly a
degreeElevate (toScaled -> b) times =
  fromScaled (degreeElevateScaled b times)
{-# INLINE degreeElevate #-}  

-- | Evaluate the bernstein polynomial using the horner rule adapted
-- for bernstein polynomials.

bernsteinEval :: (Unbox a, Fractional a)
                 => BernsteinPoly a -> a -> a
bernsteinEval (BernsteinPoly v) _
  | V.length v == 0 = 0
bernsteinEval (BernsteinPoly v) _
  | V.length v == 1 = V.unsafeHead v
bernsteinEval (BernsteinPoly v) t =
  go t (fromIntegral n) (V.unsafeIndex v 0 * u) 1
  where u = 1-t
        n = fromIntegral $ V.length v - 1
        go !tn !bc !tmp !i
          | i == n = tmp + tn*V.unsafeIndex v n
          | otherwise =
            go (tn*t) -- tn
            (bc*fromIntegral (n-i)/(fromIntegral i + 1)) -- bc
            ((tmp + tn*bc*V.unsafeIndex v i)*u) -- tmp
            (i+1) -- i
{-# SPECIALIZE bernsteinEval :: BernsteinPoly Double -> Double -> Double #-}            
            
-- | Evaluate the bernstein polynomial and first derivative
bernsteinEvalDeriv :: (Unbox t, Fractional t) => BernsteinPoly t -> t -> (t,t)
bernsteinEvalDeriv b@(BernsteinPoly v) t
  | V.length v <= 1 = (V.unsafeHead v, 0)
  | otherwise = (bernsteinEval b t, bernsteinEval (bernsteinDeriv b) t)
{-# INLINE bernsteinEvalDeriv #-}                

            
-- | Evaluate the bernstein polynomial and its derivatives.
bernsteinEvalDerivs :: (Unbox t, Fractional t) => BernsteinPoly t -> t -> [t]
bernsteinEvalDerivs b@(BernsteinPoly v) t
  | V.length v <= 1 = [V.unsafeHead v, 0]
  | otherwise = bernsteinEval b t :
                bernsteinEvalDerivs (bernsteinDeriv b) t
{-# INLINE bernsteinEvalDerivs #-}                

-- | Find the derivative of a bernstein polynomial.
bernsteinDeriv :: (Unbox a, Num a) => BernsteinPoly a -> BernsteinPoly a
bernsteinDeriv (BernsteinPoly v)
  | V.length v == 0 = zeroPoly
bernsteinDeriv (BernsteinPoly v) =
  BernsteinPoly $
  V.map (* fromIntegral (V.length v - 1)) $
  V.zipWith (-) (V.tail v) v
{-# SPECIALIZE bernsteinDeriv :: BernsteinPoly Double ->
    BernsteinPoly Double #-}  

-- | Split a bernstein polynomial
bernsteinSplit :: (Unbox a, Num a) =>
                  BernsteinPoly a -> a -> (BernsteinPoly a, BernsteinPoly a)
bernsteinSplit (BernsteinPoly v) t =
  (BernsteinPoly $ convert $
   B.map V.head interpVecs,
   BernsteinPoly $ V.reverse $ convert $
   B.map V.last $ convert interpVecs)
  where
    interp a b = (1-t)*a + t*b
    interpVecs = B.iterateN (V.length v) interpVec v
    interpVec v2 = V.zipWith interp v2 (V.tail v2)
{-# SPECIALIZE bernsteinSplit :: BernsteinPoly Double -> Double ->
    (BernsteinPoly Double, BernsteinPoly Double) #-}

addScaled :: (Unbox a, Num a) => ScaledPoly a -> ScaledPoly a -> ScaledPoly a
addScaled ba@(ScaledPoly a) bb@(ScaledPoly b)
  | la < lb = ScaledPoly $
              V.zipWith (+) (scaledCoeffs $ degreeElevateScaled ba $ lb-la) b
  | la > lb = ScaledPoly $
              V.zipWith (+) a (scaledCoeffs $ degreeElevateScaled bb $ la-lb)
  | otherwise = ScaledPoly $ V.zipWith (+) a b
  where la = V.length a
        lb = V.length b
{-# SPECIALIZE addScaled :: ScaledPoly Double -> ScaledPoly Double ->
    ScaledPoly Double #-}        

-- | Sum two bernstein polynomials.  The final degree will be the
-- maximum of either degrees.
(~+) :: (Unbox a, Fractional a) =>
        BernsteinPoly a -> BernsteinPoly a -> BernsteinPoly a
(toScaled -> a) ~+ (toScaled -> b) = fromScaled $ addScaled a b
{-# INLINE (~+) #-}

subScaled :: (Unbox a, Num a) => ScaledPoly a -> ScaledPoly a -> ScaledPoly a
subScaled ba@(ScaledPoly a) bb@(ScaledPoly b)
  | la < lb = ScaledPoly $
              V.zipWith (-) (scaledCoeffs $ degreeElevateScaled ba $ lb-la) b
  | la > lb = ScaledPoly $
              V.zipWith (-) a (scaledCoeffs $ degreeElevateScaled bb $ la-lb)
  | otherwise = ScaledPoly $ V.zipWith (-) a b
  where la = V.length a
        lb = V.length b
{-# SPECIALIZE subScaled :: ScaledPoly Double -> ScaledPoly Double ->
    ScaledPoly Double #-}        

-- | Subtract two bernstein polynomials.  The final degree will be the
-- maximum of either degrees.
(~-) :: (Unbox a, Fractional a) =>
        BernsteinPoly a -> BernsteinPoly a -> BernsteinPoly a

(toScaled -> a) ~- (toScaled -> b) = fromScaled $ subScaled a b
{-# INLINE (~-) #-}

-- | Scale a bernstein polynomial by a constant.
(*~) :: (Unbox a, Num a) => a -> BernsteinPoly a -> BernsteinPoly a
a *~ (BernsteinPoly v) = BernsteinPoly (V.map (*a) v)
{-# INLINE (*~) #-}
