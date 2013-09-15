{-# LANGUAGE BangPatterns #-}
module Math.BernsteinPoly
       (BernsteinPoly(..), bernsteinSubsegment, listToBernstein, zeroPoly, (~*), (*~), (~+),
        (~-), degreeElevate, bernsteinSplit, bernsteinEval,
        bernsteinEvalDerivs, bernsteinDeriv)
       where

import Data.List

data BernsteinPoly = BernsteinPoly {
  bernsteinDegree :: !Int,
  bernsteinCoeffs :: ![Double] }
                   deriving Show

infixl 7 ~*, *~
infixl 6 ~+, ~-

-- | Create a bernstein polynomail from a list of coëfficients.
listToBernstein :: [Double] -> BernsteinPoly
listToBernstein [] = zeroPoly
listToBernstein l = BernsteinPoly (length l - 1) l

-- | The constant zero.
zeroPoly :: BernsteinPoly
zeroPoly = BernsteinPoly 0 [0]

-- | Return the subsegment between the two parameters.
bernsteinSubsegment :: BernsteinPoly -> Double -> Double -> BernsteinPoly
bernsteinSubsegment b t1 t2 
  | t1 > t2   = bernsteinSubsegment b t2 t1
  | otherwise = snd $ flip bernsteinSplit (t1/t2) $
                fst $ bernsteinSplit b t2

-- multiply two bezier curves
-- control point i from the product of beziers P * Q
-- is sum (P_j * Q_k) where j + k = i+1

-- | Multiply two bernstein polynomials.  The final degree
-- will be the sum of either degrees.  This operation takes O((n+m)^2)
-- with n and m the degree of the beziers.

(~*) :: BernsteinPoly -> BernsteinPoly -> BernsteinPoly
(BernsteinPoly la a) ~* (BernsteinPoly lb b) =
  BernsteinPoly (la+lb) $
  zipWith (flip (/)) (binCoeff (la + lb)) $
                 init $ map sum $
                 zipWith (zipWith (*)) (repeat a') (down b') ++
                 zipWith (zipWith (*)) (tail $ tails a') (repeat $ reverse b')
  where down l = tail $ scanl (flip (:)) [] l -- [[1], [2, 1], [3, 2, 1], ...
        a' = zipWith (*) a (binCoeff la)
        b' = zipWith (*) b (binCoeff lb)


-- find the binomial coefficients of degree n.
binCoeff :: Int -> [Double]
binCoeff n = map fromIntegral $
             scanl (\x m -> x * (n-m+1) `quot` m) 1 [1..n]

-- | Degree elevate a bernstein polynomail a number of times.
degreeElevate :: BernsteinPoly -> Int -> BernsteinPoly
degreeElevate b 0 = b
degreeElevate (BernsteinPoly lp p) times =
  degreeElevate (BernsteinPoly (lp+1) (head p:inner p 1)) (times-1)
  where
    inner []  _ = error "empty bernstein coefficients"
    inner [a] _ = [a]
    inner (a:b:rest) i =
      (i*a/fromIntegral lp + b*(1 - i/fromIntegral lp))
      : inner (b:rest) (i+1)


-- | Evaluate the bernstein polynomial.
bernsteinEval :: BernsteinPoly -> Double -> Double
bernsteinEval (BernsteinPoly lp [b]) _ = b
bernsteinEval (BernsteinPoly lp (b':bs)) t = go t n (b'*u) 2 bs
  where u = 1-t
        n = fromIntegral lp
        go !tn !bc !tmp _  [b] = tmp + tn*bc*b
        go !tn !bc !tmp !i (b:rest) =
          go (tn*t)         -- tn
          (bc*(n-i+1)/i)    -- bc
          ((tmp + tn*bc*b)*u) -- tmp
          (i+1)             -- i
          rest

-- | Evaluate the bernstein polynomial and its derivatives.
bernsteinEvalDerivs :: BernsteinPoly -> Double -> [Double]
bernsteinEvalDerivs b t
  | bernsteinDegree b == 0 = [bernsteinEval b t]
  | otherwise = bernsteinEval b t :
                bernsteinEvalDerivs (bernsteinDeriv b) t

-- | Find the derivative of a bernstein polynomial.
bernsteinDeriv :: BernsteinPoly -> BernsteinPoly
bernsteinDeriv (BernsteinPoly 0 _) = zeroPoly
bernsteinDeriv (BernsteinPoly lp p) =
  BernsteinPoly (lp-1) $
  map (* fromIntegral lp) $ zipWith (-) (tail p) p

-- | Split a bernstein polynomial
bernsteinSplit :: BernsteinPoly -> Double -> (BernsteinPoly, BernsteinPoly)
bernsteinSplit (BernsteinPoly lp p) t =
  (BernsteinPoly lp $ map head controls,
   BernsteinPoly lp $ reverse $ map last controls)
  where
    interp a b = (1-t)*a + t*b
    terp [_] = []
    terp l = let ctrs = zipWith interp l (tail l)
             in ctrs : terp ctrs
    controls = p:terp p

-- | Sum two bernstein polynomials.  The final degree will be the maximum of either
-- degrees.
(~+) :: BernsteinPoly -> BernsteinPoly -> BernsteinPoly
ba@(BernsteinPoly la a) ~+ bb@(BernsteinPoly lb b)
  | la < lb = BernsteinPoly lb $
              zipWith (+) (bernsteinCoeffs $ degreeElevate ba $ lb-la) b
  | la > lb = BernsteinPoly la $
              zipWith (+) a (bernsteinCoeffs $ degreeElevate bb $ la-lb)
  | otherwise = BernsteinPoly la $
                zipWith (+) a b

-- | Subtract two bernstein polynomials.  The final degree will be the maximum of either
-- degrees.
(~-) :: BernsteinPoly -> BernsteinPoly -> BernsteinPoly
ba@(BernsteinPoly la a) ~- bb@(BernsteinPoly lb b)
  | la < lb = BernsteinPoly lb $
              zipWith (-) (bernsteinCoeffs $ degreeElevate ba (lb-la)) b
  | la > lb = BernsteinPoly la $
              zipWith (-) a (bernsteinCoeffs $ degreeElevate bb (la-lb))
  | otherwise = BernsteinPoly la $
                zipWith (-) a b

-- | Scale a bernstein polynomial by a constant.
(*~) :: Double -> BernsteinPoly -> BernsteinPoly
a *~ (BernsteinPoly lb b) = BernsteinPoly lb (map (*a) b)
