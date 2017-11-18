{-# LANGUAGE BangPatterns #-}
-- | Numerical computations used by the cubic bezier functions.  Also
-- contains functions that aren't used anymore, but might be useful on
-- its own.
module Geom2D.CubicBezier.Numeric
       (quadraticRoot, cubicRoot, cubicRootNorm, solveLinear2x2, 
        goldSearch, 
        makeSparse, SparseMatrix, sparseMulT, sparseMul,
        addMatrix, addVec, lsqMatrix, lsqSolveDist,
        decompLDL, lsqSolve,
        solveTriDiagonal, solveCyclicTriD)
       where
import qualified Data.Vector.Unboxed as V
import Data.Vector.Unboxed.Mutable as MV
import Data.Matrix.Unboxed as M
import qualified Data.Matrix.Generic as G
import qualified Data.Matrix.Unboxed.Mutable as MM
import Control.Monad.ST
import Control.Monad

sign :: (Ord a, Num a, Num t) => a -> t
sign x | x < 0 = -1
       | otherwise = 1

-- | @quadraticRoot a b c@ find the real roots of the quadratic equation
-- @a x^2 + b x + c = 0@.  It will return one, two or zero roots.

quadraticRoot :: Double -> Double -> Double -> [Double]
quadraticRoot a b c
  | a == 0 && b == 0 = []
  | a == 0 = [-c/b]
  | otherwise = result
  where
    d = b*b - 4*a*c
    q = - (b + sign b * sqrt d) / 2
    x1 = q/a
    x2 = c/q
    result | d < 0     = []
           | d == 0    = [x1]
           | otherwise = [x1, x2]

cubicRoot :: Double -> Double -> Double -> Double -> [Double]
cubicRoot a b c d
  | a == 0 = quadraticRoot b c d
  | otherwise = cubicRootNorm (b/a) (c/a) (d/a)

-- | Find real roots from the normalized cubic equation.
cubicRootNorm :: Double -> Double -> Double -> [Double]
cubicRootNorm a b c
  | r2 < q3 = [m2sqrtQ * cos(t/3) - a/3,
               m2sqrtQ * cos((t+2*pi)/3) - a/3,
               m2sqrtQ * cos((t - 2*pi)/3) - a/3]
  | otherwise = [d+e-a/3]
  where q = (a*a - 3*b) / 9
        q3 = q*q*q
        r = (2*a*a*a - 9*a*b + 27*c) / 54
        t = acos(r/sqrt q3)
        m2sqrtQ = -2*sqrt q
        r2 = r*r
        d = - sign(r)*((abs r + sqrt(r2-q3))**1/3)
        e | d == 0 = 0
          | otherwise = q/d
        

-- | @solveLinear2x2 a b c d e f@ solves the linear equation with two variables (x and y) and two systems:
-- 
-- >a x + b y + c = 0
-- >d x + e y + f = 0
-- 
-- Returns @Nothing@ if no solution is found.
solveLinear2x2 :: (Eq a, Fractional a) => a -> a -> a -> a -> a -> a -> Maybe (a, a)
solveLinear2x2 a b c d e f =
  case det of 0 -> Nothing
              _ -> Just ((c * e - b * f) / det, (a * f - c * d)  / det)
  where det = d * b - a * e
{-# SPECIALIZE solveLinear2x2 :: Double -> Double -> Double -> Double -> Double -> Double -> Maybe (Double, Double) #-}

phi :: (Floating a) => a
phi = (-1 + sqrt 5) / 2


goldSearch :: (Ord a, Floating a) => (a -> a) -> Int -> a
goldSearch f maxiter =
  goldSearch' f 0 x1 x2 1 (f 0) 
  (f x1) (f x2) (f 1) maxiter
    where x1 = 1 - phi
          x2 = phi
{-# SPECIALIZE goldSearch :: (Double -> Double) -> Int -> Double #-}          

goldSearch' :: (Ord a, Floating a) =>
               (a -> a) -> a -> a -> a ->
               a -> a -> a -> a -> a -> Int -> a
goldSearch' f x0 x1 x2 x3 y0 y1 y2 y3 maxiter
  | maxiter < 1 = snd $ maximum [(y0, x0), (y1, x1),
                                 (y2, x2), (y3, x3)]
  | y1 < y2 =
    let x25 = x1 + phi*(x3-x1)
        y25 = f x25
    in goldSearch' f x1 x2 x25 x3 y1 y2 y25 y3 (maxiter-1)
  | otherwise =
    let x05 = x2 + phi*(x0-x2)
        y05 = f x05
    in goldSearch' f x0 x05 x1 x2 y0 y05 y1 y2 (maxiter-1)


data SparseMatrix a =
  SparseMatrix (V.Vector Int)
  (V.Vector (Int, Int)) (M.Matrix a)
                      
makeSparse :: Unbox a => V.Vector Int
              -- ^ The column index of the first element of each row.
              -- Should be ascending in order.
              -> M.Matrix a
              -- ^ The adjacent coefficients in each row
              -> SparseMatrix a
              -- ^ A sparse matrix.
makeSparse v m = SparseMatrix v (sparseRanges v vars width) m
  where
    width = cols m
    vars = V.last v + width

-- give the range of (possibly) nonzero coefficients for each column.
-- The column indices are those of the dense matrix of which the
-- sparse is a representation.
sparseRanges :: V.Vector Int -> Int -> Int -> V.Vector (Int, Int)
sparseRanges v vars width = ranges 
  where
    height = V.length v
    ranges = V.scanl' nextRange (nextRange (0,0) 0) $
             V.enumFromN 1 (vars-1)
    nextRange (s,e) i = (nextStart s i, nextEnd e i)
    nextStart s i
      | s >= height = height
      | v `V.unsafeIndex` s + width <= i =
          nextStart (s+1) i
      | otherwise = s
    nextEnd e i
      | e >= height = height
      | v `V.unsafeIndex` e > i = e
      | otherwise = nextEnd (e+1) i


-- | Given a rectangular matrix M, calculate the symmetric square
-- matrix MᵀM which can be used to find a least squares solution to
-- the overconstrained system.
lsqMatrix :: (Num a, Unbox a) =>
             SparseMatrix a
             -- ^ The input system.
             -> Matrix a
             -- ^ The resulting symmetric matrix as a sparse matrix.
             -- The first element of each row is the element on the
             -- diagonal.

lsqMatrix (SparseMatrix rowStart ranges m)
  | V.length rowStart /= height =
    error "lsqMatrix: lengths don't match."
  | otherwise = M.generate (vars, width) coeff
  where
    (height, width) = dim m
    vars = V.last rowStart + width
    overlap (s1,e1) (s2, e2) =
      (max s1 s2, min e1 e2)
    realIndex (r, c) =
      (r, c - rowStart `V.unsafeIndex` r)
    coeff (r,c) = let
      (s, e) | r+c >= vars = (0, 0)
             | otherwise =
                 overlap (ranges `V.unsafeIndex` r) (ranges `V.unsafeIndex` (r+c))
      in V.foldl' (\acc i -> acc + m `M.unsafeIndex` realIndex (i, r) *
                             m `M.unsafeIndex` realIndex (i, r+c)) 0 $
         V.enumFromN s (e-s)
{-# SPECIALIZE lsqMatrix :: SparseMatrix Double -> M.Matrix Double #-}

addMatrix :: (Num a, Unbox a) => M.Matrix a -> M.Matrix a -> M.Matrix a
addMatrix = M.zipWith (+)
{-# SPECIALIZE addMatrix :: M.Matrix Double -> M.Matrix Double -> M.Matrix Double #-}

addVec :: (Num a, Unbox a) => V.Vector a -> V.Vector a -> V.Vector a
addVec = V.zipWith (+)
{-# SPECIALIZE addVec :: V.Vector Double -> V.Vector Double -> V.Vector Double #-}

-- | Multiply the vector by the transpose of the sparse matrix.
sparseMulT :: (Num a, Unbox a) =>
              V.Vector a
              -> SparseMatrix a
              -> V.Vector a
sparseMulT v (SparseMatrix rowStart ranges m)
  | V.length v /= height =
    error "sparseMulT: lengths don't match."
  | otherwise = V.generate vars coeff
  where (height, width) = dim m
        vars | V.null rowStart = 0
             | otherwise = V.unsafeLast rowStart + width
        realIndex (r, c) =
          (r, c - rowStart `V.unsafeIndex` r)
        coeff i =
          let (s, e) = ranges `V.unsafeIndex` i
          in V.foldl' (\acc j ->
                        acc + m `M.unsafeIndex` realIndex (j, i) *
                        v `V.unsafeIndex` j) 0 $
             V.enumFromN s (e-s)
{-# SPECIALIZE sparseMulT :: V.Vector Double -> SparseMatrix Double -> V.Vector Double #-}

-- | Sparse matrix * vector multiplication.
sparseMul :: (Num a, Unbox a) =>
              SparseMatrix a
              -> V.Vector a
              -> V.Vector a
sparseMul (SparseMatrix rowStart _ranges m) v
  | V.length v /= vars =
    error "sparseMulT: lengths don't match."
  | otherwise = V.generate height coeff
  where (height, width) = dim m
        vars | V.null rowStart = 0
             | otherwise = V.unsafeLast rowStart + width
        coeff i = V.sum $ V.zipWith (*)
                  (V.unsafeSlice (rowStart V.! i) width v)
                  (G.unsafeTakeRow m i)
{-# SPECIALIZE sparseMul :: SparseMatrix Double -> V.Vector Double -> V.Vector Double #-}

-- | LDL* decomposition of the sparse hermitian matrix.  The
-- first element of each row is the diagonal component of the D
-- matrix.  The following elements are the elements next to the
-- diagonal in the L* matrix (the diagonal components in L* are 1).
-- For efficiency it mutates the matrix inplace.
decompLDL :: (Fractional a, Unbox a) => M.Matrix a -> M.Matrix a
decompLDL m = runST $ do
  m2 <- M.thaw m
  let (vars, width) = dim m
  V.forM_ (V.enumFromN 0 $ vars-1) $
    \startr -> do
      pivot <- MM.unsafeRead m2 (startr, 0)
      V.forM_ (V.enumFromN 1 $ width-1) $
        \c -> do
          el <- MM.unsafeRead m2 (startr, c)
          MM.unsafeWrite m2 (startr, c) (el/pivot)
      V.forM_ (V.enumFromN 0 $ min (width-1) $ vars-startr-1) $
        \r -> do
         r0 <- MM.unsafeRead m2 (startr, r+1)
         V.forM_ (V.enumFromN 0 (width-r-1)) $
              \c -> do r1 <- MM.unsafeRead m2 (startr, r+c+1)
                       el <- MM.unsafeRead m2 (r+startr+1, c)
                       MM.unsafeWrite m2 (r+startr+1, c)
                         (el - r0*r1*pivot)
  M.unsafeFreeze m2
{-# SPECIALIZE decompLDL :: Matrix Double -> Matrix Double #-}

solveLDL :: (Fractional a, Unbox a) =>
            M.Matrix a -> V.Vector a -> V.Vector a
solveLDL m v
  | rows m /= V.length v = error "solveLDL: lengths don't match"
  | otherwise = runST $ do
      let (vars, width) = M.dim m
      sol1 <- MV.new vars
      -- forward substitution on the first (width) rows
      V.forM_ (V.enumFromN 0 $ min vars width) $
        \i -> do
          let vi = v `V.unsafeIndex` i
          s <- liftM (V.foldl' (-) vi) $
               V.forM (V.enumFromN 0 i) $
               \j -> liftM ((m `M.unsafeIndex` (j, i-j)) *)
                     (MV.unsafeRead sol1 j)
          MV.unsafeWrite sol1 i s
          
      -- forward substitution on the next (height-width) rows
      V.forM_ (V.enumFromN width $ vars - width) $
        \i -> do
          let vi = v `V.unsafeIndex` i
          s <- liftM (V.foldl' (-) vi) $
               V.forM (V.enumFromN 1 (width-1)) $
               \j -> liftM ((m `M.unsafeIndex` (i-j, j)) *)
                     (MV.unsafeRead sol1 $ i-j)
          MV.unsafeWrite sol1 i s
          
      -- backward substitution on the last (width) rows
      V.forM_ (V.enumFromN 0 $ min vars width) $
        \i -> do
          solI <- MV.unsafeRead sol1 (vars-i-1)
          let d = m `M.unsafeIndex` (vars-i-1, 0)
          s <- liftM (V.foldl' (-) (solI/d)) $
               V.forM (V.enumFromN 0 i) $
               \j -> liftM ((m `M.unsafeIndex` (vars-i-1, j+1)) *)
                     (MV.unsafeRead sol1 $ vars-i+j)
          MV.unsafeWrite sol1 (vars-i-1) s
          
      -- backward substitution on the prevous (vars-width) rows
      V.forM_ (V.enumFromN width $ vars - width) $
        \i -> do
          solI <- MV.unsafeRead sol1 (vars-i-1)
          let d = m `M.unsafeIndex` (vars-i-1, 0)
          s <- liftM (V.foldl' (-) (solI/d)) $
               V.forM (V.enumFromN 0 (width-1)) $
               \j -> liftM ((m `M.unsafeIndex` (vars-i-1, j+1)) *)
                     (MV.unsafeRead sol1 $ vars-i+j)
          MV.unsafeWrite sol1 (vars-i-1) s
          
      V.unsafeFreeze sol1
{-# SPECIALIZE solveLDL :: M.Matrix Double -> V.Vector Double -> V.Vector Double #-}
    
-- | @lsqSolve rowStart M y@ Find a least squares solution x to the
-- system xM = y.
lsqSolve :: (Fractional a, Unbox a) =>
            SparseMatrix a    -- ^ sparse matrix
         -> V.Vector a        -- ^ Right hand side vector.
         -> V.Vector a        -- ^ Solution vector
lsqSolve m@(SparseMatrix _ _ m') v
  | rows m' /= V.length v = error "lsqSolve: lengths don't match"
  | otherwise = solveLDL m2 v2
  where
    v2 = sparseMulT v m
    m2 = decompLDL $ lsqMatrix m
{-# SPECIALIZE lsqSolve :: SparseMatrix Double -> V.Vector Double -> V.Vector Double #-}

-- | @lsqSolveDist rowStart M y@ Find a least squares solution of the distance between the points.
lsqSolveDist :: (Fractional a, Unbox a) =>
                SparseMatrix (a, a) -- ^ sparse matrix
             -> V.Vector (a, a)     -- ^ Right hand side vector.
             -> V.Vector a          -- ^ Solution vector
lsqSolveDist (SparseMatrix r s m') v
  | rows m' /= V.length v = error "lsqSolve: lengths don't match"
  | otherwise = solveLDL m3 v3
  where
    v3 = sparseMulT v1 m1 `addVec` sparseMulT v2 m2
    m3 = decompLDL $ lsqMatrix m1 `addMatrix` lsqMatrix m2
    (v1, v2) = V.unzip v
    (m1', m2') = M.unzip m'
    m1 = SparseMatrix r s m1'
    m2 = SparseMatrix r s m2'
{-# SPECIALIZE lsqSolveDist :: SparseMatrix (Double, Double) -> V.Vector (Double, Double) -> V.Vector Double #-}

-- | solve a tridiagonal system.  see metafont the program: ¶ 283
solveTriDiagonal :: (Unbox a, Fractional a) =>
                    (a, a, a) -> V.Vector (a, a, a, a) -> V.Vector a
solveTriDiagonal (!b0, !c0, !d0) rows_ = solutions
  where
    twovars = V.scanl nextrow (c0/b0, d0/b0) rows_
    solutions = V.scanr nextsol vn (V.unsafeInit twovars)
    vn = snd $ V.unsafeLast twovars
    nextsol (u, v) ti = v - u*ti
    nextrow (u, v) (ai, bi, ci, di) =
      (ci/(bi - u*ai), (di - v*ai)/(bi - u*ai))
{-# SPECIALIZE solveTriDiagonal :: (Double, Double, Double) -> V.Vector (Double, Double, Double, Double) -> V.Vector Double #-}

-- | solve a cyclic tridiagonal system.  see metafont the program: ¶ 286
solveCyclicTriD :: (Unbox a, Fractional a) => V.Vector (a, a, a, a) -> V.Vector a
solveCyclicTriD rows_ = solutions
  where
    threevars = V.tail $ V.scanl nextrow (0, 0, 1) rows_
    nextrow (!u, !v, !w) (!ai, !bi, !ci, !di) =
      (ci/(bi - ai*u), (di - ai*v)/(bi - ai*u), -ai*w/(bi - ai*u))
    (totvn, totwn) = V.foldr (\(u, v, w) (v', w') ->
                               (v - u*v', w - u*w'))
                     (0, 1) (V.init threevars)
    t0 = (vn - un*totvn) / (1 - (wn - un*totwn))
    (!un, !vn, !wn) = V.last threevars
    solutions = V.scanl nextsol t0 $ V.cons (un, vn, wn) $
                V.init $ V.init $ threevars
    nextsol t (!u, !v, !w) = (v + w*t0 - t)/u
{-# SPECIALIZE solveCyclicTriD :: V.Vector (Double, Double, Double, Double) -> V.Vector Double #-}
