{-# LANGUAGE TypeFamilies, MultiParamTypeClasses, FlexibleInstances, DeriveFunctor, FunctionalDependencies #-}

-- | Basic 2 dimensional geometry functions.
module Geom2D where
import qualified Data.Vector.Generic.Mutable as M
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as V
import Control.Monad


infixl 6 ^+^, ^-^
infixl 7 *^, ^*, ^/
infixr 5 $*

data Point a = Point {
  pointX :: !a,
  pointY :: !a}
             deriving (Eq, Functor)

type DPoint = Point Double

instance Show a => Show (Point a) where
  show (Point x y) =
    "Point " ++ show x ++ " " ++ show y

-- | A transformation (x, y) -> (ax + by + c, dx + ey + d)
data Transform a = Transform {
  xformA :: !a,
  xformB :: !a,
  xformC :: !a,
  xformD :: !a,
  xformE :: !a,
  xformF :: !a }
                 deriving (Eq, Show, Functor)

data Line a = Line (Point a) (Point a)
            deriving (Show, Eq, Functor)
data Polygon a = Polygon [Point a]
               deriving (Show, Eq, Functor)

class AffineTransform a b | a -> b where
  transform :: Transform b -> a -> a

instance Num a => AffineTransform (Transform a) a where
  {-# INLINE transform #-}
  transform (Transform a' b' c' d' e' f') (Transform a b c d e f)  =
    Transform (a*a'+b'*d) (a'*b + b'*e) (a'*c+b'*f +c')
    (d'*a+e'*d) (d'*b+e'*e) (d'*c+e'*f+f')
    
instance Num a => AffineTransform (Point a) a where
  {-# INLINE transform #-}
  transform (Transform a b c d e f) (Point x y) =
    Point (a*x + b*y + c) (d*x + e*y + f)

instance Num a => AffineTransform (Polygon a) a where
  {-# INLINE transform #-}
  transform t (Polygon p) = Polygon $ map (transform t) p

newtype instance V.MVector s (Point a) = MV_Point (V.MVector s (a, a))
newtype instance V.Vector    (Point a) = V_Point  (V.Vector    (a, a))

instance V.Unbox a => V.Unbox (Point a)
instance V.Unbox a => M.MVector V.MVector (Point a) where
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicOverlaps #-}
  {-# INLINE basicUnsafeNew #-}
  {-# INLINE basicUnsafeReplicate #-}
  {-# INLINE basicUnsafeRead #-}
  {-# INLINE basicUnsafeWrite #-}
  {-# INLINE basicClear #-}
  {-# INLINE basicSet #-}
  {-# INLINE basicUnsafeCopy #-}
  {-# INLINE basicUnsafeGrow #-}
  basicLength (MV_Point v) = M.basicLength v
  basicUnsafeSlice i n (MV_Point v) = MV_Point $ M.basicUnsafeSlice i n v
  basicOverlaps (MV_Point v1) (MV_Point v2) = M.basicOverlaps v1 v2
  basicUnsafeNew n = MV_Point `liftM` M.basicUnsafeNew n
  basicUnsafeReplicate n (Point x y) = MV_Point `liftM` M.basicUnsafeReplicate n (x,y)
  basicUnsafeRead (MV_Point v) i = uncurry Point `liftM` M.basicUnsafeRead v i
  basicUnsafeWrite (MV_Point v) i (Point x y) = M.basicUnsafeWrite v i (x,y)
  basicClear (MV_Point v) = M.basicClear v
  basicSet (MV_Point v) (Point x y) = M.basicSet v (x,y)
  basicUnsafeCopy (MV_Point v1) (MV_Point v2) = M.basicUnsafeCopy v1 v2
  basicUnsafeGrow (MV_Point v) n = MV_Point `liftM` M.basicUnsafeGrow v n

instance V.Unbox a => G.Vector V.Vector (Point a) where
  {-# INLINE basicUnsafeFreeze #-}
  {-# INLINE basicUnsafeThaw #-}
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicUnsafeIndexM #-}
  {-# INLINE elemseq #-}
  basicUnsafeFreeze (MV_Point v) = V_Point `liftM` G.basicUnsafeFreeze v
  basicUnsafeThaw (V_Point v) = MV_Point `liftM` G.basicUnsafeThaw v
  basicLength (V_Point v) = G.basicLength v
  basicUnsafeSlice i n (V_Point v) = V_Point $ G.basicUnsafeSlice i n v
  basicUnsafeIndexM (V_Point v) i
                = uncurry Point `liftM` G.basicUnsafeIndexM v i
  basicUnsafeCopy (MV_Point mv) (V_Point v)
                = G.basicUnsafeCopy mv v
  elemseq _ (Point x y) z = G.elemseq (undefined :: V.Vector a) x
                       $ G.elemseq (undefined :: V.Vector a) y z

-- | Operator for applying a transformation.
($*) :: AffineTransform a b => Transform b -> a -> a
t $* p = transform t p
{-# INLINE ($*) #-}

-- | Calculate the inverse of a transformation.
inverse :: (Eq a, Num a, Fractional a) => Transform a -> Maybe (Transform a)
inverse (Transform a b c d e f) = case a*e - b*d of
  0 -> Nothing
  det -> Just $! Transform (a/det) (d/det) (-(a*c + d*f)/det) (b/det) (e/det)
         (-(b*c + e*f)/det)
{-# SPECIALIZE inverse :: Transform Double -> Maybe (Transform Double) #-}         

-- | Return the parameters (a, b, c) for the normalised equation
-- of the line: @a*x + b*y + c = 0@.
lineEquation :: Floating t => Line t -> (t, t, t)
lineEquation (Line (Point x1 y1) (Point x2 y2)) = (a, b, c)
  where a = a' / d
        b = b' / d
        c = -(y1*b' + x1*a') / d
        a' = y1 - y2
        b' = x2 - x1
        d = sqrt(a'*a' + b'*b')
{-# SPECIALIZE lineEquation :: Line Double -> (Double, Double, Double) #-}        

-- | Return the signed distance from a point to the line.  If the
-- distance is negative, the point lies to the right of the line
lineDistance :: Floating a => Line a -> Point a -> a
lineDistance l = \(Point x y) -> a*x + b*y + c
  where (a, b, c) = lineEquation l
{-# SPECIALIZE lineDistance :: Line Double -> DPoint -> Double #-}

-- | Return the point on the line closest to the given point.
closestPoint :: Fractional a => Line a -> Point a -> Point a
closestPoint (Line p1 p2) p3 = Point px py
  where
    (Point dx dy) = p2 ^-^ p1
    u = dy*pointY p3 + dx*pointX p3
    v = pointX p1*pointY p2 - pointX p2*pointY p1
    m = dx*dx + dy*dy
    px = (dx*u + dy*v) / m
    py = (dy*u - dx*v) / m
{-# specialize closestPoint :: Line Double -> Point Double -> Point Double #-}  

-- | Calculate the intersection of two lines.  If the determinant is
-- less than tolerance (parallel or coincident lines), return Nothing.
lineIntersect :: (Ord a, Floating a) => Line a -> Line a -> a -> Maybe (Point a)
lineIntersect (Line p1 p2) (Line p3 p4) eps
  | abs det <= eps = Nothing
  | otherwise = Just $ (a*^d2 ^-^ b*^d1) ^/ det
  where
    d1 = p1 ^-^ p2
    d2 = p3 ^-^ p4
    det = vectorCross d1 d2
    a = vectorCross p1 p2 
    b = vectorCross p3 p4
{-# SPECIALIZE lineIntersect :: Line Double -> Line Double -> Double -> Maybe DPoint #-}    

-- | The lenght of the vector.
vectorMag :: Floating a => Point a -> a
vectorMag (Point x y) = sqrt(x*x + y*y)
{-# INLINE vectorMag #-}

-- | The angle of the vector, in the range @(-'pi', 'pi']@.
vectorAngle :: RealFloat a => Point a -> a
vectorAngle (Point 0.0 0.0) = 0.0
vectorAngle (Point x y) = atan2 y x
{-# INLINE vectorAngle #-}

-- | The unitvector with the given angle.
dirVector :: Floating a => a -> Point a
dirVector angle = Point (cos angle) (sin angle)
{-# INLINE dirVector #-}

-- | The unit vector with the same direction.
normVector :: Floating a => Point a -> Point a
normVector p@(Point x y) = Point (x/l) (y/l)
  where l = vectorMag p
{-# INLINE normVector #-}        

-- | Scale vector by constant.
(*^) :: Num a => a -> Point a -> Point a
s *^ (Point x y) = Point (s*x) (s*y)
{-# INLINE (*^) #-}

-- | Scale vector by reciprocal of constant.
(^/) :: Fractional a => Point a -> a -> Point a
(Point x y) ^/ s = Point (x/s) (y/s)
{-# INLINE (^/) #-}

-- | Scale vector by constant, with the arguments swapped.
(^*) :: Num a => Point a -> a -> Point a
p ^* s = s *^ p
{-# INLINE (^*) #-}

-- | Add two vectors.
(^+^) :: Num a => Point a -> Point a -> Point a
(Point x1 y1) ^+^ (Point x2 y2) = Point (x1+x2) (y1+y2)
{-# INLINE (^+^) #-}

-- | Subtract two vectors.
(^-^) :: Num a => Point a -> Point a -> Point a
(Point x1 y1) ^-^ (Point x2 y2) = Point (x1-x2) (y1-y2)
{-# INLINE (^-^) #-}

-- | Dot product of two vectors.
(^.^) :: Num a => Point a -> Point a -> a
(Point x1 y1) ^.^ (Point x2 y2) = x1*x2 + y1*y2
{-# INLINE (^.^) #-}

-- | Cross product of two vectors.
vectorCross :: Num a => Point a -> Point a -> a
vectorCross (Point x1 y1) (Point x2 y2) = x1*y2 - y1*x2
{-# INLINE vectorCross #-}

-- | Distance between two vectors.
vectorDistance :: Floating a => Point a -> Point a -> a
vectorDistance p q = vectorMag (p^-^q)
{-# INLINE vectorDistance #-}

-- | Interpolate between two vectors.
interpolateVector :: (Num a) => Point a -> Point a -> a -> Point a
interpolateVector a b t = t*^b ^+^ (1-t)*^a
{-# INLINE interpolateVector #-}

-- | Create a transform that rotates by the angle of the given vector
-- and multiplies with the magnitude of the vector.
rotateScaleVec :: Num a => Point a -> Transform a
rotateScaleVec (Point x y) = Transform x (-y) 0 y x 0
{-# INLINE rotateScaleVec #-}

-- | reflect the vector over the X-axis.
flipVector :: (Num a) => Point a -> Point a
flipVector (Point x y) = Point x (-y)
{-# INLINE flipVector #-}

-- | Create a transform that rotates by the angle of the given vector
-- with the x-axis
rotateVec :: Floating a => Point a -> Transform a
rotateVec v = Transform x (-y) 0 y x 0
  where Point x y = normVector v
{-# INLINE rotateVec #-}

-- | Create a transform that rotates by the given angle (radians).
rotate :: Floating s => s -> Transform s
rotate a = Transform (cos a) (negate $ sin a) 0
           (sin a) (cos a) 0
{-# INLINE rotate #-}

-- | Rotate vector 90 degrees left.
rotate90L :: Floating s => Transform s
rotate90L = rotateVec (Point 0 1)
{-# INLINE rotate90L #-}

-- | Rotate vector 90 degrees right.
rotate90R :: Floating s => Transform s
rotate90R = rotateVec (Point 0 (-1))
{-# INLINE rotate90R #-}

-- | Create a transform that translates by the given vector.
translate :: Num a => Point a -> Transform a
translate (Point x y) = Transform 1 0 x 0 1 y
{-# INLINE translate #-}

-- | The identity transformation.
idTrans :: Num a => Transform a
idTrans = Transform 1 0 0 0 1 0
{-# INLINE idTrans #-}
