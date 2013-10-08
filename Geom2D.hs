{-# LANGUAGE BangPatterns #-}

-- | Basic 2 dimensional geometry functions.
module Geom2D where

infixl 6 ^+^, ^-^
infixl 7 *^, ^*, ^/
infixr 5 $*

data Point = Point {
  pointX :: {-# UNPACK #-} !Double,
  pointY :: {-# UNPACK #-} !Double}
           deriving Eq

instance Show Point where
  show (Point x y) =
    "Point " ++ show x ++ " " ++ show y

-- | A transformation (x, y) -> (ax + by + c, dx + ey + d)
data Transform = Transform {
  xformA :: {-# UNPACK #-} !Double,
  xformB :: {-# UNPACK #-} !Double,
  xformC :: {-# UNPACK #-} !Double,
  xformD :: {-# UNPACK #-} !Double,
  xformE :: {-# UNPACK #-} !Double,
  xformF :: {-# UNPACK #-} !Double }
               deriving Show

data Line = Line Point Point
data Polygon = Polygon [Point]

class AffineTransform a where
  transform :: Transform -> a -> a

instance AffineTransform Transform where
  transform (Transform a' b' c' d' e' f') (Transform a b c d e f)  =
    Transform (a*a'+b'*d) (a'*b + b'*e) (a'*c+b'*f +c')
    (d'*a+e'*d) (d'*b+e'*e) (d'*c+e'*f+f')
    
instance AffineTransform Point where
  transform (Transform a b c d e f) (Point x y) =
    Point (a*x + b*y + c) (d*x + e*y + f)

instance AffineTransform Polygon where
  transform t (Polygon p) = Polygon $ map (transform t) p

-- | Operator for applying a transformation.
($*) :: AffineTransform a => Transform -> a -> a
t $* p = transform t p

-- | Calculate the inverse of a transformation.
inverse :: Transform -> Maybe Transform
inverse (Transform a b c d e f) = case a*e - b*d of
  0 -> Nothing
  det -> Just $! Transform (a/det) (d/det) (-(a*c + d*f)/det) (b/det) (e/det)
         (-(b*c + e*f)/det)

-- | Return the parameters (a, b, c) for the normalised equation
-- of the line: @a*x + b*y + c = 0@.
lineEquation :: Line -> (Double, Double, Double)
lineEquation (Line (Point x1 y1) (Point x2 y2)) = (a, b, c)
  where a = a' / d
        b = b' / d
        c = -(y1*b' + x1*a') / d
        a' = y1 - y2
        b' = x2 - x1
        d = sqrt(a'*a' + b'*b')

-- | Return the the distance from a point to the line.
lineDistance :: Line -> Point -> Double
lineDistance l = \(Point x y) -> a*x + b*y + c
  where (a, b, c) = lineEquation l

-- | The lenght of the vector.
vectorMag :: Point -> Double
vectorMag (Point x y) = sqrt(x*x + y*y)
{-# INLINE vectorMag #-}

-- | The angle of the vector, in the range @(-'pi', 'pi']@.
vectorAngle :: Point -> Double
vectorAngle (Point 0.0 0.0) = 0.0
vectorAngle (Point x y) = atan2 y x
{-# INLINE vectorAngle #-}

-- | The unitvector with the given angle.
dirVector :: Double -> Point
dirVector angle = Point (cos angle) (sin angle)
{-# INLINE dirVector #-}

-- | The unit vector with the same direction.
normVector :: Point -> Point
normVector p@(Point x y) = Point (x/l) (y/l)
  where l = vectorMag p

-- | Scale vector by constant.
(*^) :: Double -> Point -> Point
s *^ (Point x y) = Point (s*x) (s*y)
{-# INLINE (*^) #-}

-- | Scale vector by reciprocal of constant.
(^/) :: Point -> Double -> Point
(Point x y) ^/ s = Point (x/s) (y/s)
{-# INLINE (^/) #-}

-- | Scale vector by constant, with the arguments swapped.
(^*) :: Point -> Double -> Point
p ^* s = s *^ p
{-# INLINE (^*) #-}

-- | Add two vectors.
(^+^) :: Point -> Point -> Point
(Point x1 y1) ^+^ (Point x2 y2) = Point (x1+x2) (y1+y2)
{-# INLINE (^+^) #-}

-- | Subtract two vectors.
(^-^) :: Point -> Point -> Point
(Point x1 y1) ^-^ (Point x2 y2) = Point (x1-x2) (y1-y2)
{-# INLINE (^-^) #-}

-- | Dot product of two vectors.
(^.^) :: Point -> Point -> Double
(Point x1 y1) ^.^ (Point x2 y2) = x1*x2 + y1*y2
{-# INLINE (^.^) #-}

-- | Cross product of two vectors.
vectorCross :: Point -> Point -> Double
vectorCross (Point x1 y1) (Point x2 y2) = x1*y2 - y1*x2
{-# INLINE vectorCross #-}

-- | Distance between two vectors.
vectorDistance :: Point -> Point -> Double
vectorDistance p q = vectorMag (p^-^q)
{-# INLINE vectorDistance #-}

-- | Interpolate between two vectors.
interpolateVector :: Point -> Point -> Double -> Point
interpolateVector a b t = t*^b ^+^ (1-t)*^a
{-# INLINE interpolateVector #-}

-- | Create a transform that rotates by the angle of the given vector
-- with the x-axis
rotateVec :: Point -> Transform
rotateVec v = Transform x (-y) 0 y x 0
  where Point x y = normVector v

-- | Create a transform that rotates by the given angle (radians).
rotate :: Double -> Transform
rotate a = Transform (cos a) (negate $ sin a) 0
           (sin a) (cos a) 0

-- | Rotate vector 90 degrees left.
rotate90L :: Transform
rotate90L = rotateVec (Point 0 1)

-- | Rotate vector 90 degrees right.
rotate90R :: Transform
rotate90R = rotateVec (Point 0 (-1))

-- | Create a transform that translates by the given vector.
translate :: Point -> Transform
translate (Point x y) = Transform 1 0 x 0 1 y

