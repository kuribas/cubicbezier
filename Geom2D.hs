-- | Basic 2 dimensional geometry functions.
module Geom2D where

infixl 6 ^+^, ^-^
infixl 7 *^, ^*

data Point = Point {
  pointX :: Double,
  pointY :: Double}
           deriving Show

data Transform = Transform {
  xformX1 :: Double,
  xformY1 :: Double,
  xformX2 :: Double,
  xformY2 :: Double,
  xformX3 :: Double,
  xformY3 :: Double }
               deriving Show

-- | The lenght of the vector.
vectorMag :: Point -> Double
vectorMag (Point x y) = sqrt(x*x + y*y)

-- | The angle of the vector, in the range @(-'pi', 'pi']@.
vectorAngle (Point 0.0 0.0) = 0.0
vectorAngle (Point x y) = atan2 y x

-- | The unitvector with the given angle.
dirVector angle = Point (cos angle) (sin angle)

-- | The unit vector with the same direction.
normVector :: Point -> Point
normVector p@(Point x y) = Point (x/l) (y/l)
  where l = vectorMag p

-- | Rotate vector 90 degrees left.
rotateVector90Left :: Point -> Point
rotateVector90Left (Point x y) = Point (-y) x

-- | Rotate vector 90 degrees right.
rotateVector90Right :: Point -> Point
rotateVector90Right (Point x y) = Point y (-x)

-- | Scale vector by constant.
(*^) :: Double -> Point -> Point
s *^ (Point x y) = Point (s*x) (s*y)

-- | Scale vector by constant, with the arguments swapped.
(^*) :: Point -> Double -> Point
p ^* s = s *^ p

-- | Add two vectors.
(^+^) :: Point -> Point -> Point
(Point x1 y1) ^+^ (Point x2 y2) = Point (x1+x2) (y1+y2)

-- | Subtract two vectors.
(^-^) :: Point -> Point -> Point
(Point x1 y1) ^-^ (Point x2 y2) = Point (x1-x2) (y1-y2)

-- | Dot product of two vectors.
dotProduct :: Point -> Point -> Double
dotProduct (Point x1 y1) (Point x2 y2) = x1*x2 + y1*y2

-- | Distance between two vectors.
vectorDistance :: Point -> Point -> Double
vectorDistance p q = vectorMag (p^-^q)

-- | Interpolate between two vectors.
interpolateVector a b t = t*^b ^+^ (1-t)*^a
