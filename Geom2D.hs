-- | Basic 2 dimensional geometry functions.

module Geom2D where

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

-- | the lenght of the vector
vectorLength :: Point -> Double
vectorLength (Point x y) = sqrt(x*x + y*y)

-- | the unit vector with the same direction
unitVector :: Point -> Point
unitVector p@(Point x y) = Point (x/l) (y/l)
  where l = vectorLength p

-- | rotate vector 90 degrees left
rotateVector90Left :: Point -> Point
rotateVector90Left (Point x y) = Point (-y) x

-- | rotate vector 90 degrees right
rotateVector90Right :: Point -> Point
rotateVector90Right (Point x y) = Point y (-x)

-- | scale vector by constant
scaleVector :: Double -> Point -> Point
scaleVector s (Point x y) = Point (s*x) (s*y)

-- | add two vectors
addVector :: Point -> Point -> Point
addVector (Point x1 y1) (Point x2 y2) = Point (x1+x2) (y1+y2)

-- | subtract two vectors
subVector :: Point -> Point -> Point
subVector (Point x1 y1) (Point x2 y2) = Point (x1-x2) (y1-y2)

-- | dot product of two vectors
dotProduct :: Point -> Point -> Point
dotProduct (Point x1 y1) (Point x2 y2) = Point (x1*x2) (y1*y2)

-- | distance between two vectors
vectorDistance :: Point -> Point -> Double
vectorDistance (Point x1 y1) (Point x2 y2) = sqrt ((x1-x2)^2 + (y1-y2)^2)
