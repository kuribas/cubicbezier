module Geom2D where

data Point = Point {
  pointX :: Double,
  pointY :: Double}

data Transform = Transform {
  xformX1 :: Double,
  xformY1 :: Double,
  xformX2 :: Double,
  xformY2 :: Double,
  xformX3 :: Double,
  xformY3 :: Double }

vectorLength (Point x y) = sqrt(x*x + y*y)

unitVector p@(Point x y) = Point (x/l) (y/l)
  where l = vectorLength p

rotateVector90Left p@(Point x y) = Point (-y, x)
rotateVector90Right p@(Point x y) = Point (y, -x)

scaleVector s (Point x y) = Point (s*x) (s*y)
addVector (Point x1 y1) (Point x2 y2) = Point (x1+x2) (y1+y2)
subVector (Point x1 y1) (Point x2 y2) = Point (x1-x2) (y1-y2)
dotProduct (Point x1 y1) (Point x2 y2) = Point (x1*x2) (y1*y2)
vectorDistance (Point x1 y1) (Point x2 y2) = sqrt ((x1-x2)^2 + (y1-y2)^2)
