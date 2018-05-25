{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
module Geom2D.CubicBezier.StrokeIntern
       where
import Geom2D
import Geom2D.CubicBezier
import Data.Monoid
import qualified Data.Vector as V

data Pen a = PenEllipse (Transform a) (Transform  a) (Transform a)

-- | A circular pen with unit radius.
penCircle :: Floating a => Pen a
penCircle = PenEllipse idTrans rotate90L rotate90R

-- | Create a pen from a path.  For predictable results the path
-- should be convex.
pathToPen :: Path Closed a -> Pen a
pathToPen = undefined

noTranslate :: Num a => Transform a -> Transform a
noTranslate (Transform a b _ c d _) =
  Transform a b 0 c d 0

instance (Floating a, Eq a) => AffineTransform (Pen a) a where
  {-# SPECIALIZE transform :: Transform Double -> Pen Double -> Pen Double #-}
  transform t (PenEllipse trans _ _) =
    let t2@(Transform a b _ d e _) = transform t trans
    in case inverse $ noTranslate t2 of
      Nothing -> pathToPen $ undefined
        -- Path [
        -- (Point c f ^+^ p, JoinLine),
        -- (Point c f ^-^ p, JoinLine)]
        where
          p | a /= 0 && b /= 0 =
                sqrt(1 + a*a/(b*b)) *^ Point a d
            | d /= 0 && e /= 0 =
                sqrt(1 + d*d/(e*e)) *^ Point a d
            | a /= 0 = Point (a+d) 0
            | b /= 0 = Point 0 (b+e)
              -- singular point: create tiny pen instead of an error
            | otherwise = Point 1e-5 1e-5
      Just inv ->
        PenEllipse t2 (transform rotate90L inv) (transform rotate90R inv)

offsetPoint :: (Floating a) =>  a -> Point a -> Point a -> Point a
offsetPoint dist start tangent =
  start ^+^ (rotate90L $* dist *^ normVector tangent)

bezierOffsetPoint :: CubicBezier Double -> Double -> Double -> (DPoint, DPoint)
bezierOffsetPoint cb dist t = (offsetPoint dist p p', p')
  where (p, p') = evalBezierDeriv cb t

-- | Calculate an offset path from the bezier curve to within
-- tolerance.  If the distance is positive offset to the left,
-- otherwise to the right. A smaller tolerance may require more bezier
-- curves in the path to approximate the offset curve
bezierOffset :: CubicBezier Double -- ^ The curve
             -> Double      -- ^ Offset distance.
             -> Maybe Int   -- ^ maximum subcurves
             -> Double      -- ^ Tolerance.
             -> Bool        -- ^ Calculate the curve faster but with
                            -- more subcurves
             -> Path Open Double        -- ^ The offset curve
bezierOffset cb dist (Just m) tol faster =
  approximatePathMax m (bezierOffsetPoint cb dist) 15 tol 0 1 faster

bezierOffset cb dist Nothing tol faster =
  approximatePath (bezierOffsetPoint cb dist) 15 tol 0 1 faster

penOffset :: (Floating a) => Pen a -> Point a -> Point a
penOffset (PenEllipse trans leftInv _) dir =
  transform trans $ normVector $ leftInv $* dir

penOffsetFun :: Pen Double -> (Double -> (DPoint, DPoint)) -> Double
             -> (Point Double, Point Double)
penOffsetFun pen f t =
  (px ^+^ penOffset pen px', px')
  where
    (px, px') = f t

strokeLine :: Floating a => Pen a -> Line a -> Line a
strokeLine pen (Line a b) = Line (a ^+^ o) (b ^+^ o)
  where
    o = penOffset pen (b ^-^ a)
    
strokeCurve :: Pen a -> CubicBezier a -> [CubicBezier a]
strokeCurve pen cb =
  undefined

ellipticArc :: Transform Double -> Transform Double
            -> Point Double -> Point Double -> CubicBezier Double
ellipticArc trans leftInv from to =
  trans $* bezierArc
  (vectorAngle $ leftInv $* from)
  (vectorAngle $ leftInv $* to)

penJoin :: Pen Double -> Point Double
        -> Point Double -> Path Open Double
penJoin pen@(PenEllipse trans leftInv _) from to =
  undefined

