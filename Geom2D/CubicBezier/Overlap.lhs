> {-# LANGUAGE MultiWayIf, PatternGuards, TemplateHaskell, BangPatterns, CPP #-}

Removing overlap from bezier paths in haskell
=============================================

This literate programming document serves as both the description 
and code for the algorithm for removing overlap and
performing set operations on bezier paths.  It can be used both for
understanding the code, and for porting the used algorithm to other 
implementations.

**Note on porting**: Porting should be fairly straightforward. 
However some care must be taken with regards to lazyness.  Often
many variables inside the `where` statement aren't
evaluated in all guards, so it's important to evaluate only those
which appear in the guards.  This library uses copying instead of 
mutation, but for the most part mutation can be used as well.
It uses the `lens` library for modifying state.  The following translation
could be made to a mutable language: 

  * reading state:
    - `view field struct`: `struct.field`
    - `get`: `state` (implicit state, usually sweepstate)
    - `use field`: `state.field`

  * writing state:
    - `set field value struct`: `struct.field = value`
    - `field .= value`: `state.field = value` (typically sweepstate)

  * modifying state:
    - `over field fun struct`: `struct.field = fun (struct.field)`
    - `modify fun`: `state = fun state`
    - `modifying field fun`: `state.field = fun (state.field)`

Let's begin with declaring the module and library imports:

> module Geom2D.CubicBezier.Overlap
>        (boolPathOp, union, intersection, difference,
>         exclusion, FillRule (..))
>        where
> import Prelude
> import Geom2D
> import Geom2D.CubicBezier.Basic
> import Geom2D.CubicBezier.Intersection
> import Geom2D.CubicBezier.Numeric
> import Math.BernsteinPoly
> import Data.Foldable (traverse_)
> import Data.Functor ((<$>))
> import Data.List (sortBy, sort, intercalate, intersperse)
> import Control.Monad.State.Strict
> import Lens.Micro
> import Lens.Micro.TH
> import Lens.Micro.Mtl
> import qualified Data.Map.Strict as M
> import qualified Data.Set as S 
> import Text.Printf
> import Data.Ratio
> import Data.Tuple
> import Data.IORef
> import Data.Maybe (isJust, isNothing, mapMaybe)

#ifdef DEBUG
> import System.IO.Unsafe (unsafePerformIO)
> import System.IO
> import Debug.Trace
#endif

The basic idea is to keep curves where one side is inside the filled
region, and the other side is outside, and discard the rest. 
Since that could be true only of a part of the curve, I also need to
split each curve when it intersects
another curve.  How to know which side is the inside, and which
side the outside?  There are two methods which are use the most: the
[*even-odd rule*](https://en.wikipedia.org/wiki/Even%E2%80%93odd_rule)
and the [*nonzero rule*](https://en.wikipedia.org/wiki/Nonzero-rule).
Instead of hardwiring it, I use higher-order functions to determine
when a turnratio is inside the region to be filled, and how the
turnratio changes with each curve.

Checking each pair of curves for intersections would work, but is
rather inefficient.  I only need to check for overlap when two curves
are adjacent.  Fortunately there exist a good method from
*computational geometry*, called the *sweep line algorithm*.  The 
idea is to sweep a vertical line over the input, starting from
leftmost point to the right (of course the opposite direction is also
possible), and to update the input dynamically.  I keep track of each
curve that intersects the sweepline by using a balanced tree of
curves.  When adding a new curve, it's only necessary to check for
intersections with the curve above and below.  Since searching on the
tree takes only `O(log n)` operations, this will save a lot of
computation.

The input is processed in horizontal order, and after splitting curves
the order must be preserved, so an ordered structure is needed.  The
standard map library from `Data.Map` is ideal, and has all the
operations needed.  This structure is called the *X-structure*,
since the data is ordered by X coordinate.:

> type XStruct = M.Map PointEvent [Curve]

I use `PointEvent` instead of just `Point`.  This way I can have a `Ord`
instance for the map, which must match the horizontal ordering.  A
newtype is ideal, since it has no extra cost, and allows me to define
a Ord instance for defining the relative order.  The value from the
map is a list, since there can be many curves starting from the same
point.

> newtype PointEvent = PointEvent DPoint
>
> instance Show PointEvent where
>    show (PointEvent (Point px py)) = printf "(%.5g, %.5g)" px py

When the x-coordinates are equal, use the y-coordinate to determine
the order.

> instance Eq PointEvent where
>   (PointEvent (Point x1 y1)) == (PointEvent (Point x2 y2)) =
>     (x1, y1) == (x2, y2)
>
> instance Ord PointEvent where
>   compare (PointEvent (Point x1 y1)) (PointEvent (Point x2 y2)) =
>     compare (x1, y2) (x2, y1)

All curves are kept left to right, so I need to remember the
curve direction for the output:

The curves intersecting the sweepline are kept in another balanced
Tree, called the *Y-structure*.  *These curves are not allowed to
overlap*, except in the endpoints, and will be ordered vertically.
The `Curve` datatype defines the ordering of the curves,
and adds additional information.  The `turnRatio` field is the
turnRatio of the area to the left for a left to right curve, and to
the right for a right to left curve.  The `changeTurn` function
determines how the turnRatio will change from up to down.  This
together with a test for the *insideness* of a certain turnratio,
allows for more flexibility.  Using this, it is possible to generalize
this algorithm to boolean operations!

The curveRank parameter is used to memoize the order in the Ystruct.
This will avoid costly comparisons for tree rebalancing etc...

The FillRule datatype is used for the exported API:

> data FillRule = EvenOdd | NonZero

> data Curve = Curve {
>   _bezier :: !(CubicBezier Double),
>   _turnRatio :: !(Int, Int),
>   _changeTurn :: !((Int, Int) -> (Int, Int)),
>   _curveRank :: Maybe (Ratio Integer)}
>
> trOne :: (Int, Int)
> trOne = (0,0)
>
> between :: Ratio Integer -> Ratio Integer -> Ratio Integer
> between a b = a +(b-a)/2
> 
> makeLenses ''Curve
>
> instance Show Curve where
>   show (Curve (CubicBezier (Point p0x p0y) (Point p1x p1y)
>               (Point p2x p2y) (Point p3x p3y)) (t1, t2) _ o) =
>     printf "Curve (%.5g, %.5g) (%.5g, %.5g) (%.5g, %.5g) (%.5g, %.5g) (%i,%i) %s" 
>     p0x p0y p1x p1y p2x p2y p3x p3y t1 t2 (show o)
> 
> type YStruct = S.Set Curve

The total state for the algorithm consists of the X-structure, the
Y-structure, and the output found so far.  I use a trick to make
access to curves above and below the current pointevent more
convenient.  I use two sets to represent a focus point into the
Y-structure, where the left set are the elements less than the
pointEvent (above), and the right set the elements greater (below):

> data SweepState = SweepState {
>   _output :: !(M.Map PointEvent [CubicBezier Double]),
>   _yStruct :: !YStruct,
>   _focusPoint :: DPoint,
>   _xStruct :: !XStruct}
>                   
> makeLenses ''SweepState

> singularC :: Point Double -> Curve
> singularC p = Curve (CubicBezier p p p p) trOne id Nothing
>

Some functions for debugging:

> showCurve (CubicBezier p0 p1 p2 p3) =
>   showPt p0 ++ showPt p1 ++ showPt p2 ++ showPt p3
>
> showPt :: DPoint -> String
> showPt (Point x y) = "(" ++ show x ++ "," ++ show y ++ ")"
>
#ifdef DEBUG
> type SweepStateM = StateT SweepState IO
>
> traceMessage :: String -> SweepStateM ()
> traceMessage msg = liftIO $ hPutStrLn stderr msg
>
> assert :: Bool -> String -> SweepStateM ()
> assert p msg = unless p $ liftIO $ hPutStrLn stderr $ "ASSERT " ++ msg
>
> assertTrace p msg e
>   | p = e
>   | otherwise = trace ("ASSERT " ++ msg) e
>
#else
> -- | output a trace of the algorithm when compiled with @-fdebug@.
> type SweepStateM  = State SweepState
>
> assertTrace _ _ e  = e
> 
> traceMessage _ = return ()
>
> assert :: Bool -> String -> SweepStateM ()
> assert _ _ = return ()
#endif

> activate, deactivate :: [Curve] -> SweepStateM ()
> activate cs = 
>   traverse_ (traceMessage . ("ACTIVATE " ++) . showCurve . view bezier) cs
> 
> deactivate cs = 
>   traverse_ (traceMessage . ("DEACTIVATE " ++) . showCurve . view bezier) cs


This handy helper function will pass the first curve above to the
given function, and if it doesn't return `Nothing`, remove it from the
state.  It does nothing when there is no curve above.

> withAbove :: (Curve -> Maybe a) -> SweepStateM (Maybe a)
> withAbove f = do
>   p <- use focusPoint
>   yStr <- use yStruct
>   let i = yStructIndex (singularC p) yStr
>   if i < 0
>     then return Nothing
>     else let c = S.elemAt i yStr
>          in case f c of
>               Nothing ->
>                 return Nothing
>               Just x -> do
>                 yStructDel i
>                 return $ Just x

The same with the curve below.

> withBelow :: (Curve -> Maybe a) -> SweepStateM (Maybe a)
> withBelow f = do
>   p <- use focusPoint
>   yStr <- use yStruct
>   let i = yStructIndex (singularC p) yStr
>       s = S.size yStr
>   if i >= s-1
>     then return Nothing
>     else let c = S.elemAt (i+1) yStr
>          in case f c of
>               Nothing ->
>                 return Nothing
>               Just x -> do
>                 yStructDel (i+1)
>                 return $ Just x

`splitYStruct` changes the focus and returns and removes any curves which end in
the current pointEvent:

> splitYStruct :: DPoint -> SweepStateM [Curve]
> splitYStruct p = do
>   yStr <- use yStruct
>   focusPoint .= p
>   traceMessage $ "CHANGEFOCUS " ++ showPt p
>   traceMessage "MSG Remove curves ending at pointevent from Y structure" 
>   let lStr = fst $ S.split (singularC p) yStr
>       rightCurves = takeWhile (\c -> cubicC3 (_bezier c) == p) $
>                     S.toDescList lStr
>       nR = length rightCurves
>       i = S.size lStr - nR
>   replicateM_ nR (yStructDel i)
>   return rightCurves
>

=== Some functions on the Sweep state:

Adding and removing curves from the X structure.

> insertX :: PointEvent -> [Curve] -> SweepStateM ()
> insertX p c =
>   modify $ over xStruct $ M.insertWith (++) p c
>
> xStructAdd :: Curve -> SweepStateM ()
> xStructAdd c = do
>   traceMessage $ "XSTRUCTADD " ++ showCurve (view bezier c)
>   insertX (PointEvent $ cubicC0 $
>            view bezier c) [c]
>
> xStructRemove :: SweepStateM (PointEvent, [Curve])
> xStructRemove = do
>   str <- use xStruct
>   return str
>   (p, c) <- zoom xStruct $ state M.deleteFindMin
>   traverse_ (traceMessage . ("XSTRUCTREM " ++) .
>              showCurve . view bezier) c
>   str <- use xStruct
>   return str
>   return (p, c)
>
> yStructIndex :: Curve -> YStruct -> Int
> yStructIndex c str = S.size (fst $ S.split c str) - 1
>
> yStructDel :: Int -> SweepStateM ()
> yStructDel i = do
>   str <- use yStruct
>   traceMessage $ "YSTRUCTREM " ++ showCurve (view bezier $ S.elemAt i str)
>   yStruct .= S.deleteAt i str
>

Insert the curve into the Y structure.  First lookup the position of
the curve, then calculate the rank of the curve, using the surrounding
elements.  Insert the curve using the rank.  This will avoid repeating
expensive operations.

> 
> yStructAdd :: Curve -> SweepStateM ()
> yStructAdd c = do
>   str <- use yStruct
>   traceMessage $ "YSTRUCTADD " ++ showCurve (view bezier c)
>   traceMessage $ "YSTRUCT: " ++
>       (concat $ intersperse "\n  " $ map (showCurve . view bezier) $ S.toList str)
>   assert (isNothing (_curveRank c))
>     "CURVE ALREADY HAS A RANK IN THE YSTRUCT" 
>   assert (yStructConsistent c str)
>    ("Y STRUCT NOT CONSISTENT WITH CURVE:" ++
>    show (map (compare c) (S.toAscList str)) )
>   let i = yStructIndex c str
>       s = S.size str
>       newC 
>         | s == 0 = set curveRank (Just 0) c
>         | i >= s-1 = set curveRank ((+1) <$> _curveRank (S.elemAt (s-1) str)) c
>         | i < 0 = set curveRank (subtract 1 <$> _curveRank (S.elemAt 0 str)) c
>         | otherwise = set curveRank (liftM2 between
>                                       (_curveRank $ S.elemAt i str)
>                                       (_curveRank $ S.elemAt (i+1) str)) c
>   assert (not $ S.member c str)
>     ("CURVE ALREADY IN YSTRUCT: " ++ showCurve (view bezier c))
>   assert (S.size str < S.size (S.insert newC str)) $
>     "CURVE NOT ADDED TO YSTRUCT" ++ show newC 
>   yStruct .= S.insert newC str
> 
> yStructConsistent :: Curve -> YStruct -> Bool
> yStructConsistent c str =
>   all (> c) $ dropWhile (< c) $
>   S.toAscList str
>
> yStructOverlap :: Curve -> YStruct -> [String]
> yStructOverlap c str =
>   mapMaybe checkOverlap $ S.toAscList str
>   where checkOverlap c2 =
>           case splitMaybe c c2 1e-5 of
>             (Nothing, Nothing) -> Nothing
>             _ -> Just (show c2 ++ "\n")

To compare curves vertically, take the the curve which starts the
rightmost, and see if it falls below or above the curve.  If the first
control points are coincident, test the last control points instead,
or the midpoint.  This works because if the first point is coincident
the curves shouldn't intersect except in the endpoints (see #splitAndOrder).
To lookup a single point, I use a singular bezier curve.

> instance Eq Curve where
>    Curve _ _ _ (Just o1) == Curve _ _ _ (Just o2) = o1 == o2
>    Curve c1 t1 ct1 _ == Curve c2 t2 ct2 _ =
>     c1 == c2 && t1 == t2 && ct1 t1 == ct2 t2
>     
> instance Ord Curve where
>   compare (Curve _ _ _ (Just o1)) (Curve _ _ _ (Just o2)) = compare o1 o2
>   compare (Curve c1@(CubicBezier p0 p1 p2 p3) tr1 _ _)
>     (Curve c2@(CubicBezier q0 q1 q2 q3) tr2 _ _)
>     | p0 == q0 = if
>         | p3 == q3 ->
>             -- compare the midpoint
>             case compVert (evalBezier c1 0.5) c2 of
>              LT -> LT
>              GT -> GT
>              EQ ->
>                -- otherwise arbitrary
>                compare (tr1, PointEvent p1, PointEvent p2)
>                (tr2, PointEvent q1, PointEvent q2)
>         | pointX p3 < pointX q3 ->
>             case compVert p3 c2 of
>             LT -> LT
>             EQ -> LT
>             GT -> GT
>         | otherwise ->
>             case compVert q3 c1 of
>              LT -> GT
>              EQ -> GT
>              GT -> LT
>     | pointX p0 < pointX q0 =
>       case compVert q0 c1 of
>        LT -> GT
>        EQ -> LT
>        GT -> LT
>     | otherwise =
>       case compVert p0 c2 of
>       LT -> LT
>       EQ -> GT
>       GT -> GT

Compare a point with a curve.  See if it falls below or above the hull
first.  Otherwise find the point on the curve with the same
X-coordinate by iterating.

> compVert :: DPoint -> CubicBezier Double -> Ordering
> compVert p c
>   | p == cubicC0 c ||
>     p == cubicC3 c = EQ
>   | compH /= EQ = compH
>   | otherwise = comparePointCurve p c
>     where
>       compH = compareHull p c

=== Test if the point is above or below the curve {#comparePC}

> comparePointCurve :: Point Double -> CubicBezier Double -> Ordering
> comparePointCurve (Point x1 y1) c1@(CubicBezier p0 p1 p2 p3)
>   | pointX p0 == x1 &&
>     pointX p0 == pointX p1 &&
>     pointX p0 == pointX p2 &&
>     pointX p0 == pointX p3 =
>     compare (pointY p0) y1
>   | otherwise = compare y2 y1
>   where
>     t = findX c1 x1 (maximum (map (abs.pointX) [p0, p1, p2, p3])*1e-14)
>     y2 = pointY $ evalBezier c1 t

=== Comparing against the hull {#hull}

Compare a point against the convex hull of the bezier.  `EQ` means the
point is inside the hull, `LT` below and `GT` above.  I am currently
only testing against the control points, some testing needs to be done
to see what is faster.

> belowLine :: DPoint -> DPoint -> DPoint -> Bool
> belowLine (Point px py) (Point lx ly) (Point rx ry)
>   | lx == rx = True
>   | (px >= lx && px <= rx) ||
>     (px <= lx && px >= rx) = py < midY
>   | otherwise = True
>   where midY = ly + (ry-ly) * (rx-lx) / (px-lx)
> 
> aboveLine :: DPoint -> DPoint -> DPoint -> Bool
> aboveLine (Point px py) (Point lx ly) (Point rx ry)
>   | lx == rx = True
>   | (px >= lx && px <= rx) ||
>     (px <= lx && px >= rx) = py > midY
>   | otherwise = True
>   where midY = ly + (ry-ly) * (rx-lx) / (px-lx)
> 
> compareHull :: DPoint -> CubicBezier Double -> Ordering
> compareHull p (CubicBezier c0 c1 c2 c3)
>   | pointY p > pointY c0 &&
>     pointY p > pointY c1 &&
>     pointY p > pointY c2 &&
>     pointY p > pointY c3 = LT
>   | pointY p < pointY c0 &&
>     pointY p < pointY c1 &&
>     pointY p < pointY c2 &&
>     pointY p < pointY c3 = GT
>   | otherwise = EQ

Preprocessing
-------------

Since the algorithm assumes curves are increasing in the horizontal
direction they have to be preprocessed first.  I split each curve
where the tangent is vertical.  If the resulting subsegment is too
small however, I just adjust the control point to make the curve
vertical at the endpoint.

I also do snaprounding to prevent points closer than the tolerance.

> makeXStruct :: ((Int, Int) -> (Int, Int)) -> ((Int, Int) -> (Int, Int)) -> Double -> [CubicBezier Double] -> XStruct
> makeXStruct chTr chTrBack tol =
>   M.fromListWith (++) .
>   concatMap (toCurve . snapRoundBezier tol) .
>   concatMap (splitVert tol)
>   where toCurve c@(CubicBezier p0 _ _ p3) =
>           case compare (pointX p0) (pointX p3) of
>            LT -> [(PointEvent p0, [Curve c trOne chTr Nothing])]
>            GT -> [(PointEvent p3, [Curve (reorient c) trOne chTrBack Nothing]),
>                   (PointEvent p0, [])]
>            -- vertical curve
>            EQ | pointY p0 > pointY p3 ->
>                 [(PointEvent p0, [Curve c trOne chTr Nothing])]
>               | otherwise ->
>                 [(PointEvent p3, [Curve (reorient c) trOne chTrBack Nothing]),
>                  (PointEvent p0, [])]
>
> splitVert :: Double -> CubicBezier Double -> [CubicBezier Double]
> splitVert tol curve@(CubicBezier c0 c1 c2 c3) = 
>   uncurry splitBezierN $
>   adjustLast $
>   adjustFirst (curve, vert)
>   where vert
>           | pointX c0 == pointX c1 &&
>             pointX c0 == pointX c2 &&
>             pointX c0 == pointX c3 = []
>           | otherwise = 
>               sort $ bezierVert curve
>         -- adjust control points to avoid small curve fragments
>         -- near the endpoints
>         adjustFirst (c@(CubicBezier p0 p1 p2 p3), t:ts)
>           | vectorDistance p0 (evalBezier c t) < tol =
>               (CubicBezier p0 (Point (pointX p0) (pointY p1)) p2 p3,
>                ts)
>         adjustFirst x = x
>         adjustLast (c@(CubicBezier p0 p1 p2 p3), ts@(_:_))
>           | vectorDistance p3 (evalBezier c $ last ts) < tol =
>               (CubicBezier p0 p1 (Point (pointX p3) (pointY p2)) p3,
>                init ts)
>         adjustLast x = x

main loop
---------

For the main loop, I remove the leftmost point from the
X-structure, and do the following steps:

  1. Split any curves which come near the current pointEvent.

  2. Send all curves to the left of the sweepline to the output, after
  filtering them based on the turning number.

  3. For each curve starting at the point, split if it intersects with
the curve above or the curve below.  Sort resulting curves vertically.
If there are no curves starting from point, test the curves above and
below instead.  Adjust the turnRatios for each curve.

  4. Insert the points in the Y structure.

  5. Loop until the X-structure is empty

> loopEvents :: ((Int, Int) -> Bool) -> Double -> SweepStateM ()
> loopEvents isInside tol = do
>   xStr <- use xStruct
>   unless (M.null xStr) $ do
>       (PointEvent p, curves) <- xStructRemove
>       activate curves
>       
>       ending <- splitYStruct p
>       activate ending
>
>       -- split near curves
>       traceMessage "MSG Split curves near the focuspoint." 
>       (ending2, rightSubCurves) <- splitNearPoints p tol
>
>       activate ending2
>       activate rightSubCurves
>       traceMessage "MSG Output curves"
>        
>       -- output curves to the left of the sweepline.
>       deactivate (ending ++ ending2)
>       filterOutput (ending ++ ending2) isInside 
>       let allCurves = rightSubCurves ++ curves
>       if null allCurves
> 
>          -- split surrounding curves
>         then do
>              traceMessage "MSG Split curves around pointevent."
>              splitSurround tol
>         else do
> 
>         -- sort curves
>         traceMessage "MSG Sort curves."
>         
>         sorted <- splitAndOrder tol allCurves
>         deactivate allCurves
>         activate sorted
> 
>         -- split curve above
>         traceMessage "MSG Split curve above sorted curves."
>         deactivate sorted
>         curves2 <- splitAbove sorted tol
>         activate curves2
> 
>         -- add curves to Y structure
>         traceMessage "MSG Add curves to Y structure."
>         deactivate curves2
>         addMidCurves curves2 tol
>
>       loopEvents isInside tol


Send curves to output
---------------------

> outputPaths :: M.Map PointEvent [CubicBezier Double] -> [ClosedPath Double]
> outputPaths m
>   | M.null m = []
>   | otherwise = outputNext m
>   where
>     lookupDelete p m =
>       case M.lookup (PointEvent p) m of
>        Nothing -> Nothing
>        Just (x:xs) -> Just (x, m')
>          where m' | null xs = M.delete (PointEvent p) m
>                   | otherwise = M.insert (PointEvent p) xs m
>        _ -> error "outputPaths: empty list inside map."
>     outputNext !m
>       | M.null m = []
>       | otherwise = 
>         let ((PointEvent p0, c0:cs), m0) =
>               M.deleteFindMin m
>             m0' | null cs = m0
>                 | otherwise = M.insert (PointEvent p0) cs m0
>         in go m0' c0 [] p0
>     go !m !next !prev !start
>       | p == start =
>           curvesToPath (reverse $ next:prev):
>           outputNext m
>       | otherwise =
>         case lookupDelete p m of
>          Nothing -> outputNext m
>          Just (x, m') -> go m' x (next:prev) start
>       where p = cubicC3 next
>
> curvesToPath :: [CubicBezier Double] -> ClosedPath Double
> curvesToPath =
>   ClosedPath .
>   map (\(CubicBezier p0 p1 p2 _) ->
>         (p0, JoinCurve p1 p2))

Filter and output the given curves.  The `isInside` function
determines the *insideness* of a give turnratio.  For example for the
nonzero-rule, this would be `(> 0)`.  This inserts the curve into the
output map.

> filterOutput :: [Curve] -> ((Int, Int) -> Bool) -> SweepStateM ()
> filterOutput curves isInside =
>   mapM_ (outputCurve isInside)  curves
>
> outputCurve :: ((Int, Int) -> Bool) -> Curve -> SweepStateM ()
> outputCurve isInside (Curve c tr op _)
>   | isInside (op tr) /= isInside tr =
>       let c' | isInside tr = reorient c
>              | otherwise = c
>       in do traceMessage $ "OUTPUT " ++ showCurve c
>             modifying output (M.insertWith (++) (PointEvent $ cubicC0 c') [c'])
>   | otherwise =
>       traceMessage $ "DISCARD " ++ showCurve c

Test for intersections and split: (#splitAndOrder)
--------------------------------------------------

Since the curves going out of the current pointEvent in the X-structure are
unordered, they need to be ordered first.  First they are ordered by
first derivative.  Since it's easier to compare two curves when they
don't overlap, remove overlap, and then sort again by comparing the
whole curve.

To do this, I implemented a monadic insertion sort.  First the curves are split
in the statemonad, then they are compared.

> splitAndOrder :: Double -> [Curve] -> SweepStateM [Curve]
> splitAndOrder tol curves =
>   sortSplit tol $
>   sortBy compDeriv curves
>
> compDeriv :: Curve -> Curve -> Ordering
> compDeriv (Curve (CubicBezier p0 p1 _ _) _ _ _)
>   (Curve (CubicBezier q0 q1 _ _) _ _ _) =
>   compare (slope (q1^-^ q0)) (slope (p1^-^p0)) 
>
> slope (Point 0 0) = 0
> slope (Point x y) = y/x


Insertion sort, by splitting and comparing.  This should be efficient
enough, since ordering by derivative should mostly order the curves.

> sortSplit :: Double -> [Curve] -> SweepStateM [Curve]
> sortSplit _ [] = return []
> sortSplit tol (x:xs) =
>   insertM x tol =<<
>   sortSplit tol xs
>
> insertM :: Curve -> Double -> [Curve] -> SweepStateM [Curve]
> insertM x _ [] = return [x]
> insertM x tol (y:ys) =
>   case curveOverlap x y tol of
>    Just (c1, c2) -> do
>      traverse_ xStructAdd c2
>      insertM c1 tol ys
>    Nothing -> do
>      (x', y') <- splitM x y tol
>      if x' < y'
>        then return (x':y':ys)
>        else (y':) <$> insertM x' tol ys
>
> splitM :: Curve -> Curve -> Double -> SweepStateM (Curve, Curve)
> splitM x y tol =
>   case splitMaybe x y tol of
>   (Just (a, b), Just (c, d)) -> do
>     xStructAdd b
>     xStructAdd d
>     return (a, c)
>   (Nothing, Just (c, d)) -> do
>     xStructAdd d
>     return (x, c)
>   (Just (a, b), Nothing) -> do
>     xStructAdd b
>     return (a, y)
>   (Nothing, Nothing) ->
>     return (x, y)

Handle intersections of the first curve at point and the curve
above. Return the curves with updated turnratios.  Some care is needed
when one of the curves is intersected at the endpoints, in order not
to create singular curves.

> updateTurnRatio :: Curve -> Curve -> Curve
> updateTurnRatio (Curve _ tr chTr _) =
>   set turnRatio (chTr tr)
>
> propagateTurnRatio :: Curve -> [Curve] -> [Curve]
> propagateTurnRatio cAbove l =
>   tail $ scanl updateTurnRatio cAbove l
>
> splitAbove :: [Curve] -> Double -> SweepStateM [Curve]
> splitAbove [] _ = return []
> splitAbove (c:cs) tol = do
>   yStr <- use yStruct
>   p <- use focusPoint
>   let i = yStructIndex (singularC p) yStr
>   if i < 0
>     then let c' = set turnRatio trOne c
>                 in return $ c':propagateTurnRatio c' cs
>     else 
>       let cAbove = S.elemAt i yStr
>       in case splitMaybe c cAbove tol of
>         (Nothing, Nothing) ->
>           return $ propagateTurnRatio cAbove $ c:cs
>         (Just (c1, c2), Nothing)
>           | cubicC3 (_bezier c1) == cubicC0 (_bezier cAbove)
>             -> do
>               xStructAdd cAbove; xStructAdd c2
>               yStructDel i
>               return $ propagateTurnRatio cAbove $ c1:cs
>           | otherwise -> do
>               xStructAdd c2
>               return $ propagateTurnRatio cAbove $ c1:cs
>         (Nothing, Just (c3, c4)) -> do
>             assert (cubicC3 (_bezier c3) /= cubicC0 (_bezier c)) $
>               "curve intersecting pointevent: cAbove " ++ show cAbove
>             xStructAdd c4
>             yStructDel i; yStructAdd c3
>             return $ propagateTurnRatio cAbove $ c:cs
>         (Just (c1, c2), Just (c3, c4)) -> do
>           xStructAdd c2; xStructAdd c4
>           yStructDel i; yStructAdd c3
>           return $ propagateTurnRatio cAbove $ c1:cs

Split curves near the point.  Return the curves starting from this point.

> splitNearPoints :: DPoint -> Double -> SweepStateM ([Curve], [Curve])
> splitNearPoints p tol = do
>   curves1 <- splitNearDir withAbove p tol
>   curves2 <- splitNearDir withBelow p tol
>   return (map fst curves1 ++ map fst curves2,
>           map snd curves1 ++ map snd curves2)
>
> splitNearDir  :: ((Curve -> Maybe (Curve, Double))
>                   -> SweepStateM (Maybe (Curve, Double)))
>               -> DPoint -> Double
>               -> SweepStateM [(Curve, Curve)]
> splitNearDir dir p tol = do
>   mbSplit <- dir $ \curve ->
>     (,) curve <$>
>     pointOnCurve tol p
>     (view bezier curve)
>   case mbSplit of
>    Nothing -> return []
>    Just (curve, t) -> do
>      let (c1, c2) = splitBezier (view bezier curve) t
>          c1' = adjust curve $ adjustC3 p $
>                snapRound tol <$> c1
>          c2' = adjust curve $ adjustC0 p $
>                snapRound tol <$> c2
>      traceMessage $ "MSG Splitting curve " ++ showCurve (view bezier curve)
>      ((c1', c2'):) <$> splitNearDir dir p tol

Add the sorted curves starting at point to the Y-structure, and test
last curve with curve below.

> addMidCurves :: [Curve] -> Double -> SweepStateM ()
> addMidCurves [] _ = return ()
> addMidCurves [c] tol =
>   splitBelow c tol
> addMidCurves (c:cs) tol = do
>   addMidCurves cs tol
>   yStructAdd c 
>   
> splitBelow :: Curve -> Double -> SweepStateM ()
> splitBelow c tol = do
>   yStr <- use yStruct
>   p <- use focusPoint
>   let i = yStructIndex (singularC p) yStr
>   if i >= S.size yStr-1
>     then yStructAdd c
>     else
>       let cBelow = S.elemAt (i+1) yStr
>       in case splitMaybe c cBelow tol of
>         (Nothing, Nothing) -> 
>           yStructAdd c
>         (Nothing, Just (c3, c4)) -> do
>           assert (cubicC3 (_bezier c3) /= cubicC0 (_bezier c)) $
>             "splitBelow: curve starting in future: c3 == " ++
>              show c3 ++ " c == " ++ show c
>           traceMessage "MSG Split lower curve only." 
>           xStructAdd c4
>           yStructDel (i+1); yStructAdd c3; yStructAdd c
>         (Just (c1, c2), Nothing) -> do
>           traceMessage "MSG split curve above lower curve"
>           assert (cubicC3 (_bezier c1) /= cubicC0 (_bezier cBelow))
>             "SPLITBELOW: CURVE INTERSECTING POINTEVENT." 
>           xStructAdd c2
>           yStructAdd c1
>         (Just (c1, c2), Just (c3, c4)) -> do
>           traceMessage "MSG split lower curve and curve above."
>           xStructAdd c2; xStructAdd c4
>           yStructDel (i+1); yStructAdd c3; yStructAdd c1

If no curves start from the point, check if the surrounding
curves overlap.

> splitSurround :: Double -> SweepStateM ()
> splitSurround tol = do
>   p <- use focusPoint
>   yStr <- use yStruct
>   let i = yStructIndex (singularC p) yStr
>       s = S.size yStr
>   when (i >= 0 && i < s-1) $
>    case splitMaybe (S.elemAt i yStr) (S.elemAt (i+1) yStr) tol of
>      (Just (c1, c2), Just (c3, c4)) -> do
>        xStructAdd c2; xStructAdd c4
>        yStructDel i; yStructDel i
>        yStructAdd c3; yStructAdd c1
>      (Just (c1, c2), Nothing) -> do
>        xStructAdd c2
>        yStructDel i; yStructAdd c1
>      (Nothing, Just (c1, c2)) -> do
>        xStructAdd c2
>        yStructDel (i+1); yStructAdd c1
>      (Nothing, Nothing) ->
>        return ()


=== Find curve intersections

Test if both curves intersect.  Split one or both of the curves when
they intersect.  Also snapround each point, and make sure the point
of overlap is the same in both curves.

> splitMaybe :: Curve -> Curve -> Double ->
>               (Maybe (Curve, Curve),
>                Maybe (Curve, Curve))
> splitMaybe c1 c2 tol =
>   assertTrace (null errMsg) errMsg
>   (adjustSplit c1 <$> fst n,
>    adjustSplit c2 <$> snd n)
>   where
>     errMsg = checkSplitCurve b1 b2 n
>     n = nextIntersection b1 b2 tol $
>         bezierIntersection b1 b2 pTol
>     pTol = min (bezierParamTolerance b1 tol)
>            (bezierParamTolerance b2 tol)
>     b1 = view bezier c1
>     b2 = view bezier c2
>
> checkOverlap c1 c2 =
>   any (/=p) ps
>   where x0 = pointX $ cubicC0 c2
>         x3 = pointX $ cubicC3 c2
>         t0 | pointX (cubicC0 c1) >= pointX (cubicC0 c2) = 0
>            | otherwise = findX c1 (pointX (cubicC0 c2)) 1e-7
>         t1 | pointX (cubicC3 c1) <= pointX (cubicC3 c2) = 1
>            | otherwise = findX c1 (pointX (cubicC3 c2)) 1e-7
>         comp t = let (Point bx by) = evalBezier c1 (t0 + (t1-t0)*t/10)
>                  in compare (pointY (evalBezier c2 (findX c2 bx 1e-7))) by
>         (p:ps) = map comp [1..9]
> 
> checkDirection :: CubicBezier Double -> CubicBezier Double -> CubicBezier Double -> String
> checkDirection c1 c2 c@(CubicBezier p1 _ _ p2)
>   | PointEvent p1 < PointEvent p2 = ""
>   | otherwise = "Curve has wrong direction: " ++ showCurve c ++
>       "after splitting " ++ showCurve c1 ++ " " ++ showCurve c2
>   
> checkSplitCurve :: CubicBezier Double -> CubicBezier Double ->
>                    (Maybe (CubicBezier Double, CubicBezier Double),
>                      Maybe (CubicBezier Double, CubicBezier Double)) -> String
> checkSplitCurve c1 c2 (Nothing, Nothing) =
>   if checkOverlap c1 c2
>   then "Curves overlap: " ++ showCurve c1 ++ " " ++ showCurve c2
>   else ""
>
> checkSplitCurve c1 c2 (Just (c3, c4), Just (c5, c6)) =
>   checkDirection c1 c2 c3 ++ checkDirection c1 c2 c4 ++
>   checkDirection c1 c2 c5 ++ checkDirection c1 c2 c6
>
> 
> checkSplitCurve c1 c2 (Just (c3, c4), Nothing) =
>   checkDirection c1 c2 c3 ++ checkDirection c1 c2 c4 ++
>   (if cubicC3 c2 == cubicC3 c3 then "" else
>    "second curve doesn't split first curve: " ++ showCurve c1 ++ " " ++ showCurve c2)
>
> checkSplitCurve c1 c2 (Nothing, Just (c3, c4)) =
>   checkDirection c1 c2 c3 ++ checkDirection c1 c2 c4 ++
>   (if cubicC3 c1 == cubicC3 c3 then "" else
>    "first curve doesn't split second curve: " ++ showCurve c1 ++ " " ++ showCurve c2)
>
> adjustSplit :: Curve -> (CubicBezier Double, CubicBezier Double) -> (Curve, Curve)
> adjustSplit curve (b1, b2) = (adjust1 b1, adjust1 b2)
>   where
>     adjust1 b = (if PointEvent (cubicC0 b) > PointEvent (cubicC3 b)
>                 then revertCurve else id) $
>                set curveRank Nothing $ set bezier b curve
>
> revertCurve :: Curve -> Curve
> revertCurve (Curve bez tr chtr rank) =
>   Curve (reorient bez) (swap tr) (swap . chtr . swap) rank
>
> adjust :: Curve -> CubicBezier Double -> Curve
> adjust curve curve2 = set curveRank Nothing $ set bezier curve2 curve
>
> snapRoundBezier :: Double -> CubicBezier Double -> CubicBezier Double
> snapRoundBezier tol = fmap (snapRound tol)
>

Given a list of intersection parameters, split at the next
intersection, but don't split at the first or last control point, or
when the two curves are (nearly) coincident.  Note that list of
intersections is read lazily, in order not to evaluate more
intersections that necessary.

> nextIntersection :: CubicBezier Double -> CubicBezier Double -> Double -> [(Double, Double)]
>                  -> (Maybe (CubicBezier Double, CubicBezier Double),
>                      Maybe (CubicBezier Double, CubicBezier Double))
> nextIntersection _ _ _ [] = (Nothing, Nothing)
> nextIntersection b1@(CubicBezier p0 _ _ p3) b2@(CubicBezier q0 _ _ q3) tol ((t1, t2): ts)
>   | atStart1 && atStart2 =
>       nextIntersection b1 b2 tol ts
>   | bezierEqual b1l b2l tol =
>       nextIntersection b1 b2 tol ts
>   | otherwise =
>       assertTrace (vectorDistance x1 x2 < tol)
>       ("ASSERT: DISTANCE IS LARGER THAN TOLERANCE: " ++ show t1 ++ " " ++ show t2 ++ " " ++ showCurve b1 ++ " " ++ showCurve b2)
>       (bs1, bs2)
>   where
>     bs1 | atStart1 || atEnd1 = Nothing 
>         | otherwise = Just (adjustC3 pMid2 $ snapRoundBezier tol b1l,
>                             adjustC0 pMid2 $ snapRoundBezier tol b1r)
>     bs2 | atStart2 || atEnd2 = Nothing
>         | otherwise = Just (adjustC3 pMid2 $ snapRoundBezier tol b2l,
>                             adjustC0 pMid2 $ snapRoundBezier tol b2r)
>     pMid | atStart1 = cubicC0 b1
>          | atEnd1 = cubicC3 b1
>          | atStart2 = cubicC0 b2
>          | atEnd2 = cubicC3 b2
>          | otherwise = snapRound tol <$> cubicC3 b1l

The intersection point can be in the past (if the curve is nearly vertical), so if that happens move the intersection point a tiny bit aside.  We may need to reorient the subcurve after the intersection point as well (see adjustSplit).

>     pMid2 | PointEvent pMid <= PointEvent (cubicC0 b1) ||
>             PointEvent pMid <= PointEvent (cubicC0 b2) =
>             Point (max (pointX x1) (pointX x2) + tol) (pointY pMid)
>           | otherwise = pMid
>     x1 = evalBezier b1 t1
>     x2 = evalBezier b2 t2
>     atStart1 = vectorDistance (cubicC0 b1) x1 < 3*tol
>     atStart2 = vectorDistance (cubicC0 b2) x2 < 3*tol
>     atEnd1 = vectorDistance (cubicC3 b1) x1 < 3*tol
>     atEnd2 = vectorDistance (cubicC3 b2) x2 < 3*tol
>     (b1l, b1r) = splitBezier b1 t1
>     (b2l, b2r) = splitBezier b2 t2
> 
> adjustC0 :: Point a -> CubicBezier a -> CubicBezier a
> adjustC0 p (CubicBezier _ p1 p2 p3) = CubicBezier p p1 p2 p3
>
> adjustC3 :: Point a -> CubicBezier a -> CubicBezier a
> adjustC3 p (CubicBezier p0 p1 p2 _) = CubicBezier p0 p1 p2 p

=== Check if curves overlap.

If the curves overlap, combine the overlapping part into one curve.
To compare the curves, I first split the longest curve so that the
velocities in the first control point match, then compare those curves
for equality.

> curveOverlap :: Curve -> Curve -> Double
>              -> Maybe (Curve, Maybe Curve)
> curveOverlap c1 c2 tol
>   -- starting in the same point
>   | p0 /= q0 = Nothing
>   | colinear (view bezier c1) tol = if
>       | not $ colinear (view bezier c2) tol ->
>           Nothing
>       | vectorDistance (p3^-^p0)
>         ((q3^-^q0) ^* (d1/d2)) > tol ->
>           Nothing
>       | p3 == q3 -> 
>           Just (combineCurves c2 c1,
>                 Nothing)
>       | d1 > d2 ->
>           Just (combineCurves c2 c1,
>                 Just $ adjust c1 $
>                 CubicBezier q3
>                 (snapRound tol <$> interpolateVector q3 p3 (1/3))
>                 (snapRound tol <$> interpolateVector q3 p3 (2/3))
>                 p3)
>       | otherwise ->
>           Just (combineCurves c1 c2,
>                 Just $ adjust c2 $
>                 CubicBezier p3
>                 (snapRound tol <$> interpolateVector p3 q3 (1/3))
>                 (snapRound tol <$> interpolateVector p3 q3 (2/3))
>                 q3)
>   -- equalize velocities, and compare           
>   | v1 == 0 ||
>     v2 == 0 = Nothing
>   | v1 > v2 = if bezierEqual b2 b1l tol
>               then Just (combineCurves c2 c1,
>                          if checkEmpty b1r tol
>                          then Nothing
>                          else Just $ adjust c1 $
>                               adjustC0 (cubicC3 b2) $
>                               snapRoundBezier tol b1r)
>               else Nothing
>         
>   | otherwise =
>       if bezierEqual b1 b2l tol
>               then Just (combineCurves c1 c2,
>                          if checkEmpty b2r tol
>                          then Nothing
>                          else Just $ adjust c2 $
>                               adjustC0 (cubicC3 b1) $
>                               snapRoundBezier tol b2r)
>               else Nothing
>   where
>     (b1l, b1r) = splitBezier b1 (v2/v1)
>     (b2l, b2r) = splitBezier b2 (v1/v2)
>     b1@(CubicBezier p0 p1 _ p3) = view bezier c1
>     b2@(CubicBezier q0 q1 _ q3) = view bezier c2
>     d1 = vectorDistance p0 p3
>     d2 = vectorDistance q0 q3
>     v1 = vectorDistance p0 p1
>     v2 = vectorDistance q0 q1
>
> checkEmpty :: CubicBezier Double -> Double -> Bool
> checkEmpty (CubicBezier p0 p1 p2 p3) tol = 
>   p0 == p3 &&
>   vectorDistance p0 p1 < tol &&
>   vectorDistance p0 p2 < tol

Curves can be combined if they are equal, just by composing their
changeTurn functions.

> combineCurves :: Curve -> Curve -> Curve
> combineCurves c1 c2 =
>   over changeTurn (view changeTurn c2 .) c1

=== Snaprounding

> snapRound :: Double -> Double -> Double
> snapRound tol v =
>   fromInteger (round (v/tol)) * tol

=== Test if the point is on the curve (within tolerance) {#oncurve}

> pointOnCurve :: Double -> DPoint -> CubicBezier Double -> Maybe Double
> pointOnCurve tol p c1
>   | t <- closest c1 p tol,
>     p2 <- evalBezier c1 t,
>     vectorDistance p p2 < tol = Just t
>   | otherwise = Nothing

=== Testing beziers for approximate equality {#eq}

If the control points of two bezier curves are within a distance `eps`
from each other, then both curves will all so be at least within
distance `eps` from each other.  This can be proven easily:
subtracting both curves gives the distance curve.  Since each control
point of this curve lies within a circle of radius `eps`, by the
convex hull property, the curve will also be inside the circle, so the
distances between each point will never exceed `eps`.

> bezierEqual :: CubicBezier Double -> CubicBezier Double -> Double -> Bool
> bezierEqual cb1@(CubicBezier a0 a1 a2 a3) cb2@(CubicBezier b0 b1 b2 b3) tol
>   -- controlpoints equal within tol
>   | vectorDistance a1 b1 < tol &&
>     vectorDistance a2 b2 < tol &&
>     vectorDistance a3 b3 < tol &&
>     vectorDistance a0 b0 < tol = True
>   -- compare if both are colinear and close together
>   | dist < tol &&
>     colinear cb1 ((tol-dist)/2) &&
>     colinear cb2 ((tol-dist)/2) = True
>   | otherwise = False
>   where dist = max (abs $ ld b0) (abs $ ld b3)
>         ld = lineDistance (Line a0 a3)

Higher level functions
----------------------

> fillFunction :: FillRule -> Int -> Bool
> fillFunction NonZero = (/=0)
> fillFunction EvenOdd = odd
>
> newSweep :: XStruct -> SweepState
> newSweep xStr = SweepState M.empty S.empty undefined xStr
>
> runSweep :: SweepState -> SweepStateM () -> SweepState
#if DEBUG
> runSweep sweep m = 
>   unsafePerformIO $ do
>   hPutStrLn stderr "XSTRUCTBEGIN" 
>   mapM_ (hPutStrLn stderr . showCurve) $ map (view bezier) $
>    concat $ M.elems $ view xStruct sweep
>   hPutStrLn stderr "XSTRUCTEND" 
>   execStateT m sweep
#else
> runSweep sweep m =
>   execState m sweep
#endif

> -- | `O((n+m)*log(n+m))`, for `n` segments and `m` intersections.
> -- Union of paths, removing overlap and rounding to the given
> -- tolerance.
> union :: [ClosedPath Double] -- ^ Paths
>          -> FillRule         -- ^ input fillrule
>          -> Double           -- ^ Tolerance
>          -> [ClosedPath Double]
> union p fill tol =
>   outputPaths $ view output $ runSweep sweep $ 
>   loopEvents (fillFunction fill . fst) tol
>   where
>     sweep = newSweep xStr
>     xStr = makeXStruct (over _1 $ subtract 1) (over _1 (+1)) tol $
>            concatMap closedPathCurves p
>
> union' p fill tol =
>   outputPaths $ view output $ runSweep sweep $ 
>   loopEvents (fillFunction fill . fst) tol
>   where
>     sweep = newSweep xStr
>     xStr = makeXStruct (over _1 $ subtract 1) (over _1 (+1)) tol p
>
> 
> -- | `O((n+m)*log(n+m))`, for `n` segments and `m` intersections.
> -- Combine paths using the given boolean operation
> boolPathOp :: (Bool -> Bool -> Bool) -- ^ operation
>           -> [ClosedPath Double]     -- ^ first path (merged with union)
>           -> [ClosedPath Double]     -- ^ second path (merged with union)
>           -> FillRule                -- ^ input fillrule
>           -> Double                  -- ^ tolerance 
>           -> [ClosedPath Double]
> boolPathOp op p1 p2 fill tol =
>   outputPaths $ view output $ runSweep sweep $
>   loopEvents isInside tol
>   where
>     isInside (a, b) = fillFunction fill a `op`
>                       fillFunction fill b
>     sweep = newSweep xStr
>
>     xStr = M.unionWith (++)
>            (makeXStruct 
>             (over _1 (subtract 1))
>             (over _1 (+1)) tol $
>             concatMap closedPathCurves p1)
>            (makeXStruct
>             (over _2 (subtract 1))
>             (over _2 (+1)) tol $
>             concatMap closedPathCurves p2)
>
> intersection, difference, exclusion ::
>   [ClosedPath Double] -> [ClosedPath Double] ->
>   FillRule -> Double -> [ClosedPath Double]
>
> -- | path intersection  
> intersection = boolPathOp (&&)
>
> -- | path difference
> difference = boolPathOp (\a b -> a && not b)
>
> -- | path exclusion
> exclusion = boolPathOp (\a b -> if a then not b else b)
>

handy for debugging: 

>
> -- mkBezier (a, b) (c, d) (e, f) (g, h) = CubicBezier (Point a b) (Point c d) (Point e f) (Point g h)
> --x = 
