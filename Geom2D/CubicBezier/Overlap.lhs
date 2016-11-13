> {-# LANGUAGE MultiWayIf, PatternGuards, TemplateHaskell, BangPatterns #-}

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
    - `field %= fun`: `state.field = fun (state.field)`

Let's begin with declaring the module and library imports:

> module Geom2D.CubicBezier.Overlap
>        (boolPathOp, union, intersection, difference,
>         exclusion, FillRule (..))
>        where
> import Prelude hiding (mapM)
> import Geom2D
> import Geom2D.CubicBezier.Basic
> import Geom2D.CubicBezier.Intersection
> import Math.BernsteinPoly
> import Data.Traversable (mapM)
> import Data.Functor ((<$>))
> import Data.List (sortBy, sort)
> import Control.Monad.State hiding (mapM)
> import Lens.Micro
> import Lens.Micro.TH
> import Lens.Micro.Mtl
> import qualified Data.Map.Strict as M
> import qualified Data.Set as S

The basic idea is to keep curves where one side is inside the filled
region, and the other side is outside, and discard the rest. 
Since that could be true only of a part of the curve, we also need to
split each curve when it intersects
another curve.  How to know which side is the inside, and which
side the outside?  There are two methods which are use the most: the
[*even-odd rule*](https://en.wikipedia.org/wiki/Even%E2%80%93odd_rule)
and the [*nonzero rule*](https://en.wikipedia.org/wiki/Nonzero-rule).
Instead of hardwiring it, I use higher-order functions to determine
when a turnratio is inside the region to be filled, and how the
turnratio changes with each curve.

Checking each pair of curves for intersections would work, but is
rather inefficient.  We only need to check for overlap when two curves
are adjacent.  Fortunately there exist a good method from
*computational geometry*, called the *sweep line algorithm*.  The 
idea is to sweep a vertical line over the input, starting from
leftmost point to the right (of course the opposite direction is also
possible), and to update the input dynamically.  We keep track of each
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
>                    deriving Show

When the x-coordinates are equal, use the y-coordinate to determine
the order.

> instance Eq PointEvent where
>   (PointEvent (Point x1 y1)) == (PointEvent (Point x2 y2)) =
>     (x1, y1) == (x2, y2)
>
> instance Ord PointEvent where
>   compare (PointEvent (Point x1 y1)) (PointEvent (Point x2 y2)) =
>     compare (x1, y2) (x2, y1)

All curves are kept left to right, so we need to remember the
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

The FillRule datatype is used for the exported API:

> data FillRule = EvenOdd | NonZero

> data Curve = Curve {
>   _bezier :: !(CubicBezier Double),
>   _turnRatio :: !(Int, Int),
>   _changeTurn :: !((Int, Int) -> (Int, Int))}
>
> trOne :: (Int, Int)
> trOne = (0,0)
> 
> makeLenses ''Curve
>
> instance Show Curve where
>   show (Curve b a _) =
>     "Curve " ++ show b ++ " " ++ show a
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
>   _yStructLeft :: !YStruct,
>   _yStructRight :: !YStruct,
>   _xStruct :: !XStruct}
>                   deriving Show
>                   
> makeLenses ''SweepState

Changing the focus point can be done efficiently in `O(log n)` by
mering and splitting again:

> changeFocus :: DPoint -> SweepState -> SweepState
> changeFocus p sweep =
>   let (lStr, rStr) =
>         S.split (Curve (CubicBezier p p p p) trOne id) $
>         S.union (view yStructLeft sweep) (view yStructRight sweep)
>   in set yStructLeft lStr $
>      set yStructRight rStr
>      sweep

This handy helper function will pass the first curve above to the
given function, and if it doesn't return `Nothing`, remove it from the
state.  It does nothing when there is no curve above.

> withAbove :: (Curve -> Maybe a) -> State SweepState (Maybe a)
> withAbove f = do
>   lStr <- use yStructLeft
>   if S.null lStr
>     then return Nothing
>     else let (c, lStr') = S.deleteFindMax lStr
>          in case f c of
>              Nothing ->
>                return Nothing
>              Just x -> do
>                yStructLeft .= lStr'
>                return $ Just x

The same with the curve below.

> withBelow :: (Curve -> Maybe a) -> State SweepState (Maybe a)
> withBelow f = do
>   rStr <- use yStructRight
>   if S.null rStr
>     then return Nothing
>     else let (c, rStr') = S.deleteFindMin rStr
>          in case f c of
>              Nothing ->
>                return Nothing
>              Just x -> do
>                yStructRight .= rStr'
>                return $ Just x

`splitYStruct` changes the focus and returns and removes any curves which end in
the current pointEvent:

> splitYStruct :: DPoint -> State SweepState [Curve]
> splitYStruct p = do
>   modify $ changeFocus p
>   let go = do
>         mbC <- withAbove $ \c ->
>           -- remove and return c if it ends in point p
>           
>           guard (cubicC3 (_bezier c) == p) >> Just c
>         case mbC of
>          Just c ->
>            (c:) <$> go
>          Nothing -> return []
>   go


=== Some functions on the Sweep state:

Adding and removing curves from the X structure.

> insertX :: PointEvent -> [Curve] -> SweepState -> SweepState
> insertX p c =
>   over xStruct $ M.insertWith (++) p c
>
> xStructAdd :: Curve -> SweepState -> SweepState
> xStructAdd c =
>   insertX (PointEvent $ cubicC0 $
>                         view bezier c) [c]
>
> xStructRemove :: State SweepState (PointEvent, [Curve])
> xStructRemove = zoom xStruct $ state M.deleteFindMin

To compare curves vertically, take the the curve which starts the
rightmost, and see if it falls below or above the curve.  If the first control points are coincident, test the
last control points instead. The curves in the Y-structure shouldn't
intersect (except in the endpoints), so these cases don't have to be
handled.  To lookup a single point, I use a singular bezier curve.

> instance Eq Curve where
>   Curve c1 t1 ct1 == Curve c2 t2 ct2 =
>     c1 == c2 && t1 == t2 && ct1 (ct2 t1) == t1
>     
> instance Ord Curve where
>   compare (Curve c1@(CubicBezier p0 p1 p2 p3) tr1 _)
>     (Curve c2@(CubicBezier q0 q1 q2 q3) tr2 _)
>     | p0 == q0 = if
>         | p3 == q3 ->
>             -- compare the midpoint
>             case (compVert (evalBezier c1 0.5) c2) of
>              LT -> LT
>              GT -> GT
>              EQ ->
>                -- otherwise arbitrary
>                compare (tr1, PointEvent p1, PointEvent p2)
>                (tr2, PointEvent q1, PointEvent q2)
>         | pointX p3 < pointX q3 ->
>             case (compVert p3 c2) of
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
>       case (compVert p0 c2) of
>       LT -> LT
>       EQ -> GT
>       GT -> GT

Compare a point with a curve.  See if it falls below or above the hull
first.  Otherwise find the point on the curve with the same
X-coordinate by solving a cubic equation.

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
>     t = findX x1 c1 $
>         maximum (map maxp [p0, p1, p2, p3])*1e-12
>     maxp (Point x y) = max (abs x) (abs y)
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
>            LT -> [(PointEvent p0, [Curve c trOne chTr])]
>            GT -> [(PointEvent p3, [Curve (reorient c) trOne chTrBack]),
>                   (PointEvent p0, [])]
>            -- vertical curve
>            EQ | pointY p0 > pointY p3 ->
>                 [(PointEvent p0, [Curve c trOne chTr])]
>               | otherwise ->
>                 [(PointEvent p3, [Curve (reorient c) trOne chTrBack]),
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

For the main loop, we remove the leftmost point from the
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

> loopEvents :: ((Int, Int) -> Bool) -> Double -> SweepState -> SweepState
> loopEvents isInside tol sweep 
>   | M.null $ view xStruct sweep = sweep
>   | otherwise =
>     loopEvents isInside tol $!
>     flip execState sweep $ do
>       -- remove leftmost point from X structure
>       (PointEvent p, curves) <- xStructRemove
>       -- change focus, and remove curves ending at current
>       -- pointevent from Y structure
>       ending <- splitYStruct p
>       -- split near curves
>       (ending2, rightSubCurves) <- splitNearPoints p tol
>       -- output curves to the left of the sweepline.
>       modify $ filterOutput (ending ++ ending2) isInside 
>       let allCurves = rightSubCurves ++ curves
>       if null allCurves
>          -- split surrounding curves
>         then splitSurround tol
>         else do
>         -- sort curves
>         sorted <- splitAndOrder tol allCurves
>         -- split curve above
>         curves2 <- splitAbove sorted tol
>         -- add curves to Y structure
>         addMidCurves curves2 tol

Send curves to output
---------------------

> outputPaths :: (M.Map PointEvent [CubicBezier Double]) -> [ClosedPath Double]
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
>         let ((PointEvent p0, (c0:cs)), m0) =
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

> filterOutput :: [Curve] -> ((Int, Int) -> Bool) -> SweepState -> SweepState
> filterOutput curves isInside sweep =
>   foldl (flip $ outputCurve isInside) sweep curves
>
> outputCurve :: ((Int, Int) -> Bool) -> Curve -> SweepState -> SweepState
> outputCurve isInside (Curve c tr op)
>   | isInside (op tr) /= isInside tr =
>       let c' | isInside tr = reorient c
>              | otherwise = c
>       in over output (M.insertWith (++) (PointEvent $ cubicC0 c') [c'])
>   | otherwise = id

Test for intersections and split:
---------------------------------

Since the curves going out of the current pointEvent in the X-structure are
unordered, they need to be ordered first.  First they are ordered by
first derivative.  Since it's easier to compare two curves when they
don't overlap, remove overlap, and then sort again by comparing the
whole curve.

To do this, I implemented a monadic insertion sort.  First the curves are split
in the statemonad, then they are compared.

> splitAndOrder :: Double -> [Curve] -> State SweepState [Curve]
> splitAndOrder tol curves =
>   sortSplit tol $
>   sortBy compDeriv curves
>
> compDeriv :: Curve -> Curve -> Ordering
> compDeriv (Curve (CubicBezier p0 p1 _ _) _ _)
>   (Curve (CubicBezier q0 q1 _ _) _ _) =
>   compare (vectorCross (p1^-^p0) (q1^-^ q0)) 0

Insertion sort, by splitting and comparing.  This should be efficient
enough, since ordering by derivative should mostly order the curves.

> sortSplit :: Double -> [Curve] -> State SweepState [Curve]
> sortSplit _ [] = return []
> sortSplit tol (x:xs) =
>   insertM x tol =<<
>   sortSplit tol xs
>
> insertM :: Curve -> Double -> [Curve] -> State SweepState [Curve]
> insertM x _ [] = return [x]
> insertM x tol (y:ys) =
>   case curveOverlap x y tol of
>    Just (c1, c2) -> do
>      mapM (modify . xStructAdd) c2
>      insertM c1 tol ys
>    Nothing -> do
>      (x', y') <- splitM x y tol
>      if x' < y'
>        then return (x':y':ys)
>        else (y':) <$> insertM x' tol ys
>
> splitM :: Curve -> Curve -> Double -> State SweepState (Curve, Curve)
> splitM x y tol =
>   case splitMaybe x y tol of
>   (Just (a, b), Just (c, d)) -> do
>     modify $ insertX (PointEvent $ cubicC0 $ view bezier b) [b, d]
>     return (a, c)
>   (Nothing, Just (c, d)) -> do
>     modify $ insertX (PointEvent $ cubicC0 $ view bezier d) [d]
>     return (x, c)
>   (Just (a, b), Nothing) -> do
>     modify $ insertX (PointEvent $ cubicC0 $ view bezier b) [b]
>     return (a, y)
>   (Nothing, Nothing) ->
>     return (x, y)

Handle intersections of the first curve at point and the curve
above. Return the curves with updated turnratios.  Some care is needed
when one of the curves is intersected at the endpoints, in order not
to create singular curves.

> updateTurnRatio :: Curve -> Curve -> Curve
> updateTurnRatio (Curve _ tr chTr) =
>   set turnRatio (chTr tr)
>
> propagateTurnRatio :: Curve -> [Curve] -> [Curve]
> propagateTurnRatio cAbove l =
>   tail $ scanl updateTurnRatio cAbove l
>
> splitAbove :: [Curve] -> Double -> State SweepState [Curve]
> splitAbove [] _ = return []
> splitAbove (c:cs) tol = do
>   lStr <- use yStructLeft
>   if S.null lStr
>     then let c' = set turnRatio trOne c
>          in return $ c':propagateTurnRatio c' cs
>     else do
>     let (cAbove, lStr') = S.deleteFindMax lStr
>     case splitMaybe c cAbove tol of
>      (Nothing, Nothing) ->
>        return $ propagateTurnRatio cAbove $ c:cs
>      (Just (c1, c2), Nothing) ->
>        if cubicC3 (_bezier c1) == cubicC0 (_bezier cAbove)
>        then do
>          modify $ xStructAdd cAbove . xStructAdd c2
>          yStructLeft .= lStr'
>          return $ propagateTurnRatio cAbove $ c1:cs
>        else do
>          modify $ xStructAdd c2
>          return $ propagateTurnRatio cAbove $ c1:cs
>      (Nothing, Just (c3, c4)) ->
>        if cubicC3 (_bezier c3) == cubicC0 (_bezier c)
>        then error "curve intersecting pointevent"
>        else do
>          modify $ xStructAdd c4
>          yStructLeft .= S.insert c3 lStr'
>          return $ propagateTurnRatio cAbove $ c:cs
>      (Just (c1, c2), Just (c3, c4)) -> do
>        modify $ xStructAdd c2 . xStructAdd c4
>        yStructLeft .= S.insert c3 lStr'
>        return $ propagateTurnRatio cAbove $ c1:cs

Split curves near the point.  Return the curves starting from this point, and the index of the last split point

> splitNearPoints :: DPoint -> Double -> State SweepState ([Curve], [Curve])
> splitNearPoints p tol = do
>   curves1 <- splitNearDir withAbove p tol
>   curves2 <- splitNearDir withBelow p tol
>   return (map fst curves1 ++ map fst curves2,
>           map snd curves1 ++ map snd curves2)
>
> splitNearDir  :: ((Curve -> Maybe (Curve, Double))
>                   -> State SweepState (Maybe (Curve, Double)))
>               -> DPoint -> Double
>               -> State SweepState [(Curve, Curve)]
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
>      ((c1', c2'):) <$> splitNearDir dir p tol

Add the sorted curves starting at point to the Y-structure, and test
last curve with curve below.

> addMidCurves :: [Curve] -> Double -> State SweepState ()
> addMidCurves [] _ = return ()
> addMidCurves [c] tol =
>   splitBelow c tol
> addMidCurves (c:cs) tol = do
>   yStructLeft %= S.insert c 
>   addMidCurves cs tol
>   
> splitBelow :: Curve -> Double -> State SweepState ()
> splitBelow c tol = do
>   rStr <- use yStructRight
>   let (cBelow, rStr') = S.deleteFindMin rStr
>   if S.null rStr
>     then yStructLeft %= S.insert c
>     else
>     case splitMaybe c cBelow tol of
>      (Nothing, Nothing) ->
>        yStructLeft %= S.insert c
>      (Nothing, Just (c3, c4)) ->
>        if cubicC3 (_bezier c3) == cubicC0 (_bezier c)
>        then error "internal error: splitBelow: curve starting in future"
>        else do
>          modify $ xStructAdd c4
>          yStructLeft %= S.insert c . S.insert c3
>          yStructRight .= rStr'
>      (Just (c1, c2), Nothing) ->
>        if cubicC3 (_bezier c1) == cubicC0 (_bezier cBelow)
>        then error "internal error: splitBelow: curve intersecting pointevent."
>        else do
>          modify $ xStructAdd c2
>          yStructLeft %= S.insert c1
>      (Just (c1, c2), Just (c3, c4)) -> do
>        modify $ xStructAdd c2 . xStructAdd c4
>        yStructLeft %= S.insert c1 . S.insert c3
>        yStructRight .= rStr'

If no curves start from the point, we have to check if the surrounding
curves overlap.

> splitSurround :: Double -> State SweepState ()
> splitSurround tol = do
>   lStr <- use yStructLeft
>   rStr <- use yStructRight
>   if S.null lStr || S.null rStr
>     then return ()
>     else do
>     let (cAbove, lStr') = S.deleteFindMax lStr
>         (cBelow, rStr') = S.deleteFindMin rStr
>     case splitMaybe cAbove cBelow tol of
>      (Just (c1, c2), Just (c3, c4)) -> do
>        modify $ xStructAdd c2 .
>          xStructAdd c4
>        yStructLeft .= S.insert c1 lStr'
>        yStructRight .= S.insert c3 rStr'
>      (Just (c1, c2), Nothing) -> do
>        modify $ xStructAdd c2
>        yStructLeft .= S.insert c1 lStr'
>      (Nothing, Just (c1, c2)) -> do
>        modify $ xStructAdd c2
>        yStructRight .= S.insert c1 rStr'
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
>   (adjustSplit c1 <$> fst n,
>    adjustSplit c2 <$> snd n)
>   where
>     n = nextIntersection b1 b2 tol $
>         bezierIntersection b1 b2 pTol
>     pTol = min (bezierParamTolerance b1 tol)
>            (bezierParamTolerance b2 tol)
>     b1 = view bezier c1
>     b2 = view bezier c2
>
> adjustSplit :: Curve -> (CubicBezier Double, CubicBezier Double) -> (Curve, Curve)
> adjustSplit curve (b1, b2)   =
>   (set bezier b1 curve,
>    set bezier b2 curve)
>
> adjust :: Curve -> CubicBezier Double -> Curve
> adjust curve curve2 = set bezier curve2 curve
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
>   | atStart1 =
>     (Nothing,
>      Just (adjustC3 p0 $ snapRoundBezier tol b2l,
>            adjustC0 p0 $ snapRoundBezier tol b2r))
>   | atStart2 =
>     (Just (adjustC3 q0 $ snapRoundBezier tol b1l,
>            adjustC0 q0 $ snapRoundBezier tol b1r),
>      Nothing)
>   | atEnd1 && atEnd2 = (Nothing, Nothing)
>   | atEnd1 =
>     (Nothing,
>      Just (adjustC3 p3 $ snapRoundBezier tol b2l,
>            adjustC0 p3 $ snapRoundBezier tol b2r))
>   | atEnd2 =
>     (Just (adjustC3 q3 $ snapRoundBezier tol b1l,
>            adjustC0 q3 $ snapRoundBezier tol b1r),
>      Nothing)
>   | bezierEqual b1l b2l tol =
>       nextIntersection b1 b2 tol ts
>   | otherwise =
>       let pMid = snapRound tol <$> cubicC3 b1l
>       in (Just (snapRoundBezier tol b1l,
>                 snapRoundBezier tol b1r),
>           Just (adjustC3 pMid $ snapRoundBezier tol b2l,
>                 adjustC0 pMid $ snapRoundBezier tol b2r))
>    where
>      x1 = evalBezier b1 t1
>      x2 = evalBezier b2 t2
>      atStart1 = vectorDistance (cubicC0 b1) x1 < tol
>      atStart2 = vectorDistance (cubicC0 b2) x2 < tol
>      atEnd1 = vectorDistance (cubicC3 b1) x1 < tol
>      atEnd2 = vectorDistance (cubicC3 b2) x2 < tol
>      (b1l, b1r) = splitBezier b1 t1
>      (b2l, b2r) = splitBezier b2 t2
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
>   | (t:_) <-
>     closest c1 p tol,
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

=== Finding the on curve point at the X coordinate {#findx}

Solve a cubic equation to find the X coordinate.  This should be
converted to a closed form solver for efficiency.

> findX :: Double -> CubicBezier Double -> Double -> Double
> findX x c@(CubicBezier p0 p1 p2 p3) eps =
>   head $ bezierFindRoot bez 0 1 $
>   bezierParamTolerance c (eps/10)
>   where bez = listToBernstein 
>               (map pointX [p0, p1, p2, p3]) ~-
>               listToBernstein [x, x, x, x]

Higher level functions
----------------------

> fillFunction :: FillRule -> Int -> Bool
> fillFunction NonZero = (>0)
> fillFunction EvenOdd = odd

> -- | `O((n+m)*log(n+m))`, for `n` segments and `m` intersections.
> -- Union of paths, removing overlap and rounding to the given
> -- tolerance.
> union :: [ClosedPath Double] -- ^ Paths
>          -> FillRule         -- ^ input fillrule
>          -> Double           -- ^ Tolerance
>          -> [ClosedPath Double]
> union p fill tol =
>   outputPaths out
>   where
>     out = view output $
>           loopEvents (fillFunction fill . fst) tol sweep
>     sweep = SweepState M.empty S.empty S.empty xStr
>     xStr = makeXStruct (over _1 $ subtract 1) (over _1 (+1)) tol $
>            concatMap closedPathCurves p
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
>   outputPaths $ view output $
>   loopEvents isInside tol sweep
>   where
>     isInside (a, b) = fillFunction fill a `op`
>                       fillFunction fill b
>     sweep = SweepState M.empty S.empty S.empty xStr
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
