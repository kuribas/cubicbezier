{-# LANGUAGE BangPatterns, DeriveFunctor #-}
-- | This module implements an extension to paths as used in
-- D.E.Knuth's /Metafont/.  Metafont gives an alternate way
-- to specify paths using bezier curves.  I'll give a brief overview of
-- the metafont curves.  A more in depth explanation can be found in 
-- /The MetafontBook/.
-- 
-- Each spline has a tension parameter, which is a relative measure of
-- the length of the curve.  You can specify the tension for the left
-- side and the right side of the spline separately.  By default
-- metafont gives a tension of 1, which gives a good looking curve.
-- Tensions shouldn't be less than 3/4, but this implementation
-- doesn't check for it.  If you want to avoid points of inflection on
-- the spline, you can use @TensionAtLeast@ instead of @Tension@,
-- which will adjust the length of the control points so they fall
-- into the /bounding triangle/, if such a triangle exist.
--
-- You can either give directions for each node, or let metafont find
-- them.  Metafont will solve a set of equations to find the
-- directions.  You can also let metafont find directions at corner
-- points by setting the /curl/, which is how much the point /curls/
-- at that point.  At endpoints a curl of 1 is implied when it is not
-- given.
--
-- Metafont will then find the control points from the path for you.
-- You can also specify the control points explicitly.
--
-- Here is an example path from the metafont program text:
-- 
-- @
-- z0..z1..tension atleast 1..{curl 2}z2..z3{-1,-2}..tension 3 and 4..z4..controls z45 and z54..z5
-- @
-- 
-- This path is equivalent to:
--
-- @
-- z0{curl 1}..tension atleast 1 and atleast 1..{curl 2}z2{curl 2}..tension 1 and 1..
-- {-1,-2}z3{-1,-2}..tension 3 and 4..z4..controls z45 and z54..z5
-- @
--
-- This path can be used with the following datatype:
-- 
-- @
-- OpenMetaPath [ (z0, MetaJoin Open (Tension 1) (Tension 1) Open)
--              , (z1, MetaJoin Open (TensionAtLeast 1) (TensionAtLeast 1) (Curl 2))
--              , (z2, MetaJoin Open (Tension 1) (Tension 1) Open)
--              , (z3, MetaJoin (Direction (Point (-1) (-2))) (Tension 3) (Tension 4) Open)
--              , (z4, Controls z45 z54)
--              ] z5
-- @
--
-- Cyclic paths are similar, but use the @CyclicMetaPath@ contructor.
-- There is no ending point, since the ending point will be the same
-- as the first point.

module Geom2D.CubicBezier.MetaPath
       (unmeta, MetaPath (..), MetaJoin (..), MetaNodeType (..), Tension (..))
where
import Geom2D
import Geom2D.CubicBezier.Basic
import Data.List
import Text.Printf
import qualified Data.Vector as V

data MetaPath a = OpenMetaPath [(Point a, MetaJoin a)] (Point a)
                -- ^ A metapath with endpoints
              | CyclicMetaPath [(Point a, MetaJoin a)]
                -- ^ A metapath with cycles.  The last join joins the
                -- last point with the first.
              deriving (Eq, Functor)

data MetaJoin a = MetaJoin { metaTypeL :: MetaNodeType a
                           -- ^ The nodetype going out of the
                           -- previous point.  The metafont default is
                           -- @Open@.
                         , tensionL :: Tension
                           -- ^ The tension going out of the previous point.
                           -- The metafont default is 1.
                         , tensionR :: Tension
                           -- ^ The tension going into the next point.
                           -- The metafont default is 1.
                         , metaTypeR :: MetaNodeType a
                           -- ^ The nodetype going into the next point.
                           -- The metafont default is @Open@.
                         }
              | Controls (Point a) (Point a)
                -- ^ Specify the control points explicitly.
              deriving (Show, Eq, Functor)

data MetaNodeType a = Open
                    -- ^ An open node has no direction specified.  If
                    -- it is an internal node, the curve will keep the
                    -- same direction going into and going out from
                    -- the node.  If it is an endpoint or corner
                    -- point, it will have curl of 1.
                  | Curl {curlgamma :: Double}
                    -- ^ The node becomes and endpoint or a corner
                    -- point.  The curl specifies how much the segment
                    -- `curves`.  A curl of `gamma` means that the
                    -- curvature is `gamma` times that of the
                    -- following node.
                  | Direction {nodedir :: Point a}
                    -- ^ The node has a given direction.
                  deriving (Eq, Show, Functor)

data Tension = Tension {tensionValue :: Double}
               -- ^ The tension value specifies how /tense/ the curve is.
               -- A higher value means the curve approaches a line
               -- segment, while a lower value means the curve is more
               -- round.  Metafont doesn't allow values below 3/4.
             | TensionAtLeast {tensionValue :: Double}
               -- ^ Like @Tension@, but keep the segment inside the
               -- bounding triangle defined by the control points, if
               -- there is one.
             deriving (Eq, Show)

instance Show a => Show (MetaPath a) where
  show (CyclicMetaPath nodes) =
    showPath nodes ++ "cycle"
  show (OpenMetaPath nodes lastpoint) =
    showPath nodes ++ showPoint lastpoint

showPath :: Show a => [(Point a, MetaJoin a)] -> String
showPath = concatMap showNodes
  where
    showNodes (p, Controls u v) =
      showPoint p ++ "..controls " ++ showPoint u ++ "and " ++ showPoint v ++ ".."
    showNodes (p, MetaJoin m1 t1 t2 m2) =
      showPoint p ++ typename m1 ++ ".." ++ tensions ++ typename m2
      where
        tensions
          | t1 == t2 && t1 == Tension 1 = ""
          | t1 == t2 = printf "tension %s.." (showTension t1)
          | otherwise = printf "tension %s and %s.."
                        (showTension t1) (showTension t2)
    showTension (TensionAtLeast t) = printf "atleast %.3f" t :: String
    showTension (Tension t) = printf "%.3f" t :: String
    typename Open = ""
    typename (Curl g) = printf "{curl %.3f}" g :: String
    typename (Direction dir) = printf "{%s}" (showPoint dir) :: String
    
showPoint :: Show a => Point a -> String
showPoint (Point x y) = "(" ++ show x ++ ", " ++ show y ++ ")"

-- | Create a normal path from a metapath.
unmeta :: MetaPath Double -> Path Double
unmeta (OpenMetaPath nodes endpoint) =
  unmetaOpen (flip sanitize endpoint $ removeEmptyDirs nodes) endpoint

unmeta (CyclicMetaPath nodes) =
  case spanList bothOpen (removeEmptyDirs nodes) of
    ([], []) -> error "empty metapath"
    (l, []) -> if fst (last l) == fst (head l)
               then unmetaAsOpen l []
               else unmetaCyclic l
    (l, m:n) ->
      if leftOpen (m:n)
      then unmetaAsOpen (l++[m]) n
      else unmetaAsOpen l (m:n)

-- solve a cyclic metapath as an open path if possible.
-- rotate to the defined node, and rotate back after
-- solving the path.
unmetaAsOpen :: [(DPoint, MetaJoin Double)] -> [(DPoint, MetaJoin Double)] -> Path Double
unmetaAsOpen l m = ClosedPath (l'++m') 
  where n = length m
        OpenPath o _ =
          unmetaOpen (sanitizeCycle (m++l)) (fst $ head (m ++l))
        (m',l') = splitAt n o

unmetaOpen :: [(DPoint, MetaJoin Double)] -> DPoint -> Path Double
unmetaOpen nodes endpoint =
  let subsegs = openSubSegments nodes endpoint
      path = joinSegments $ map unmetaSubSegment subsegs
  in OpenPath path endpoint

-- decompose into a list of subsegments that need to be solved.
openSubSegments :: [(DPoint, MetaJoin Double)] -> DPoint -> [MetaPath Double]
openSubSegments [] _   = []
openSubSegments l lastPoint =
  case spanList (not . breakPoint) l of
    (m, n:o) ->
      let point = case o of
            ((p,_):_) -> p
            _ -> lastPoint
      in OpenMetaPath (m ++ [n]) point :
         openSubSegments o lastPoint
    _ -> error "openSubSegments': unexpected end of segments"

-- join subsegments into one segment
joinSegments :: [Path Double] -> [(DPoint, PathJoin Double)]
joinSegments = concatMap nodes
  where nodes (OpenPath n _) = n
        nodes (ClosedPath n) = n

-- solve a cyclic metapath where all angles depend on the each other.
unmetaCyclic :: [(DPoint, MetaJoin Double)] -> Path Double
unmetaCyclic nodes =
  let points = map fst nodes
      chords = zipWith (^-^) (tail $ cycle points) points
      tensionsA = map (tensionL . snd) nodes
      tensionsB = map (tensionR . snd) nodes
      turnAngles = zipWith turnAngle chords (tail $ cycle chords)
      thetas = solveCyclicTriD $
               eqsCycle tensionsA
               points
               tensionsB
               turnAngles
      phis = zipWith (\x y -> -(x+y)) turnAngles (tail $ cycle thetas)
  in ClosedPath $ zip points $
     zipWith6 unmetaJoin points (tail $ cycle points)
     thetas phis tensionsA tensionsB

-- solve a subsegment
unmetaSubSegment :: MetaPath Double -> Path Double

-- the simple case where the control points are already given.
unmetaSubSegment (OpenMetaPath [(p, Controls u v)] q) =
  OpenPath [(p, JoinCurve u v)] q

-- otherwise solve the angles, and find the control points
unmetaSubSegment (OpenMetaPath nodes lastpoint) =
  let points = map fst nodes ++ [lastpoint]
      joins = map snd nodes
      chords = zipWith (^-^) (tail points) points
      tensionsA = map tensionL joins
      tensionsB = map tensionR joins
      turnAngles = zipWith turnAngle chords (tail chords) ++ [0]
      thetas = solveTriDiagonal $
               eqsOpen points joins chords turnAngles
               (map tensionValue tensionsA)
               (map tensionValue tensionsB)
      phis = zipWith (\x y -> -(x+y)) turnAngles (tail thetas)
      pathjoins =
        zipWith6 unmetaJoin points (tail points) thetas phis tensionsA tensionsB
  in OpenPath (zip points pathjoins) lastpoint

unmetaSubSegment _ = error "unmetaSubSegment: subsegment should not be cyclic"

removeEmptyDirs :: [(DPoint, MetaJoin Double)] -> [(DPoint, MetaJoin Double)]
removeEmptyDirs = map remove
  where remove (p, MetaJoin (Direction (Point 0 0)) tl tr jr) = remove (p, MetaJoin Open tl tr jr)
        remove (p, MetaJoin jl tl tr (Direction (Point 0 0))) = (p, MetaJoin jl tl tr Open)
        remove j = j

-- if p == q, it will become a control point
bothOpen :: [(DPoint, MetaJoin Double)] -> Bool
bothOpen ((p, MetaJoin Open _ _ Open):(q, _):_) = p /= q  
bothOpen [(_, MetaJoin Open _ _ Open)] = True
bothOpen _ = False

leftOpen :: [(DPoint, MetaJoin Double)] -> Bool
leftOpen ((p, MetaJoin Open _ _ _):(q, _):_) = p /= q  
leftOpen [(_, MetaJoin Open _ _ _)] = True
leftOpen _ = False

sanitizeCycle :: [(DPoint, MetaJoin Double)] -> [(DPoint, MetaJoin Double)]
sanitizeCycle [] = []
sanitizeCycle l = take n $ tail $
                  sanitize (drop (n-1) $ cycle l) (fst $ head l)
  where n = length l

sanitize :: [(DPoint, MetaJoin Double)] -> DPoint -> [(DPoint, MetaJoin Double)]
sanitize [] _ = []

-- ending open => curl
sanitize [(p, MetaJoin m t1 t2 Open)] r =
  if p == r
  then [(p, Controls p p)]
  else [(p, MetaJoin m t1 t2 (Curl 1))]

sanitize ((p, MetaJoin m1 tl tr Open): rest@(node2:node3:_)) r
  | (fst node2 == fst node3) && (metaTypeL (snd node2) == Open) =
    (p, MetaJoin m1 tl tr (Curl 1)) : sanitize rest r
    
sanitize (node1@(p, MetaJoin m1 tl tr m2): node2@(q, MetaJoin n1 sl sr n2): rest) r
  | p == q =
    -- if two consecutive points are the same, just make a curve with all control points the same
    -- we still have to propagate a curl or given direction.
    let newnode = (p, Controls p p)
    in case (m2, n1) of
      (Curl g, Open) -> -- curl, open => explicit, curl
        newnode : sanitize ((q, MetaJoin (Curl g) sl sr n2):rest) r
      (Direction dir, Open) ->   -- given, open => explicit, given
        newnode : sanitize ((q, MetaJoin (Direction dir) sl sr n2) : rest) r
      (Open, Open) ->   -- open, open => explicit, curl
        newnode : sanitize ((q, MetaJoin (Curl 1) sl sr n2) : rest) r
      _ -> newnode : sanitize (node2:rest) r
  | otherwise =
    case (m2, n1) of
      (Curl g, Open) -> -- curl, open => curl, curl
        node1 : sanitize ((q, MetaJoin (Curl g) sl sr n2):rest) r
      (Open, Curl g) -> -- open, curl => curl, curl
        (p, MetaJoin m1 tl tr (Curl g)) : sanitize (node2:rest) r
      (Direction dir, Open) ->   -- given, open => given, given
        node1 : sanitize ((q, MetaJoin (Direction dir) sl sr n2) : rest) r
      (Open, Direction dir) ->   -- open, given => given, given
        (p, MetaJoin m1 tl tr (Direction dir)) : sanitize (node2:rest) r
      _ -> node1 : sanitize (node2:rest) r

sanitize ((p, m): (q, n): rest) r =
  case (m, n) of
    (Controls _u v, MetaJoin Open t1 t2 mt2) -- explicit, open => explicit, given
      | q == v    -> (p, m) : sanitize ((q, MetaJoin (Curl 1) t1 t2 mt2): rest) r
      | otherwise -> (p, m) : sanitize ((q, MetaJoin (Direction (q^-^v)) t1 t2 mt2): rest) r
    (MetaJoin mt1 tl tr Open, Controls u _v) -- open, explicit => given, explicit
      | u == p    -> (p, MetaJoin mt1 tl tr (Curl 1)) : sanitize ((q, n): rest) r 
      | otherwise -> (p, MetaJoin mt1 tl tr (Direction (u^-^p))) : sanitize ((q, n): rest) r
    _ -> (p, m) : sanitize ((q, n) : rest) r

sanitize (n:l) r = n:sanitize l r

spanList :: ([a] -> Bool) -> [a] -> ([a], [a])
spanList _ xs@[] =  (xs, xs)
spanList p xs@(x:xs')
  | p xs =  let (ys,zs) = spanList p xs' in (x:ys,zs)
  | otherwise    =  ([],xs)

-- break the subsegment if the angle to the left or the right is defined or a curl.
breakPoint :: [(DPoint, MetaJoin Double)] -> Bool
breakPoint ((_,  MetaJoin _ _ _ Open):(_, MetaJoin Open _ _ _):_) = False
breakPoint _ = True

-- solve the tridiagonal system for t[i]:
-- a[n] t[i-1] + b[n] t[i] + c[n] t[i+1] = d[i]
-- where a[0] = c[n] = 0
-- by first rewriting it into
-- the system t[i] + u[i] t[i+1] = v[i]
-- where u[n] = 0
-- then solving for t[n]
-- see metafont the program: ¶ 283
solveTriDiagonal :: [(Double, Double, Double, Double)] -> [Double]
solveTriDiagonal [] = error "solveTriDiagonal: not enough equations"
solveTriDiagonal ((_, b0, c0, d0): rows) =
  V.toList $ solveTriDiagonal2 (b0, c0, d0) (V.fromList rows)

solveTriDiagonal2 :: (Double, Double, Double) -> V.Vector (Double, Double, Double, Double) -> V.Vector Double
solveTriDiagonal2 (!b0, !c0, !d0) rows = solutions
  where
    twovars = V.scanl nextrow (c0/b0, d0/b0) rows
    solutions = V.scanr nextsol vn (V.unsafeInit twovars)
    vn = snd $ V.unsafeLast twovars
    nextsol (u, v) ti = v - u*ti
    nextrow (u, v) (ai, bi, ci, di) =
      (ci/(bi - u*ai), (di - v*ai)/(bi - u*ai))

-- test = ((80.0,58.0,51.0),[(-432.0,78.0,102.0,503.0),(71.0,-82.0,20.0,2130.0),(52.39,-10.43,4.0,56.0),(34.0,38.0,0.0,257.0)])
-- [-15.726940528143576,22.571642107784243,-78.93751365259996,-297.27313545829384,272.74438435742667]
      
     
-- solve the cyclic tridiagonal system.
-- see metafont the program: ¶ 286
solveCyclicTriD :: [(Double, Double, Double, Double)] -> [Double]
solveCyclicTriD rows = solutions
  where
    (!un, !vn, !wn): threevars =
      reverse $ tail $ scanl nextrow (0, 0, 1) rows
    nextrow (!u, !v, !w) (!ai, !bi, !ci, !di) =
      (ci/(bi - ai*u), (di - ai*v)/(bi - ai*u), -ai*w/(bi - ai*u))
    (totvn, totwn) = foldl (\(v', w') (u, v, w) ->
                             (v - u*v', w - u*w'))
                     (0, 1) threevars
    t0 = (vn - un*totvn) / (1 - (wn - un*totwn))
    solutions = scanl nextsol t0
                ((un, vn, wn) : reverse (tail threevars))
    nextsol t (!u, !v, !w) = (v + w*t0 - t)/u

turnAngle :: DPoint -> DPoint -> Double
turnAngle (Point 0 0) _ = 0
turnAngle (Point x y) q = vectorAngle $ rotateVec p $* q
  where p = Point x (-y)

zipNext :: [b] -> [(b, b)]
zipNext [] = []
zipNext l = zip l (tail $ cycle l)

-- find the equations for a cycle containing only open points
eqsCycle :: [Tension] -> [DPoint] -> [Tension]
         -> [Double] -> [(Double, Double, Double, Double)]
eqsCycle tensionsA points tensionsB turnAngles = 
  zipWith4 eqTension
  (zipNext (map tensionValue tensionsA))
  (zipNext dists)
  (zipNext turnAngles)
  (zipNext (map tensionValue tensionsB))
  where 
    dists = zipWith vectorDistance points (tail $ cycle points)

-- find the equations for an path with open points.
-- The first and last node should be a curl or a given angle

eqsOpen :: [DPoint] -> [MetaJoin Double] -> [DPoint] -> [Double]
        -> [Double] -> [Double] -> [(Double, Double, Double, Double)]
eqsOpen _ [MetaJoin mt1 t1 t2 mt2] [delta] _ _ _ =
  let replaceType Open = Curl 1
      replaceType t = t
  in case (replaceType mt1, replaceType mt2) of
    (Curl g, Direction dir) ->
      [eqCurl0 g (tensionValue t1) (tensionValue t2) 0,
       (0, 1, 0, turnAngle delta dir)]
    (Direction dir, Curl g) ->
      [(0, 1, 0, turnAngle delta dir),
       eqCurlN g (tensionValue t1) (tensionValue t2)]
    (Direction dir, Direction dir2) ->
      [(0, 1, 0, turnAngle delta dir),
       (0, 1, 0, turnAngle delta dir2)]
    (Curl _, Curl _) ->
      [(0, 1, 0, 0), (0, 1, 0, 0)]
    _ -> undefined

eqsOpen points joins chords turnAngles tensionsA tensionsB =
  eq0 : restEquations joins tensionsA dists turnAngles tensionsB
  where
    dists = zipWith vectorDistance points (tail points)      
    eq0 = case head joins of
      (MetaJoin (Curl g) _ _ _) -> eqCurl0 g (head tensionsA) (head tensionsB) (head turnAngles)
      (MetaJoin (Direction dir) _ _ _) -> (0, 1, 0, turnAngle (head chords) dir)
      (MetaJoin Open _ _ _) -> eqCurl0 1 (head tensionsA) (head tensionsB) (head turnAngles)
      (Controls _ _) -> error "eqsOpen: illegal join"

    restEquations [lastnode] (tensionA:_) _ _ (tensionB:_) =
      case lastnode of
        MetaJoin _ _ _ (Curl g) -> [eqCurlN g tensionA tensionB]
        MetaJoin _ _ _ Open -> [eqCurlN 1 tensionA tensionB]
        MetaJoin _ _ _ (Direction dir) -> [(0, 1, 0, turnAngle (last chords) dir)]
        (Controls _ _) -> error "eqsOpen: illegal join"

    restEquations (_:othernodes) (tensionA:restTA) (d:restD) (turn:restTurn) (tensionB:restTB) =
      eqTension (tensionA, head restTA) (d, head restD) (turn, head restTurn) (tensionB, head restTB) :
      restEquations othernodes restTA restD restTurn restTB

    restEquations _ _ _ _ _ = error "eqsOpen: illegal rest equations"

-- the equation for an open node
eqTension :: (Double, Double) -> (Double, Double)
          -> (Double, Double) -> (Double, Double)
          -> (Double, Double, Double, Double)
eqTension (tensionA', tensionA) (dist', dist) (psi', psi) (tensionB', tensionB) =
  (a, b+c, d, -b*psi' - d*psi)
  where
    a = tensionB' * tensionB' / (tensionA' * dist')
    b = (3 - 1/tensionA') * tensionB' * tensionB' / dist'
    c = (3 - 1/tensionB) * tensionA * tensionA / dist
    d = tensionA * tensionA / (tensionB * dist)

-- the equation for a starting curl
eqCurl0 :: Double -> Double -> Double -> Double -> (Double, Double, Double, Double)
eqCurl0 gamma tensionA tensionB psi = (0, c, d, r)
  where
    c = chi/tensionA + 3 - 1/tensionB
    d = (3 - 1/tensionA)*chi + 1/tensionB
    chi = gamma*tensionB*tensionB / (tensionA*tensionA)
    r = -d*psi

-- the equation for an ending curl
eqCurlN :: Double -> Double -> Double -> (Double, Double, Double, Double)
eqCurlN gamma tensionA tensionB = (a, b, 0, 0)
  where
    a = (3 - 1/tensionB)*chi + 1/tensionA
    b = chi/tensionB + 3 - 1/tensionA
    chi = gamma*tensionA*tensionA / (tensionB*tensionB)

-- getting the control points
unmetaJoin :: DPoint -> DPoint -> Double -> Double -> Tension -> Tension -> PathJoin Double
unmetaJoin !z0 !z1 !theta !phi !alpha !beta
  | abs phi < 1e-4 && abs theta < 1e-4 = JoinLine z1
  | otherwise = JoinCurve u v
  where Point dx dy = z1^-^z0
        bounded = (sf <= 0 && st <= 0 && sf <= 0) ||
                  (sf >= 0 && st >= 0 && sf >= 0)
        rr' = velocity st sf ct cf alpha
        ss' = velocity sf st cf ct beta
        stf = st*cf + sf*ct -- sin (theta + phi)
        st = sin theta
        sf = sin phi
        ct = cos theta
        cf = cos phi
        rr = case alpha of
          TensionAtLeast _ | bounded ->
            min rr' (sf/stf)
          _ -> rr'
        ss = case beta of
          TensionAtLeast _ | bounded ->
            min ss' (st/stf)
          _ -> ss'
        -- u = z0 + rr * (rotate theta chord)
        u = z0 ^+^ rr *^ Point (dx*ct - dy*st) (dy*ct + dx*st)
        -- v = z1 - ss * (rotate (-phi) chord)
        v = z1 ^-^ ss *^ Point (dx*cf + dy*sf) (dy*cf - dx*sf)

constant1, constant2, sqrt2 :: Double
constant1 = 3*(sqrt 5 - 1)/2
constant2 = 3*(3 - sqrt 5)/2
sqrt2 = sqrt 2

-- another magic formula by John Hobby.
velocity :: Double -> Double -> Double
         -> Double -> Tension -> Double
velocity st sf ct cf t =
  min 4 $ 
  (2 + sqrt2 * (st - sf/16)*(sf - st/16)*(ct - cf)) /
  ((3 + constant1*ct + constant2*cf) * tensionValue t)
