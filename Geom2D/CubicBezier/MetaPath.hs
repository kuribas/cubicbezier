{-# LANGUAGE DeriveFoldable #-}
{-# LANGUAGE DeriveTraversable #-}
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
       (unmetaOpen, unmetaClosed, ClosedMetaPath(..), OpenMetaPath (..),
        MetaJoin (..), MetaNodeType (..), Tension (..))
where
import Geom2D
import Geom2D.CubicBezier.Basic
import Text.Printf
import Data.Monoid
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector as BV
import Geom2D.CubicBezier.Numeric

data OpenMetaPath a = OpenMetaPath [(Point a, MetaJoin a)] (Point a)
                        -- ^ A metapath with endpoints
                    deriving (Functor, Traversable, Foldable)


data ClosedMetaPath a = ClosedMetaPath [(Point a, MetaJoin a)]
                        -- ^ A metapath with cycles.  The last join
                        -- joins the last point with the first.
                      deriving (Eq, Functor, Traversable, Foldable)

data MetaJoin a = MetaJoin { metaTypeL :: MetaNodeType a
                           -- ^ The nodetype going out of the
                           -- previous point.  The metafont default is
                           -- @Open@.
                           , tensionL :: Tension a
                             -- ^ The tension going out of the previous point.
                             -- The metafont default is 1.
                           , tensionR :: Tension a
                             -- ^ The tension going into the next point.
                             -- The metafont default is 1.
                           , metaTypeR :: MetaNodeType a
                             -- ^ The nodetype going into the next point.
                             -- The metafont default is @Open@.
                           }
                | Controls (Point a) (Point a)
                  -- ^ Specify the control points explicitly.
                deriving (Show, Eq, Functor, Traversable, Foldable)
                         
data MetaNodeType a = Open
                    -- ^ An open node has no direction specified.  If
                    -- it is an internal node, the curve will keep the
                    -- same direction going into and going out from
                    -- the node.  If it is an endpoint or corner
                    -- point, it will have curl of 1.
                  | Curl {curlgamma :: a}
                    -- ^ The node becomes and endpoint or a corner
                    -- point.  The curl specifies how much the segment
                    -- `curves`.  A curl of `gamma` means that the
                    -- curvature is `gamma` times that of the
                    -- following node.
                  | Direction {nodedir :: Point a}
                    -- ^ The node has a given direction.
                  deriving (Eq, Show, Functor, Foldable, Traversable)

data Tension a = Tension {tensionValue :: a}
               -- ^ The tension value specifies how /tense/ the curve is.
               -- A higher value means the curve approaches a line
               -- segment, while a lower value means the curve is more
               -- round.  Metafont doesn't allow values below 3/4.
             | TensionAtLeast {tensionValue :: a}
               -- ^ Like @Tension@, but keep the segment inside the
               -- bounding triangle defined by the control points, if
               -- there is one.
             deriving (Eq, Show, Functor, Foldable, Traversable)

instance (Show a, Real a) => Show (ClosedMetaPath a) where
  show (ClosedMetaPath nodes) =
    showPath nodes ++ "cycle"

instance (Show a, Real a) => Show (OpenMetaPath a) where
  show (OpenMetaPath nodes lastpoint) =
    showPath nodes ++ showPoint lastpoint

showPath :: (Show a, Real a) => [(Point a, MetaJoin a)] -> String
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
    showTension (TensionAtLeast t) = printf "atleast %.3f" (realToFrac t :: Double) :: String
    showTension (Tension t) = printf "%.3f" (realToFrac t :: Double) :: String
    typename Open = ""
    typename (Curl g) = printf "{curl %.3f}" (realToFrac g :: Double) :: String
    typename (Direction dir) = printf "{%s}" (showPoint dir) :: String

zipWithNext :: (V.Unbox a, V.Unbox b) => (a -> a -> b) -> V.Vector a -> V.Vector a -> V.Vector b
zipWithNext f v1 v2 = V.generate (V.length v1) val
  where val i | i < V.length v1-1 = f (v1 V.! i) (v2 V.! (i+1))
              | otherwise = f (v1 V.! i) (v2 V.! 0)

zipWithNext4 :: V.Unbox b => (Double -> Double -> Double -> Double ->
                 Double -> Double -> Double -> Double -> b)
             -> V.Vector Double -> V.Vector Double
             -> V.Vector Double -> V.Vector Double
             -> V.Vector b
zipWithNext4 f v1 v2 v3 v4 = V.generate (V.length v1) val
  where val i | i < V.length v1-1 =
                  f (v1 V.! i) (v1 V.! (i+1))
                  (v2 V.! i) (v2 V.! (i+1))
                  (v3 V.! i) (v3 V.! (i+1))
                  (v4 V.! i) (v4 V.! (i+1) )
              | otherwise =
                  f (v1 V.! i) (v1 V.! 0)
                  (v2 V.! i) (v2 V.! 0)
                  (v3 V.! i) (v3 V.! 0)
                  (v4 V.! i) (v4 V.! 0) 
                            
showPoint :: Show a => Point a -> String
showPoint (Point x y) = "(" ++ show x ++ ", " ++ show y ++ ")"

-- | Create a normal path from a metapath.
unmetaOpen :: OpenMetaPath Double -> Path Open Double
unmetaOpen (OpenMetaPath nodes endpoint) =
  unmetaOpen'
  (flip sanitize endpoint $ removeEmptyDirs nodes)
  endpoint

unmetaOpen' :: [(DPoint, MetaJoin Double)] -> DPoint -> Path Open Double
unmetaOpen' nodes endpoint =
  joinSegments $ map unmetaSubSegment subsegs
  where
    subsegs = openSubSegments nodes endpoint
    
unmetaClosed :: ClosedMetaPath Double -> Path Closed Double
unmetaClosed (ClosedMetaPath nodes) =
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
unmetaAsOpen :: [(DPoint, MetaJoin Double)] -> [(DPoint, MetaJoin Double)] -> Path Closed Double
unmetaAsOpen l [] = closeOpenPath $ unmetaOpen' (sanitizeCycle l) (fst $ head l)
unmetaAsOpen l m = closeOpenPath $ Path l'' <> Path m''
  where n = length m
        Path o = unmetaOpen' (sanitizeCycle (m++l)) (fst $ head (m ++l))
        (m',ct:l') = splitAt n o
        (l'', m'') =
          case ct of
            CurveTo _ _ p -> (LineTo p: l', m'++[ct])
            LineTo _ -> (ct: l', m')
               

-- decompose into a list of subsegments that need to be solved.
openSubSegments :: [(DPoint, MetaJoin Double)] -> DPoint -> [OpenMetaPath Double]
openSubSegments [] _ =  []
openSubSegments l lastPoint = 
  case spanList (not . breakPoint) l of
    (m, n:o) ->
      let point = case o of
            ((p,_):_) -> p
            _ -> lastPoint
      in OpenMetaPath (m ++ [n]) point :
         openSubSegments o lastPoint
    _ -> error "openSubSegments': unexpected end of segments"

spanList :: ([a] -> Bool) -> [a] -> ([a], [a])
spanList _ xs@[] =  (xs, xs)
spanList p xs@(x:xs')
  | p xs =  let (ys,zs) = spanList p xs' in (x:ys,zs)
  | otherwise    =  ([],xs)

-- break the subsegment if the angle to the left or the right is defined or a curl.
-- break the subsegment if the angle to the left or the right is defined or a curl.
breakPoint :: [(DPoint, MetaJoin Double)] -> Bool
breakPoint ((_,  MetaJoin _ _ _ Open):(_, MetaJoin Open _ _ _):_) = False
breakPoint _ = True

-- join subsegments into one segment
joinSegments :: [Path Open Double] -> Path Open Double
joinSegments = mconcat

-- solve a cyclic metapath where all angles depend on the each other.
unmetaCyclic :: [(DPoint, MetaJoin Double)] -> Path Closed Double
unmetaCyclic nodes =
  let points = V.fromList $ map fst nodes
      chords = zipWithNext (flip (^-^)) points points
      tensionsA = BV.fromList $ map (tensionL . snd) nodes
      tensionsB = BV.fromList $ map (tensionR . snd) nodes
      turnAngles = zipWithNext turnAngle chords chords
      thetas = solveCyclicTriD $
               eqsCycle tensionsA points tensionsB turnAngles
      phis = zipWithNext (\x y -> -(x+y)) turnAngles thetas
  in Path $ LineTo (V.head points) :
     BV.toList (BV.zipWith6 unmetaJoin (V.convert points)
                (V.convert $ V.tail points V.++ V.singleton (V.head points))
     (V.convert thetas) (V.convert phis) tensionsA tensionsB)

-- solve a subsegment
unmetaSubSegment :: OpenMetaPath Double -> Path Open Double

-- the simple case where the control points are already given.
unmetaSubSegment (OpenMetaPath [(p, Controls u v)] q) =
  Path [LineTo p, CurveTo u v q]

-- otherwise solve the angles, and find the control points
unmetaSubSegment (OpenMetaPath nodes lastpoint) =
  let points = V.fromList $ map fst nodes ++ [lastpoint]
      joins = BV.fromList $ map snd nodes
      chords = V.zipWith (^-^) (V.tail points) points
      tensionsA = BV.map tensionL joins
      tensionsB = BV.map tensionR joins
      turnAngles = V.zipWith turnAngle chords (V.tail chords) V.++ V.singleton 0
      thetas = solveTriDiagonal2 $
               eqsOpen points joins chords turnAngles
               (V.convert $ BV.map tensionValue tensionsA)
               (V.convert $ BV.map tensionValue tensionsB)
      phis = V.zipWith (\x y -> -(x+y)) turnAngles (V.tail thetas)
      pathjoins =
        BV.zipWith6 unmetaJoin (V.convert points) (BV.tail $ V.convert points)
        (V.convert thetas) (V.convert phis) tensionsA tensionsB
  in Path $ LineTo (V.head points) : BV.toList pathjoins

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

-- solve the tridiagonal system for t[i]:
-- a[n] t[i-1] + b[n] t[i] + c[n] t[i+1] = d[i]
-- where a[0] = c[n] = 0
-- by first rewriting it into
-- the system t[i] + u[i] t[i+1] = v[i]
-- where u[n] = 0
-- then solving for t[n]
-- see metafont the program: ¶ 283
solveTriDiagonal2 :: V.Vector (Double, Double, Double, Double) -> V.Vector Double
solveTriDiagonal2 v
  | V.length v == 0 = error "solveTriDiagonal: not enough equations"
  | otherwise = solveTriDiagonal (b0, c0, d0) rows
    where (_, b0, c0, d0) = V.head v
          rows = V.tail v

-- test = ((80.0,58.0,51.0),[(-432.0,78.0,102.0,503.0),(71.0,-82.0,20.0,2130.0),(52.39,-10.43,4.0,56.0),(34.0,38.0,0.0,257.0)])
-- [-15.726940528143576,22.571642107784243,-78.93751365259996,-297.27313545829384,272.74438435742667]
      
turnAngle :: DPoint -> DPoint -> Double
turnAngle (Point 0 0) _ = 0
turnAngle (Point x y) q = vectorAngle $ rotateVec p $* q
  where p = Point x (-y)

-- find the equations for a cycle containing only open points
eqsCycle :: BV.Vector (Tension Double) -> V.Vector DPoint
         -> BV.Vector (Tension Double) -> V.Vector Double
         -> V.Vector (Double, Double, Double, Double)
eqsCycle tensionsA points tensionsB turnAngles = 
  zipWithNext4 eqTension
  (V.convert $ BV.map tensionValue tensionsA)
  dists turnAngles
  (V.convert $ BV.map tensionValue tensionsB)
  where 
    dists = V.zipWith vectorDistance points (V.tail points V.++ points)

-- find the equations for an path with open points.
-- The first and last node should be a curl or a given angle

replaceOpen :: Num a => MetaNodeType a -> MetaNodeType a
replaceOpen Open = Curl 1
replaceOpen t = t

eqsOpen :: V.Vector DPoint -> BV.Vector (MetaJoin Double) -> V.Vector DPoint
        -> V.Vector Double -> V.Vector Double -> V.Vector Double
        -> V.Vector (Double, Double, Double, Double)
eqsOpen points joins chords turnAngles tensionsA tensionsB
  | BV.length joins == 0 = error "empty segment"
  | BV.length joins == 1 =
      let MetaJoin mt1 t1 t2 mt2 = BV.head joins
          delta = V.head chords
      in case (replaceOpen mt1, replaceOpen mt2) of
        (Curl g, Direction dir) ->
          V.fromList [
          (eqCurl0 g (tensionValue t1) (tensionValue t2) 0),
          (0, 1, 0, turnAngle delta dir)]
        (Direction dir, Curl g) ->
          V.fromList [
          (0, 1, 0, turnAngle delta dir),
          eqCurlN g (tensionValue t1) (tensionValue t2)]
        (Direction dir, Direction dir2) ->
          V.fromList [
          (0, 1, 0, turnAngle delta dir),
          (0, 1, 0, turnAngle delta dir2)]
        (Curl _, Curl _) ->
          V.fromList [
          (0, 1, 0, 0),
          (0, 1, 0, 0)]
        _ -> error "Illegal open path"
  | otherwise = V.generate (BV.length joins+1) eqs
  where
    dists = V.zipWith vectorDistance points (V.tail points)      
    eqs i | i == 0 = case BV.head joins of
              (MetaJoin (Curl g) _ _ _) ->
                eqCurl0 g (V.head tensionsA) (V.head tensionsB) (V.head turnAngles)
              (MetaJoin (Direction dir) _ _ _) ->
                (0, 1, 0, turnAngle (V.head chords) dir)
              (MetaJoin Open _ _ _) ->
                eqCurl0 1 (V.head tensionsA) (V.head tensionsB) (V.head turnAngles)
              _ -> error "eqsOpen: illegal join"
          | i < BV.length joins =
              eqTension (tensionsA V.! (i-1)) (tensionsA V.! i) (dists V.! (i-1)) (dists V.! i)
              (turnAngles V.! (i-1)) (turnAngles V.! i) (tensionsB V.! (i-1)) (tensionsB V.! i)
          | otherwise = case joins BV.! (i-1) of
              MetaJoin _ _ _ (Curl g) ->
                eqCurlN g (tensionsA V.! (i-1)) (tensionsB V.! (i-1))
              MetaJoin _ _ _ (Direction dir) ->
                (0, 1, 0, turnAngle (V.last chords) dir)
              MetaJoin _ _ _ Open ->
                eqCurlN 1 (tensionsA V.! (i-1)) (tensionsB V.! (i-1))
              _ -> error "eqsOpen: illegal join"

-- the equation for an open node
eqTension :: Double -> Double -> Double -> Double 
          -> Double -> Double -> Double -> Double
          -> (Double, Double, Double, Double)
eqTension tensionA' tensionA dist' dist psi' psi tensionB' tensionB =
  (a, b+c, d, -b*psi' - d*psi)
  where
    a = tensionB' * tensionB' / (tensionA' * dist')
    b = (3 - 1/tensionA') * tensionB' * tensionB' / dist'
    c = (3 - 1/tensionB) * tensionA * tensionA / dist
    d = tensionA * tensionA / (tensionB * dist)

-- the equation for a starting curl
eqCurl0 :: Double -> Double -> Double -> Double
        -> (Double, Double, Double, Double)
eqCurl0 gamma tensionA tensionB psi = (0, c, d, r)
  where
    c = chi/tensionA + 3 - 1/tensionB
    d = (3 - 1/tensionA)*chi + 1/tensionB
    chi = gamma*tensionB*tensionB / (tensionA*tensionA)
    r = -d*psi

-- the equation for an ending curl
eqCurlN :: Double -> Double -> Double
        -> (Double, Double, Double, Double)
eqCurlN gamma tensionA tensionB = (a, b, 0, 0)
  where
    a = (3 - 1/tensionB)*chi + 1/tensionA
    b = chi/tensionB + 3 - 1/tensionA
    chi = gamma*tensionA*tensionA / (tensionB*tensionB)

-- getting the control points
unmetaJoin :: DPoint -> DPoint -> Double -> Double -> Tension Double
           -> Tension Double -> PathJoin Double
unmetaJoin !z0 !z1 !theta !phi !alpha !beta
  | abs phi < 1e-4 && abs theta < 1e-4 = LineTo z1
  | otherwise = CurveTo u v z1
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
         -> Double -> Tension Double -> Double
velocity st sf ct cf t =
  min 4 $ 
  (2 + sqrt2 * (st - sf/16)*(sf - st/16)*(ct - cf)) /
  ((3 + constant1*ct + constant2*cf) * tensionValue t)
