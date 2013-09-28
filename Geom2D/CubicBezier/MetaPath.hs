{-# LANGUAGE BangPatterns #-}
module Geom2D.CubicBezier.MetaPath
where
import Geom2D.CubicBezier
import Data.List
import Text.Printf
import Debug.Trace

data MetaPath = OpenMetaPath [(Point, MetaJoin)] Point
              | CyclicMetaPath [(Point, MetaJoin)]

data MetaJoin = Implicit {metaType1 :: MetaNodeType,
                          metaType2 :: MetaNodeType}
              | Explicit Point Point
              deriving Show

data MetaNodeType = Open {getTension :: Tension }
                  | Curl {curlgamma :: Double, getTension :: Tension}
                  | Given {nodedir :: Point, getTension :: Tension}
                  deriving Show

data Tension = Tension {tensionValue :: Double}
             | TensionAtLeast {tensionValue :: Double}
             deriving (Eq, Show)

instance Show MetaPath where
  show (CyclicMetaPath nodes) =
    showPath nodes ++ "cyclic"
  show (OpenMetaPath nodes lastpoint) =
    showPath nodes ++ showPoint lastpoint

showPath :: [(Point, MetaJoin)] -> [Char]
showPath = concatMap showNodes
  where
    showNodes (p, Explicit u v) =
      showPoint p ++ "..controls " ++ showPoint u ++ "and " ++ showPoint v ++ ".."
    showNodes (p, Implicit m1 m2) =
      showPoint p ++ typename m1 ++ ".." ++ tensions ++ typename m2
      where
        tensions
          | getTension m1 == getTension m2 && getTension m1 == Tension 1 = ""
          | getTension m1 == getTension m2 = printf "tension %s .." (showTension $ getTension m1)
          | otherwise = printf "tensions %s and %s.."
                        (showTension $ getTension m1) (showTension $ getTension m2)
    showTension (TensionAtLeast t) = printf "atleast %.3f" t :: String
    showTension (Tension t) = printf "%.3f" t :: String
    typename (Open _) = ""
    typename (Curl g _) = printf "{Curl %.3f}" g :: String
    typename (Given dir _) = printf "{%s}" (showPoint dir) :: String
    
showPoint :: Point -> String
showPoint (Point x y) = printf "(%.3f, %.3f)" x y

-- | Create a normal path from a metapath.
unmeta :: MetaPath -> Path
unmeta (OpenMetaPath nodes endpoint) =
  let subsegs = openSubSegments (sanitizeOpen nodes) endpoint
      path = joinSegments $ map unmetaSubSegment subsegs
  in OpenPath path endpoint

unmeta (CyclicMetaPath nodes) =
  case span (bothOpen . snd) nodes of
    (l, []) -> unmetaCyclic l
    (l, (m:n)) ->
      if leftOpen $ snd m
      then unmetaAsOpen (l++[m]) n
      else unmetaAsOpen l (m:n)

-- decompose into a list of subsegments that need to be solved.
openSubSegments :: [(Point, MetaJoin)] -> Point -> [MetaPath]
openSubSegments l p = openSubSegments' (tails l) p

openSubSegments' :: [[(Point, MetaJoin)]] -> Point -> [MetaPath]
openSubSegments' [[]] _ = []
openSubSegments' [] _   = []
openSubSegments' l lastPoint = case break breakPoint l of
  (m, n:o) ->
    let point = case o of
          (((p,_):_):_) -> p
          _ -> lastPoint
    in OpenMetaPath (map head (m ++ [n])) point :
       openSubSegments' o lastPoint
  _ -> error "unexpected end of segments"

-- join subsegments into one segment
joinSegments :: [Path] -> [(Point, PathJoin)]
joinSegments = concatMap nodes
  where nodes (OpenPath n _) = n
        nodes (ClosedPath n) = n

-- solve a cyclic metapath where all angles depend on the each other.
unmetaCyclic :: [(Point, MetaJoin)] -> Path
unmetaCyclic nodes =
  let points = map fst nodes
      chords = zipWith (^-^) points (last points : points)
      tensionsA = (map (getTension . metaType1 . snd) nodes)
      tensionsB = (map (getTension . metaType2 . snd) nodes)
      turnAngles = zipWith turnAngle chords (tail $ cycle chords)
      thetas = solveCyclicTriD $
               eqsCycle tensionsA
               points
               tensionsB
               turnAngles
      phis = zipWith (\x y -> -(x+y)) turnAngles (tail thetas)
  in ClosedPath $ zip points $
     zipWith6 unmetaJoin points (tail points)
     thetas phis tensionsA tensionsB

-- solve a cyclic metapath as an open path if possible.
-- rotate to the defined node, and rotate back after
-- solving the path.
unmetaAsOpen :: [(Point, MetaJoin)] -> [(Point, MetaJoin)] -> Path
unmetaAsOpen l m = ClosedPath (l'++m') 
  where n = length m
        OpenPath o _ = unmeta $ OpenMetaPath (m++l) (fst $ head m)
        (m',l') = splitAt n o

-- solve a subsegment
unmetaSubSegment :: MetaPath -> Path

-- the simple case where the control points are already given.
unmetaSubSegment (OpenMetaPath [(p, Explicit u v)] q) =
  OpenPath [(p, JoinCurve u v)] q

-- otherwise solve the angles, and find the control points
unmetaSubSegment mp@(OpenMetaPath nodes lastpoint) =
  let points = map fst nodes ++ [lastpoint]
      joins = map snd nodes
      chords = zipWith (^-^) (tail points) points
      tensionsA = map (getTension . metaType1) joins
      tensionsB = map (getTension . metaType2) joins
      turnAngles = zipWith turnAngle chords (tail chords) ++ [0]
      thetas = solveTriDiagonal $
               eqsOpen points joins chords turnAngles
               (map tensionValue tensionsA)
               (map tensionValue tensionsB)
      phis = zipWith (\x y -> -x-y) turnAngles (tail thetas)
      pathjoins = zipWith6 unmetaJoin points (tail points) thetas phis tensionsA tensionsB
  in OpenPath (zip points pathjoins) lastpoint

unmetaTest mp@(OpenMetaPath nodes' lastpoint) =
  let nodes = sanitizeOpen nodes'
      points = map fst nodes ++ [lastpoint]
      joins = map snd nodes
      chords = zipWith (^-^) (tail points) points
      tensionsA = map (getTension . metaType1) joins
      tensionsB = map (getTension . metaType2) joins
      turnAngles = zipWith turnAngle chords (tail chords) ++ [0]
      eqs = --solveTriDiagonal $
               eqsOpen points joins chords turnAngles
               (map tensionValue tensionsA)
               (map tensionValue tensionsB)
      --pathjoins = zipWith6 unmetaJoin points (tail points) thetas turnAngles tensionsA tensionsB
  in eqs


bothOpen :: MetaJoin -> Bool
bothOpen (Implicit (Open _) (Open _)) = True
bothOpen _ = False

leftOpen :: MetaJoin -> Bool
leftOpen (Implicit (Open _) _) = True
leftOpen _ = False

-- replace open nodetypes with more defined nodetypes if possible
sanitizeOpen :: [(Point, MetaJoin)] -> [(Point, MetaJoin)]
sanitizeOpen [] = []

-- starting open => curl
sanitizeOpen ((p, Implicit (Open t) m):rest) =
  sanitizeRest ((p, Implicit (Curl 1 t) m):rest)
sanitizeOpen l = sanitizeRest l
   
sanitizeRest :: [(Point, MetaJoin)] -> [(Point, MetaJoin)]
sanitizeRest [] = []

-- ending open => curl
sanitizeRest [(p, Implicit m (Open t))] =
  [(p, Implicit m (Curl 1 t))]

sanitizeRest (node1@(p, Implicit m1 m2): node2@(q, Implicit n1 n2): rest) =
  case (m2, n1) of
    (Curl _ _, Open t) -> -- curl, open => curl, curl
      node1 : sanitizeRest ((q, Implicit (Curl 1 t) n2):rest)
    (Open t, Curl _ _) -> -- open, curl => curl, curl
      (p, Implicit m1 (Curl 1 t)) : sanitizeRest (node2:rest)
    (Given dir _, Open t) ->   -- given, open => given, given
      node1 : sanitizeRest ((q, (Implicit (Given dir t) n2)) : rest)
    (Open t, Given dir _) ->   -- open, given => given, given
      (p, Implicit m1 (Given dir t)) : sanitizeRest (node2:rest)
    _ -> node1 : sanitizeRest (node2:rest)

sanitizeRest ((p, m): (q, n): rest) =
  case (m, n) of
    (Explicit _u v, Implicit (Open t) mt2) ->  -- explicit, open => explicit, given
      (p, m) : sanitizeRest ((q, Implicit (Given (q^-^v) t) mt2): rest)
    (Implicit mt1 (Open t), Explicit u _v) ->  -- open, explicit => given, explicit
      (p, Implicit mt1 (Given (u^-^p) t)) : sanitizeRest ((q, n): rest)
    _ -> (p, m) : sanitizeRest ((q, n) : rest)

sanitizeRest (n:l) = n:sanitizeRest l

-- break the subsegment if the angle to the left or the right is defined or a curl.
breakPoint :: [(Point, MetaJoin)] -> Bool
breakPoint ((_,  Implicit _ (Open _)):(_, Implicit (Open _) _):_) = False
breakPoint _ = True

-- solve the tridiagonal system for t[i]:
-- a[n] t[i-1] + b[i] t[i] + c[b] t[i+1] = d[i]
-- where a[0] = c[n] = 0
-- by first rewriting it into
-- the system t[i] + u[i] t[i+1] = v[i]
-- where u[n] = 0
-- then solving for t[n]
-- see metafont the program: ¶ 283
solveTriDiagonal :: [(Double, Double, Double, Double)] -> [Double]
solveTriDiagonal [] = error "tridiagonal: not enough equations"
solveTriDiagonal ((_, b0, c0, d0): rows) = solutions
  where
    ((_, vn): twovars) =
      reverse $ scanl nextrow (c0/b0, d0/b0) rows
    nextrow (u, v) (ai, bi, ci, di) =
      (ci/(bi - u*ai), (di - v*ai)/(bi - u*ai))
    solutions = reverse $ scanl nextsol vn twovars
    nextsol ti (u, v) = v - u*ti

solveTriDiagonal2 ((_, b0, c0, d0): rows) = reverse $ scanl nextrow (c0/b0, d0/b0) rows
  where
    nextrow (u, v) (ai, bi, ci, di) =
      (ci/(bi - u*ai), (di - v*ai)/(bi - u*ai))
    

-- test = ((80.0,58.0,51.0),[(-432.0,78.0,102.0,503.0),(71.0,-82.0,20.0,2130.0),(52.39,-10.43,4.0,56.0),(34.0,38.0,0.0,257.0)])

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

turnAngle :: Point -> Point -> Double
turnAngle (Point x y) q = vectorAngle $ rotateVec p $* q
  where p = Point x (-y)

zipPrev :: [a] -> [(a, a)]
zipPrev [] = []
zipPrev l = zip (last l : l) l

-- find the equations for a cycle containing only open points
eqsCycle tensionsA points tensionsB turnAngles = 
  zipWith4 eqTension
  (zipPrev (map tensionValue tensionsA))
  (zipPrev dists)
  (zipPrev turnAngles)
  (zipPrev (map tensionValue tensionsB))
  where 
    dists = zipWith vectorDistance points (tail $ cycle points)

-- find the equations for an path with open points.
-- The first and last node should be a curl or a given angle

eqsOpen [p, q] [join] [delta] _ _ _ =
  case join of
    Implicit (Curl _ _) (Curl _ _) ->
      [(0, 1, 0, 0), (0, 1, 0, 0)]
    Implicit (Curl g t) (Given dir t2) ->
      [eqCurl0 g (tensionValue t) (tensionValue t2) 0,
       (0, 1, 0, turnAngle delta dir)]
    Implicit (Given dir t) (Curl g t2) ->
      [(0, 1, 0, turnAngle delta dir),
       eqCurlN g (tensionValue t) (tensionValue t2)]
    _ -> error "illegal nodetype in subsegment"

eqsOpen points joins chords turnAngles tensionsA tensionsB =
  eq0 : restEquations joins tensionsA dists turnAngles tensionsB
  where
    dists = zipWith vectorDistance points (tail points)      
    eq0 = case head joins of
      (Implicit (Curl g _) _) -> eqCurl0 g (head tensionsA) (head tensionsB) (head turnAngles)
      (Implicit (Given dir _) _) -> (0, 1, 0, turnAngle (head chords) dir)
      _ -> error "illegal subsegment first nodetype"

    restEquations [lastnode] (tensionA:_) _ _ (tensionB:_) =
      case lastnode of
        Implicit _ (Curl g _) -> [eqCurlN g tensionA tensionB]
        Implicit _ (Given dir _) -> [(0, 1, 0, turnAngle (last chords) dir)]
        _  -> error "illegal subsegment last nodetype"

    restEquations (_:othernodes) (tensionA:restTA) (d:restD) (turn:restTurn) (tensionB:restTB) =
      eqTension (tensionA, head restTA) (d, head restD) (turn, head restTurn) (tensionB, head restTB) :
      restEquations othernodes restTA restD restTurn restTB

    restEquations _ _ _ _ _ = error "illegal rest equations"

-- the equation for an open node
eqTension :: (Double, Double) -> (Double, Double)
          -> (Double, Double) -> (Double, Double)
          -> (Double, Double, Double, Double)
eqTension (tensionA', tensionA) (dist', dist) (psi', psi) (tensionB', tensionB) =
  (a, b+c, d, -b*psi' - d*psi)
  where
    a = (tensionB' * tensionB' / (tensionA' * dist'))
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

tension t u = Implicit (Open (Tension t)) (Open (Tension u))
open = Implicit (Open (Tension 1)) (Open (Tension 1))

-- the equation for an ending curl
eqCurlN :: Double -> Double -> Double -> (Double, Double, Double, Double)
eqCurlN gamma tensionA tensionB = (a, b, 0, 0)
  where
    a = (3 - 1/tensionB)*chi + 1/tensionA
    b = chi/tensionB + 3 - 1/tensionA
    chi = gamma*tensionA*tensionA / (tensionB*tensionB)

-- magic formula for getting the control points by John Hobby
unmetaJoin !z0 !z1 !theta !phi !alpha !beta
  --- | abs (phi - theta) < 1e-5 = JoinLine
  | otherwise = JoinCurve u v
  where chord@(Point dx dy) = z1^-^z0
        bounded = (signum sf == signum st) &&
                  (signum st == signum stf)
        rr' = velocity st sf ct cf alpha
        ss' = velocity sf st cf ct beta
        stf = st*cf + sf*ct -- sin (theta + phi)
        st = sin theta
        sf = sin phi
        ct = cos theta
        cf = cos phi
        rr = case alpha of
          TensionAtLeast _ | bounded ->
            min rr' (st/stf)
          _ -> rr'
        ss = case beta of
          TensionAtLeast _ | bounded ->
            min ss' (sf/stf)
          _ -> ss'
        u = z0 ^+^ rr *^ Point (dx*ct - dy*st) (dy*ct + dx*st)  -- z0 + rr * (rotate theta chord)
        v = z1 ^-^ ss *^ Point (dx*cf + dy*sf) (dy*cf - dx*sf)  -- z1 - ss * (rotate (-phi) chord)

constant1, constant2, sqrt2 :: Double
constant1 = 3/2*(sqrt 5 - 1)
constant2 = 3/2*(3 - sqrt 5)
sqrt2 = sqrt 2

-- another magic formula by John Hobby.
velocity :: Double -> Double -> Double
         -> Double -> Tension -> Double
velocity st sf ct cf t =
  (2 + sqrt2 * (st - sf/16)*(sf - st/16)*(ct - cf)) /
  ((3 + constant1*ct + constant2*cf) * tensionValue t)

