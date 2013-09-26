{-# LANGUAGE BangPatterns #-}
module Geom2D.CubicBezier.MetaPath
where
import Geom2D.CubicBezier
import Test.QuickCheck
import Data.List


data MetaPath = OpenMetaPath [(Point, MetaJoin)] Point
              | CyclicMetaPath [(Point, MetaJoin)]

data MetaJoin = Implicit {metaType1 :: MetaNodeType,
                          metaType2 :: MetaNodeType}
              | Explicit Point Point

data MetaNodeType = Open {getTension :: Tension }
                  | Curl {gamma :: Double, getTension :: Tension}
                  | Given {dir :: Double, getTension :: Tension}

data Tension = Tension {tensionValue :: Double}
             | TensionAtLeast {tensionValue :: Double}

-- | Create a normal path from a metapath.
unmeta :: MetaPath -> Path
unmeta (OpenMetaPath nodes endpoint) =
  let subsegs = openSubSegments (sanitizeOpen nodes) endpoint
      path = joinSegments $ map unmetaSubsegment subsegs
  in OpenPath path endpoint

unmeta (CyclicMetaPath nodes) =
  case span (bothOpen . snd) nodes of
    (l, []) -> unmetaCyclic l
    (l, (m:n)) ->
      if leftOpen $ snd m
      then unmetaAsOpen (l++[m]) n
      else unmetaAsOpen l (m:n)

unmetaAsOpen :: [(Point, MetaJoin)] -> [(Point, MetaJoin)] -> Path
unmetaAsOpen l m = ClosedPath (l'++m') 
  where n = length m
        OpenPath o p = unmeta $ OpenMetaPath (m++l) (fst $ head m)
        (m',l') = splitAt n o
        
unmetaCyclic :: [(Point, MetaJoin)] -> Path
unmetaCyclic nodes =
  let angles = solveCyclicTriD $
               eqsCycle (map (tensionValue . getTension . metaType1 . snd) nodes)
               points
               (map (tensionValue . getTension . metaType2 . snd) nodes)
      points = map fst nodes
  in ClosedPath $ zip points $
     zipWith3 unmetaJoin points 
     (insertAngles (map snd nodes) angles)
     (tail points ++ points)

joinSegments :: [Path] -> [(Point, PathJoin)]
joinSegments = concatMap nodes
  where nodes (OpenPath n _) = n

unmetaSubsegment :: MetaPath -> Path
unmetaSubsegment (OpenMetaPath [(p, Explicit u v)] q) =
  OpenPath [(p, JoinCurve u v)] q

unmetaSubsegment mp@(OpenMetaPath nodes lastpoint) =
  let angles = solveTriDiagonal $ eqsOpen mp
      newnodes = insertAngles (map snd nodes) angles
      points = map fst nodes ++ [lastpoint]
      pathjoins = zipWith3 unmetaJoin points newnodes (tail points)
  in OpenPath (zip points pathjoins) lastpoint

bothOpen :: MetaJoin -> Bool
bothOpen (Implicit (Open _) (Open _)) = True
bothOpen _ = False

leftOpen :: MetaJoin -> Bool
leftOpen (Implicit (Open _) _) = True
leftOpen _ = False

-- replace open nodetypes if possible
sanitizeOpen :: [(Point, MetaJoin)] -> [(Point, MetaJoin)]
sanitizeOpen [] = []

-- starting open => curl
sanitizeOpen l@((p, Implicit (Open t) m):_) =
  (p, Implicit (Curl 1 t) m): sanitizeRest l
sanitizeOpen l = sanitizeRest l
   
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

sanitizeRest ((p, m): (q, n): rest) =
  case (m, n) of
    (Explicit u v, Implicit (Open t) mt2) ->  -- explicit, open => explicit, given
      (p, m) : sanitizeRest ((q, Implicit (Given dir t) mt2): rest)
      where dir = vectorAngle (q^-^v)
    (Implicit mt1 (Open t), Explicit u v) ->  -- open, explicit => given, explicit
      (p, Implicit mt1 (Given dir t)) : sanitizeRest ((q, n): rest)
      where dir = vectorAngle (u^-^p)

sanitizeRest (n:l) = n:sanitizeRest l

breakPoint :: [(Point, MetaJoin)] -> Bool
breakPoint ((_,  Implicit _ (Open _)):(_, Implicit (Open _) _):_) = False
breakPoint _ = True

-- decompose into a list of subsegments that need to be solved.
openSubSegments :: [(Point, MetaJoin)] -> Point -> [MetaPath]
openSubSegments l p = openSubSegments' (tails l) p

openSubSegments' :: [[(Point, MetaJoin)]] -> Point -> [MetaPath]
openSubSegments' [[]] _ = []
openSubSegments' [] _   = []
openSubSegments' l lastPoint = case break breakPoint l of
  (m, n:o) ->
    let point = case o of
          ([(p,_)]:_) -> p
          otherwise -> lastPoint
    in OpenMetaPath (map head (m ++ [n])) point :
       openSubSegments' o lastPoint
  otherwise -> error "unexpected end of segments"

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
turnAngle p q = vectorAngle q - vectorAngle p

zipPrev :: [a] -> [(a, a)]
zipPrev [] = []
zipPrev l = zip (last l : l) l

zipNext :: [a] -> [(a, a)]
zipNext [] = []
zipNext l = zip l (tail l)

-- neighbours [] = []
-- neighbours [a] = [(a, a, a)]
-- zipNeighbours f ls@(l : ls'@(l' : ls'')) =
--   ((last ls'), l, l') : zipWith3 f ls ls' (ls''++ls)

-- find the equations for a cycle containing only open points
eqsCycle :: [Double] -> [Point] -> [Double]
         -> [(Double, Double, Double, Double)]
eqsCycle tensionsA points tensionsB = 
  zipWith4 eqTension
  (zipPrev tensionsA)
  (zipPrev dists)
  (zipPrev turnAngles)
  (zipPrev tensionsB)
  where
    turnAngles = zipWith turnAngle chords (tail $ cycle chords)
    dists = zipWith vectorDistance points (tail $ cycle points)
    chords = zipWith (^-^) points (last points : points)

-- find the equations for a subsegment consisting of a number of open
-- points starting and ending with either a curl or a given direction.
eqsOpen :: MetaPath -> [(Double, Double, Double, Double)]
eqsOpen (OpenMetaPath [(p, Implicit n m)] q) =
  case (n, m) of
    (Curl g t, Curl g2 t2)     -> [(0, 1, 0, a), (0, 1, 0, a)]
      where a = vectorAngle (q^-^p)
    (Curl g t, Given dir t2) -> [eqCurl0 g (tensionValue t) (tensionValue t2) 0,
                                 (0, 1, 0, dir)]
    (Given dir t, Curl g t2) -> [(0, 1, 0, dir),
                                 eqCurlN g (tensionValue t) (tensionValue t2)]
    otherwise                -> error "illegal nodetype in subsegment"

eqsOpen (OpenMetaPath nodes lastpoint) =
  eq0 : restEquations joins tensionsA dists turnAngles tensionsB
  where
    points = map fst nodes ++ [lastpoint]
    joins = map snd nodes
    dists = zipWith vectorDistance chords (tail chords)
    chords = zipWith (^-^) (tail points) points
    tensionsA = map (tensionValue . getTension . metaType1) joins
    tensionsB = map (tensionValue . getTension . metaType2) joins
    turnAngles = zipWith turnAngle chords (tail chords)

    eq0 = case head joins of
      (Implicit (Curl g t) _) -> eqCurl0 g (head tensionsA) (head tensionsB) (head turnAngles)
      (Implicit (Given angle _) _) -> (0, 1, 0, angle)
      _ -> error "illegal subsegment first nodetype"

    restEquations [lastnode] (tensionA:_) _ _ (tensionB:_) =
      case lastnode of
        Implicit _ (Curl g _) -> [eqCurlN g tensionA tensionB]
        Implicit _ (Given angle _) -> [(0, 1, 0, angle)]
        otherwise -> error "illegal subsegment last nodetype"

    restEquations (node:othernodes) (tensionA:restTA) (d:restD) (turn:restTurn) (tensionB:restTB) =
      eqTension (tensionA, head restTA) (d, head restD) (turn, head restTurn) (tensionB, head restTB) :
      restEquations othernodes restTA restD restTurn restTB

eqTension :: (Double, Double) -> (Double, Double)
          -> (Double, Double) -> (Double, Double)
          -> (Double, Double, Double, Double)
eqTension (tensionA', tensionA) (d', d) (psi', psi) (tensionB', tensionB) =
  (a, b+c, d, -b*psi' - d*psi)
  where
    a = (tensionB' * tensionB' / (tensionA' * d'))
    b = (3 - 1/tensionA') * tensionB' * tensionB' / d'
    c = (3 - 1/tensionB) * tensionA * tensionA / d
    d = tensionA * tensionA / (tensionB * d)
    
eqCurl0 :: Double -> Double -> Double -> Double -> (Double, Double, Double, Double)
eqCurl0 gamma tensionA tensionB psi = (0, c, d, r)
  where
    c = chi/tensionA + 3 - 1/tensionB
    d = (3 - 1/tensionA)*chi + 1/tensionB
    chi = gamma*tensionB*tensionB / tensionA*tensionA
    r = -d*psi

eqCurlN :: Double -> Double -> Double -> (Double, Double, Double, Double)
eqCurlN gamma tensionA tensionB = (a, b, 0, 0)
  where
    a = (3 - 1/tensionB)*chi + 1/tensionA
    b = chi/tensionB + 3 - 1/tensionA
    chi = gamma*tensionA*tensionA / tensionB*tensionB

insertAngles :: [MetaJoin] -> [Double] -> [MetaJoin]
insertAngles _ [] = []
insertAngles ((Implicit m n): rest) (a:b:angles) =
  Implicit (Given a (getTension m)) (Given b (getTension n)) :
  insertAngles rest (b:angles)

-- turn a metajoin into a simple join
unmetaJoin :: Point -> MetaJoin -> Point -> PathJoin
unmetaJoin _ (Explicit u v) _  = JoinCurve u v

-- magic formula for getting the control points by John Hobby
unmetaJoin !z0 (Implicit (Given !w0 !alpha) (Given !w1 !beta)) !z1
  | abs theta < 1e-5 && abs phi < 1e-5 = JoinLine
  | otherwise = JoinCurve u v
  where chord@(Point dx dy) = z1^-^z0
        chordAngle = vectorAngle chord
        theta = chordAngle - w0
        phi   = w1 - chordAngle
        bounded = (signum sf == signum st) &&
                  (signum st == signum stf)
        rr' = velocity theta phi st sf ct cf alpha
        ss' = velocity phi theta sf st cf ct beta
        stf = st*cf + sf*ct -- sin (theta + phi)
        st = sin theta
        sf = sin phi
        ct = cos theta
        cf = cos phi
        rr = case alpha of
          TensionAtLeast t | bounded ->
            min rr' (st/stf)
          _ -> rr'
        ss = case beta of
          TensionAtLeast t | bounded ->
            min ss' (sf/stf)
          _ -> ss'
        u = z0 ^+^ rr *^ Point (dx*ct - dy*st) (dy*ct + dx*st)  -- (rotate theta $* chord)
        v = z1 ^-^ ss *^ Point (dx*cf + dy*sf) (dy*cf - dx*sf)  -- (rotate (-phi) $* chord)

unmetaJoin _ _ _ = error "Cannot convert join."

constant1 :: Double
constant1 = 3*(1 + (sqrt 5 - 1))/2

constant2 :: Double
constant2 = 3*(3 - sqrt 5)/2

-- another magic formula by John Hobby.
velocity :: Double -> Double -> Double -> Double
         -> Double -> Double -> Tension -> Double
velocity theta phi st sf ct cf t =
  (2 + sqrt 2 * (st - sf/16)*(sf - st/16)*(ct - cf)) /
  (3 + constant1*ct + constant2*cf) * tensionValue t

