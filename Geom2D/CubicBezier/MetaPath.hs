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

data Tension = Tension Double
             | TensionAtLeast Double

-- | Create a normal path from a metapath.
unmeta (OpenMetaPath nodes endpoint) =
  let subsegs = openSubSegments (sanitizeOpen nodes) endpoint
      angles = map (solveTriDiagonal . eqsOpen) subsegs
      newnodes = zipWith insertAngles nodes angles
      points = map fst newnodes ++ endpoint
      path = zipWith (uncurry unmetaJoin) newnodes (tail points)
  in Path (head points) (zip path (tail points))



-- replace open nodetypes if possible
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
      (p, Given dir t) : sanitizeRest (node2:rest)

sanitizeRest ((p, m): (q, n): rest) =
  case (m, n) of
    (Explicit u v, Implicit (Open t) mt2) ->  -- explicit, open => explicit, given
      (p, m) : sanitizeRest ((q, Implicit (Given dir t) n): rest)
      where dir = vectorAngle (q^-^v)
    (Implicit mt1 (Open t), Explicit u v) ->  -- open, explicit => given, explicit
      (p, Implicit m (Given dir t)) : sanitizeRest ((q, n): rest)
      where dir = vectorAngle (u^-^p)

sanitizeRest (n:l) = n:sanitizeRest l

breakPoint ((_,  Implicit _ (Open _)):(_, Implicit (Open _) _):_) = False
breakPoint _ = True

-- decompose into a list of subsegments that need to be solved.
openSubSegments l = openSubSegments' (tails l)

openSubSegments' [[]] _ = []
openSubSegments' [] _   = []
openSubSegments' l lastPoint = case break breakPoint l of
  (m, n:o) ->
    let point = case o of
          ((p,_):_) -> p
          otherwise -> lastPoint
    in OpenMetaPath (map head (m ++ n)) point :
       openSubSegments o
  otherwise -> "unexpected end of segments"

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

turnAngle p q = vectorAngle q - vectorAngle p
turnAngle3 p q r = turnAngle (q^-^p) (r^-^q)

zipPrev [] = []
zipPrev l = zip (last l : l) l

zipNext [] = []
zipNext l = zip l (tail l)

-- neighbours [] = []
-- neighbours [a] = [(a, a, a)]
-- zipNeighbours f ls@(l : ls'@(l' : ls'')) =
--   ((last ls'), l, l') : zipWith3 f ls ls' (ls''++ls)

-- find the equations for a cycle containing only open points
eqsCycle tensionsA points tensionsB = 
  zipWith4 findEqTension
  (zipPrev tensionsA)
  (zipPrev dists)
  (zipPrev turnAngles)
  (zipPrev tensionsB)
  where
    turnAngles = zipWith turnAngle chords (tail $ cycle chords)
    dists = zipWith vectorDistance points (tail $ cycle points)
    chords = zipWith (^-^) points (withLast points)

-- find the equations for a subsegment consisting of a number of open
-- points starting and ending with either a curl or a given direction.
eqsOpen [] = error "empty subsegment list"
    
eqsOpen [MetaNode p0 _ rn0, MetaNode p1 ln1 _] =
  case (rn0, ln1) of
    (Curl g t, Curl g2 t2)     -> [(0, 1, 0, a), (0, 1, 0, a)]
      where a = vectorAngle (p1^-^p0)
    (Curl g t, Given dir t2) -> [eqCurl0 g t t2 0, (0, 1, 0, dir)]
    (Given dir t, Curl g t2) -> [(0, 1, 0, dir), eqCurlN g t t2]
    otherwise                -> error "illegal nodetype in subsegment"

eqsOpen points joins =
  eq0 : restEquations nodeRest tensionsA dists turnAngles tensionsB
  where
    dists = zipWith vectorDistance chords (tail chords)
    chords = zipWith (^-^) (tail points) points
    tensionsA = map (getTension . metaJoin1) joins
    tensionsB = map (getTension . metaJoin2) joins
    turnAngles = zipWith turnAngle chords (tail chords)

    eq0 = case head joins of
      (_, Implicit (Curl g t) _) -> eqCurl0 g (head tensionsA) (head tensionsB) (head turnAngles)
      (_, Implicit MetaNode _ _ (Given angle _) -> (0, 1, 0, angle)
      otherwise -> error "illegal subsegment first nodetype"

    restEquations [lastnode] (tensionA:_) _ _ (tensionB:_) =
      case lastnode of
        MetaNode _ (Curl _ _) _ -> eqCurlN g tensionA tensionB
        MetaNode _ (Given angle _) -> (0, 1, 0, angle)
        otherwise -> error "illegal subsegment last nodetype"

    restEquations (node:othernodes) (tensionA:restTA) (d:restD) (turn:restTurn) (tensionB:restTB) =
      eqTension (tensionA, head restTA) (d, head restD) (turn, head restTurn) (tensionB, head restTB) :
      restEquations othernodes restTA restD retTurn restTB

eqTension (tensionA', tensionA) (d', d) (psi', psi) (tensionB', tensionB) =
  (a, b+c, d, -b*psi' - d*psi)
  where
    a = (tensionB' * tensionB' / (tensionA' * d'))
    b = (3 - 1/tensionA') * tensionB' * tensionB' / d'
    c = (3 - 1/tensionB) * tensionA * tensionA / d
    d = tensionA * tensionA / (tensionB * d))
    
eqCurl0 (Curl gamma) tensionA tensionB psi = (0, c, d, r)
  where
    c = chi/tensionA + 3 - 1/tensionB
    d = (3 - 1/tensionA)*chi + 1/tensionB
    chi = gamma*tensionB*tensionB / tensionA*tensionA
    r = -d*psi

eqCurlN (Curl gamma) tensionA tensionB = (a, b, 0, 0)
  where
    a = (3 - 1/tensionB)*chi + 1/tensionA
    b = chi/tensionB + 3 - 1/tensionA
    chi = gamma*tensionA*tensionA / tensionB*tensionB

insertAngles _ [] = []
insertAngles ((Implicit m n): rest) (a:b:angles) =
  Implicit (Given a (getTension m)) (Given b (getTension n)) :
  insertAngles rest (b:angles)

-- turn a metajoin into a simple join
unmetaJoin _ (Explicit u v) _  = JoinCurve u v

-- magic formula for getting the control points by John Hobby
unmetaJoin z0 (Implicit (Given w0 alpha) (Given w1 beta)) z1
  | abs theta < 1e-5 && abs phi < 1e-5 = (z0, JoinLine)
  | otherwise = JoinCurve u v
  where chord = z1^-^z0
        chordAngle = vectorAngle chord
        theta = chordAngle - w0
        phi   = w1 - chordAngle
        u = z0 ^+^ hobbyMagic theta phi / alpha *^ (rotate theta $* chord)
        v = z1 ^-^ hobbyMagic phi theta / beta *^ (rotate (-phi) $* chord)

unmetaJoin _ _ _ = error "Cannot convert join."

-- another magic formula by John Hobby.
hobbyMagic theta phi = (2 + sqrt2 * (sinTheta - sinPhi/16)*(sinPhi - sinTheta/16)*(cosTheta - cosPhi)) /
                       3*(1 + (sqrt5 - 1)*cosTheta/2 + (3 - sqrt5)*cosPhi/2)
  where sqrt2 = sqrt 2
        sqrt5 = sqrt 5
        sinTheta = sin theta
        sinPhi = sin phi
        cosTheta = cos theta
        cosPhi = cos phi
