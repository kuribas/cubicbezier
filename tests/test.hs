import Test.Tasty
import Test.Tasty.HUnit
import Geom2D.CubicBezier
import Control.Monad
import Text.Parsec
import Text.Parsec.String
import Text.Parsec.Error
import MPTest
import NumTest

tests :: TestTree
tests = testGroup "Tests" [mfTests, numTests]

num :: Parser Double
num = 
  liftM (read.concat) $ sequence
  [ option "" $ string "-"
  , many1 digit
  , option "" $ string "."
  , option "0" $ many digit]

pointP :: Parser Point
pointP = do
  char '('; spaces
  n <- num; spaces
  char ','; spaces
  m <- num ; spaces
  char ')'
  return (Point n m)

nodeP :: Parser MetaNodeType
nodeP = option Open specialNode
  where specialNode = do
          char '{'; spaces
          node <- choice [
            do string "curl"
               spaces
               n <- num
               return (Curl n),
            do p <- pointP
               return $ Direction p]
          spaces; char '}'
          return node

tensionAmount :: Parser Tension
tensionAmount = do
  cons <- option Tension
          (string "atleast" >>
           return TensionAtLeast)
  spaces
  n <- num
  return $ cons n
         
tensionP :: Parser (Tension, Tension)
tensionP =
  option (Tension 1, Tension 1) $
  do string "tension";
     spaces;
     t1 <- tensionAmount
     spaces
     t2 <- option t1 (do string "and"
                         spaces
                         tensionAmount)
     spaces
     string ".."
     return (t1, t2)

mpRest :: Point -> Parser MetaPath
mpRest p = do
  leftNode <- nodeP; spaces
  string ".."; spaces
  (tl, tr) <- tensionP; spaces
  rightNode <- nodeP; spaces
  mp <- mpP
  return $ case mp of
    OpenMetaPath joins q ->
      (OpenMetaPath ((p, MetaJoin leftNode tl
                         tr rightNode):joins) q)
    CyclicMetaPath joins ->
      CyclicMetaPath ((p, MetaJoin leftNode tl
                         tr rightNode):joins)

mpP :: Parser MetaPath
mpP =
  do p <- pointP
     spaces
     option (OpenMetaPath [] p) (mpRest p)
  <|> do
    string "cycle"
    return (CyclicMetaPath [])

pathRest :: Point -> Parser Path
pathRest p = do
  string ".."; spaces
  string "controls"; spaces
  n <- pointP; spaces
  string "and"; spaces
  m <- pointP; spaces
  string ".."; spaces
  path <- pathP
  return $ case path of
    OpenPath joins q ->
      (OpenPath ((p, JoinCurve n m):joins) q)
    ClosedPath joins ->
      ClosedPath ((p, JoinCurve n m):joins)
  
pathP :: Parser Path
pathP =
  do p <- pointP
     spaces
     option (OpenPath [] p) (pathRest p)
  <|> do
    string "cycle"
    return (ClosedPath [])

tryParse :: Parser a -> String -> a
tryParse p s =
  case parse p "" s of
   Left err -> error $ concatMap messageString $
               errorMessages err
   Right res -> res
  
  
doubleEq :: (Ord a, Fractional a) => a -> a -> Bool
doubleEq a b =
  abs (a - b) < 0.01

pointEq :: Point -> Point -> Bool
pointEq (Point a b) (Point c d) =
  doubleEq a c && doubleEq b d

joinEq :: PathJoin -> PathJoin -> Bool
joinEq JoinLine JoinLine = True
joinEq (JoinCurve a b) (JoinCurve c d) =
  pointEq a c && pointEq b d
joinEq _ _ = True

pathEq :: Path -> Path -> Bool
pathEq (OpenPath joins p) (OpenPath joins2 q) =
  pointEq p q && length joins == length joins2 &&
  and (zipWith
   (\(p1, j1) (p2, j2) ->
     pointEq p1 p2 && joinEq j1 j2)
   joins joins2)

pathEq (ClosedPath joins) (ClosedPath joins2) =
  and (zipWith
   (\(p1, j1) (p2, j2) ->
     pointEq p1 p2 && joinEq j1 j2)
   joins joins2)

thetas :: Path -> [Double]
thetas (OpenPath j p) =
  zipWith3 theta
  (map fst j)
  (tail (map fst j) ++ [p])
  (map snd j)
  where
    theta q r (JoinLine) = 0
    theta q r (JoinCurve c1 _) =
      vectorAngle (c1^-^q) - vectorAngle (r^-^q)
    
thetas (ClosedPath j) =
  thetas (OpenPath j (fst $ head j))

phis :: Path -> [Double]
phis (OpenPath j p) =
  zipWith3 phi
  (map fst j)
  (tail (map fst j) ++ [p])
  (map snd j)
  where
    phi q r (JoinLine) = 0
    phi q r (JoinCurve _ c2) =
      vectorAngle (q^-^r) - vectorAngle (c2^-^r)

phis (ClosedPath j) =
  phis (OpenPath j (fst $ head j))

testPath :: TestName -> String -> TestTree
testPath p1 p2 =
  testCase p1 $ 
  assertBool "Incorrect metapath." $
  unmeta (tryParse mpP p1) `pathEq`
  tryParse pathP p2

-- These tests were created by running mf, typing expr after the
-- prompt, and entering the metapaths.
mfTests :: TestTree
mfTests = testGroup "Metafont" [
  testPath "(0,0)..(4,3)"
  "(0,0)..controls (1.33333,1) and (2.66667,2) ..(4,3)",

  testPath "(0,0){(1,-2)}..(4,3)"
  "(0,0)..controls (1.81548,-3.63095) and (6.97739,0.24046)..(4,3)",

  testPath "(0,0)..{(1,-2)}(4,3)"
  "(0,0)..controls (-2.97739,2.75954) and (2.18452,6.63095)..(4,3)",

  testPath "(0,0){curl 2}..(4,3)"
  "(0,0)..controls (1.33333,1) and (2.66667,2)..(4,3)",

  testPath "(0,0){(2, 3)}..{(1, 2)}(4,3)"
  "(0,0)..controls (0.95523,1.43285) and (3.21622,1.43243)..(4,3)",

  testPath "(0,0)..(4,3)..(-2, 1)"
  "(0,0)..controls (2.08194,-1.42896) and (4.78885,0.60123)..(4,3)..controls (2.67747,7.02158) and (-3.35492,5.01077)..(-2,1)",

  testPath "(1,1)..tension 0.8 and 1.2..(3,4)..tension 10 ..(-10,-10)"
  "(1,1)..controls (-2.7088,-12.93713) and (13.27118,14.12433)..(3,4)..controls (2.54623,3.55272) and (-9.58751,-9.5144)..(-10,-10)",
  
  testPath "(0,0){curl 2}..(4,3)..(-2, 1)"
  "(0,0)..controls (1.14464,-2.66646) and (6.04007,-0.56508)..(4,3)..controls (2.2501,6.05801) and (-2.43489,4.49635)..(-2,1)",

  testPath "(0,0){(-3, -2)}..(4,3)..(-2, 1)"
  "(0,0)..controls (-3.65675,-2.43784) and (1.35551,2.07506)..(4,3)..controls (27.8797,11.35223) and (-26.11505,-6.64606)..(-2,1)",

  testPath "(0,0)..(2,3)..(4,4)..cycle"
  "(0,0)..controls (-0.27211,1.267) and (0.9676,2.15346)..(2,3)..controls (2.60509,3.49615) and (3.2241,4.08679)..(4,4)..controls (12.90535,3.00386) and (1.91997,-8.93997)..cycle",

  testPath "(0,0)..tension 0.9 and 1.1 ..(2,3)..(4,4)..cycle"
  "(0,0)..controls (-0.39941,1.39384) and (0.99234,2.26094)..(2,3)..controls (2.62666,3.45963) and (3.22433,4.07909)..(4,4)..controls (12.2955,3.15413) and (2.40324,-8.38663)..cycle",

    testPath "(0,0)..(2,3){(1,1)}..(4,4)..cycle"
  "(0,0)..controls (-0.24208,1.27483) and (1.07744,2.07744)..(2,3)..controls (2.56248,3.56248) and (3.22197,4.11229)..(4,4)..controls (12.86206,2.72092) and (1.68616,-8.87949)..cycle"
  ]

main :: IO ()
main = defaultMain tests
