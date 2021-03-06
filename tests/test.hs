import Test.Tasty
import Test.Tasty.HUnit
import Geom2D.CubicBezier
import Control.Monad
import Text.Parsec
import Text.Parsec.String
import Text.Parsec.Error

num :: Parser Double
num = 
  liftM (read.concat) $ sequence
  [ option "" $ string "-"
  , many1 digit
  , option "" $ string "."
  , option "0" $ many digit]

pointP :: Parser DPoint
pointP = do
  char '('; spaces
  n <- num; spaces
  char ','; spaces
  m <- num ; spaces
  char ')'
  return (Point n m)

nodeP :: Parser (MetaNodeType Double)
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

tensionAmount :: Parser (Tension Double)
tensionAmount = do
  cons <- option Tension
          (string "atleast" >>
           return TensionAtLeast)
  spaces
  n <- num
  return $ cons n
         
tensionP :: Parser (Tension Double, Tension Double)
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

openMetaRest :: DPoint -> Parser (OpenMetaPath Double)
openMetaRest p = do
  leftNode <- nodeP; spaces
  string ".."; spaces
  (tl, tr) <- tensionP; spaces
  rightNode <- nodeP; spaces
  OpenMetaPath joins q <- openMetaP
  return $ OpenMetaPath ((p, MetaJoin leftNode tl
                          tr rightNode):joins) q

closedMetaRest :: DPoint -> Parser (ClosedMetaPath Double)
closedMetaRest p = do
  leftNode <- nodeP; spaces
  string ".."; spaces
  (tl, tr) <- tensionP; spaces
  rightNode <- nodeP; spaces
  ClosedMetaPath joins <- closedMetaP
  return $ ClosedMetaPath ((p, MetaJoin leftNode tl
                             tr rightNode):joins)

openRest :: DPoint -> Parser (OpenPath Double)
openRest p = do
  string ".."; spaces
  string "controls"; spaces
  p1 <- pointP; spaces
  string "and"; spaces
  p2 <- pointP; spaces
  string ".."; spaces
  OpenPath joins q <- openP
  return $ OpenPath ((p, JoinCurve p1 p2):joins) q

closedRest :: DPoint -> Parser (ClosedPath Double)
closedRest p = do
  string ".."; spaces
  string "controls"; spaces
  p1 <- pointP; spaces
  string "and"; spaces
  p2 <- pointP; spaces
  string ".."; spaces
  ClosedPath joins <- closedP
  return $ ClosedPath ((p, JoinCurve p1 p2):joins)

openP :: Parser (OpenPath Double)
openP = 
  do p <- pointP
     spaces
     option (OpenPath [] p) (openRest p)

closedP :: Parser (ClosedPath Double)
closedP =
  do p <- pointP
     spaces
     closedRest p
  <|> do
    string "cycle"
    return (ClosedPath [])

openMetaP :: Parser (OpenMetaPath Double)
openMetaP =
  do p <- pointP
     spaces
     option (OpenMetaPath [] p) (openMetaRest p)

closedMetaP :: Parser (ClosedMetaPath Double)
closedMetaP =
  do p <- pointP
     spaces
     closedMetaRest p
  <|> do
    string "cycle"
    return (ClosedMetaPath [])

tryParse :: Parser a -> String -> a
tryParse p s =
  case parse p "" s of
   Left err -> error $ concatMap messageString $
               errorMessages err
   Right res -> res
  
  
doubleEq :: (Ord a, Fractional a) => a -> a -> Bool
doubleEq a b =
  abs (a - b) < 0.01

pointEq :: DPoint -> DPoint -> Bool
pointEq (Point a b) (Point c d) =
  doubleEq a c && doubleEq b d

joinEq :: PathJoin Double -> PathJoin Double -> Bool
joinEq JoinLine JoinLine  = True
joinEq (JoinCurve a b) (JoinCurve c d) =
  pointEq a c && pointEq b d
joinEq _ _ = True

openPathEq (OpenPath joins p) (OpenPath joins2 q) =
  pointEq p q && length joins == length joins2 &&
  and (zipWith
   (\(p1, j1) (p2, j2) ->
     pointEq p1 p2 && joinEq j1 j2)
   joins joins2)

closedPathEq (ClosedPath joins) (ClosedPath joins2) =
  and (zipWith
   (\(p1, j1) (p2, j2) ->
     pointEq p1 p2 && joinEq j1 j2)
   joins joins2)

openThetas :: OpenPath Double -> [Double]
openThetas (OpenPath j p) =
  zipWith3 theta
  (map fst j)
  (tail (map fst j) ++ [p])
  (map snd j)
  where
    theta q r JoinLine = 0
    theta q r (JoinCurve c1 _) =
      vectorAngle (c1^-^q) - vectorAngle (r^-^q)
    
closedThetas (ClosedPath j) =
  openThetas (OpenPath j (fst $ head j))

openPhis :: OpenPath Double -> [Double]
openPhis (OpenPath j p) =
  zipWith3 phi
  (map fst j)
  (tail (map fst j) ++ [p])
  (map snd j)
  where
    phi q r JoinLine = 0
    phi q r (JoinCurve _ c2) =
      vectorAngle (q^-^r) - vectorAngle (c2^-^r)

closedPhis (ClosedPath j) =
  openPhis (OpenPath j (fst $ head j))

testOpen :: TestName -> String -> TestTree
testOpen p1 p2 =
  testCase p1 $ 
  assertBool "Incorrect metapath." $
  unmetaOpen (tryParse openMetaP p1) `openPathEq`
  tryParse openP p2

testClosed :: TestName -> String -> TestTree
testClosed p1 p2 =
  testCase p1 $ 
  assertBool "Incorrect metapath." $
  unmetaClosed (tryParse closedMetaP p1) `closedPathEq`
  tryParse closedP p2

-- These tests were created by running mf, typing expr after the
-- prompt, and entering the metapaths.
mfTests :: TestTree
mfTests = testGroup "Metafont" [
  testOpen "(0,0)..(4,3)"
  "(0,0)..controls (1.33333,1) and (2.66667,2) ..(4,3)",

  testOpen "(0,0){(1,-2)}..(4,3)"
  "(0,0)..controls (1.81548,-3.63095) and (6.97739,0.24046)..(4,3)",

  testOpen "(0,0)..{(1,-2)}(4,3)"
  "(0,0)..controls (-2.97739,2.75954) and (2.18452,6.63095)..(4,3)",

  testOpen "(0,0){curl 2}..(4,3)"
  "(0,0)..controls (1.33333,1) and (2.66667,2)..(4,3)",

  testOpen "(0,0){(2, 3)}..{(1, 2)}(4,3)"
  "(0,0)..controls (0.95523,1.43285) and (3.21622,1.43243)..(4,3)",

  testOpen "(0,0)..(4,3)..(-2, 1)"
  "(0,0)..controls (2.08194,-1.42896) and (4.78885,0.60123)..(4,3)..controls (2.67747,7.02158) and (-3.35492,5.01077)..(-2,1)",

  testOpen "(1,1)..tension 0.8 and 1.2..(3,4)..tension 10 ..(-10,-10)"
  "(1,1)..controls (-2.7088,-12.93713) and (13.27118,14.12433)..(3,4)..controls (2.54623,3.55272) and (-9.58751,-9.5144)..(-10,-10)",
  
  testOpen "(0,0){curl 2}..(4,3)..(-2, 1)"
  "(0,0)..controls (1.14464,-2.66646) and (6.04007,-0.56508)..(4,3)..controls (2.2501,6.05801) and (-2.43489,4.49635)..(-2,1)",

  testOpen "(0,0){(-3, -2)}..(4,3)..(-2, 1)"
  "(0,0)..controls (-3.65675,-2.43784) and (1.35551,2.07506)..(4,3)..controls (27.8797,11.35223) and (-26.11505,-6.64606)..(-2,1)",

  testClosed "(0,0)..(2,3)..(4,4)..cycle"
  "(0,0)..controls (-0.27211,1.267) and (0.9676,2.15346)..(2,3)..controls (2.60509,3.49615) and (3.2241,4.08679)..(4,4)..controls (12.90535,3.00386) and (1.91997,-8.93997)..cycle",

  testClosed "(0,0)..tension 0.9 and 1.1 ..(2,3)..(4,4)..cycle"
  "(0,0)..controls (-0.39941,1.39384) and (0.99234,2.26094)..(2,3)..controls (2.62666,3.45963) and (3.22433,4.07909)..(4,4)..controls (12.2955,3.15413) and (2.40324,-8.38663)..cycle",

  testClosed "(0,0)..(2,3){(1,1)}..(4,4)..cycle"
  "(0,0)..controls (-0.24208,1.27483) and (1.07744,2.07744)..(2,3)..controls (2.56248,3.56248) and (3.22197,4.11229)..(4,4)..controls (12.86206,2.72092) and (1.68616,-8.87949)..cycle"
  ]

tests :: TestTree
tests = testGroup "Tests" [mfTests]


main :: IO ()
main = defaultMain tests
