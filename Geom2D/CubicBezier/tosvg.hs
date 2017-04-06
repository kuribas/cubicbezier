{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE QuasiQuotes #-}
import Control.Applicative
import Control.Monad.State
import Data.Char
import Data.Maybe
import Data.List
import Data.String.Here.Interpolated
import Data.Hashable
import GHC.Generics
import qualified Data.HashMap.Strict as H
import qualified Data.HashSet as S
import System.Environment

data Curve = Curve Point Point Point Point
  deriving (Eq, Generic)
data Point = Point Double Double
  deriving (Eq, Generic)

instance Hashable Curve
instance Hashable Point

data Trace =
  XStruct [Curve] |
  MSG String |
  ChangeFocus Point |
  XStructAdd Curve |
  XStructRem Curve |
  YStructAdd Curve |
  YStructRem Curve |
  Activate Curve |
  DeActivate Curve |
  Output Curve |
  Discard Curve

skipChar :: Char -> StateT String Maybe ()
skipChar c = StateT $ \case
  (c1:r) | c == c1 -> Just ((), r)
  _ -> Nothing

readMb :: Read a => StateT String Maybe a
readMb = StateT $ listToMaybe . reads

readCurve :: String -> Maybe Curve
readCurve = evalStateT $ do
  p0 <- readPoint
  p1 <- readPoint
  p2 <- readPoint
  p3 <- readPoint
  return $ Curve p0 p1 p2 p3

readPoint :: StateT String Maybe Point
readPoint = do
  skipChar '('
  x <- readMb
  skipChar ','
  y <- readMb
  skipChar ')'
  return $ Point x y

parseTrace :: [String] -> [Trace]
parseTrace [] = []
parseTrace (str:rest) = case pref of
  "MSG" -> MSG arg:parseTrace rest
  "XSTRUCTADD" -> ret XStructAdd
  "XSTRUCTREM" -> ret XStructRem
  "YSTRUCTADD" -> ret YStructAdd
  "YSTRUCTREM" -> ret YStructRem
  "ACTIVATE"   -> ret Activate
  "DEACTIVATE" -> ret DeActivate
  "OUTPUT"     -> ret Output
  "DISCARD"    -> ret Discard
  "CHANGEFOCUS" -> maybe id (:)
                   (ChangeFocus <$> evalStateT readPoint arg)
                   (parseTrace rest)
  "XSTRUCTBEGIN" ->
    let (curves, rest2) = span (/= "XSTRUCTEND") rest
    in (XStruct $ mapMaybe readCurve curves) : 
       (parseTrace $ case rest2 of [] -> []; (_:r) -> r)
  _ -> parseTrace rest
  where
    ret cons = maybe id (:) (cons <$> readCurve arg) $ parseTrace rest
    (pref, str2) = span isAlpha str
    arg = dropWhile isSpace str2

getCurve :: Trace -> [Curve]
getCurve (XStruct cs) = cs
getCurve (XStructAdd c) = [c]
getCurve (XStructRem c) = [c]
getCurve (YStructAdd c) = [c]
getCurve (YStructRem c) = [c]
getCurve (Activate c) = [c]
getCurve (DeActivate c) = [c]
getCurve (Output c) = [c]
getCurve (Discard c) = [c]
getCurve _ = []

getXStruct (XStruct cs) = Just cs
getXStruct _ = Nothing

makeObjects :: Double -> H.HashMap Curve Int -> [Trace] -> String
makeObjects scale hm tr = unlines $ concatMap makeObject $ H.toList hm
  where xStruct :: S.HashSet Curve
        xStruct = case listToMaybe $ mapMaybe getXStruct tr of
          Just x -> S.fromList x
          _ -> S.empty
        makeObject (c@(Curve (Point p0x p0y) (Point p1x p1y) (Point p2x p2y) (Point p3x p3y)), id) =
          [[i|<path id="p${id}" d="M${p0x} ${p0y} C${p1x} ${p1y} ${p2x} ${p2y} ${p3x} ${p3y}" class="${if S.member c xStruct then "xstruct" else "hidden"}"/>|]
          ,[i|<circle id="lp${id}" cx="${p0x}" cy="${p0y}" class="xstruct" r="${1.5/scale}"/>|]
          ,[i|<circle id="rp${id}" cx="${p3x}" cy="${p3y}" class="xstruct" r="${1.5/scale}"/>|]]
          
makeAction :: H.HashMap Curve Int -> [Trace] -> [String]
makeAction _ [] = []
makeAction hm (XStruct cs:r) = makeAction hm r
makeAction hm (MSG s: r)
  | null steps = makeAction hm r
  | otherwise = msg : makeAction hm r2
  where
    (actions, r2) = span (isJust.setAction hm) r
    steps = concat $ intersperse "\n," $ mapMaybe (setAction hm) actions
    msg = [i|{message: "${s}", action: "set", steps:
[${steps}]}|]
makeAction hm (ChangeFocus (Point x y):r) =
  ([i|{message: "change focuspoint to (${x}, ${y}).", action: "focus", next: [${x}, ${y}], prev: [${x}, ${y}]}|]
    ++ "\n") : makeAction hm r
makeAction hm r = makeAction hm (MSG "": r)

setAction hm (XStructAdd c) =
  Just [i|{id: ${fromJust $ H.lookup c hm}, set: true, class: "xstruct"}|]
setAction hm (XStructRem c) =
  Just [i|{id: ${fromJust $ H.lookup c hm}, set: false, class: "xstruct"}|]
setAction hm (YStructAdd c) = 
  Just [i|{id: ${fromJust $ H.lookup c hm}, set: true, class: "ystruct"}|]
setAction hm (YStructRem c) =
  Just [i|{id: ${fromJust $ H.lookup c hm}, set: false, class: "ystruct"}|]
setAction hm (Activate c) = 
  Just [i|{id: ${fromJust $ H.lookup c hm}, set: true, class: "active"}|]
setAction hm (DeActivate c) =
  Just [i|{id: ${fromJust $ H.lookup c hm}, set: false, class: "active"}|]
setAction hm (Output c) = 
  Just [i|{id: ${fromJust $ H.lookup c hm}, set: true, class: "out"}|]
setAction hm (Discard c) =
  Just [i|{id: ${fromJust $ H.lookup c hm}, set: true, class: "discarded"}|]
setAction _ _ = Nothing

outputString :: Int -> Int -> Double -> Double -> Double -> String -> String -> (Double, Double) -> String  
outputString width height scale transX transY actions objects focus =
  [template|./tosvg.svg|]

curveMap :: [Trace] -> H.HashMap Curve Int
curveMap trace = curves
  where
    curves = fst $ foldl addCurve (H.empty, 1) $
             concatMap getCurve trace
    addCurve (mp, i) c
      | H.member c mp = (mp, i)
      | otherwise = (H.insert c i mp, i+1)

toSvg :: Int -> Int -> String -> String
toSvg width height input = outputString width height scale transX transY actions objects (findFocus trace)
  where trace = parseTrace $ lines input
        cm = curveMap trace
        actions = concat $ intersperse ",\n " $ makeAction cm trace
        objects = makeObjects scale cm trace
        (minX, minY, maxX, maxY) = limits cm
        scale = min (fromIntegral (width-20) / (maxX - minX))
                (fromIntegral (height-40) / (maxY - minY))
        transX = 10 - minX * scale
        transY = 10 - minY * scale
        findFocus (ChangeFocus (Point x y):_) = (x, y)
        findFocus (_:r) = findFocus r
        findFocus [] = error "no focuspoint found"

limits :: H.HashMap Curve Int -> (Double, Double, Double, Double)
limits = foldl getLimit (1/0, 1/0, -1/0, -1/0) . H.keys 
  where
    getLimit (minX, minY, maxX, maxY)
      (Curve (Point p0x p0y) (Point p1x p1y) (Point p2x p2y) (Point p3x p3y)) =
      (minimum [p0x, p1x, p2x, p3x, minX],
       minimum [p0y, p1y, p2y, p3y, minY],
       maximum [p0x, p1x, p2x, p3x, maxX],
       maximum [p0y, p1y, p2y, p3y, maxY]
      )

usage = error "usage: tosvg [-w width] [-h height]"

parseArgs :: (Int, Int) -> [String] -> IO (Int, Int)
parseArgs (w, h) ("-w":sw:r) =
  case reads sw of
    [] -> usage
    ((i,_):_) -> parseArgs (i, h) r

parseArgs (w, h) ("-h":sh:r) =
  case reads sh of
    [] -> usage
    ((i,_):_) -> parseArgs (w, i) r
parseArgs def _ = return def

main = do
  (width, height) <- parseArgs (1200, 700) =<< getArgs
  interact (toSvg width height)
