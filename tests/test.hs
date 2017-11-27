import Test.Tasty
import Test.Tasty.HUnit
import Geom2D.CubicBezier
import Control.Monad
import Text.Parsec
import Text.Parsec.String
import Text.Parsec.Error
import MFTest (mftests) where

tests :: TestTree
tests = testGroup "Tests" [mfTests]

main :: IO ()
main = defaultMain tests
