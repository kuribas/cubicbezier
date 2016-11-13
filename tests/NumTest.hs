module NumTest (numTests) where
import Test.Tasty
import Test.Tasty.HUnit
import Geom2D.CubicBezier.Numeric
import Data.Matrix.Unboxed as M
import Data.Vector.Unboxed as V

m1 :: Matrix Double
m1 = matrix 5  [
  7.31151,    8.22469,    14.8986,  17.4782,     2.26534,
 17.1487,     6.1471,     16.269,    4.72402,    6.44706,
 10.8667,     0.309989,   8.19124,   5.56535 ,   5.75257,
 10.4157,    12.2949,     7.37355,   8.25757,   15.8631,
  1.62244,    5.93015,    2.60174,  12.1159,     2.76559,
  0.370118,   3.71538,    3.95067,  13.7553,    11.1942,
 15.5189,     6.16186,   17.444,     0.259559,  19.8304,
  5.85888,   10.5329,     1.5679,    4.2007,    14.9874,
 13.1124,    18.1901,    16.2589,   10.1622,    10.966,
  2.01525,    0.663806,  11.305,     7.15144,    2.9295,
  32, 1, 2, 5, 3,
  5, 6, 22, 4, 5,
  6, 23, 54, 1, 4]

v1:: Vector Int
v1 = V.fromList [0, 0, 0, 0, 3, 3, 3, 3, 6, 6, 6, 6, 6]

v2 :: Vector Double
v2 = V.fromList [8.12236, 9.67724, 5.88729, 6.37359, 4.91577,
                 2.60334, 2.88849, 7.80152, 4.32252, 0.813729, 2, 5, 1]

numTests = undefined

mpTests :: TestTree
mpTests = testGroup "Numeric" []

