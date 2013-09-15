module Geom2D.CubicBezier.MetaPath
where
import Geom2D.CubicBezier
import Test.QuickCheck
import Data.List

  
-- data MetaPath = OpenMetaPath MetaNode Curl [(MetaJoin, MetaNode)] MetaJoin MetaNode Curl
--               | ClosedMetaPath MetaNode [(MetaJoin, MetaNode)] MetaJoin

-- newtype Curl = Curl Double

-- data MetaNode = MetaNode Point (Maybe Point) -- node and direction
-- data MetaJoin = JoinControls Point Point
--               | JoinTension Tension Tension

-- unmeta (MetaPath) = undefined

-- solve the tridiagonal system for t[i]:
-- a[n] t[i-1] + b[i] t[i] + c[b] t[i+1] = d[i]
-- where a[0] = c[n] = 0
-- by first rewriting it into
-- the system t[i] + u[i] t[i+1] = v[i]
-- where u[n] = 0
-- then solving for t[n]
-- see metafont the program: ¶ 283
solveTriDiagonal :: (Double, Double, Double)
                 -> [(Double, Double, Double, Double)]
                 -> [Double]
solveTriDiagonal _ [] = error "tridiagonal: not enough equations"
solveTriDiagonal (b0, c0, d0) rows = solutions
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
solveCyclicTriD = 