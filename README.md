cubicbezier
===========

Haskell library for manipulating cubic bezier curves.  The original goal
is to support typography, but it may be useful for general graphics.

Features implemented:

 * least squares approximation of a curve through some points
 * finding tangents parallel to a vector
 * evaluating bezier curves and derivatives
 * curvature and radius of curvature
 * intersections between two curves
 * intersections between a curve and a line
 * finding inflection points and cusps
 * affine transformations on bezier curves
 * creating paths from meta paths (as in D.E.Knuth's _metafont_)
 
Features todo:

 * self intersections
 * creating an approximate outline for a curve
 * adding and subtracting closed paths
 * removing overlap from paths

