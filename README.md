cubicbezier
===========

Haskell library for manipulating cubic bezier curves.  The original goal
is to support typography, but it may be useful for general graphics.

Features implemented:
 * finding least squares approximation of a curve through some points
 * finding tangents parallel to a vector
 * evaluating bezier curves and derivatives
 * finding curvature and radius of curvature

Features todo:

 * finding intersection points
 * approximating an outline
 * creating paths from meta paths (as in D.E.Knuth's _metafont_)
 * adding and subtracting closed paths
 * removing overlap from paths
 * affine transformations on bezier curves

 other TODO:
  * eliminate large dependencies (hmatrix)