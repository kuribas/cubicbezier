cubicbezier
===========

Haskell library for manipulating cubic and quadratic bezier curves.
It is meant as a low level library to support 2D graphics and
typography applications.

Features implemented:

 * evaluating bezier curves and derivatives
 * approximation of a curve through some points
 * removing overlap and boolean operations on paths
 * finding tangents parallel to a vector
 * curvature and radius of curvature
 * intersections between two curves
 * intersections between a curve and a line
 * finding inflection points and cusps
 * affine transformations on bezier curves
 * creating paths from meta paths (as in D.E.Knuth's _metafont_)
 
Features todo:

 * calligraphic strokes
