Introduction
============

This is a library containing a set of basic basis functions
which show up frequently in finite element methods and meshfree
methods commonly.

Basis Functions
===============

The following basis functions are currently supported:

* Multiquadrics
* Triangular center-of-mass polynomials
* Reproducing kernel particles
* Partition of unity basis (via Shepard method)

The goal of including these functions involves providing
methods for manipulating them in any way required to
generate inner product matrices for common operators.
Derivatives of all basis functions are included among
the basis function sets listed. As of the writing of this
documentation, for some bases (like multquadrics) the
generic derivative functions are incomplete.

Targeted Changes
================

Ideally, we would like to be able generate derivatives
of any basis function to any order and in any dimension.
To do this, algorithmic differentiation is the best route
to take. There are several methods which could be used.
But several libraries for automatical differentiation
are available as Python modules which could be called
from any language.
