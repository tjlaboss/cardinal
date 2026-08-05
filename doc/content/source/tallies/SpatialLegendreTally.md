# SpatialLegendreTally

## Description

The `SpatialLegendreTally` class wraps an OpenMC tally with three
`SpatialLegendreFilter`s (one each for the $x$, $y$, and $z$ axes) to compute a
Cartesian [!ac](FET) of a score over a box-shaped region. Unlike a
[CellTally](CellTally.md) or [MeshTally](MeshTally.md), a `SpatialLegendreTally`
does not create any `ElementalAuxVariable`s. Instead, it populates a MOOSE
`FunctionSeries` (from the functional expansion tools module) named
`<score>_<function_suffix>` with the expansion coefficients. That continuous
function is then projected onto the finite element solution with a
`FunctionSeriesToAux` auxkernel, reusing MOOSE's own [!ac](FET) reconstruction
machinery rather than writing piecewise-constant values directly.

The expansion orders in each direction are set with `orders`, and the physical
bounds of the expansion are set with `lower_left` and `upper_right`. Only the
collision estimator is supported.

## Example Input File Syntax

!listing /test/tests/neutronics/filters/fet/legendre_integral.i
  block=Problem

!syntax parameters /Problem/Tallies/SpatialLegendreTally

!syntax inputs /Problem/Tallies/SpatialLegendreTally
