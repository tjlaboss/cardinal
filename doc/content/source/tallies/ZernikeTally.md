# ZernikeTally

## Description

The `ZernikeTally` class wraps an OpenMC tally with a `SpatialLegendreFilter`
(axial, in $z$) and a `ZernikeFilter` (radial, over a disc) to compute a
cylindrical [!ac](FET) of a score, using the `CylindricalDuo` series type. Like
the [SpatialLegendreTally](SpatialLegendreTally.md), a `ZernikeTally` does not
create any `ElementalAuxVariable`s; it populates a MOOSE `FunctionSeries` named
`<score>_<function_suffix>` with the expansion coefficients, which is then
projected onto the finite element solution with a `FunctionSeriesToAux`
auxkernel.

The expansion `orders` are given as the axial (Legendre) order followed by the
radial (Zernike) order. The disc is defined by `radius` and `centroid`, and the
axial extent runs from the centroid up to `max_z` (the lower $z$ bound is
inferred symmetrically about the centroid). Only the collision estimator is
supported.

!syntax parameters /Problem/Tallies/ZernikeTally

!syntax inputs /Problem/Tallies/ZernikeTally
