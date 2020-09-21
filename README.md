# Grundmann-Möller n-dimensional Simplex Quadrature

* [GitHub](https://github.com/eschnett/GrundmannMoeller.jl): Source
  code repository
* [![GitHub
  CI](https://github.com/eschnett/GrundmannMoeller.jl/workflows/CI/badge.svg)](https://github.com/eschnett/GrundmannMoeller.jl/actions)

## Description

This package calculates integrations points and weights for numerical
quadrature (i.e. integrating a function) over an `n`-dimensional
simplex. It supports arbitrary (odd) degrees of accuracy, and supports
arbitrary floating-point types.

## Sample usage

```Julia
julia> using GrundmannMoeller

julia> # Obtain quadrature scheme (2 dimensions, degree 5; degree must be odd)

julia> scheme = grundmann_moeller(Float64, Val(2), 5);

julia> length(scheme.weights)
10

julia> # Apply scheme

julia> vertices = [[0,0], [1,0], [0,1]];

julia> f(x) = x[1]*x[2]
f (generic function with 1 method)

julia> q = integrate(f, scheme, vertices)
0.04166666666666668

julia> q ≈ 1/24
true
```

## Reference

### Create a quadrature scheme

```Julia
    grundmann_moeller(::Type{T}, ::Val{D}, degree::Int)
```
* `T`: desired floating point type
* `D`: dimension
* `degree`: desired polynomial degree of accuracy (must be odd)

### Evaluate an integral (apply the scheme)

```Julia
    integrate(fun, scheme, vertices::AbstractVector)
    integrate(fun, scheme, vertices::AbstractMatrix)
```
* `fun`: integrand, should accept an
  [`SVector`](https://github.com/JuliaArrays/StaticArrays.jl) as
  argument
* `scheme`: quadrature scheme
* `vertices`: vertices of the simplex

The vertices need to be passed either as a vector-of-vectors or as a
matrix. In the first case, there need to be `D+1` points with `D`
coordinates each. In the second case, the matrix needs to have size
`D`×`D+1`.

## Provenance of this package

This package is a Julia translation of one of the algorithms of the
Python [`quadpy`](https://github.com/nschloe/quadpy) package by Nico
Schlömer.

The algorithm itself was published by
* A. Grundmann, H. M. Möller, "Invariant integration formulas for the
  n-simplex by combinatorial methods", SIAM J. Numer. Anal. 15,
  282-290 (1978), [DOI
  10.1137/0715019](https://doi.org/10.1137/0715019).

## Related work

See [`SimplexQuad.jl`](https://github.com/eschnett/SimplexQuad.jl) for
a simplex quadrature package that uses a different algorithm.
