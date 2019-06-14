[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.com/PetrKryslUCSD/FinEtoolsDeforNonlinear.jl.svg?branch=master)](https://travis-ci.com/PetrKryslUCSD/FinEtoolsDeforNonlinear.jl)

[![][docs-latest-img]][docs-latest-url]

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: http://petrkryslucsd.github.io/FinEtoolsDeforNonlinear.jl/latest/

# FinEtoolsDeforNonlinear: Nonlinear stress analysis application

`FinEtools` is a package for basic operations on finite element meshes.
`FinEtoolsDeforNonlinear` is a package using `FinEtools` to solve nonlinear stress analysis problems.
At the moment,  statics and hyper elastic materials are included.

## News

- 06/11/2019: Applications are now separated  out from the `FinEtools` package.

[Past news](oldnews.md)

## How to run

The [FinEtools](https://github.com/PetrKryslUCSD/FinEtools.jl) package is
needed. The entire setup of `FinEtoolsDeforNonlinear` can be performed with
```julia
] activate .; instantiate
```

The package `FinEtoolsDeforNonlinear` can be tested as
```julia
] activate .; instantiate; test
```

There are a number of examples covering statics and dynamics. The examples may
be executed as described in the  [conceptual guide to
`FinEtools`](https://petrkryslucsd.github.io/FinEtools.jl/latest).
