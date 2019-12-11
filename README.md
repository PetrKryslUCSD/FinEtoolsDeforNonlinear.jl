[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.com/PetrKryslUCSD/FinEtoolsDeforNonlinear.jl.svg?branch=master)](https://travis-ci.com/PetrKryslUCSD/FinEtoolsDeforNonlinear.jl)

[![][docs-latest-img]][docs-latest-url]

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: http://petrkryslucsd.github.io/FinEtoolsDeforNonlinear.jl/latest/

# FinEtoolsDeforNonlinear: Nonlinear stress analysis application

[`FinEtools`](https://github.com/PetrKryslUCSD/FinEtools.jl.git) is a package
for basic operations on finite element meshes. `FinEtoolsDeforNonlinear` is a
package using `FinEtools` to solve nonlinear stress analysis problems. At the
moment,  statics and dynamics with hyperelastic materials are included.

## News

- 12/13/2019: Instrumented an example of transient (explicit) dynamics so that runs in parallel on multiple threads.
- 12/09/2019: Added an example of transient (explicit) dynamics.

[Past news](oldnews.md)

## How to test the package

Here is a record of a session to install this package and test it. You should
see something similar. The git bash running on Windows 10 was used in this
example.

Clone the repo:
```
$ git clone https://github.com/PetrKryslUCSD/FinEtoolsDeforNonlinear.jl.git
Cloning into 'FinEtoolsDeforNonlinear.jl'...
remote: Enumerating objects: 70, done.
remote: Counting objects: 100% (70/70), done.
remote: Compressing objects: 100% (47/47), done.
remote: Total 70 (delta 18), reused 66 (delta 17), pack-reused 0
Unpacking objects: 100% (70/70), done.
```
Change your working directory, and run Julia:
```
$ cd FinEtoolsDeforNonlinear.jl/

PetrKrysl@Spectre MINGW64 /tmp/exp/FinEtoolsDeforNonlinear.jl (master)
$ ~/AppData/Local/Julia-1.2.0-rc1/bin/julia.exe
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.2.0-rc1.0 (2019-05-30)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |
```
Activate and instantiate the environment:
```
(v1.2) pkg> activate .; instantiate
[ Info: activating environment at `C:\Users\PETRKR~1\AppData\Local\Temp\exp\FinEtoolsDeforNonlinear.jl\Project.toml`.
   Cloning default registries into `C:\Users\PetrKrysl\.julia`
   Cloning registry from "https://github.com/JuliaRegistries/General.git"
     Added registry `General` to `C:\Users\PetrKrysl\.julia\registries\General`
   Cloning git-repo `https://github.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl.git`
  Updating git-repo `https://github.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl.git`
   Cloning git-repo `https://github.com/PetrKryslUCSD/FinEtools.jl`
  Updating git-repo `https://github.com/PetrKryslUCSD/FinEtools.jl`
 Installed DefaultApplication ── v0.1.3
 Installed Crayons ───────────── v4.0.0
 Installed OrderedCollections ── v1.1.0
 Installed Arpack ────────────── v0.3.1
 Installed BinaryProvider ────── v0.5.4
 Installed StaticArrays ──────── v0.11.0
 Installed UnicodePlots ──────── v1.1.0
 Installed Compat ────────────── v2.1.0
 Installed Requires ──────────── v0.5.2
 Installed WoodburyMatrices ──── v0.4.1
 Installed Missings ──────────── v0.4.1
 Installed Ratios ────────────── v0.3.1
 Installed Interpolations ────── v0.12.2
 Installed DocStringExtensions ─ v0.7.0
 Installed DataStructures ────── v0.15.0
 Installed Tokenize ──────────── v0.5.4
 Installed MacroTools ────────── v0.5.0
 Installed SortingAlgorithms ─── v0.3.1
 Installed OffsetArrays ──────── v0.11.0
 Installed ArgCheck ──────────── v1.0.1
 Installed AxisAlgorithms ────── v1.0.0
 Installed CSTParser ─────────── v0.6.0
 Installed StatsBase ─────────── v0.30.0
 Installed Parameters ────────── v0.10.3
 Installed PGFPlotsX ─────────── v0.3.8
  Building Arpack ───→ `C:\Users\PetrKrysl\.julia\packages\Arpack\cu5By\deps\build.log`
  Building PGFPlotsX → `C:\Users\PetrKrysl\.julia\packages\PGFPlotsX\PZlVQ\deps\build.log`
```
Test the package:
```
(FinEtoolsDeforNonlinear) pkg> test
   Testing FinEtoolsDeforNonlinear
 Resolving package versions...
Test Summary: | Pass  Total
Materials     |    4      4
 20.068692 seconds (13.37 M allocations: 675.316 MiB, 1.52% gc time)
Test Summary: | Pass  Total
Operations    |   13     13
  8.833191 seconds (24.65 M allocations: 1.213 GiB, 6.02% gc time)
   Testing FinEtoolsDeforNonlinear tests passed
```

## Examples


There are a number of examples covering statics and dynamics. The examples may
be executed as described in the  [conceptual guide to
`FinEtools`](https://petrkryslucsd.github.io/FinEtools.jl/latest).
