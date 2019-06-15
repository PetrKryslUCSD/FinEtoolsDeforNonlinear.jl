# Guide

## Modules

- Nonlinear deformation:  `AlgoDeforNonlinearModule` (algorithms),  `FEMMDeforNonlinearModule`,  `MatDeforNonlinearModule`, `MatDeforNeohookeanModule` (hyperelastic material models).

## Nonlinear deformation FEM  machines

For  the base machine for linear deformation, `FEMMDeforNonlinearBase`, assumes standard isoparametric  finite elements. It evaluates  the interior integrals:

- The stiffness matrix, the geometric stiffness matrix, the loading corresponding to prescribed nonzero displacement.

- The load vector corresponding to resisting forces due to deformation of the material.

## Material and Material Orientation

### Materials for nonlinear deformation analysis

The module `MatDeforNonlinearModule` defines the interface  for nonlinear materials (based upon the model-reduction type, 3-D, 2-D, or 1-D) are provided.

Currently  there are material types for hyper elastic materials. The user may add  additional material types by deriving from `AbstractMatDeforNonlinear` and equipping them with two methods: (1) compute the tangent moduli, and (2) update the material state.

For full generality, material types  should implement these methods for fully three-dimensional, plane strain and plane stress, 2D axially symmetric, and one-dimensional deformation models.

## Algorithms

### Nonlinear deformation algorithms

There are algorithms for

- Nonlinear static analysis.
