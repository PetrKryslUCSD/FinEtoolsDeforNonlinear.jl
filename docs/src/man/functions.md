# Functions

## Helpers


```@autodocs
Modules = [FinEtools, FinEtoolsDeforNonlinear.AssemblyModule,  ]
Private = true
Order = [:function]
```


## FEM machines

### Nonlinear deformation

#### Simple FE models

```@autodocs
Modules = [FinEtools, FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule, FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule ]
Private = true
Order = [:function]
```

## Algorithms

### Nonlinear deformation

```@autodocs
Modules = [FinEtools, FinEtoolsDeforNonlinear.AlgoDeforNonlinearModule]
Private = true
Order = [:function]
```

## Material models

### Material models for nonlinear deformation

```@autodocs
Modules = [FinEtools, FinEtoolsDeforNonlinear.MatDeforNonlinearModule, ]
Private = true
Order = [:function]
```

### Material models for neohookean hyperelasticity

```@autodocs
Modules = [FinEtools,  FinEtoolsDeforNonlinear.MatDeforNeohookeanModule, FinEtoolsDeforNonlinear.MatDeforNeohookeanADModule, FinEtoolsDeforNonlinear.MatDeforNeohookeanNaiveModule]
Private = true
Order = [:function]
```

### Material models for St Venant-Kirchhoff hyperelasticity

```@autodocs
Modules = [FinEtools,  FinEtoolsDeforNonlinear.MatDeforStVKModule, FinEtoolsDeforNonlinear.MatDeforStVKADModule]
Private = true
Order = [:function]
```

### Material models for Mooney-Rivlin hyperelasticity

```@autodocs
Modules = [FinEtools,  FinEtoolsDeforNonlinear.MatDeforMooneyRivlinADModule, FinEtoolsDeforNonlinear.MatDeforI1RivlinADModule]
Private = true
Order = [:function]
```
