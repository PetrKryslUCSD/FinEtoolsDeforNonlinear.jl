var documenterSearchIndex = {"docs":
[{"location":"guide/guide.html#Guide-1","page":"Guide","title":"Guide","text":"","category":"section"},{"location":"guide/guide.html#Modules-1","page":"Guide","title":"Modules","text":"","category":"section"},{"location":"guide/guide.html#","page":"Guide","title":"Guide","text":"Nonlinear deformation:  AlgoDeforNonlinearModule (algorithms),  FEMMDeforNonlinearModule,  MatDeforNonlinearModule, MatDeforNeohookeanModule (hyperelastic material models).","category":"page"},{"location":"guide/guide.html#Nonlinear-deformation-FEM-machines-1","page":"Guide","title":"Nonlinear deformation FEM  machines","text":"","category":"section"},{"location":"guide/guide.html#","page":"Guide","title":"Guide","text":"For  the base machine for nonlinear deformation, FEMMDeforNonlinearBase, assumes standard isoparametric  finite elements. It evaluates  the interior integrals:","category":"page"},{"location":"guide/guide.html#","page":"Guide","title":"Guide","text":"The stiffness matrix, the geometric stiffness matrix, the loading corresponding to prescribed nonzero displacement.\nThe load vector corresponding to restoring forces due to deformation of the material.","category":"page"},{"location":"guide/guide.html#","page":"Guide","title":"Guide","text":"The nonlinear deformation package FinEtoolsDeforNonlinear is dependent on the linear deformation package, FinEtoolsDeforLinear.","category":"page"},{"location":"guide/guide.html#Material-and-Material-Orientation-1","page":"Guide","title":"Material and Material Orientation","text":"","category":"section"},{"location":"guide/guide.html#Materials-for-nonlinear-deformation-analysis-1","page":"Guide","title":"Materials for nonlinear deformation analysis","text":"","category":"section"},{"location":"guide/guide.html#","page":"Guide","title":"Guide","text":"The module MatDeforNonlinearModule defines the interface  for nonlinear materials (based upon the model-reduction type, 3-D, 2-D, or 1-D) are provided.","category":"page"},{"location":"guide/guide.html#","page":"Guide","title":"Guide","text":"Currently  there are material types for hyper elastic materials. The user may add  additional material types by deriving from AbstractMatDeforNonlinear and equipping them with two methods: (1) compute the tangent moduli, and (2) update the material state.","category":"page"},{"location":"guide/guide.html#","page":"Guide","title":"Guide","text":"The material model assumes that the strains are provided in the local, material-attached coordinate system. The stresses (and the tangent moduli) also need to be output in this local coordinate system.","category":"page"},{"location":"guide/guide.html#","page":"Guide","title":"Guide","text":"For full generality, material types  should implement these methods for fully three-dimensional, plane strain and plane stress, 2D axially symmetric, and one-dimensional deformation models.","category":"page"},{"location":"guide/guide.html#Algorithms-1","page":"Guide","title":"Algorithms","text":"","category":"section"},{"location":"guide/guide.html#Nonlinear-deformation-algorithms-1","page":"Guide","title":"Nonlinear deformation algorithms","text":"","category":"section"},{"location":"guide/guide.html#","page":"Guide","title":"Guide","text":"There are algorithms for","category":"page"},{"location":"guide/guide.html#","page":"Guide","title":"Guide","text":"Nonlinear static analysis.","category":"page"},{"location":"index.html#FinEtools-(Finite-Element-tools)-Documentation-1","page":"Home","title":"FinEtools (Finite Element tools) Documentation","text":"","category":"section"},{"location":"index.html#Conceptual-guide-1","page":"Home","title":"Conceptual guide","text":"","category":"section"},{"location":"index.html#","page":"Home","title":"Home","text":"The construction of the toolkit is described: the composition of modules, the basic data structures, the methodology of computing quantities required in the finite element methodology, and more.","category":"page"},{"location":"index.html#","page":"Home","title":"Home","text":"Pages = [\n    \"guide/guide.md\",\n]\nDepth = 1","category":"page"},{"location":"index.html#Manual-1","page":"Home","title":"Manual","text":"","category":"section"},{"location":"index.html#","page":"Home","title":"Home","text":"The description of the types and the functions, organized by module and/or other logical principle.","category":"page"},{"location":"index.html#","page":"Home","title":"Home","text":"Pages = [\n    \"man/types.md\",\n    \"man/functions.md\",\n]\nDepth = 2","category":"page"},{"location":"man/types.html#Types-1","page":"Types","title":"Types","text":"","category":"section"},{"location":"man/types.html#FEM-machines-1","page":"Types","title":"FEM machines","text":"","category":"section"},{"location":"man/types.html#Nonlinear-deformation-1","page":"Types","title":"Nonlinear deformation","text":"","category":"section"},{"location":"man/types.html#","page":"Types","title":"Types","text":"Modules = [FinEtools, FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule, FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule]\nPrivate = true\nOrder = [:type]","category":"page"},{"location":"man/types.html#FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule.AbstractFEMMDeforNonlinear","page":"Types","title":"FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule.AbstractFEMMDeforNonlinear","text":"AbstractFEMMDeforNonlinear <: AbstractFEMMDeforLinear\n\nAbstract type of FEMM for nonlinear deformation.\n\n\n\n\n\n","category":"type"},{"location":"man/types.html#FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule.FEMMDeforNonlinear","page":"Types","title":"FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule.FEMMDeforNonlinear","text":"FEMMDeforNonlinear{MR<:AbstractDeforModelRed,  S<:AbstractFESet, F<:Function, M<:AbstractMatDeforNonlinear} <: AbstractFEMMDeforNonlinear\n\nClass for nonlinear deformation finite element modeling machine.\n\n\n\n\n\n","category":"type"},{"location":"man/types.html#FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule.FEMMDeforNonlinear-Union{Tuple{M}, Tuple{F}, Tuple{S}, Tuple{MR}, Tuple{Type{MR},IntegDomain{S,F},CSys,M}} where M<:FinEtoolsDeforNonlinear.MatDeforNonlinearModule.AbstractMatDeforNonlinear where F<:Function where S<:AbstractFESet where MR<:AbstractDeforModelRed","page":"Types","title":"FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule.FEMMDeforNonlinear","text":"FEMMDeforNonlinear(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M) where {MR<:AbstractDeforModelRed, S<:AbstractFESet, F<:Function, M<:AbstractMatDeforNonlinear}\n\nConstructor of nonlinear deformation finite element modeling machine.\n\n\n\n\n\n","category":"method"},{"location":"man/types.html#FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule.FEMMDeforNonlinear-Union{Tuple{M}, Tuple{F}, Tuple{S}, Tuple{MR}, Tuple{Type{MR},IntegDomain{S,F},M}} where M<:FinEtoolsDeforNonlinear.MatDeforNonlinearModule.AbstractMatDeforNonlinear where F<:Function where S<:AbstractFESet where MR<:AbstractDeforModelRed","page":"Types","title":"FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule.FEMMDeforNonlinear","text":"FEMMDeforNonlinear(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M) where {MR<:AbstractDeforModelRed, S<:AbstractFESet, F<:Function, M<:AbstractMatDeforNonlinear}\n\nConstructor of nonlinear deformation finite element modeling machine.\n\nThe material coordinate system is the default.\n\n\n\n\n\n","category":"method"},{"location":"man/types.html#Material-models-1","page":"Types","title":"Material models","text":"","category":"section"},{"location":"man/types.html#Material-models-for-nonlinear-deformation-1","page":"Types","title":"Material models for nonlinear deformation","text":"","category":"section"},{"location":"man/types.html#","page":"Types","title":"Types","text":"Modules = [FinEtools, FinEtoolsDeforNonlinear.MatDeforNonlinearModule, ]\nPrivate = true\nOrder = [:type]","category":"page"},{"location":"man/types.html#FinEtoolsDeforNonlinear.MatDeforNonlinearModule.AbstractMatDeforNonlinear","page":"Types","title":"FinEtoolsDeforNonlinear.MatDeforNonlinearModule.AbstractMatDeforNonlinear","text":"AbstractMatDeforNonlinear <: AbstractMatDefor\n\nAbstract nonlinear material.\n\n\n\n\n\n","category":"type"},{"location":"man/types.html#Material-models-for-neohookean-hyperelasticity-1","page":"Types","title":"Material models for neohookean hyperelasticity","text":"","category":"section"},{"location":"man/types.html#","page":"Types","title":"Types","text":"Modules = [FinEtools,  FinEtoolsDeforNonlinear.MatDeforNeohookeanModule,]\nPrivate = true\nOrder = [:type]","category":"page"},{"location":"man/types.html#FinEtoolsDeforNonlinear.MatDeforNeohookeanModule.MatDeforNeohookean","page":"Types","title":"FinEtoolsDeforNonlinear.MatDeforNeohookeanModule.MatDeforNeohookean","text":"MatDeforNeohookean{MR<:AbstractDeforModelRed, MTAN<:Function, MUPD<:Function} <: AbstractMatDeforNonlinear\n\nType for triaxial neohookean hyperelastic material.\n\nnote: Note\n\n\nThe material object is not thread safe. It holds temporary arrays. If the object is to be used in a multi-threaded environment, each thread must have its own private copy.\n\nnote: Note\n\n\nJulia threads do not cooperate with BLAS threads. All the BLAS calls are eliminated from the implementation of this material in order to obtain good parallel efficiency.\n\n\n\n\n\n","category":"type"},{"location":"man/types.html#FinEtoolsDeforNonlinear.MatDeforNeohookeanModule.MatDeforNeohookean-Tuple{Type{DeforModelRed3D},Float64,Float64,Float64}","page":"Types","title":"FinEtoolsDeforNonlinear.MatDeforNeohookeanModule.MatDeforNeohookean","text":"MatDeforNeohookean(mr::Type{DeforModelRed3D}, mass_density::FFlt, E::FFlt, nu::FFlt)\n\nCreate triaxial neohookean hyperelastic material.\n\n\n\n\n\n","category":"method"},{"location":"man/types.html#FinEtoolsDeforNonlinear.MatDeforNeohookeanModule.MatDeforNeohookean-Union{Tuple{MR}, Tuple{Type{MR},Float64,Float64}} where MR","page":"Types","title":"FinEtoolsDeforNonlinear.MatDeforNeohookeanModule.MatDeforNeohookean","text":"    MatDeforNeohookean(mr::Type{MR}, E::FFlt, nu::FFlt) where {MR}\n\nCreate neohookean isotropic elastic material.\n\nThe mass density is the default value.\n\n\n\n\n\n","category":"method"},{"location":"man/types.html#Material-models-for-St-Venant-Kirchhoff-hyperelasticity-1","page":"Types","title":"Material models for St Venant-Kirchhoff hyperelasticity","text":"","category":"section"},{"location":"man/types.html#","page":"Types","title":"Types","text":"Modules = [FinEtools,  FinEtoolsDeforNonlinear.MatDeforStVKModule,]\nPrivate = true\nOrder = [:type]","category":"page"},{"location":"man/types.html#FinEtoolsDeforNonlinear.MatDeforStVKModule.MatDeforStVK","page":"Types","title":"FinEtoolsDeforNonlinear.MatDeforStVKModule.MatDeforStVK","text":"MatDeforStVK{MR<:AbstractDeforModelRed, MTAN<:Function, MUPD<:Function} <: AbstractMatDeforNonlinear\n\nType for triaxial St Venant-Kirchhoff hyperelastic material.\n\n\n\n\n\n","category":"type"},{"location":"man/types.html#FinEtoolsDeforNonlinear.MatDeforStVKModule.MatDeforStVK-Tuple{Type{DeforModelRed3D},Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64}","page":"Types","title":"FinEtoolsDeforNonlinear.MatDeforStVKModule.MatDeforStVK","text":"MatDeforStVK(mr::Type{DeforModelRed3D}, mass_density::FFlt, E1::FFlt, E2::FFlt, E3::FFlt, G12::FFlt, G13::FFlt, G23::FFlt, nu12::FFlt, nu13::FFlt, nu23::FFlt)\n\nCreate triaxial St Venant-Kirchhoff hyperelastic material.\n\nIn general, the material is assumed to be orthotropic. There is a specialized constructor for an isotropic version.\n\n\n\n\n\n","category":"method"},{"location":"man/types.html#FinEtoolsDeforNonlinear.MatDeforStVKModule.MatDeforStVK-Union{Tuple{MR}, Tuple{Type{MR},Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64}} where MR","page":"Types","title":"FinEtoolsDeforNonlinear.MatDeforStVKModule.MatDeforStVK","text":"MatDeforStVK(mr::Type{MR}, E::FFlt, nu::FFlt) where {MR}\n\nCreate St Venant-Kirchhoff isotropic elastic material.\n\nThe mass density is the default value.\n\n\n\n\n\n","category":"method"},{"location":"man/types.html#FinEtoolsDeforNonlinear.MatDeforStVKModule.MatDeforStVK-Union{Tuple{MR}, Tuple{Type{MR},Float64,Float64}} where MR","page":"Types","title":"FinEtoolsDeforNonlinear.MatDeforStVKModule.MatDeforStVK","text":"MatDeforStVK(mr::Type{MR}, E::FFlt, nu::FFlt) where {MR}\n\nCreate St Venant-Kirchhoff isotropic elastic material.\n\nThe mass density is the default value.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#Functions-1","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"man/functions.html#Helpers-1","page":"Functions","title":"Helpers","text":"","category":"section"},{"location":"man/functions.html#","page":"Functions","title":"Functions","text":"Modules = [FinEtools, FinEtoolsDeforNonlinear.AssemblyModule,  ]\nPrivate = true\nOrder = [:function]","category":"page"},{"location":"man/functions.html#FinEtools.AssemblyModule.assemble!-Union{Tuple{D}, Tuple{MV}, Tuple{T}, Tuple{FinEtoolsDeforNonlinear.AssemblyModule.SysvecAssemblerOpt{T},MV,D}} where D<:(AbstractArray{Int64,N} where N) where MV<:(AbstractArray{T,N} where N) where T<:Number","page":"Functions","title":"FinEtools.AssemblyModule.assemble!","text":"assemble!(self::SysvecAssemblerOpt{T}, vec::MV,\n  dofnums::D) where {T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}\n\nAssemble an elementwise vector.\n\nThe method assembles a column element vector using the vector of degree of freedom numbers for the rows.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtools.AssemblyModule.makevector!-Tuple{FinEtoolsDeforNonlinear.AssemblyModule.SysvecAssemblerOpt}","page":"Functions","title":"FinEtools.AssemblyModule.makevector!","text":"makevector!(self::SysvecAssemblerOpt)\n\nMake the global vector. Beware: the buffer itself is returned as the vector, not a copy.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtools.AssemblyModule.startassembly!-Union{Tuple{T}, Tuple{FinEtoolsDeforNonlinear.AssemblyModule.SysvecAssemblerOpt{T},Int64}} where T<:Number","page":"Functions","title":"FinEtools.AssemblyModule.startassembly!","text":"startassembly!(self::SysvecAssemblerOpt{T},\n  ndofs_row::FInt) where {T<:Number}\n\nStart assembly.\n\nThe method makes the buffer for the vector assembly. It must be called before the first call to the method assemble.\n\nndofs_row= Total number of degrees of freedom.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FEM-machines-1","page":"Functions","title":"FEM machines","text":"","category":"section"},{"location":"man/functions.html#Nonlinear-deformation-1","page":"Functions","title":"Nonlinear deformation","text":"","category":"section"},{"location":"man/functions.html#Simple-FE-models-1","page":"Functions","title":"Simple FE models","text":"","category":"section"},{"location":"man/functions.html#","page":"Functions","title":"Functions","text":"Modules = [FinEtools, FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule, FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule ]\nPrivate = true\nOrder = [:function]","category":"page"},{"location":"man/functions.html#FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule.geostiffness-Union{Tuple{T}, Tuple{A}, Tuple{FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule.AbstractFEMMDeforNonlinear,A,NodalField{Float64},NodalField{T},NodalField{T},Float64,Float64}} where T<:Number where A<:AbstractSysmatAssembler","page":"Functions","title":"FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule.geostiffness","text":"geostiffness(self::AbstractFEMMDeforNonlinear, assembler::A, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, tn::FFlt, dtn::FFlt) where {A<:AbstractSysmatAssembler, T<:Number}\n\nCompute and assemble geometric stiffness matrix.\n\nArguments\n\nassembler = matrix assembler,\ngeom = geometry field: reference coordinates of the nodes,\nun1 = displacement field at time tn1 = tn + dtn,\nun = displacement field at time tn,\ntn = time in step n\ndtn = increment of time\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule.nzebcloads-Union{Tuple{T}, Tuple{A}, Tuple{FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule.AbstractFEMMDeforNonlinear,A,NodalField{Float64},NodalField{T},NodalField{T},NodalField{T},Float64,Float64}} where T<:Number where A<:SysvecAssembler","page":"Functions","title":"FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule.nzebcloads","text":"nzebcloads(self::AbstractFEMMDeforNonlinear, assembler::A, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, du::NodalField{T}, tn::FFlt, dtn::FFlt) where {A<:SysvecAssembler, T<:Number}\n\nCompute and assemble load vector due to prescribed increment of displacements.\n\nArguments\n\nassembler = matrix assembler,\ngeom = geometry field: reference coordinates of the nodes,\nun1 = displacement field at time tn1 = tn + dtn,\nun = displacement field at time tn,\ndu = field of displacement increment prescribed at time tn1 = tn + dtn,   The increment is stored in the fixed values of this field.\ntn = time in step n\ndtn = increment of time\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule.restoringforce-Union{Tuple{T}, Tuple{A}, Tuple{FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule.AbstractFEMMDeforNonlinear,A,NodalField{Float64},NodalField{T},NodalField{T},Float64,Float64}, Tuple{FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule.AbstractFEMMDeforNonlinear,A,NodalField{Float64},NodalField{T},NodalField{T},Float64,Float64,Any}} where T<:Number where A<:AbstractSysvecAssembler","page":"Functions","title":"FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule.restoringforce","text":"restoringforce(self::AbstractFEMMDeforNonlinear, assembler::A, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, tn::FFlt, dtn::FFlt, savestate = false) where {A<:AbstractSysvecAssembler, T<:Number}\n\nCompute the restoring force vector.\n\nNote: This method UPDATES the state of the FEMM object.  In particular, the material state gets updated.\n\nArguments\n\nassembler = vector assembler,\ngeom = geometry field: reference coordinates of the nodes,\nun1 = displacement field at time tn1 = tn + dtn,\nun = displacement field at time tn,\ntn = time in step n\ndtn = increment of time\nsavestate = bool flag: should we modify the material states (savestate = true)? Otherwise work with a copy of the material state.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule.stiffness-Union{Tuple{T}, Tuple{A}, Tuple{FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule.AbstractFEMMDeforNonlinear,A,NodalField{Float64},NodalField{T},NodalField{T},Float64,Float64}} where T<:Number where A<:AbstractSysmatAssembler","page":"Functions","title":"FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule.stiffness","text":"stiffness(self::AbstractFEMMDeforNonlinear, assembler::A, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, tn::FFlt, dtn::FFlt) where {A<:AbstractSysmatAssembler, T<:Number}\n\nCompute and assemble  stiffness matrix.\n\nArguments\n\nassembler = matrix assembler,\ngeom = geometry field: reference coordinates of the nodes,\nun1 = displacement field at time tn1 = tn + dtn,\nun = displacement field at time tn,\ntn = time in step n\ndtn = increment of time\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#Algorithms-1","page":"Functions","title":"Algorithms","text":"","category":"section"},{"location":"man/functions.html#Nonlinear-deformation-2","page":"Functions","title":"Nonlinear deformation","text":"","category":"section"},{"location":"man/functions.html#","page":"Functions","title":"Functions","text":"Modules = [FinEtools, FinEtoolsDeforNonlinear.AlgoDeforNonlinearModule]\nPrivate = true\nOrder = [:function]","category":"page"},{"location":"man/functions.html#FinEtoolsDeforNonlinear.AlgoDeforNonlinearModule.nonlinearstatics-Tuple{Dict{String,Any}}","page":"Functions","title":"FinEtoolsDeforNonlinear.AlgoDeforNonlinearModule.nonlinearstatics","text":"AlgoDeforNonlinearModule.nonlinearstatics(modeldata::FDataDict)\n\nAlgorithm for static nonlinear deformation (stress) analysis.\n\nThe algorithm chooses steps from the array of load multipliers: the step takes it precisely from the preceding step to the next step in one go.\n\nArgument\n\nmodeldata = dictionary with values for keys\n\n\"fens\"  = finite element node set\n\"regions\"  = array of region dictionaries\n\"essential_bcs\" = array of essential boundary condition dictionaries\n\"traction_bcs\" = array of traction boundary condition dictionaries\n\"temperature_change\" = dictionary of data for temperature change\n\nFor each region (connected piece of the domain made of a particular material), mandatory, the  region dictionary  contains values for keys:\n\n\"femm\" = finite element model machine (mandatory);\n\nFor essential boundary conditions (optional) each dictionary would hold\n\n\"displacement\" = when this key is not present, the assumption is     that the displacement is fixed at zero (0). Otherwise, this needs to     be set to a function with signature f(x, lambda), where x is the     location of the node in the reference configuration, and lambda is the     load factor. In other words, whenever this quantity is supplied, it is     implied that the displacement depends on the load factor.\n\"component\" = which component is prescribed  (1, 2, 3)?\n\"node_list\" = list of nodes on the boundary to which the condition     applies (mandatory)\n\nFor traction boundary conditions (optional) each dictionary would hold\n\n\"femm\" = finite element model machine (mandatory);\n\"traction_vector\" = traction vector, a force-intensity (ForceIntensity) object.\n\nControl parameters The following attributes  may be supplied:\n\n\"load_multipliers\" = For what load multipliers should the solution be         calculated? Array of monotonically increasing numbers.\n\"line_search\" = Should we use line search? Boolean.  Default = true.\n\"maxdu_tol\" = Tolerance on the magnitude  of the largest incremental         displacement component.\n\"maxbal_tol\" = Tolerance on the magnitude of the largest out-of-balance         force component.\n\"iteration_observer\" = observer function to be called     after each iteration.  Default is to do nothing.     The observer function has a signature               iteration_observer(lambda,iter,du,modeldata)     where lambda is the current load factor, iter is the iteration     number, du is the nodal field of current displacement increments.\n\"increment_observer\" = observer function     to be called after convergence is reached in each step (optional)     The observer function has a signature               output(lambda, modeldata)     where lambda is the current load factor. Default is to do nothing.     The increment observer can refer to the following key-value pairs in modeldata:     - \"un1\": converged displacement in current step     - \"un\": converged displacement in the last step     - \"t\": current value of load factor     - \"dt\": increment of the load factor\n\nOutput\n\nmodeldata = the dictionary on input is augmented with the keys\n\n\"geom\" = the nodal field that is the geometry\n\"u\" = the nodal field that is the computed displacement\n\"reactions\" = computed reaction nodal field\n\"timing\" = dictionary with timing results\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#Material-models-1","page":"Functions","title":"Material models","text":"","category":"section"},{"location":"man/functions.html#Material-models-for-nonlinear-deformation-1","page":"Functions","title":"Material models for nonlinear deformation","text":"","category":"section"},{"location":"man/functions.html#","page":"Functions","title":"Functions","text":"Modules = [FinEtools, FinEtoolsDeforNonlinear.MatDeforNonlinearModule, ]\nPrivate = true\nOrder = [:function]","category":"page"},{"location":"man/functions.html#FinEtoolsDeforNonlinear.MatDeforNonlinearModule.estimatesoundspeed-Union{Tuple{M}, Tuple{M}} where M<:FinEtoolsDeforNonlinear.MatDeforNonlinearModule.AbstractMatDeforNonlinear","page":"Functions","title":"FinEtoolsDeforNonlinear.MatDeforNonlinearModule.estimatesoundspeed","text":"estimatesoundspeed(self::M)  where {M<:AbstractMatDeforNonlinear}\n\nEstimate sound speed in the undeformed state.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsDeforNonlinear.MatDeforNonlinearModule.newstate-Union{Tuple{M}, Tuple{M}} where M<:FinEtoolsDeforNonlinear.MatDeforNonlinearModule.AbstractMatDeforNonlinear","page":"Functions","title":"FinEtoolsDeforNonlinear.MatDeforNonlinearModule.newstate","text":"newstate(self::M) where {M<:AbstractMatDeforNonlinear}\n\nCreate an initial material state at an integration point.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsDeforNonlinear.MatDeforNonlinearModule.tangentmoduli!-Union{Tuple{M}, Tuple{M,Array{Float64,2},Array{Float64,1},Array{Float64,2},Array{Float64,2},Float64,Float64,Array{Float64,2},Int64}} where M<:FinEtoolsDeforNonlinear.MatDeforNonlinearModule.AbstractMatDeforNonlinear","page":"Functions","title":"FinEtoolsDeforNonlinear.MatDeforNonlinearModule.tangentmoduli!","text":"tangentmoduli!(self::M, D::FFltMat, statev::FFltVec, Fn1::FFltMat,\n    Fn::FFltMat, tn::FFlt, dtn::FFlt, loc::FFltMat, label::FInt)\n    where {M<:AbstractMatDeforNonlinear}\n\nCalculate the material stiffness matrix.\n\nArguments\n\nself = material\nD = matrix of tangent moduli, supplied as a buffer and overwritten. Returned as output.\nstatev = material state vector, the content of this vector must not change   inside this function.\nFn1 = deformation gradient at time tn1 = tn + dtn,\nFn = deformation gradient at time tn,\ntn = time in step n\ndtn = increment of time\nloc = location of the integration point in the reference coordinates (time t0),\nlabel = label of the element containing the integration point\n\nnote: Note\n\n\nThe deformation gradients and the matrix of the tangent moduli are expressed with respect to the local material coordinate system.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsDeforNonlinear.MatDeforNonlinearModule.totlag2curr!-Tuple{Any,Any,Any}","page":"Functions","title":"FinEtoolsDeforNonlinear.MatDeforNonlinearModule.totlag2curr!","text":"totlag2curr!(c, C, F)\n\nConvert a total Lagrangean constitutive matrix to a current Lagrangean one (sometimes known as \"Eulerian\").\n\nC    = Lagrangean constitutive matrix, 6x6, symmetric\nF    = current deformation gradient, FiJ = partial xi / partial X_J\n\nThe transformation is cijkl = 1/J CIJKL FiI FjJ FkK FlL. In the present case the fourth-order tensor is represented with a 6 x 6 matrix.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsDeforNonlinear.MatDeforNonlinearModule.totlag2curr4th!-Tuple{Any,Any,Any}","page":"Functions","title":"FinEtoolsDeforNonlinear.MatDeforNonlinearModule.totlag2curr4th!","text":"totlag2curr4th!(c, C, F)\n\nConvert a total Lagrangean constitutive matrix to a current Lagrangean one (sometimes known as \"Eulerian\").\n\nC    = Lagrangean constitutive matrix, fourth-order tensor\nF    = current deformation gradient, FiJ = partial xi / partial X_J\n\nThe transformation is cijkl = 1/J CIJKL FiI FjJ FkK FlL. Both the input and the output are fourth-order tensors.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsDeforNonlinear.MatDeforNonlinearModule.totlag2currsymm!-Tuple{Any,Any,Any}","page":"Functions","title":"FinEtoolsDeforNonlinear.MatDeforNonlinearModule.totlag2currsymm!","text":"totlag2currsymm!(c, C, F)\n\nConvert a total Lagrangean constitutive matrix to a current Lagrangean one (sometimes known as \"Eulerian\").\n\nC    = Lagrangean constitutive matrix, 6x6, symmetric\nF    = current deformation gradient, FiJ = partial xi / partial X_J\n\nThe transformation is cijkl = 1/J CIJKL FiI FjJ FkK FlL. In the present case the fourth-order tensor is represented with a 6 x 6 matrix.\n\nnote: Note\n\n\nThe Lagrangean material stiffness matrices, both input and output, are presumed symmetric.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsDeforNonlinear.MatDeforNonlinearModule.update!-Union{Tuple{M}, Tuple{M,Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,2},Array{Float64,2},Float64,Float64}, Tuple{M,Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,2},Array{Float64,2},Float64,Float64,Array{Float64,2}}, Tuple{M,Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,2},Array{Float64,2},Float64,Float64,Array{Float64,2},Int64}, Tuple{M,Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,2},Array{Float64,2},Float64,Float64,Array{Float64,2},Int64,Any}} where M<:FinEtoolsDeforNonlinear.MatDeforNonlinearModule.AbstractMatDeforNonlinear","page":"Functions","title":"FinEtoolsDeforNonlinear.MatDeforNonlinearModule.update!","text":"update!(self::M, statev::FFltVec, stress::FFltVec, output::FFltVec,\n    Fn1::FFltMat, Fn::FFltMat, tn::FFlt, dtn::FFlt,\n    loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)\n    where {M<:AbstractMatDeforNonlinear}\n\nUpdate material state.\n\nArguments\n\nself = material\nstatev = state variables: array which is (if necessary) allocated  in an appropriate    size, filled with the state variables, and returned. The contents of this    vector may change as the state of the material is updated by the logic inside    this function. If this change is to be saved, it must happen outside of this    function.\ncauchy = Cauchy stress vector, allocated by the caller with a size of the number of stress and strain components, nstressstrain. The components of the stress vector are calculated and stored in the stress vector.\noutput =  array which is (if necessary) allocated  in an appropriate size,    filled with the output quantity, and returned.\nFn1 = deformation gradient at time tn1 = tn + dtn,\nFn = deformation gradient at time tn,\ntn = time in step n\ndtn = increment of time\nloc = location of the integration point in the reference coordinates (time t0),\nlabel = label of the element containing the integration point\n\nOutput\n\ncauchy = Cauchy stress vector\noutput = output array\n\nnote: Note\n\n\nThe deformation gradients and the stress vector are expressed with respect to the local material coordinate system.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#Material-models-for-neohookean-hyperelasticity-1","page":"Functions","title":"Material models for neohookean hyperelasticity","text":"","category":"section"},{"location":"man/functions.html#","page":"Functions","title":"Functions","text":"Modules = [FinEtools,  FinEtoolsDeforNonlinear.MatDeforNeohookeanModule, FinEtoolsDeforNonlinear.MatDeforNeohookeanADModule, FinEtoolsDeforNonlinear.MatDeforNeohookeanNaiveModule]\nPrivate = true\nOrder = [:function]","category":"page"},{"location":"man/functions.html#FinEtoolsDeforNonlinear.MatDeforNonlinearModule.estimatesoundspeed-Union{Tuple{M}, Tuple{M}} where M<:FinEtoolsDeforNonlinear.MatDeforNeohookeanModule.MatDeforNeohookean","page":"Functions","title":"FinEtoolsDeforNonlinear.MatDeforNonlinearModule.estimatesoundspeed","text":"estimatesoundspeed(self::M)  where {M<:AbstractMatDeforNonlinear}\n\nEstimate sound speed in the undeformed state.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsDeforNonlinear.MatDeforNonlinearModule.estimatesoundspeed-Union{Tuple{M}, Tuple{M}} where M<:FinEtoolsDeforNonlinear.MatDeforNeohookeanNaiveModule.MatDeforNeohookeanNaive","page":"Functions","title":"FinEtoolsDeforNonlinear.MatDeforNonlinearModule.estimatesoundspeed","text":"estimatesoundspeed(self::M)  where {M<:AbstractMatDeforNonlinear}\n\nEstimate sound speed in the undeformed state.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#Material-models-for-St-Venant-Kirchhoff-hyperelasticity-1","page":"Functions","title":"Material models for St Venant-Kirchhoff hyperelasticity","text":"","category":"section"},{"location":"man/functions.html#","page":"Functions","title":"Functions","text":"Modules = [FinEtools,  FinEtoolsDeforNonlinear.MatDeforStVKModule, FinEtoolsDeforNonlinear.MatDeforStVKADModule]\nPrivate = true\nOrder = [:function]","category":"page"},{"location":"man/functions.html#Material-models-for-Mooney-Rivlin-hyperelasticity-1","page":"Functions","title":"Material models for Mooney-Rivlin hyperelasticity","text":"","category":"section"},{"location":"man/functions.html#","page":"Functions","title":"Functions","text":"Modules = [FinEtools,  FinEtoolsDeforNonlinear.MatDeforMooneyRivlinADModule, FinEtoolsDeforNonlinear.MatDeforI1RivlinADModule]\nPrivate = true\nOrder = [:function]","category":"page"}]
}
