"""
Test of the nonlinear FEM machine
"""
module m1testop1
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule: stiffness
using LinearAlgebra: norm
using SparseArrays
using DelimitedFiles
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookean(mr, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    fens, fes = H8block(L, W, H, 1, 1, 1)
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)

    geom = NodalField(fens.xyz)
    un1 = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    un = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    numberdofs!(un1)
    tn, dtn = 0.0, 1.0
    K = stiffness(femm, geom, un1, un, tn, dtn)
    ijs = readdlm("m1testop1.csv", ',', Float64)
    matK  = sparse(Int.(ijs[:, 1]), Int.(ijs[:, 2]), ijs[:, 3], size(K, 1), size(K, 2))
    # @show norm(K - matK) / maximum(matK[:])
    @test norm(K - matK) / maximum(matK[:]) < 1.0e-4
end
end
using .m1testop1
m1testop1.test()

module m1testop2
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule: stiffness
using LinearAlgebra: norm
using SparseArrays
using DelimitedFiles
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookean(mr, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    fens, fes = H8hexahedron([0 0 0; 1.2 0 -0.1; 1.0 0.9 0.1; -0.05 1.1 0; 0 0 1.03; 1.2 0 1.1; 1.0 0.95 0.81; -0.05 1.01 0.9], 1, 1, 1)
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)

    geom = NodalField(fens.xyz)
    un1 = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    un = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    numberdofs!(un1)
    tn, dtn = 0.0, 1.0
    K = stiffness(femm, geom, un1, un, tn, dtn)
    ijs = readdlm("m1testop2.csv", ',', Float64)
    matK  = sparse(Int.(ijs[:, 1]), Int.(ijs[:, 2]), ijs[:, 3], size(K, 1), size(K, 2))
    # @show norm(K - matK) / maximum(matK[:])
    @test norm(K - matK) / maximum(matK[:]) < 1.0e-4
end
end
using .m1testop2
m1testop2.test()

module m1testop3
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule: stiffness
using LinearAlgebra: norm
using SparseArrays
using DelimitedFiles
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookean(mr, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    fens, fes = H8hexahedron([0 0 0; 1.2 0 -0.1; 1.0 0.9 0.1; -0.05 1.1 0; 0 0 1.03; 1.2 0 1.1; 1.0 0.95 0.81; -0.05 1.01 0.9], 1, 1, 1)
    Rmout = fill(0.0, 3, 3)
    rv = vec([1 1.1 0.3])
    rotmat3!(Rmout, rv)
    mcsys = CSys(Rmout)
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), mcsys, m)

    geom = NodalField(fens.xyz)
    un1 = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    un = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    numberdofs!(un1)
    tn, dtn = 0.0, 1.0
    K = stiffness(femm, geom, un1, un, tn, dtn)
    ijs = readdlm("m1testop3.csv", ',', Float64)
    matK  = sparse(Int.(ijs[:, 1]), Int.(ijs[:, 2]), ijs[:, 3], size(K, 1), size(K, 2))
    # @show norm(K - matK) / maximum(matK[:])
    @test norm(K - matK) / maximum(matK[:]) < 1.0e-4
end
end
using .m1testop3
m1testop3.test()

module m1testop4
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule: stiffness
using LinearAlgebra: norm, cross
using SparseArrays
using DelimitedFiles
# using Debugger
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookean(mr, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    fens, fes = H8hexahedron([0 0 0; 1.2 0 -0.1; 1.0 0.9 0.1; -0.05 1.1 0; 0 0 1.03; 1.2 0 1.1; 1.0 0.95 0.81; -0.05 1.01 0.9], 1, 1, 1)
    function  mcsys(csmatout::FFltMat, XYZ, tangents, fe_label)
        # @bp
        n1 = vec(XYZ); n1[1] = 0; n1  = n1 / norm(n1);
        n3 = vec([1, 0, 0]);
        n2 = cross(n3, n1);
        csmatout .= hcat(n1, n2, n3);
    end
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), CSys(3, 3, mcsys), m)

    geom = NodalField(fens.xyz)
    un1 = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    un = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    numberdofs!(un1)
    tn, dtn = 0.0, 1.0
    K = stiffness(femm, geom, un1, un, tn, dtn)
    ijs = readdlm("m1testop4.csv", ',', Float64)
    matK  = sparse(Int.(ijs[:, 1]), Int.(ijs[:, 2]), ijs[:, 3], size(K, 1), size(K, 2))
    # @show norm(K - matK) / maximum(matK[:])
    @test norm(K - matK) / maximum(matK[:]) < 1.0e-4
end
end
using .m1testop4
m1testop4.test()

module m1testop5
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule: stiffness
using LinearAlgebra: norm
using SparseArrays
using DelimitedFiles
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookean(mr, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    fens, fes = H8hexahedron([0 0 0; 1.2 0 -0.1; 1.0 0.9 0.1; -0.05 1.1 0; 0 0 1.03; 1.2 0 1.1; 1.0 0.95 0.81; -0.05 1.01 0.9], 1, 1, 1)
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)

    geom = NodalField(fens.xyz)
    un1 = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    un = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    numberdofs!(un1)
    d(XYZ) = 0.1 * [XYZ[1]^2*XYZ[3], XYZ[2]^2*XYZ[3], XYZ[2]*XYZ[3]-XYZ[3]^3];
    for i = 1:count(fens)
        un1.values[i, :]  .+= d(fens.xyz[i, :]);
    end
    #println("un1.values = $(un1.values)")
    tn, dtn = 0.0, 1.0
    K = stiffness(femm, geom, un1, un, tn, dtn)
    ijs = readdlm("m1testop5.csv", ',', Float64)
    matK  = sparse(Int.(ijs[:, 1]), Int.(ijs[:, 2]), ijs[:, 3], size(K, 1), size(K, 2))
    # @show norm(K - matK) / maximum(matK[:])
    @test norm(K - matK) / maximum(matK[:]) < 1.0e-4
end
end
using .m1testop5
m1testop5.test()

module m1testop6
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule: stiffness
using LinearAlgebra: norm, cross
using SparseArrays
using DelimitedFiles
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookean(mr, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    fens, fes = H8hexahedron([0 0 0; 1.2 0 -0.1; 1.0 0.9 0.1; -0.05 1.1 0; 0 0 1.03; 1.2 0 1.1; 1.0 0.95 0.81; -0.05 1.01 0.9], 1, 1, 1)
    function  mcsys(csmatout::FFltMat, XYZ, tangents, fe_label)
        # @bp
        n1 = vec(XYZ); n1[1] = 0; n1  = n1 / norm(n1);
        n3 = vec([1, 0, 0]);
        n2 = cross(n3, n1);
        csmatout .= hcat(n1, n2, n3);
    end
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), CSys(3, 3, mcsys), m)

    geom = NodalField(fens.xyz)
    un1 = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    un = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    numberdofs!(un1)
    d(XYZ) = 0.1 * [XYZ[1]^2*XYZ[3], XYZ[2]^2*XYZ[3], XYZ[2]*XYZ[3]-XYZ[3]^3];
    for i = 1:count(fens)
        un1.values[i, :]  .+= d(fens.xyz[i, :]);
    end
    #println("un1.values = $(un1.values)")
    tn, dtn = 0.0, 1.0
    K = stiffness(femm, geom, un1, un, tn, dtn)
    ijs = readdlm("m1testop5.csv", ',', Float64)
    matK  = sparse(Int.(ijs[:, 1]), Int.(ijs[:, 2]), ijs[:, 3], size(K, 1), size(K, 2))
    # @show norm(K - matK) / maximum(matK[:])
    @test norm(K - matK) / maximum(matK[:]) < 1.0e-4
end
end
using .m1testop6
m1testop6.test()

module m1testop7
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule: stiffness
using LinearAlgebra: norm, cross
using SparseArrays
using DelimitedFiles
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookean(mr, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    fens, fes = H8hexahedron([0 0 0; 1.2 0 -0.1; 1.0 0.9 0.1; -0.05 1.1 0; 0 0 1.03; 1.2 0 1.1; 1.0 0.95 0.81; -0.05 1.01 0.9], 1, 1, 1)
    Rmout = fill(0.0, 3, 3)
    rv = vec([-0.56 -0.1361 0.35])
    rotmat3!(Rmout, rv)
    mcsys = CSys(Rmout)
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), mcsys, m)

    geom = NodalField(fens.xyz)
    un1 = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    un = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    numberdofs!(un1)
    d(XYZ) = 0.1 * [XYZ[1]^2*XYZ[3], XYZ[2]^2*XYZ[3], XYZ[2]*XYZ[3]-XYZ[3]^3];
    for i = 1:count(fens)
        un1.values[i, :]  .+= d(fens.xyz[i, :]);
    end
    #println("un1.values = $(un1.values)")
    tn, dtn = 0.0, 1.0
    K = stiffness(femm, geom, un1, un, tn, dtn)
    ijs = readdlm("m1testop5.csv", ',', Float64)
    matK  = sparse(Int.(ijs[:, 1]), Int.(ijs[:, 2]), ijs[:, 3], size(K, 1), size(K, 2))
    # @show norm(K - matK) / maximum(matK[:])
    @test norm(K - matK) / maximum(matK[:]) < 1.0e-4
end
end
using .m1testop7
m1testop7.test()

module m1testop7a
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule: restoringforce
using LinearAlgebra: norm, cross
using SparseArrays
using DelimitedFiles
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookean(mr, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    fens, fes = H8hexahedron([0 0 0; 1.2 0 -0.1; 1.0 0.9 0.1; -0.05 1.1 0; 0 0 1.03; 1.2 0 1.1; 1.0 0.95 0.81; -0.05 1.01 0.9], 1, 1, 1)
    Rmout = fill(0.0, 3, 3)
    rv = vec([-0.56 -0.1361 0.35])
    rotmat3!(Rmout, rv)
    mcsys = CSys(Rmout)
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), mcsys, m)

    geom = NodalField(fens.xyz)
    un1 = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    un = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    numberdofs!(un1)
    d(XYZ) = 0.1 * [XYZ[1]^2*XYZ[3], XYZ[2]^2*XYZ[3], XYZ[2]*XYZ[3]-XYZ[3]^3];
    for i = 1:count(fens)
        un1.values[i, :]  .+= d(fens.xyz[i, :]);
    end
    # println("un1.values = $(un1.values)")
    tn, dtn = 0.0, 1.0
    F = restoringforce(femm, geom, un1, un, tn, dtn)
    matF = readdlm("m1testop7.csv", ',', Float64)
    # @show norm(F - matF) / maximum(matF)
    @test norm(F - matF) / maximum(matF) < 1.0e-4
end
end
using .m1testop7a
m1testop7a.test()

module m1testop8
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule: restoringforce
using LinearAlgebra: norm, cross
using SparseArrays
using DelimitedFiles
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookean(mr, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    fens, fes = H8hexahedron([0 0 0; 1.2 0 -0.1; 1.0 0.9 0.1; -0.05 1.1 0; 0 0 1.03; 1.2 0 1.1; 1.0 0.95 0.81; -0.05 1.01 0.9], 1, 1, 1)
    function  mcsys(csmatout::FFltMat, XYZ, tangents, fe_label)
        # @bp
        n1 = vec(XYZ); n1[1] = 0; n1  = n1 / norm(n1);
        n3 = vec([1, 0, 0]);
        n2 = cross(n3, n1);
        csmatout .= hcat(n1, n2, n3);
    end
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), CSys(3, 3, mcsys), m)

    geom = NodalField(fens.xyz)
    un1 = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    un = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    numberdofs!(un1)
    d(XYZ) = 0.1 * [XYZ[1]^2*XYZ[3], XYZ[2]^2*XYZ[3], XYZ[2]*XYZ[3]-XYZ[3]^3];
    for i = 1:count(fens)
        un1.values[i, :]  .+= d(fens.xyz[i, :]);
    end
    # println("un1.values = $(un1.values)")
    tn, dtn = 0.0, 1.0
    F = restoringforce(femm, geom, un1, un, tn, dtn)
    matF = readdlm("m1testop7.csv", ',', Float64)
    # @show norm(F - matF) / maximum(matF)
    @test norm(F - matF) / maximum(matF) < 1.0e-4
end
end
using .m1testop8
m1testop8.test()

module m1testop9
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule: restoringforce
using LinearAlgebra: norm, cross
using SparseArrays
using DelimitedFiles
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookean(mr, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    fens, fes = H8hexahedron([0 0 0; 1.2 0 -0.1; 1.0 0.9 0.1; -0.05 1.1 0; 0 0 1.03; 1.2 0 1.1; 1.0 0.95 0.81; -0.05 1.01 0.9], 1, 1, 1)
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)

    geom = NodalField(fens.xyz)
    un1 = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    un = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    numberdofs!(un1)
    d(XYZ) = 0.1 * [XYZ[1]^2*XYZ[3], XYZ[2]^2*XYZ[3], XYZ[2]*XYZ[3]-XYZ[3]^3];
    for i = 1:count(fens)
        un1.values[i, :]  .+= d(fens.xyz[i, :]);
    end
    # println("un1.values = $(un1.values)")
    tn, dtn = 0.0, 1.0
    F = restoringforce(femm, geom, un1, un, tn, dtn)
    matF = readdlm("m1testop7.csv", ',', Float64)
    # @show norm(F - matF) / maximum(matF)
    @test norm(F - matF) / maximum(matF) < 1.0e-4
end
end
using .m1testop9
m1testop9.test()

module m1testop10
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule: nzebcloads
using LinearAlgebra: norm
using SparseArrays
using DelimitedFiles
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookean(mr, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    fens, fes = H8block(L, W, H, 1, 1, 1)
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)

    geom = NodalField(fens.xyz)
    un1 = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    un = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    for i in [3, 7, 8]
        for j in 1:3
            setebc!(un1, [i], true, j, 0.0)
        end
    end
    numberdofs!(un1)
    du = deepcopy(un1)
    tn, dtn = 0.0, 1.0
    d(XYZ) = 100.1 * [XYZ[1]^2*XYZ[3], XYZ[2]^2*XYZ[3], XYZ[2]*XYZ[3]-XYZ[3]^3];
    for i = 1:count(fens)
        du.fixed_values[i, :]  .+= d(fens.xyz[i, :]);
    end
    # @show du
    F = nzebcloads(femm, geom, un1, un, du, tn, dtn)
    matF = readdlm("m1testop10.csv", ',', Float64)
    # @show norm(K - matK) / maximum(matK[:])
    @test norm(F - matF) / maximum(matF[:]) < 1.0e-4
end
end
using .m1testop10
m1testop10.test()

module m1testop11
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule: geostiffness
using LinearAlgebra: norm
using SparseArrays
using DelimitedFiles
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookean(mr, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    fens, fes = H8hexahedron([0 0 0; 1.2 0 -0.1; 1.0 0.9 0.1; -0.05 1.1 0; 0 0 1.03; 1.2 0 1.1; 1.0 0.95 0.81; -0.05 1.01 0.9], 1, 1, 1)
    Rmout = fill(0.0, 3, 3)
    rv = vec([-0.56 -0.1361 0.35])
    rotmat3!(Rmout, rv)
    mcsys = CSys(Rmout)
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), mcsys, m)

    geom = NodalField(fens.xyz)
    un1 = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    un = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    numberdofs!(un1)
    d(XYZ) = 0.1 * [XYZ[1]^2*XYZ[3], XYZ[2]^2*XYZ[3], XYZ[2]*XYZ[3]-XYZ[3]^3];
    for i = 1:count(fens)
        un1.values[i, :]  .+= d(fens.xyz[i, :]);
    end
    #println("un1.values = $(un1.values)")
    tn, dtn = 0.0, 1.0
    Kgeo = geostiffness(femm, geom, un1, un, tn, dtn)
    # @show Kgeo[5,5]
    ijs = readdlm("m1testop11.csv", ',', Float64)
    matK  = sparse(Int.(ijs[:, 1]), Int.(ijs[:, 2]), ijs[:, 3], size(Kgeo, 1), size(Kgeo, 2))
    # @show norm(K - matK) / maximum(matK[:])
    @test norm(Kgeo - matK) / maximum(matK[:]) < 1.0e-4
end
end
using .m1testop11
m1testop11.test()

module m1testop12
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule: geostiffness
using LinearAlgebra: norm
using SparseArrays
using DelimitedFiles
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookean(mr, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    fens, fes = H8hexahedron([0 0 0; 1.2 0 -0.1; 1.0 0.9 0.1; -0.05 1.1 0; 0 0 1.03; 1.2 0 1.1; 1.0 0.95 0.81; -0.05 1.01 0.9], 1, 1, 1)
    # Rmout = fill(0.0, 3, 3)
    # rv = vec([-0.56 -0.1361 0.35])
    # rotmat3!(Rmout, rv)
    # mcsys = CSys(Rmout)
    # femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), mcsys, m)
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)

    geom = NodalField(fens.xyz)
    un1 = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    un = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    numberdofs!(un1)
    d(XYZ) = 0.1 * [XYZ[1]^2*XYZ[3], XYZ[2]^2*XYZ[3], XYZ[2]*XYZ[3]-XYZ[3]^3];
    for i = 1:count(fens)
        un1.values[i, :]  .+= d(fens.xyz[i, :]);
    end
    #println("un1.values = $(un1.values)")
    tn, dtn = 0.0, 1.0
    Kgeo = geostiffness(femm, geom, un1, un, tn, dtn)
    # @show Kgeo[5,5]
    ijs = readdlm("m1testop11.csv", ',', Float64)
    matK  = sparse(Int.(ijs[:, 1]), Int.(ijs[:, 2]), ijs[:, 3], size(Kgeo, 1), size(Kgeo, 2))
    # @show norm(K - matK) / maximum(matK[:])
    @test norm(Kgeo - matK) / maximum(matK[:]) < 1.0e-4
end
end
using .m1testop12
m1testop12.test()
