using Test, LinearAlgebra, SparseArrays, Random
using QDLDL
rng = Random.MersenneTwister(1312)

@testset "update values" begin

    # create random KKT system
    nz = 100
    nc = 70
    H = sprand(nz,nz,0.05);
    H = H'*H + I
    b = randn(nz + nc)
    A1 = sprand(nc,nz,0.8)
    K1 = [H A1';A1 -1e-3*I]

    # get triu of K1 and create QDLDL struct
    triuK1 = triu(K1);
    F = qdldl(triuK1)

    # compare with backslash
    @test norm(K1\b - F\b,Inf) < 1e-12

    # create a new KKT system with the same sparsity pattern
    A2 = copy(A1)
    A2.nzval .= randn(length(A2.nzval))
    K2 = [H A2';A2 -1e-7*I]
    triuK2 = triu(K2)

    # update factorization of F in place (non allocating)
    update_values!(F,1:length(triuK2.nzval),triuK2.nzval)
    refactor!(F)

    # compare with backslash using the refactored F
    @test norm(K2\b - F\b,Inf)<1e-12

end
