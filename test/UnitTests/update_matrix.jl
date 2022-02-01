using Test, LinearAlgebra, SparseArrays, Random
using QDLDL
rng = Random.MersenneTwister(2401)

@testset "matrix update checks" begin

  # random KKT system
  m = 20
  n = 30
  Q = random_psd(n)
  G = sprandn(rng, m, n, 0.3)
  A1,b1 = build_kkt(Q, G)

  # create factorisation
  F = qdldl(A1)

  # compare qdldl solve to \
  @test norm(F\b1 - A1\b1) < 1e-10

  # update kkt with new data in the same sparsity pattern
  Q2 = copy(Q)
  Q2.nzval .= randn(rng, nnz(Q2))
  Q2 = Q2' + Q2
  G2 = copy(G)
  G2.nzval .= randn(rng, nnz(G2))
  A2,b2 = build_kkt(Q2, G2)

  # update factorisation in-place and test solution
  update_A!(F,A2)
  @test norm(F\b2 - A2\b2) < 1e-10

end
