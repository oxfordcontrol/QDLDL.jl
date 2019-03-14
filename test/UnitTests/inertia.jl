
using Test, LinearAlgebra, SparseArrays, Random
using QDLDL
rng = Random.MersenneTwister(2701)

@testset "inertia checks" begin

  m = 20
  n = 30

  A = random_psd(n) + Diagonal(rand(n))
  B = sprandn(m,n,0.2)
  C = -Diagonal(rand(m))
  M = [A B' ; B C]

  b = randn(m+n)
  F = qdldl(M)

  @test positive_inertia(F) == n

end

nothing
