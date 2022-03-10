using Test, LinearAlgebra, SparseArrays, Random
using QDLDL
rng = Random.MersenneTwister(2401)

@testset "refactoring checks" begin

  m = 20
  n = 30

  A  = random_psd(n) + Diagonal(rand(n))
  B  = sprandn(m,n,0.2)
  C1 = -Diagonal(rand(m))
  C2 = -Diagonal(rand(m))
  M1 = [A B' ; B C1]
  M2 = [A B' ; B C2]
  M3 = [A B' ; B -pi*I]

  b = randn(m+n)
  F = qdldl(M1,perm = randperm(m+n))

  @test norm(M1\b - F\b) <= 1e-10

  #update the diagonal of M and F and try again
  update_diagonal!(F,(n+1):(m+n),diag(C2))
  refactor!(F)
  @test norm(M2\b - F\b) <= 1e-10

  #update to a scalar and try again
  update_diagonal!(F,(n+1):(m+n),-pi)
  refactor!(F)
  @test norm(M3\b - F\b) <= 1e-10


end
