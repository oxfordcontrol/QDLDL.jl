
using Test, LinearAlgebra, SparseArrays, Random
using QDLDL
rng = Random.MersenneTwister(1312)

@testset "non-quasidefinite" begin

  m = 20

  A = random_psd(m)
  A[:,10] .= 0
  A[10,:] .= 0

  @test_throws ErrorException qdldl(A)

end

nothing
