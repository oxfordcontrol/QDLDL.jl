
using Test, LinearAlgebra, SparseArrays, Random
using QDLDL
rng = Random.MersenneTwister(0706)

@testset "linear solves" begin

  m = 20
  n = 30

  A = random_psd(n)
  B = sprandn(m,n,0.2)
  C = -random_psd(m)
  M = [A B' ; B C]

  b = randn(m+n)
  F = qdldl(M)

  #basic solve
  @test norm(F\b -  M\b, Inf) <= 1e-10

  #solve function
  @test norm(solve(F,b) -  M\b, Inf) <= 1e-10

  #solve in place
  x = copy(b)
  solve!(F,x)
  @test norm(x -  M\b, Inf) <= 1e-10

  #scalar system
  A = sprandn(1,1,1.)
  b = randn(1)
  @test norm(qdldl(A)\b -  A\b, Inf) <= 1e-10

  #32 bit floats, 32 bit ints
  b = randn(m+n)
  M32 = SparseMatrixCSC{Float32,Int32}(M)
  b32 = Vector{Float32}(b)
  @test norm(qdldl(M32)\b32 -  M\b, Inf) <= 1e-6

  #32 bit floats, 64 bit ints
  Mx = SparseMatrixCSC{Float32,Int64}(M)
  @test norm(qdldl(Mx)\b32 -  M\b, Inf) <= 1e-6

  #64 bit floats, 32 bit ints
  Mx = SparseMatrixCSC{Float64,Int32}(M)
  @test norm(qdldl(Mx)\b -  M\b, Inf) <= 1e-6


end

nothing
