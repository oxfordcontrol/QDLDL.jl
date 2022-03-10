
using Test, LinearAlgebra, SparseArrays, Random
using QDLDL

@testset "regularization" begin

  a = 1e-20
  A = sparse([a 0;0 a])

  signs = [1,1]
  f = qdldl(A,Dsigns = signs)
  @test f.workspace.D == [1e-7,1e-7]
  @test regularized_entries(f) == 2

  signs = [1,-1]
  f = qdldl(A,Dsigns = signs)
  @test f.workspace.D == [1e-7,-1e-7]
  @test regularized_entries(f) == 2

  δ = 1e-8
  f = qdldl(A,Dsigns = signs,regularize_delta = δ)
  @test f.workspace.D == [δ,-δ]
  @test regularized_entries(f) == 2

  ϵ = 1e-30
  f = qdldl(A,Dsigns = signs,regularize_eps = ϵ, regularize_delta = δ)
  @test f.workspace.D == [a,-δ]
  @test regularized_entries(f) == 1

  #just setting A[1,1] = 0 should still allow
  #regularization, since A[1,1] will be structurally
  #entry on the diagonal still
  A[1,1] = 0.
  f = qdldl(A,Dsigns = signs)
  @test f.workspace.D == [1e-7,-1e-7]

  #but the same matrix should fail on non regularized
  #factorization since its D[1,1] entry will be exactly zero
  @test_throws ErrorException qdldl(A)

  #regularization should also fail when there is
  #no structural entry on the diagonal, regardless
  #of whether regularization is enabled
  A = sparse([0. 0;0 a])
  @test_throws ErrorException qdldl(A)
  @test_throws ErrorException qdldl(A,Dsigns = signs)


  #make sure that signs are being properly permuted
  A = sparse(I(3).*a)
  signs = [1,1,-1]
  perm  = [2,3,1]
  f = qdldl(A;perm=perm,Dsigns = signs)
  @test f.workspace.D == [1e-7,-1e-7,1e-7]
  @test regularized_entries(f) == 3

  #try a biggish problem
  rng = Random.MersenneTwister(2401)
  m = 20
  n = 30
  A = random_psd(n) + Diagonal(rand(n))
  B = sprandn(m,n,0.2)
  C = -Diagonal(rand(m))
  M = [A B' ; B C]
  b = randn(m+n)
  s = [ones(Int64,n);-ones(Int64,m)]
  F = qdldl(M,perm=randperm(m+n),Dsigns = s)
  #basic solve
  @test norm(F\b -  M\b, Inf) <= 1e-10



end

nothing
