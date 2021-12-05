using Random, Test
rng = Random.MersenneTwister(12345)

function random_psd(n)

    A = sprandn(rng,n,n,0.2)
    A = A+A';
    A = A + Diagonal(((sum(abs.(A),dims=1)[:]))) #make diagonally dominant

end

function build_kkt(Q,G)
  n = size(Q,1)
  m = size(G,1)

  A = [Q G';G zeros(m,m)] + 1e-2*Diagonal([ones(n);-ones(m)])
  A = sparse(A)

  b = randn(rng, n + m)

  return A,b
end

@testset "All Unit Tests" begin

  include("./UnitTests/basic.jl")
  include("./UnitTests/refactoring.jl")
  include("./UnitTests/non-quasidef.jl")
  include("./UnitTests/inertia.jl")
  include("./UnitTests/update_matrix.jl")

end
nothing
