# QDLDL.jl - A free LDL factorisation routine 
[![Build Status](https://travis-ci.com/oxfordcontrol/QDLDL.jl.svg?branch=master)](https://travis-ci.com/oxfordcontrol/QDLDL.jl)
[![codecov](https://codecov.io/gh/oxfordcontrol/QDLDL.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/oxfordcontrol/QDLDL.jl)

QDLDL is a factorisation routine for quasi-definite linear systems `Ax=b`. This is a pure Julia implementation of the C language QDLDL solver (https://github.com/oxfordcontrol/qdldl) with some additional functionality implemented to support refactorisations. 


## Getting Started

QDLDL requirers Julia v1.0 and can be added via the Julia package manager (type `]`): `pkg> add QDLDL`. Make the package available in you project with `using QDLDL`.

## Using QDLDL
QDLDL can be used to solve linear systems `Ax = b`.
Given a quasidefinite matrix `A` and right-hand side vector `b`, compute the factorisation `F` with:
```julia
F = qdldl(A)
```
Solve the linear system for `x` with:
```julia
x = solve(F, b)
```
This will allocate new memory for `x`. To solve in-place and overwrite `b` with `x` use:
```julia
solve!(F, b)
```

## Authors

* [Paul Goulart](http://users.ox.ac.uk/~engs1373/)


## License

This project is licensed under the Apache License - see the [LICENSE.md](LICENSE) file for details.
