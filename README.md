# CompactBasisFunctions

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGNI.github.io/CompactBasisFunctions.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGNI.github.io/CompactBasisFunctions.jl/dev)
[![PkgEval Status](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/C/CompactBasisFunctions.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/C/CompactBasisFunctions.html)
[![Build Status](https://github.com/JuliaGNI/CompactBasisFunctions.jl/workflows/CI/badge.svg)](https://github.com/JuliaGNI/CompactBasisFunctions.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaGNI/CompactBasisFunctions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGNI/CompactBasisFunctions.jl)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.4317806.svg)](https://doi.org/10.5281/zenodo.4317806)

This package provides a set of basis functions, mostly compactly supported, which are implemented as [Continuum Arrays](https://github.com/JuliaApproximation/ContinuumArrays.jl). Bases are accessed like arrays with continuous dimensions, e.g. as `b[0.1, 2]` to evaluate the second basis function in the point `0.1`. Operations such as derivatives, inner products and mass matrices are implemented as [Lazy Array](https://github.com/JuliaArrays/LazyArrays.jl) operations, providing a high-level linear algebra interface. Functions in a basis are represented by a lazy multiplication of the basis and a vector of coefficients, which materializes only upon evaluation of the product.

## References

If you use CompactBasisFunctions.jl in your work, please consider citing it by

```
@misc{Kraus:2020:CompactBasisFunctions,
  title={CompactBasisFunctions.jl: Compactly supported basis functions in Julia},
  author={Kraus, Michael},
  year={2020},
  howpublished={\url{https://github.com/JuliaGNI/CompactBasisFunctions.jl}},
  doi={10.5281/zenodo.4317806}
}
```
