# What the doctests of this package need in scope.
#
# Included both by `docs/make.jl` and by the `doctest` job of `.github/workflows/CI.yml`, so that a
# documentation build and a doctest run cannot disagree. One definition, two callers: the CI
# workflow is byte-identical in every repository and cannot carry per-package knowledge, and a
# second copy of this list is exactly the thing that goes stale.
#
# The Legendre and Lagrange examples integrate against a quadrature rule, and the accessors are
# shared with QuadratureRules, so QuadratureRules has to be in scope. `vandermonde_matrix` and
# `vandermonde_matrix_inverse` are not exported, so they are imported by name.

using Documenter: DocMeta

using CompactBasisFunctions
using QuadratureRules

DocMeta.setdocmeta!(
    CompactBasisFunctions,
    :DocTestSetup,
    :(using CompactBasisFunctions;
    using QuadratureRules;
    import CompactBasisFunctions: vandermonde_matrix, vandermonde_matrix_inverse);
    recursive = true
)
