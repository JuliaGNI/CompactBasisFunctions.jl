# Changelog

All notable changes to CompactBasisFunctions.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Releases before v0.3.0 are
not covered here; see the git history for those.


## [0.3.0]

**This release is breaking.** The Chebyshev basis now lives on `[0,1]` instead of `[-1,+1]`, so
every Chebyshev node, basis function and derivative value changes. The `FastTransforms` dependency
is gone, and the lower bounds on Julia, `ContinuumArrays` and `QuadratureRules` all move up.

### Changed

- **The Chebyshev basis is defined on `[0,1]`.** `Base.axes(C::Chebyshev)` already returned
  `(Inclusion(0..1), eachindex(C))`, while the nodes and the basis functions lived on `[-1,+1]`,
  so `grid(C)` did not lie inside `axes(C,1)`. The declaration was the correct one and the values
  were wrong. The basis functions are now the Chebyshev polynomials evaluated at `2x-1`, and both
  derivatives carry the resulting chain-rule factor, i.e. `d/dx Tᵢ(2x-1) = 2i·U_{i-1}(2x-1)`. This
  follows the pattern `Legendre` already used, and makes all four bases agree on their domain.

  Any code that evaluates a `Chebyshev` basis, its derivative, or its nodes therefore gets
  different numbers. Where the old values are wanted, pull the argument back onto `[0,1]`:

  ```julia
  C = ChebyshevT(4)
  C[(x + 1) / 2, j]           # the value the old basis returned at x ∈ [-1,+1]
  2 .* nodes(C) .- 1          # the old node locations, up to their ordering
  ```

- **`nodes(C)` and `grid(C)` are in ascending order.** `FastTransforms.chebyshevpoints` returned
  the points in descending order, so this reverses them independently of the rescaling above. Code
  that indexes into the node vector, rather than iterating over it, is affected even after
  accounting for the change of domain.

- **`nodes(C)` and `grid(C)` return a `Vector{T}`** instead of the lazy
  `FastTransforms.ChebyshevGrid{kind,T}`, and the `XT` type parameter of `Chebyshev` changes with
  it. The node values themselves are slightly more accurate: `QuadratureRules.chebyshev_nodes`
  evaluates the closed form in `BigFloat` and rounds once to `T`, where `chebyshevpoints` computed
  in `T` throughout.

- **A `Chebyshev{kind,T}` whose nodes are not representable in `T` throws `InexactError` at
  construction** rather than on element access, because the conversion is no longer deferred by a
  lazy grid. In practice `T` must be a floating-point type; the exceptions are the cases where the
  nodes happen to be exact, such as `ChebyshevU(Integer, 2)`, whose nodes are `0` and `1`. This
  requirement is now documented on `Chebyshev`.

- `ChebyshevU(1)` still errors — Chebyshev points of the second kind need at least two points —
  but with an `ErrorException` from `QuadratureRules` in place of the old `ArgumentError` from
  `FastTransforms`.

- **Compat bounds.** Julia moves from `1.6` to `1.10`, and the CI matrix now covers `1.10`, `1.12`
  and `^1.13.0-0` plus nightly. The accreted `ContinuumArrays` list `0.8, 0.9, …, 0.20` is trimmed
  to `0.18, 0.19, 0.20`. `QuadratureRules` requires `0.1.10`. The missing `LinearAlgebra = "1"` and
  `Test = "1"` bounds are added, as is `Documenter = "1"` in `docs/Project.toml`.

- All node generation goes through `QuadratureRules`, whose 0.1.10 release added the `*_points`
  (on `[-1,+1]`) and `*_nodes` (on `[0,1]`) accessors:

  ```
  chebyshevpoints(T, n, Val(kind))    ->  chebyshev_nodes(T, n, Val(kind))
  GaussLegendreQuadrature(n).nodes    ->  gauss_legendre_nodes(n)
  LobattoLegendreQuadrature(n).nodes  ->  lobatto_legendre_nodes(n)
  ```

  For `LagrangeGauß` and `LagrangeLobatto` this is a cleanup with no change in behaviour: the same
  values of the same type, without constructing the quadrature weights only to discard them.

- Docstrings added for `Chebyshev`, for both `_chebyshev` helpers and for both Chebyshev derivative
  evaluators, each stating the interval it applies to.

### Fixed

- **`Lagrange` no longer loses precision above `Float64`.** The constructor allocated its buffers
  with `zeros(n)` and `zeros(n,n)`, i.e. in `Float64` regardless of `T`, so
  `diffs[i,j] = x[i] - x[j]` rounded every difference on assignment and `denom[i] = 1/p` rounded
  again; the `denom::XT` and `diffs::Matrix{T}` fields then widened those `Float64` values back up.
  A `Lagrange{BigFloat}` reported `BigFloat` while carrying only `Float64` precision, and since
  `L.diffs` also feeds the derivative evaluation, derivatives were affected too. The buffers are
  now allocated in `T`.

  On 256-bit Gauss-Legendre nodes the maximum deviation of the basis from the exact cardinal
  function drops from `2.7e-17` to below `1e-70`. `Float64` and `Float32` bases are unchanged in
  value; `eltype(L.denom)` and `eltype(L.diffs)` now follow `T`.

- **Precompilation on Julia 1.13.** `QuadratureRules` 0.1.8 dropped `GenericLinearAlgebra`, whose
  unconditional definition of `LinearAlgebra.eigencopy_oftype(::UpperHessenberg, S)` collides with
  the one Julia 1.13 provides, causing a method overwriting error. The `0.1.10` lower bound picks
  that up.

- Regression tests for both numerical changes above, neither of which was covered by a single
  assertion before: dropping a chain-rule factor or reverting the `zeros(T, …)` buffers left the
  suite green. Chebyshev now asserts that the nodes lie inside `axes(C,1)`, are ascending and `n`
  in number, pins the two node vectors that changed, and checks basis and derivative values against
  their closed forms, with the `i=2` branches pinning the chain-rule factor. Lagrange asserts the
  partition of unity and its derivative to full `BigFloat` precision.

### Removed

- **`FastTransforms` is no longer a dependency.** It was pulled in for the single function
  `chebyshevpoints`. Dropping it also removes `FastTransforms_jll`, `FFTW`, `ToeplitzMatrices`,
  `GenericFFT` and `DSP` from the dependency tree.


[0.3.0]: https://github.com/JuliaGNI/CompactBasisFunctions.jl/compare/v0.2.15...main
