# Changelog

All notable changes to CompactBasisFunctions.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Releases before v0.3.0 are
not covered here; see the git history for those.


## [0.3.1]

### Fixed

- **An evaluation now runs in the wider of the basis's element type and the argument's.** A basis
  of extended precision reported a precision it did not carry:

  ```julia
  b = Legendre(BigFloat, 6)
  b[0.3, 5] == b[BigFloat(0.3), 5]        # false before, true now
  ```

  `Legendre` widened only the trailing `√(2j+1)` factor and ran Bonnet's recurrence in the type of
  the argument, so a `BigFloat` basis evaluated at a `Float64` point returned a `BigFloat` correct
  to `1.3e-16` relative and no further. `Bernstein` and `Chebyshev` did not widen at all and
  returned a `Float64` outright — although their derivative evaluators already promoted, so a
  basis and its derivative disagreed about their own element type. This is the same class of bug
  as the `Lagrange` buffer allocation fixed in 0.3.0.

  All six evaluators now convert the argument first. The three that already promoted —
  `Legendre`'s derivative and both of `Chebyshev`'s — did so only after the shift `2x-1` onto
  ``[-1,1]``, which left the shift itself running at the argument's precision and then recorded
  its rounding in more digits; the conversion now precedes the shift everywhere. That omission
  is invisible for `x ≥ 0.25`, where `2x-1` is exact in `Float64`, so it takes a point such as
  `0.1` to see it. `Lagrange` needed no change; it promotes elementwise against its own nodes.

  The promotion only ever widens, so an argument more precise than the basis is untouched and
  every value at matching precision is unchanged. What changes is the return type of `Bernstein`
  and `Chebyshev` where `T` is wider than the argument, and the accuracy of `Legendre` there.

- **`Lagrange` accepted node sets whose differences are degenerate.** The check used `allunique`,
  which compares with `isequal`, and so asked a different question from the one that matters:

  ```julia
  allunique([0.0, -0.0])      # true — yet 0.0 - (-0.0) is 0.0
  Lagrange([0.0, -0.0])       # accepted, denominators Inf
  Lagrange([0.0, NaN, 1.0])   # accepted, every value NaN
  ```

  It now tests the differences that each denominator is built from, rejecting a vanishing one as
  a repeated node and a non-finite node as such. Their product has to be representable too, and
  distinct finite nodes can still overflow or underflow it — `Lagrange([0.0, 1e160, 2e160,
  3e160])`, which used to come out with every denominator `0`, now says that the nodes are sound
  and the product is not, and suggests rescaling them. Repeated nodes still throw as before,
  repeated `NaN` still throws, and `-0.0` is still a perfectly good node among distinct ones.

### Changed

- `docs/Project.toml` drops its `[sources]` block: the documentation workflow already develops the
  package from the checkout, and `[sources]` is understood only by Pkg 1.11 and later while this
  package supports Julia 1.10. The root `Project.toml` dropped its own in 0.3.0 for the same
  reason.


## [0.3.0]

**This release is breaking.** Three numerical results change: the Chebyshev basis now lives on
`[0,1]` instead of `[-1,+1]`, so every Chebyshev node, basis function and derivative value is
different; the Legendre derivative gains the normalisation factor it was missing; and the
`ChebyshevU` derivative returns finite values at the endpoints instead of `NaN`. The accessors
`basis`, `degree`, `nodes`, `nnodes` and `order` are now extended from `GeometricBase` rather than
defined here, the `FastTransforms` dependency is gone, and the lower bounds on Julia,
`ContinuumArrays` and `QuadratureRules` all move up.

### Added

- **A manual.** The documentation was a single `index.md` containing `@autodocs`, which produced
  almost nothing, because only the four basis structs carried docstrings — one line each — and the
  accessors carried none. There are now a conceptual page on polynomial approximation, a page per
  basis family, a usage page, a generated API reference and a bibliography, with every example a
  `jldoctest` so the build verifies it. Docstrings are written for everything exported, for the six
  derivative types and for the Vandermonde functions.

- **This changelog.**

### Fixed

- **The Legendre derivative was missing the `√(2j+1)` normalisation** that its basis functions
  carry, so `(d*Legendre(n))[x,j]` was not the derivative of `Legendre(n)[x,j]` — it was wrong by
  exactly that factor for every `j > 0`:

  ```julia
  b = Legendre(4);  d = Derivative(axes(b,1))
  (d*b)[0.37, 1]                             # 2.0        — before
  (b[0.37+h, 1] - b[0.37-h, 1]) / 2h         # 3.4641… = 2√3
  ```

  The normalisation was introduced in v0.2.0 and the derivative was never updated with it.
  Nothing in the test suite related the two — the basis assertions used `√(2j+1)` and the
  derivative assertions used the unnormalised values, so each agreed with itself. The suite now
  checks every basis against a central difference taken in 256-bit `BigFloat`; the other three
  were already consistent.

  This changes results for anyone differentiating a Legendre basis, which includes the CGVI and
  DGVI integrators of GeometricIntegrators and NonlinearIntegrators.

- **`(d*ChebyshevU(n))[x, j]` returned `NaN` at `x = 0` and `x = 1`** — the whole boundary of the
  basis's own domain — for every `j`. The closed form
  `U'ᵢ = ((i+1)Tᵢ₊₁ - x̃ Uᵢ) / (x̃² - 1)` is `0/0` at `x̃ = ±1`, although the derivative is a
  polynomial and finite there. It is replaced by the differentiated recurrence
  `U'ᵢ = 2Uᵢ₋₁ + 2x̃ U'ᵢ₋₁ - U'ᵢ₋₂`, which has no division and so no singularity. Interior values
  are unchanged bit for bit, apart from cases that returned `-0.0` and now return `0.0`. The
  existing tests checked interior points only, with a comment saying why.

- **Every polynomial evaluator was exponential in the degree.** `_bernstein`, `_chebyshev`,
  `_legendre` and `_legendre_derivative` each recursed into two subproblems per step and
  recomputed shared subtrees. One evaluation, measured:

  | n | Bernstein (mid index) | Legendre | ChebyshevT |
  |---|---|---|---|
  | 10 | 2.3 µs | 0.36 µs | 0.21 µs |
  | 20 | 1.5 ms | 51 µs | 29 µs |
  | 25 | **45 ms** | 571 µs | 323 µs |

  Iterating the same recurrences upwards gives 120 ns, 196 ns and 119 ns at `n = 25`, growing
  linearly thereafter. For Chebyshev and Legendre the per-step arithmetic and its order are
  unchanged, so results are **bitwise identical**, verified over 840 values. Bernstein's
  recurrence has two indices, so it is evaluated from the closed form `C(p,j) xʲ (1-x)^(p-j)`
  instead; 241 of 1170 sampled values differ, by at most `4.7e-16` relative, and are on balance
  closer to exact than the recursion they replace. Accuracy was never the problem — only the cost.

- **`Lagrange` accepted duplicate nodes** and silently produced `Inf` in its denominators and
  `NaN` on evaluation. It now throws an `ArgumentError` naming the nodes.

- **`grid` on a modal basis** raised `MethodError: no method matching grid_axis(…)` from
  ContinuumArrays' internals, naming nothing a caller could act on. `nodes`, `nnodes` and `grid`
  now report that `Bernstein` and `Legendre` are modal bases and therefore have no nodes. The
  generic stubs in `basis.jl` name the offending type instead of saying `"Not implemented!"`.

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
  the one Julia 1.13 provides, causing a method overwriting error. The `0.2` lower bound below
  picks that up.

- Regression tests for the numerical changes above, several of which no assertion covered: dropping
  a chain-rule factor, reverting the `zeros(T, …)` buffers, or omitting the Legendre normalisation
  all left the suite green. Chebyshev now asserts that the nodes lie inside `axes(C,1)`, are
  ascending and `n` in number, pins the two node vectors that changed, and checks basis and
  derivative values against their closed forms at the endpoints as well as inside. Lagrange asserts
  the cardinal property and the partition of unity, the latter to full `BigFloat` precision.
  Bernstein is checked against the recurrence it no longer uses, evaluated in `BigFloat`. Every
  basis is checked against a central difference, and each family has a cost guard at high degree.
  `runtests.jl` seeds the RNG, so the `rand()` draws several files rely on are reproducible and the
  Vandermonde tests cannot pass or fail by luck of the draw.

### Changed

- **The accessors are extended from `GeometricBase`.** `basis`, `degree`, `nodes`, `nnodes` and
  `order` were defined here, and independently in QuadratureRules, RungeKutta and the integrator
  packages, so they were distinct functions that happened to share a name. Loading two of them
  together made the name resolve to nothing:

  ```julia
  using QuadratureRules, CompactBasisFunctions
  order(Legendre(3))    # UndefVarError: `order` not defined
  ```

  which is the normal combination, since this package depends on QuadratureRules. All five are now
  imported from `GeometricBase`, which declares them method-free for exactly this purpose — the
  pattern `grid` already followed with `ContinuumArrays.grid`. Requires `GeometricBase` 0.14.8 and
  `QuadratureRules` 0.2.

  `nodes` and `nnodes` are exported now as well; `nbasis` remains this package's own, as nothing
  else in the ecosystem defines it.

  Note that `basis` still collides with `ContinuumArrays`, which exports its own and means
  something different by it: for a basis object ContinuumArrays returns the basis itself, where
  this package returns the vector of basis functions. Qualify or import explicitly when loading
  both.

- **`Lagrange` no longer stores `vdminv`.** Every constructor filled it with
  `vandermonde_matrix_inverse(x)` and nothing ever read it — not here, not in GeometricIntegrators,
  RungeKutta, NonlinearIntegrators or MultiSymplectic — so each basis inverted a Vandermonde matrix
  for nothing. The function itself stays.

- Out-of-range basis indices throw `BoundsError` rather than `AssertionError`, so the derivative and
  the basis agree on what an out-of-range index does.

- `Bernstein` and `Legendre` gain the `isapprox` that `Chebyshev` and `Lagrange` already had.

- `vandermonde_matrix_inverse` allocates in `T` rather than allocating `Float64` and converting —
  the same class of bug as the `Lagrange` element-type fix above — so integer node vectors now work
  instead of throwing `InexactError`.

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
  to `0.18, 0.19, 0.20`. `GeometricBase` is a new dependency at `0.14.8`, and `QuadratureRules`
  requires `0.2`, both for the shared accessors above. The missing `LinearAlgebra = "1"` and
  `Test = "1"` bounds are added, `Random` joins the test target, and `docs/Project.toml` gains
  `DocumenterCitations = "1"` and `QuadratureRules = "0.2"` alongside the `Documenter = "1"` it
  already carried.

- All node generation goes through the `*_nodes` functions of `QuadratureRules`, which as of its
  0.2 release return the nodes on `[0,1]` — the interval of every basis here — unless the
  `interval` keyword asks for `[-1,+1]` instead:

  ```
  chebyshevpoints(T, n, Val(kind))    ->  chebyshev_nodes(T, n, Val(kind))
  GaussLegendreQuadrature(n).nodes    ->  gauss_legendre_nodes(n)
  LobattoLegendreQuadrature(n).nodes  ->  lobatto_legendre_nodes(n)
  ```

  For `LagrangeGauß` and `LagrangeLobatto` this is a cleanup with no change in behaviour: the same
  values of the same type, without constructing the quadrature weights only to discard them.

- Docstrings added for `Chebyshev`, for both `_chebyshev` helpers and for both Chebyshev derivative
  evaluators, each stating the interval it applies to.

### Removed

- **`FastTransforms` is no longer a dependency.** It was pulled in for the single function
  `chebyshevpoints`. Dropping it also removes `FastTransforms_jll`, `FFTW`, `ToeplitzMatrices`,
  `GenericFFT` and `DSP` from the dependency tree.


[0.3.1]: https://github.com/JuliaGNI/CompactBasisFunctions.jl/compare/v0.3.0...main
[0.3.0]: https://github.com/JuliaGNI/CompactBasisFunctions.jl/compare/v0.2.15...v0.3.0
