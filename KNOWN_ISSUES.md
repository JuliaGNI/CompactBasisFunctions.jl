# Known issues

What is known to be broken or incomplete and is not fixed yet. Delete an entry when its issue is
fixed; the fix goes in `CHANGELOG.md`.

## Upstream

### K1 · `fatou lint` reports false `parse-error` findings in three test files

- location: `test/bernstein.jl:156`
- evidence: `fatou lint --output json .` with fatou 0.20.0 reports 17 `parse-error` findings,
  severity `error`: 3 in `test/bernstein.jl` (real site `test/bernstein.jl:156`), 8 in
  `test/chebyshev.jl` (real sites `test/chebyshev.jl:257` and `:260`), 6 in `test/legendre.jl`
  (real sites `test/legendre.jl:146` and `:149`). Each real site is a
  `(@elapsed for _ in 1:100 … end) < 1.0` construct. fatou's parser reads the `(` before a macro
  call with a `for` body as an unclosed generator, and one site then cascades into several
  findings ("unclosed comprehension", "trailing tokens after statement", "unexpected block
  keyword", "expected `end`"). `Meta.parseall` on each of the three files gives no error nodes,
  and the test suite runs them.
- kind: upstream
- found: 2026-09-13

## Tests

### K2 · `quality/jet.jl` does not see an unstable return value in an entry point's own frame

- location: `test/quality/jet.jl:33`
- evidence: `report_opt` reports runtime dispatch, not a non-concrete return type. The mutants
  `order(b::PolynomialBasis) = Base.inferencebarrier(nbasis(b))` and
  `Base.inferencebarrier(_eval(b, x, j))` in `getindex(b::PolynomialBasis, x::Number, j::Integer)`
  both survive `quality/jet.jl`; the `order` mutant also survives `basis.jl`. An `@inferred` on the
  counting accessors in `test/basis.jl` would catch it.
- kind: missing test
- found: 2026-09-28

### K3 · `quality/jet.jl` covers only the `Float64` types of the `@allocated` test

- location: `test/quality/jet.jl:15`
- evidence: other tests of the same entry points use `BigFloat` bases, a `BigFloat` point, and
  `isapprox(ChebyshevT(Float32, 3), ChebyshevT(3))`. `report_opt` gives 0 reports at each of those
  types for all 11 entry points, so lines for them would pass today; they are not in the file.
- kind: missing test
- found: 2026-09-28
