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
