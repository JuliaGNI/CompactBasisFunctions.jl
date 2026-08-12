```@meta
CurrentModule = CompactBasisFunctions
DocTestSetup = quote
    using CompactBasisFunctions
    using QuadratureRules
    import CompactBasisFunctions: vandermonde_matrix, vandermonde_matrix_inverse
end
```

# Polynomial Approximation

The problem this package addresses is to represent a function ``f`` on an interval by a
finite combination of simpler functions,

```math
f(x) \approx \sum_{j} c_j \, \phi_j(x) ,
```

where the ``\phi_j`` are the *basis functions* and the ``c_j`` the *coefficients*. Once such
a representation is fixed, the operations one actually wants — differentiating, integrating,
taking inner products, solving a variational problem — become linear algebra on the
coefficients, and the only thing that varies from one method to the next is which basis was
chosen and why.

Everything in this package is polynomial: each ``\phi_j`` is a polynomial of degree at most
``p``, so the four bases all span the same space, the polynomials of degree ``\le p``. They
differ only in *which* basis of that space they provide. That sounds like a distinction
without a difference, and mathematically it is; numerically and practically it is not. A
basis determines what the coefficients mean, how expensive and how stable it is to evaluate
them, and how easily boundary conditions, positivity, or orthogonality can be imposed.

For the general theory see [Trefethen:2019:ApproximationTheory](@cite) and, for the
spectral-methods perspective that motivates most of the choices here,
[Boyd:2001:ChebyshevFourier](@cite) and [Canuto:2006:SpectralMethods](@cite).

## Interpolation and projection

There are two natural ways to pick the coefficients.

**Interpolation** demands agreement at a finite set of *nodes* ``x_i``:

```math
\sum_j c_j \, \phi_j(x_i) = f(x_i) \qquad \text{for all } i .
```

This is a linear system with the matrix ``\Phi_{ij} = \phi_j(x_i)``. With as many nodes as
basis functions it has a unique solution whenever the nodes are distinct, and the quality of
the result depends on the nodes far more than on the basis.

**Projection** demands agreement in an integral sense, minimising ``\| f - \sum_j c_j
\phi_j \|`` in the ``L^2`` norm, which gives the normal equations

```math
\sum_j M_{ij} \, c_j = \int f(x) \, \phi_i(x) \, dx , \qquad
M_{ij} = \int \phi_i(x) \, \phi_j(x) \, dx ,
```

with ``M`` the *mass matrix* of the basis. Projection needs no nodes but it needs integrals,
which in practice are evaluated by quadrature — which is why this package is built on
[QuadratureRules.jl](https://github.com/JuliaGNI/QuadratureRules.jl) and why the two share
their [`nodes`](@ref), `weights` and [`order`](@ref) accessors.

Interpolation is cheap and local; projection is optimal in the norm it minimises. Variational
integrators, the application this package was written for, use both: a Lagrange or Chebyshev
basis to interpolate a trajectory through nodes, and a quadrature rule to approximate the
action integral [Marsden:2001:DiscreteMechanics](@cite).

## Nodal and modal bases

The distinction that organises this package is whether a basis comes with nodes.

A **nodal** basis is built from a set of points, and is *cardinal* with respect to them:
``\phi_j(x_i) = \delta_{ij}``. The interpolation matrix ``\Phi`` is then the identity, so the
coefficients simply *are* the function values, ``c_j = f(x_j)``, and interpolation costs
nothing. [`Lagrange`](@ref) is the nodal basis in its purest form, and [`Chebyshev`](@ref)
carries the Chebyshev points as its nodes.

A **modal** basis has no nodes. Its coefficients are the coefficients of an expansion —
amplitudes of modes — and recovering them from function values requires solving a system or
computing integrals. [`Legendre`](@ref) and [`Bernstein`](@ref) are modal.

This is why [`nodes`](@ref), [`nnodes`](@ref) and `grid` are defined for the first group and
throw an informative error for the second:

```jldoctest
julia> nnodes(LagrangeLobatto(3))
3

julia> nodes(Legendre(3))
ERROR: Legendre is a modal basis and has no nodes, so nodes is not defined for it.
[...]
```

[`nbasis`](@ref), [`order`](@ref) and [`degree`](@ref) are defined for all four.

## The reference interval

Every basis in this package is defined on the interval ``[0,1]``, and `axes(b, 1)` is
`Inclusion(0..1)` for all of them. `QuadratureRules` follows the same convention, so nodes,
weights and basis functions compose without rescaling.

The classical theory, on the other hand, lives on the symmetric interval ``[-1,+1]``, which
is where the Chebyshev and Legendre polynomials are defined and where their orthogonality
relations take their usual form. The two views are related by the affine map

```math
\tilde{x} = 2x - 1 , \qquad x \in [0,1] \Longleftrightarrow \tilde{x} \in [-1,+1] ,
```

and this package evaluates the classical polynomial at ``\tilde{x}``. Differentiating brings
in the chain-rule factor

```math
\frac{d}{dx} \, \phi(2x-1) = 2 \, \phi'(2x-1) ,
```

which is why every derivative of a shifted basis carries a factor `2`. [`Bernstein`](@ref)
and [`Lagrange`](@ref) need no shift: the former is defined on ``[0,1]`` to begin with, and
the latter is defined by whatever nodes it is given.

An interval ``[a,b]`` other than the reference one is handled by the caller, by composing
with ``x = a + (b-a)\xi``; a derivative then picks up ``1/(b-a)``.

## Choice of nodes

For a nodal basis, the nodes matter more than anything else. Interpolating at equidistant
nodes is unstable: as the degree grows, the interpolant develops oscillations near the ends of
the interval that grow without bound even for perfectly smooth functions. This is the *Runge
phenomenon* [Runge:1901:Interpolation](@cite), and the underlying reason is that the
Lebesgue constant of equidistant nodes grows exponentially in the degree.

Nodes clustered towards the endpoints avoid this. Two families are used here:

- the **Chebyshev points**, the roots or extrema of a Chebyshev polynomial, whose Lebesgue
  constant grows only logarithmically — these are the nodes of [`Chebyshev`](@ref);
- the **Gauß-** and **Lobatto-Legendre** nodes, the points of the corresponding quadrature
  rules, available as [`LagrangeGauß`](@ref) and [`LagrangeLobatto`](@ref).

The choice between the last two is usually about the boundary. Gauß nodes lie strictly
inside the interval and integrate to higher order; Lobatto nodes include both endpoints, so
the coefficients of an expansion include the boundary values themselves — what a method needs
when it has to impose or read off conditions there.

## The Vandermonde matrix

For the monomial basis ``\phi_j(x) = x^{j-1}`` the interpolation matrix is the *Vandermonde
matrix* ``V_{ij} = x_i^{j-1}``, available as [`vandermonde_matrix`](@ref). It maps monomial
coefficients to values at the nodes, so its inverse maps values to monomial coefficients, and
composing the two directions converts between the nodal and monomial representations of the
same polynomial.

```jldoctest
julia> V = vandermonde_matrix([0.0, 0.5, 1.0]);

julia> V * [1.0, 0.0, 1.0]                # values of 1 + x² at the nodes
3-element Vector{Float64}:
 1.0
 1.25
 2.0

julia> vandermonde_matrix_inverse([0.0, 0.5, 1.0]) * [1.0, 1.25, 2.0]
3-element Vector{Float64}:
 1.0
 0.0
 1.0
```

[`vandermonde_matrix_inverse`](@ref) computes the inverse in closed form, as a product of a
lower and an upper triangular factor whose entries are known explicitly from the nodes, rather
than by a general factorisation.

That said, the Vandermonde matrix is the standard example of an ill-conditioned matrix: its
condition number grows exponentially in the number of nodes, whatever the nodes are. This is
not a defect of the algorithm but of the monomial basis, and it is the practical argument for
working in one of the four bases here instead of in coefficients of ``x^k``.

## Orthogonality

A basis is orthogonal when its mass matrix is diagonal, and orthonormal when the mass matrix
is the identity. Then projection needs no linear solve at all: the coefficients are just the
inner products against the basis functions.

[`Legendre`](@ref) is orthonormal on ``[0,1]`` by construction. Since
``\int_0^1 P_i(2x-1) P_j(2x-1) \, dx = \delta_{ij}/(2j+1)``, scaling the ``j``-th basis
function by ``\sqrt{2j+1}`` gives exactly ``\int_0^1 L_i L_j \, dx = \delta_{ij}``, which is
what that factor in the definition is for.

```jldoctest
julia> l = Legendre(4);  quad = GaussLegendreQuadrature(8);

julia> M = [sum(weights(quad)[k] * l[nodes(quad)[k], i] * l[nodes(quad)[k], j]
                for k in eachindex(nodes(quad))) for i in 0:3, j in 0:3];

julia> maximum(abs, M - [i == j for i in 0:3, j in 0:3]) < 1e-14
true
```

The Chebyshev polynomials are orthogonal too, but with respect to a weighted inner product,
with weight ``(1-\tilde{x}^2)^{-1/2}`` for the first kind and ``(1-\tilde{x}^2)^{+1/2}`` for
the second [Mason:2003:ChebyshevPolynomials](@cite); they are not orthogonal in the
unweighted ``L^2`` inner product that the mass matrix above uses.

## Partition of unity

A basis forms a *partition of unity* when ``\sum_j \phi_j(x) = 1`` for every ``x``. Then a
constant function is represented by constant coefficients, an expansion reproduces constants
exactly, and — if the basis functions are also non-negative — the expansion is a convex
combination of its coefficients and therefore bounded by them.

[`Bernstein`](@ref) has both properties on ``[0,1]``, which is what makes it the basis of
Bézier curves and gives the convex-hull property that computer-aided design relies on
[Farouki:2012:BernsteinPolynomialBasis](@cite). [`Lagrange`](@ref) forms a partition of
unity as well, though its cardinal functions do take negative values.

A useful consequence: differentiating the identity gives ``\sum_j \phi_j'(x) = 0``, so the
rows of any derivative matrix sum to zero. Both identities make good tests, precisely because
they hold exactly rather than approximately, and the size of the deviation from them measures
the accumulated rounding error of the evaluation. This is how the arbitrary-precision
behaviour of [`Lagrange`](@ref) is checked in the test suite.

## Evaluation by recurrence

Classical orthogonal polynomials are not evaluated from their explicit coefficients, which
would be both slow and unstable, but from three-term recurrences
[Gautschi:1967:ComputationalRecurrence](@cite):

```math
T_j = 2\tilde{x} T_{j-1} - T_{j-2} , \qquad
U_j = 2\tilde{x} U_{j-1} - U_{j-2} , \qquad
j P_j = (2j-1) \tilde{x} P_{j-1} - (j-1) P_{j-2} .
```

These are stable in the upward direction, and this package iterates them upwards, carrying
the two previous values. Derivatives are obtained by differentiating the same recurrences term
by term and carrying the derivative alongside the value:

```math
P_j' = \frac{(2j-1) P_{j-1} + (2j-1) \tilde{x} P_{j-1}' - (j-1) P_{j-2}'}{j} , \qquad
U_j' = 2 U_{j-1} + 2\tilde{x} U_{j-1}' - U_{j-2}' .
```

Doing this by *recursion* rather than iteration is a trap. Each step then descends into two
subproblems and recomputes shared subtrees, so the cost of one value grows like ``\varphi^j``
rather than ``j``. Earlier versions of this package did exactly that, at 45 ms for a single
degree-24 Bernstein value; the iterative form takes about 120 ns and grows linearly.

The Bernstein recurrence

```math
B_{j,p} = (1-x) \, B_{j,p-1} + x \, B_{j-1,p-1}
```

has two indices, so it cannot be collapsed into a single upward sweep for one ``j``; that
basis is evaluated from the closed form
``B_{j,p}(x) = \binom{p}{j} x^j (1-x)^{p-j}`` instead, with the binomial coefficient
accumulated multiplicatively over integer partial products; see
[Evaluation in closed form](@ref) for where that stops being exact.

## Bases as quasi-matrices

The interface follows [ContinuumArrays.jl](https://github.com/JuliaApproximation/ContinuumArrays.jl)
[Olver:2020:ContinuumArrays](@cite): a basis is a matrix with one continuous and one discrete
axis, a *quasi-matrix*. The continuous axis is the domain, the discrete one indexes the basis
functions, so that `b[x, j]` evaluates ``\phi_j(x)`` and reads exactly like indexing.

The payoff is that the linear algebra of the approximation is expressible directly:
`Derivative(axes(b,1)) * b` is the derivative operator applied to the basis, and `b * c` is
the function with coefficients `c` — both as lazy products, materialising only when indexed.
See [Usage](@ref) for what that looks like in practice.
