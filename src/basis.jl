
_not_implemented(f, b) = error("$(f) is not implemented for $(typeof(b)).")

"""
Error for the accessors that a modal basis cannot answer: its coefficients are those of an
expansion, not values at points, so it has no nodes and no grid. Saying so beats the
`MethodError` about ContinuumArrays' `grid_axis` that `grid` would otherwise produce.
"""
_no_nodes(b, f) =
    error("$(nameof(typeof(b))) is a modal basis and has no nodes, so $(f) is not defined for it.")

# `nodes` and `nnodes` are deliberately left unimplemented for the modal bases, which have
# no nodes at all; those carry their own methods with a message saying so. A generic
# `grid(::Basis)` is not defined here on purpose: ContinuumArrays dispatches `grid` on
# `Any`, so a method on its `Basis` supertype would also capture its own spline bases.

basis(b::Basis) = _not_implemented("basis", b)
nodes(b::Basis) = _not_implemented("nodes", b)

nbasis(b::Basis) = _not_implemented("nbasis", b)
nnodes(b::Basis) = _not_implemented("nnodes", b)

order(b::Basis) = _not_implemented("order", b)
degree(b::Basis) = _not_implemented("degree", b)
