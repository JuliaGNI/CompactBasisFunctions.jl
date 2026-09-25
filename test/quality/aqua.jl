using Aqua
using CompactBasisFunctions
using Test

# Package-level quality assurance: type piracy, method ambiguities, stale and duplicated
# dependencies, undefined exports, unbound type parameters, `Project.toml` validity. These are
# the faults the rest of the suite is structurally unable to see — it exercises behaviour, and
# every one of these is a property of the package as a whole.
#
# Piracy is the one that bites here: the bases subtype ContinuumArrays' `Basis` and extend
# GeometricBase's accessors, so a method on `Basis` itself would own neither the function nor
# the argument type. Every method belongs to this package's own supertypes instead.
Aqua.test_all(CompactBasisFunctions)
