module ScriForwardDiffExt

import Scri
import ForwardDiff

# Peel dual numbers to their value type (recursively, for nested duals), so that
# `transform!` builds its parameter-independent grid and analysis factorizations at the
# primal precision.  See the docstring of `Scri.primal_float`.
Scri.primal_float(::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N} = Scri.primal_float(V)

end  # module ScriForwardDiffExt
