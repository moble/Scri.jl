module Scri

import Quaternionic
import Quaternionic:
    QuatVec,
    Rotor,
    Lorentz,
    Boost,
    absvec,
    abs2vec,
    𝐤,
    components,
    basetype,
    from_spherical_coordinates,
    ℂreal
import SphericalFunctions:
    SphericalFunctions,
    sYlm_matrix,
    sYlmCalculator,
    set_R!,
    ð,
    golden_ratio_spiral_rotors,
    D
import ForwardDiff
import LinearAlgebra: LinearAlgebra, mul!, ldiv!, lu, qr
import Logging
import Hwloc
import OhMyThreads
import Base.Threads: nthreads
import TestItems: @testitem, @testmodule

# These are just for precompilation
using PrecompileTools: @setup_workload, @compile_workload
using Random: Xoshiro

# The L2 cache size of the machine running the code, in bytes, which sets how many time
# samples `transform!` handles in each block.  It is read when the package loads, in
# `__init__`, not when it is precompiled, so that a depot shared between machines does not
# carry the build machine's value.  The precompile workload runs before `__init__`, so
# `transform!` also reads it if it is still unset.
const cachesize_L2 = Ref(0)

"""
    read_cache_size!()

Read the L2 cache size of this machine from Hwloc into `cachesize_L2[]`.
"""
function read_cache_size!()
    Logging.with_logger(Logging.NullLogger()) do
        return cachesize_L2[] = Hwloc.cachesize(:L2)
    end
    return nothing
end

__init__() = read_cache_size!()

include("cubic_spline.jl")
include("signs.jl")
include("conventions.jl")
include("data_components.jl")
include("utilities.jl")
include("quadrature.jl")
include("objective_plan.jl")
include("transform_plan.jl")
include("aberration.jl")
include("bms.jl")
include("transform.jl")
include("gradient.jl")

export transform!,
    transform_objective,
    objective_and_gradient!,
    bms_from_parameters,
    parameters_from_bms,
    ObjectivePlan,
    TransformPlan,
    compute_t′,
    pixel_waveform,
    diagnostics,
    BMS,
    DataComponents,
    Conventions,
    represent!,
    compose,
    conformal_factor,
    lorentz,
    supertranslation,
    frame_rotation,
    boost_velocity,
    translation,
    time_translation,
    space_translation,
    proper_supertranslation
# Note: the single-character accessors `ℐ` and `c_α`, and `ℓₘₐₓ`/`ell_max`, are
# deliberately NOT exported (collision-prone names); use them qualified, as `Scri.ℐ(g)`.

include("precompilation.jl")

end
