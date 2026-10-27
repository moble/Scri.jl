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
    SphericalFunctions, sYlm_matrix, sYlmCalculator, set_R!, ð, leja_rotors, SSHTMatrix, D
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
# samples `transform!` handles in each block.  Calling `cachesize_L2()` reads it from Hwloc
# once per process, on the first call; concurrent first calls wait for that one reading.
# A value read while precompiling is not stored in the package image, so a depot shared
# between machines does not carry the build machine's value.
const cachesize_L2 = OncePerProcess{Int}() do
    Logging.with_logger(Logging.NullLogger()) do
        return Hwloc.cachesize(:L2)
    end
end

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
