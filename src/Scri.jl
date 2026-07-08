module Scri

import Quaternionic
import Quaternionic:
    QuatVec,
    Rotor,
    Lorentz,
    Boost,
    absvec,
    𝐤,
    components,
    basetype,
    from_spherical_coordinates,
    ℂreal
import SphericalFunctions: ₛ𝐘, ð, golden_ratio_spiral_rotors
import LinearAlgebra: LinearAlgebra, mul!, ldiv!, lu, qr
import OffsetArrays: OffsetVector
import Logging
import Hwloc
import Polyester
import OhMyThreads
import Base.Threads: nthreads
import TestItems: @testitem, @testmodule

# These are just for precompilation
using PrecompileTools: @setup_workload, @compile_workload
using Random: Xoshiro

const cachesize_L2 = Logging.with_logger(Logging.NullLogger()) do
    return Hwloc.cachesize(:L2)
end
const cachesize_L3 = Logging.with_logger(Logging.NullLogger()) do
    try
        try
            Hwloc.cachesize(:L3)
        catch
            Hwloc.gettopology().mem
        end
    catch
        cachesize_L2
    end
end

include("cubic_spline.jl")
include("signs.jl")
include("conventions.jl")
include("data_components.jl")
include("utilities.jl")
include("aberration.jl")
include("bms.jl")
include("transform.jl")

export transform!, diagnostics, BMS

include("precompilation.jl")

end
