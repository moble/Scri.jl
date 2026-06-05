module Scri

import Quaternionic: QuatVec, Rotor, absvec, 𝐤, value, components, basetype
import SphericalFunctions: ₛ𝐘, ð, golden_ratio_spiral_rotors
import LinearAlgebra: mul!, ldiv!, lu, I, qr
import OffsetArrays: OffsetVector
import Logging
import Hwloc
import Polyester
import OhMyThreads
import Base.Threads: nthreads

# These are just for precompilation
using PrecompileTools: @setup_workload, @compile_workload
using Random: Xoshiro

const cachesize_L2 = Logging.with_logger(Logging.NullLogger()) do
    Hwloc.cachesize(:L2)
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
include("data_components.jl")
include("utilities.jl")
include("aberration.jl")
include("transform.jl")

export transform!, diagnostics

include("precompilation.jl")

end
