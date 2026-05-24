module Scri

import Quaternionic: QuatVec, Rotor, absvec, abs2vec, 𝐤, ×, ×̂, value, components, basetype
import SphericalFunctions: SSHTDirect, ₛ𝐘, ð, ð̄, golden_ratio_spiral_rotors
import LinearAlgebra: mul!, ldiv!, lu, I, BLAS, qr
import OffsetArrays: OffsetVector
import Logging
import Hwloc
import Polyester
import OhMyThreads
import Base.Threads: nthreads

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

include("utilities.jl")
include("cubic_spline.jl")
include("data_components.jl")
include("aberration.jl")
include("transform.jl")

export transform!

end
