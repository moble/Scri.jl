module ScriPlotsExt

import Scri
import Plots
import Colors

const red = Colors.colorant"red"
const yellow = Colors.colorant"yellow"
const green = Colors.colorant"green"

function Scri.diagnostics(t, data, data_components)
    d = Scri.diagnostics(data, data_components)
    plots = Dict{Symbol, Plots.Plot}()
    for (comp, values) ∈ d
        s = Scri.spin_weight(comp)
        ℓₘₐₓ = size(values, 2) - 1
        colors = reshape(
            [
                [red for _ ∈ 0:abs(s)-1];
                range(green, stop=yellow, length=ℓₘₐₓ-abs(s)+1)
            ],
            1, :,
        )
        labels = reshape([string(ℓ) for ℓ ∈ 0:ℓₘₐₓ], 1, :)
        p = Plots.plot(
            t, sqrt.(values), linecolors=colors, label=labels,
            yscale=:log10, title=string(comp), xlabel="Time", ylabel="Energy"
        )
        plots[comp] = p
    end
    return plots
end

end  # module ScriPlotsExt
