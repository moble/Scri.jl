# One big problem with the `transform!` function is that it's important to specialize on the
# `data_components` type.  But that's so flexible, there are many possibilities.  Here, I've
# selected a few of the most likely ones I expect to be used in practice.  But it's already
# a lot, and presumably leads to slower load times for the package.  I don't actually know
# that precompilation is *all that* necessary, but this is the tradeoff I've chosen for now.

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the
    # size of the precompile file and potentially make loading faster.  Nothing defined in
    # this section will be defined outside of this block.
    ℓₘₐₓ = 4
    Nᵐ = (ℓₘₐₓ+1)^2
    Nᵗ = 17
    t₁, t₂ = -50.0, 50.0
    rng = Xoshiro(1234)
    t = collect(LinRange(t₁, t₂, Nᵗ))
    v⃗ = 1e-8 * randn(rng, QuatVec{Float64})
    R = randn(rng, Rotor{Float64})
    αᵢₙ = 1e-8 * randn(rng, Complex{Float64}, Nᵐ)

    @compile_workload begin
        # All calls in this block will be precompiled, regardless of whether they belong to
        # this package or not (on Julia 1.8 and higher).

        for data_components ∈ (
            (:ψ₄,), (:ψ₄, :ψ₃), (:ψ₄, :ψ₃, :ψ₂), (:ψ₄, :ψ₃, :ψ₂, :ψ₁),
            (:ψ₄, :ψ₃, :ψ₂, :ψ₁, :ψ₀),
            (:σ,), (:σ, :ψ₄), (:σ, :ψ₄, :ψ₃), (:σ, :ψ₄, :ψ₃, :ψ₂), (:σ, :ψ₄, :ψ₃, :ψ₂, :ψ₁),
            (:σ, :ψ₄, :ψ₃, :ψ₂, :ψ₁, :ψ₀),
            (:h,), (:h, :ψ₄), (:h, :ψ₄, :ψ₃), (:h, :ψ₄, :ψ₃, :ψ₂), (:h, :ψ₄, :ψ₃, :ψ₂, :ψ₁),
            (:h, :ψ₄, :ψ₃, :ψ₂, :ψ₁, :ψ₀),
            (:News,),
        )
            Nᵈ = length(data_components)
            data = randn(rng, Complex{Float64}, Nᵐ, Nᵗ, Nᵈ)
            Scri.transform!(data, t, v⃗, R, αᵢₙ; data_components)
        end
    end
end
