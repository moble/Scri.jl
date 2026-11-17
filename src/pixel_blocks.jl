# Blocked synthesis of pixel time series from mode data.
#
# The objective and its gradient need the data's time series at many pixels.  Synthesizing
# each pixel on its own, with one dot product per time sample, reads the whole mode array
# once per pixel, and that dominates their cost.  Instead the pixels are taken in blocks.
# For each block, parallel tasks first evaluate the harmonics at every pixel; then one
# matrix product per component synthesizes the time series of the whole block, using BLAS's
# own threads; then parallel tasks process each pixel's series.  The stages run one after
# another and are never nested, so the number of BLAS threads never needs to be changed.
# The threading notes in the documentation explain why that matters.

"""
    pixel_blocks(Nᵖ, bytes_per_pixel, budget)

Split the pixels `1:Nᵖ` into consecutive blocks, returned as a vector of ranges.  Each
block holds as many pixels as fit in `budget` bytes, at `bytes_per_pixel` bytes each, but
at least four for each thread, so that the parallel stages have enough work to share.

The budget only limits memory.  Larger blocks are never slower: BLAS does its own cache
blocking within each product, and fewer, larger products spread better over BLAS's
threads, with fewer switches between the stages.  So a budget that fits a block in cache
gains nothing.
"""
function pixel_blocks(Nᵖ::Integer, bytes_per_pixel::Integer, budget::Integer)
    B = max(1, min(Nᵖ, max(4 * nthreads(), budget ÷ max(1, bytes_per_pixel))))
    return [b:min(b + B - 1, Nᵖ) for b ∈ 1:B:Nᵖ]
end

"""
    foreach_chunk(f, n, buffers)

Split `1:n` into at most `length(buffers)` consecutive chunks, and call `f(buffers[c],
chunk)` for each chunk, numbered `c`, in a parallel task of its own.  The buffers are
indexed by chunk, never by thread, so no two tasks share one, and they can be allocated
once and reused for every block.
"""
function foreach_chunk(f, n::Integer, buffers::AbstractVector)
    n > 0 || return nothing
    chunks = OhMyThreads.index_chunks(1:n; n=min(n, length(buffers)))
    OhMyThreads.tforeach(enumerate(chunks); chunking=false) do (c, chunk)
        return f(buffers[c], chunk)
    end
    return nothing
end

"""
    harmonics_block!(Y, calcs, Rₚ, block)

Evaluate the spin-weighted spherical harmonics at the rotors `Rₚ[block]`, in parallel tasks
that each use one of the calculators `calcs`.  The harmonics at `Rₚ[block[q]]` are stored
in `Y[:, q, :]`, arranged as [`sYlm_rows!`](@ref) arranges them, with one slice of the last
dimension for each spin weight that the calculators serve.
"""
function harmonics_block!(Y::AbstractArray{<:Any,3}, calcs, Rₚ, block)
    foreach_chunk(length(block), calcs) do calc, chunk
        for q ∈ chunk
            sYlm_rows!(view(Y, :, q, :), calc, Rₚ[block[q]])
        end
    end
    return Y
end

@doc raw"""
    synthesize_block!(σ, data, Y, spins, N₀)

Synthesize the time series of a block of pixels from the modes in `data`, which is only
read:
```math
σ[h, q, k] = ∑_{m=s_k^2+1}^{N₀} Y[m, q, s_k - s₀ + 1]\, \mathrm{data}[m, h, k],
```
where ``s_k`` is `spins[k]`, the spin weight of component `k`, and ``s₀`` is the least of
them.  So `Y[:, q, :]` holds the harmonics at pixel `q` as [`harmonics_block!`](@ref)
leaves them, and `N₀` is the number of modes within the data's band limit.  The modes with
``ℓ < |s_k|`` are skipped, as everywhere in [`transform!`](@ref).

When `data` and `Y` hold the same BLAS element type, there is one matrix product for each
component, which uses BLAS's own threads.  So there is when `Y` holds ForwardDiff dual
numbers on the element type of `data`, by way of [`stack_partials`](@ref).  Otherwise each
sum is a [`pixel_dot`](@ref), and the pixels are shared among parallel tasks.
"""
function synthesize_block!(
    σ::StridedArray{T,3}, data::StridedArray{T,3}, Y::StridedArray{T,3}, spins, N₀
) where {T<:LinearAlgebra.BlasComplex}
    s₀ = minimum(spins)
    for k ∈ axes(σ, 3)
        s = spins[k]
        m = (s ^ 2 + 1):N₀
        mul!(view(σ, :, :, k), transpose(view(data, m, :, k)), view(Y, m, :, s - s₀ + 1))
    end
    return σ
end
function synthesize_block!(
    σ::StridedArray{Complex{D},3},
    data::StridedArray{Complex{V},3},
    Y::StridedArray{Complex{D},3},
    spins,
    N₀,
) where {V<:LinearAlgebra.BlasReal,D<:ForwardDiff.Dual{<:Any,V}}
    # Dual harmonics with primal data: the values and partials of the harmonics of each
    # spin weight are stacked once and shared by the components of that spin weight.
    s₀ = minimum(spins)
    for s ∈ unique(spins)
        m = (s ^ 2 + 1):N₀
        Yˢ = stack_partials(view(Y, m, :, s - s₀ + 1))
        for k ∈ findall(==(s), spins)
            unstack_partials!(view(σ, :, :, k), transpose(view(data, m, :, k)) * Yˢ, false)
        end
    end
    return σ
end
function synthesize_block!(σ, data, Y, spins, N₀)
    s₀ = minimum(spins)
    OhMyThreads.tforeach(axes(σ, 2)) do q
        for k ∈ axes(σ, 3)
            s = spins[k]
            Yˢ = view(Y, :, q, s - s₀ + 1)
            @inbounds for h ∈ axes(σ, 1)
                σ[h, q, k] = pixel_dot(Yˢ, view(data, :, h, k), s ^ 2 + 1, N₀)
            end
        end
    end
    return σ
end

"""
    synthesis_blocks(data, Rₚ, spins; extra_bytes_per_pixel=0,
                     budget=synthesis_budget(data))

Split the pixels at the rotors `Rₚ` into blocks with [`pixel_blocks`](@ref), counting the
memory that [`harmonics_block!`](@ref) and [`synthesize_block!`](@ref) need for each pixel
of a block, plus `extra_bytes_per_pixel` for whatever else the caller stores per pixel,
within `budget` bytes.  The default budget is the size of `data` itself, or 128 MiB if
that is larger.  So for large data the blocks scale with the problem and take no more
memory than the mode array the caller already holds, while for short time series, whose
pixels' harmonics take as much room as their series, the blocks are not made needlessly
small.  Tests pass a tiny budget to get several blocks from a small problem.
"""
function synthesis_blocks(
    data,
    Rₚ,
    spins;
    extra_bytes_per_pixel::Integer=0,
    budget::Integer=synthesis_budget(data),
)
    Nᵐ, Nᵗ, Nᵈ = size(data)
    Tʸ = Complex{basetype(eltype(Rₚ))}
    Tˣ = promote_type(eltype(data), Tʸ)
    Nˢ = maximum(spins) - minimum(spins) + 1
    bytes = Nᵐ * Nˢ * sizeof(Tʸ) + Nᵗ * Nᵈ * sizeof(Tˣ) + extra_bytes_per_pixel
    return pixel_blocks(length(Rₚ), bytes, budget)
end

"""
    synthesis_budget(data)

The default memory budget of [`synthesis_blocks`](@ref), in bytes: the size of `data`, or
128 MiB if that is larger.
"""
synthesis_budget(data) = max(sizeof(data), 2^27)

"""
    block_workspace(blocks, data, t′, Rₚ, spins, TR)

Allocate the memory that [`transform_pixel_blocks`](@ref) needs to transform the data at
the pixels whose rotors are `Rₚ`, a block at a time, for the blocks `blocks`, as the
named tuple `(; Y, σ, calcs, buffers)`.  `Y` holds the harmonics of a block, and `σ` the
time series synthesized from them.  `calcs` holds one harmonics calculator for each
parallel task, and `buffers` one tuple `(σⱼ, σ̈ⱼ, μⱼ)` for each parallel task: one
pixel's time series, its spline's second derivatives, and its transformed series, of
element type `Complex{TR}`.  A caller that already holds memory of these shapes and types
may pass its own named tuple instead.
"""
function block_workspace(blocks, data, t′, Rₚ, spins, ::Type{TR}) where {TR}
    Nᵐ, Nᵗ, Nᵈ = size(data)
    Nᵗ′ = length(t′)
    ℓₘₐₓ = isqrt(Nᵐ) - 1
    spin_range = minimum(spins):maximum(spins)
    B = maximum(length, blocks)
    Y = Array{Complex{basetype(eltype(Rₚ))}}(undef, Nᵐ, B, length(spin_range))
    σ = Array{promote_type(eltype(data), eltype(Y))}(undef, Nᵗ, B, Nᵈ)
    calcs = [sYlmCalculator(one(eltype(Rₚ)), ℓₘₐₓ, spin_range) for _ ∈ 1:nthreads()]
    buffers = [
        (
            Matrix{Complex{TR}}(undef, Nᵈ, Nᵗ),
            Matrix{Complex{TR}}(undef, Nᵈ, Nᵗ),
            Matrix{Complex{TR}}(undef, Nᵈ, Nᵗ′),
        ) for _ ∈ 1:nthreads()
    ]
    return (; Y, σ, calcs, buffers)
end

"""
    transform_pixel_blocks(on_pixel!, on_block!, blocks, data, t, t′, cubic_spline_cache,
                           geometry, dc, spins, N₀, TR)
    transform_pixel_blocks(on_pixel!, on_block!, blocks, data, t, t′, cubic_spline_cache,
                           geometry, dc, spins, N₀, workspace)

Transform the data at every pixel, one block of pixels at a time.  The named tuple
`geometry = (; Rₚ, κ⁻¹ₚ, αₚ, ðt′╱2κₚ, ð²αₚ)` holds each pixel's rotor in the rest frame
and the scalars that [`transform_pixel!`](@ref) needs there, as `transform_objective`
computes them.

For each block in `blocks`, the time series of its pixels are synthesized with
[`harmonics_block!`](@ref) and [`synthesize_block!`](@ref).  Then parallel tasks transform
each pixel `j = block[q]` with `transform_pixel!`, into a buffer `μⱼ` of element type
`Complex{TR}`, and call `on_pixel!(μⱼ, j, q)`.  The buffer is reused for the next pixel,
so `on_pixel!` must copy anything it keeps.  When every pixel of the block is done,
`on_block!(block)` is called.

The memory for all of this is allocated once, by [`block_workspace`](@ref), and used for
every block.  Pass a `workspace` of that form in place of `TR` to use memory the caller
already holds.  When the call returns, the workspace still holds the harmonics and the
synthesized series of the last block.
"""
function transform_pixel_blocks(
    on_pixel!,
    on_block!,
    blocks,
    data,
    t,
    t′,
    cubic_spline_cache::CubicSplineCache,
    geometry,
    dc::DataComponents,
    spins,
    N₀,
    ::Type{TR},
) where {TR}
    workspace = block_workspace(blocks, data, t′, geometry.Rₚ, spins, TR)
    return transform_pixel_blocks(
        on_pixel!,
        on_block!,
        blocks,
        data,
        t,
        t′,
        cubic_spline_cache,
        geometry,
        dc,
        spins,
        N₀,
        workspace,
    )
end
function transform_pixel_blocks(
    on_pixel!,
    on_block!,
    blocks,
    data,
    t,
    t′,
    cubic_spline_cache::CubicSplineCache,
    geometry,
    dc::DataComponents,
    spins,
    N₀,
    workspace::NamedTuple,
)
    (; Rₚ, κ⁻¹ₚ, αₚ, ðt′╱2κₚ, ð²αₚ) = geometry
    (; Y, σ, calcs, buffers) = workspace
    for block ∈ blocks
        nb = length(block)
        Yb, σb = view(Y, :, 1:nb, :), view(σ, :, 1:nb, :)
        harmonics_block!(Yb, calcs, Rₚ, block)
        synthesize_block!(σb, data, Yb, spins, N₀)
        foreach_chunk(nb, buffers) do (σⱼ, σ̈ⱼ, μⱼ), chunk
            for q ∈ chunk
                j = block[q]
                permutedims!(σⱼ, view(σb, :, q, :), (2, 1))
                transform_pixel!(
                    μⱼ,
                    σⱼ,
                    σ̈ⱼ,
                    cubic_spline_cache,
                    t,
                    t′,
                    κ⁻¹ₚ[j],
                    αₚ[j],
                    ðt′╱2κₚ[1, j],
                    ðt′╱2κₚ[2, j],
                    ð²αₚ[j],
                    dc,
                )
                on_pixel!(μⱼ, j, q)
            end
        end
        on_block!(block)
    end
    return nothing
end

"""
    transformed_projection!(μ̂, μᵇ, blocks, data, t, t′, cubic_spline_cache, geometry, dc,
                            spins, N₀, Ȳ, workspace)

Project the transformed waveform onto the modes up to an `ObjectivePlan`'s `ℓ_fit` by its
quadrature.  The result is stored in `μ̂`, which holds one matrix for each component `k`:
```
μ̂[k][n, i] = Σⱼ Ȳ[k][n, j] μⱼ[k, i],
```
where `μⱼ` is pixel `j`'s transformed time series, from
[`transform_pixel_blocks`](@ref), and `Ȳ` is the plan's weighted conjugate harmonics.  The
other arguments are as for `transform_pixel_blocks`, with the blocks `blocks` and the
`workspace` that it takes.  The transformed series of each block are gathered into
`μᵇ[i, q, k]`, with time first, so that each pixel's series fills contiguous memory, and
one matrix product for each component adds the block's share to `μ̂[k]`.  So `μᵇ` must
have at least as many columns as the largest block.  Returns `μ̂`.
"""
function transformed_projection!(
    μ̂, μᵇ, blocks, data, t, t′, cubic_spline_cache, geometry, dc, spins, N₀, Ȳ, workspace
)
    Nᵈ = size(data, 3)
    for k ∈ 1:Nᵈ
        fill!(μ̂[k], 0)
    end
    gather!(μⱼ, j, q) = permutedims!(view(μᵇ, :, q, :), μⱼ, (2, 1))
    function project!(block)
        for k ∈ 1:Nᵈ
            μᵇₖ = transpose(view(μᵇ, :, 1:length(block), k))
            add_product!(μ̂[k], view(Ȳ[k], :, block), μᵇₖ)
        end
    end
    transform_pixel_blocks(
        gather!,
        project!,
        blocks,
        data,
        t,
        t′,
        cubic_spline_cache,
        geometry,
        dc,
        spins,
        N₀,
        workspace,
    )
    return μ̂
end

"""
    transformed_projection(data, t, t′, cubic_spline_cache, geometry, dc, spins, N₀, Ȳ, TR)

Return the projection that [`transformed_projection!`](@ref) computes, as a vector of
matrices of element type `Complex{TR}`, in new memory and with blocks of its own choosing.
"""
function transformed_projection(
    data, t, t′, cubic_spline_cache, geometry, dc, spins, N₀, Ȳ, ::Type{TR}
) where {TR}
    Nᵈ = size(data, 3)
    Nᵗ′ = length(t′)
    μ̂ = [Matrix{Complex{TR}}(undef, size(Ȳ[k], 1), Nᵗ′) for k ∈ 1:Nᵈ]
    blocks = synthesis_blocks(
        data, geometry.Rₚ, spins; extra_bytes_per_pixel=Nᵗ′ * Nᵈ * sizeof(Complex{TR})
    )
    μᵇ = Array{Complex{TR}}(undef, Nᵗ′, maximum(length, blocks), Nᵈ)
    workspace = block_workspace(blocks, data, t′, geometry.Rₚ, spins, TR)
    return transformed_projection!(
        μ̂,
        μᵇ,
        blocks,
        data,
        t,
        t′,
        cubic_spline_cache,
        geometry,
        dc,
        spins,
        N₀,
        Ȳ,
        workspace,
    )
end

"""
    stack_partials(B)

Return the matrix `[value(B) ∂₁B ∂₂B … ∂ₙB]` of plain complex numbers, for a matrix `B` of
complex numbers whose real and imaginary parts are ForwardDiff dual numbers with `N`
partials.  The block `∂ₚB` holds the `p`-th partials of the real and imaginary parts as the
real and imaginary parts of one complex number.  A product `A * B` with `A` free of dual
numbers is linear in `B`, so `A * stack_partials(B)` holds the value and every partial of
the product, side by side, and [`unstack_partials!`](@ref) reassembles them.  That way one
BLAS product replaces the generic product of dual numbers.
"""
function stack_partials(
    B::AbstractMatrix{Complex{D}}
) where {V,N,D<:ForwardDiff.Dual{<:Any,V,N}}
    n = size(B, 2)
    Bˢ = Matrix{Complex{V}}(undef, size(B, 1), n * (N + 1))
    @inbounds for c ∈ 1:n, r ∈ axes(B, 1)
        x, y = reim(B[r, c])
        Bˢ[r, c] = complex(ForwardDiff.value(x), ForwardDiff.value(y))
        ∂x, ∂y = ForwardDiff.partials(x), ForwardDiff.partials(y)
        for p ∈ 1:N
            Bˢ[r, c + n * p] = complex(∂x[p], ∂y[p])
        end
    end
    return Bˢ
end

"""
    unstack_partials!(C, Cˢ, add)

Store in `C` the complex dual numbers whose values and partials are the blocks of `Cˢ`,
laid out as by [`stack_partials`](@ref), or add them to `C` if `add` is true.
"""
function unstack_partials!(
    C::AbstractMatrix{Complex{D}}, Cˢ::AbstractMatrix, add::Bool
) where {V,N,D<:ForwardDiff.Dual{<:Any,V,N}}
    n = size(C, 2)
    @inbounds for c ∈ 1:n, r ∈ axes(C, 1)
        z = Cˢ[r, c]
        x = D(real(z), ForwardDiff.Partials(ntuple(p -> real(Cˢ[r, c + n * p]), Val(N))))
        y = D(imag(z), ForwardDiff.Partials(ntuple(p -> imag(Cˢ[r, c + n * p]), Val(N))))
        C[r, c] = add ? C[r, c] + complex(x, y) : complex(x, y)
    end
    return C
end

"""
    add_product!(C, A, B)

Add the matrix product `A * B` to `C`.  When `A` holds plain complex numbers of a BLAS type
and `B` complex dual numbers on that type, the product is one BLAS product, through
[`stack_partials`](@ref); otherwise it is `mul!(C, A, B, true, true)`.
"""
add_product!(C, A, B) = mul!(C, A, B, true, true)
function add_product!(
    C::AbstractMatrix{Complex{D}},
    A::StridedMatrix{Complex{V}},
    B::AbstractMatrix{Complex{D}},
) where {V<:LinearAlgebra.BlasReal,D<:ForwardDiff.Dual{<:Any,V}}
    return unstack_partials!(C, A * stack_partials(B), true)
end
