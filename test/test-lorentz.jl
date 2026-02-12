# Tests for Lorentz{T}
#
# These are spec tests: they define the interface and invariants that the
# Lorentz{T} type must satisfy, and will fail until that type is implemented.
#
# Strategy: metamorphic testing.  We factor mathematical properties into
# predicate functions in LorentzProperties, then call those predicates for
# many random inputs in separate @testitem blocks.  This separates the
# statement of a property from the test harness.
#
# Assumed public API for Lorentz{T}:
#   Scri.Rotation(R::Quaternionic.Rotor{T})  — construct from a pure rotation
#   Λ₁ * Λ₂                                  — composition (right-to-left)
#   inv(Λ)                                   — group inverse
#   one(Λ)  /  one(Lorentz{T})               — identity element
#   Λ(v::AbstractVector)                     — apply to a 4-vector [t,x,y,z]

# ============================================================================
# Property predicates
# ============================================================================

@testmodule LorentzProperties begin
    import LinearAlgebra: norm
    import Quaternionic

    """Minkowski inner product with signature −+++."""
    function minkowski(v, w)
        return -v[1] * w[1] + v[2] * w[2] + v[3] * w[3] + v[4] * w[4]
    end

    """Minkowski norm squared (signed; negative for timelike vectors)."""
    minkowski_norm²(v) = minkowski(v, v)

    """True iff Λ preserves the Minkowski inner product ⟨v, w⟩."""
    function metric_preserved(Λ, v, w; atol=1e-12)
        return abs(minkowski(Λ(v), Λ(w)) - minkowski(v, w)) ≤ atol
    end

    """True iff (Λ₁ * Λ₂)(v) ≈ Λ₁(Λ₂(v)) — composition is sequential application."""
    function composition_consistent(Λ₁, Λ₂, v; atol=1e-12)
        return norm((Λ₁ * Λ₂)(v) - Λ₁(Λ₂(v))) ≤ atol
    end

    """True iff Λ * inv(Λ) and inv(Λ) * Λ both act as the identity on v."""
    function inverse_law(Λ, v; atol=1e-12)
        return norm((Λ * inv(Λ))(v) - v) ≤ atol && norm((inv(Λ) * Λ)(v) - v) ≤ atol
    end

    """True iff one(Λ) maps every v to itself."""
    function identity_law(Λ, v; atol=1e-12)
        return norm(one(Λ)(v) - v) ≤ atol
    end

    """True iff (Λ₁ * Λ₂) * Λ₃ and Λ₁ * (Λ₂ * Λ₃) agree on v."""
    function associativity_law(Λ₁, Λ₂, Λ₃, v; atol=1e-12)
        return norm(((Λ₁ * Λ₂) * Λ₃)(v) - (Λ₁ * (Λ₂ * Λ₃))(v)) ≤ atol
    end

    """True iff the rotation Λ leaves the time component (index 1) of v unchanged."""
    function preserves_time(Λ, v; atol=1e-12)
        return abs(Λ(v)[1] - v[1]) ≤ atol
    end

    """True iff the rotation Λ preserves the Euclidean norm of the spatial part of v."""
    function preserves_spatial_norm(Λ, v; atol=1e-12)
        return abs(norm(Λ(v)[2:4]) - norm(v[2:4])) ≤ atol
    end

    """
    True iff Rotation(R) and Rotation(−R) produce the same map on v.

    This checks that the double cover Spin⁺(3,1) → SO⁺(3,1) has kernel {±1}:
    both sheets of the cover must give the same Lorentz transformation.
    """
    function double_cover(Λ_pos, Λ_neg, v; atol=1e-12)
        return norm(Λ_pos(v) - Λ_neg(v)) ≤ atol
    end

    """
    True iff the map Spin(3) ∋ R ↦ Rotation(R) ∈ Lorentz is a group homomorphism,
    i.e. Rotation(R₁ * R₂)(v) ≈ (Rotation(R₁) * Rotation(R₂))(v).
    """
    function rotor_homomorphism(Λ₁₂, Λ₁, Λ₂, v; atol=1e-12)
        return norm(Λ₁₂(v) - (Λ₁ * Λ₂)(v)) ≤ atol
    end

    """
    Both real GA norm conditions from the `Quaternionic.components(::Lorentz)` docstring:
      (i)  (R¹)²+(Rˣʸ)²+(Rˣᶻ)²+(Rʸᶻ)² − (Rᵗˣ)²−(Rᵗʸ)²−(Rᵗᶻ)²−(Rᵗˣʸᶻ)² = 1
      (ii) R¹·Rᵗˣʸᶻ − Rʸᶻ·Rᵗˣ + Rˣᶻ·Rᵗʸ − Rˣʸ·Rᵗᶻ = 0
    These are necessary and sufficient for R·R̃ = 1 in the even subalgebra of Cl(3,1).
    """
    function ga_norm_conditions(Λ; atol=1e-12)
        R¹, Rᵗˣ, Rᵗʸ, Rᵗᶻ, Rˣʸ, Rˣᶻ, Rʸᶻ, Rᵗˣʸᶻ = Quaternionic.components(Λ)
        quad = R¹^2 + Rˣʸ^2 + Rˣᶻ^2 + Rʸᶻ^2 - Rᵗˣ^2 - Rᵗʸ^2 - Rᵗᶻ^2 - Rᵗˣʸᶻ^2
        cross = R¹ * Rᵗˣʸᶻ - Rʸᶻ * Rᵗˣ + Rˣᶻ * Rᵗʸ - Rˣʸ * Rᵗᶻ
        return abs(quad - 1) ≤ atol && abs(cross) ≤ atol
    end

    """
    The GA reverse (= group inverse for unit rotors) negates all grade-2 GA components
    and leaves the grade-0 component R¹ and grade-4 component Rᵗˣʸᶻ unchanged.
    """
    function ga_reverse_components(Λ; atol=1e-12)
        c = Quaternionic.components(Λ)
        ci = Quaternionic.components(inv(Λ))
        return abs(ci[1] - c[1]) ≤ atol &&  # R¹     unchanged  (grade 0)
               abs(ci[2] + c[2]) ≤ atol &&  # Rᵗˣ    negated    (grade 2)
               abs(ci[3] + c[3]) ≤ atol &&  # Rᵗʸ    negated    (grade 2)
               abs(ci[4] + c[4]) ≤ atol &&  # Rᵗᶻ    negated    (grade 2)
               abs(ci[5] + c[5]) ≤ atol &&  # Rˣʸ    negated    (grade 2)
               abs(ci[6] + c[6]) ≤ atol &&  # Rˣᶻ    negated    (grade 2)
               abs(ci[7] + c[7]) ≤ atol &&  # Rʸᶻ    negated    (grade 2)
               abs(ci[8] - c[8]) ≤ atol     # Rᵗˣʸᶻ  unchanged  (grade 4)
    end
end

# ============================================================================
# Shared data
# ============================================================================

@testsnippet RotationData begin
    import Quaternionic
    import Random: seed!

    seed!(42)

    const T = Float64
    const n = 20

    # Uniformly random unit rotors from Spin(3), via Haar measure on S³
    rotors = [randn(Quaternionic.Rotor{T}) for _ ∈ 1:n]

    # Corresponding Lorentz{T} elements (rotation subgroup only)
    Λs = [Scri.Rotation(R) for R ∈ rotors]

    # Random general 4-vectors [t, x, y, z]
    general_vecs = [randn(T, 4) for _ ∈ 1:n]

    # Random purely spatial 4-vectors [0, x, y, z]
    spatial_vecs = [[zero(T); randn(T, 3)] for _ ∈ 1:n]
end

@testsnippet LorentzData begin
    import Quaternionic
    import LinearAlgebra: normalize
    import Random: seed!

    seed!(123)

    const T = Float64
    const n = 20

    # Random unit rotors and the corresponding Lorentz elements
    rotors = [randn(Quaternionic.Rotor{T}) for _ ∈ 1:n]
    rot_Λs = Scri.Rotation.(rotors)

    # Random rapidities (bounded away from zero) and random boost directions
    rapidities = abs.(randn(T, n)) .+ T(0.1)
    directions = [normalize(randn(T, 3)) for _ ∈ 1:n]
    boost_Λs = [Scri.Boost(η, n̂) for (η, n̂) ∈ zip(rapidities, directions)]

    # Interleaved sequence rot₁, boost₁, rot₂, boost₂, … and its cumulative products.
    # composed[k] = rot₁ * boost₁ * … * (k-th element): tests norm preservation at every
    # intermediate step.
    mixed_seq = [x for pair ∈ zip(rot_Λs, boost_Λs) for x ∈ pair]
    composed = accumulate(*, mixed_seq)

    # Random general 4-vectors
    general_vecs = [randn(T, 4) for _ ∈ 1:n]
end

# ============================================================================
# Group structure
# ============================================================================

@testitem "Rotation: identity element" tags = [:unit, :fast] setup = [
    LorentzProperties, RotationData
] begin
    for v ∈ general_vecs, Λ ∈ Λs
        @test LorentzProperties.identity_law(Λ, v)
    end
end

@testitem "Rotation: composition matches sequential application" tags = [:unit, :fast] setup = [
    LorentzProperties, RotationData
] begin
    for i ∈ 1:(n - 1)
        Λ₁, Λ₂ = Λs[i], Λs[i + 1]
        for v ∈ general_vecs
            @test LorentzProperties.composition_consistent(Λ₁, Λ₂, v)
        end
    end
end

@testitem "Rotation: inverse" tags = [:unit, :fast] setup = [
    LorentzProperties, RotationData
] begin
    for Λ ∈ Λs, v ∈ general_vecs
        @test LorentzProperties.inverse_law(Λ, v)
    end
end

@testitem "Rotation: associativity" tags = [:unit, :fast] setup = [
    LorentzProperties, RotationData
] begin
    for i ∈ 1:(n - 2)
        Λ₁, Λ₂, Λ₃ = Λs[i], Λs[i + 1], Λs[i + 2]
        for v ∈ general_vecs
            @test LorentzProperties.associativity_law(Λ₁, Λ₂, Λ₃, v)
        end
    end
end

# ============================================================================
# Minkowski isometry
# ============================================================================

@testitem "Rotation: preserves Minkowski metric" tags = [:unit, :fast] setup = [
    LorentzProperties, RotationData
] begin
    for Λ ∈ Λs, i ∈ 1:(n - 1)
        v, w = general_vecs[i], general_vecs[i + 1]
        @test LorentzProperties.metric_preserved(Λ, v, w)
    end
end

@testitem "Rotation: null vectors remain null" tags = [:unit, :fast] setup = [
    LorentzProperties, RotationData
] begin
    for Λ ∈ Λs, v ∈ spatial_vecs
        # Construct a null vector ℓ = [‖v_spatial‖, v_spatial] (future-pointing)
        import LinearAlgebra: norm
        v_sp = v[2:4]
        ℓ = [norm(v_sp); v_sp]
        ℓ′ = Λ(ℓ)
        @test abs(LorentzProperties.minkowski_norm²(ℓ′)) ≤ 1e-12
    end
end

# ============================================================================
# Rotation-specific invariants
# ============================================================================

@testitem "Rotation: preserves time component" tags = [:unit, :fast] setup = [
    LorentzProperties, RotationData
] begin
    for Λ ∈ Λs, v ∈ general_vecs
        @test LorentzProperties.preserves_time(Λ, v)
    end
end

@testitem "Rotation: preserves spatial norm" tags = [:unit, :fast] setup = [
    LorentzProperties, RotationData
] begin
    for Λ ∈ Λs, v ∈ spatial_vecs
        @test LorentzProperties.preserves_spatial_norm(Λ, v)
    end
end

# ============================================================================
# Double cover: Spin⁺(3,1) → SO⁺(3,1)
# ============================================================================

@testitem "Rotation: double cover — R and −R give the same transformation" tags = [
    :unit, :fast
] setup = [LorentzProperties, RotationData] begin
    for (R, v) ∈ zip(rotors, general_vecs)
        Λ_pos = Scri.Rotation(R)
        Λ_neg = Scri.Rotation(-R)
        @test LorentzProperties.double_cover(Λ_pos, Λ_neg, v)
    end
end

# ============================================================================
# Spin(3) → Lorentz is a group homomorphism
# ============================================================================

@testitem "Rotation: Spin(3) → Lorentz{T} is a group homomorphism" tags = [:unit, :fast] setup = [
    LorentzProperties, RotationData
] begin
    for i ∈ 1:(n - 1)
        R₁, R₂ = rotors[i], rotors[i + 1]
        Λ₁₂ = Scri.Rotation(R₁ * R₂)
        Λ₁ = Scri.Rotation(R₁)
        Λ₂ = Scri.Rotation(R₂)
        for v ∈ general_vecs
            @test LorentzProperties.rotor_homomorphism(Λ₁₂, Λ₁, Λ₂, v)
        end
    end
end

# ============================================================================
# Axis-fixing: a rotation about a given axis leaves that axis invariant
# ============================================================================

@testitem "Rotation: rotation about an axis fixes that axis direction" tags = [:unit, :fast] begin
    import Quaternionic
    import LinearAlgebra: norm, normalize

    # Construct an explicit rotation about the z-axis by angle θ.
    # Rotor formula: R = cos(θ/2) + sin(θ/2) k  (k = unit z-quaternion)
    for θ ∈ [0.3, 1.0, π / 2, π, 2π - 0.1]
        w, z_coeff = cos(θ / 2), sin(θ / 2)
        R = Quaternionic.Rotor(w, 0.0, 0.0, z_coeff)
        Λ = Scri.Rotation(R)

        # The z-axis spatial 4-vector [0, 0, 0, 1] must be fixed
        z_vec = [0.0, 0.0, 0.0, 1.0]
        @test Λ(z_vec) ≈ z_vec atol = 1e-12

        # Vectors in the xy-plane should be rotated but stay in the xy-plane
        v_xy = [0.0, 1.0, 0.5, 0.0]
        v_xy′ = Λ(v_xy)
        @test abs(v_xy′[4]) ≤ 1e-12   # z-component stays zero
        @test abs(v_xy′[1]) ≤ 1e-12   # time stays zero
    end

    # Repeat for rotations about the x-axis
    for θ ∈ [0.7, π / 3]
        w, x_coeff = cos(θ / 2), sin(θ / 2)
        R = Quaternionic.Rotor(w, x_coeff, 0.0, 0.0)
        Λ = Scri.Rotation(R)

        x_vec = [0.0, 1.0, 0.0, 0.0]
        @test Λ(x_vec) ≈ x_vec atol = 1e-12
    end

    # Repeat for rotations about the y-axis
    for θ ∈ [1.2, π / 4]
        w, y_coeff = cos(θ / 2), sin(θ / 2)
        R = Quaternionic.Rotor(w, 0.0, y_coeff, 0.0)
        Λ = Scri.Rotation(R)

        y_vec = [0.0, 0.0, 1.0, 0.0]
        @test Λ(y_vec) ≈ y_vec atol = 1e-12
    end
end

# ============================================================================
# Angle composition (metamorphic)
# ============================================================================

@testitem "Rotation: composing with itself doubles the rotation angle" tags = [:unit, :fast] begin
    import Quaternionic

    # R(2θ, n̂) = R(θ, n̂) * R(θ, n̂): rotation by 2θ about the same axis
    # equals composing the rotation by θ with itself.
    for θ ∈ [0.3, 0.9, 1.5, π / 3]
        # Rotation by θ about z-axis
        R_θ = Quaternionic.Rotor(Quaternionic.Quaternion(cos(θ / 2), 0.0, 0.0, sin(θ / 2)))
        R_2θ = Quaternionic.Rotor(Quaternionic.Quaternion(cos(θ), 0.0, 0.0, sin(θ)))  # cos(2θ/2) = cos(θ), sin(2θ/2) = sin(θ)

        Λ_θ = Scri.Rotation(R_θ)
        Λ_2θ = Scri.Rotation(R_2θ)

        v = [0.0, 1.0, 0.0, 0.0]
        @test Λ_2θ(v) ≈ (Λ_θ * Λ_θ)(v) atol = 1e-12
    end
end

# ============================================================================
# 2π rotation is the identity on vectors (but NOT on spinors — Spin group)
# ============================================================================

@testitem "Rotation: 2π rotation acts as the identity on 4-vectors" tags = [:unit, :fast] begin
    import Quaternionic
    import Random: seed!, randn

    seed!(7)

    # R for 2π about z-axis = Rotor(-1, 0, 0, 0): this is the non-trivial
    # element of the kernel of Spin(3) → SO(3), so it acts as identity on vectors.
    R_2π = Quaternionic.Rotor(-1.0, 0.0, 0.0, 0.0)
    Λ_2π = Scri.Rotation(R_2π)

    for _ ∈ 1:10
        v = randn(4)
        @test Λ_2π(v) ≈ v atol = 1e-12
    end
end

# ============================================================================
# Boost constructors and GA algebraic properties
# ============================================================================

@testitem "Boost: explicit GA component values" tags = [:unit, :fast] begin
    import Quaternionic

    # For a pure boost along 𝐧̂ with rapidity η the rotor is
    #   𝐑 = cosh(η/2)·𝟏 + sinh(η/2)·(nˣ𝐭𝐱 + nʸ𝐭𝐲 + nᶻ𝐭𝐳)
    # so only R¹ and the single Rᵗ* component for that axis are nonzero.
    for η ∈ [0.3, 0.7, 1.2, 1.8, 2.5]
        ch, sh = cosh(η / 2), sinh(η / 2)

        # --- boost along z ---
        c = Quaternionic.components(Scri.Boost(η, [0.0, 0.0, 1.0]))
        @test c[1] ≈ ch atol = 1e-14  # R¹
        @test c[2] ≈ 0.0 atol = 1e-14  # Rᵗˣ
        @test c[3] ≈ 0.0 atol = 1e-14  # Rᵗʸ
        @test c[4] ≈ sh atol = 1e-14  # Rᵗᶻ
        @test c[5] ≈ 0.0 atol = 1e-14  # Rˣʸ
        @test c[6] ≈ 0.0 atol = 1e-14  # Rˣᶻ
        @test c[7] ≈ 0.0 atol = 1e-14  # Rʸᶻ
        @test c[8] ≈ 0.0 atol = 1e-14  # Rᵗˣʸᶻ

        # --- boost along x ---
        c = Quaternionic.components(Scri.Boost(η, [1.0, 0.0, 0.0]))
        @test c[1] ≈ ch atol = 1e-14  # R¹
        @test c[2] ≈ sh atol = 1e-14  # Rᵗˣ
        @test c[3] ≈ 0.0 atol = 1e-14  # Rᵗʸ
        @test c[4] ≈ 0.0 atol = 1e-14  # Rᵗᶻ
        @test c[5] ≈ 0.0 atol = 1e-14  # Rˣʸ
        @test c[6] ≈ 0.0 atol = 1e-14  # Rˣᶻ
        @test c[7] ≈ 0.0 atol = 1e-14  # Rʸᶻ
        @test c[8] ≈ 0.0 atol = 1e-14  # Rᵗˣʸᶻ

        # --- boost along y ---
        c = Quaternionic.components(Scri.Boost(η, [0.0, 1.0, 0.0]))
        @test c[1] ≈ ch atol = 1e-14  # R¹
        @test c[2] ≈ 0.0 atol = 1e-14  # Rᵗˣ
        @test c[3] ≈ sh atol = 1e-14  # Rᵗʸ
        @test c[4] ≈ 0.0 atol = 1e-14  # Rᵗᶻ
        @test c[5] ≈ 0.0 atol = 1e-14  # Rˣʸ
        @test c[6] ≈ 0.0 atol = 1e-14  # Rˣᶻ
        @test c[7] ≈ 0.0 atol = 1e-14  # Rʸᶻ
        @test c[8] ≈ 0.0 atol = 1e-14  # Rᵗˣʸᶻ
    end
end

@testitem "Boost: known action on 4-vectors" tags = [:unit, :fast] begin
    # A boost along z by rapidity η maps [t, x, y, z] →
    #   [t·cosh(η) + z·sinh(η),  x,  y,  t·sinh(η) + z·cosh(η)]
    for η ∈ [0.3, 0.7, 1.2, 1.8, 2.5]
        Λ = Scri.Boost(η, [0.0, 0.0, 1.0])
        ch, sh = cosh(η), sinh(η)

        # timelike unit vector: [1, 0, 0, 0] → [cosh(η), 0, 0, sinh(η)]
        @test Λ([1.0, 0.0, 0.0, 0.0]) ≈ [ch, 0.0, 0.0, sh] atol = 1e-13

        # spacelike z-unit: [0, 0, 0, 1] → [sinh(η), 0, 0, cosh(η)]
        @test Λ([0.0, 0.0, 0.0, 1.0]) ≈ [sh, 0.0, 0.0, ch] atol = 1e-13

        # transverse x-unit unchanged: no mixing with t or z
        @test Λ([0.0, 1.0, 0.0, 0.0]) ≈ [0.0, 1.0, 0.0, 0.0] atol = 1e-13

        # future-null vector along z → rescaled by exp(η)
        eη = exp(η)
        @test Λ([1.0, 0.0, 0.0, 1.0]) ≈ eη .* [1.0, 0.0, 0.0, 1.0] atol = 1e-13

        # past-null vector against z → rescaled by exp(−η)
        eη⁻¹ = exp(-η)
        @test Λ([1.0, 0.0, 0.0, -1.0]) ≈ eη⁻¹ .* [1.0, 0.0, 0.0, -1.0] atol = 1e-13
    end
end

@testitem "Boost: collinear rapidities add (metamorphic)" tags = [:unit, :fast] begin
    import LinearAlgebra: norm

    # For boosts along the same axis, rapidities are additive:
    #   Boost(η₁, n̂) * Boost(η₂, n̂) acts identically to Boost(η₁+η₂, n̂).
    for (η₁, η₂) ∈ [(0.3, 0.5), (1.0, 1.5), (0.1, 2.0), (0.7, 0.7)]
        for n̂ ∈ [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            Λ_prod = Scri.Boost(η₁, n̂) * Scri.Boost(η₂, n̂)
            Λ_sum = Scri.Boost(η₁ + η₂, n̂)
            for v ∈ [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.5, -0.3], [1.0, 0.2, 0.3, 0.9]]
                @test norm(Λ_prod(v) - Λ_sum(v)) ≤ 1e-12
            end
        end
    end
end

@testitem "Lorentz: GA norm conditions hold for pure transformations" tags = [:unit, :fast] setup = [
    LorentzProperties, LorentzData
] begin
    for Λ ∈ rot_Λs
        @test LorentzProperties.ga_norm_conditions(Λ)
    end
    for Λ ∈ boost_Λs
        @test LorentzProperties.ga_norm_conditions(Λ)
    end
end

@testitem "Lorentz: GA norm conditions preserved under composition (metamorphic)" tags = [
    :unit, :fast
] setup = [LorentzProperties, LorentzData] begin
    # Each element of `composed` is a prefix product of alternating random rotations and
    # boosts.  Since each factor satisfies R·R̃ = 1, so must every prefix product — this
    # is the "norms multiply" property of the quaternion algebra over ℂ.
    for Λ ∈ composed
        @test LorentzProperties.ga_norm_conditions(Λ)
    end
end

@testitem "Lorentz: GA reverse is the group inverse (component test)" tags = [:unit, :fast] setup = [
    LorentzProperties, LorentzData
] begin
    # The GA reverse ~R negates all grade-2 components and fixes grades 0 and 4.
    # Since inv(Λ) is implemented via conj(rotor(Λ)) = ~R, this must hold for all Λ.
    for Λ ∈ rot_Λs
        @test LorentzProperties.ga_reverse_components(Λ)
    end
    for Λ ∈ boost_Λs
        @test LorentzProperties.ga_reverse_components(Λ)
    end
    for Λ ∈ composed
        @test LorentzProperties.ga_reverse_components(Λ)
    end
end

@testitem "Boost: group properties" tags = [:unit, :fast] setup = [
    LorentzProperties, LorentzData
] begin
    # identity
    for Λ ∈ boost_Λs, v ∈ general_vecs
        @test LorentzProperties.identity_law(Λ, v)
    end
    # composition = sequential application
    for i ∈ 1:(n - 1), v ∈ general_vecs
        @test LorentzProperties.composition_consistent(boost_Λs[i], boost_Λs[i + 1], v)
    end
    # inverse
    for Λ ∈ boost_Λs, v ∈ general_vecs
        @test LorentzProperties.inverse_law(Λ, v)
    end
    # associativity
    for i ∈ 1:(n - 2), v ∈ general_vecs
        @test LorentzProperties.associativity_law(
            boost_Λs[i], boost_Λs[i + 1], boost_Λs[i + 2], v
        )
    end
    # metric preservation
    for Λ ∈ boost_Λs, i ∈ 1:(n - 1)
        v, w = general_vecs[i], general_vecs[i + 1]
        @test LorentzProperties.metric_preserved(Λ, v, w)
    end
end
