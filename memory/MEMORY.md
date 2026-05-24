# Scri.jl Project Memory

## Status

`Lorentz{T}` (rotation subgroup) is implemented and fully tested (3558 tests
passing).  Stores `Quaternionic.Rotor{Complex{T}}` internally — ready for the
boost sector.  Boost constructors not yet implemented.

## Architecture (from notes/plan.md and docs/)

- `AbstractBMS{T}` — root
  - `BMS{T}` — full group element: `α::Supertranslation{T}`, `Λ::Lorentz{T}`
  - `Supertranslation{T}` — αₗₘ modes; includes `SpacetimeTranslation{T}`
  - `Lorentz{T}` — spinor/rotor `R`; constructed via `Rotation(...)` or `Boost(...)`
- `Lorentz{T}` may internally use `Quaternion{Complex{T}}` (biquaternion) — TBD
- Only runtime dep: `Quaternionic = "3.1.1"` (unit quaternion arithmetic)

## Lorentz{T} implementation notes

- Internal rep: `Quaternionic.Rotor{Complex{T}}` (biquaternion); group condition
  `w²+x²+y²+z²=1` over ℂ.  Rotations: all real.  Boosts: imaginary parts.
- SL(2ℂ) matrix: `M = [[w+iz, y+ix], [-y+ix, w-iz]]`
- 4-vector as 2×2 Hermitian: `H = [[t+vz, vx-ivy], [vx+ivy, t-vz]]`
- **Action is `H' = M†HM`** (NOT `M H M†`).  The map `R → M` is an
  anti-homomorphism (`M(R₁R₂) = M(R₂)M(R₁)`), so `M†` conjugates the order
  back to give `Λ_{R₁R₂}(v) = Λ_{R₁}(Λ_{R₂}(v))`.
- `inv(Λ)` = `conj(rotor(Λ))` (quaternionic conjugate).

- `Scri.Rotation(R::Quaternionic.Rotor{T})` — construct rotation
- `Λ₁ * Λ₂` — composition (right-to-left, standard)
- `inv(Λ)` — group inverse
- `one(Λ)` / `one(Lorentz{T})` — identity
- `Λ(v::AbstractVector)` — apply to 4-vector `[t,x,y,z]`

## Key files

- `src/Scri.jl` — main source (stub)
- `test/test-lorentz.jl` — Lorentz{T} spec tests (metamorphic, property-based)
- `test/test-basic-test.jl` — placeholder tests
- `docs/src/lorentz_group.md` — null tetrad, null rotation formulas in GA
- `docs/src/bms_group.md` — BMS group law, conformal factor K, decomposition
- `notes/plan.md` — full type hierarchy and BMS algorithm sketch

## Physics conventions (CLAUDE.md)

- Metric: −+++ (Minkowski `ℝ^{3,1}`)
- `𝐭² = −1`, `𝐱² = 𝐲² = 𝐳² = +1`
- Null tetrad: ℓ = (𝐭+𝐳)/√2, 𝐧 = (𝐭−𝐳)/√2, 𝐦 = (𝐱+𝐈₃𝐲)/√2
- Rotors act by conjugation: `𝐯′ = R𝐯R̃`
- Use `Spin⁺(3,1)` — never SU(2), SL(2,C)
- NP quantities are "Weyl *components*", not "Weyl scalars"

## Code style

- Unicode everywhere; `for x ∈` not `for x in`; `JuliaFormatter` blue style
- Two spaces between sentences in prose/comments
- Tests: `@testitem` blocks, TestItemRunner, run with `julia --project=test test/runtests.jl`
