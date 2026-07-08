#!/usr/bin/env julia
# examples/rotation_sphere.jl
#
# Visualizes the effect of a passive frame rotation on the celestial sphere.
# An "F"-shaped ψ₂ field (spin-0) is expanded in spherical harmonics, then
# transformed with transform! (pure rotation, zero boost, zero supertranslation),
# and the result is plotted on 3-D spheres.
#
# Panel 1: observer A's natural view (camera at −y, looking toward +y).
#   x → right, y → into screen, z → up.  F reads correctly from outside.
# Panel 2: observer B's natural view (camera at −x, looking toward +x).
#   x′ → right (= A's −y), y′ → into screen (= A's +x), z′ → up.
#   The rotated F is face-on from B's perspective too.
#
# Run from the Scri.jl root:
#   julia --project=examples examples/rotation_sphere.jl

using LinearAlgebra: I, lu
using PythonPlot
using PythonCall: pyimport, pylist, pyexec, pydict
import Quaternionic:
    Quaternionic,
    rotor,
    components,
    Rotor,
    QuatVec,
    normalize,
    from_spherical_coordinates,
    𝐢,
    𝐣,
    𝐤
import SphericalFunctions: ₛ𝐘, golden_ratio_spiral_rotors
import Scri
import Scri: transform!

# ─────────────────────────────────────────────────────────────────────────────
# F-shaped mask
# ─────────────────────────────────────────────────────────────────────────────

function in_F(x, y, z)
    y < 0.12 && return false
    in_vert = -0.32 < -x < -0.04 && -0.10 < z < 0.68
    in_top = -0.32 < -x < 0.42 && 0.52 < z < 0.68
    in_mid = -0.32 < -x < 0.22 && 0.15 < z < 0.30
    return in_vert || in_top || in_mid
end
function in_F(R)
    θ, ϕ = Quaternionic.to_spherical_coordinates(R)
    in_vert = 1π/12 < θ < 6π/12 && 0 < ϕ < 1π/8
    in_top = 1π/12 < θ < 2π/12 && 0 < ϕ < 4π/8
    in_mid = 3π/12 < θ < 4π/12 && 0 < ϕ < 3π/8
    return in_vert || in_top || in_mid
end

# ─────────────────────────────────────────────────────────────────────────────
# Build ψ₂ mode weights via SSHT
# ─────────────────────────────────────────────────────────────────────────────

const ℓₘₐₓ = 32
const N_modes = (ℓₘₐₓ + 1)^2

# For a step-function mask (discontinuous, not band-limited), a square LU solve
# does exact interpolation that oscillates wildly off the sample grid (Runge
# phenomenon).  Use an overdetermined grid (ℓ_samp ≈ 2ℓₘₐₓ → ~4× more points)
# so that \ does a least-squares projection, giving proper Gibbs-bounded modes.
const ℓ_samp = 2 * ℓₘₐₓ   # (ℓ_samp+1)² ≈ 4 × N_modes sample points

println("Computing SSHT sample grid…")
const R_f = exp(Quaternionic.QuatVecF64(0.0, 0.0, 1.0) * (π/4))
const Rₛ = golden_ratio_spiral_rotors(0, ℓ_samp, Float64)
const mask_vals = [in_F(vec(R(𝐤))...) for R ∈ Rₛ]
#const mask_vals = [in_F(R) for R ∈ Rₛ]
const Y_samp = ₛ𝐘(0, ℓₘₐₓ, Float64, Rₛ)         # rectangular: (4(ℓₘₐₓ+1)² × N_modes)
const α_orig = Y_samp \ complex(mask_vals)         # least-squares → ψ₂ modes

# ─────────────────────────────────────────────────────────────────────────────
# Rotate via transform!
# ─────────────────────────────────────────────────────────────────────────────
# CW 90° about z sends +y hemisphere (A's "front") to +x hemisphere (B's "front").
# Panel 2's camera sits at −x looking toward +x, where the rotated F is face-on.

#const R_rot = rotor(-π/2, 0.0, 0.0, 1.0)   # CW 90° about z
const R_rot = exp(normalize(𝐢/3+𝐣+𝐤) * (π/3) / 2)   # rotate by π/3 about (y+z)/√2
const M = Quaternionic.to_rotation_matrix(R_rot) # active rotation matrix

# DataComponents(:ψ₂, :ψ₃, :ψ₄): ψ₂ = spin-0, ψ₃ = −1, ψ₄ = −2.
# ψ₂ requires ψ₃ and ψ₄ on ℐ⁺; set ψ₃ = ψ₄ = 0.
const dc = Scri.DataComponents(:ψ₂, :ψ₃, :ψ₄)

const N_time = 4                                    # ≥ 4 required by cubic spline
const t_arr = [0.0, 1.0, 2.0, 3.0]

data = zeros(ComplexF64, N_modes, N_time, 3)        # (modes, time, components)
for j ∈ 1:N_time
    data[:, j, 1] .= α_orig                         # ψ₂ = F pattern, constant in time
end

println("Running transform!…")
transform!(
    data,
    t_arr,
    QuatVec(0.0, 0.0, 0.0),            # zero boost
    R_rot,                             # frame rotation
    zeros(ComplexF64, 1),              # zero supertranslation
    dc,
)

const α_rot = data[:, 1, 1]            # rotated ψ₂ modes (constant → any slice)

# ─────────────────────────────────────────────────────────────────────────────
# Evaluate on a dense grid for plotting
# ─────────────────────────────────────────────────────────────────────────────

const Nθ, Nφ = 120, 240
const θg = range(0.0, Float64(π), Nθ)
const φg = range(0.0, 2Float64(π), Nφ)

println("Building dense evaluation grid…")
const Rs_mat = [R_f*from_spherical_coordinates(Float64(t), Float64(p)) for t ∈ θg, p ∈ φg]
const Y_dense = ₛ𝐘(0, ℓₘₐₓ, Float64, vec(Rs_mat))

const f_orig = clamp.(reshape(real.(Y_dense * α_orig), Nθ, Nφ), 0.0, 1.0)
const f_rot = clamp.(reshape(real.(Y_dense * α_rot), Nθ, Nφ), 0.0, 1.0)

# Sphere Cartesian coordinates (Nθ × Nφ)
const Xs = [sin(t) * cos(p) for t ∈ θg, p ∈ φg]
const Ys = [sin(t) * sin(p) for t ∈ θg, p ∈ φg]
const Zs = [cos(t) for t ∈ θg, p ∈ φg]

# ─────────────────────────────────────────────────────────────────────────────
# Face colors helper
# ─────────────────────────────────────────────────────────────────────────────
# matplotlib stores the (Nθ−1)×(Nφ−1) faces in row-major order:
#   face k (0-indexed) = (row=k÷(Nφ−1), col=k%(Nφ−1)).
# vec(fc') produces this ordering from a Julia (Nθ−1)×(Nφ−1) matrix.
# Passing 1-D Julia arrays to np.asarray avoids axis-reversal ambiguity.

function face_colors_np(np, f::Matrix{Float64})
    fc =
        (
            f[1:(end - 1), 1:(end - 1)] .+ f[2:end, 1:(end - 1)] .+ f[1:(end - 1), 2:end] .+
            f[2:end, 2:end]
        ) .* 0.25
    clamp!(fc, 0.0, 1.0)

    v = np.asarray(vec(fc'))        # 1-D numpy array, row-major ordering

    # Blend: grey (v=0) → crimson (v=1)
    return np.column_stack(
        pylist([
            0.62 + 0.24 * v,            # R: 0.62 → 0.86
            0.62 - 0.54 * v,            # G: 0.62 → 0.08
            0.62 - 0.38 * v,            # B: 0.62 → 0.24
            0.45 + 0.55 * v,            # alpha: 0.45 bg → 1.0 F
        ])
    )                              # shape (N, 4)
end

# ─────────────────────────────────────────────────────────────────────────────
# Axes arrows helper
# ─────────────────────────────────────────────────────────────────────────────
# mplot3d depth-sorts each artist by a single average z, so an arrow spanning
# the sphere lands arbitrarily in front of or behind it.  We disable depth
# sorting (computed_zorder=false) and layer manually:
#   Z_FAR    (1): outer arrow segments pointing away from the camera
#   Z_INSIDE (2): inner segments, origin → sphere surface
#   Z_SPHERE (3): the translucent sphere (tints everything below)
#   Z_NEAR   (4): outer segments pointing toward the camera, plus their labels
# Each arrow is split at the sphere surface; the outer piece carries the head.

const Z_FAR, Z_INSIDE, Z_SPHERE, Z_NEAR = 1, 2, 3, 4

# 3-D quiver heads are just two line "wings".  For proper filled arrowheads,
# subclass FancyArrowPatch to project 3-D endpoints to 2-D at draw time —
# the head is then rendered as a clean 2-D "-|>" from any viewing angle.
const arrow3d_ns = pydict()
pyexec(
    """
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return min(zs)
""",
    arrow3d_ns,
)
const Arrow3D = arrow3d_ns["Arrow3D"]

# Unit vector from scene origin toward the camera, given view_init angles
function camera_direction(elev, azim)
    return [cosd(elev) * cosd(azim), cosd(elev) * sind(azim), sind(elev)]
end

function draw_axes!(ax, M_axes, labels, cam; scale=1.55, alpha=1.0, linestyle="-")
    colors = ("tab:red", "tab:green", "tab:blue")
    for (i, (lbl, col)) ∈ enumerate(zip(labels, colors))
        u = M_axes[:, i]                       # unit axis direction
        v = u .* scale                          # arrow tip
        z_outer = sum(u .* cam) ≥ 0 ? Z_NEAR : Z_FAR

        # Inner segment: origin → sphere surface, always under the sphere
        ax.plot(
            [0.0, u[1]],
            [0.0, u[2]],
            [0.0, u[3]];
            color=col,
            alpha=alpha,
            linewidth=2.4,
            linestyle=linestyle,
            zorder=Z_INSIDE,
        )

        # Outer segment with a filled 2-D arrowhead: surface → tip
        arrow = Arrow3D(
            pylist([u[1], v[1]]),
            pylist([u[2], v[2]]),
            pylist([u[3], v[3]]);
            arrowstyle="-|>",
            mutation_scale=18,
            color=col,
            alpha=alpha,
            linewidth=2.4,
            linestyle=linestyle,
            shrinkA=0,
            shrinkB=0,
            zorder=z_outer,
        )
        ax.add_artist(arrow)

        ax.text(
            v[1] * 1.13,
            v[2] * 1.13,
            v[3] * 1.13,
            lbl;
            color=col,
            alpha=alpha,
            fontsize=16,
            fontweight="bold",
            ha="center",
            va="center",
            zorder=z_outer,
        )
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Figure
# ─────────────────────────────────────────────────────────────────────────────

println("Rendering figure…")
const np = pyimport("numpy")
fig = figure(; figsize=(14, 7))

const I3 = Matrix(1.0I, 3, 3)

panels = [
    (f_orig, 20, I3, ("x", "y", "z"), M, ("x′", "y′", "z′"), "Observer A — original frame"),
    (
        f_rot,
        20,
        I3,
        ("x′", "y′", "z′"),
        inv(M),
        ("x", "y", "z"),
        "Observer B — rotated frame",
    ),
]

for (k, (f_field, azim_val, M_solid, lbl_solid, M_ghost, lbl_ghost, title)) ∈
    enumerate(panels)
    ax = fig.add_subplot(120 + k; projection="3d", computed_zorder=false)

    # Solid translucent sphere with F pattern
    rgba_np = face_colors_np(np, f_field)
    surf = ax.plot_surface(
        Xs,
        Ys,
        Zs;
        rstride=1,
        cstride=1,
        shade=false,
        linewidth=0,
        antialiased=false,
        zorder=Z_SPHERE,
    )
    surf.set_facecolor(rgba_np)

    # Solid primary axes and ghost secondary axes, split at the sphere surface
    cam = camera_direction(22, azim_val)
    draw_axes!(ax, M_solid, lbl_solid, cam; alpha=1.0)
    draw_axes!(ax, M_ghost, lbl_ghost, cam; alpha=0.45, linestyle="--")

    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()
    ax.view_init(; elev=22, azim=azim_val)
    ax.set_title(title; fontsize=12, pad=10)
end

# fig.suptitle("Passive frame rotation on the celestial sphere  (ℓmax = $ℓₘₐₓ)",
#              fontsize=14, y=1.02)
tight_layout()

outfile = joinpath(@__DIR__, "rotation_sphere.png")
savefig(outfile; dpi=150, bbox_inches="tight")
println("Saved: ", outfile)
