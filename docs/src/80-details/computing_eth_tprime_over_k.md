# Computing ``ðt'/k``

The combination ``ðt'/k`` appears in the BMS transformation of the
Weyl components.  Thus, efficiently computing it is a key part of the
BMS pipeline.

We have the coordinate transformation

```math
\begin{gather}
u' = k(u + σᵅ α),
\\
\frac{1}{k} = γ(1 - σⁱv⃗⋅n̂),
\end{gather}
```

where ``σᵅ`` accounts for differences in the definition of the
supertranslation parameter ``α``, and ``σⁱ`` accounts for whether
we're dealing with ``ℐ⁺`` (``σⁱ = +1``) or ``ℐ⁻`` (``σⁱ = -1``).

Note that ``1/k`` has a simple form — and specifically, it decomposes
into a pure ``ℓ = 0`` plus a pure ``ℓ = 1`` function on the
sphere.  This allows us to compute ``ð(1/k)`` in closed form, without
needing to evaluate the synthesis sum at all.  The ``ð`` operator is a
derivation, meaning that it obeys the Leibniz rule, so as usual we
have ``ð(k) = -ð(1/k) k^2``.

```math
\begin{align}
\frac{ðt'}{k} &= -ð(1/k) k (u + σᵅ α) + σᵅ ðα \\
&= \frac{σⁱ ð(v⃗⋅n̂)}{1 - σⁱv⃗⋅n̂} (u + σᵅ α) + σᵅ ðα.
\end{align}
```

Since ``α`` and ``ðα`` may contain — in principle — arbitrarily high
``ℓ`` modes, it will be best to compute them via
`SphericalFunctions.jl`.  However, the ``v⃗⋅n̂`` and ``ð(v⃗⋅n̂)``
terms are simple enough that it will be more efficient to compute them
directly in closed form, using rotor components.  That is the
objective of what follows.

## Computing ``ð(v⃗⋅n̂)``

So we have simplified the problem to computing ``ð(v⃗⋅n̂)``.  The
product ``v⃗⋅n̂`` is an ``s=0``, ``ℓ = 1`` function.  The ``ð``
operator has a very simple action in mode space, so one approach is to
just evaluate the mode weights of ``v⃗⋅n̂`` and apply the known
mode-space formula for ``ð``.

Given the standard spherical-harmonic formulas

```math
\begin{aligned}
Y_{1,0}(θ, φ) &= \sqrt{\frac{3}{4π}}\,\cos θ,
\\
Y_{1,\pm 1}(θ, φ) &= \mp\sqrt{\frac{3}{8π}}\,\sin θ\,e^{\pm iφ},
\end{aligned}
```

we have

```math
\begin{aligned}
\cos θ &= \sqrt{\frac{4π}{3}}\,Y_{1,0},
\\
\sin θ \cos φ &= -\frac{1}{2}\sqrt{\frac{8π}{3}}\,\left(Y_{1,1} - Y_{1,-1}\right),
\\
\sin θ \sin φ &= \frac{i}{2}\sqrt{\frac{8π}{3}}\,\left(Y_{1,1} + Y_{1,-1}\right).
\end{aligned}
```

We can re-express the dot product as

```math
\begin{aligned}
v⃗⋅n̂
&= v_x\sin θ\cos φ + v_y\sin θ\sin φ + v_z\cos θ \\
&= -v_x\frac{1}{2}\sqrt{\frac{8π}{3}}\,\left(Y_{1,1} - Y_{1,-1}\right)
   + v_y\frac{i}{2}\sqrt{\frac{8π}{3}}\,\left(Y_{1,1} + Y_{1,-1}\right)
   + v_z \sqrt{\frac{4π}{3}} Y_{1,0} \\
\end{aligned}
```

Now, we also know that the ``ð`` operator acts on spin-weighted spherical harmonics by

```math
ð {}_{s}Y_{ℓ,m} = \sqrt{(ℓ-s)(ℓ+s+1)}\, {}_{s+1}Y_{ℓ,m}.
```

Of course, with ``s=0`` and ``ℓ=1``, the factor is just ``\sqrt{2}``, so we have

```math
\begin{aligned}
ð(v⃗⋅n̂)
&= -v_x\sqrt{\frac{4π}{3}}\,\left({}_{1}Y_{1,1} - {}_{1}Y_{1,-1}\right)
   + v_y i\sqrt{\frac{4π}{3}}\,\left({}_{1}Y_{1,1} + {}_{1}Y_{1,-1}\right)
   + v_z \sqrt{\frac{8π}{3}} {}_{1}Y_{1,0} \\
&= (-v_x+ v_y i)\sqrt{\frac{4π}{3}}\,{}_{1}Y_{1,1}
   + (v_x+ v_y i)\sqrt{\frac{4π}{3}}\, {}_{1}Y_{1,-1}
   + v_z \sqrt{\frac{8π}{3}} {}_{1}Y_{1,0}
\end{aligned}
```

Another approach is to use the original form of ``ð`` in terms of
partial derivatives with respect to spherical coordinates.  For
``s=0``, it is just

```math
ð f = -\left(\partial_θ f + \frac{i}{\sin θ}\,\partial_φ f\right).
```

Applying this to the explicit formula for ``v⃗⋅n̂`` given above, we have

```math
\begin{aligned}
ð(v⃗⋅n̂)
&= -v_x\cos θ\cos φ - v_y\cos θ\sin φ + v_z\sin θ
 - \frac{i}{\sin θ} \left(-v_x\sin θ\sin φ + v_y\sin θ\cos φ\right) \\
&= -v_x(\cos θ\cos φ - i\sin φ) - v_y(\cos θ\sin φ + i\cos φ) + v_z\sin θ
\end{aligned}
```

We have a pair of expressions for Wikipedia:

```math
\begin{aligned}
{}_1 Y_{10}(\theta,\phi)     &=  \sqrt{\frac{3}{8\pi}}\,\sin\theta, \\
{}_1 Y_{1\pm 1}(\theta,\phi) &= -\sqrt{\frac{3}{16\pi}}(1 \mp \cos\theta)\,e^{\pm i\phi}.
\end{aligned}
```

Plugging these in, we can verify that the two approaches give the same
result.

However, it will be better to compute this directly in terms of the
rotor components, which is more efficient and numerically well-behaved
(no need to evaluate any trigonometric functions, and no spurious
poles at the poles).  This is the objective of the next section.

## SWSHs and quaternions

The SWSH is related to Wigner's ``D`` matrix by

```math
{}_{s}Y_{\ell,m}(R) = (-1)^{s} \sqrt{ \frac{2\ell + 1} {4\pi}} \mathfrak{D}^{(\ell)}_{m, -s} (R),
```

and

```math
  \mathfrak{D}^{(\ell)}_{m',m}(R) = \sum_{\rho}
  \binom{\ell+m'} {\rho}\, \binom{\ell-m'} {\ell-\rho-m}\,
  (-1)^{\rho}\, R_s^{\ell+m'-\rho}\, \bar{R}_s^{\ell-\rho-m}\,
  R_a^{\rho-m'+m}\, \bar{R}_a^{\rho}\, \sqrt{ \frac{ (\ell+m)!\, (\ell-m)! } { (\ell+m')!\, (\ell-m')! } },
```

where ``R_s = R_1 + i R_z`` and ``R_a = R_y + i R_x``.  This gives a
direct way to evaluate the SWSHs in terms of the rotor components,
without needing to compute any trigonometric functions.

```math
\begin{gather}
R = \cos(θ/2) \cos(ϕ/2) - 𝐢 \sin(θ/2) \sin(ϕ/2) + 𝐣 \sin(θ/2) \cos(ϕ/2) + 𝐤 \cos(θ/2) \sin(ϕ/2),
\\
R_s = R_1 + i R_z = \cos(θ/2) \cos(ϕ/2) + i \cos(θ/2) \sin(ϕ/2) = \cos(θ/2) e^{iϕ/2},
\\
R_a = R_y + i R_x = \sin(θ/2) \cos(ϕ/2) + i (-\sin(θ/2) \sin(ϕ/2)) = \sin(θ/2) e^{-iϕ/2}.
\end{gather}
```

```math
\begin{aligned}
{}_{1}Y_{1,m}(R) &= -\sqrt{\frac{3} {4\pi}} \mathfrak{D}^{(1)}_{m, -1} (R) \\
&= -\sqrt{\frac{3} {4\pi}} \sum_{\rho}
  \binom{1+m} {\rho}\, \binom{1-m} {2-\rho}\,
  (-1)^{\rho}\, R_s^{1+m-\rho}\, \bar{R}_s^{2-\rho}\,
  R_a^{\rho-m-1}\, \bar{R}_a^{\rho}\, \sqrt{ \frac{ 2 } { (1+m)!\, (1-m)! } }
\end{aligned}
```

For ``m=0``, we must have ``\rho=1``:

```math
\begin{aligned}
{}_{1}Y_{1,0}(R)
&= \sqrt{\frac{3} {2\pi}} \bar{R}_s\, \bar{R}_a \\
&= \sqrt{\frac{3} {2\pi}} \sin(θ/2) \cos(θ/2) \\
&= \sqrt{\frac{3} {8\pi}} \sin(θ)
\end{aligned}
```

For ``m=1`` we have ``\rho=2``:

```math
\begin{aligned}
{}_{1}Y_{1,1}(R)
&= -\sqrt{\frac{3} {4\pi}} \bar{R}_a^{2} \\
&= -\sqrt{\frac{3} {4\pi}} \sin^2(θ/2) e^{iϕ} \\
&= -\sqrt{\frac{3} {16\pi}} (1 - \cos θ) e^{iϕ}
\end{aligned}
```

For ``m=-1``, we have ``\rho=0``:

```math
\begin{aligned}
{}_{1}Y_{1,-1}(R)
&= -\sqrt{\frac{3} {4\pi}} \bar{R}_s^{2} \\
&= -\sqrt{\frac{3} {4\pi}} \cos^2(θ/2) e^{-iϕ} \\
&= -\sqrt{\frac{3} {16\pi}} (1 + \cos θ) e^{-iϕ}
\end{aligned}
```

We have derived the same expressions in terms of spherical coordinates
as are found on Wikipedia, but we have also found their simple
expressions in terms of the rotor components, which is what we will
actually use.  We now express the full result

## Expressing ``ðt'/k`` with rotor components

A significant simplification will come from the fact that the rotation
of `v` by a rotor `R` is expressed as `R*v*conj(R) = R(v)` with this
result:

```julia
function (R::Rotor)(v::QuatVec)
    quatvec(SA[
        false,
        ((R[1]^2 + R[2]^2 - R[3]^2 - R[4]^2)*v[2]
            + (R[1]*R[3] + R[2]*R[4])*2v[4] + (R[2]*R[3] - R[1]*R[4])*2v[3]),
        ((R[1]^2 - R[2]^2 + R[3]^2 - R[4]^2)*v[3]
            + (R[2]*R[3] + R[1]*R[4])*2v[2] + (R[3]*R[4] - R[1]*R[2])*2v[4]),
        ((R[1]^2 + R[4]^2 - R[2]^2 - R[3]^2)*v[4]
            + (R[1]*R[2] + R[3]*R[4])*2v[3] + (R[2]*R[4] - R[1]*R[3])*2v[2])
    ])
end
```

Now, we expand the SWHSs in rotor components and find

```math
\begin{aligned}
ð(v⃗⋅n̂)
&= (-v_x + v_y i)\sqrt{\frac{4π}{3}}\,{}_{1}Y_{1,1}
   + (v_x + v_y i)\sqrt{\frac{4π}{3}}\, {}_{1}Y_{1,-1}
   + v_z \sqrt{\frac{8π}{3}} {}_{1}Y_{1,0} \\
&= (v_x - v_y i) \bar{R}_a^{2}
   - (v_x + v_y i) \bar{R}_s^{2}
   + 2v_z \bar{R}_s\, \bar{R}_a \\
&= (v_x - v_y i) (R_y - i R_x)^{2}
   - (v_x + v_y i) (R_1 - i R_z)^{2}
   + 2v_z (R_1 - i R_z)\, (R_y - i R_x) \\
&= (v_x - v_y i) (R_y^2 - R_x^2 - 2i R_x R_y) \\&\quad
   - (v_x + v_y i) (R_1^2 - R_z^2 - 2i R_1 R_z) \\&\quad
   + 2v_z (R_1 R_y - i R_y R_z - i R_1 R_x - R_x R_z) \\
&= v_x (- R_1^2 - R_x^2 + R_y^2 + R_z^2 + 2i (R_1 R_z - R_x R_y)) \\&\quad
   - i v_y (R_1^2 - R_x^2 + R_y^2 - R_z^2 - 2i (R_x R_y + R_1 R_z)) \\&\quad
   + 2v_z (R_1 R_y - i R_y R_z - i R_1 R_x - R_x R_z) \\
&= -\left( \tilde{R} \vec{v} R \right) \cdot \left( \hat{x} + i \hat{y} \right).
\end{aligned}
```

Similarly, we have

```math
\begin{aligned}
v⃗⋅n̂
&= (R_1^2 + R_z^2 - R_x^2 - R_y^2)v_z + 2(-R_1 R_x + R_y R_z)v_y + 2(R_x R_z + R_1 R_y)v_x \\
&= \left( \tilde{R} \vec{v} R \right) \cdot \hat{z}
\end{aligned}
```

These are very easy to calculate.  If we define ``Λ = R̃v⃗R`` then we
have

```math
\begin{align}
\frac{ðt'}{k}
&= \frac{σⁱ ð(v⃗⋅n̂)}{1 - σⁱv⃗⋅n̂} (u + σᵅ α) + σᵅ ðα \\
&= \frac{-σⁱ (Λ_x + i Λ_y)}{1 - σⁱ Λ_z} (u + σᵅ α) + σᵅ ðα.
\end{align}
```

The term proportional to ``u`` will vary between time steps, so we
factor out a term constant in time and one proportional to time:

```math
\frac{ðt'}{k}
= \left(\frac{ðt'}{k}\right)_0 + \left(\frac{ðt'}{k}\right)_1 u,
```

where

```math
\begin{aligned}
\left(\frac{ðt'}{k}\right)_0 &= σᵅ \left( \frac{Λ_x + i Λ_y}{Λ_z - σⁱ} α + ðα \right), \\
\left(\frac{ðt'}{k}\right)_1 &= \frac{Λ_x + i Λ_y}{Λ_z - σⁱ}.
\end{aligned}
```
