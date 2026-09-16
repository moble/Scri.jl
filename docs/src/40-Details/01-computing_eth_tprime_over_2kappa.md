# [Computing ``ðt'/2κ``](@id computing_eth_tprime_over_2kappa)

The combination ``ðt'/2κ`` (where ``t`` can represent either ``u`` or
``v``) appears in the BMS transformation of the Weyl components.
Thus, computing it accurately and efficiently is a key part of the BMS
pipeline.

We have the coordinate transformation

```math
\begin{gather}
t' = κ(t - c_α α),
\\
\frac{1}{κ} = γ(1 - ℐv⃗⋅k̂),
\end{gather}
```

where ``c_α`` accounts for differences in the definition of the
supertranslation parameter ``α``, and ``ℐ`` accounts for whether
we're dealing with ``ℐ⁺`` (``ℐ = +1``) or ``ℐ⁻`` (``ℐ = -1``).

Note that ``1/κ`` has a simple form — and specifically, it decomposes
into a pure ``ℓ = 0`` plus a pure ``ℓ = 1`` function on the
sphere.  This allows us to compute ``ð(1/κ)`` in closed form, without
needing to evaluate the synthesis sum at all.  The ``ð`` operator is a
derivation, meaning that it obeys the Leibniz rule, so as usual we
have ``ð(κ) = -ð(1/κ) κ^2``.

```math
\begin{align}
\frac{ðt'}{2κ}
&= \frac{ð(κ)}{2κ} (t - c_α α) - \frac{c_α ðα}{2} \\
&= \frac{-ð(1/κ)}{2}κ (t - c_α α) - \frac{c_α ðα}{2} \\
&= \frac{ℐ ð(v⃗⋅k̂)}{2(1 - ℐv⃗⋅k̂)} (t - c_α α) - \frac{c_α ðα}{2} \\
&= \frac{-ð(v⃗⋅k̂)}{2(v⃗⋅k̂-ℐ)} (t - c_α α) - \frac{c_α ðα}{2}.
\end{align}
```

Since ``α`` and ``ðα`` may contain — in principle — arbitrarily high
``ℓ`` modes, it will be best to compute them via
`SphericalFunctions.jl`.  However, the ``v⃗⋅k̂`` and ``ð(v⃗⋅k̂)``
terms are simple enough that it will be more efficient to compute them
directly in closed form, using rotor components.  That is the
objective of what follows.

## ``v⃗⋅k̂`` as a function on ``\mathrm{Spin}(3)``

As mentioned when we [introduced the ``ð`` operator](@ref
the-operator-eth), ``ð`` is best understood as a derivative operator
on the group ``\mathrm{Spin}(3)``.  It is not immediately obvious that
``v⃗⋅k̂`` needs to be expressed as a function on ``\mathrm{Spin}(3)``,
but in order to obtain an expression for ``ðt'/2κ`` that is a function
on ``\mathrm{Spin}(3)``, we need to express all of the terms in that
form.  The solution is simple:

```math
v⃗⋅k̂ : Q ↦ v⃗ ⋅ (Q ẑ Q̄) = (Q̄ v⃗ Q) ⋅ ẑ.
```

That is, ``Q`` is interpreted as a rotation that takes the ``ẑ``
basis vector to the direction ``k̂``, and the inner product is taken
with ``v⃗``.  Note the equivalent form in the final expression, which
will be make the final result more efficient to calculate.  In fact,
we will define

```math
\vec{λ} = Q̄ v⃗ Q,
```

which can be calculated once, very efficiently.  Then we just need to
take the ``ẑ`` component of ``\vec{λ}`` to get ``v⃗⋅k̂ = λᶻ``, and we
will see that its other components are exactly what we need to compute
``ð(v⃗⋅k̂)``.

## Computing ``ð(v⃗⋅k̂)``

We can now immediately apply the definition of ``ð`` to compute
``ð(v⃗⋅k̂)``.  We have

```math
\begin{aligned}
ð(v⃗⋅k̂)
&= c_ð R_{x̂+iŷ} v⃗⋅k̂ \\
&= -c_ð i \left.\frac{d}{dϵ}\right|_{ϵ=0} \left\{
  v⃗ ⋅ \left[Q e^{-ϵ(x̂+iŷ)/2} ẑ e^{ϵ(x̂+iŷ)/2} Q̄\right] \right\} \\
&= -c_ð i v⃗ ⋅ \left\{Q \left[ \left.\frac{d}{dϵ}\right|_{ϵ=0} \left(
  e^{-ϵ(x̂+iŷ)/2} ẑ e^{ϵ(x̂+iŷ)/2}\right) \right] Q̄ \right\} \\
&= -c_ð i v⃗ ⋅ \left\{Q \left[ -\frac{x̂+iŷ}{2} ẑ + ẑ \frac{x̂+iŷ}{2} \right] Q̄ \right\} \\
&= -c_ð v⃗ ⋅ \left[Q \left( iŷ+x̂ \right) Q̄\right] \\
&= -c_ð \left( Q̄ v⃗ Q \right) ⋅ \left( x̂+iŷ \right) \\
&= -c_ð \left( λˣ + i λʸ \right).
\end{aligned}
```

## Putting it together

The final result is simple to calculate:

```math
\begin{aligned}
\frac{ðt'}{2κ}
&= \frac{-ð(v⃗⋅k̂)}{2(v⃗⋅k̂-ℐ)} (t - c_α α) - \frac{c_α ðα}{2} \\
&= c_ð \frac{λˣ + i λʸ}{2(λᶻ-ℐ)} (t - c_α α) - \frac{c_α ðα}{2}.
\end{aligned}
```

The term proportional to ``t`` will vary between time steps, so we
factor out a term constant in time and one proportional to time:

```math
\frac{ðt'}{2κ}
= \left(\frac{ðt'}{2κ}\right)_0 + \left(\frac{ðt'}{2κ}\right)_1 t,
```

where

```math
\begin{aligned}
\left(\frac{ðt'}{2κ}\right)_0 &= -c_α \left( c_ð \frac{λˣ + i λʸ}{2(λᶻ-ℐ)} α + \frac{ðα}{2} \right), \\
\left(\frac{ðt'}{2κ}\right)_1 &= c_ð \frac{λˣ + i λʸ}{2(λᶻ-ℐ)}.
\end{aligned}
```
