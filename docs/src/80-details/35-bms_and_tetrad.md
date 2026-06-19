# BMS and the tetrad

The next question is how [the standard tetrad](@ref "The standard
tetrad") transforms under coordinate transformations — and
specifically BMS transformations.  One *very* important point is that
the tetrad is *defined* in terms of the coordinates.  In particular,
if we have an unprimed coordinate system and a primed coordinate
system, we actually have *two distinct* tetrads.  For example, we have
these tetrads that are regular near ``ℐ⁺``:

```math
\begin{aligned}
\tilde{l} &= -\frac{1}{\sqrt{2}} ∂_ω &\qquad\qquad
\tilde{l}' &= -\frac{1}{\sqrt{2}} ∂_{ω'} \\
\tilde{m} &= \frac{1}{\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_φ\right) &
\tilde{m}' &= \frac{1}{\sqrt{2}} \left(∂_{θ'} + \frac{i}{\sin θ'} ∂_{φ'}\right) \\
\tilde{m̄} &= \frac{1}{\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_φ\right) &
\tilde{m̄}' &= \frac{1}{\sqrt{2}} \left(∂_{θ'} - \frac{i}{\sin θ'} ∂_{φ'}\right) \\
\tilde{n} &= \sqrt{2} ∂ᵤ &
\tilde{n}' &= \sqrt{2} ∂_{u'}
\end{aligned}
```

It will be instructive to see how these two tetrads are related under
a supertranslation ``u' = u - εᵅα``.  On a surface of constant ``u'``,
we have ``u = u' + εᵅα``, but ``u'`` is constant with respect to other
primed coordinates.  Note that all other coordinates are unchanged.
The chain rule immediately shows that

```math
\begin{aligned}
∂_{u'} &= ∂_u, \\
∂_{ω'} &= ∂_ω, \\
∂_{θ'} &= ∂_θ + εᵅ \frac{∂α}{∂{θ'}} ∂_u, \\
∂_{φ'} &= ∂_φ + εᵅ \frac{∂α}{∂{φ'}} ∂_u.
\end{aligned}
```

Plugging these into the primed tetrad, we find

```math
\begin{aligned}
\tilde{l}' &= \tilde{l}, \\
\tilde{m}' &= \tilde{m} + εᵅ \frac{\tilde{m}(α)}{\sqrt{2}} \tilde{n}, \\
\tilde{m̄}' &= \tilde{m̄} + εᵅ \frac{\tilde{m̄}(α)}{\sqrt{2}} \tilde{n}, \\
\tilde{n}' &= \tilde{n} + εᵅ \left( \frac{∂α}{∂{θ'}} \tilde{m} + \frac{∂α}{∂{φ'}} \tilde{m̄} \right) + \frac{1}{2} εᵅ εᵝ \left( \frac{∂α}{∂{θ'}} \frac{∂β}{∂{θ'}} + \frac{∂α}{∂{φ'}} \frac{∂β}{∂{φ'}} \right) \tilde{l}.
\end{aligned}
```
