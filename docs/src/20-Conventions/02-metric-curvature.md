# Metric and Curvature

The default signature used in this package is ``{-}{+}{+}{+}``.  A
`Conventions` parameter with ``c_s = -1`` flips the signature to
``{+}{-}{-}{-}``.

Our definitions of the curvature quantities generally follow the
conventions of Misner, Thorne, and Wheeler
[MisnerThorneWheeler_1973](@cite).  Their Eq. (14.36) defines the
Christoffel symbols by

```math
{\Gamma^a}_{bc} = \frac{1}{2} g^{ad} \bigl(
  \partial_b g_{cd} + \partial_c g_{bd} - \partial_d g_{bc}
\bigr).
```

The Riemann tensor is the only one of these that will frequently be
seen with a different sign.  We default to the expression given in
MTW's Eq. (8.44), though we include the `Conventions` parameter
``c_R`` in the definition:

```math
{R^a}_{bcd} = c_R \left(
  \partial_c {\Gamma^a}_{bd}
  - \partial_d {\Gamma^a}_{bc}
  + {\Gamma^a}_{ce} {\Gamma^e}_{bd}
  - {\Gamma^a}_{de} {\Gamma^e}_{bc}
\right).
```

It's worth pausing a moment to consider what another author's choices
might lead to.  Suppose author ``A`` uses the opposite signature,
``c_s = -1``.  We'll use prefixed superscripts ``[A]`` on all
quantities used by the author.  Then their metric, ``{}^{[A]}g_{ab}``
is related to our metric — which we always write without the
superscripts as ``g_{ab}`` — according to

```math
\begin{gathered}
{}^{[A]}g_{ab} = c_s g_{ab} \\
{}^{[A]}g^{ab} = c_s g^{ab},
\end{gathered}
```

(remembering that ``c_s ∈ \{-1, 1\}``).  We assume that there is
universal agreement on the *formula* used to define the Christoffel
symbols, so

```math
{}^{[A]}{\Gamma^a}_{bc} = \frac{1}{2} {}^{[A]}g^{ad} \bigl(
  \partial_b {}^{[A]}g_{cd} + \partial_c {}^{[A]}g_{bd} - \partial_d {}^{[A]}g_{bc}
\bigr),
```

but due to cancellation of two factors of ``c_s`` when we substitute
the expressions for the metrics, we have ``{}^{[A]}{\Gamma^a}_{bc} =
{\Gamma^a}_{bc}``.  Now suppose this author defines

```math
{}^{[A]}{R^a}_{bcd} = - \left(
  \partial_c {}^{[A]}{\Gamma^a}_{bd}
  - \partial_d {}^{[A]}{\Gamma^a}_{bc}
  + {}^{[A]}{\Gamma^a}_{ce} {}^{[A]}{\Gamma^e}_{bd}
  - {}^{[A]}{\Gamma^a}_{de} {}^{[A]}{\Gamma^e}_{bc}
\right).
```

This author uses ``c_R = -1``, and we have ``{}^{[A]}{R^a}_{bcd} = -
{R^a}_{bcd}`` — or generally ``{}^{[A]}{R^a}_{bcd} = c_R
{R^a}_{bcd}``.  But now consider the *lowered-index* Riemann tensor:

```math
{}^{[A]}R_{abcd}
= {}^{[A]}g_{ae} {}^{[A]}{R^e}_{bcd}
= c_s c_R g_{ae} {R^e}_{bcd}
= c_s c_R {R}_{abcd}.
```

The first equation above is what the author would write down in their
own conventions (though they wouldn't have written the superscripts,
obviously).

The Ricci tensor is almost universally defined as

```math
R_{ab} = {R^c}_{acb},
```

as in Eq. (8.47) of MTW.  So we have ``{}^{[A]}R_{ab} = c_R R_{ab}``.
The Ricci scalar is then

```math
{}^{[A]}R = {}^{[A]}g^{ab} {}^{[A]}R_{ab} = c_s c_R R.
```

And the Weyl tensor is universally defined to track the sign of the
Riemann tensor, as in Eq. (13.50) of MTW:

```math
C_{abcd} =
R_{abcd}
- \frac{1}{2} (g_{ac} R_{bd} - g_{ad} R_{bc} + g_{bd} R_{ac} - g_{bc} R_{ad})
+ \frac{1}{6} R (g_{ac} g_{bd} - g_{ad} g_{bc}).
```

We can see that

```math
{}^{[A]}C_{abcd} = c_s c_R C_{abcd},
```

and each term in the definition is consistent individually with this
factor.
