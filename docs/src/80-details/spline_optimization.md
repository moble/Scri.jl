# Cubic spline optimization

Naively, the single most expensive part of the numerical BMS
transformation is the interpolation of the data components from one
time grid to another.  However, there are significant opportunities
for optimization in our particular use case, so we implement a custom
cubic spline interpolation method that is tailored to our needs.  The
main idea is to interleave the construction of the spline coefficients
with the interpolation and transformation, so that we can avoid
redundant passes over the data and keep the working set small enough
to fit in L1 cache — or possibly even in registers — giving a
significant speedup over a more general-purpose implementation.

We use a natural cubic spline for this interpolation, which is a
piecewise cubic polynomial that is twice continuously differentiable
and has zero second derivatives at the endpoints.  The spline
coefficients are computed using the Thomas algorithm for solving the
tridiagonal system that arises from the natural cubic spline
conditions.

We suppose that we have a set of input time points ``tⱼ`` and
corresponding data values ``dₖⱼ`` for each component ``k`` at each
time point ``j``.  We want to evaluate the spline at a set of query
time points ``t'_{j'}`` that are contained within the range of the
``tⱼ`` (no extrapolation beyond machine-precision errors).  We abuse
notation to write the spline as a function ``dₖ``, which is given by
```math
\begin{aligned}
dₖ(t) &=
\frac{d̈ₖⱼ₊₁}{6} \frac{(t - tⱼ)³}{tⱼ₊₁-tⱼ}
+ \frac{d̈ₖⱼ}{6} \frac{(tⱼ₊₁ - t)³}{tⱼ₊₁-tⱼ}
\\ &\qquad
+ \frac{dₖⱼ₊₁ - d̈ₖⱼ₊₁ (tⱼ₊₁-tⱼ)²}{6} \frac{t - tⱼ}{tⱼ₊₁-tⱼ}
\\ &\qquad
+ \frac{dₖⱼ - d̈ₖⱼ (tⱼ₊₁-tⱼ)²}{6} \frac{tⱼ₊₁ - t}{tⱼ₊₁-tⱼ},
\end{aligned}
```
where the ``d̈ₖⱼ`` are found by imposing ``C^2`` continuity between
adjacent segments, giving rise to the tridiagonal system
```math
(tᵢ-tᵢ₋₁) d̈ₖⱼ₋₁ + 2 (tᵢ₊₁ - tᵢ₋₁) d̈ₖⱼ + (tᵢ₊₁ - tᵢ) d̈ₖⱼ₊₁
= 6 \left( \frac{dₖⱼ₊₁ - dₖⱼ}{tᵢ₊₁ - tᵢ} - \frac{dₖⱼ - dₖⱼ₋₁}{tᵢ - tᵢ₋₁} \right).
```
The latter may be written in matrix form as
```math
A\, d̈ₖ = rₖ,
```
where ``A`` is a tridiagonal matrix that depends only on the input
time points, and ``rₖ`` is a vector that depends on the input time
points and data values for component ``k``.  The goal is to solve for
the second derivatives ``d̈ₖ`` given ``rₖ``, so that we may evaluate
the spline at the query points.  To do so, we first decompose the
matrix ``A`` into its LU factors
```math
A = L U
=
\begin{pmatrix}
1  &    &   &  \\
l₁ & 1  &   &  \\
   & l₂ & 1 &  \\
   &    & \ddots & \ddots
\end{pmatrix}
\begin{pmatrix}
u₁ & s₁ &    &    &  \\
   & u₂ & s₂ &    &  \\
   &    & u₃ & s₃ &  \\
   &    &    & \ddots & \ddots
\end{pmatrix},
```
This decomposition into ``L`` and ``U`` only needs to be computed once
for a given set of input time points, and can be reused for all data
components.  As usual, the solution proceeds in two steps: first we
solve ``L zₖ = rₖ`` for the intermediate vector ``zₖ = Ud̈ₖ`` by
forward substitution with ``L``, and then we solve for ``d̈ₖ`` by
backward substitution with ``U``.  Because ``L`` and ``U`` are simple
bidiagonal matrices, these steps are very efficient and can be
implemented with minimal overhead.  But they must be done recursively:
```math
zₖⱼ = rₖⱼ - lⱼ zₖⱼ₋₁, \quad j ∈ 2:N-1,
```
and
```math
d̈ₖⱼ = \frac{zₖⱼ - sⱼ d̈ₖⱼ₊₁}{uⱼ}, \quad j ∈ N-1:-1:2
```

Now, given the data, we can compute the right-hand side vector ``rₖ``
for each component and proceed with the forward and backward
substitutions to find the second derivatives.  Finally, we can
evaluate the spline at the query points using the formula given above.
