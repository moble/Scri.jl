# Augmented direct SSHT

The spin-spherical-harmonic transform (SSHT) is a key component of the
BMS pipeline, and its performance is critical to the overall
efficiency of the transformation.  A great deal of literature exists
on fast algorithms for computing the SSHT, primarily because the
Cosmic-Microwave-Background (CMB) community has a strong interest in
analyzing large datasets giving polarization information on the entire
sky.  However, it is important to note that the CMB community uses
very high angular resolution — values of `ℓₘₐₓ` into the thousands.
The gravitational-wave modeling community, on the other hand,
typically uses much lower angular resolution, with `ℓₘₐₓ` values
typically not beyond 8.  When performing BMS transformations of such
data, we need to account for the fact that there will be mode mixing
up to significantly higher values of `ℓ` due to the supertranslation
and boost effects, but still only typically up to `ℓₘₐₓ` of several
dozen as long as the boost is not too extreme.  This opens up the
possibility of using a much more direct algorithm for computing the
SSHT, which has the advantage of being simpler to implement, and
enabling *in-place* transforms that avoid the need for large temporary
arrays.

However, there is an important complication when using the "direct"
SSHT in BMS transformations that involve data components with multiple
spin weights.  The natural set of pixels on which to evaluate a field
*depends on the spin weight*.  Specifically, the number of pixels
would naturally be the same as the number of (nonzero) modes for that
spin weight up to the given ``ℓₘₐₓ``, which is ``Nᵖ = (ℓₘₐₓ + 1)² -
s²``; higher-spin fields would have fewer pixels than lower-spin
fields.  But if the data components need to be mixed together, we need
to evaluate them on the same set of pixels.

## Direct SSHT
For a given spin weight ``s``, we suppose that we have a set of mode
weights ``𝐦``.  Then it is easy to evaluate those weights on a set of
pixels by matrix multiplication:
```math
{}_{s}𝐘 · 𝐦 = 𝐩,
```
where ``{}_{s}𝐘_{ij} = {}{}_{s}Y_{ℓ,m}(R_j)`` is the synthesis matrix of the
SSHT, the index ``(ℓ,m)`` corresponds to the ``i``-th entry of the
vector of mode weights, and ``R_j`` is the position of the ``j``-th
pixel.  The pixel values are then modified as needed, and the analysis step is given by the inverse operation:
```math
{}_{s}𝐘^{-1} · 𝐩 = 𝐦.
```
Obviously, rather than computing the inverse of the synthesis matrix,
we can perform the analysis step by solving the linear system as usual
using an LU or QR factorization that is precomputed and stored in the
`{}_{s}𝐘` object.

## Augmenting the matrix

When all spin-weight components share a common pixel grid and there
are different spin weights, the synthesis matrix ``{}_{s}𝐘`` has ``Nᵖ``
rows and ``Nᵐ`` columns, with ``Nᵖ > Nᵐ``: there are more pixels than
modes.  The analysis step is then an overdetermined linear system, and
a naive least-squares approach is both more expensive and less clean
than what the structure of the problem already provides.

### Column space and null space

The ``Nᵐ`` columns of ``{}_{s}𝐘`` span an ``Nᵐ``-dimensional subspace of
``ℂ^{Nᵖ}`` — the *column space* of ``{}_{s}𝐘``.  Any physically meaningful
pixel vector must lie exactly there: the synthesis equation ``𝐩 =
{{}_{s}𝐘} \cdot 𝐦`` admits no solution unless ``𝐩`` is already in
col(``{}_{s}𝐘``).

The orthogonal complement of col(``{}_{s}𝐘``) has dimension ``Nᵖ - Nᵐ``
and is called the *null space of ``{}_{s}𝐘ᴴ``*: the set of pixel vectors
``𝐯 ∈ ℂ^{Nᵖ}`` satisfying ``{}_{s}𝐘ᴴ\, 𝐯 = 0``.  Such a vector is
orthogonal to every mode shape on the pixel grid; no linear
combination of harmonics up to ``ℓₘₐₓ`` can produce any component in
this direction.  This is not a numerical artifact — it is the
geometric consequence of evaluating a finite mode basis on a pixel
grid that is larger than strictly necessary.

### QR decomposition identifies both subspaces at once

The QR decomposition of ``{}_{s}𝐘`` produces a unitary ``Q ∈ ℂ^{Nᵖ×Nᵖ}``
such that the first ``Nᵐ`` columns of ``Q`` form an orthonormal basis
for col(``{}_{s}𝐘``).  Because ``Q`` is unitary, its remaining ``Nᵖ - Nᵐ``
columns are automatically orthogonal to the first ``Nᵐ`` and therefore
span *exactly* the null space of ``{}_{s}𝐘ᴴ``.

It so happens that `LinearAlgebra.qr` in Julia stores the ``Q`` factor
as Householder reflections.  To obtain the orthogonal columns, we just
need to apply the reflections to the identity matrix, after which we
obtain ``Q_⊥`` by slicing.

### Augmenting to a square, prefactorable system

Prepending ``Q_⊥`` to ``{}_{s}𝐘`` yields a square ``Nᵖ × Nᵖ`` matrix whose
columns span all of ``ℂ^{Nᵖ}``.  Every pixel vector decomposes
uniquely into a component in col(``{}_{s}𝐘``) and a component in
col(``Q_⊥``), so the augmented system
```math
\bigl[Q_⊥ \big| {{}_{s}𝐘}\bigr]
\begin{bmatrix} ξ \\ 𝐦 \end{bmatrix}
= 𝐩
```
always has a unique solution.  The mode weights ``𝐦`` appear in the
lower ``Nᵐ`` entries of the solution vector, and the auxiliary
variable ``ξ`` absorbs whatever null-space content is present in
``𝐩``.

The length of ``ξ`` is ``Nᵖ - Nᵐ``.  By putting it at the front, we
happen to be adding storage corresponding to modes for which ``ℓ <
|s|``.  These are trivially zero analytically, so we frequently omit
them.  But when processing multiple data components together, it is
more convenient to keep them around, so that we have a uniform
three-dimensional array.  By performing the linear solve of this
augmented system, we can use that extra storage to make the forward
and reverse SSHT transformations both in-place and allocation-free.

Analytically, assuming no numerical errors and infinite bandwidth, the
null-space component should be absent and ``ξ ≈ 0``.  A large ``ξ``
therefore signals that the data contain power that the mode basis
cannot represent: a truncation artifact, or a sufficiently extreme
transformation that has pushed power above ``ℓₘₐₓ``.  Checking the
power in the ``ξ`` component therefore provides a consistency check on
the transformation at small extra cost.

Because ``Q_⊥`` has orthonormal columns, the augmented matrix inherits
the conditioning of ``{}_{s}𝐘`` itself.  Its LU factorization is
computed once during setup; every subsequent analysis step reduces to
a single in-place triangular solve, with no temporary allocations and
no regularization choices.  When ``Nᵖ = Nᵐ`` — that is, when the spin
weight is low enough that the mode basis already fills the pixel grid
— no augmentation is necessary and a direct LU of ``{}_{s}𝐘``
suffices.
