# BMS Transformations

Though we go into detail about the *structure* of the BMS group
[elsewhere](@ref bms_group), we need to specify some basic conventions
about the coordinate transformation.  For anything not stated on this
page, read the theory behind it in the "BMS Transformations" section.

```@raw html
<!-- markdownlint-disable MD049 -->
```

We use the *__passive__* convention as found in most of the early
literature, whereby spacetime events remain fixed, but the coordinates
used to describe them change.  This contrasts with much of the more
recent literature, which tends to focus on charges and currents, for
which the active convention is more natural.

```@raw html
<!-- markdownlint-enable MD049 -->
```

We assume two sets (primed and unprimed) of asymptotic coordinates as
introduced on the [tetrad-conventions page](@ref "Tetrad").  The BMS
supertranslation affects only the time coordinate:

```math
u' = u - c_α α(θ,ϕ).
```

Here, ``c_α`` is a convention parameter, and ``α`` is a function only
of the angular coordinates.

The retarded time ``u`` is, of course, the coordinate on ``ℐ⁺``.
There is technically a *distinct* group of BMS transformations on
``ℐ⁻``, for which the advanced time ``v`` is the coordinate, which
transforms as

```math
v' = v - c_α α(θ,ϕ).
```

These two groups are frequently labeled BMS⁺ and BMS⁻, respectively.
Technically, it makes no sense to ask what the effect of a BMS⁺
transformation is on ``ℐ⁻``, or vice versa.  Nonetheless, it turns out
that the two groups are isomorphic — though in a slightly complicated
way.  So through this identification, we can treat the two groups as a
single BMS group.

## [Future and past null infinity](@id scri_pm_conventions)

The data that we manipulate *and* the various transformations exist on
a chosen piece of ``ℐ`` — the two pieces of which are future null
infinity ``ℐ⁺`` and past null infinity ``ℐ⁻``.  Knowing which piece we
are working on is important, and we will frequently need to retain a
sign to track that information.  By mild abuse of notation, we will
denote this variable simply as ``ℐ`` — or `ℐ` in the code.  The sign
is defined as

```math
ℐ = \begin{cases}
+1 & \text{future null infinity } ℐ⁺ \text{ (outgoing radiation)}, \\
-1 & \text{past null infinity } ℐ⁻ \text{ (incoming radiation)}.
\end{cases}
```

This information is crucial to understand — for example — the meaning
of the data described by a [`DataComponents`](@ref
Scri.DataComponents) object.  That is, it is used to specify which
actual points in the spacetime are represented by the data.  However,
there is a second, distinct, and fairly subtle piece of information
that is also important.

The subtlety is how the celestial sphere ``𝕊²`` is *labeled* in each
case.  At ``ℐ⁺`` we label a point by the direction in which outgoing
radiation *propagates*, so we use the future-pointing null ray ``𝐤 =
(1, k̂)``.  At ``ℐ⁻`` we work on the observer's past light cone — the
sphere of directions from which light *arrives* (the astronomer's
sky).  To see a source you look *opposite* to the light's direction of
travel, so the ``ℐ⁻`` ray is ``𝐤 = (1, -k̂)``.  That is, the labeling
of ``ℐ⁻`` is __antipodal__ to the ``ℐ⁺`` labeling.  Equivalently, a
single free null geodesic of Minkowski space joins a point ``k̂`` of
``ℐ⁺`` to the antipodal point ``-k̂`` of ``ℐ⁻``.  Identifying the two
spheres this way is the *antipodal matching* used, e.g., in the
analysis of gravitational scattering and soft theorems
[Strominger_2014, Strominger_2017](@cite).  Both choices are captured
by writing the null ray as ``𝐤 = (1, ℐ k̂)``.

A supertranslation must be represented as a function on the sphere, so
it is important to know which sphere is being used.  For example, it
would be reasonable to define a supertranslation acting on ``ℐ⁻`` via
its isomorphic BMS⁺ element.  Therefore each `BMS` object must track a
sign ``ℐ`` to indicate which piece of ``ℐ`` the supertranslation
should be interpreted as acting on.

So we have two independent places where a sign ``ℐ`` is used: in
`DataComponents` and in `BMS`.  They are conceptually distinct, but
can be related directly.
