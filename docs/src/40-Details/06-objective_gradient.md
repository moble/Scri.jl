```@meta
CurrentModule = Scri
```

# [The gradient of the objective](@id objective_gradient)

[`objective_and_gradient!`](@ref) returns the objective of
[`transform_objective`](@ref) together with its gradient with respect
to every parameter of the BMS transformation.  It costs about two
evaluations of the objective, however many parameters there are,
whereas forward-mode automatic differentiation costs about one
evaluation per parameter.  This page sets up the notation, writes out
the objective, shows the structure of the computation in a diagram,
and then derives each piece.

## Notation

### Indices

Every quantity is written with all of its indices, and every sum is
written out; repeated indices are never summed implicitly.

| Index | Runs over | Number |
| :--- | :--- | :--- |
| ``i`` | directions ``R'_i`` of the grid in the new frame | ``N_\mathrm{dir}`` |
| ``j`` | output times ``t'_j`` | ``N_{t'}`` |
| ``a, b, b'`` | input times ``t_a`` | ``N_t`` |
| ``k, k'`` | components of the data, such as ``h``, ``ψ_4``, and ``ψ_3``, of spin weight ``s_k`` | ``N_\mathrm{comp}`` |
| ``ℓ, m`` | modes | ``N_\mathrm{modes}`` |
| ``p`` | parameters ``θ_p`` | ``P`` |
| ``r`` | the eight real numbers that describe the geometry of one direction | ``8`` |

### The parameters

The parameter vector is

```math
θ = \left(𝐮,\ 𝐛,\ a_1, …, a_{(ℓ_α+1)^2}\right),
```

made of a rotation vector ``𝐮``, a rapidity vector ``𝐛``, and the real
degrees of freedom ``a`` of the supertranslation's modes up to
``ℓ_α``.  [`bms_from_parameters`](@ref) maps ``θ`` to the boost
velocity, the frame rotation, and the supertranslation:

```math
v⃗ = \tanh|𝐛|\, \frac{𝐛}{|𝐛|},
\qquad
R = \exp(𝐮/2),
\qquad
α_{ℓ,m} \text{ from the } a,
```

where the docstring of `bms_from_parameters` gives the layout of the
``a``.

### The data, the target, and the weights

The input is the modes ``d_{ℓ,m;ak}`` of component ``k`` at input time
``t_a``.  The [`ObjectivePlan`](@ref) holds everything else, none of
which depends on ``θ``:

- the grid's directions ``R'_i``, with quadrature weights ``w_i``;
- the output times ``t'_j``, with weights ``v_j``;
- a weight ``c_k`` for each component;
- the target ``d_{\mathrm{t};ijk}``: the value of component ``k`` of
  the target waveform in direction ``R'_i`` at time ``t'_j``, computed
  once from the target's modes;
- the normalization ``𝒩``, which is either ``1`` or the target's norm.

### The geometry of each direction

For each direction ``R'_i`` of the grid, the transformation determines

- the rotor ``R_i`` of the corresponding direction in the original
  frame;
- the conformal factor ``κ_i``;
- the supertranslation and its derivatives there, ``α_i``, ``ðα_i``,
  and ``ð^2α_i``;
- the null-rotation parameter, which is affine in the original
  retarded time ``u``:

  ```math
  \left(\frac{ðu'}{2κ}\right)_i = ν^{(0)}_i + u\, ν^{(1)}_i
  ```

  (see [Computing ``ðt'/2κ``](@ref computing_eth_tprime_over_2kappa);
  ``ðα_i`` enters through ``ν^{(0)}_i``).

Everything after this is computed from two sets of numbers per
direction.  The first is the harmonics at the original direction, for
each spin weight ``s`` among the components:

```math
{}_sY_{i;ℓ,m} = {}_sY_{ℓ,m}(R_i).
```

The second is eight real numbers:

```math
\left(σ_{i1}, …, σ_{i8}\right)
= \left( κ_i^{-1},\ α_i,\
  \operatorname{Re} ν^{(0)}_i,\ \operatorname{Im} ν^{(0)}_i,\
  \operatorname{Re} ν^{(1)}_i,\ \operatorname{Im} ν^{(1)}_i,\
  \operatorname{Re} ð^2α_i,\ \operatorname{Im} ð^2α_i \right).
```

In direction ``i``, the output time ``t'_j`` corresponds to the
original time

```math
τ_{ij} = κ_i^{-1}\, t'_j + α_i.
```

### The objective

The objective is

```math
J = \frac{1}{𝒩} \sum_i J_i,
\qquad
J_i = w_i \sum_j v_j \sum_k c_k
      \left| d'_{ijk} - d_{\mathrm{t};ijk} \right|^2,
```

where ``d'_{ijk}`` is the transformed waveform in direction ``R'_i`` at
output time ``t'_j``.  It is built in three steps.

**Step 1, the mode sum,** gives the value of component ``k`` in the
original direction ``R_i`` at each input time:

```math
x_{iak} = \sum_{ℓ,m} d_{ℓ,m;ak}\, {}_{s_k}Y_{i;ℓ,m}.
```

**Step 2, interpolation,** evaluates the natural cubic spline through
those samples at the times ``τ_{ij}``:

```math
y_{ijk} = \sum_a S_{ija}\, x_{iak}.
```

The weights ``S_{ija}`` are real.  They depend on ``τ_{ij}`` and the
input times, but not on the data.

**Step 3, the transformation laws,** mix the components:

```math
d'_{ijk} = \sum_{k'} M_{ijkk'}\, y_{ijk'} + b_{ijk}.
```

The coefficients ``M_{ijkk'}`` hold the conformal factors and the
mixing by the null rotation, so they depend on ``κ_i^{-1}`` and on
``ν^{(0)}_i + τ_{ij}\, ν^{(1)}_i``.  On ``ℐ^+``, for example, the
transformed ``ψ_3`` is ``κ^{-3}`` times ``ψ_3`` plus a multiple of
``ψ_4`` proportional to the null-rotation parameter.  The offsets
``b_{ijk}`` are nonzero only for the shear and the strain, and are
proportional to ``κ_i^{-1}`` and to ``ð^2α_i`` or its conjugate.
[`mix_components!`](@ref) applies both, and
[`mixing_matrix!`](@ref Scri.mixing_matrix!) extracts ``M_{ijkk'}``
from it.

All together,

```math
d'_{ijk} = \sum_{k'} M_{ijkk'} \sum_a S_{ija} \sum_{ℓ,m}
           d_{ℓ,m;ak'}\, {}_{s_{k'}}Y_{i;ℓ,m} + b_{ijk}.
```

Step 1 costs ``N_t\, N_\mathrm{comp}\, N_\mathrm{modes}`` complex
multiply-adds per direction, which dominates everything else except
at the smallest band limits.

### Derivatives with respect to complex quantities

``J`` is real, but most of the quantities it depends on are complex.
For a complex quantity ``z``, treat ``z`` and ``z^*`` as independent
variables, so that

```math
\frac{∂J}{∂z^*}
= \frac{1}{2}\left( \frac{∂J}{∂(\operatorname{Re} z)}
  + \mathrm{i}\, \frac{∂J}{∂(\operatorname{Im} z)} \right),
\qquad
δJ = 2 \operatorname{Re} \sum_z \left(\frac{∂J}{∂z^*}\right)^* δz.
```

Every backward step below is an instance of one rule.  If
``z_k = \sum_{k'} A_{kk'}\, x_{k'}`` with coefficients that do not
depend on ``x``, then ``z^*_k = \sum_{k'} A^*_{kk'}\, x^*_{k'}``, so

```math
\frac{∂J}{∂x^*_{k'}} = \sum_k A^*_{kk'}\, \frac{∂J}{∂z^*_k}.
```

When the ``A_{kk'}`` are real, the conjugate does nothing.  The code
stores ``2\, ∂J/∂z^*``, which is the gradient with respect to the real
and imaginary parts packed into one complex number; that only moves
the factors of 2.  In the code, a name such as `∇x` holds
``2\, ∂J/∂x^*``, or ``∂J/∂x`` when ``x`` is real.

## The structure of the computation

The objective is a sum over directions.  In each direction, the
expensive work is the chain of steps 1 to 3, which is linear in the
input data, and the parameters enter it only through ``{}_sY_{i;ℓ,m}`` and
``σ_{ir}``.  So the chain is differentiated in reverse, applying the
rule above to each step, which costs one more pass of the same size.
The geometry, which is cheap, is differentiated forward with
ForwardDiff's dual numbers.  The diagram shows the flow for one
direction.

```@raw html
<div class="gradient-flow-diagram" style="text-align: center; margin: 1.5em 0;">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1070 790" width="100%" style="max-width: 900px;" role="img"
     aria-label="The flow of the objective and its gradient through each direction of the grid">
  <style>
    .grad-flow text { fill: currentColor; text-anchor: middle; font-size: 17px; white-space: pre; }
    .grad-flow text.cap { font-size: 14.5px; }
    .grad-flow text.lab { font-size: 16px; }
    .grad-flow .mi { font-family: KaTeX_Math, "Times New Roman", serif; font-style: italic; }
    .grad-flow .mu { font-family: KaTeX_Main, "Times New Roman", serif; }
    .grad-flow .pr { font-family: "Times New Roman", "Liberation Serif", serif; }
    .grad-flow .box { stroke: currentColor; stroke-width: 1.2; fill: none; }
    .grad-flow .geom { fill: rgba(56, 132, 255, 0.17); }
    .grad-flow .rev { fill: rgba(255, 140, 0, 0.22); }
    .grad-flow .plain { fill: rgba(128, 128, 128, 0.10); }
    .grad-flow .heavy { stroke-width: 3; }
    .grad-flow .edge { stroke: currentColor; stroke-width: 1.3; fill: none; marker-end: url(#grad-flow-arrow); }
  </style>
  <defs>
    <marker id="grad-flow-arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7"
            orient="auto-start-reverse">
      <path d="M0,0 L10,5 L0,10 z" fill="currentColor"/>
    </marker>
  </defs>
  <g class="grad-flow">
    <path class="edge " d="M 52,54 L 70,54"/>
  <path class="edge " d="M 166,54 L 195,54"/>
  <path class="edge " d="M 222,96 L 222,139"/>
  <path class="edge " d="M 567,96 L 567,139"/>
  <path class="edge " d="M 222,197 L 222,273"/>
  <path class="edge " d="M 527,197 L 476,273"/>
  <path class="edge " d="M 607,197 L 658,273"/>
  <path class="edge " d="M 94,318 L 117,318"/>
  <path class="edge " d="M 327,318 L 347,318"/>
  <path class="edge " d="M 557,318 L 577,318"/>
  <path class="edge " d="M 787,318 L 805,318"/>
  <path class="edge " d="M 1019,318 L 1046,318"/>
  <text x="1052" y="323" class="lab" style="text-anchor: start"><tspan
      class="mi">J</tspan></text>
  <path class="edge " d="M 912,363 L 912,431"/>
  <text x="921" y="402.0" class="lab" style="text-anchor: start"><tspan
      class="mi">d</tspan><tspan
      class="pr">′</tspan><tspan
      class="mi" dy="5" font-size="72%">ijk</tspan><tspan
      class="mu" dy="-5"> − </tspan><tspan
      class="mi">d</tspan><tspan
      class="mu" dy="5" font-size="72%">t</tspan><tspan
      class="mu" font-size="72%">;</tspan><tspan
      class="mi" font-size="72%">ijk</tspan></text>
  <path class="edge rev" d="M 805,476 L 787,476"/>
  <path class="edge rev" d="M 577,476 L 557,476"/>
  <path class="edge rev" d="M 347,476 L 327,476"/>
  <path class="edge " d="M 892,521 L 862,562"/>
  <path class="edge " d="M 222,521 L 222,673"/>
  <path class="edge " d="M 800,638 L 800,673"/>
    <g>
    <rect class="box io" x="8" y="32" width="44" height="44" rx="7"/>
    <text x="30" y="59"><tspan
      class="mi">θ</tspan></text>
  </g>
  <g>
    <rect class="box geom" x="70" y="32" width="96" height="44" rx="7"/>
    <text x="118" y="59"><tspan
      class="mi">v</tspan><tspan
      class="mu">⃗, </tspan><tspan
      class="mi">R</tspan><tspan
      class="mu">, </tspan><tspan
      class="mi">α</tspan></text>
  </g>
  <g>
    <rect class="box geom" x="195" y="12" width="830" height="84" rx="7"/>
    <text x="610" y="38" class="cap">for each direction of the new frame’s grid: the original direction, the conformal factor,</text>
    <text x="610" y="59" class="cap">the supertranslation and its derivatives there, and the null-rotation parameter</text>
    <text x="610" y="80"><tspan
      class="mi">R</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5">,   </tspan><tspan
      class="mi">κ</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5">,   </tspan><tspan
      class="mi">α</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5">,   ð</tspan><tspan
      class="mi">α</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5">,   ð</tspan><tspan
      class="mu" dy="-7" font-size="72%">2</tspan><tspan
      class="mi" dy="7">α</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5">,   </tspan><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-12" dx="-4.9" font-size="72%">(0)</tspan><tspan
      class="mu" dy="7"> + </tspan><tspan
      class="mi">u</tspan><tspan
      class="mu"> </tspan><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-12" dx="-4.9" font-size="72%">(1)</tspan></text>
  </g>
  <g>
    <rect class="box geom" x="117" y="139" width="210" height="58" rx="7"/>
    <text x="222" y="162" class="cap">harmonics</text>
    <text x="222" y="184"><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="-5">Y</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" font-size="72%">;</tspan><tspan
      class="mi" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" dy="-5">,   ∂</tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="-5">Y</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" font-size="72%">;</tspan><tspan
      class="mi" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="5" font-size="72%">p</tspan></text>
  </g>
  <g>
    <rect class="box geom" x="457" y="139" width="220" height="58" rx="7"/>
    <text x="567" y="162" class="cap">eight real numbers</text>
    <text x="567" y="184"><tspan
      class="mi">σ</tspan><tspan
      class="mi" dy="5" font-size="72%">ir</tspan><tspan
      class="mu" dy="-5">,   ∂</tspan><tspan
      class="mi">σ</tspan><tspan
      class="mi" dy="5" font-size="72%">ir</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="5" font-size="72%">p</tspan></text>
  </g>
  <g>
    <rect class="box io" x="10" y="289" width="84" height="58" rx="7"/>
    <text x="52" y="312" class="cap">input</text>
    <text x="52" y="334"><tspan
      class="mi">d</tspan><tspan
      class="mi" dy="5" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" font-size="72%">;</tspan><tspan
      class="mi" font-size="72%">ak</tspan></text>
  </g>
  <g>
    <rect class="box plain heavy" x="117" y="273" width="210" height="90" rx="7"/>
    <text x="222" y="312" class="cap">1. mode sum</text>
    <text x="222" y="334"><tspan
      class="mi">x</tspan><tspan
      class="mi" dy="5" font-size="72%">iak</tspan><tspan
      class="mu" dy="-5"> = Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi">d</tspan><tspan
      class="mi" dy="5" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" font-size="72%">;</tspan><tspan
      class="mi" font-size="72%">ak</tspan><tspan
      class="mu" dy="-5">    </tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="3" font-size="58%">k</tspan><tspan
      class="mi" dy="-8">Y</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" font-size="72%">;</tspan><tspan
      class="mi" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan></text>
  </g>
  <g>
    <rect class="box plain" x="347" y="273" width="210" height="90" rx="7"/>
    <text x="452" y="302" class="cap">2. interpolate in time</text>
    <text x="452" y="323"><tspan
      class="mi">y</tspan><tspan
      class="mi" dy="5" font-size="72%">ijk</tspan><tspan
      class="mu" dy="-5"> = Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">a</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi">S</tspan><tspan
      class="mi" dy="5" font-size="72%">ija</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi">x</tspan><tspan
      class="mi" dy="5" font-size="72%">iak</tspan></text>
    <text x="452" y="344"><tspan
      class="mi">τ</tspan><tspan
      class="mi" dy="5" font-size="72%">ij</tspan><tspan
      class="mu" dy="-5"> = </tspan><tspan
      class="mi">κ</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-12" dx="-4.9" font-size="72%">−1</tspan><tspan
      class="mi" dy="7">t</tspan><tspan
      class="pr">′</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5"> + </tspan><tspan
      class="mi">α</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan></text>
  </g>
  <g>
    <rect class="box plain" x="577" y="273" width="210" height="90" rx="7"/>
    <text x="682" y="302" class="cap">3. transformation laws</text>
    <text x="682" y="323"><tspan
      class="mi">d</tspan><tspan
      class="pr">′</tspan><tspan
      class="mi" dy="5" font-size="72%">ijk</tspan><tspan
      class="mu" dy="-5"> = Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">k</tspan><tspan
      class="pr" font-size="72%">′</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi">M</tspan><tspan
      class="mi" dy="5" font-size="72%">ijkk</tspan><tspan
      class="pr" font-size="72%">′</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi">y</tspan><tspan
      class="mi" dy="5" font-size="72%">ijk</tspan><tspan
      class="pr" font-size="72%">′</tspan></text>
    <text x="682" y="344"><tspan
      class="mu">+ </tspan><tspan
      class="mi">b</tspan><tspan
      class="mi" dy="5" font-size="72%">ijk</tspan></text>
  </g>
  <g>
    <rect class="box plain" x="805" y="273" width="214" height="90" rx="7"/>
    <text x="912" y="302" class="cap">4. compare with target</text>
    <text x="912" y="323"><tspan
      class="mi">J</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5"> = </tspan><tspan
      class="mi">w</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5"> Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">jk</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi">v</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mi" dy="-5">c</tspan><tspan
      class="mi" dy="5" font-size="72%">k</tspan></text>
    <text x="912" y="344"><tspan
      class="mu">× |</tspan><tspan
      class="mi">d</tspan><tspan
      class="pr">′</tspan><tspan
      class="mi" dy="5" font-size="72%">ijk</tspan><tspan
      class="mu" dy="-5"> − </tspan><tspan
      class="mi">d</tspan><tspan
      class="mu" dy="5" font-size="72%">t</tspan><tspan
      class="mu" font-size="72%">;</tspan><tspan
      class="mi" font-size="72%">ijk</tspan><tspan
      class="mu" dy="-5">|</tspan><tspan
      class="mu" dy="-7" font-size="72%">2</tspan></text>
  </g>
  <g>
    <rect class="box rev" x="805" y="431" width="214" height="90" rx="7"/>
    <text x="912" y="460" class="cap">sensitivity to the result</text>
    <text x="912" y="481"><tspan
      class="mu">∂</tspan><tspan
      class="mi">J</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">d</tspan><tspan
      class="pr">′</tspan><tspan
      class="mu" dy="-7" font-size="72%">*</tspan><tspan
      class="mi" dy="12" dx="-6.4" font-size="72%">ijk</tspan><tspan
      class="mu" dy="-5"> =</tspan></text>
    <text x="912" y="502"><tspan
      class="mi">w</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mi" dy="-5">v</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mi" dy="-5">c</tspan><tspan
      class="mi" dy="5" font-size="72%">k</tspan><tspan
      class="mu" dy="-5">(</tspan><tspan
      class="mi">d</tspan><tspan
      class="pr">′</tspan><tspan
      class="mi" dy="5" font-size="72%">ijk</tspan><tspan
      class="mu" dy="-5"> − </tspan><tspan
      class="mi">d</tspan><tspan
      class="mu" dy="5" font-size="72%">t</tspan><tspan
      class="mu" font-size="72%">;</tspan><tspan
      class="mi" font-size="72%">ijk</tspan><tspan
      class="mu" dy="-5">)</tspan></text>
  </g>
  <g>
    <rect class="box rev" x="577" y="431" width="210" height="90" rx="7"/>
    <text x="682" y="460" class="cap">back through the laws</text>
    <text x="682" y="481"><tspan
      class="mu">∂</tspan><tspan
      class="mi">J</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">y</tspan><tspan
      class="mu" dy="-7" dx="1.5" font-size="72%">*</tspan><tspan
      class="mi" dy="12" dx="-6.4" font-size="72%">ijk</tspan><tspan
      class="pr" font-size="72%">′</tspan><tspan
      class="mu" dy="-5"> =</tspan></text>
    <text x="682" y="502"><tspan
      class="mu">Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">k</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi">M</tspan><tspan
      class="mu" dy="-7" dx="1.5" font-size="72%">*</tspan><tspan
      class="mi" dy="12" dx="-6.4" font-size="72%">ijkk</tspan><tspan
      class="pr" font-size="72%">′</tspan><tspan
      class="mu" dy="-5"> ∂</tspan><tspan
      class="mi">J</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">d</tspan><tspan
      class="pr">′</tspan><tspan
      class="mu" dy="-7" font-size="72%">*</tspan><tspan
      class="mi" dy="12" dx="-6.4" font-size="72%">ijk</tspan></text>
  </g>
  <g>
    <rect class="box rev" x="347" y="431" width="210" height="90" rx="7"/>
    <text x="452" y="460" class="cap">back through interpolation</text>
    <text x="452" y="481"><tspan
      class="mu">∂</tspan><tspan
      class="mi">J</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">x</tspan><tspan
      class="mu" dy="-7" dx="1.5" font-size="72%">*</tspan><tspan
      class="mi" dy="12" dx="-6.4" font-size="72%">iak</tspan><tspan
      class="mu" dy="-5"> =</tspan></text>
    <text x="452" y="502"><tspan
      class="mu">Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi">S</tspan><tspan
      class="mi" dy="5" font-size="72%">ija</tspan><tspan
      class="mu" dy="-5"> ∂</tspan><tspan
      class="mi">J</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">y</tspan><tspan
      class="mu" dy="-7" dx="1.5" font-size="72%">*</tspan><tspan
      class="mi" dy="12" dx="-6.4" font-size="72%">ijk</tspan></text>
  </g>
  <g>
    <rect class="box rev heavy" x="117" y="431" width="210" height="90" rx="7"/>
    <text x="222" y="460" class="cap">back through the mode sum</text>
    <text x="222" y="481"><tspan
      class="mu">∂</tspan><tspan
      class="mi">J</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="-5">Y</tspan><tspan
      class="mu" dy="-7" dx="1.5" font-size="72%">*</tspan><tspan
      class="mi" dy="12" dx="-6.4" font-size="72%">i</tspan><tspan
      class="mu" font-size="72%">;</tspan><tspan
      class="mi" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" dy="-5"> =</tspan></text>
    <text x="222" y="502"><tspan
      class="mu">Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">k</tspan><tspan
      class="mu" font-size="72%">:</tspan><tspan
      class="mi" font-size="72%">s</tspan><tspan
      class="mi" dy="3" font-size="58%">k</tspan><tspan
      class="mu" dy="-3" font-size="72%">=</tspan><tspan
      class="mi" font-size="72%">s</tspan><tspan
      class="mu" dy="-5"> Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">a</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi">d</tspan><tspan
      class="mu" dy="-7" dx="1.5" font-size="72%">*</tspan><tspan
      class="mi" dy="12" dx="-6.4" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" font-size="72%">;</tspan><tspan
      class="mi" font-size="72%">ak</tspan><tspan
      class="mu" dy="-5"> ∂</tspan><tspan
      class="mi">J</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">x</tspan><tspan
      class="mu" dy="-7" dx="1.5" font-size="72%">*</tspan><tspan
      class="mi" dy="12" dx="-6.4" font-size="72%">iak</tspan></text>
  </g>
  <g>
    <rect class="box geom" x="580" y="562" width="440" height="76" rx="7"/>
    <text x="800" y="594" class="cap">steps 2 and 3 again, with derivatives in the eight numbers</text>
    <text x="800" y="616"><tspan
      class="mu">∂</tspan><tspan
      class="mi">J</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">σ</tspan><tspan
      class="mi" dy="5" font-size="72%">ir</tspan><tspan
      class="mu" dy="-5"> = 2 </tspan><tspan
      class="mu">Re</tspan><tspan
      class="mu"> Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">jk</tspan><tspan
      class="mu" dy="-5"> (∂</tspan><tspan
      class="mi">J</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">d</tspan><tspan
      class="pr">′</tspan><tspan
      class="mu" dy="-7" font-size="72%">*</tspan><tspan
      class="mi" dy="12" dx="-6.4" font-size="72%">ijk</tspan><tspan
      class="mu" dy="-5">)</tspan><tspan
      class="mu" dy="-7" font-size="72%">*</tspan><tspan
      class="mu" dy="7"> ∂</tspan><tspan
      class="mi">d</tspan><tspan
      class="pr">′</tspan><tspan
      class="mi" dy="5" font-size="72%">ijk</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">σ</tspan><tspan
      class="mi" dy="5" font-size="72%">ir</tspan></text>
  </g>
  <g>
    <rect class="box plain" x="70" y="673" width="900" height="62" rx="7"/>
    <text x="520" y="698" class="cap">chain rule, summed over directions</text>
    <text x="520" y="720"><tspan
      class="mu">∂</tspan><tspan
      class="mi">J</tspan><tspan
      class="mu">/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="5" font-size="72%">p</tspan><tspan
      class="mu" dy="-5"> = (1/𝒩) Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5"> [ 2 </tspan><tspan
      class="mu">Re</tspan><tspan
      class="mu"> Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" dy="-5"> (∂</tspan><tspan
      class="mi">J</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="-5">Y</tspan><tspan
      class="mu" dy="-7" dx="1.5" font-size="72%">*</tspan><tspan
      class="mi" dy="12" dx="-6.4" font-size="72%">i</tspan><tspan
      class="mu" font-size="72%">;</tspan><tspan
      class="mi" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" dy="-5">)</tspan><tspan
      class="mu" dy="-7" font-size="72%">*</tspan><tspan
      class="mu" dy="7"> ∂</tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="-5">Y</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" font-size="72%">;</tspan><tspan
      class="mi" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="5" font-size="72%">p</tspan><tspan
      class="mu" dy="-5"> + Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">r</tspan><tspan
      class="mu" dy="-5"> ∂</tspan><tspan
      class="mi">J</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">σ</tspan><tspan
      class="mi" dy="5" font-size="72%">ir</tspan><tspan
      class="mu" dy="-5"> ∂</tspan><tspan
      class="mi">σ</tspan><tspan
      class="mi" dy="5" font-size="72%">ir</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="5" font-size="72%">p</tspan><tspan
      class="mu" dy="-5"> ]</tspan></text>
  </g>

  <g class="legend">
    <rect class="box geom" x="110" y="758" width="22" height="16" rx="3"/>
    <text x="140" y="771" class="cap" style="text-anchor: start">forward mode: dual numbers, cheap</text>
    <rect class="box rev" x="420" y="758" width="22" height="16" rx="3"/>
    <text x="450" y="771" class="cap" style="text-anchor: start">reverse pass, by hand</text>
    <rect class="box plain heavy" x="650" y="758" width="22" height="16" rx="3"/>
    <text x="680" y="771" class="cap" style="text-anchor: start">a sum over every mode and time</text>
  </g>
  </g>
</svg>
</div>
```

The top row computes the geometry of the direction, with derivatives
with respect to all the parameters.  The middle row is the objective
itself, and the row below it runs the same steps backwards.  The boxes
with heavy borders are the two sums over every mode and every time,
which dominate the cost.

## Where the parameters enter

Neither the data nor the target depends on ``θ``.  Within direction
``i``, ``θ`` enters only through ``{}_sY_{i;ℓ,m}`` and ``σ_{ir}``, so the
chain rule gives

```math
\frac{∂J}{∂θ_p}
= \frac{1}{𝒩} \sum_i \left[
    2 \operatorname{Re} \sum_{s,ℓ,m}
    \left(\frac{∂J_i}{∂\,{}_sY^*_{i;ℓ,m}}\right)^*
    \frac{∂\,{}_sY_{i;ℓ,m}}{∂θ_p}
  + \sum_r \frac{∂J_i}{∂σ_{ir}}\, \frac{∂σ_{ir}}{∂θ_p}
\right].
```

There are three kinds of factor here, computed in three different
ways:

- The derivatives of the geometry, ``∂\,{}_sY_{i;ℓ,m}/∂θ_p`` and
  ``∂σ_{ir}/∂θ_p``, come from forward-mode differentiation in all the
  parameters.  They are cheap: a few operations per number and
  parameter.
- The sensitivity to the harmonics, ``∂J_i/∂\,{}_sY^*_{i;ℓ,m}``, comes from
  the reverse pass through steps 3, 2, and 1, holding ``σ_{ir}`` fixed.
- The sensitivity to the eight numbers, ``∂J_i/∂σ_{ir}``, comes from
  forward-mode differentiation of steps 2 and 3 alone, in those eight
  numbers.

Each ``J_i`` depends only on the quantities of direction ``i``.  So
each direction is computed and then differentiated at once, and nothing
has to be stored from one direction to the next.

## The reverse pass

**From the objective:**

```math
\frac{∂J_i}{∂d'^*_{ijk}}
= w_i\, v_j\, c_k \left( d'_{ijk} - d_{\mathrm{t};ijk} \right).
```

**Back through the transformation laws:**

```math
\frac{∂J_i}{∂y^*_{ijk'}}
= \sum_k M^*_{ijkk'}\, \frac{∂J_i}{∂d'^*_{ijk}}.
```

The offsets ``b_{ijk}`` do not depend on the data, so they contribute
nothing here.  Their dependence on the geometry is not lost: it enters
through ``∂J_i/∂σ_{ir}``, below.

**Back through the interpolation:**

```math
\frac{∂J_i}{∂x^*_{iak}}
= \sum_j S_{ija}\, \frac{∂J_i}{∂y^*_{ijk}}.
```

The sum over ``j`` runs the interpolation backwards, so it needs the
structure of ``S_{ija}``.  If ``t_n ≤ τ_{ij} < t_{n+1}``, the spline's
value is a combination of the two samples ``x_{ink}`` and
``x_{i,n+1,k}`` and of the spline's second derivatives
``\ddot{x}_{ink}`` and ``\ddot{x}_{i,n+1,k}`` there:

```math
y_{ijk} = \sum_a E_{ija}\, x_{iak} + \sum_b F_{ijb}\, \ddot{x}_{ibk},
```

where ``E_{ija}`` and ``F_{ijb}`` are nonzero only for
``a, b ∈ \{n, n+1\}``.  The second derivatives at the interior knots
solve a tridiagonal system, and vanish at the ends:

```math
\sum_{b'} Q_{bb'}\, \ddot{x}_{ib'k} = \sum_a B_{ba}\, x_{iak},
```

where ``Q_{bb'}`` is symmetric and ``B_{ba}`` takes second
differences.  So

```math
S_{ija} = E_{ija} + \sum_{b,b'} F_{ijb}\, \left(Q^{-1}\right)_{bb'}\, B_{b'a},
```

and, writing ``g_{ijk}`` for ``∂J_i/∂y^*_{ijk}`` and using the
symmetry of ``Q``,

```math
\sum_j S_{ija}\, g_{ijk}
= \sum_j E_{ija}\, g_{ijk} + \sum_{b'} B_{b'a}\, z_{ib'k},
\qquad
\sum_{b'} Q_{bb'}\, z_{ib'k} = \sum_j F_{ijb}\, g_{ijk}.
```

That is, each ``g_{ijk}`` is spread onto the samples and second
derivatives on either side of ``τ_{ij}``, the second-derivative parts
go through the same tridiagonal solve, and the result is spread by the
second differences.  [`spline_transpose!`](@ref Scri.spline_transpose!)
takes ``E_{ija}`` and ``F_{ijb}`` from the same evaluation function,
and chooses each time's interval in the same way, as the forward
interpolation, so the two agree to roundoff.

**Back through the mode sum:**  ``x_{iak}`` is linear in the
harmonics, with the modes as coefficients, so

```math
\frac{∂J_i}{∂\,{}_sY^*_{i;ℓ,m}}
= \sum_{k:\, s_k = s}\ \sum_a d^*_{ℓ,m;ak}\, \frac{∂J_i}{∂x^*_{iak}}.
```

The sum over ``k`` collects the components that share a spin weight,
and so share harmonics.  This is a sum of the same size as step 1, and
it is the second of the two expensive operations.

## The eight numbers of the geometry

The sensitivity of ``J_i`` to the eight numbers of direction ``i`` is

```math
\frac{∂J_i}{∂σ_{ir}}
= 2 \operatorname{Re} \sum_{j,k}
  \left(\frac{∂J_i}{∂d'^*_{ijk}}\right)^*
  \left.\frac{∂d'_{ijk}}{∂σ_{ir}}\right|_{x},
```

where the derivative holds the original time series ``x_{iak}``
fixed.  The eight numbers enter steps 2 and 3 in three ways:
``κ_i^{-1}`` and ``α_i`` move the interpolation times ``τ_{ij}``;
``κ_i^{-1}`` and the null-rotation parameter set ``M_{ijkk'}``, which
also depends on ``τ_{ij}``; and ``κ_i^{-1}`` and ``ð^2α_i`` set the
offsets ``b_{ijk}``.  Rather than differentiating each of these by
hand, steps 2 and 3 are run once more for the direction, with dual
numbers that hold derivatives with respect to the eight numbers.  The
spline's second derivatives depend only on ``x_{iak}``, so they stay
plain numbers.  This costs a small multiple of steps 2 and 3, with no
mode sum, and it includes the offsets of the shear and the strain.

## Assembling the gradient

For the derivatives of the geometry, ``θ`` is turned into dual numbers
holding derivatives with respect to all the parameters.  The same
functions that compute the objective's geometry then return it with
its derivatives: [`bms_from_parameters`](@ref), the boost and rotation
of the grid, [`supertranslation_on_grid`](@ref
Scri.supertranslation_on_grid), [`compute_ðt′╱2κ`](@ref
Scri.compute_ðt′╱2κ), the conformal factor, and the recursion for the
harmonics at ``R_i``.  These are combined with ``∂J_i/∂\,{}_sY^*_{i;ℓ,m}``
and ``∂J_i/∂σ_{ir}`` as in the formula above, and summed over the
directions, in parallel over threads.

The harmonics depend only on the boost and the rotation, not on the
``a``, but they hold derivatives with respect to all the parameters.
That extra work is a few percent of the total.

## Comparing only the modes up to ``ℓ_\mathrm{fit}``

With ``ℓ_\mathrm{fit} = L``, the objective compares projections onto
the modes with ``ℓ ≤ L``, using the harmonics at the grid's own
directions, ``{}_sY'_{i;ℓ,m} = {}_sY_{ℓ,m}(R'_i)``:

```math
J = \frac{1}{𝒩} \sum_j v_j \sum_k c_k \sum_{ℓ=|s_k|}^{L} \sum_m
    \left| Π_{ℓ,m;jk} - Π_{\mathrm{t};ℓ,m;jk} \right|^2,
\qquad
Π_{ℓ,m;jk} = \sum_i w_i\, {}_{s_k}Y'^*_{i;ℓ,m}\, d'_{ijk},
```

where ``Π_{\mathrm{t};ℓ,m;jk}`` is the target's projection.  This is
not a sum of terms that each depend on one direction.  But ``J``
depends on ``d'_{ijk}`` only through ``Π_{ℓ,m;jk}``, which is linear in
it, so the rule above gives

```math
\frac{∂(𝒩 J)}{∂d'^*_{ijk}}
= w_i\, v_j\, c_k \sum_{ℓ=|s_k|}^{L} \sum_m
  \left( Π_{ℓ,m;jk} - Π_{\mathrm{t};ℓ,m;jk} \right) {}_{s_k}Y'_{i;ℓ,m}:
```

the residual of the projection, evaluated back in the direction
``R'_i``.  Each direction's reverse pass starts from this in place of
``∂J_i/∂d'^*_{ijk}``.  Since it needs ``Π_{ℓ,m;jk}`` before any
direction can start, a first pass over the directions computes the
projection, at the cost of one evaluation of the objective.

## Cost and memory

Per direction, the objective costs one mode sum and the gradient two,
plus work that involves no mode sum: the second run of steps 2 and 3
with eight derivatives, the harmonics' derivatives, and the geometry.
Forward-mode differentiation of the objective instead carries a
derivative for every parameter through the mode sum.  Measured with
strain and ``ψ_4`` data, ``N_t = 400``, and 8 threads, as multiples
of the time for one evaluation of the objective:

| ``ℓ_\mathrm{max}`` | ``ℓ_\mathrm{fit}`` | ``P`` | gradient | ForwardDiff |
| :---: | :---: | :---: | :---: | :---: |
| 8 | — | 7 | 2.4 | 9.6 |
| 8 | — | 22 | 2.4 | 21.7 |
| 8 | 4 | 7 | 3.0 | 10.2 |
| 8 | 4 | 22 | 3.1 | 29.8 |
| 16 | — | 7 | 1.6 | 6.6 |
| 16 | — | 22 | 1.6 | 14.4 |
| 16 | 8 | 7 | 2.5 | 7.1 |
| 16 | 8 | 22 | 2.7 | 18.0 |

The gradient's multiple is smaller at the larger band limit, where the
mode sums take a larger share of the time.  With ``ℓ_\mathrm{fit}``,
the first pass adds about one evaluation.

No result of one direction is needed by another, so each thread
reuses a few buffers the size of one direction's time series and
harmonics.  A general reverse-mode tool would instead record every
operation of every mode sum, of which there are
``N_\mathrm{dir}\, N_t\, N_\mathrm{comp}\, N_\mathrm{modes}``.

## How it is checked

The tests in `test/gradient.jl` check each piece and the whole:

- ``\sum_j \sum_a X_j^*\, S_{ija}\, Z_a
  = \sum_a \left(\sum_j S_{ija}\, X_j\right)^* Z_a`` for the spline's
  backward sum, for random ``X_j`` and ``Z_a`` on uneven grids, from
  four knots up, with times on the knots and at the ends;
- ``M_{ijkk'}`` and ``b_{ijk}`` against [`mix_components!`](@ref), for
  components that mix, on both null infinities and in more than one
  convention;
- the value and the gradient against ForwardDiff through
  [`transform_objective`](@ref), at the identity and away from it, for
  components that mix, on both null infinities, with weights,
  normalization, and ``ℓ_\mathrm{fit}``; they agree to about
  ``10^{-15}``, and to about ``10^{-25}`` in `Double64`;
- the gradient against central differences, independently of
  ForwardDiff.
