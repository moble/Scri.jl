# Claude Code — Project Guidance for Scri.jl

## Physics conventions

- Metric signature is **-+++** (Minkowski: `ℝ^{3,1}`).
- Units where `G = c = 1` throughout.
- We use **Geometric Algebra** (= Clifford Algebra, but with geometric
  emphasis and developed over `ℝ`).  Complex structures arise naturally
  within the algebra; do not introduce them _ad hoc_.
- Always speak in terms of **Spin groups**.  Never use the accidental
  isomorphisms SU(2) or SL(2,C) — even for the double cover of the
  proper orthochronous Lorentz group, use `Spin⁺(3,1)` (or just "the
  Spin group") and the language of rotors.
- Rotations and boosts are represented as **rotors** acting by
  conjugation: `𝐯' = R𝐯R̃`.  Composition of transformations is rotor
  multiplication.
- The proper orthochronous Lorentz group is the quotient of its Spin
  double cover; prefer to work at the level of the cover.

## Writing style

- Two spaces between sentences in all prose (documentation, docstrings,
  comments).  We write in monospace contexts and need the visual
  separation.
- Documentation lives in Documenter.jl Markdown files under `docs/src/`.

## Julia code style

- **Strongly prefer Unicode** in both source code and documentation.
  Use literal `α`, `β`, `γ`, `𝐯`, `𝐰`, `∈`, `∧`, `⊗`, … everywhere
  rather than ASCII spellings.
- Formatter: `JuliaFormatter` with the `blue` style (see
  `.JuliaFormatter.toml`).  Notable: `for … ∈` (not `in`).
- In LaTeX inside docstrings or Markdown, write literal Unicode rather
  than `\alpha`, `\beta`, etc.  MathJax/KaTeX handle Unicode directly,
  and it makes the raw docstring far more readable.  E.g., write
  `` `α` `` not `` `\alpha` ``.
- Follow existing patterns: exported names are `PascalCase` types and
  `snake_case` functions; internal helpers are prefixed with `_`.

## Running tests

Tests use `TestItemRunner.jl`.  Individual tests are `@testitem` blocks
scattered through the source; `test/runtests.jl` discovers and runs them.

**While iterating on a specific piece of code**, run only the relevant
subset to keep the feedback loop tight:

```sh
# by name substring
julia --project=test test/runtests.jl --name "rotor composition"

# by filename substring
julia --project=test test/runtests.jl --file lorentz

# by pattern (matches name OR filename)
julia --project=test test/runtests.jl --pattern "boost"

# by tag (AND logic — all listed tags must be present)
julia --project=test test/runtests.jl --tags unit,fast

# add --verbose / -v to see individual test results
julia --project=test test/runtests.jl --pattern "boost" --verbose
```

To see what tags are defined: `julia --project=test test/runtests.jl --list-tags`

**Run the full test suite** before any significant milestone: completing a
feature, declaring a bug fixed, or making a commit that touches more than
one logical unit.  Full suite:

```sh
julia --project=test test/runtests.jl
```

## Key dependencies

- `Quaternionic.jl` — quaternion / rotor arithmetic.
- `Grassmann.jl` / Geometric Algebra primitives live under `Grassmann/`.
- Documentation built with Documenter.jl (`docs/make.jl`).

## What to avoid

- Do **not** introduce `SU(2)`, `SL(2,C)`, or any accidental-isomorphism
  language when a Spin-group description is available.
- Do **not** call the Newman-Penrose quantities `Ψ₀`…`Ψ₄` the **Weyl
  scalars**.  They are **Weyl components**.  Calling them scalars misleads
  readers into thinking they transform as scalars under Lorentz
  transformations; they do not.
- Do **not** use `\alpha`-style LaTeX escapes when the Unicode character
  can be written directly.
- Do **not** use `for x in` — always `for x ∈`.
