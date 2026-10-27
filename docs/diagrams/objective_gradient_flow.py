# Generate the inline SVG for the diagram on the documentation page "The gradient of the
# objective", docs/src/40-Details/06-objective_gradient.md.
#
#     python3 docs/diagrams/objective_gradient_flow.py           # update the page in place
#     python3 docs/diagrams/objective_gradient_flow.py out.svg   # or write the SVG to a file
#
# Math is written in a tiny markup: `_{...}` or `_x` is a subscript, `^{...}` or `^x` a
# superscript; Latin and Greek letters are set in the math italic font, everything else
# upright.  The text is placed by estimated character widths, so check the page in a
# browser after changing a label.
import html, os, re, sys

ITALIC = re.compile(r"[A-Za-zα-ωϑϕΑ-ΩΣ∇ℓ]")
UPRIGHT_LETTERS = set("ΣΔ∇ð∂")


SUB = (5, 3, 2)  # baseline drop of a subscript, by nesting level
SUP = (-7, -4, -3)  # baseline rise of a superscript, by nesting level
SIZE = (None, "72%", "58%")


def runs(text):
    """Split `text` into italic and upright runs; a prime is a run of its own."""
    buf, cur = "", None
    for ch in text:
        if ch == "′":
            if buf:
                yield buf, cur
            yield ch, False
            buf, cur = "", None
            continue
        it = bool(ITALIC.match(ch)) and ch not in UPRIGHT_LETTERS
        if cur is None or it == cur:
            buf += ch
            cur = it
        else:
            yield buf, cur
            buf, cur = ch, it
    if buf:
        yield buf, cur


CHAR_W = (8.5, 6.4, 5.2)  # rough advance of one character, by nesting level, in pixels


def math(s):
    """Return tspans for the markup `s`, which may nest scripts.  A subscript and a
    superscript on the same base are stacked, using a rough estimate of their widths."""
    segs = []  # [text, italic, level, baseline offset, dx]

    def width(lo, hi):
        return sum(len(t) * CHAR_W[min(lv, 2)] for t, _, lv, _, _ in segs[lo:hi])

    def parse(s, level, offset):
        i = 0
        prev = None  # (kind, first segment, end segment) of a script just closed
        while i < len(s):
            ch = s[i]
            if ch in "_^":
                if s[i + 1] == "{":
                    depth, j = 0, i + 1
                    while True:
                        depth += {"{": 1, "}": -1}.get(s[j], 0)
                        if depth == 0:
                            break
                        j += 1
                    body, i = s[i + 2 : j], j + 1
                else:
                    body, i = s[i + 1], i + 2
                shift = (SUB if ch == "_" else SUP)[min(level, 2)]
                first = len(segs)
                parse(body, level + 1, offset + shift)
                if ch == "^" and first > 0 and len(segs) > first and segs[first - 1][1]:
                    segs[first][4] += 1.5  # italic correction
                if prev and prev[0] != ch and prev[2] == first and len(segs) > first:
                    # Stack this script under or over the previous one.
                    w_prev = width(prev[1], prev[2])
                    w_this = width(first, len(segs))
                    segs[first][4] -= w_prev
                    pending_dx[0] = max(0.0, w_prev - w_this)
                    prev = None
                else:
                    prev = (ch, first, len(segs))
            else:
                prev = None
                if ch == "«":
                    j = s.index("»", i)
                    new = [(s[i + 1 : j], False)]
                    i = j + 1
                else:
                    j = i
                    while j < len(s) and s[j] not in "_^«":
                        j += 1
                    new = list(runs(s[i:j]))
                    i = j
                for text, it in new:
                    segs.append([text, it, level, offset, pending_dx[0]])
                    pending_dx[0] = 0.0

    pending_dx = [0.0]  # owed to the next segment after a stacked pair
    parse(s, 0, 0)
    out, cur = [], 0
    for text, it, level, off, dx in segs:
        cls = "pr" if text == "′" else ("mi" if it else "mu")
        attrs = f'class="{cls}"'
        if off != cur:
            attrs += f' dy="{off - cur}"'
            cur = off
        if dx:
            attrs += f' dx="{dx:.1f}"'
        if level:
            attrs += f' font-size="{SIZE[min(level, 2)]}"'
        # The line break sits inside the tag, where it is not text, to keep lines short.
        out.append(f"<tspan\n      {attrs}>{html.escape(text)}</tspan>")
    return "".join(out)


def text_lines(cx, cy, lines):
    """Lines of a box, centered on (cx, cy).  Each line is ('cap', str) or ('math', str)."""
    gap = 21
    y0 = cy - gap * (len(lines) - 1) / 2 + 5
    parts = []
    for n, (kind, s) in enumerate(lines):
        y = y0 + n * gap
        if kind == "cap":
            # A Greek function name in a caption is set in the math font, where ν is not v.
            cap = re.sub(r"([σιμν])", r'<tspan class="mi">\1</tspan>', html.escape(s))
            parts.append(f'<text x="{cx}" y="{y:.0f}" class="cap">{cap}</text>')
        else:
            parts.append(f'<text x="{cx}" y="{y:.0f}">{math(s)}</text>')
    return "\n    ".join(parts)


boxes = {}


def box(name, cx, cy, w, h, cls, lines):
    boxes[name] = (cx, cy, w, h)
    x, y = cx - w / 2, cy - h / 2
    return (
        f'<g>\n    <rect class="box {cls}" x="{x:.0f}" y="{y:.0f}" width="{w}" height="{h}" rx="7"/>\n    '
        + text_lines(cx, cy, lines)
        + "\n  </g>"
    )


def edge(points, cls="", label=None, label_at=None, anchor="start"):
    d = "M " + " L ".join(f"{x:.0f},{y:.0f}" for x, y in points)
    s = f'<path class="edge {cls}" d="{d}"/>'
    if label:
        lx, ly = label_at
        s += f'\n  <text x="{lx}" y="{ly}" class="lab" style="text-anchor: {anchor}">{math(label)}</text>'
    return s


C1, C2, C3, C4, C5 = 52, 222, 452, 682, 912
W, W5 = 210, 214
R0, R1, R2, R3, R4, R5 = 54, 168, 318, 476, 600, 704
CS = (C3 + C4) // 2  # the eight numbers feed steps 2 and 3
CB = 800  # the eight numbers' sensitivity
H2 = 90  # height of the boxes of the middle rows

parts = []
# Geometry: forward mode in all parameters.
parts.append(box("θ", 30, R0, 44, 44, "io", [("math", "θ")]))
parts.append(box("vRα", 118, R0, 96, 44, "geom", [("math", "v⃗, R, α")]))
parts.append(
    box(
        "geom",
        610,
        R0,
        830,
        84,
        "geom",
        [
            ("cap", "for each pixel of the new frame’s grid: the original direction, the conformal factor,"),
            ("cap", "the supertranslation and its derivatives there, and the null-rotation parameter"),
            ("math", "R_j,   κ_j,   α_j,   ðα_j,   ð^2α_j,   b_j^{(0)} + t b_j^{(1)}"),
        ],
    )
)
parts.append(box("Y", C2, R1, W, 58, "geom", [("cap", "harmonics"), ("math", "_{s}Y^{(ℓ,m)},   ∂(_{s}Y^{(ℓ,m)})/∂θ^p")]))
parts.append(box("σ", CS, R1, 220, 58, "geom", [("cap", "eight real numbers"), ("math", "g^n,   ∂g^n/∂θ^p")]))

# Forward pass on plain numbers.
parts.append(box("D", C1, R2, 84, 58, "io", [("cap", "input"), ("math", "^{k}d^h_{(ℓ,m)}")]))
parts.append(
    box("sum", C2, R2, W, H2, "plain heavy", [("cap", "1. synthesis, σ"), ("math", "^{k}x^h = ^{k}d^h_{(ℓ,m)} _{s_k}Y^{(ℓ,m)}")])
)
parts.append(
    box(
        "interp",
        C3,
        R2,
        W,
        H2,
        "plain",
        [("cap", "2. interpolation, ι"), ("math", "^{k}y^i = S^i_h ^{k}x^h"), ("math", "τ_i = κ^{−1}t′_i + α")],
    )
)
parts.append(
    box(
        "laws",
        C4,
        R2,
        W,
        H2,
        "plain",
        [("cap", "3. mixing, μ"), ("math", "^{k}d^{′i} = M^k_{k′} ^{k′}y^i"), ("math", "+ q^k")],
    )
)
parts.append(
    box(
        "cmp",
        C5,
        R2,
        W5,
        H2,
        "plain",
        [("cap", "4. norm, ν"), ("math", "Φ_j = w_j Σ_{ik} v_iu_k"), ("math", "× |^{k}d^{′i} − ^{k}c^i|^2")],
    )
)

# Reverse pass on plain numbers.
parts.append(
    box(
        "start",
        C5,
        R3,
        W5,
        H2,
        "rev",
        [("cap", "back through ν"), ("math", "(Ω_{d′})_{ki} ="), ("math", "2w_jv_iu_k(^{k}d^{′i} − ^{k}c^i)^*")],
    )
)
parts.append(
    box(
        "MH",
        C4,
        R3,
        W,
        H2,
        "rev",
        [("cap", "back through μ"), ("math", "(Ω_y)_{k′i} ="), ("math", "(Ω_{d′})_{ki} M^k_{k′}")],
    )
)
parts.append(
    box(
        "ST",
        C3,
        R3,
        W,
        H2,
        "rev",
        [("cap", "back through ι"), ("math", "(Ω_x)_{kh} ="), ("math", "(Ω_y)_{ki} S^i_h")],
    )
)
parts.append(
    box(
        "sum2",
        C2,
        R3,
        W,
        H2,
        "rev heavy",
        [
            ("cap", "back through σ"),
            ("math", "_{s}(Ω_Y)_{(ℓ,m)} ="),
            ("math", "Σ_{k:s_k=s} (Ω_x)_{kh} ^{k}d^h_{(ℓ,m)}"),
        ],
    )
)

# The eight numbers, and the chain rule.
parts.append(
    box(
        "σbar",
        CB,
        R4,
        440,
        76,
        "geom",
        [
            ("cap", "ι and μ again, with derivatives in the eight numbers"),
            ("math", "(Ω_g)_n = «Re» [(Ω_{d′})_{ki} ∂(^{k}d^{′i})/∂g^n]"),
        ],
    )
)
parts.append(
    box(
        "chain",
        520,
        R5,
        900,
        62,
        "plain",
        [
            ("cap", "chain rule, summed over pixels"),
            ("math", "∂Φ/∂θ^p = (1/𝒩) Σ_j [ Σ_s «Re» [_{s}(Ω_Y)_{(ℓ,m)} ∂(_{s}Y^{(ℓ,m)})/∂θ^p] + (Ω_g)_n ∂g^n/∂θ^p ]"),
        ],
    )
)

h0, h1, h2 = 42, 29, H2 // 2  # half-heights of rows 0, 1, and 2 and 3
edges = []
edges.append(edge([(52, R0), (70, R0)]))
edges.append(edge([(166, R0), (195, R0)]))
edges.append(edge([(C2, R0 + h0), (C2, R1 - h1)]))
edges.append(edge([(CS, R0 + h0), (CS, R1 - h1)]))
edges.append(edge([(C2, R1 + h1), (C2, R2 - h2)]))
edges.append(edge([(CS - 40, R1 + h1), (C3 + 24, R2 - h2)]))
edges.append(edge([(CS + 40, R1 + h1), (C4 - 24, R2 - h2)]))
edges.append(edge([(C1 + 42, R2), (C2 - W // 2, R2)]))
edges.append(edge([(C2 + W // 2, R2), (C3 - W // 2, R2)]))
edges.append(edge([(C3 + W // 2, R2), (C4 - W // 2, R2)]))
edges.append(edge([(C4 + W // 2, R2), (C5 - W5 // 2, R2)]))
edges.append(edge([(C5 + W5 // 2, R2), (1046, R2)], label="Φ", label_at=(1052, R2 + 5)))
edges.append(
    edge([(C5, R2 + h2), (C5, R3 - h2)], label="^{k}d^{′i} − ^{k}c^i", label_at=(C5 + 9, (R2 + R3) / 2 + 5))
)
edges.append(edge([(C5 - W5 // 2, R3), (C4 + W // 2, R3)], "rev"))
edges.append(edge([(C4 - W // 2, R3), (C3 + W // 2, R3)], "rev"))
edges.append(edge([(C3 - W // 2, R3), (C2 + W // 2, R3)], "rev"))
edges.append(edge([(C5 - 20, R3 + h2), (C5 - 50, R4 - 38)]))
edges.append(edge([(C2, R3 + h2), (C2, R5 - 31)]))
edges.append(edge([(CB, R4 + 38), (CB, R5 - 31)]))

legend_y = 770
legend = f"""
  <g class="legend">
    <rect class="box geom" x="110" y="{legend_y - 12}" width="22" height="16" rx="3"/>
    <text x="140" y="{legend_y + 1}" class="cap" style="text-anchor: start">forward mode: dual numbers, cheap</text>
    <rect class="box rev" x="420" y="{legend_y - 12}" width="22" height="16" rx="3"/>
    <text x="450" y="{legend_y + 1}" class="cap" style="text-anchor: start">reverse pass, by hand</text>
    <rect class="box plain heavy" x="650" y="{legend_y - 12}" width="22" height="16" rx="3"/>
    <text x="680" y="{legend_y + 1}" class="cap" style="text-anchor: start">a sum over every mode and time</text>
  </g>"""

svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1070 790" width="100%" style="max-width: 900px;" role="img"
     aria-label="The flow of the objective and its gradient through each direction of the grid">
  <style>
    .grad-flow text {{ fill: currentColor; text-anchor: middle; font-size: 17px; white-space: pre; }}
    .grad-flow text.cap {{ font-size: 14.5px; }}
    .grad-flow text.lab {{ font-size: 16px; }}
    .grad-flow .mi {{ font-family: KaTeX_Math, "Times New Roman", serif; font-style: italic; }}
    .grad-flow .mu {{ font-family: KaTeX_Main, "Times New Roman", serif; }}
    .grad-flow .pr {{ font-family: "Times New Roman", "Liberation Serif", serif; }}
    .grad-flow .box {{ stroke: currentColor; stroke-width: 1.2; fill: none; }}
    .grad-flow .geom {{ fill: rgba(56, 132, 255, 0.17); }}
    .grad-flow .rev {{ fill: rgba(255, 140, 0, 0.22); }}
    .grad-flow .plain {{ fill: rgba(128, 128, 128, 0.10); }}
    .grad-flow .heavy {{ stroke-width: 3; }}
    .grad-flow .edge {{ stroke: currentColor; stroke-width: 1.3; fill: none; marker-end: url(#grad-flow-arrow); }}
  </style>
  <defs>
    <marker id="grad-flow-arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7"
            orient="auto-start-reverse">
      <path d="M0,0 L10,5 L0,10 z" fill="currentColor"/>
    </marker>
  </defs>
  <g class="grad-flow">
  {chr(10).join("  " + e for e in edges)}
  {chr(10).join("  " + p for p in parts)}
  {legend}
  </g>
</svg>"""
svg = "\n".join(line.rstrip() for line in svg.split("\n"))
if len(sys.argv) > 1:
    open(sys.argv[1], "w").write(svg)
else:
    page = os.path.join(os.path.dirname(__file__), "..", "src", "40-Details", "06-objective_gradient.md")
    text = open(page).read()
    i, j = text.index("<svg "), text.index("</svg>") + len("</svg>")
    open(page, "w").write(text[:i] + svg + text[j:])
