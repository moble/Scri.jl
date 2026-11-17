# Generate the inline SVG for the diagram on the documentation page "The gradient of the
# objective", docs/src/40-Details/06-objective_gradient.md.
#
#     python3 docs/diagrams/objective_gradient_flow.py           # update the page in place
#     python3 docs/diagrams/objective_gradient_flow.py out.svg   # or write the SVG to a file
#     python3 docs/diagrams/objective_gradient_flow.py page.md   # or update another copy
#
# Math is written in a tiny markup: `_{...}` or `_x` is a subscript, `^{...}`
# or `^x` a superscript; Latin and Greek letters are set in the math italic font,
# everything else upright; `«...»` is upright text, and an empty `«»` keeps two scripts
# from being stacked.  A combining macron after a letter puts a bar over it, drawn as a
# spacing macron moved back over the letter, because the math fonts have no combining
# marks.  A subscript is placed under a superscript by estimated character widths, so
# check the page in a browser after changing a label.
import html, os, re, sys, unicodedata

ITALIC = re.compile(r"[A-Za-zα-ωϑϕΑ-ΩΣ∇ℓ]")
UPRIGHT_LETTERS = set("ΣΔΦ∇ð∂")


SUB = (5, 3, 2)  # baseline drop of a subscript, by nesting level
SUP = (-7, -4, -3)  # baseline rise of a superscript, by nesting level
SIZE = (None, "72%", "58%")


def runs(text):
    """Split `text` into italic and upright runs; a prime is a run of its own, and a
    combining mark stays with the letter before it."""
    buf, cur = "", None
    for ch in text:
        if unicodedata.combining(ch):
            # The letter before the mark becomes a run of its own, with the mark.
            if len(buf) > 1:
                yield buf[:-1], cur
            yield buf[-1:] + ch, cur
            buf, cur = "", None
            continue
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


FONT_PX = (17.0, 12.24, 9.86)  # font size by nesting level: 17 px, then 72% and 58% of it
# Advance widths in ems, close to those of the KaTeX fonts; anything else counts as half an em.
# They are only used to place a subscript under a superscript, so they need not be exact.
EM = dict(
    a=0.53, b=0.43, c=0.43, d=0.52, e=0.47, f=0.49, g=0.48, h=0.58, i=0.35, j=0.41, k=0.52,
    l=0.30, m=0.88, n=0.60, o=0.49, p=0.50, q=0.45, r=0.45, s=0.47, t=0.36, u=0.57, v=0.49,
    w=0.72, x=0.57, y=0.49, z=0.47, M=0.97, S=0.61, W=0.94, Y=0.58, P=0.64,
)
EM.update({"ℓ": 0.42, "θ": 0.47, "ι": 0.35, "κ": 0.58, "μ": 0.60, "ν": 0.49, "σ": 0.57, "τ": 0.44,
           "′": 0.28, ":": 0.28, ",": 0.28, " ": 0.25, "=": 0.78, "(": 0.39, ")": 0.39, "∂": 0.53})
BAR = "ˉ"  # a spacing macron, about 8.5 pixels wide
BAR_BACK = {"μ": 9.5, "c": 8.0}  # how far back the bar moves to sit over each letter


def advance(text, level):
    """The estimated width of `text` in pixels at nesting level `level`."""
    if text == BAR:
        return 0.0
    return FONT_PX[min(level, 2)] * sum(EM.get(ch, 0.5) for ch in text if not unicodedata.combining(ch))


def math(s):
    """Return tspans for the markup `s`, which may nest scripts.  A subscript and a
    superscript on the same base are stacked, using a rough estimate of their widths."""
    segs = []  # [text, italic, level, baseline offset, dx]

    def width(lo, hi):
        return sum(advance(t, lv) for t, _, lv, _, _ in segs[lo:hi])

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
                before = first - 2 if first > 1 and segs[first - 1][0] == BAR else first - 1
                if ch == "^" and first > 0 and len(segs) > first and segs[before][1]:
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
                    barred = len(text) > 1 and unicodedata.combining(text[-1])
                    base = text[:-1] if barred else text
                    segs.append([base, it, level, offset, pending_dx[0]])
                    pending_dx[0] = 0.0
                    if barred:
                        back = BAR_BACK[base]
                        segs.append([BAR, False, level, offset - 1.5, -back])
                        pending_dx[0] = back - 7.5  # return to the end of the letter

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


def caption(s):
    # A Greek function name in a caption is set in the math font, where ν is not v.
    return re.sub(r"([σιμνθ])", r'<tspan class="mi">\1</tspan>', html.escape(s))


def text_lines(cx, cy, lines):
    """Lines of a box, centered on (cx, cy).  Each line is ('cap', str) or ('math', str)."""
    gap = 22
    y0 = cy - gap * (len(lines) - 1) / 2 + 5
    parts = []
    for n, (kind, s) in enumerate(lines):
        y = y0 + n * gap
        if kind == "cap":
            parts.append(f'<text x="{cx}" y="{y:.0f}" class="cap">{caption(s)}</text>')
        else:
            parts.append(f'<text x="{cx}" y="{y:.0f}">{math(s)}</text>')
    return "\n    ".join(parts)


def box(cx, cy, w, h, cls, lines):
    x, y = cx - w / 2, cy - h / 2
    return (
        f'<g>\n    <rect class="box {cls}" x="{x:.0f}" y="{y:.0f}" width="{w}" height="{h}" rx="7"/>\n    '
        + text_lines(cx, cy, lines)
        + "\n  </g>"
    )


def edge(points, label=None, label_at=None, anchor="start", arrow=True):
    d = "M " + " L ".join(f"{x:.0f},{y:.0f}" for x, y in points)
    s = f'<path class="edge{"" if arrow else " bare"}" d="{d}"/>'
    if label:
        lx, ly = label_at
        s += f'\n  <text x="{lx}" y="{ly}" class="lab" style="text-anchor: {anchor}">{math(label)}</text>'
    return s


def label(x, y, s, anchor="start", cls="lab"):
    body = caption(s) if cls == "cap" else math(s)
    return f'<text x="{x}" y="{y}" class="{cls}" style="text-anchor: {anchor}">{body}</text>'


# Columns: the geometry, then the four stages.  Rows: the evaluation, the data, the
# backward pass, the eight numbers, and the contraction.
CG, CS, CI, CM, CN = 95, 341, 553, 761, 989  # centers
WG, WS, WI, WM, WN = 170, 218, 174, 210, 214  # widths
RB, RD, RC, RE, RF = 104, 213, 322, 448, 574  # centers of the rows
H = 84  # height of the boxes of the evaluation and of the backward pass
RAIL = 22  # height of the line that takes g to ι and μ
left = lambda c, w: c - w / 2
right = lambda c, w: c + w / 2

parts = []

# The geometry: forward mode in all the parameters.
parts.append(
    box(
        CG,
        RB,
        WG,
        H + 24,
        "fwd",
        [
            ("cap", "geometry of the pixel,"),
            ("cap", "with its P derivatives"),
            ("math", "_{s}Y^{ℓm}_j,   g^n_j"),
        ],
    )
)

# The evaluation.
parts.append(
    box(
        CS,
        RB,
        WS,
        H,
        "plain heavy",
        [("cap", "synthesis, σ"), ("math", "σ^{kh}_j = d^{kh}_{ℓm} _{s_k}Y^{ℓm}_j")],
    )
)
parts.append(box(CI, RB, WI, H, "plain", [("cap", "interpolation, ι"), ("math", "ι^{ki}_j = S_j«»^i«»_h σ^{kh}_j")]))
parts.append(
    box(CM, RB, WM, H, "plain", [("cap", "mixing, μ"), ("math", "μ^{ki}_j = M_j«»^k«»_{k′} ι^{k′i}_j + q^k_j")])
)
parts.append(
    box(
        CN,
        RB,
        WN,
        H,
        "plain",
        [("cap", "norm, ν"), ("math", "ν_j = Σ_{ik} W_{jki} |μ^{ki}_j − c^{ki}_j|^2")],
    )
)

# The data, which the synthesis and its backward step both use.
parts.append(box(CS, RD, 150, 50, "io", [("math", "«data»  d^{kh}_{ℓm}")]))

# The backward pass.
parts.append(
    box(
        CN,
        RC,
        WN,
        H,
        "rev",
        [("cap", "back through ν"), ("math", "∂ν_j/∂μ^{ki}_j ="), ("math", "W_{jki} (μ̄^{ki}_j − c̄^{ki}_j)")],
    )
)
parts.append(
    box(
        CM,
        RC,
        WM,
        H,
        "rev",
        [("cap", "back through μ"), ("math", "∂ν_j/∂ι^{ki}_j ="), ("math", "(∂ν_j/∂μ^{k′i}_j) M_j«»^{k′}«»_k")],
    )
)
parts.append(
    box(
        CI,
        RC,
        WI,
        H,
        "rev",
        [("cap", "back through ι"), ("math", "∂ν_j/∂σ^{kh}_j ="), ("math", "(∂ν_j/∂ι^{ki}_j) S_j«»^i«»_h")],
    )
)
parts.append(
    box(
        CS,
        RC,
        WS,
        H,
        "rev heavy",
        [
            ("cap", "back through σ"),
            ("math", "∂ν_j/∂_{s}Y^{ℓm}_j ="),
            ("math", "Σ_{k:s_k=s} (∂ν_j/∂σ^{kh}_j) d^{kh}_{ℓm}"),
        ],
    )
)

# The eight numbers, and the contraction.
CE, WE = (CI + CM) // 2 + 10, 480
parts.append(
    box(
        CE,
        RE,
        WE,
        78,
        "fwd",
        [
            ("cap", "ι and μ, with derivatives in the eight numbers, give ∂μ/∂g"),
            ("math", "∂ν_j/∂g^n_j = 2 «Re» [(∂ν_j/∂μ^{ki}_j) ∂μ^{ki}_j/∂g^n_j]"),
        ],
    )
)
CF, WF = 568, 1056
parts.append(
    box(
        CF,
        RF,
        WF,
        96,
        "plain",
        [
            ("cap", "contract the covectors with the P vectors, and sum over the pixels"),
            ("math", "∂ν_j/∂θ^p = 2 «Re» [ Σ_s (∂ν_j/∂_{s}Y^{ℓm}_j) ∂_{s}Y^{ℓm}_j/∂θ^p ] + (∂ν_j/∂g^n_j) ∂g^n_j/∂θ^p"),
            ("math", "∂Φ/∂θ^p = (1/𝒩) Σ_j ∂ν_j/∂θ^p"),
        ],
    )
)

h2 = H / 2
hg = (H + 24) / 2
edges = []
# The parameters into the geometry.
edges.append(edge([(CG, 8), (CG, RB - hg)], label="θ^p", label_at=(CG + 9, 26)))
# The harmonics into the synthesis.
edges.append(
    edge([(right(CG, WG), RB), (left(CS, WS), RB)], label="_{s}Y^{ℓm}_j", label_at=(right(CG, WG) + 25, RB - 11), anchor="middle")
)
# The eight numbers into the interpolation and the mixing.
xg = CG + 60
edges.append(edge([(xg, RB - hg), (xg, RAIL), (CM, RAIL), (CM, RB - h2)], label="g^n_j", label_at=(xg + 30, RAIL - 8)))
edges.append(edge([(CI, RAIL), (CI, RB - h2)]))
# Along the evaluation.
edges.append(edge([(right(CS, WS), RB), (left(CI, WI), RB)]))
edges.append(edge([(right(CI, WI), RB), (left(CM, WM), RB)]))
edges.append(edge([(right(CM, WM), RB), (left(CN, WN), RB)]))
# The data into the synthesis and into its backward step.
edges.append(edge([(CS, RD - 25), (CS, RB + h2)]))
edges.append(edge([(CS, RD + 25), (CS, RC - h2)]))
# The residual starts the backward pass.
edges.append(
    edge([(CN, RB + h2), (CN, RC - h2)], label="μ^{ki}_j − c^{ki}_j", label_at=(CN - 10, RD + 5), anchor="end")
)
# Along the backward pass.
edges.append(edge([(left(CN, WN), RC), (right(CM, WM), RC)]))
edges.append(edge([(left(CM, WM), RC), (right(CI, WI), RC)]))
edges.append(edge([(left(CI, WI), RC), (right(CS, WS), RC)]))
# To the eight numbers, and to the contraction.
edges.append(edge([(CN, RC + h2), (CN, RE), (right(CE, WE), RE)]))
edges.append(edge([(CE, RE + 39), (CE, RF - 48)]))
edges.append(edge([(CS, RC + h2), (CS, RF - 48)], label="∂ν_j/∂_{s}Y^{ℓm}_j", label_at=(CS - 10, RE + 5), anchor="end"))
edges.append(edge([(CG - 40, RB + hg), (CG - 40, RF - 48)]))
edges.append(label(CG - 30, RD + 28, "∂_{s}Y^{ℓm}_j/∂θ^p,"))
edges.append(label(CG - 30, RD + 54, "∂g^n_j/∂θ^p"))

legend_y = 660
legend = f"""
  <g class="legend">
    <rect class="box fwd" x="80" y="{legend_y - 12}" width="22" height="16" rx="3"/>
    <text x="110" y="{legend_y + 1}" class="cap" style="text-anchor: start">forward mode: dual numbers</text>
    <rect class="box rev" x="356" y="{legend_y - 12}" width="22" height="16" rx="3"/>
    <text x="386" y="{legend_y + 1}" class="cap" style="text-anchor: start">backward pass, by hand</text>
    <rect class="box plain heavy" x="606" y="{legend_y - 12}" width="22" height="16" rx="3"/>
    <text x="636" y="{legend_y + 1}" class="cap" style="text-anchor: start">a sum over every mode, time, and component</text>
  </g>"""

svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1106 684" width="100%" style="max-width: 1000px;" role="img"
     aria-label="The flow of the objective and its gradient for one pixel">
  <style>
    .grad-flow text {{ fill: currentColor; text-anchor: middle; font-size: 17px; white-space: pre; }}
    .grad-flow text.cap {{ font-size: 14.5px; }}
    .grad-flow text.lab {{ font-size: 16px; }}
    .grad-flow .mi {{ font-family: KaTeX_Math, "Times New Roman", serif; font-style: italic; }}
    .grad-flow .mu {{ font-family: KaTeX_Main, "Times New Roman", serif; }}
    .grad-flow .pr {{ font-family: "Times New Roman", "Liberation Serif", serif; }}
    .grad-flow .box {{ stroke: currentColor; stroke-width: 1.2; fill: none; }}
    .grad-flow .fwd {{ fill: rgba(56, 132, 255, 0.17); }}
    .grad-flow .rev {{ fill: rgba(255, 140, 0, 0.22); }}
    .grad-flow .plain {{ fill: rgba(128, 128, 128, 0.10); }}
    .grad-flow .heavy {{ stroke-width: 3; }}
    .grad-flow .edge {{ stroke: currentColor; stroke-width: 1.3; fill: none; marker-end: url(#grad-flow-arrow); }}
    .grad-flow .edge.bare {{ marker-end: none; }}
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

target = sys.argv[1] if len(sys.argv) > 1 else None
if target and target.endswith(".svg"):
    open(target, "w").write(svg)
else:
    page = target or os.path.join(
        os.path.dirname(__file__), "..", "src", "40-Details", "06-objective_gradient.md"
    )
    text = open(page).read()
    i, j = text.index("<svg "), text.index("</svg>") + len("</svg>")
    open(page, "w").write(text[:i] + svg + text[j:])
