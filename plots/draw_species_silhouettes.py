"""Silhouettes of the eleven species, dark grey on white, one PDF per species.

Each drawing is a hand-written SVG path set on a 600 x 600 canvas, rendered to
PDF (and a PNG preview) with rsvg-convert.
"""
import subprocess
from pathlib import Path

OUT = Path(__file__).parent / "species_silhouettes"
INK = "#2E2E2E"

HEAD = ('<svg xmlns="http://www.w3.org/2000/svg" width="600" height="600" '
        'viewBox="0 0 600 600"><rect width="600" height="600" fill="#ffffff"/>')


def fill(d):
    return f'<path fill="{INK}" d="{d}"/>'


def line(d, w=6):
    return (f'<path fill="none" stroke="{INK}" stroke-width="{w}" '
            f'stroke-linecap="round" d="{d}"/>')


def cut(d):
    return f'<path fill="#ffffff" d="{d}"/>'


def leafy(cx, cy, rx, ry, rot=0):
    """A pointed leaf (two arcs) centred on cx, cy."""
    d = (f"M{cx - rx} {cy} Q{cx} {cy - ry} {cx + rx} {cy} "
         f"Q{cx} {cy + ry} {cx - rx} {cy} Z")
    return f'<path fill="{INK}" transform="rotate({rot} {cx} {cy})" d="{d}"/>'


SPECIES = {}

# --- Chlamydomonas reinhardtii: ovoid cell, two apical flagella, eyespot ----
SPECIES["Chlamydomonas_reinhardtii"] = "".join([
    fill("M300 250 C371 250 415 305 415 378 C415 452 369 505 300 505 "
         "C231 505 185 452 185 378 C185 305 229 250 300 250 Z"),
    line("M283 252 C258 196 222 152 176 122", 9),
    line("M317 252 C342 196 378 152 424 122", 9),
    cut("M300 320 m-26 0 a26 22 0 1 0 52 0 a26 22 0 1 0 -52 0 Z"),
])

# --- Marchantia polymorpha: thallus, gemma cups, umbrella archegoniophore ---
def _archegoniophore(cx, cy, r, stalk_to):
    """Stalked receptacle with digitate rays."""
    import math
    rays = []
    for i in range(9):
        a = math.pi * (0.04 + 0.92 * i / 8)
        x = cx - math.cos(a) * r * 0.92
        y = cy + math.sin(a) * r * 0.26 + 12
        rays.append(f'<ellipse fill="{INK}" cx="{x:.0f}" cy="{y:.0f}" rx="{r*0.10:.0f}" '
                    f'ry="{r*0.24:.0f}" transform="rotate({(i - 4) * 12} {x:.0f} {y:.0f})"/>')
    return "".join([
        line(f"M{cx} {cy} L{cx} {stalk_to}", 10),
        "".join(rays),
        f'<path fill="{INK}" d="M{cx - r} {cy + 8} '
        f'C{cx - r*0.88:.0f} {cy - r*0.52:.0f} {cx + r*0.88:.0f} {cy - r*0.52:.0f} '
        f'{cx + r} {cy + 8} C{cx + r*0.5:.0f} {cy + r*0.22:.0f} '
        f'{cx - r*0.5:.0f} {cy + r*0.22:.0f} {cx - r} {cy + 8} Z"/>',
    ])


SPECIES["Marchantia_polymorpha"] = "".join([
    _archegoniophore(300, 214, 92, 430),
    fill("M300 424 C256 424 204 436 164 454 C124 472 110 496 122 510 "
         "C138 528 198 534 252 526 C276 522 292 514 300 504 "
         "C308 514 324 522 348 526 C402 534 462 528 478 510 "
         "C490 496 476 472 436 454 C396 436 344 424 300 424 Z"),
    cut("M204 476 m-26 0 a26 14 0 1 0 52 0 a26 14 0 1 0 -52 0 Z"),
    line("M172 520 L160 564 M226 530 L220 574 M280 534 L278 578 "
         "M330 532 L336 576 M386 524 L398 566 M434 512 L450 552", 5),
])

# --- Lunularia cruciata: forked thallus with crescent gemma cups ------------
SPECIES["Lunularia_cruciata"] = "".join([
    fill("M300 464 C206 464 132 436 124 382 C116 326 174 288 228 302 "
         "C260 310 284 332 300 360 C316 332 340 310 372 302 "
         "C426 288 484 326 476 382 C468 436 394 464 300 464 Z"),
    cut("M196 372 a28 28 0 0 0 56 0 a28 18 0 0 1 -56 0 Z"),
    line("M186 456 L174 502 M238 466 L232 514 M290 470 L290 520 "
         "M342 468 L350 516 M396 456 L410 500", 5),
])

# --- Physcomitrium patens: leafy gametophyte with seta and capsule ----------
_leaves = "".join(
    leafy(300 + dx, 452 + dy, 78, 30, rot)
    for dx, dy, rot in [(-72, 10, 16), (72, 10, -16), (-60, -26, 34),
                        (60, -26, -34), (-28, 34, 6), (28, 34, -6)])
SPECIES["Physcomitrium_patens"] = "".join([
    line("M300 430 C294 350 300 274 302 214", 8),
    fill("M302 214 C342 214 364 190 364 158 C364 130 340 98 302 94 "
         "C264 98 240 130 240 158 C240 190 262 214 302 214 Z"),
    _leaves,
    line("M266 512 L248 556 M300 518 L300 564 M336 512 L356 556", 5),
])

# --- Selaginella moellendorffii: branched frond of scale leaves -------------
def _branch(x0, y0, x1, y1, n, size):
    out = [line(f"M{x0} {y0} L{x1} {y1}", 7)]
    for i in range(n):
        t = (i + 1) / (n + 1)
        x = x0 + (x1 - x0) * t
        y = y0 + (y1 - y0) * t
        out.append(leafy(x - size, y, size, size * 0.55, -20))
        out.append(leafy(x + size, y, size, size * 0.55, 20))
    return "".join(out)


def _branch(x0, y0, x1, y1, n, size, w=7):
    out = [line(f"M{x0} {y0} L{x1} {y1}", w)]
    for i in range(n):
        t = (i + 1) / (n + 1)
        x = x0 + (x1 - x0) * t
        y = y0 + (y1 - y0) * t
        out.append(leafy(x - size * 0.95, y, size, size * 0.62, -24))
        out.append(leafy(x + size * 0.95, y, size, size * 0.62, 24))
    return "".join(out)


SPECIES["Selaginella_moellendorffii"] = "".join([
    line("M300 556 C302 490 300 430 300 300", 9),
    _branch(300, 306, 300, 132, 7, 19),
    _branch(300, 344, 186, 210, 6, 17),
    _branch(300, 344, 414, 210, 6, 17),
    _branch(300, 410, 156, 322, 5, 16),
    _branch(300, 410, 444, 322, 5, 16),
    _branch(300, 470, 172, 436, 4, 14),
    _branch(300, 470, 428, 436, 4, 14),
])

# --- Azolla filiculoides: overlapping bilobed leaves on a branching axis ----
def _azolla(x0, y0, x1, y1, n, r):
    out = [line(f"M{x0} {y0} L{x1} {y1}", 5)]
    for i in range(n):
        t = (i + 0.5) / n
        x = x0 + (x1 - x0) * t
        y = y0 + (y1 - y0) * t
        out.append(f'<ellipse fill="{INK}" cx="{x - r*1.1:.0f}" cy="{y:.0f}" '
                   f'rx="{r}" ry="{r*0.74:.0f}" transform="rotate(-24 {x:.0f} {y:.0f})"/>')
        out.append(f'<ellipse fill="{INK}" cx="{x + r*1.1:.0f}" cy="{y:.0f}" '
                   f'rx="{r}" ry="{r*0.74:.0f}" transform="rotate(24 {x:.0f} {y:.0f})"/>')
    return "".join(out)


def _azolla(x0, y0, x1, y1, n, r, w=5):
    out = [line(f"M{x0} {y0} L{x1} {y1}", w)]
    for i in range(n):
        t = (i + 0.5) / n
        x = x0 + (x1 - x0) * t
        y = y0 + (y1 - y0) * t
        for sign in (-1, 1):
            out.append(
                f'<ellipse fill="{INK}" cx="{x + sign*r*1.05:.0f}" cy="{y:.0f}" '
                f'rx="{r}" ry="{r*0.7:.0f}" '
                f'transform="rotate({sign*26} {x:.0f} {y:.0f})"/>')
    return "".join(out)


SPECIES["Azolla_filiculoides"] = "".join([
    _azolla(300, 468, 300, 168, 11, 17),
    _azolla(296, 400, 178, 300, 7, 14),
    _azolla(304, 362, 422, 262, 7, 14),
    _azolla(296, 452, 190, 392, 5, 13),
    _azolla(304, 424, 412, 372, 5, 13),
    line("M300 470 C294 508 288 536 280 562 M300 470 C310 506 318 534 328 558", 5),
])

# --- Salvinia cucullata: floating cupped leaves and submerged root-leaf -----
def _cup(cx, cy, w, h, rot):
    """Hooded floating leaf: a thick crescent that opens upwards."""
    d = (f"M{cx - w} {cy} C{cx - w} {cy + h} {cx + w} {cy + h} {cx + w} {cy} "
         f"C{cx + w*0.6:.0f} {cy + h*0.26:.0f} {cx - w*0.6:.0f} {cy + h*0.26:.0f} "
         f"{cx - w} {cy} Z")
    return f'<path fill="{INK}" transform="rotate({rot} {cx} {cy})" d="{d}"/>'


SPECIES["Salvinia_cucullata"] = "".join([
    fill("M104 350 C206 340 394 340 496 350 C394 362 206 362 104 350 Z"),
    _cup(172, 282, 74, 104, -10),
    _cup(428, 282, 74, 104, 10),
    _cup(300, 250, 82, 112, 0),
    line("M176 330 L176 306 M424 330 L424 306 M300 330 L300 286", 7),
    line("M300 354 C297 412 293 464 288 514", 6),
    line("M297 380 C262 394 238 418 224 450 M295 416 C258 432 234 462 222 496 "
         "M292 452 C262 474 244 502 236 536", 5),
    line("M303 380 C338 394 362 418 376 450 M305 416 C342 432 366 462 378 496 "
         "M308 452 C338 474 356 502 364 536", 5),
])

# --- Picea abies: conical spruce -------------------------------------------
SPECIES["Picea_abies"] = "".join([
    fill("M300 60 L330 128 L312 124 L352 200 L330 196 L376 280 L348 274 "
         "L400 364 L364 356 L424 456 L322 442 L322 544 L278 544 L278 442 "
         "L176 456 L236 356 L200 364 L252 274 L224 280 L270 196 L248 200 "
         "L288 124 L270 128 Z"),
])

# --- Arabidopsis thaliana: rosette, bolting stem, siliques -----------------
def _siliques(pts):
    out = []
    for x, y, rot in pts:
        sign = 1 if x > 300 else -1
        out.append(line(f"M{300 + sign*3} {y + 12} L{x - sign*22} {y + 7}", 4))
        out.append(f'<ellipse fill="{INK}" cx="{x}" cy="{y}" rx="27" ry="7" '
                   f'transform="rotate({rot} {x} {y})"/>')
    return "".join(out)


_rosette = "".join(
    leafy(300 + dx, 486 + dy, 80, 28, rot)
    for dx, dy, rot in [(-76, 6, 16), (76, 6, -16), (-62, -26, 34),
                        (62, -26, -34), (-28, 30, 6), (28, 30, -6)])
SPECIES["Arabidopsis_thaliana"] = "".join([
    line("M300 470 C296 380 302 300 300 132", 10),
    _siliques([(350, 410, -40), (250, 374, 40), (354, 334, -40),
               (248, 296, 40), (348, 256, -40), (254, 220, 40),
               (340, 188, -46)]),
    "".join(leafy(300 + dx, 120 + dy, 18, 10, rot) for dx, dy, rot in
            [(-24, -6, 34), (24, -6, -34), (0, -24, 0)]),
    _rosette,
    line("M280 528 L266 570 M300 532 L300 574 M322 528 L338 570", 5),
])

# --- Solanum lycopersicum: vine with compound leaves and a fruit ------------
def _pinna(x, y, rot):
    return "".join([
        line(f"M{x} {y} L{x + 78} {y - 26}", 5),
        leafy(x + 20, y - 18, 22, 13, rot - 18),
        leafy(x + 46, y - 26, 20, 12, rot - 18),
        leafy(x + 22, y + 6, 22, 13, rot + 18),
        leafy(x + 48, y - 2, 20, 12, rot + 18),
        leafy(x + 80, y - 26, 18, 11, rot),
    ])


def _pinna(x, y, flip=False):
    out = [line(f"M{x} {y} L{x + 104} {y - 34}", 6)]
    for dx, dy, rot in [(28, -22, -26), (62, -34, -26),
                        (30, 8, 22), (64, -4, 22)]:
        out.append(leafy(x + dx, y + dy, 30, 17, rot))
    out.append(leafy(x + 106, y - 34, 26, 15, -22))
    body = "".join(out)
    if flip:
        return f'<g transform="scale(-1 1) translate(-600 0)">{body}</g>'
    return body


SPECIES["Solanum_lycopersicum"] = "".join([
    line("M286 566 C282 466 302 374 296 282 C292 214 314 162 346 120", 11),
    _pinna(296, 486),
    _pinna(304, 400, flip=True),
    _pinna(298, 310),
    _pinna(312, 220, flip=True),
    line("M294 350 C254 346 226 362 212 392", 6),
    f'<circle fill="{INK}" cx="206" cy="418" r="34"/>',
    f'<circle fill="{INK}" cx="258" cy="446" r="26"/>',
])

# --- Persea americana: canopy tree with a fruit -----------------------------
SPECIES["Persea_americana"] = "".join([
    fill("M286 556 L286 372 C266 344 238 336 214 342 C228 324 228 302 218 288 "
         "C242 300 266 300 286 290 L286 240 L314 240 L314 320 "
         "C332 306 358 302 378 310 C362 324 354 344 356 362 "
         "C338 350 324 354 314 370 L314 556 Z"),
    fill("M300 84 C398 84 480 146 480 214 C480 286 402 330 300 330 "
         "C198 330 120 286 120 214 C120 146 202 84 300 84 Z"),
    line("M404 324 C406 340 404 352 402 362", 5),
    fill("M402 360 C418 360 428 378 428 396 C440 408 446 424 446 440 "
         "C446 470 426 490 402 490 C378 490 358 470 358 440 "
         "C358 424 364 408 376 396 C376 378 386 360 402 360 Z"),
])

OUT.mkdir(exist_ok=True)
for name, body in SPECIES.items():
    svg = OUT / f"{name}.svg"
    svg.write_text(HEAD + body + "</svg>")
    subprocess.run(["rsvg-convert", "-f", "pdf", "-o", str(OUT / f"{name}.pdf"),
                    str(svg)], check=True)
    subprocess.run(["rsvg-convert", "-w", "600", "-h", "600", "-o",
                    str(OUT / f"{name}.png"), str(svg)], check=True)
print(f"{len(SPECIES)} species written to {OUT}")
