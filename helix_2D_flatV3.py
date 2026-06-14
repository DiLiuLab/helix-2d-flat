#!/usr/bin/env python3
"""
helix_2D_flatV3.py

Draw a flat 2D nucleic-acid duplex schematic as an SVG. The input is one
DNA/RNA strand written 5' to 3'. The second strand is computed as the
antiparallel complement. The drawing uses phosphate circles, pentose
pentagons, rectangular bases, and explicit line-segment bonds in a
NUCPLOT/PDBsum-like flat style.

Inputs:
  - A sequence on the command line, or a text/FASTA file with --input.
  - Optional drawing parameters such as chain colors, strand direction,
    element sizes, and bond lengths.

Output:
  - An SVG file. If no output path is supplied, the script derives one by
    adding _helix2d before the extension of the input file, or writes
    helix_2D_flatV3_out.svg for direct sequence input.

Example:
  python helix_2D_flatV3.py CGCGCGCGCGCG -o duplex.svg --left-color '#d842c2' --right-color '#f3c0e8'

GUI:
  Run with no arguments, or with --gui, to open the Tkinter GUI.
"""

from __future__ import annotations

import argparse
import html
import math
import os
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


SCRIPT_VERSION = "2026-06-13e-V3"
APP_NAME = "Helix 2D Flat"


ALLOWED_BASES = set("ACGTURYSWKMBDHVN")
IGNORE_CHARS = set(" \t\r\n-_.0123456789'`\"/\\")

DEFAULT_SUGAR_RADIUS = 27.0
DEFAULT_PAIR_GAP = 1.5 * 2.0 * DEFAULT_SUGAR_RADIUS * math.sin(math.pi / 5.0)


DNA_COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A", "U": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "D": "H", "H": "D", "V": "B", "N": "N",
}
RNA_COMPLEMENT = dict(DNA_COMPLEMENT)
RNA_COMPLEMENT["A"] = "U"
RNA_COMPLEMENT["T"] = "A"


def resource_path(relative_path: str) -> str:
    """Return a bundled or source-tree path for an application resource."""
    bundle_root = getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(bundle_root, relative_path)


@dataclass
class DrawingOptions:
    left_color: str = "#d842c2"
    right_color: str = "#f3c0e8"
    outline_color: str = "#111111"
    text_color: str = "#111111"
    oxygen_color: str = "#e83e52"
    background: str = "#ffffff"
    transparent: bool = False
    input_chain: str = "left"
    left_direction: str = "up-to-bottom"
    nucleic_acid: str = "auto"
    show_numbers: bool = True
    terminal_labels: bool = True
    terminal_style: str = "phosphate"  # phosphate, oxygen, none
    show_base_pair_lines: bool = False
    base_pair_line_color: str = "#111111"
    row_spacing: float = 104.0
    purine_width: float = 104.0
    pyrimidine_width: float = 86.0
    ambiguous_width: float = 94.0
    base_height: float = 42.0
    sugar_radius: float = DEFAULT_SUGAR_RADIUS
    phosphate_radius: float = 27.0
    pair_gap: float = DEFAULT_PAIR_GAP
    base_bond_length: float = 30.0
    phosphate_kink_length: float = 30.0
    margin_x: float = 44.0
    margin_top: float = 150.0
    margin_bottom: float = 112.0
    left_sugar_x: float = 170.0
    left_phosphate_x: float = 0.0
    font_family: str = "Arial, Helvetica, sans-serif"
    font_size: float = 28.0
    base_label_padding: float = 10.0
    rotate_bottom_up_labels: bool = False
    terminal_font_size: float = 30.0
    stroke_width: float = 3.0
    title: str = "Flat 2D nucleic acid helix V3"


@dataclass
class ResidueRow:
    y: float
    left_base: str
    left_number: int
    right_base: str
    right_number: int


def clean_sequence_text(text: str) -> str:
    """Return an uppercase sequence after removing FASTA headers and spacing."""
    kept: List[str] = []
    bad_letters: List[str] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith(">") or line.startswith(";"):
            continue
        for ch in line.upper():
            if ch in ALLOWED_BASES:
                kept.append(ch)
            elif ch in IGNORE_CHARS:
                continue
            elif ch.isalpha():
                bad_letters.append(ch)
            else:
                continue

    if bad_letters:
        unique = "".join(sorted(set(bad_letters)))
        raise ValueError(
            "Unsupported sequence letter(s): {}. Allowed letters are {}.".format(
                unique, "".join(sorted(ALLOWED_BASES))
            )
        )

    seq = "".join(kept)
    if not seq:
        raise ValueError("No sequence letters were found.")
    return seq


def load_sequence_from_args(args: argparse.Namespace) -> Tuple[str, Optional[str]]:
    """Load sequence from --seq, --input, positional sequence, or file path."""
    sources = [bool(args.sequence_opt), bool(args.input_file), bool(args.sequence)]
    if sum(sources) == 0:
        raise ValueError("Please provide a sequence, --seq, --input, or run --gui.")
    if sum(sources) > 1:
        raise ValueError("Provide only one of positional sequence, --seq, or --input.")

    if args.input_file:
        with open(args.input_file, "r", encoding="utf-8") as handle:
            return clean_sequence_text(handle.read()), args.input_file

    if args.sequence_opt:
        return clean_sequence_text(args.sequence_opt), None

    # Positional argument can be either a literal sequence or an existing file.
    if args.sequence and os.path.exists(args.sequence) and os.path.isfile(args.sequence):
        with open(args.sequence, "r", encoding="utf-8") as handle:
            return clean_sequence_text(handle.read()), args.sequence
    return clean_sequence_text(args.sequence or ""), None


def detect_nucleic_acid(seq: str, requested: str) -> str:
    """Return dna or rna, respecting the user request when supplied."""
    requested = requested.lower()
    if requested in ("dna", "rna"):
        return requested
    if requested != "auto":
        raise ValueError("nucleic_acid must be auto, dna, or rna.")
    has_t = "T" in seq
    has_u = "U" in seq
    if has_t and has_u:
        raise ValueError("The sequence contains both T and U. Use --nucleic-acid dna or rna to choose a complement alphabet.")
    if has_u:
        return "rna"
    return "dna"


def complement_base(base: str, nucleic_acid: str) -> str:
    mapping = RNA_COMPLEMENT if nucleic_acid == "rna" else DNA_COMPLEMENT
    return mapping.get(base.upper(), "N")


def normalize_left_direction(direction: str) -> str:
    """Normalize user spellings of the left strand's 5'->3' direction."""
    value = (direction or "").strip().lower().replace("_", "-").replace(" ", "-")
    aliases = {
        "bottom-to-up": "bottom-to-up",
        "bottom-to-top": "bottom-to-up",
        "bottom-up": "bottom-to-up",
        "up": "bottom-to-up",
        "upward": "bottom-to-up",
        "up-to-bottom": "up-to-bottom",
        "top-to-bottom": "up-to-bottom",
        "top-down": "up-to-bottom",
        "down": "up-to-bottom",
        "downward": "up-to-bottom",
    }
    if value not in aliases:
        raise ValueError("left_direction must be bottom-to-up or up-to-bottom.")
    return aliases[value]


def opposite_direction(direction: str) -> str:
    direction = normalize_left_direction(direction)
    if direction == "bottom-to-up":
        return "up-to-bottom"
    return "bottom-to-up"


def direction_y_sign(direction: str) -> int:
    """Return +1 when 5'->3' proceeds downward in SVG coordinates, else -1."""
    direction = normalize_left_direction(direction)
    if direction == "up-to-bottom":
        return 1
    return -1


def position_for_row(row_index: int, direction: str, n: int) -> int:
    """Return the 5'->3' residue number for a top-to-bottom row index."""
    if direction_y_sign(direction) == 1:
        return row_index + 1
    return n - row_index


def row_index_for_position(position: int, direction: str, n: int) -> int:
    """Return the top-to-bottom row index for a 5'->3' residue number."""
    if position < 1 or position > n:
        raise ValueError("Residue position {} is outside 1..{}.".format(position, n))
    if direction_y_sign(direction) == 1:
        return position - 1
    return n - position


def side_direction(side: str, options: DrawingOptions) -> str:
    """Return the 5'->3' direction for a drawn side."""
    left = normalize_left_direction(options.left_direction)
    if side == "left":
        return left
    if side == "right":
        return opposite_direction(left)
    raise ValueError("side must be left or right.")


def make_rows(seq: str, options: DrawingOptions) -> List[ResidueRow]:
    """Create top-to-bottom residue rows for the drawing."""
    n = len(seq)
    na = detect_nucleic_acid(seq, options.nucleic_acid)
    left_dir = side_direction("left", options)
    right_dir = side_direction("right", options)
    rows: List[ResidueRow] = []

    for row_index in range(n):
        y = options.margin_top + row_index * options.row_spacing
        left_number = position_for_row(row_index, left_dir, n)
        right_number = position_for_row(row_index, right_dir, n)

        if options.input_chain == "left":
            left_base = seq[left_number - 1]
            right_base = complement_base(left_base, na)
        else:
            right_base = seq[right_number - 1]
            left_base = complement_base(right_base, na)

        rows.append(ResidueRow(y, left_base, left_number, right_base, right_number))
    return rows


def base_box_width(base: str, options: DrawingOptions, label: str) -> float:
    """Return base rectangle width, keeping purines wider than pyrimidines."""
    base = base.upper()
    if base in ("A", "G", "R"):
        width = options.purine_width
    elif base in ("C", "T", "U", "Y"):
        width = options.pyrimidine_width
    else:
        width = options.ambiguous_width

    estimated_label_width = len(label) * options.font_size * 0.58 + 20.0
    return max(width, estimated_label_width)


def derive_output_path(output: Optional[str], input_file: Optional[str]) -> str:
    """Choose an SVG output path, adding a suffix before the extension when possible."""
    if output:
        root, ext = os.path.splitext(output)
        if not ext:
            return output + ".svg"
        if ext.lower() != ".svg":
            return root + ".svg"
        return output

    if input_file:
        directory, filename = os.path.split(input_file)
        root, _ext = os.path.splitext(filename)
        return os.path.join(directory or ".", root + "_helix2d.svg")

    return "helix_2D_flatV3_out.svg"


def points_to_string(points: Iterable[Tuple[float, float]]) -> str:
    return " ".join("{:.2f},{:.2f}".format(x, y) for x, y in points)


def regular_pentagon(cx: float, cy: float, r: float, point_angle_degrees: float) -> List[Tuple[float, float]]:
    """Return a pentagon with one vertex at point_angle_degrees."""
    angle0 = math.radians(point_angle_degrees)
    return [
        (cx + r * math.cos(angle0 + 2.0 * math.pi * k / 5.0),
         cy + r * math.sin(angle0 + 2.0 * math.pi * k / 5.0))
        for k in range(5)
    ]


def svg_line(x1: float, y1: float, x2: float, y2: float, color: str, width: float, extra: str = "") -> str:
    return (
        '<line x1="{:.2f}" y1="{:.2f}" x2="{:.2f}" y2="{:.2f}" '
        'stroke="{}" stroke-width="{:.2f}" stroke-linecap="round" {}/>'
    ).format(x1, y1, x2, y2, html.escape(color), width, extra)


def svg_text(x: float, y: float, text: str, size: float, color: str, family: str,
             weight: str = "700", anchor: str = "middle", rotate: float = 0.0) -> str:
    transform = ""
    if abs(rotate) > 1e-9:
        transform = ' transform="rotate({:.2f} {:.2f} {:.2f})"'.format(rotate, x, y)
    return (
        '<text x="{:.2f}" y="{:.2f}" text-anchor="{}" dominant-baseline="middle" '
        'font-family="{}" font-size="{:.2f}" font-weight="{}" fill="{}"{}>{}</text>'
    ).format(
        x, y, html.escape(anchor), html.escape(family), size, html.escape(weight),
        html.escape(color), transform, html.escape(text)
    )


def svg_rect(x: float, y: float, w: float, h: float, fill: str, stroke: str, stroke_width: float) -> str:
    return (
        '<rect x="{:.2f}" y="{:.2f}" width="{:.2f}" height="{:.2f}" '
        'fill="{}" stroke="{}" stroke-width="{:.2f}"/>'
    ).format(x, y, w, h, html.escape(fill), html.escape(stroke), stroke_width)


def svg_circle(cx: float, cy: float, r: float, fill: str, stroke: str, stroke_width: float) -> str:
    return (
        '<circle cx="{:.2f}" cy="{:.2f}" r="{:.2f}" fill="{}" '
        'stroke="{}" stroke-width="{:.2f}"/>'
    ).format(cx, cy, r, html.escape(fill), html.escape(stroke), stroke_width)


def svg_polygon(points: List[Tuple[float, float]], fill: str, stroke: str, stroke_width: float) -> str:
    return '<polygon points="{}" fill="{}" stroke="{}" stroke-width="{:.2f}"/>'.format(
        points_to_string(points), html.escape(fill), html.escape(stroke), stroke_width
    )


def side_outer_sign(side: str) -> int:
    if side == "left":
        return -1
    if side == "right":
        return 1
    raise ValueError("side must be left or right.")


def side_fill_color(side: str, options: DrawingOptions) -> str:
    return options.left_color if side == "left" else options.right_color


def sugar_inward_x_offset(options: DrawingOptions) -> float:
    """Horizontal offset from sugar center to C1' when C2' points down."""
    return options.sugar_radius * math.cos(math.radians(18.0))


def side_sugar_x(side: str, options: DrawingOptions) -> float:
    """Return sugar x-position, using pair_gap to keep fixed sugar-base bond lengths."""
    if side == "left":
        return options.left_sugar_x

    canonical_pair_width = options.purine_width + options.pyrimidine_width
    inward_vertex_distance = 2.0 * options.base_bond_length + canonical_pair_width + options.pair_gap
    return options.left_sugar_x + 2.0 * sugar_inward_x_offset(options) + inward_vertex_distance


def side_phosphate_x(side: str, options: DrawingOptions) -> float:
    """Return phosphate x-position.

    When left_phosphate_x is positive, it is treated as a legacy explicit
    left-side x-coordinate and mirrored onto the right side. The default value
    is 0, which enables automatic placement: the phosphate center is placed
    just outside the average of the C3'/C4' outer sugar vertices, offset by
    phosphate_kink_length. This keeps the phosphate bonds short, balanced, and
    visually close to the sugar-pentagon edge length.
    """
    if side not in ("left", "right"):
        raise ValueError("side must be left or right.")

    if options.left_phosphate_x > 0.0:
        offset = options.left_sugar_x - options.left_phosphate_x
        if side == "left":
            return options.left_phosphate_x
        return side_sugar_x("right", options) + offset

    # Average horizontal distance from sugar center to the two outer vertices
    # C3' and C4' in the oriented pentagon. The phosphate radius is included
    # because bonds stop at the circle edge rather than at the hidden center.
    outer_avg = options.sugar_radius * (math.cos(math.radians(18.0)) + math.cos(math.radians(54.0))) / 2.0
    center_offset = outer_avg + options.phosphate_kink_length + options.phosphate_radius
    if side == "left":
        return side_sugar_x("left", options) - center_offset
    return side_sugar_x("right", options) + center_offset


def sugar_center_y(row_y: float, options: DrawingOptions) -> float:
    """Return the pentagon center y so C1' lies on the base-pair row."""
    return row_y - options.sugar_radius * math.sin(math.radians(18.0))


def polar_point(cx: float, cy: float, r: float, angle_degrees: float) -> Tuple[float, float]:
    angle = math.radians(angle_degrees)
    return cx + r * math.cos(angle), cy + r * math.sin(angle)


def point_on_circle_toward(center: Tuple[float, float], toward: Tuple[float, float], radius: float) -> Tuple[float, float]:
    """Return the point on a circle nearest to another point.

    Bonds to phosphates are drawn to the circle edge, not to the hidden center,
    so the visible line lengths stay comparable to the sugar-ring edge length.
    """
    if radius <= 0.0:
        return center
    dx = toward[0] - center[0]
    dy = toward[1] - center[1]
    dist = math.hypot(dx, dy)
    if dist < 1e-9:
        return center
    return center[0] + radius * dx / dist, center[1] + radius * dy / dist


def phosphate_bond_endpoint(center: Tuple[float, float], toward: Tuple[float, float],
                            position: int, options: DrawingOptions) -> Tuple[float, float]:
    """Endpoint for a bond drawn from a sugar/kink point to a phosphate.

    Internal phosphates and terminal 5' phosphates are circles, so lines stop at
    the circle circumference. The legacy terminal oxygen style has no circle, so
    that special case keeps the old center endpoint.
    """
    radius = options.phosphate_radius
    if position == 1 and options.terminal_style != "phosphate":
        radius = 0.0
    return point_on_circle_toward(center, toward, radius)


def sugar_atom_points(side: str, y: float, options: DrawingOptions) -> Dict[str, Tuple[float, float]]:
    """Return named sugar-ring vertices.

    V3 draws the two strands as C2-rotation partners. The left sugar keeps the
    requested orientation: C4'-O4' is horizontal and C2' is exactly at the
    bottom. The right sugar is the 180-degree rotation of that shape, so its
    C4'-O4' edge is also horizontal but appears at the bottom, and its C2'
    vertex appears at the top. This makes the two backbones look like true
    antiparallel C2-related strands rather than unrelated mirror drawings.

    C1' is kept on the base-pair row on both sides, so sugar-base bonds remain
    horizontal.
    """
    if side not in ("left", "right"):
        raise ValueError("side must be left or right.")

    cx = side_sugar_x(side, options)
    cy = y
    r = options.sugar_radius

    # Offsets for a regular pentagon whose left-strand C4'-O4' edge is
    # horizontal and whose C2' vertex points downward. These are expressed
    # relative to the base-pair row, with C1' at y = row_y.
    c1_dx = r * math.cos(math.radians(18.0))
    c4_dx = r * math.cos(math.radians(54.0))
    c2_dy = r * (1.0 - math.sin(math.radians(18.0)))
    c4_dy = r * (math.sin(math.radians(54.0)) + math.sin(math.radians(18.0)))

    left_offsets = {
        "C1": (c1_dx, 0.0),
        "C2": (0.0, c2_dy),
        "C3": (-c1_dx, 0.0),
        "C4": (-c4_dx, -c4_dy),
        "O4": (c4_dx, -c4_dy),
    }

    if side == "left":
        return {name: (cx + dx, cy + dy) for name, (dx, dy) in left_offsets.items()}

    # 180-degree rotation of the left sugar local geometry.
    return {name: (cx - dx, cy - dy) for name, (dx, dy) in left_offsets.items()}


def sugar_points(side: str, y: float, options: DrawingOptions) -> List[Tuple[float, float]]:
    """Return sugar pentagon points in ring order: C1'-C2'-C3'-C4'-O4'."""
    atoms = sugar_atom_points(side, y, options)
    return [atoms["C1"], atoms["C2"], atoms["C3"], atoms["C4"], atoms["O4"]]


def sugar_anchors(side: str, y: float, options: DrawingOptions) -> Dict[str, Tuple[float, float]]:
    """Return the key sugar vertices used for bonds.

    The generic top_outer/bottom_outer names are historical names used by the
    backbone layout code. In V3, they are assigned so the right-side backbone is
    the C2-rotated counterpart of the left-side backbone.
    """
    atoms = sugar_atom_points(side, y, options)
    if side == "left":
        top_outer = atoms["C4"]
        bottom_outer = atoms["C3"]
    elif side == "right":
        top_outer = atoms["C3"]
        bottom_outer = atoms["C4"]
    else:
        raise ValueError("side must be left or right.")

    return {
        "inward": atoms["C1"],
        "bottom_outer": bottom_outer,
        "top_outer": top_outer,
        "C1": atoms["C1"],
        "C2": atoms["C2"],
        "C3": atoms["C3"],
        "C4": atoms["C4"],
        "O4": atoms["O4"],
    }


def sugar_anchor(side: str, y: float, anchor_name: str, options: DrawingOptions) -> Tuple[float, float]:
    anchors = sugar_anchors(side, y, options)
    if anchor_name not in anchors:
        raise ValueError("Unknown sugar anchor: {}".format(anchor_name))
    return anchors[anchor_name]


def five_prime_anchor_name(direction: str) -> str:
    if direction_y_sign(direction) == 1:
        return "top_outer"
    return "bottom_outer"


def three_prime_anchor_name(direction: str) -> str:
    if direction_y_sign(direction) == 1:
        return "bottom_outer"
    return "top_outer"


def circle_intersections_equal_radius(c1: Tuple[float, float], c2: Tuple[float, float],
                                      radius: float) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """Return two intersection points of equal-radius circles.

    If the requested radius is too small, it is expanded just enough to produce
    a stable solution. This keeps the two visible phosphate bonds as equal as
    possible even after users change row spacing or element sizes.
    """
    dx = c2[0] - c1[0]
    dy = c2[1] - c1[1]
    d = math.hypot(dx, dy)
    if d < 1e-9:
        return (c1[0], c1[1] + radius), (c1[0], c1[1] - radius)

    use_radius = max(radius, d / 2.0 + 0.01)
    mid_x = (c1[0] + c2[0]) / 2.0
    mid_y = (c1[1] + c2[1]) / 2.0
    half = d / 2.0
    height = math.sqrt(max(0.0, use_radius * use_radius - half * half))

    ux = dx / d
    uy = dy / d
    px = -uy
    py = ux

    return (mid_x + px * height, mid_y + py * height), (mid_x - px * height, mid_y - py * height)


def choose_outer_point(side: str, p1: Tuple[float, float], p2: Tuple[float, float]) -> Tuple[float, float]:
    """Choose the point farther away from the base-pair center line."""
    outer = side_outer_sign(side)
    return p1 if outer * p1[0] >= outer * p2[0] else p2


def phosphate_center_for_residue(side: str, position: int, rows: List[ResidueRow], options: DrawingOptions) -> Tuple[float, float]:
    """Return the phosphate center at the 5' side of a residue position.

    The center is placed geometrically so the visible bond lengths from the
    sugar/kink point to the phosphate circle edge and from the previous sugar
    vertex to the phosphate circle edge are nearly equal to the short
    sugar-to-kink segment.
    """
    n = len(rows)
    direction = side_direction(side, options)
    row_index = row_index_for_position(position, direction, n)
    y = rows[row_index].y
    outer = side_outer_sign(side)

    five_anchor = sugar_anchor(side, y, five_prime_anchor_name(direction), options)
    kink = (five_anchor[0] + outer * options.phosphate_kink_length, five_anchor[1])

    circle_radius = options.phosphate_radius
    if position == 1 and options.terminal_style != "phosphate":
        circle_radius = 0.0

    target_center_distance = options.phosphate_kink_length + circle_radius

    if position == 1:
        # The terminal 5' phosphate has only the two-segment connection, so put
        # its center one visible bond length plus one radius away from the kink.
        return (kink[0], kink[1] - direction_y_sign(direction) * target_center_distance)

    previous_row_index = row_index_for_position(position - 1, direction, n)
    previous_y = rows[previous_row_index].y
    three_anchor = sugar_anchor(side, previous_y, three_prime_anchor_name(direction), options)

    p1, p2 = circle_intersections_equal_radius(kink, three_anchor, target_center_distance)
    return choose_outer_point(side, p1, p2)


def row_labels(row: ResidueRow, show_numbers: bool) -> Tuple[str, str]:
    if show_numbers:
        return "{}{}".format(row.left_base, row.left_number), "{}{}".format(row.right_base, row.right_number)
    return row.left_base, row.right_base


def make_base_geometry(row: ResidueRow, options: DrawingOptions) -> Dict[str, float]:
    """Return base rectangles. Sugar-base bonds have fixed length, so pair gaps shift."""
    left_label, right_label = row_labels(row, options.show_numbers)
    left_w = base_box_width(row.left_base, options, left_label)
    right_w = base_box_width(row.right_base, options, right_label)

    left_inward = sugar_anchor("left", row.y, "inward", options)
    right_inward = sugar_anchor("right", row.y, "inward", options)

    left_x = left_inward[0] + options.base_bond_length
    right_x = right_inward[0] - options.base_bond_length - right_w

    return {
        "left_x": left_x,
        "left_w": left_w,
        "right_x": right_x,
        "right_w": right_w,
        "y": row.y - options.base_height / 2.0,
        "h": options.base_height,
    }


def base_pair_bond_count(left_base: str, right_base: str) -> int:
    """Return the schematic H-bond count for canonical Watson-Crick pairs."""
    pair = {left_base.upper(), right_base.upper()}
    if pair == {"G", "C"}:
        return 3
    if pair == {"A", "T"} or pair == {"A", "U"}:
        return 2
    return 1


def dashed_base_pair_offsets(count: int, options: DrawingOptions) -> List[float]:
    """Return y-offsets for multiple dashed base-pair lines.

    The offsets are intentionally more separated than in V3 so the two
    A:T/A:U lines and three G:C lines read as distinct H-bond equivalents.
    """
    spread = max(8.0, min(options.base_height * 0.34, options.font_size * 0.50))
    if count <= 1:
        return [0.0]
    if count == 2:
        return [-0.70 * spread, 0.70 * spread]
    return [-spread, 0.0, spread]


def add_dashed_base_pair_lines(parts: List[str], row: ResidueRow, geom: Dict[str, float], options: DrawingOptions) -> None:
    """Draw 2 dashed lines for A:T/A:U and 3 dashed lines for G:C."""
    count = base_pair_bond_count(row.left_base, row.right_base)
    width = max(1.0, options.stroke_width * 0.75)
    extra = 'stroke-dasharray="5 5"'
    x1 = geom["left_x"] + geom["left_w"]
    x2 = geom["right_x"]
    for offset in dashed_base_pair_offsets(count, options):
        y = row.y + offset
        parts.append(svg_line(x1, y, x2, y, options.base_pair_line_color, width, extra))


def add_backbone_bonds(parts: List[str], rows: List[ResidueRow], options: DrawingOptions) -> None:
    """Add all line-segment bonds behind the nucleotide shapes."""
    sw = options.stroke_width
    color = options.outline_color

    # Glycosidic-style sugar-base bonds. These are fixed length, so the inner
    # gap between bases shifts horizontally depending on which side has a purine.
    for row in rows:
        geom = make_base_geometry(row, options)
        y = row.y
        left_inward = sugar_anchor("left", y, "inward", options)
        right_inward = sugar_anchor("right", y, "inward", options)
        parts.append(svg_line(left_inward[0], left_inward[1], geom["left_x"], left_inward[1], color, sw))
        parts.append(svg_line(geom["right_x"] + geom["right_w"], right_inward[1], right_inward[0], right_inward[1], color, sw))

        if options.show_base_pair_lines:
            add_dashed_base_pair_lines(parts, row, geom, options)

    # Phosphodiester-style backbone bonds. Each phosphate is placed at the 5'
    # side of one residue. It connects to that residue's 4'/5' side with two
    # line segments and, if internal, to the previous residue's 3' side with one
    # line segment.
    n = len(rows)
    for side in ("left", "right"):
        direction = side_direction(side, options)
        outer = side_outer_sign(side)
        for position in range(1, n + 1):
            if position == 1 and options.terminal_style == "none":
                continue

            row_index = row_index_for_position(position, direction, n)
            y = rows[row_index].y
            p_x, p_y = phosphate_center_for_residue(side, position, rows, options)

            five_anchor = sugar_anchor(side, y, five_prime_anchor_name(direction), options)
            phosphate_center = (p_x, p_y)
            kink = (five_anchor[0] + outer * options.phosphate_kink_length, five_anchor[1])
            p_edge_from_kink = phosphate_bond_endpoint(phosphate_center, kink, position, options)
            parts.append(svg_line(five_anchor[0], five_anchor[1], kink[0], kink[1], color, sw))
            parts.append(svg_line(kink[0], kink[1], p_edge_from_kink[0], p_edge_from_kink[1], color, sw))

            if position > 1:
                previous_row_index = row_index_for_position(position - 1, direction, n)
                previous_y = rows[previous_row_index].y
                three_anchor = sugar_anchor(side, previous_y, three_prime_anchor_name(direction), options)
                p_edge_from_three = phosphate_bond_endpoint(phosphate_center, three_anchor, position, options)
                parts.append(svg_line(three_anchor[0], three_anchor[1], p_edge_from_three[0], p_edge_from_three[1], color, sw))


def add_terminal_3prime_stubs(parts: List[str], rows: List[ResidueRow], options: DrawingOptions) -> None:
    """Add short line stubs at the two 3' ends."""
    if not rows:
        return
    n = len(rows)
    sw = options.stroke_width
    color = options.outline_color

    for side in ("left", "right"):
        direction = side_direction(side, options)
        outer = side_outer_sign(side)
        vertical = direction_y_sign(direction)
        row_index = row_index_for_position(n, direction, n)
        y = rows[row_index].y
        start = sugar_anchor(side, y, three_prime_anchor_name(direction), options)
        # Keep the terminal 3' stub close to the same visible bond length.
        end = (start[0] + outer * options.phosphate_kink_length * 0.80,
               start[1] + vertical * options.phosphate_kink_length * 0.60)
        parts.append(svg_line(start[0], start[1], end[0], end[1], color, sw))


def visual_text_anchor_for_rotated_label(side: str, rotate_label: bool) -> str:
    """Return SVG text-anchor that gives the requested visual edge alignment.

    Left-strand base labels are visually left-aligned; right-strand labels are
    visually right-aligned. When a label is rotated 180 degrees, the SVG anchor
    has to be inverted to preserve that visual alignment.
    """
    visual_anchor = "start" if side == "left" else "end"
    if not rotate_label:
        return visual_anchor
    if visual_anchor == "start":
        return "end"
    if visual_anchor == "end":
        return "start"
    return visual_anchor


def base_label_position(side: str, geom: Dict[str, float], row_y: float, options: DrawingOptions) -> Tuple[float, float, str, float]:
    """Return x, y, anchor, and rotation for a base label."""
    pad = max(0.0, options.base_label_padding)
    direction = side_direction(side, options)
    rotate_label = options.rotate_bottom_up_labels and normalize_left_direction(direction) == "bottom-to-up"
    rotate = 180.0 if rotate_label else 0.0

    if side == "left":
        x = geom["left_x"] + pad
    elif side == "right":
        x = geom["right_x"] + geom["right_w"] - pad
    else:
        raise ValueError("side must be left or right.")
    anchor = visual_text_anchor_for_rotated_label(side, rotate_label)
    return x, row_y + 1.0, anchor, rotate


def add_shapes_and_base_labels(parts: List[str], rows: List[ResidueRow], options: DrawingOptions) -> None:
    sw = options.stroke_width
    outline = options.outline_color
    font = options.font_family
    n = len(rows)

    # Phosphate circles at the 5' side of each residue. The terminal 5' group is
    # a normal phosphate circle by default, but can be drawn as a legacy red O or hidden.
    for side in ("left", "right"):
        fill = side_fill_color(side, options)
        for position in range(1, n + 1):
            p_x, p_y = phosphate_center_for_residue(side, position, rows, options)
            if position == 1:
                if options.terminal_style == "phosphate":
                    parts.append(svg_circle(p_x, p_y, options.phosphate_radius, fill, outline, sw))
                elif options.terminal_style == "oxygen":
                    parts.append(svg_text(p_x, p_y, "O", options.terminal_font_size,
                                          options.oxygen_color, font, weight="700"))
                else:
                    continue
            else:
                parts.append(svg_circle(p_x, p_y, options.phosphate_radius, fill, outline, sw))

    for row in rows:
        geom = make_base_geometry(row, options)
        y = row.y
        left_label, right_label = row_labels(row, options.show_numbers)

        # Sugars.
        parts.append(svg_polygon(sugar_points("left", y, options), options.left_color, outline, sw))
        parts.append(svg_polygon(sugar_points("right", y, options), options.right_color, outline, sw))

        # Bases.
        parts.append(svg_rect(geom["left_x"], geom["y"], geom["left_w"], geom["h"],
                              options.left_color, outline, sw))
        parts.append(svg_rect(geom["right_x"], geom["y"], geom["right_w"], geom["h"],
                              options.right_color, outline, sw))

        # Base labels: visually left-aligned on the left strand and
        # right-aligned on the right strand. Optional 180-degree rotation can be
        # applied to labels on strands whose 5'->3' direction is bottom-to-up.
        lx, ly, l_anchor, l_rotate = base_label_position("left", geom, y, options)
        rx, ry, r_anchor, r_rotate = base_label_position("right", geom, y, options)
        parts.append(svg_text(lx, ly, left_label, options.font_size, options.text_color,
                              font, weight="700", anchor=l_anchor, rotate=l_rotate))
        parts.append(svg_text(rx, ry, right_label, options.font_size, options.text_color,
                              font, weight="700", anchor=r_anchor, rotate=r_rotate))


def add_terminal_labels(parts: List[str], rows: List[ResidueRow], options: DrawingOptions) -> None:
    """Add 5' and 3' text labels after the nucleotide shapes."""
    if not rows or not options.terminal_labels:
        return

    n = len(rows)
    font = options.font_family
    size = options.terminal_font_size

    for side in ("left", "right"):
        direction = side_direction(side, options)
        outer = side_outer_sign(side)
        anchor = "end" if side == "left" else "start"

        # 5' label near the terminal phosphate/oxygen position.
        p_x, p_y = phosphate_center_for_residue(side, 1, rows, options)
        label_x = p_x + outer * (options.phosphate_radius + 18.0)
        parts.append(svg_text(label_x, p_y, "5'", size, options.text_color, font, weight="700", anchor=anchor))

        # 3' label near the terminal stub.
        row_index = row_index_for_position(n, direction, n)
        y = rows[row_index].y
        three_anchor = sugar_anchor(side, y, three_prime_anchor_name(direction), options)
        vertical = direction_y_sign(direction)
        stub_end = (three_anchor[0] + outer * options.phosphate_kink_length * 0.80,
                    three_anchor[1] + vertical * options.phosphate_kink_length * 0.60)
        label_x = stub_end[0] + outer * 22.0
        label_y = stub_end[1] + vertical * 1.0
        parts.append(svg_text(label_x, label_y, "3'", size, options.text_color, font, weight="700", anchor=anchor))


def svg_dimensions(rows: List[ResidueRow], options: DrawingOptions) -> Tuple[float, float]:
    if rows:
        content_bottom = rows[-1].y + options.margin_bottom + options.row_spacing / 2.0
    else:
        content_bottom = options.margin_top + options.margin_bottom

    right_edge = side_phosphate_x("right", options) + options.phosphate_radius + options.margin_x
    width = right_edge + options.terminal_font_size * 1.6
    height = content_bottom
    return width, height


def create_helix_svg(seq: str, options: DrawingOptions) -> str:
    """Create the SVG text for the duplex drawing."""
    if not seq:
        raise ValueError("Sequence is empty.")
    if options.input_chain not in ("left", "right"):
        raise ValueError("input_chain must be left or right.")
    options.left_direction = normalize_left_direction(options.left_direction)
    if options.terminal_style not in ("phosphate", "oxygen", "none"):
        raise ValueError("terminal_style must be phosphate, oxygen, or none.")

    rows = make_rows(seq, options)
    width, height = svg_dimensions(rows, options)

    parts: List[str] = []
    parts.append('<?xml version="1.0" encoding="UTF-8"?>')
    parts.append('<svg xmlns="http://www.w3.org/2000/svg" width="{:.0f}" height="{:.0f}" viewBox="0 0 {:.2f} {:.2f}">'.format(
        width, height, width, height
    ))
    parts.append('<title>{}</title>'.format(html.escape(options.title)))
    if not options.transparent:
        parts.append('<rect x="0" y="0" width="100%" height="100%" fill="{}"/>'.format(html.escape(options.background)))

    parts.append('<g id="bonds">')
    add_backbone_bonds(parts, rows, options)
    add_terminal_3prime_stubs(parts, rows, options)
    parts.append('</g>')

    parts.append('<g id="nucleotides">')
    add_shapes_and_base_labels(parts, rows, options)
    parts.append('</g>')

    parts.append('<g id="terminal-labels">')
    add_terminal_labels(parts, rows, options)
    parts.append('</g>')

    parts.append('</svg>')
    return "\n".join(parts) + "\n"


def write_svg_file(seq: str, output_path: str, options: DrawingOptions) -> None:
    svg = create_helix_svg(seq, options)
    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write(svg)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Draw a flat 2D nucleic-acid duplex schematic as SVG.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("sequence", nargs="?", help="Sequence string, or an existing text/FASTA file path.")
    parser.add_argument("--seq", dest="sequence_opt", help="Sequence string, 5' to 3'.")
    parser.add_argument("-i", "--input", dest="input_file", help="Text or FASTA file containing one sequence.")
    parser.add_argument("-o", "--output", help="Output SVG file.")
    parser.add_argument("--gui", action="store_true", help="Open the GUI.")
    parser.add_argument("--version", action="version", version="{} {}".format(APP_NAME, SCRIPT_VERSION))
    parser.add_argument("--input-chain", choices=["left", "right"], default=DrawingOptions.input_chain,
                        help="Which drawn chain is the input 5'->3' strand.")
    parser.add_argument("--left-direction", "--left-strand-direction",
                        choices=["bottom-to-up", "up-to-bottom", "bottom-to-top", "top-to-bottom", "bottom-up", "top-down", "up", "down"],
                        default=DrawingOptions.left_direction,
                        help="Direction of the left strand in the 5'->3' sense. The right strand is antiparallel.")
    parser.add_argument("--nucleic-acid", choices=["auto", "dna", "rna"], default="auto",
                        help="Complement alphabet to use.")

    parser.add_argument("--left-color", default=DrawingOptions.left_color, help="Fill color for the left chain.")
    parser.add_argument("--right-color", default=DrawingOptions.right_color, help="Fill color for the right chain.")
    parser.add_argument("--outline-color", default=DrawingOptions.outline_color, help="Outline/bond color.")
    parser.add_argument("--text-color", default=DrawingOptions.text_color, help="Base and terminal label color.")
    parser.add_argument("--oxygen-color", default=DrawingOptions.oxygen_color, help="Legacy 5' oxygen label color.")
    parser.add_argument("--background", default=DrawingOptions.background, help="Background color.")
    parser.add_argument("--transparent", action="store_true", help="Do not draw a background rectangle.")

    parser.add_argument("--row-spacing", type=float, default=DrawingOptions.row_spacing, help="Vertical distance between base pairs.")
    parser.add_argument("--purine-width", type=float, default=DrawingOptions.purine_width, help="Rectangle width for A/G.")
    parser.add_argument("--pyrimidine-width", type=float, default=DrawingOptions.pyrimidine_width, help="Rectangle width for C/T/U.")
    parser.add_argument("--ambiguous-width", type=float, default=DrawingOptions.ambiguous_width, help="Rectangle width for ambiguous bases.")
    parser.add_argument("--base-height", type=float, default=DrawingOptions.base_height, help="Base rectangle height.")
    parser.add_argument("--sugar-radius", type=float, default=DrawingOptions.sugar_radius, help="Pentagon radius for sugar.")
    parser.add_argument("--phosphate-radius", type=float, default=DrawingOptions.phosphate_radius, help="Circle radius for phosphate; default diameter is close to the sugar-pentagon edge length.")
    parser.add_argument("--pair-gap", type=float, default=None,
                        help="Target gap for a canonical purine-pyrimidine base pair. If omitted, this is 1.5x the sugar-pentagon edge length.")
    parser.add_argument("--base-bond-length", type=float, default=DrawingOptions.base_bond_length,
                        help="Length of each sugar-to-base bond; default is matched to the backbone bond scale.")
    parser.add_argument("--phosphate-kink-length", type=float, default=DrawingOptions.phosphate_kink_length,
                        help="Length of the short sugar-to-kink segment in each two-segment phosphate bond; default is matched to the base-bond scale.")
    parser.add_argument("--font-size", type=float, default=DrawingOptions.font_size, help="Base label font size.")
    parser.add_argument("--terminal-font-size", type=float, default=DrawingOptions.terminal_font_size, help="5'/3'/O label font size.")
    parser.add_argument("--stroke-width", type=float, default=DrawingOptions.stroke_width, help="Stroke width for outlines and bonds.")
    parser.add_argument("--base-label-padding", type=float, default=DrawingOptions.base_label_padding,
                        help="Horizontal padding used for left/right aligned base labels inside their rectangles.")
    parser.add_argument("--rotate-bottom-up-labels", action="store_true",
                        help="Rotate base labels 180 degrees on any strand whose 5'->3' direction is bottom-to-up.")

    parser.add_argument("--no-residue-numbers", action="store_true", help="Show only base letters, not residue numbers.")
    parser.add_argument("--no-terminal-labels", action="store_true", help="Hide 5' and 3' labels.")
    parser.add_argument("--terminal-style", choices=["phosphate", "oxygen", "none"], default=DrawingOptions.terminal_style,
                        help="How to draw the two terminal 5' positions.")
    parser.add_argument("--show-base-pair-lines", action="store_true", help="Draw dashed base-pair lines: two for A:T/A:U and three for G:C.")
    parser.add_argument("--base-pair-line-color", default=DrawingOptions.base_pair_line_color, help="Color for optional base-pair lines.")
    parser.add_argument("--title", default=DrawingOptions.title, help="SVG title metadata.")
    return parser


def options_from_args(args: argparse.Namespace) -> DrawingOptions:
    return DrawingOptions(
        left_color=args.left_color,
        right_color=args.right_color,
        outline_color=args.outline_color,
        text_color=args.text_color,
        oxygen_color=args.oxygen_color,
        background=args.background,
        transparent=args.transparent,
        input_chain=args.input_chain,
        left_direction=normalize_left_direction(args.left_direction),
        nucleic_acid=args.nucleic_acid,
        show_numbers=not args.no_residue_numbers,
        terminal_labels=not args.no_terminal_labels,
        terminal_style=args.terminal_style,
        show_base_pair_lines=args.show_base_pair_lines,
        base_pair_line_color=args.base_pair_line_color,
        row_spacing=args.row_spacing,
        purine_width=args.purine_width,
        pyrimidine_width=args.pyrimidine_width,
        ambiguous_width=args.ambiguous_width,
        base_height=args.base_height,
        sugar_radius=args.sugar_radius,
        phosphate_radius=args.phosphate_radius,
        pair_gap=args.pair_gap if args.pair_gap is not None else 1.5 * 2.0 * args.sugar_radius * math.sin(math.pi / 5.0),
        base_bond_length=args.base_bond_length,
        phosphate_kink_length=args.phosphate_kink_length,
        font_size=args.font_size,
        base_label_padding=args.base_label_padding,
        rotate_bottom_up_labels=args.rotate_bottom_up_labels,
        terminal_font_size=args.terminal_font_size,
        stroke_width=args.stroke_width,
        title=args.title,
    )


def run_gui() -> int:
    try:
        import tkinter as tk
        from tkinter import colorchooser, filedialog, messagebox, ttk
    except Exception as exc:  # pragma: no cover - depends on local Python build
        sys.stderr.write("Could not import Tkinter for GUI mode: {}\n".format(exc))
        return 2

    try:
        root = tk.Tk()
    except Exception as exc:  # pragma: no cover - depends on display availability
        sys.stderr.write("Could not open Tkinter GUI window: {}\n".format(exc))
        return 2

    root.title("{} SVG generator - {}".format(APP_NAME, SCRIPT_VERSION))
    root.geometry("980x820")

    icon_path = resource_path(os.path.join("assets", "helix_2d_flat_icon.png"))
    try:
        app_icon = tk.PhotoImage(file=icon_path)
        root.iconphoto(True, app_icon)
        root._helix_app_icon = app_icon
    except Exception:
        # The icon is optional when the script is copied without its assets.
        pass

    defaults = DrawingOptions()
    main = ttk.Frame(root, padding=12)
    main.pack(fill="both", expand=True)
    main.columnconfigure(1, weight=1)

    seq_var = tk.StringVar(value="CGCGCGCGCGCG")
    output_var = tk.StringVar(value=os.path.abspath("helix_2D_flatV3_out.svg"))
    input_chain_var = tk.StringVar(value=defaults.input_chain)
    left_direction_var = tk.StringVar(value=defaults.left_direction)
    nucleic_acid_var = tk.StringVar(value=defaults.nucleic_acid)
    left_color_var = tk.StringVar(value=defaults.left_color)
    right_color_var = tk.StringVar(value=defaults.right_color)
    outline_color_var = tk.StringVar(value=defaults.outline_color)
    text_color_var = tk.StringVar(value=defaults.text_color)
    oxygen_color_var = tk.StringVar(value=defaults.oxygen_color)
    background_var = tk.StringVar(value=defaults.background)
    terminal_style_var = tk.StringVar(value=defaults.terminal_style)

    row_spacing_var = tk.StringVar(value=str(defaults.row_spacing))
    purine_width_var = tk.StringVar(value=str(defaults.purine_width))
    pyrimidine_width_var = tk.StringVar(value=str(defaults.pyrimidine_width))
    base_height_var = tk.StringVar(value=str(defaults.base_height))
    sugar_radius_var = tk.StringVar(value=str(defaults.sugar_radius))
    phosphate_radius_var = tk.StringVar(value=str(defaults.phosphate_radius))
    pair_gap_var = tk.StringVar(value="{:.1f}".format(defaults.pair_gap))
    base_bond_length_var = tk.StringVar(value=str(defaults.base_bond_length))
    phosphate_kink_length_var = tk.StringVar(value=str(defaults.phosphate_kink_length))
    font_size_var = tk.StringVar(value=str(defaults.font_size))
    stroke_width_var = tk.StringVar(value=str(defaults.stroke_width))
    base_label_padding_var = tk.StringVar(value=str(defaults.base_label_padding))

    show_numbers_var = tk.BooleanVar(value=defaults.show_numbers)
    terminal_labels_var = tk.BooleanVar(value=defaults.terminal_labels)
    transparent_var = tk.BooleanVar(value=defaults.transparent)
    show_base_pair_lines_var = tk.BooleanVar(value=defaults.show_base_pair_lines)
    rotate_bottom_up_labels_var = tk.BooleanVar(value=defaults.rotate_bottom_up_labels)

    def choose_color(var: tk.StringVar) -> None:
        chosen = colorchooser.askcolor(color=var.get(), parent=root)
        if chosen and chosen[1]:
            var.set(chosen[1])

    def browse_output() -> None:
        path = filedialog.asksaveasfilename(
            parent=root,
            title="Save SVG as",
            defaultextension=".svg",
            filetypes=[("SVG files", "*.svg"), ("All files", "*.*")],
            initialfile=os.path.basename(output_var.get()),
        )
        if path:
            output_var.set(path)

    row_idx = 0
    ttk.Label(
        main,
        text="Script version: {}".format(SCRIPT_VERSION),
        foreground="#555555",
    ).grid(row=row_idx, column=0, columnspan=3, sticky="w", pady=(0, 6))
    row_idx += 1

    ttk.Label(main, text="Sequence, 5' to 3'").grid(row=row_idx, column=0, sticky="nw", padx=(0, 8), pady=3)
    seq_entry = ttk.Entry(main, textvariable=seq_var)
    seq_entry.grid(row=row_idx, column=1, columnspan=2, sticky="ew", pady=3)
    row_idx += 1

    ttk.Label(main, text="Output SVG").grid(row=row_idx, column=0, sticky="w", padx=(0, 8), pady=3)
    ttk.Entry(main, textvariable=output_var).grid(row=row_idx, column=1, sticky="ew", pady=3)
    ttk.Button(main, text="Browse...", command=browse_output).grid(row=row_idx, column=2, sticky="ew", padx=(8, 0), pady=3)
    row_idx += 1

    def add_readonly_combo(label: str, var: tk.StringVar, values: Sequence[str], note: str = "") -> None:
        nonlocal row_idx
        ttk.Label(main, text=label).grid(row=row_idx, column=0, sticky="w", padx=(0, 8), pady=3)
        combo = ttk.Combobox(main, textvariable=var, values=list(values), state="readonly", width=18)
        combo.grid(row=row_idx, column=1, sticky="w", pady=3)
        if note:
            ttk.Label(main, text=note).grid(row=row_idx, column=2, sticky="w", padx=(8, 0), pady=3)
        row_idx += 1

    add_readonly_combo("Input chain", input_chain_var, ("left", "right"), "which side gets the input sequence")

    # Some Tk themes render ttk.OptionMenu/dropdown arrows poorly. Radio buttons
    # are clearer here and make the strand-direction setting visible at a glance.
    ttk.Label(main, text="Left 5'->3' direction").grid(row=row_idx, column=0, sticky="w", padx=(0, 8), pady=3)
    direction_frame = ttk.Frame(main)
    direction_frame.grid(row=row_idx, column=1, columnspan=2, sticky="w", pady=3)
    ttk.Radiobutton(
        direction_frame,
        text="bottom-to-up (5' bottom, 3' top)",
        variable=left_direction_var,
        value="bottom-to-up",
    ).grid(row=0, column=0, sticky="w", padx=(0, 18))
    ttk.Radiobutton(
        direction_frame,
        text="up-to-bottom (5' top, 3' bottom)",
        variable=left_direction_var,
        value="up-to-bottom",
    ).grid(row=0, column=1, sticky="w")
    ttk.Label(
        direction_frame,
        text="Right strand is always antiparallel.",
        foreground="#555555",
    ).grid(row=1, column=0, columnspan=2, sticky="w", pady=(2, 0))
    row_idx += 1

    add_readonly_combo("Nucleic acid", nucleic_acid_var, ("auto", "dna", "rna"))

    color_rows = [
        ("Left chain color", left_color_var),
        ("Right chain color", right_color_var),
        ("Outline/bond color", outline_color_var),
        ("Text color", text_color_var),
        ("Legacy 5' oxygen color", oxygen_color_var),
        ("Background", background_var),
    ]
    for label, var in color_rows:
        ttk.Label(main, text=label).grid(row=row_idx, column=0, sticky="w", padx=(0, 8), pady=3)
        ttk.Entry(main, textvariable=var, width=14).grid(row=row_idx, column=1, sticky="w", pady=3)
        ttk.Button(main, text="Color...", command=lambda v=var: choose_color(v)).grid(row=row_idx, column=2, sticky="w", padx=(8, 0), pady=3)
        row_idx += 1

    ttk.Separator(main, orient="horizontal").grid(row=row_idx, column=0, columnspan=3, sticky="ew", pady=8)
    row_idx += 1

    numeric_frame = ttk.Frame(main)
    numeric_frame.grid(row=row_idx, column=0, columnspan=3, sticky="ew")
    for c in range(4):
        numeric_frame.columnconfigure(c, weight=1)

    def add_numeric(parent: ttk.Frame, r: int, c: int, label: str, var: tk.StringVar) -> None:
        ttk.Label(parent, text=label).grid(row=r, column=c * 2, sticky="w", padx=(0, 6), pady=3)
        ttk.Entry(parent, textvariable=var, width=9).grid(row=r, column=c * 2 + 1, sticky="w", padx=(0, 16), pady=3)

    numeric_items = [
        ("Row spacing", row_spacing_var),
        ("Purine width", purine_width_var),
        ("Pyrimidine width", pyrimidine_width_var),
        ("Base height", base_height_var),
        ("Sugar radius", sugar_radius_var),
        ("Phosphate radius", phosphate_radius_var),
        ("Pair gap", pair_gap_var),
        ("Base bond length", base_bond_length_var),
        ("Backbone bond length", phosphate_kink_length_var),
        ("Font size", font_size_var),
        ("Label padding", base_label_padding_var),
        ("Stroke width", stroke_width_var),
    ]
    for idx, (label, var) in enumerate(numeric_items):
        add_numeric(numeric_frame, idx // 2, idx % 2, label, var)
    row_idx += 1

    ttk.Separator(main, orient="horizontal").grid(row=row_idx, column=0, columnspan=3, sticky="ew", pady=8)
    row_idx += 1

    add_readonly_combo("Terminal 5' style", terminal_style_var, ("phosphate", "oxygen", "none"))

    ttk.Checkbutton(main, text="Show residue numbers", variable=show_numbers_var).grid(row=row_idx, column=0, columnspan=2, sticky="w", pady=3)
    row_idx += 1
    ttk.Checkbutton(main, text="Show 5' and 3' labels", variable=terminal_labels_var).grid(row=row_idx, column=0, columnspan=2, sticky="w", pady=3)
    row_idx += 1
    ttk.Checkbutton(main, text="Transparent background", variable=transparent_var).grid(row=row_idx, column=0, columnspan=2, sticky="w", pady=3)
    row_idx += 1
    ttk.Checkbutton(main, text="Show dashed base-pair lines (A:T=2, G:C=3)", variable=show_base_pair_lines_var).grid(row=row_idx, column=0, columnspan=2, sticky="w", pady=3)
    row_idx += 1
    ttk.Checkbutton(main, text="Rotate labels on bottom-to-up strands", variable=rotate_bottom_up_labels_var).grid(row=row_idx, column=0, columnspan=2, sticky="w", pady=3)
    row_idx += 1

    def parse_float(var: tk.StringVar, name: str) -> float:
        try:
            return float(var.get())
        except ValueError:
            raise ValueError("{} must be a number.".format(name))

    def generate_svg_from_gui() -> None:
        try:
            seq = clean_sequence_text(seq_var.get())
            options = DrawingOptions(
                left_color=left_color_var.get(),
                right_color=right_color_var.get(),
                outline_color=outline_color_var.get(),
                text_color=text_color_var.get(),
                oxygen_color=oxygen_color_var.get(),
                background=background_var.get(),
                transparent=transparent_var.get(),
                input_chain=input_chain_var.get(),
                left_direction=normalize_left_direction(left_direction_var.get()),
                nucleic_acid=nucleic_acid_var.get(),
                show_numbers=show_numbers_var.get(),
                terminal_labels=terminal_labels_var.get(),
                terminal_style=terminal_style_var.get(),
                show_base_pair_lines=show_base_pair_lines_var.get(),
                row_spacing=parse_float(row_spacing_var, "Row spacing"),
                purine_width=parse_float(purine_width_var, "Purine width"),
                pyrimidine_width=parse_float(pyrimidine_width_var, "Pyrimidine width"),
                base_height=parse_float(base_height_var, "Base height"),
                sugar_radius=parse_float(sugar_radius_var, "Sugar radius"),
                phosphate_radius=parse_float(phosphate_radius_var, "Phosphate radius"),
                pair_gap=parse_float(pair_gap_var, "Pair gap"),
                base_bond_length=parse_float(base_bond_length_var, "Base bond length"),
                phosphate_kink_length=parse_float(phosphate_kink_length_var, "Backbone bond length"),
                font_size=parse_float(font_size_var, "Font size"),
                base_label_padding=parse_float(base_label_padding_var, "Label padding"),
                rotate_bottom_up_labels=rotate_bottom_up_labels_var.get(),
                stroke_width=parse_float(stroke_width_var, "Stroke width"),
            )
            output_path = derive_output_path(output_var.get().strip() or None, None)
            write_svg_file(seq, output_path, options)
            messagebox.showinfo("Done", "SVG written to:\n{}".format(os.path.abspath(output_path)), parent=root)
        except Exception as exc:
            messagebox.showerror("Error", str(exc), parent=root)

    ttk.Button(main, text="Generate SVG", command=generate_svg_from_gui).grid(row=row_idx, column=0, columnspan=3, sticky="ew", pady=(12, 0))

    seq_entry.focus_set()
    root.mainloop()
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    if len(argv) == 0 or "--gui" in argv:
        return run_gui()

    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        seq, input_file = load_sequence_from_args(args)
        output_path = derive_output_path(args.output, input_file)
        options = options_from_args(args)
        write_svg_file(seq, output_path, options)
        print("Wrote {}".format(os.path.abspath(output_path)))
        return 0
    except Exception as exc:
        sys.stderr.write("Error: {}\n".format(exc))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
