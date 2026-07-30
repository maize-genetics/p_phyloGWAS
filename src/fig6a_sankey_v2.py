"""
Figure 6A (Sankey version): OG filtering across three analysis stages
with envPC overlap structure shown at the phylogenetic step.

19,613 OGs are tested. The phylogenetic association test partitions the
significant subset into envPC1 / envPC2 / envPC3 groups, with overlaps
between them. Selection shift and differential expression further filter
each group independently.
"""

import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
import numpy as np

PHYLOGWAS_ROOT = os.environ.get("PHYLOGWAS_ROOT", "/workdir/sh2246/p_phyloGWAS")
OUT_DIR = os.path.join(PHYLOGWAS_ROOT, "output/figure")

# ════════════════════════════════════════════════════════════════════════
# DATA
# ════════════════════════════════════════════════════════════════════════

TOTAL_TESTED = 19_613

# Overlap structure at the phylogenetic stage
# (these are the *exclusive* counts in each Venn region)
PHYLO = {
    "envPC1_only":           156,   # 160 - 3 - 1 - 0
    "envPC2_only":           103,   # 109 - 3 - 3 - 0
    "envPC3_only":            64,   #  68 - 1 - 3 - 0
    "envPC1_envPC2":           3,
    "envPC1_envPC3":           1,
    "envPC2_envPC3":           3,
    "envPC1_envPC2_envPC3":    0,
}

# Group totals at each downstream stage (no overlap data shown beyond phylo)
# envPC2/envPC3 corrected 2026-07-30: the RELAX-significance layer (hyphyCandidate2/3 in
# 11_candidateOGInvestigation.ipynb) referenced a nonexistent column on the wrong dataframe
# (hyphy_drought$wet/hyphy_clay$wet instead of hyphy_wet$OG/hyphy_clay$OG), silently dropping
# every wet/clay-significant OG. envPC1 (cold/warm) was unaffected and stays the same.
SELECTION = {"envPC1": 30, "envPC2": 28, "envPC3": 13}
EXPRESSION = {"envPC1": 10, "envPC2": 10, "envPC3": 7}

# ════════════════════════════════════════════════════════════════════════
# STYLE
# ════════════════════════════════════════════════════════════════════════

COLORS = {
    "envPC1":  "#185FA5",
    "envPC2":  "#3B6D11",
    "envPC3":  "#BA7517",
    "neutral": "#5F5E5A",
    "filtered":"#B4B2A9",
}

RIBBON_ALPHA = 0.45      # pure-color ribbons
OVERLAP_ALPHA = 0.55     # overlap ribbons (gradient-like via blend)
FILTERED_ALPHA = 0.10    # gray "filtered out" ribbons

# Visual minimums so tiny counts (1, 3) don't vanish
MIN_RIBBON_HEIGHT = 4.0
MIN_BAR_SEGMENT  = 4.0

# ════════════════════════════════════════════════════════════════════════
# LAYOUT
# ════════════════════════════════════════════════════════════════════════

FIG_W, FIG_H = 6, 2.8           # inches
COLS = {                       # column x-positions (left edge of bar)
    "tested":     0.5,
    "phylo":      4.5,
    "selection":  8.5,
    "expression":12.5,
}
BAR_WIDTH = 0.35

# y-extent: 0 (bottom) to ~6 (top of plotting area)
# Larger y values = higher on screen
Y_TOP = 6.0
Y_BOTTOM = -17.5   # leave room for legend and filtered ribbon below the bars
GAP_BETWEEN_GROUPS = 0.2      # vertical gap between envPC bars at col 2+

# Height scaling: sqrt so small numbers stay visible
def scale_height(count, factor=0.10):
    """Convert a count to a bar height using square-root scaling."""
    if count <= 0:
        return 0
    h = np.sqrt(count) * factor
    return max(h, MIN_BAR_SEGMENT * 0.01)  # in axis units


# ════════════════════════════════════════════════════════════════════════
# RIBBON DRAWING
# ════════════════════════════════════════════════════════════════════════

def draw_ribbon(ax, x0, y0_top, y0_bot, x1, y1_top, y1_bot,
                color, alpha=0.5, zorder=1):
    """
    Draw a Sankey ribbon as a filled shape with two cubic Bezier curves
    (top edge and bottom edge) connecting two vertical slices.

    The control points are placed at the horizontal midpoint with the
    same y as their respective endpoint — this gives the classic
    smooth S-curve flow.
    """
    # Horizontal midpoint for control points
    cx0 = x0 + 0.5 * (x1 - x0)
    cx1 = x0 + 0.5 * (x1 - x0)

    verts = [
        (x0, y0_top),                # start top
        (cx0, y0_top),               # control 1
        (cx1, y1_top),               # control 2
        (x1, y1_top),                # end top
        (x1, y1_bot),                # down right edge
        (cx1, y1_bot),               # control 1 (bottom curve, reversed)
        (cx0, y0_bot),               # control 2
        (x0, y0_bot),                # end bottom
        (x0, y0_top),                # close
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.CLOSEPOLY,
    ]
    path = Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor=color, edgecolor="none",
                                alpha=alpha, zorder=zorder)
    ax.add_patch(patch)


def blend_colors(c1, c2, w=0.5):
    """Linear blend of two hex colors → rgb tuple."""
    def hex2rgb(h):
        h = h.lstrip("#")
        return tuple(int(h[i:i+2], 16) / 255 for i in (0, 2, 4))
    r1 = hex2rgb(c1); r2 = hex2rgb(c2)
    return tuple(r1[i] * (1 - w) + r2[i] * w for i in range(3))


# ════════════════════════════════════════════════════════════════════════
# BUILD FIGURE
# ════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H), dpi=300)
ax.set_xlim(0, 14)
ax.set_ylim(Y_BOTTOM, Y_TOP + 1.5)   # extra space top for labels, bottom for legend
ax.set_aspect("auto")
ax.axis("off")

# ── Stage labels (between columns) ──────────────────────────────────────
stage_label_y = Y_TOP + 1.8
stage_labels = [
    ((COLS["tested"]    + COLS["phylo"])      / 2 + BAR_WIDTH/2,
     "Phylogenetic", "association"),
    ((COLS["phylo"]     + COLS["selection"])  / 2 + BAR_WIDTH/2,
     "Selection",    "shift"),
    ((COLS["selection"] + COLS["expression"]) / 2 + BAR_WIDTH/2,
     "Differential", "expression"),
]
for x, l1, l2 in stage_labels:
    ax.text(x, stage_label_y,      l1, ha="center", va="bottom",
            fontsize=9, fontweight="medium")
    ax.text(x, stage_label_y - 0.18, l2, ha="center", va="top",
            fontsize=9, fontweight="medium")
    # Underline
#    ax.plot([x - 0.9, x + 0.9], [stage_label_y - 0.42, stage_label_y - 0.42],
#            color="#D1D1CC", linewidth=0.5)

# Column subtitles
sub_y = Y_TOP + 0.25
for label, key in [("19,613 OGs", "tested")]:
    ax.text(COLS[key] + BAR_WIDTH/2, sub_y, label,
            ha="center", va="bottom", fontsize=8)

# ════════════════════════════════════════════════════════════════════════
# COLUMN 1: Tested pool (single bar, partitioned virtually for ribbon sources)
# ════════════════════════════════════════════════════════════════════════

# Partition heights for ribbon source slices (sqrt-scaled)
parts = [
    ("envPC1_only",         PHYLO["envPC1_only"],        COLORS["envPC1"], None),
    ("envPC1_envPC2",       PHYLO["envPC1_envPC2"],      None,             ("envPC1","envPC2")),
    ("envPC2_only",         PHYLO["envPC2_only"],        COLORS["envPC2"], None),
    ("envPC2_envPC3",       PHYLO["envPC2_envPC3"],      None,             ("envPC2","envPC3")),
    ("envPC3_only",         PHYLO["envPC3_only"],        COLORS["envPC3"], None),
    ("envPC1_envPC3",       PHYLO["envPC1_envPC3"],      None,             ("envPC1","envPC3")),
]
part_heights = {p[0]: max(scale_height(p[1], factor=0.18), MIN_BAR_SEGMENT * 0.01) for p in parts}

total_top = Y_TOP - 0.2
significant_total_h = sum(part_heights.values())

# Filtered-out bar (gray, capped visually)
filtered_h = 14
filtered_top = total_top - significant_total_h - 0.1
filtered_bot = filtered_top - filtered_h

# Draw "tested" column as a single neutral bar + filtered gray below
ax.add_patch(mpatches.Rectangle(
    (COLS["tested"], total_top - significant_total_h),
    BAR_WIDTH, significant_total_h,
    facecolor=COLORS["neutral"], edgecolor="none", alpha=0.8))
ax.add_patch(mpatches.Rectangle(
    (COLS["tested"], filtered_bot),
    BAR_WIDTH, filtered_h,
    facecolor=COLORS["filtered"], edgecolor="none", alpha=0.7))

# Labels
ax.text(COLS["tested"] - BAR_WIDTH, total_top - significant_total_h/2,
        f"330", ha="center", va="center",
        fontsize=8, color="#2C2C2A")
#ax.text(COLS["tested"] - 0.08, filtered_bot + filtered_h/2,
#        f"filtered out\n{TOTAL_TESTED - 330:,}", ha="right", va="center",
#        fontsize=9, color="#5F5E5A")

# Compute y ranges for each source partition slice in column 1
src_y = {}
y_cursor = total_top
for key, count, _, _ in parts:
    h = part_heights[key]
    src_y[key] = (y_cursor, y_cursor - h)   # (top, bottom)
    y_cursor -= h

ax.text(COLS["tested"] + BAR_WIDTH+0.08, src_y["envPC1_envPC2"][1]+0.01,
        f"3", ha="center", va="center",
        fontsize=6, color="#2C2C2A")

ax.text(COLS["tested"] + BAR_WIDTH+0.08, src_y["envPC1_envPC3"][1]+0.01,
        f"1", ha="center", va="center",
        fontsize=6, color="#2C2C2A")

ax.text(COLS["tested"] + BAR_WIDTH+0.08, src_y["envPC2_envPC3"][1]+0.01,
        f"3", ha="center", va="center",
        fontsize=6, color="#2C2C2A")
# ════════════════════════════════════════════════════════════════════════
# COLUMN 2: envPC1, envPC2, envPC3 bars
# ════════════════════════════════════════════════════════════════════════

# Heights for each envPC group (proportional to total in group, sqrt-scaled)
group_h = {
    "envPC1": scale_height(160, factor=0.5),
    "envPC2": scale_height(109, factor=0.5),
    "envPC3": scale_height( 68, factor=0.5),
}

# Stack the three group bars vertically at column 2
col2_top = Y_TOP - 0.2
col2_y = {}
y_cursor = col2_top
for grp in ["envPC1", "envPC2", "envPC3"]:
    h = group_h[grp]
    col2_y[grp] = (y_cursor, y_cursor - h)
    y_cursor -= h + GAP_BETWEEN_GROUPS

# Within each envPC bar, allocate sub-slots for unique + overlap segments
# The unique segment goes in the middle; overlap segments at the edges
# facing the adjacent group bar.

def split_bar(grp, top, bottom):
    """
    Subdivide an envPC bar into y-ranges for its component parts.
    Returns dict mapping part_key → (y_top, y_bottom).
    """
    total = bottom_count = {
        "envPC1": 160, "envPC2": 109, "envPC3": 68,
    }[grp]
    H = top - bottom
    slots = {}

    if grp == "envPC1":
        # Order top→bottom: unique, ∩envPC2 (at bottom, faces envPC2 below)
        # Actually let's put: unique (top), ∩envPC3 (middle), ∩envPC2 (bottom)
        # so ∩envPC2 touches the envPC2 bar below
        h_u  = H * PHYLO["envPC1_only"] / total
        h_13 = H * PHYLO["envPC1_envPC3"] / total
        h_12 = H * PHYLO["envPC1_envPC2"] / total
        # Bump tiny ones to a visible minimum
        h_13 = max(h_13, 0.04)
        h_12 = max(h_12, 0.04)
        h_u  = H - h_13 - h_12
        slots["envPC1_only"]      = (top,                top - h_u)
        slots["envPC1_envPC3"]    = (top - h_u,          top - h_u - h_13)
        slots["envPC1_envPC2"]    = (top - h_u - h_13,   bottom)

    elif grp == "envPC2":
        # ∩envPC1 at top (faces envPC1 above), unique middle, ∩envPC3 at bottom
        h_u  = H * PHYLO["envPC2_only"] / total
        h_12 = H * PHYLO["envPC1_envPC2"] / total
        h_23 = H * PHYLO["envPC2_envPC3"] / total
        h_12 = max(h_12, 0.04)
        h_23 = max(h_23, 0.04)
        h_u  = H - h_12 - h_23
        slots["envPC2_only"]    = (top,                  top - h_u)
        slots["envPC1_envPC2"]      = (top - h_u,           top - h_12 - h_u)
        slots["envPC2_envPC3"]    = (top - h_12 - h_u,     bottom)

    else:  # envPC3
        # ∩envPC2 at top (faces envPC2 above), ∩envPC1 small, unique middle/bottom
        h_u  = H * PHYLO["envPC3_only"] / total
        h_23 = H * PHYLO["envPC2_envPC3"] / total
        h_13 = H * PHYLO["envPC1_envPC3"] / total
        h_23 = max(h_23, 0.04)
        h_13 = max(h_13, 0.03)
        h_u  = H - h_23 - h_13
        slots["envPC3_only"]    = (top,                  top - h_u)
        slots["envPC1_envPC3"]    = (top - h_u,           top - h_u - h_13)
        slots["envPC2_envPC3"]      = (top - h_u - h_13,    bottom)

    return slots

col2_slots = {grp: split_bar(grp, *col2_y[grp]) for grp in ["envPC1", "envPC2", "envPC3"]}

# Draw column 2 main bars (colored solidly per group)
for grp, (top, bot) in col2_y.items():
    ax.add_patch(mpatches.Rectangle(
        (COLS["phylo"], bot), BAR_WIDTH, top - bot,
        facecolor=COLORS[grp], edgecolor="none"))
    # Group label to right
    ax.text(COLS["phylo"] - BAR_WIDTH, (top + bot) / 2,
            f"{ {'envPC1':160,'envPC2':109,'envPC3':68}[grp] }",
            ha="center", va="center", fontsize=8, color=COLORS[grp])

# ════════════════════════════════════════════════════════════════════════
# RIBBONS COL 1 → COL 2
# ════════════════════════════════════════════════════════════════════════

x0 = COLS["tested"] + BAR_WIDTH
x1 = COLS["phylo"]

# Unique groups: single ribbon each
unique_routes = [
    ("envPC1_only", "envPC1"),
    ("envPC2_only", "envPC2"),
    ("envPC3_only", "envPC3"),
]
for src_key, tgt_grp in unique_routes:
    y0_top, y0_bot = src_y[src_key]
    y1_top, y1_bot = col2_slots[tgt_grp][src_key]
    draw_ribbon(ax, x0, y0_top, y0_bot, x1, y1_top, y1_bot,
                color=COLORS[tgt_grp], alpha=RIBBON_ALPHA, zorder=2)

# Overlap groups: ribbon splits in two, fanning to both target bars
overlap_routes = [
    ("envPC1_envPC2", "envPC1", "envPC2"),
    ("envPC1_envPC3", "envPC1", "envPC3"),
    ("envPC2_envPC3", "envPC2", "envPC3"),
]
for src_key, tgt_a, tgt_b in overlap_routes:
    y0_top, y0_bot = src_y[src_key]
    blend_color = blend_colors(COLORS[tgt_a], COLORS[tgt_b])
    # Each branch carries the full source thickness (it's the SAME OGs in both)
    y1_top_a, y1_bot_a = col2_slots[tgt_a][src_key]
    y1_top_b, y1_bot_b = col2_slots[tgt_b][src_key]
    draw_ribbon(ax, x0, y0_top, y0_bot, x1, y1_top_a, y1_bot_a,
                color=blend_color, alpha=OVERLAP_ALPHA, zorder=3)
    draw_ribbon(ax, x0, y0_top, y0_bot, x1, y1_top_b, y1_bot_b,
                color=blend_color, alpha=OVERLAP_ALPHA, zorder=3)

# Filtered out: drop downward from below the significant slice
y0_top = filtered_top
y0_bot = filtered_bot
y1_top = y0_bot-3
y1_bot = y0_bot-3
draw_ribbon(ax, x0, y0_top, y0_bot,
            (COLS["tested"] + COLS["phylo"]) / 2, y1_top, y1_bot,
            color=COLORS["filtered"], alpha=FILTERED_ALPHA, zorder=1)

# ════════════════════════════════════════════════════════════════════════
# COLUMN 3: Selection shift
# ════════════════════════════════════════════════════════════════════════

col3_h = {grp: scale_height(SELECTION[grp], factor=0.70) for grp in SELECTION}
col3_y = {}
y_cursor = col2_top
for grp in ["envPC1", "envPC2", "envPC3"]:
    h = col3_h[grp]
    col3_y[grp] = (y_cursor, y_cursor - h)
    y_cursor -= h + GAP_BETWEEN_GROUPS * 1.5

for grp, (top, bot) in col3_y.items():
    ax.add_patch(mpatches.Rectangle(
        (COLS["selection"], bot), BAR_WIDTH, top - bot,
        facecolor=COLORS[grp], edgecolor="none"))
    ax.text(COLS["selection"] - BAR_WIDTH , (top + bot) / 2,
            f"{SELECTION[grp]}", ha="center", va="center",
            fontsize=8, color=COLORS[grp])

# Ribbons col 2 → col 3
x0 = COLS["phylo"] + BAR_WIDTH
x1 = COLS["selection"]
for grp in ["envPC1", "envPC2", "envPC3"]:
    src_top, src_bot = col2_y[grp]
    tgt_top, tgt_bot = col3_y[grp]
    src_h = src_top - src_bot
    tgt_h = (tgt_top - tgt_bot)/49*25
    # Passing ribbon (sized to target)
    draw_ribbon(ax, x0, src_top, src_top - tgt_h,
                x1, tgt_top, tgt_bot,
                color=COLORS[grp], alpha=RIBBON_ALPHA, zorder=2)
    # Filtered drop
    drop_top = src_top - tgt_h
    drop_bot = src_bot
    if drop_top > drop_bot:
        mid_x = (x0 + x1) / 2
        mid_y_top = drop_bot - 5
        mid_y_bot = drop_bot - 5
        draw_ribbon(ax, x0, drop_top, drop_bot,
                    mid_x, mid_y_top, mid_y_bot,
                    color=COLORS["filtered"], alpha=FILTERED_ALPHA, zorder=1)

# ════════════════════════════════════════════════════════════════════════
# COLUMN 4: Differential expression
# ════════════════════════════════════════════════════════════════════════

col4_h = {grp: scale_height(EXPRESSION[grp], factor=0.9) for grp in EXPRESSION}
col4_y = {}
y_cursor = col2_top
for grp in ["envPC1", "envPC2", "envPC3"]:
    h = col4_h[grp]
    col4_y[grp] = (y_cursor, y_cursor - h)
    y_cursor -= h + GAP_BETWEEN_GROUPS * 1.5

for grp, (top, bot) in col4_y.items():
    ax.add_patch(mpatches.Rectangle(
        (COLS["expression"], bot), BAR_WIDTH, top - bot,
        facecolor=COLORS[grp], edgecolor="none"))
    ax.text(COLS["expression"] - BAR_WIDTH, (top + bot) / 2,
            f"{EXPRESSION[grp]}", ha="center", va="center",
            fontsize=8, color=COLORS[grp], fontweight="medium")

# Ribbons col 3 → col 4
x0 = COLS["selection"] + BAR_WIDTH
x1 = COLS["expression"]
for grp in ["envPC1", "envPC2", "envPC3"]:
    src_top, src_bot = col3_y[grp]
    tgt_top, tgt_bot = col4_y[grp]
    tgt_h = (tgt_top - tgt_bot)/81*49
    draw_ribbon(ax, x0, src_top, src_top - tgt_h,
                x1, tgt_top, tgt_bot,
                color=COLORS[grp], alpha=RIBBON_ALPHA, zorder=2)
    drop_top = src_top - tgt_h
    drop_bot = src_bot
    if drop_top > drop_bot:
        mid_x = (x0 + x1) / 2
        mid_y_top = drop_bot - 5
        mid_y_bot = drop_bot - 5
        draw_ribbon(ax, x0, drop_top, drop_bot,
                    mid_x, mid_y_top, mid_y_bot,
                    color=COLORS["filtered"], alpha=FILTERED_ALPHA, zorder=1)

# Final group labels far right
for grp, (top, bot) in col4_y.items():
    ax.text(COLS["expression"] + BAR_WIDTH + 0.1, (top + bot) / 2,
            grp, ha="left", va="center",
            fontsize=10, color=COLORS[grp], fontweight="medium")

# ════════════════════════════════════════════════════════════════════════
# OVERLAP LEGEND
# ════════════════════════════════════════════════════════════════════════


# ════════════════════════════════════════════════════════════════════════
# SAVE
# ════════════════════════════════════════════════════════════════════════

plt.tight_layout(pad=0.5)
os.makedirs(OUT_DIR, exist_ok=True)
png_path = os.path.join(OUT_DIR, "fig6a_sankey_v2.png")
svg_path = os.path.join(OUT_DIR, "fig6a_sankey_v2.svg")
plt.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
plt.savefig(svg_path, bbox_inches="tight", facecolor="white")
plt.show()
print(f"Saved {png_path} and {svg_path}")
