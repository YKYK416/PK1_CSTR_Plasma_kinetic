#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create a prospective blocked factorial design for Fe(110) pathway tests.

This is a pre-registered design artifact, not an experimental result.  It
separates independent reactor/catalyst runs from within-run time traces and
keeps the calibrated local-field condition explicit.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import random
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from _paths import LEGACY_HONG_ROOT

PROJECT = LEGACY_HONG_ROOT
DEFAULT_OUTPUT = PROJECT / "analysis" / "p16_fe110_discriminating_design_20260901_r2"
SEED = 20260901

SURFACES = ("Fe(110)-representative catalyst", "Matched inert surface")
FIELDS = ("F_s,low (independently calibrated)", "F_s,high (independently calibrated)")
ISOTOPES = ("14N2/H2", "15N2/H2")
BLOCKS = ("Block 1", "Block 2", "Block 3")


def build_schedule() -> list[dict]:
    arms = list(itertools.product(SURFACES, FIELDS, ISOTOPES))
    records = []
    for block_index, block in enumerate(BLOCKS, start=1):
        ordered = arms.copy()
        random.Random(SEED + block_index).shuffle(ordered)
        for run_order, (surface, field, isotope) in enumerate(ordered, start=1):
            records.append({
                "block": block,
                "block_run_order": run_order,
                "surface": surface,
                "local_field_target": field,
                "nitrogen_feed": isotope,
                "independent_unit": f"{block.replace(' ', '')}_run{run_order:02d}",
                "mandatory_acceptance": "matched p,T,xN2,tau; field calibration; V-I-Q closure; blank/isotope checks",
            })
    return records


def box(ax, x, y, title, body, color, width=0.23, height=0.16):
    ax.add_patch(FancyBboxPatch((x - width / 2, y - height / 2), width, height, boxstyle="round,pad=0.012",
                                facecolor=color, edgecolor="#666666", linewidth=0.55))
    ax.text(x, y + 0.031, title, ha="center", va="center", fontsize=6.8, fontweight="bold")
    ax.text(x, y - 0.031, body, ha="center", va="center", fontsize=5.95)


def arrow(ax, start, end, label=""):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=8.5, linewidth=0.85, color="#626262"))
    if label:
        ax.text((start[0] + end[0]) / 2, (start[1] + end[1]) / 2 + 0.03, label, ha="center", va="center", fontsize=5.8)


def draw_figure(output: Path, records: list[dict]) -> None:
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.2, "axes.titlesize": 9.1,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(12.2, 8.4))
    grid = fig.add_gridspec(2, 2, left=0.075, right=0.98, bottom=0.085, top=0.93, hspace=0.47, wspace=0.28)

    ax = fig.add_subplot(grid[0, 0]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title("Prospective 2$^3$ mechanism-discrimination factorial", loc="left", pad=8)
    box(ax, 0.16, 0.70, "Surface", "Fe(110)-representative\nvs matched inert", "#E2F0D9")
    box(ax, 0.16, 0.38, "Local field", "F$_s$,low vs F$_s$,high\nfrom independent calibration", "#FFF2CC")
    box(ax, 0.16, 0.10, "N source", "$^{14}$N$_2$/H$_2$ vs\n$^{15}$N$_2$/H$_2$", "#DDEBF7")
    ax.text(0.45, 0.54, "2 × 2 × 2\n= 8 treatment arms", ha="center", va="center", fontsize=10, fontweight="bold")
    arrow(ax, (0.28, 0.70), (0.37, 0.57)); arrow(ax, (0.28, 0.38), (0.37, 0.52)); arrow(ax, (0.28, 0.10), (0.37, 0.47))
    box(ax, 0.76, 0.54, "Independent block", "reconditioned catalyst /\nreactor-day unit", "#FCE4D6", width=0.27)
    arrow(ax, (0.54, 0.54), (0.62, 0.54))
    ax.text(0.58, 0.65, "randomize 8 arms", ha="center", fontsize=5.9)
    ax.text(0.51, 0.19, "Keep pressure, gas temperature, x$_{N2}$, residence time,\nflow and measurement protocol matched inside each block.", ha="center", fontsize=6.35, style="italic")
    ax.text(-0.10, 1.06, "A", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    ax = fig.add_subplot(grid[0, 1]); ax.axis("off")
    ax.set_title("Independent replication and run-order guardrail", loc="left", pad=8)
    rows = [[record["block"], record["block_run_order"],
             "Fe" if record["surface"].startswith("Fe") else "inert",
             "low" if "low" in record["local_field_target"] else "high",
             record["nitrogen_feed"].replace("/H2", "")]
            for record in records]
    table = ax.table(cellText=rows, colLabels=["Block", "Order", "Surface", "F$_s$", "N feed"],
                     cellLoc="center", colLoc="center", colWidths=[0.24, 0.14, 0.22, 0.18, 0.22], bbox=[0.03, 0.12, 0.94, 0.74])
    table.auto_set_font_size(False); table.set_fontsize(5.7)
    for (row, _col), cell in table.get_celld().items():
        cell.set_linewidth(0.30)
        if row == 0: cell.set_facecolor("#DDEBF7"); cell.set_text_props(weight="bold")
        elif row % 2 == 0: cell.set_facecolor("#F7F7F7")
    ax.text(0.03, 0.035, "Seed = 20260901. Three blocks shown are a reproducible pilot / variance-estimation layout,\nnot a formal power calculation. Do not count time samples from one run as independent replicates.", fontsize=6.15, style="italic")
    ax.text(-0.13, 1.06, "B", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    ax = fig.add_subplot(grid[1, 0]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title("Pre-registered measurement panel and acceptance gates", loc="left", pad=8)
    measurements = [
        ("Electrical", "V–I–Q and applied\nwaveform", "#DDEBF7"),
        ("Gas state", "NH / NH$_x$ diagnostic;\nN, N(2D), N(2P) proxy", "#FCE4D6"),
        ("N source", "$^{15}$NH$_3$ plus blank\nand memory checks", "#E2F0D9"),
        ("Surface", "operando / before-after\ncoverage & material check", "#FFF2CC"),
    ]
    for (title, body, color), x in zip(measurements, (0.14, 0.38, 0.62, 0.86)):
        box(ax, x, 0.63, title, body, color, width=0.20, height=0.18)
        arrow(ax, (x, 0.53), (x, 0.36))
    box(ax, 0.50, 0.23, "Acceptance before route interpretation", "field-transfer calibration; matched reactor state; isotope closure;\nblank pass; independent block replication", "#F4CCCC", width=0.56, height=0.16)
    ax.text(0.50, 0.055, "NH$_3$ alone is intentionally excluded as a route-identification observable.", ha="center", fontsize=6.5, style="italic")
    ax.text(-0.10, 1.06, "C", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    ax = fig.add_subplot(grid[1, 1]); ax.axis("off")
    ax.set_title("What each contrast can—and cannot—falsify", loc="left", pad=8)
    rows = [
        ("Fe vs inert", "Tests material boundary after\nmatching electrical input.", "Does not identify\nER or LH."),
        ("F$_s$,low vs F$_s$,high", "Tests field response after\nindependent F$_s$ calibration.", "Cannot use E/N\nas F$_s$ surrogate."),
        ("$^{15}$N$_2$ vs $^{14}$N$_2$", "Tests feed-N provenance\nwith blank acceptance.", "Does not select a\nfirst-N–H route."),
        ("State + surface panel", "Constrains ER–LH with\nstate + surface signatures.", "Requires hybrid ROP\nfor flux ranking."),
    ]
    table = ax.table(cellText=rows, colLabels=["Contrast", "Can test", "Cannot establish alone"], cellLoc="left", colLoc="left",
                     colWidths=[0.20, 0.42, 0.38], bbox=[0.0, 0.18, 1.0, 0.65])
    table.auto_set_font_size(False); table.set_fontsize(6.0)
    for (row, _col), cell in table.get_celld().items():
        cell.set_linewidth(0.35)
        if row == 0: cell.set_facecolor("#DDEBF7"); cell.set_text_props(weight="bold")
        elif row % 2 == 0: cell.set_facecolor("#F7F7F7")
    ax.text(0.0, 0.055, "Analysis unit: independent block/run. Model blocks and repeated time samples explicitly;\npredefine the route-ranking analysis only after Figure 22 G0–G3 are passed.", fontsize=6.25, style="italic")
    ax.text(-0.13, 1.06, "D", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    fig.text(0.5, 0.017, "Figure 25. Prospective blocked factorial validation design for distinguishing gas–surface pathway hypotheses.\nIt is a pre-registered protocol, not an experimental data figure or a claimed surface-flux result.", ha="center", fontsize=7.45)
    fig.savefig(output / "Figure_25_fe110_pathway_discrimination_design.png", dpi=600, bbox_inches="tight")
    fig.savefig(output / "Figure_25_fe110_pathway_discrimination_design.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(output: Path, records: list[dict]) -> None:
    lines = [
        "# Prospective Fe(110) pathway-discrimination validation design",
        "",
        "**Status:** pre-registered proposed experiment. No experimental data, surface flux or pathway ranking is reported.",
        "",
        "## Question",
        "",
        "Can a material-dependent and field-calibrated response, together with nitrogen-source closure and state/surface diagnostics, distinguish surface first-N–H alternatives from a gas-phase-only explanation?",
        "",
        "## Design",
        "",
        "A blocked 2^3 factorial crosses: (i) Fe(110)-representative versus matched inert surface, (ii) independently calibrated low versus high local surface field, and (iii) 14N2/H2 versus 15N2/H2. The eight arms are randomized within each independent reactor/catalyst reconditioning block. The supplied schedule contains three blocks (24 runs) with seed 20260901.",
        "",
        "The independent unit is a completed, independently reconditioned reactor/catalyst run. Time-resolved signals acquired within that run are repeated measurements, not independent replicates. Three blocks are a transparent pilot and variance-estimation layout; an effect-size and variance-informed power calculation is still required before a confirmatory campaign.",
        "",
        "## Mandatory controls",
        "",
        "- Match pressure, temperature, inlet composition, residence time, flow and acquisition protocol within every block.",
        "- Calibrate F_s independently; do not identify it with gas E/N or fit it only to NH3.",
        "- Record V–I–Q and waveform for every arm so material effects are not confused with electrical mismatch.",
        "- Include isotope blanks, carry-over/cleaning checks and nitrogen-balance acceptance before interpreting 15NH3.",
        "- Pre-register the hierarchical analysis with block as a grouping factor and within-run traces as repeated measures.",
        "",
        "## Inference boundary",
        "",
        "The factorial contrasts can test a material-conditioned response, a field-conditioned response and feed-nitrogen provenance. They cannot alone establish the dominant ER or LH first-N–H flux. Such a ranking remains conditional on the field-transfer and shared-boundary/hybrid-ROP gates stated in Figure 22.",
        "",
        "## Deliverables",
        "",
        "- `Figure_25_fe110_pathway_discrimination_design.png/.pdf`",
        "- `P16_fe110_pathway_factorial_schedule.csv`",
        "- `P16_fe110_pathway_discrimination_design.md`",
    ]
    (output / "P16_fe110_pathway_discrimination_design.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Prospective blocked factorial design for Fe(110) pathway discrimination")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    records = build_schedule()
    with (output / "P16_fe110_pathway_factorial_schedule.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=records[0].keys()); writer.writeheader(); writer.writerows(records)
    draw_figure(output, records)
    write_report(output, records)
    print(f"independent_runs={len(records)}")
    print(f"blocks={len(BLOCKS)}; arms_per_block={len(records) // len(BLOCKS)}; seed={SEED}")
    print(f"output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
