#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create a graph-informed, pre-experimental model-discrimination atlas.

This is deliberately a *design* figure, not an experimental-result figure.
It reads only accepted continuous-CSTR ROP/sensitivity outputs and makes the
condition-selection logic auditable.  It never reads rejected pulsed outputs.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parent.parent
FIGURES = PROJECT / "analysis" / "p10_model_discrimination_atlas_20260831"
SOURCE = PROJECT / "analysis" / "figures_20260829" / "Figure_3_dense_EN_N2_source_data.csv"
SENSITIVITY = PROJECT / "analysis" / "p6_local_reaction_control_20260830" / "P6_local_sensitivities.csv"

CONDITIONS = [
    (20, 0.1, "L", "Low-field\nentry-limited"),
    (140, 0.1, "P", "Productive\nreference"),
    (240, 0.9, "H", "High-field, N$_2$-rich\nloss-dominated"),
]
CHANNELS = [
    ("named_h2star_pct", "Named electronic H$_2$*", "#0072B2"),
    ("association_pct", "H + N + M", "#D55E00"),
    ("n2p_h2_pct", "N($^2$P) + H$_2$", "#E69F00"),
    ("rydberg_h2star_pct", "Rydberg-H$_2$*", "#CC79A7"),
    ("n2d_h2_pct", "N($^2$D) + H$_2$", "#56B4E9"),
]


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7.2,
        "axes.labelsize": 8.2, "axes.titlesize": 8.8,
        "xtick.labelsize": 7.0, "ytick.labelsize": 7.0,
        "legend.fontsize": 6.2, "axes.linewidth": 0.7,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def panel(ax: plt.Axes, label: str) -> None:
    ax.text(-0.12, 1.08, label, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top")


def selected_frame() -> pd.DataFrame:
    dense = pd.read_csv(SOURCE)
    rows = []
    for en_td, n2_fraction, code, role in CONDITIONS:
        matched = dense[(dense["en_td"] == en_td) & np.isclose(dense["n2_fraction"], n2_fraction)]
        if len(matched) != 1:
            raise ValueError(f"Expected exactly one source-data row for {en_td} Td, xN2={n2_fraction}")
        row = matched.iloc[0].to_dict()
        row.update({"code": code, "role": role, "condition": f"{en_td} Td, x$_{{N_2}}$={n2_fraction:.1f}"})
        secondary_columns = [key for key, _label, _colour in CHANNELS[1:]]
        row["largest_secondary_pct"] = max(float(row[key]) for key in secondary_columns)
        row["dominance_margin"] = float(row["named_h2star_pct"]) / max(row["largest_secondary_pct"], 1.0e-20)
        rows.append(row)
    return pd.DataFrame(rows)


def draw_entry_hierarchy(ax: plt.Axes, frame: pd.DataFrame) -> None:
    x = np.arange(len(frame))
    bottom = np.zeros(len(frame))
    for key, label, colour in CHANNELS:
        values = frame[key].to_numpy(float)
        ax.bar(x, values, bottom=bottom, width=0.64, color=colour, edgecolor="white", linewidth=0.35, label=label)
        bottom += values
    ax.set_ylim(0, 100)
    ax.set_xticks(x, [f"{row.code}\n{int(row.en_td)} Td, {row.n2_fraction:.1f}" for row in frame.itertuples()])
    ax.set_ylabel("Direct NH-entry flux (%)")
    ax.set_title("Accepted ROP hierarchy selects contrasting test states")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=2)
    for xpos, value in zip(x, frame["named_h2star_pct"]):
        ax.text(xpos, value / 2, f"{value:.1f}%", ha="center", va="center", fontsize=6.5, color="white", fontweight="bold")
    panel(ax, "A")


def draw_logic(ax: plt.Axes, frame: pd.DataFrame) -> None:
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.set_title("Graph-informed discrimination logic", loc="left", pad=8)
    positions = {"L": (0.16, 0.78), "P": (0.16, 0.48), "H": (0.16, 0.18)}
    colors = {"L": "#B3CDE3", "P": "#0072B2", "H": "#E69F00"}
    descriptions = {
        "L": "weak activation\nentry-limited control",
        "P": "H$_2$* direct entry\nstrong NH$_3$ formation",
        "H": "secondary-entry increase\nand turnover/loss stress",
    }
    for row in frame.itertuples():
        x, y = positions[row.code]
        ax.add_patch(FancyBboxPatch((x - 0.12, y - 0.095), 0.24, 0.19, boxstyle="round,pad=0.012",
                                    facecolor=colors[row.code], edgecolor="#3F3F3F", linewidth=0.65))
        text_colour = "white" if row.code == "P" else "black"
        ax.text(x, y + 0.038, row.code, ha="center", va="center", fontsize=10, fontweight="bold", color=text_colour)
        ax.text(x, y - 0.032, descriptions[row.code], ha="center", va="center", fontsize=5.8, color=text_colour)
    ax.add_patch(FancyBboxPatch((0.59, 0.61), 0.32, 0.22, boxstyle="round,pad=0.015", facecolor="#EAF3F8", edgecolor="#3F3F3F", linewidth=0.65))
    ax.text(0.75, 0.76, "Direct-entry layer", ha="center", fontsize=7.0, fontweight="bold")
    ax.text(0.75, 0.67, "H$_2$* → NH versus\nsecondary NH-entry routes", ha="center", fontsize=6.2)
    ax.add_patch(FancyBboxPatch((0.59, 0.20), 0.32, 0.22, boxstyle="round,pad=0.015", facecolor="#E8F3E5", edgecolor="#3F3F3F", linewidth=0.65))
    ax.text(0.75, 0.35, "Turnover layer", ha="center", fontsize=7.0, fontweight="bold")
    ax.text(0.75, 0.26, "NHx hydrogenation, ionic loss,\ntransport and wall/surface loss", ha="center", fontsize=6.2)
    for source, target, rad in (("L", (0.59, 0.72), 0.0), ("P", (0.59, 0.72), 0.12), ("P", (0.59, 0.31), -0.12), ("H", (0.59, 0.31), 0.0)):
        start = positions[source]
        ax.add_patch(FancyArrowPatch((start[0] + 0.125, start[1]), target, arrowstyle="-|>",
                                     mutation_scale=10, lw=1.15, color="#4C4C4C",
                                     connectionstyle=f"arc3,rad={rad}"))
    ax.text(0.52, 0.055, "Prediction is falsified if independent observables require\na different ordering of entry or turnover controls.",
            ha="center", va="center", fontsize=6.1, style="italic")
    panel(ax, "B")


def draw_measurement_matrix(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.set_title("Pre-specified observables and interpretation boundary", loc="left", pad=8)
    rows = [
        ["NH$_3$ outlet", "FTIR / calibrated GC", "Mandatory product and N-balance anchor"],
        ["V–I–Q waveform", "time-resolved electrical", "Constrains power and E/N(t); not n$_e$ by itself"],
        ["NH or NHx", "state-sensitive optical diagnostic", "Separates entry proxy from final NH$_3$"],
        ["$^{15}$N label", "MS / isotope-resolved product", "Closes the nitrogen-source attribution"],
        ["H$_2$* / n$_e$", "only if technically validated", "Constrain model; do not infer from a single OES line"],
    ]
    table = ax.table(cellText=rows, colLabels=["Observable", "Method class", "Use / limitation"],
                     cellLoc="left", colLoc="left", colWidths=[0.19, 0.27, 0.54], bbox=[0.0, 0.03, 1.0, 0.82])
    table.auto_set_font_size(False); table.set_fontsize(6.05)
    for (row, _column), cell in table.get_celld().items():
        cell.set_linewidth(0.35)
        if row == 0:
            cell.set_facecolor("#DDEBF7"); cell.set_text_props(weight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#F7F7F7")
    panel(ax, "C")


def draw_decision_ledger(ax: plt.Axes, frame: pd.DataFrame) -> None:
    ax.axis("off")
    ax.set_title("Selection ledger from accepted continuous-CSTR evidence", loc="left", pad=8)
    rows = []
    for row in frame.itertuples():
        rows.append([
            row.code,
            f"{int(row.en_td)} Td, {row.n2_fraction:.1f}",
            f"{row.named_h2star_pct:.2f}",
            f"{row.largest_secondary_pct:.2f}",
            f"{row.dominance_margin:.2f}",
            {"L": "entry contrast", "P": "turnover reference", "H": "secondary/loss contrast"}[row.code],
        ])
    table = ax.table(cellText=rows,
                     colLabels=["State", "E/N, x$_{N_2}$", "Named H$_2$*\n(%)", "Largest secondary\n(%)", "Margin\n(×)", "Use"],
                     cellLoc="center", colLoc="center", colWidths=[0.08, 0.19, 0.16, 0.18, 0.10, 0.29],
                     bbox=[0.0, 0.36, 1.0, 0.48])
    table.auto_set_font_size(False); table.set_fontsize(6.15)
    for (row, _column), cell in table.get_celld().items():
        cell.set_linewidth(0.35)
        if row == 0:
            cell.set_facecolor("#DDEBF7"); cell.set_text_props(weight="bold")
        elif row == 2:
            cell.set_facecolor("#F7F7F7")
    ax.text(0.0, 0.24, "Data status: accepted terminal P0-window CW/CSTR ROP only.\nNo pulse calculation or experimental result is shown.", fontsize=6.2, weight="bold")
    ax.text(0.0, 0.12, "Margin = named electronic H$_2$* share / largest secondary direct-NH share.\nIt ranks contrast for testing; it is not a confidence interval.", fontsize=6.05)
    ax.text(0.0, 0.005, "Figure role: candidate design figure for the next manuscript version;\ndo not insert into the current result section until independent tests exist.", fontsize=6.0, style="italic")
    panel(ax, "D")


def make_figure(frame: pd.DataFrame, output: Path) -> None:
    # A double-column design figure: fixed margins keep the two tabular panels
    # readable at the final manuscript scale.  ``constrained_layout`` is not
    # used because Matplotlib can collapse table-containing axes.
    fig = plt.figure(figsize=(12.0, 8.55))
    grid = fig.add_gridspec(2, 2, left=0.06, right=0.985, top=0.94, bottom=0.075,
                            hspace=0.62, wspace=0.32,
                            height_ratios=(1.0, 1.03), width_ratios=(1.0, 1.04))
    draw_entry_hierarchy(fig.add_subplot(grid[0, 0]), frame)
    draw_logic(fig.add_subplot(grid[0, 1]), frame)
    draw_measurement_matrix(fig.add_subplot(grid[1, 0]))
    draw_decision_ledger(fig.add_subplot(grid[1, 1]), frame)
    fig.text(0.5, 0.018, "Graph-informed model-discrimination atlas. Proposed verification design; not experimental evidence.", ha="center", fontsize=7.0)
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / "Figure_19_graph_informed_model_discrimination_atlas.png", dpi=600)
    fig.savefig(output / "Figure_19_graph_informed_model_discrimination_atlas.pdf")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=FIGURES)
    args = parser.parse_args()
    configure()
    frame = selected_frame()
    args.output.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output / "Figure_19_condition_selection_source.csv", index=False)
    make_figure(frame, args.output)
    print(frame[["code", "en_td", "n2_fraction", "named_h2star_pct", "largest_secondary_pct", "dominance_margin"]].to_csv(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
