#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Draw a source-controlled graphical abstract for the Hong–Fe(110) paper.

This is a conceptual scope diagram. It reports no new numerical quantity and
does not depict either gas or surface route as dominant.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle


PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = PROJECT / "analysis" / "graphical_abstract_20260901_r1"


def molecule(ax, x: float, y: float, labels: tuple[str, ...], colors: tuple[str, ...], scale: float = 1.0) -> None:
    offsets = ((0, 0), (0.055 * scale, 0.03 * scale), (0.055 * scale, -0.03 * scale), (-0.055 * scale, 0))
    for index, (label, color) in enumerate(zip(labels, colors)):
        dx, dy = offsets[index]
        ax.add_patch(Circle((x + dx, y + dy), 0.022 * scale, facecolor=color, edgecolor="#4D4D4D", lw=0.5, zorder=4))
        ax.text(x + dx, y + dy, label, ha="center", va="center", fontsize=5.2 * scale, color="white", zorder=5)
    for index in range(1, len(labels)):
        dx, dy = offsets[index]
        ax.plot([x, x + dx], [y, y + dy], color="#606060", lw=0.8, zorder=3)


def arrow(ax, start, end, color="#5A5A5A", text="", bend=0.0):
    ax.add_patch(FancyArrowPatch(start, end, connectionstyle=f"arc3,rad={bend}", arrowstyle="-|>", mutation_scale=11,
                                 linewidth=1.25, color=color, zorder=2))
    if text:
        ax.text((start[0] + end[0]) / 2, (start[1] + end[1]) / 2 + (0.035 if bend >= 0 else -0.035), text,
                ha="center", va="center", fontsize=7.1, color="#303030")


def draw(output: Path) -> None:
    mpl.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, ax = plt.subplots(figsize=(12.0, 5.4))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.5, 0.94, "Pathway-resolved plasma ammonia synthesis: from gas-phase baseline to testable gas–surface competition", ha="center", va="center", fontsize=14, fontweight="bold")

    ax.add_patch(FancyBboxPatch((0.035, 0.20), 0.16, 0.54, boxstyle="round,pad=0.015", facecolor="#F5F8FB", edgecolor="#5E7E9E", lw=1.0))
    ax.text(0.115, 0.68, "Feed", ha="center", va="center", fontsize=10.5, fontweight="bold")
    molecule(ax, 0.090, 0.52, ("N", "N"), ("#1F77B4", "#1F77B4"), 1.2)
    molecule(ax, 0.090, 0.34, ("H", "H"), ("#D55E00", "#D55E00"), 1.2)
    ax.text(0.115, 0.25, "controlled N$_2$/H$_2$ inlet", ha="center", fontsize=7.1)

    ax.add_patch(FancyBboxPatch((0.255, 0.20), 0.22, 0.54, boxstyle="round,pad=0.015", facecolor="#EAF3FA", edgecolor="#3C7DA6", lw=1.0))
    ax.text(0.365, 0.68, "0D plasma-CSTR baseline", ha="center", va="center", fontsize=10.5, fontweight="bold")
    ax.add_patch(Rectangle((0.285, 0.38), 0.16, 0.17, facecolor="white", edgecolor="#4F81A8", lw=0.9))
    for xpos in (0.30, 0.34, 0.38, 0.42):
        ax.plot([xpos, xpos + 0.012], [0.395, 0.535], color="#72B7D2", lw=1.4)
    ax.text(0.365, 0.46, "continuous E/N\nCSTR + P0", ha="center", va="center", fontsize=8.2)
    ax.text(0.365, 0.27, "gas NH-entry fractions\nand downstream NH$_x$ turnover", ha="center", fontsize=7.3)

    ax.add_patch(FancyBboxPatch((0.535, 0.48), 0.19, 0.26, boxstyle="round,pad=0.015", facecolor="#FFF4D6", edgecolor="#B8872B", lw=1.0))
    ax.text(0.63, 0.68, "Gas route", ha="center", va="center", fontsize=10, fontweight="bold")
    molecule(ax, 0.585, 0.57, ("N",), ("#1F77B4",), 1.0)
    molecule(ax, 0.650, 0.57, ("H*", "H*"), ("#56B4E9", "#56B4E9"), 0.9)
    arrow(ax, (0.605, 0.565), (0.675, 0.565), color="#C57B00")
    ax.text(0.63, 0.505, "state-resolved gas NH entry", ha="center", fontsize=7.2)

    ax.add_patch(FancyBboxPatch((0.535, 0.20), 0.19, 0.20, boxstyle="round,pad=0.015", facecolor="#EAF5E5", edgecolor="#5D9155", lw=1.0))
    ax.text(0.63, 0.35, "Fe(110) structural branch", ha="center", va="center", fontsize=9.4, fontweight="bold")
    for xpos in (0.565, 0.605, 0.645, 0.685):
        ax.add_patch(Circle((xpos, 0.265), 0.016, facecolor="#666666", edgecolor="#333333", lw=0.4))
    ax.text(0.63, 0.215, "ER and LH first N–H alternatives", ha="center", fontsize=7.0)

    ax.add_patch(FancyBboxPatch((0.765, 0.20), 0.13, 0.54, boxstyle="round,pad=0.015", facecolor="#FCEAEA", edgecolor="#B45B5B", lw=1.0))
    ax.text(0.83, 0.67, "Comparison gate", ha="center", va="center", fontsize=9.7, fontweight="bold")
    ax.text(0.83, 0.54, r"$F_s=\beta E_{bulk}$", ha="center", fontsize=10.2, fontweight="bold")
    ax.text(0.83, 0.41, "shared reactor state\nsite conservation\nROP fractions", ha="center", fontsize=7.5)
    ax.text(0.83, 0.26, "no cross-boundary\nflux claim before pass", ha="center", fontsize=6.9, style="italic")

    ax.add_patch(FancyBboxPatch((0.925, 0.20), 0.055, 0.54, boxstyle="round,pad=0.015", facecolor="#F1E9FB", edgecolor="#8064A2", lw=1.0))
    ax.text(0.952, 0.68, "NH$_3$", ha="center", va="center", fontsize=10.5, fontweight="bold")
    molecule(ax, 0.946, 0.48, ("N", "H", "H", "H"), ("#1F77B4", "#D55E00", "#D55E00", "#D55E00"), 0.65)
    ax.text(0.952, 0.27, "product +\nindependent observables", ha="center", fontsize=6.6)

    arrow(ax, (0.195, 0.47), (0.255, 0.47), color="#477DA3")
    arrow(ax, (0.475, 0.54), (0.535, 0.60), color="#C57B00")
    arrow(ax, (0.475, 0.38), (0.535, 0.30), color="#568C53")
    arrow(ax, (0.725, 0.60), (0.765, 0.57), color="#A65A5A")
    arrow(ax, (0.725, 0.30), (0.765, 0.37), color="#A65A5A")
    arrow(ax, (0.895, 0.48), (0.925, 0.48), color="#8064A2")
    ax.text(0.5, 0.07, "Graphical abstract. The gas-phase Hong conclusion is quantitatively reported within its CW/CSTR boundary; the Fe(110) branch supplies structural alternatives and a falsifiable harmonization protocol, not an unverified surface-flux ranking.", ha="center", va="center", fontsize=7.5, style="italic")
    fig.savefig(output / "Graphical_abstract_scope.png", dpi=600, bbox_inches="tight")
    fig.savefig(output / "Graphical_abstract_scope.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Draw graphical abstract for Hong–Fe(110) manuscript")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(); output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    draw(output)
    (output / "README.md").write_text(
        "# Graphical abstract\n\n"
        "Conceptual scope diagram for the English manuscript. It contains no new data and does not rank gas and surface fluxes. "
        "The attempted AI schematic generation was unavailable because its connection failed; this reproducible vector-like fallback is drawn from the audited model boundary.\n",
        encoding="utf-8")
    print(f"output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
