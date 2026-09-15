#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared visual language for the manuscript publication figures.

The module deliberately contains presentation rules only.  It neither changes
source data nor infers new kinetics.  Plot scripts may retain their local
layouts while obtaining a common, colour-blind-safe visual hierarchy.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


# Stable semantic colours used throughout the Hong/Fe(110) manuscript.
PALETTE = {
    "named": "#0072B2",       # named electronic-H2* direct-NH entry
    "secondary": "#CC79A7",   # Rydberg/special excited competitor
    "excited_n": "#56B4E9",   # excited-N competitor
    "turnover": "#009E73",    # productive NHx hydrogenation
    "loss": "#D55E00",        # ion-mediated product loss
    "transition": "#E69F00",  # transition/uncertain regime
    "reference": "#6B7280",   # neutral reference/background
    "soft_blue": "#DCEAF7",
    "soft_green": "#E7F2DF",
    "soft_orange": "#FCE7D6",
    "soft_red": "#F8DFDF",
}


def apply_theme() -> None:
    """Apply an accessible, compact double-column publication theme."""
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "font.size": 8.0,
        "axes.labelsize": 8.5,
        "axes.titlesize": 8.8,
        "xtick.labelsize": 7.0,
        "ytick.labelsize": 7.0,
        "legend.fontsize": 6.5,
        "axes.linewidth": 0.65,
        "lines.linewidth": 1.35,
        "lines.markersize": 4.6,
        "axes.titlepad": 5.0,
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.transparent": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "mathtext.default": "regular",
    })


def refine_figure(fig: plt.Figure, *, keep_boxed: bool = False) -> None:
    """Reduce chart junk without modifying plotted values or annotations."""
    for ax in fig.axes:
        # Colourbars have no conventional data spines.
        if ax.get_label() == "<colorbar>":
            ax.tick_params(length=2.2, width=0.6)
            continue
        ax.tick_params(direction="out", length=3.0, width=0.65, pad=2.4)
        if not keep_boxed:
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        # Preserve grids deliberately requested by an individual script, but
        # make them supporting structure rather than a competing visual layer.
        for line in list(ax.get_xgridlines()) + list(ax.get_ygridlines()):
            line.set_color("#CBD5E1")
            line.set_alpha(0.36)
            line.set_linewidth(0.45)
        legend = ax.get_legend()
        if legend is not None:
            legend.set_frame_on(False)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.11, 1.07, label, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="right", clip_on=False)


def export_figure(fig: plt.Figure, output: Path, stem: str, *, dpi: int = 600,
                  keep_boxed: bool = False) -> tuple[Path, Path]:
    """Apply final styling and export matching vector and high-DPI files."""
    output.mkdir(parents=True, exist_ok=True)
    refine_figure(fig, keep_boxed=keep_boxed)
    pdf = output / f"{stem}.pdf"
    png = output / f"{stem}.png"
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.035)
    fig.savefig(png, dpi=dpi, bbox_inches="tight", pad_inches=0.035)
    return pdf, png
