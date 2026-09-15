#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Render main-text Figure 10 from the P6 local kinetic-control table."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


TOOL_DIR = Path(__file__).resolve().parent
PROJECT = TOOL_DIR.parent
DEFAULT_INPUT = PROJECT / "analysis" / "p6_local_reaction_control_20260830" / "P6_local_sensitivities.csv"
DEFAULT_CASES = PROJECT / "analysis" / "p6_local_reaction_control_20260830" / "P6_control_case_results.csv"
DEFAULT_OUTPUT = PROJECT / "analysis" / "figures_20260829"

ROLE_COLORS = {
    "primary NH source": "#0072B2", "secondary NH source": "#CC79A7",
    "associative NH source": "#E69F00", "NH3 formation": "#009E73",
    "NH3 ionization loss": "#D55E00", "ion-mediated NH3 loss": "#8C564B",
}


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8.5,
        "axes.titlesize": 8.8, "xtick.labelsize": 7, "ytick.labelsize": 6.8,
        "legend.fontsize": 6.3, "axes.linewidth": 0.7, "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def panel(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.15, 1.07, letter, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top")


def short_label(label: str) -> str:
    return (label.replace("Named electronic H2*", "named electronic H$_2$*")
            .replace("Rydberg-H2*", "Rydberg-H$_2$*")
            .replace("N(2D,2P)", "N($^2$D,$^2$P)")
            .replace("H + N + M", "H + N + M")
            .replace("NH2", "NH$_2$").replace("NH3+", "NH$_3^+$")
            .replace("NH3", "NH$_3$").replace("H2", "H$_2$"))


def plot_tornado(ax: plt.Axes, data: pd.DataFrame, cases: pd.DataFrame, representative: str, output: str,
                 title: str, letter: str) -> None:
    sub = data.loc[(data.representative.eq(representative)) & data.output.eq(output)].copy()
    sub = sub.reindex(sub.sensitivity.abs().sort_values().index)
    colors = sub.role.map(ROLE_COLORS).fillna("#777777")
    y = np.arange(len(sub))
    ax.barh(y, sub.sensitivity, color=colors, edgecolor="none")
    ax.axvline(0, color="black", lw=0.7)
    ax.set_yticks(y, [short_label(value) for value in sub.control_label])
    ax.set_xlabel("Local sensitivity, S = d ln(y) / d ln(k)")
    ax.set_title(title)
    ax.grid(axis="x", color="#DDDDDD", lw=0.45, zorder=0)
    ax.set_axisbelow(True)
    baseline = cases.loc[(cases.representative.eq(representative)) & cases.control.eq("baseline")].iloc[0]
    ax.text(0.02, 0.03,
            rf"baseline $\dot n_{{NH_3}}={baseline['NH3_productivity_cm-3s-1']:.2e}$" + "\n" +
            rf"$D/P={baseline['NH3_loss_to_formation']:.4f}$",
            transform=ax.transAxes, va="bottom", fontsize=6.1,
            bbox={"facecolor": "white", "edgecolor": "#AAAAAA", "alpha": 0.92, "pad": 1.2})
    panel(ax, letter)


def make_figure(data: pd.DataFrame, cases: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.7, 7.35), constrained_layout=True)
    plot_tornado(axes[0, 0], data, cases, "productivity_maximum", "NH3_productivity_cm-3s-1",
                 r"Production control at the NH$_3$ maximum", "A")
    plot_tornado(axes[0, 1], data, cases, "loss_dominated", "NH3_productivity_cm-3s-1",
                 r"Production control in the loss-dominated state", "B")
    plot_tornado(axes[1, 0], data, cases, "loss_dominated", "NH3_loss_to_formation",
                 r"Control of NH$_3$ destruction / formation", "C")
    axes[1, 0].text(0.52, 0.92, "Near $D/P=1$, formation and loss\nco-vary; productivity is more diagnostic.",
                     transform=axes[1, 0].transAxes, fontsize=6.0, va="top")

    ax = axes[1, 1]
    sub = data.loc[data.output.eq("named_H2star_source_share")].copy()
    labels = sub.drop_duplicates("control").set_index("control").loc[:, "control_label"]
    order = (sub.groupby("control").sensitivity.apply(lambda x: x.abs().max())
             .sort_values().index.tolist())
    codes = {"productivity_maximum": ("NH$_3$ maximum", "#0072B2", "o"),
             "loss_dominated": ("loss-dominated", "#D55E00", "s")}
    y = np.arange(len(order))
    for representative, (label, color, marker) in codes.items():
        values = sub.loc[sub.representative.eq(representative)].set_index("control").loc[order, "sensitivity"]
        ax.plot(values, y, marker=marker, ms=4.2, lw=1.05, color=color, label=label)
    ax.axvline(0, color="black", lw=0.7)
    ax.set_yticks(y, [short_label(labels[item]) for item in order])
    ax.set_xlabel(r"$S$ of named-H$_2$* NH-source share")
    ax.set_title("Robustness of named-H$_2$* source dominance")
    ax.grid(axis="x", color="#DDDDDD", lw=0.45); ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower right")
    panel(ax, "D")

    handles = [mpl.patches.Patch(color=color, label=role) for role, color in ROLE_COLORS.items()]
    fig.legend(handles=handles, ncol=3, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.095))
    fig.text(0.5, -0.155,
             "Local sensitivities use symmetric 1/1.10 and 1.10 rate multipliers. Finite pathway disabling (Figure 2) and differential local control answer distinct causal questions.",
             ha="center", fontsize=7)
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / "Figure_10_local_kinetic_control.pdf", bbox_inches="tight")
    fig.savefig(output / "Figure_10_local_kinetic_control.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot Hong Figure 10 local kinetic-control analysis")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    data = pd.read_csv(args.input)
    cases = pd.read_csv(args.cases)
    expected = 2 * 8 * 3
    if len(data) != expected:
        raise ValueError(f"Expected {expected} sensitivity rows, found {len(data)}")
    configure()
    if len(cases) != 34:
        raise ValueError(f"Expected 34 P6 case rows, found {len(cases)}")
    make_figure(data, cases, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
