#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate the two additional, main-text-only figures for the Hong study.

The script intentionally reads only the validated CSTR/CW calculations and
the closed-system regression.  It neither reads nor represents the rejected
pulsed-afterglow runs.  Figure 1 documents the P0 numerical basis; Figure 6
provides a symmetry-aware, quantitative view of NH-source topology.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd

from _paths import LEGACY_HONG_ROOT

TOOL_DIR = Path(__file__).resolve().parent
PROJECT = LEGACY_HONG_ROOT
RUNS = PROJECT / "Reproduction" / "2017Hong" / "方案2_脉冲能效地图" / "build_pulse" / "runs"
FIGURES = PROJECT / "analysis" / "figures_20260829"

COLORS = {
    "blue": "#0072B2", "sky": "#56B4E9", "green": "#009E73",
    "orange": "#E69F00", "vermillion": "#D55E00", "purple": "#CC79A7",
    "gray": "#777777", "black": "#000000",
}


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 9,
        "axes.titlesize": 9, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 6.7, "axes.linewidth": 0.7, "pdf.fonttype": 42,
        "ps.fonttype": 42, "savefig.transparent": False,
    })


def panel(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.15, 1.08, letter, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top")


def save(fig: plt.Figure, output: Path, name: str) -> None:
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(output / f"{name}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def read_fortran_csv(path: Path) -> pd.DataFrame:
    """Read a fixed-width-friendly CSV emitted by the Fortran driver."""
    return pd.read_csv(path, skipinitialspace=True)


def make_figure1(output: Path) -> None:
    closed = read_fortran_csv(RUNS / "p0_closed_rhs_regression" /
                              "pulse_cycles_p0_closed_rhs_regression.csv")
    # The RHS-patched CSTR regression is the P0 reference: it reaches three
    # consecutive cycles below both 1e-3 residual thresholds at cycle 55.
    cstr_path = RUNS / "p0_cstr_rhs_tau1e-2" / "pulse_cycles_p0_cstr_rhs_tau1e-2.csv"
    cstr = read_fortran_csv(cstr_path)
    dense = pd.read_csv(output / "Figure_3_dense_EN_N2_source_data.csv")

    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.15), constrained_layout=True)

    # A. A deliberately compact reactor schematic that fixes the numerical scope.
    ax = axes[0, 0]
    ax.axis("off")
    vessel = FancyBboxPatch((0.28, 0.22), 0.44, 0.56, boxstyle="round,pad=0.035",
                            facecolor="#EAF3F8", edgecolor=COLORS["blue"], lw=1.4)
    ax.add_patch(vessel)
    ax.add_patch(FancyArrowPatch((0.03, 0.50), (0.28, 0.50), arrowstyle="-|>",
                                 mutation_scale=12, lw=1.2, color=COLORS["green"]))
    ax.add_patch(FancyArrowPatch((0.72, 0.50), (0.97, 0.50), arrowstyle="-|>",
                                 mutation_scale=12, lw=1.2, color=COLORS["vermillion"]))
    ax.add_patch(FancyArrowPatch((0.50, 0.98), (0.50, 0.78), arrowstyle="-|>",
                                 mutation_scale=12, lw=1.2, color=COLORS["orange"]))
    ax.text(0.50, 0.62, "0D plasma CSTR", ha="center", va="center", weight="bold")
    ax.text(0.50, 0.48, r"Hong 2017 kinetics", ha="center", va="center")
    ax.text(0.50, 0.36, r"$\tau$ = 10 ms; $T_g$ = 300 K", ha="center", va="center", fontsize=7)
    ax.text(0.02, 0.59, r"feed: N$_2$ + H$_2$", color=COLORS["green"], fontsize=7)
    ax.text(0.76, 0.59, r"outlet", color=COLORS["vermillion"], fontsize=7)
    ax.text(0.53, 0.93, r"CW field; fixed $n_e$", color=COLORS["orange"], fontsize=7)
    ax.text(0.50, 0.07, "P0: terminal three-cycle steady window", ha="center", fontsize=7)
    ax.set_title("Validated steady-state calculation")
    panel(ax, "A")

    ax = axes[0, 1]
    ax.plot(closed["cycle"], closed["NH3_end"], color=COLORS["gray"], lw=1.1,
            label="closed-system regression")
    ax.plot(cstr["cycle"], cstr["NH3_end"], color=COLORS["blue"], lw=1.5,
            label=r"CSTR, $\tau$ = 10 ms")
    ax.set_yscale("log")
    ax.set_xlabel("P0 cycle")
    ax.set_ylabel(r"Cycle-end [NH$_3$] (cm$^{-3}$)")
    ax.set_title("Open flow establishes a periodic steady state")
    ax.legend(frameon=False, loc="best")
    panel(ax, "B")

    ax = axes[1, 0]
    c = cstr.loc[cstr["cycle"] > 1]
    ax.semilogy(c["cycle"], c["rel_state_max"], color=COLORS["blue"], marker="o", ms=2.1,
                lw=1.0, label="state residual")
    ax.semilogy(c["cycle"], c["rel_dNH3_cycle"], color=COLORS["vermillion"], marker="s", ms=2.0,
                lw=1.0, label=r"NH$_3$ residual")
    ax.axhline(1e-3, color=COLORS["black"], ls="--", lw=0.8, label="P0 threshold")
    ax.set_xlabel("P0 cycle")
    ax.set_ylabel("Relative cycle-to-cycle change")
    ax.set_ylim(1e-5, 1.5)
    ax.set_title("CSTR P0 convergence at 120 Td, $x_{N_2}$ = 0.3333")
    ax.legend(frameon=False, ncol=1, loc="upper right")
    panel(ax, "C")

    ax = axes[1, 1]
    ok = dense["run_status"].isin(["success", "retry_success"])
    if not bool(ok.all()):
        raise ValueError("Dense map contains a non-validated run.")
    retry = dense["run_status"].eq("retry_success")
    quality = np.maximum(dense["rel_state_max"], dense["rel_dnh3_cycle"])
    dots = ax.scatter(dense.loc[~retry, "cycles_used"], quality[~retry],
                      c=np.log10(dense.loc[~retry, "nh3_final_cm-3"].clip(lower=1e-20)),
                      cmap="cividis", s=20, edgecolor="none", label="first-pass success")
    ax.scatter(dense.loc[retry, "cycles_used"], quality[retry], marker="x", s=38,
               linewidth=1.1, color=COLORS["vermillion"], label="retry success")
    ax.axhline(1e-3, color=COLORS["black"], ls="--", lw=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("Cycles used")
    ax.set_ylabel("max(P0 residuals)")
    ax.set_title("All 108 dense-scan cases pass P0")
    ax.legend(frameon=False, loc="lower right")
    cb = fig.colorbar(dots, ax=ax, fraction=0.047, pad=0.04)
    cb.set_label(r"log$_{10}$[NH$_3$] (cm$^{-3}$)")
    panel(ax, "D")
    save(fig, output, "Figure_1_CSTR_P0_validation")


def source_entropy(data: pd.DataFrame) -> pd.Series:
    cols = ["named_h2star_pct", "rydberg_h2star_pct", "n2d_h2_pct", "n2p_h2_pct",
            "vib_h2_pct", "association_pct"]
    fractions = data[cols].to_numpy(float) / 100.0
    fractions = np.clip(fractions, 1e-300, 1.0)
    return pd.Series(-(fractions * np.log2(fractions)).sum(axis=1), index=data.index)


def make_figure6(output: Path) -> None:
    data = pd.read_csv(output / "Figure_3_dense_EN_N2_source_data.csv")
    n2s, ens = sorted(data.n2_fraction.unique()), sorted(data.en_td.unique())
    pvt = lambda col: data.pivot(index="en_td", columns="n2_fraction", values=col).loc[ens, n2s]
    entropy = source_entropy(data)
    data = data.assign(source_entropy_bits=entropy)
    channels = ["rydberg_h2star_pct", "n2d_h2_pct", "n2p_h2_pct", "vib_h2_pct", "association_pct"]
    short = {"rydberg_h2star_pct": "Rydberg-H2*", "n2d_h2_pct": "N(2D)+H2",
             "n2p_h2_pct": "N(2P)+H2", "vib_h2_pct": "vibrational H2",
             "association_pct": "H+N+M"}
    runner = data[channels].idxmax(axis=1).map(short)
    data = data.assign(runner_up=runner)

    fig = plt.figure(figsize=(11.6, 5.25), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, width_ratios=(1.60, 0.93, 0.93), height_ratios=(1, 1))

    # A. State-ledger: safe group-theoretical content, without inventing g/u parity.
    ax = fig.add_subplot(grid[:, 0]); ax.axis("off")
    ax.set_title("Symmetry ledger for the NH-source manifold", loc="left", pad=8)
    ax.text(0.00, 0.94, r"D$_{\infty h}$ is applicable to isolated H$_2$ states;", fontsize=7.4)
    ax.text(0.00, 0.89, "electron-impact kinetics require collision dynamics in addition", fontsize=7.4)
    ax.text(0.00, 0.84, "to state symmetry.  The mechanism abbreviations omit g/u and +/-.", fontsize=7.4)
    headers = ["family", r"$S,\Lambda$", "symmetry status", "NH role"]
    rows = [
        [r"H$_2$(v=0–3)", r"X $^1\Sigma_g^+$", "complete term", "vibrational"],
        [r"B3SIG / A3SIG", r"$^3\Sigma$", "g/u, +/- omitted", "named H$_2$*"],
        [r"B1SIG", r"$^1\Sigma$", "g/u, +/- omitted", "named H$_2$*"],
        [r"C3PI", r"$^3\Pi$", "g/u, +/- omitted", "named H$_2$*"],
        [r"Rydberg-H$_2$*", "unspecified", "term unresolved", "secondary"],
        [r"N($^2$P), N($^2$D)", "atomic", r"not molecular D$_{\infty h}$", "secondary"],
    ]
    table = ax.table(cellText=rows, colLabels=headers, cellLoc="left", colLoc="left",
                     colWidths=[0.28, 0.16, 0.34, 0.22], bbox=[0.00, 0.19, 1.0, 0.60])
    table.auto_set_font_size(False); table.set_fontsize(6.6)
    for (r, c), cell in table.get_celld().items():
        cell.set_linewidth(0.25)
        if r == 0:
            cell.set_facecolor("#DDEBF7"); cell.set_text_props(weight="bold")
        elif r in (2, 3, 4, 5):
            cell.set_facecolor("#F8F8F8")
    ax.text(0.00, 0.11, "Interpretation rule: labels support classification by spin and", fontsize=7.2, weight="bold")
    ax.text(0.00, 0.06, "orbital projection only; they cannot alone impose a selection rule", fontsize=7.2)
    ax.text(0.00, 0.01, "or prove a rate ordering.  Flux integration supplies that ranking.", fontsize=7.2)
    panel(ax, "A")

    # B. Entropy makes continuous redistribution mathematically visible.
    ax = fig.add_subplot(grid[0, 1])
    im = ax.imshow(pvt("source_entropy_bits"), origin="lower", aspect="auto", cmap="magma", vmin=0, vmax=0.85)
    ax.set_xticks(range(len(n2s)), [f"{x:.1f}" for x in n2s])
    ax.set_yticks(range(0, len(ens), 2), [str(ens[i]) for i in range(0, len(ens), 2)])
    ax.set_xlabel(r"N$_2$ inlet fraction"); ax.set_ylabel("E/N (Td)")
    ax.set_title("Secondary-source diversity")
    cb = fig.colorbar(im, ax=ax, fraction=0.047, pad=0.04); cb.set_label("Shannon entropy (bit)")
    panel(ax, "B")

    # C. Runner-up channel field, a categorical pathway-switch map.
    ax = fig.add_subplot(grid[0, 2])
    cats = ["Rydberg-H2*", "N(2D)+H2", "N(2P)+H2", "vibrational H2", "H+N+M"]
    cmap = mpl.colors.ListedColormap([COLORS["purple"], COLORS["sky"], COLORS["orange"], COLORS["green"], COLORS["vermillion"]])
    arr = data.assign(code=pd.Categorical(data.runner_up, categories=cats).codes).pivot(index="en_td", columns="n2_fraction", values="code").loc[ens, n2s]
    ax.imshow(arr, origin="lower", aspect="auto", cmap=cmap, vmin=-0.5, vmax=4.5)
    ax.set_xticks(range(len(n2s)), [f"{x:.1f}" for x in n2s])
    ax.set_yticks(range(0, len(ens), 2), [str(ens[i]) for i in range(0, len(ens), 2)])
    ax.set_xlabel(r"N$_2$ inlet fraction"); ax.set_ylabel("E/N (Td)")
    ax.set_title("Identity of the largest secondary source")
    handles = [mpl.patches.Patch(color=cmap(i), label=cats[i]) for i in range(len(cats))]
    ax.legend(handles=handles, fontsize=5.65, frameon=False, loc="upper center", ncol=2,
              bbox_to_anchor=(0.5, -0.27), columnspacing=0.75, handlelength=1.0)
    panel(ax, "C")

    # D. A flux-topology plot for the two composition extremes, normalized to total NH production.
    ax = fig.add_subplot(grid[1, 1:])
    use_channels = [("named_h2star_pct", "named electronic H$_2$*", COLORS["blue"]),
                    ("rydberg_h2star_pct", "Rydberg-H$_2$*", COLORS["purple"]),
                    ("n2p_h2_pct", "N($^2$P)+H$_2$", COLORS["orange"]),
                    ("association_pct", "H+N+M", COLORS["vermillion"])]
    for n2, ls in [(0.1, "-"), (0.9, "--")]:
        sub = data.loc[np.isclose(data.n2_fraction, n2)].sort_values("en_td")
        for col, name, color in use_channels:
            ax.plot(sub.en_td, sub[col], color=color, ls=ls, lw=1.2, marker="o", ms=2.0,
                    label=f"{name}; $x_{{N_2}}$={n2:.1f}")
    ax.set_yscale("log"); ax.set_ylim(1e-5, 120); ax.set_xlim(min(ens), max(ens))
    ax.set_xticks(ens[::2]); ax.set_xlabel("E/N (Td)")
    ax.set_ylabel("Integrated direct NH source (%)")
    ax.set_title("Continuous redistribution without loss of named-H$_2$* dominance")
    # Explain visual encoding once, avoiding a twelve-entry legend.
    handles = [mpl.lines.Line2D([], [], color=c, lw=1.5, label=n) for _, n, c in use_channels]
    handles += [mpl.lines.Line2D([], [], color="black", lw=1.4, ls="-", label=r"$x_{N_2}=0.1$"),
                mpl.lines.Line2D([], [], color="black", lw=1.4, ls="--", label=r"$x_{N_2}=0.9$")]
    ax.legend(handles=handles, frameon=False, ncol=2, loc="upper right")
    panel(ax, "D")
    save(fig, output, "Figure_6_symmetry_aware_NH_path_topology")

    # Reproducible source table for the new quantitative panels.
    data[["en_td", "n2_fraction", "source_entropy_bits", "runner_up", "named_h2star_pct",
          "rydberg_h2star_pct", "n2p_h2_pct", "association_pct"]].to_csv(
              output / "Figure_6_source_data.csv", index=False)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=FIGURES)
    args = parser.parse_args()
    configure()
    make_figure1(args.output)
    make_figure6(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
