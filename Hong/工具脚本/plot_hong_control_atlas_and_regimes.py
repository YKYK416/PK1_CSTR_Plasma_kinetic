#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Render manuscript Figures 11 and 12 from P7 and the accepted 108-point data."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


TOOL_DIR = Path(__file__).resolve().parent
PROJECT = TOOL_DIR.parent
P7 = PROJECT / "analysis" / "p7_control_switching_atlas_20260830" / "P7_control_sensitivities.csv"
P5 = PROJECT / "analysis" / "figures_20260829" / "Figure_8_9_NHx_CSTR_summary.csv"
DENSE = PROJECT / "analysis" / "figures_20260829" / "Figure_3_dense_EN_N2_source_data.csv"
OUTPUT = PROJECT / "analysis" / "figures_20260829"

CONTROLS = ("NH_to_NH3_association", "NH2_hydrogenation", "proton_NH3_ionization")
CONTROL_LABELS = {
    "NH_to_NH3_association": r"NH + H$_2$ + M $\rightarrow$ NH$_3$ + M",
    "NH2_hydrogenation": r"H + NH$_2$ + M $\rightarrow$ NH$_3$ + M",
    "proton_NH3_ionization": r"H$^+$ + NH$_3$ $\rightarrow$ NH$_3^+$ + H",
}
CONTROL_COLORS = {"NH_to_NH3_association": "#0072B2", "NH2_hydrogenation": "#009E73",
                  "proton_NH3_ionization": "#D55E00"}
REGIMES = ("formation-starved", "productive hydrogenation", "transition turnover", "ion-loss cancellation")
REGIME_COLORS = {"formation-starved": "#9E9E9E", "productive hydrogenation": "#009E73",
                 "transition turnover": "#E69F00", "ion-loss cancellation": "#D55E00"}


def configure() -> None:
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8.5,
                         "axes.titlesize": 9, "xtick.labelsize": 7, "ytick.labelsize": 7,
                         "legend.fontsize": 6.6, "axes.linewidth": 0.7, "pdf.fonttype": 42,
                         "ps.fonttype": 42})


def panel(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.14, 1.15, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top", clip_on=False)


def setup_axes(ax: plt.Axes, n2s: list[float], ens: list[int]) -> None:
    ax.set_xticks(range(len(n2s)), [f"{x:.1f}" for x in n2s])
    ax.set_yticks(range(len(ens)), [str(value) for value in ens])
    ax.set_xlabel(r"N$_2$ inlet fraction")
    ax.set_ylabel("E/N (Td)")
    ax.set_xticks(np.arange(-0.5, len(n2s), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(ens), 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.42, alpha=0.82)
    ax.tick_params(which="minor", bottom=False, left=False)


def sensitivity_matrix(data: pd.DataFrame, control: str, ens: list[int], n2s: list[float]) -> np.ndarray:
    return (data.loc[data.control.eq(control)].pivot(index="en_td", columns="n2_fraction",
            values="sensitivity_productivity").reindex(index=ens, columns=n2s).to_numpy(float))


def draw_high_loss_frames(ax: plt.Axes, turnover: pd.DataFrame, ens: list[int], n2s: list[float]) -> None:
    sampled = turnover.loc[turnover.en_td.isin(ens) & turnover.n2_fraction.isin(n2s)]
    for row in sampled.loc[sampled.NH3_loss_to_formation.ge(0.9)].itertuples():
        ax.add_patch(Rectangle((n2s.index(row.n2_fraction) - 0.5, ens.index(row.en_td) - 0.5), 1, 1,
                               fill=False, edgecolor="black", lw=1.0, zorder=6))


def make_figure11(sensitivity: pd.DataFrame, turnover: pd.DataFrame, output: Path) -> None:
    ens, n2s = sorted(sensitivity.en_td.unique()), sorted(sensitivity.n2_fraction.unique())
    max_abs = float(np.max(np.abs(sensitivity.sensitivity_productivity)))
    fig, axes = plt.subplots(2, 2, figsize=(9.7, 6.65), constrained_layout=True)
    for ax, control, letter in zip(axes.flat[:3], CONTROLS, "ABC"):
        image = ax.imshow(sensitivity_matrix(sensitivity, control, ens, n2s), origin="lower", aspect="auto",
                          cmap="RdBu_r", vmin=-max_abs, vmax=max_abs)
        setup_axes(ax, n2s, ens); ax.set_title(CONTROL_LABELS[control])
        matrix = sensitivity_matrix(sensitivity, control, ens, n2s)
        for row, en in enumerate(ens):
            for col, n2 in enumerate(n2s):
                value = matrix[row, col]
                text_color = "white" if abs(value) > 0.55 * max_abs else "black"
                ax.text(col, row, f"{value:.2f}", ha="center", va="center", fontsize=5.6, color=text_color)
        draw_high_loss_frames(ax, turnover, ens, n2s)
        colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        colorbar.set_label(r"$S$ of NH$_3$ productivity")
        panel(ax, letter)

    ax = axes[1, 1]
    largest = (sensitivity.assign(abs_s=lambda x: x.sensitivity_productivity.abs())
               .sort_values("abs_s").groupby(["en_td", "n2_fraction"], as_index=False).tail(1))
    codes = {control: number for number, control in enumerate(CONTROLS)}
    matrix = (largest.assign(code=largest.control.map(codes)).pivot(index="en_td", columns="n2_fraction",
              values="code").reindex(index=ens, columns=n2s).to_numpy(float))
    cmap = mpl.colors.ListedColormap([CONTROL_COLORS[control] for control in CONTROLS])
    ax.imshow(matrix, origin="lower", aspect="auto", cmap=cmap, vmin=-0.5, vmax=2.5)
    negative = largest.loc[largest.sensitivity_productivity.lt(0)]
    for _, row in negative.iterrows():
        ax.text(n2s.index(row.n2_fraction), ens.index(row.en_td), "−", color="white", ha="center", va="center",
                fontsize=12, fontweight="bold")
    setup_axes(ax, n2s, ens); ax.set_title("Largest absolute local productivity control")
    draw_high_loss_frames(ax, turnover, ens, n2s)
    handles = [mpl.patches.Patch(color=CONTROL_COLORS[control], label=CONTROL_LABELS[control]) for control in CONTROLS]
    handles.append(mpl.lines.Line2D([], [], color="black", marker="$-$", linestyle="None", label="negative S"))
    ax.legend(handles=handles, frameon=False, loc="upper center", bbox_to_anchor=(0.50, -0.24), ncol=2)
    panel(ax, "D")
    fig.text(0.5, 0.005, "Cell labels are central local sensitivities from explicit symmetric 1/1.10 and 1.10 perturbations; no interpolation. "
             "Black frames identify states with $D/P \\geq 0.9$ in the independent 108-point turnover calculation.",
             ha="center", fontsize=7)
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / "Figure_11_kinetic_control_switching_atlas.pdf", bbox_inches="tight")
    fig.savefig(output / "Figure_11_kinetic_control_switching_atlas.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    sensitivity.to_csv(output / "Figure_11_control_switching_source_data.csv", index=False)


def classify_regimes(turnover: pd.DataFrame) -> pd.DataFrame:
    out = turnover.copy()
    peak = float(out["nh3_productivity_cm-3s-1"].max())
    ratio, product = out["NH3_loss_to_formation"], out["nh3_productivity_cm-3s-1"]
    out["regime"] = np.select([ratio.ge(0.9), product.ge(0.1 * peak) & ratio.lt(0.5),
                                product.lt(0.1 * peak) & ratio.lt(0.1)],
                               ["ion-loss cancellation", "productive hydrogenation", "formation-starved"],
                               default="transition turnover")
    out["log10_productivity"] = np.log10(out["nh3_productivity_cm-3s-1"].clip(lower=1e-300))
    return out


def make_figure12(sensitivity: pd.DataFrame, turnover: pd.DataFrame, dense: pd.DataFrame, output: Path) -> None:
    regimes = classify_regimes(turnover)
    ens_full, n2s_full = sorted(regimes.en_td.unique()), sorted(regimes.n2_fraction.unique())
    codes = {name: number for number, name in enumerate(REGIMES)}
    regime_matrix = (regimes.assign(code=regimes.regime.map(codes)).pivot(index="en_td", columns="n2_fraction",
                     values="code").reindex(index=ens_full, columns=n2s_full).to_numpy(float))
    local = (sensitivity.assign(abs_s=lambda x: x.sensitivity_productivity.abs()).sort_values("abs_s")
             .groupby(["en_td", "n2_fraction"], as_index=False).tail(1))
    direct = dense[["en_td", "n2_fraction", "named_h2star_pct"]]
    local = local.merge(direct, on=["en_td", "n2_fraction"], how="left")

    fig = plt.figure(figsize=(10.3, 5.5), constrained_layout=True)
    grid = fig.add_gridspec(1, 3, width_ratios=(1.05, 1.0, 1.45))
    ax = fig.add_subplot(grid[0, 0])
    cmap = mpl.colors.ListedColormap([REGIME_COLORS[name] for name in REGIMES])
    ax.imshow(regime_matrix, origin="lower", aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5)
    setup_axes(ax, n2s_full, ens_full)
    ax.set_yticks(range(0, len(ens_full), 2), [str(value) for value in ens_full[::2]])
    ax.set_title(r"NH$_3$ turnover regime, 108 accepted points")
    for name in REGIMES:
        ax.text(1.02, 0.95 - 0.09 * REGIMES.index(name), f"{(regimes.regime == name).sum()}: {name}",
                transform=ax.transAxes, fontsize=5.8, va="top")
    handles = [mpl.patches.Patch(color=REGIME_COLORS[name], label=name) for name in REGIMES]
    ax.legend(handles=handles, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=1)
    panel(ax, "A")

    ax = fig.add_subplot(grid[0, 1])
    ens, n2s = sorted(local.en_td.unique()), sorted(local.n2_fraction.unique())
    control_codes = {control: number for number, control in enumerate(CONTROLS)}
    matrix = (local.assign(code=local.control.map(control_codes)).pivot(index="en_td", columns="n2_fraction",
              values="code").reindex(index=ens, columns=n2s).to_numpy(float))
    cmap_control = mpl.colors.ListedColormap([CONTROL_COLORS[control] for control in CONTROLS])
    ax.imshow(matrix, origin="lower", aspect="auto", cmap=cmap_control, vmin=-0.5, vmax=2.5)
    negative = local.loc[local.sensitivity_productivity.lt(0)]
    for _, row in negative.iterrows():
        ax.text(n2s.index(row.n2_fraction), ens.index(row.en_td), "−", color="white", ha="center", va="center",
                fontsize=12, fontweight="bold")
    setup_axes(ax, n2s, ens); ax.set_title("Dominant sampled local control")
    ax.text(0.5, -0.25, "Targeted 30-point atlas\nminus sign = suppressive control", transform=ax.transAxes,
            ha="center", va="top", fontsize=7)
    panel(ax, "B")

    ax = fig.add_subplot(grid[0, 2]); ax.axis("off")
    ax.set_title("Two-layer kinetic design rule", loc="left", pad=10)
    boxes = [(0.04, 0.72, 0.28, 0.15, "Electronic H$_2$*\nNH-entry gate", "#DDEBF7"),
             (0.37, 0.72, 0.24, 0.15, "NH radical\nnetwork", "#E2F0D9"),
             (0.66, 0.72, 0.28, 0.15, "NH$_x$ hydrogenation\nproduct formation", "#D9EAD3"),
             (0.66, 0.36, 0.28, 0.15, "Ion-mediated\nNH$_3$ loss", "#F4CCCC")]
    for x, y, w, h, label, face in boxes:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02", facecolor=face, edgecolor="#555555", lw=0.8))
        ax.text(x+w/2, y+h/2, label, ha="center", va="center", fontsize=7.5)
    for start, end, color in [((0.32, 0.795), (0.37, 0.795), "#0072B2"), ((0.61, 0.795), (0.66, 0.795), "#009E73"),
                              ((0.80, 0.72), (0.80, 0.51), "#D55E00")]:
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=11, color=color, lw=1.4))
    share_min, share_max = local.named_h2star_pct.min(), local.named_h2star_pct.max()
    counts = local.control.value_counts()
    ax.text(0.04, 0.59, f"Direct NH entry remains named-H$_2$* dominated\nacross sampled states: {share_min:.1f}–{share_max:.1f}%.", fontsize=7.6)
    ax.text(0.04, 0.23, "Design rule:\noperate where NH$_x$ hydrogenation controls productivity\nand avoid the $D/P \\geq 0.9$ ion-loss cancellation regime.",
            fontsize=8.0, weight="bold")
    ax.text(0.04, 0.08, "Sampled dominant controls: " + ", ".join(f"{CONTROL_LABELS[key]} ({counts.get(key, 0)})" for key in CONTROLS),
            fontsize=6.4, wrap=True)
    ax.text(0.04, 0.01,
            "Regime definitions: formation-starved = $P<0.1P_{max}$ and $D/P<0.1$; productive = $P\\geq0.1P_{max}$ and $D/P<0.5$; "
            "loss cancellation = $D/P\\geq0.9$; remainder = transition.", fontsize=5.6, wrap=True)
    panel(ax, "C")
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / "Figure_12_two_layer_mechanistic_regime_map.pdf", bbox_inches="tight")
    fig.savefig(output / "Figure_12_two_layer_mechanistic_regime_map.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    regimes.to_csv(output / "Figure_12_mechanistic_regime_source_data.csv", index=False)
    local.to_csv(output / "Figure_12_local_control_overlay_source_data.csv", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot Hong Figures 11 and 12")
    parser.add_argument("--p7", type=Path, default=P7)
    parser.add_argument("--p5", type=Path, default=P5)
    parser.add_argument("--dense", type=Path, default=DENSE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    sensitivity, turnover, dense = pd.read_csv(args.p7), pd.read_csv(args.p5), pd.read_csv(args.dense)
    if len(sensitivity) != 90 or not np.isfinite(sensitivity.sensitivity_productivity).all():
        raise ValueError("Figures 11/12 require 90 finite P7 sensitivity rows")
    if len(turnover) != 108 or len(dense) != 108:
        raise ValueError("Figure 12 requires the complete accepted 108-point turnover and source tables")
    configure(); make_figure11(sensitivity, turnover, args.output); make_figure12(sensitivity, turnover, dense, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
