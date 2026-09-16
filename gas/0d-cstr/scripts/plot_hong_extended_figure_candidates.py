#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build standalone candidate Figures 13--16 from accepted Hong CSTR outputs.

These figures deliberately reuse only the accepted 108-point terminal windows
and the explicit 30-point P7 local-control atlas.  They are candidate figures,
not manuscript inserts: no kinetic calculation, interpolation, or manuscript
source file is modified by this script.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd

from _paths import LEGACY_HONG_ROOT

TOOL_DIR = Path(__file__).resolve().parent
PROJECT = LEGACY_HONG_ROOT
FIG = PROJECT / "analysis" / "figures_20260829"
DENSE = FIG / "Figure_3_dense_EN_N2_source_data.csv"
TURNOVER = FIG / "Figure_8_9_NHx_CSTR_summary.csv"
REACTIONS = FIG / "Figure_8_NHx_reaction_turnover.csv"
P7 = PROJECT / "analysis" / "p7_control_switching_atlas_20260830" / "P7_control_sensitivities.csv"

COLORS = {
    "named": "#0072B2", "rydberg": "#CC79A7", "n2d": "#56B4E9",
    "n2p": "#E69F00", "vib": "#009E73", "association": "#D55E00",
    "formation": "#009E73", "loss": "#D55E00", "net": "#0072B2",
    "starved": "#9E9E9E", "productive": "#009E73", "transition": "#E69F00",
    "cancellation": "#D55E00",
}
SOURCE_COLUMNS = {
    "named_h2star_pct": ("Named electronic H$_2$*", COLORS["named"]),
    "rydberg_h2star_pct": ("Rydberg-H$_2$*", COLORS["rydberg"]),
    "n2d_h2_pct": ("N($^2$D)+H$_2$", COLORS["n2d"]),
    "n2p_h2_pct": ("N($^2$P)+H$_2$", COLORS["n2p"]),
    "vib_h2_pct": ("vibrational H$_2$", COLORS["vib"]),
    "association_pct": ("H+N+M", COLORS["association"]),
}
CONTROL_COLORS = {
    "NH_to_NH3_association": COLORS["named"],
    "NH2_hydrogenation": COLORS["formation"],
    "proton_NH3_ionization": COLORS["loss"],
}
CONTROL_SHORT = {
    "NH_to_NH3_association": r"NH+H$_2$+M $\rightarrow$ NH$_3$+M",
    "NH2_hydrogenation": r"H+NH$_2$+M $\rightarrow$ NH$_3$+M",
    "proton_NH3_ionization": r"H$^+$+NH$_3$ $\rightarrow$ NH$_3^+$+H",
}


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 9,
        "axes.titlesize": 9.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 6.5, "axes.linewidth": 0.7, "pdf.fonttype": 42,
        "ps.fonttype": 42, "savefig.transparent": False,
    })


def panel(ax: plt.Axes, label: str) -> None:
    ax.text(-0.14, 1.08, label, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top", clip_on=False)


def save(fig: plt.Figure, output: Path, stem: str) -> None:
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(output / f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def decorate_grid(ax: plt.Axes, n2s: list[float], ens: list[int]) -> None:
    ax.set_xticks(range(len(n2s)), [f"{x:.1f}" for x in n2s])
    ax.set_yticks(range(0, len(ens), 2), [str(ens[i]) for i in range(0, len(ens), 2)])
    ax.set_xticks(np.arange(-0.5, len(n2s), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(ens), 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.35, alpha=0.75)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_xlabel(r"N$_2$ inlet fraction")
    ax.set_ylabel("E/N (Td)")


def add_margin(dense: pd.DataFrame) -> pd.DataFrame:
    out = dense.copy()
    secondary_cols = [key for key in SOURCE_COLUMNS if key != "named_h2star_pct"]
    out["largest_secondary_pct"] = out[secondary_cols].max(axis=1)
    out["largest_secondary_column"] = out[secondary_cols].idxmax(axis=1)
    out["named_to_secondary_ratio"] = out["named_h2star_pct"] / out["largest_secondary_pct"].clip(lower=1e-300)
    out["log10_named_to_secondary"] = np.log10(out["named_to_secondary_ratio"])
    out["dominance_margin_pct_point"] = out["named_h2star_pct"] - out["largest_secondary_pct"]
    return out


def matrix(table: pd.DataFrame, column: str) -> tuple[np.ndarray, list[int], list[float]]:
    ens = sorted(table.en_td.unique().tolist())
    n2s = sorted(table.n2_fraction.unique().tolist())
    data = table.pivot(index="en_td", columns="n2_fraction", values=column).reindex(index=ens, columns=n2s)
    return data.to_numpy(float), ens, n2s


def figure13(dense: pd.DataFrame, turnover: pd.DataFrame, output: Path) -> None:
    table = add_margin(dense)
    values, ens, n2s = matrix(table, "log10_named_to_secondary")
    fig, axes = plt.subplots(2, 2, figsize=(10.7, 7.1), constrained_layout=True)
    ax = axes[0, 0]
    im = ax.imshow(values, origin="lower", aspect="auto", cmap="cividis")
    decorate_grid(ax, n2s, ens)
    ax.set_title("Named-H$_2$* dominance over the largest secondary source")
    cb = fig.colorbar(im, ax=ax, fraction=0.047, pad=0.04)
    cb.set_label(r"log$_{10}(\Phi_{named}/\Phi_{largest\ secondary})$")
    panel(ax, "A")

    ax = axes[0, 1]
    categories = list(SOURCE_COLUMNS)[1:]
    names = [SOURCE_COLUMNS[key][0] for key in categories]
    cmap = mpl.colors.ListedColormap([SOURCE_COLUMNS[key][1] for key in categories])
    codes = pd.Categorical(table.largest_secondary_column, categories=categories).codes
    payload = (table.assign(code=codes).pivot(index="en_td", columns="n2_fraction", values="code")
               .reindex(index=ens, columns=n2s).to_numpy(float))
    ax.imshow(payload, origin="lower", aspect="auto", cmap=cmap, vmin=-0.5, vmax=len(categories) - 0.5)
    decorate_grid(ax, n2s, ens)
    ax.set_title("Identity of the largest secondary direct-NH source")
    handles = [mpl.patches.Patch(color=SOURCE_COLUMNS[key][1], label=names[i]) for i, key in enumerate(categories)]
    ax.legend(handles=handles, frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.22))
    panel(ax, "B")

    ax = axes[1, 0]
    ax.hist(table.log10_named_to_secondary, bins=np.linspace(table.log10_named_to_secondary.min(),
            table.log10_named_to_secondary.max(), 15), color=COLORS["named"], edgecolor="white")
    ax.axvline(table.log10_named_to_secondary.median(), color="black", ls="--", lw=1.0,
               label=f"median = {table.log10_named_to_secondary.median():.2f}")
    ax.set_xlabel(r"log$_{10}(\Phi_{named}/\Phi_{largest\ secondary})$")
    ax.set_ylabel("Accepted grid points")
    ax.set_title("Distribution of the dominance margin (n = 108)")
    ax.legend(frameon=False)
    panel(ax, "C")

    ax = axes[1, 1]
    merged = table.merge(turnover[["en_td", "n2_fraction", "nh3_productivity_cm-3s-1"]],
                         on=["en_td", "n2_fraction"], validate="one_to_one")
    points = ax.scatter(merged.log10_named_to_secondary, merged["nh3_productivity_cm-3s-1"],
                        c=merged.en_td, s=25 + 45 * merged.n2_fraction, cmap="cividis",
                        edgecolor="black", linewidth=0.25)
    ax.set_yscale("log")
    ax.set_xlabel(r"log$_{10}(\Phi_{named}/\Phi_{largest\ secondary})$")
    ax.set_ylabel(r"NH$_3$ outlet productivity (cm$^{-3}$ s$^{-1}$)")
    ax.set_title("Route dominance and product level are distinct coordinates")
    cb = fig.colorbar(points, ax=ax, fraction=0.047, pad=0.04); cb.set_label("E/N (Td)")
    panel(ax, "D")
    fig.text(0.5, -0.035, "All panels use the accepted 108-point CSTR/CW grid. The margin compares the named family with the single largest secondary direct-NH family at each state.", ha="center", fontsize=7)
    save(fig, output, "Figure_13_direct_NH_dominance_margin")
    table.to_csv(output / "Figure_13_direct_NH_dominance_margin_source_data.csv", index=False)


def representative_rows(table: pd.DataFrame) -> pd.DataFrame:
    keys = [(20, 0.1, "Low field"), (140, 0.1, r"NH$_3$ maximum"), (240, 0.9, "High field, N$_2$-rich")]
    rows = []
    for en, n2, label in keys:
        item = table.loc[(table.en_td.eq(en)) & np.isclose(table.n2_fraction, n2)].iloc[0].copy()
        item["label"] = label
        rows.append(item)
    return pd.DataFrame(rows)


def figure14(dense: pd.DataFrame, turnover: pd.DataFrame, reactions: pd.DataFrame, output: Path) -> None:
    reps = representative_rows(add_margin(dense).merge(turnover, on=["en_td", "n2_fraction", "tag"], validate="one_to_one"))
    fig = plt.figure(figsize=(11.5, 6.6), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=(1.2, 1.0), height_ratios=(1.0, 1.05))
    ax = fig.add_subplot(grid[:, 0])
    source_keys = list(SOURCE_COLUMNS)
    bottoms = np.zeros(len(reps))
    x = np.arange(len(reps))
    for key in source_keys:
        label, color = SOURCE_COLUMNS[key]
        ax.bar(x, reps[key], bottom=bottoms, color=color, width=0.62, label=label)
        bottoms += reps[key].to_numpy(float)
    ax.set_ylim(0, 100); ax.set_ylabel("Direct NH-source contribution (%)")
    ax.set_xticks(x, [f"{row.label}\n{int(row.en_td)} Td, $x_{{N_2}}$={row.n2_fraction:.1f}" for row in reps.itertuples()])
    ax.set_title("Upstream direct-NH entry remains electronically excited-H$_2$ dominated")
    ax.legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.13))
    panel(ax, "A")

    ax = fig.add_subplot(grid[0, 1])
    labels = [str(item) for item in reps.label]
    for name, column, color, marker in [(r"formation $P$", "NH3_formation_cm-3s-1", COLORS["formation"], "o"),
                                        (r"destruction $D$", "NH3_loss_cm-3s-1", COLORS["loss"], "s"),
                                        (r"net outlet", "nh3_productivity_cm-3s-1", COLORS["net"], "D")]:
        ax.plot(x, reps[column], marker=marker, ms=6, lw=1.5, color=color, label=name)
    ax.set_yscale("log"); ax.set_xticks(x, labels, rotation=10, ha="right")
    ax.set_ylabel(r"NH$_3$ flux (cm$^{-3}$ s$^{-1}$)")
    ax.set_title("Downstream product turnover")
    ax.legend(frameon=False, loc="lower left")
    panel(ax, "B")

    ax = fig.add_subplot(grid[1, 1]); ax.axis("off")
    ax.set_title("Dominant formation/loss steps at the product level", loc="left", pad=8)
    y = 0.91
    for _, row in reps.iterrows():
        subset = reactions.loc[reactions.tag.eq(row["tag"]) & reactions.target.eq("NH3")]
        source = subset.loc[subset.direction.eq("formation")].nlargest(1, "flux_cm-3s-1").iloc[0]
        sink = subset.loc[subset.direction.eq("loss")].nlargest(1, "flux_cm-3s-1").iloc[0]
        ax.text(0.01, y, f"{row['label']}: {int(row['en_td'])} Td, $x_{{N_2}}$={row['n2_fraction']:.1f}", fontsize=7.2, weight="bold")
        ax.text(0.04, y - 0.09, "formation  " + source.reaction_short, color=COLORS["formation"], fontsize=6.2)
        ax.text(0.04, y - 0.17, "loss          " + sink.reaction_short, color=COLORS["loss"], fontsize=6.2)
        outlet = row["nh3_productivity_cm-3s-1"]
        ax.text(0.04, y - 0.25, f"$D/P={row['NH3_loss_to_formation']:.4f}$; outlet = {outlet:.2e} cm$^{{-3}}$ s$^{{-1}}$", fontsize=6.4)
        y -= 0.32
    panel(ax, "C")
    fig.text(0.5, -0.035, "Candidate causal cascade: source-family fractions, absolute NH$_3$ turnover and the largest elementary product-level formation/loss term are reported on their native scales; widths are not cross-scale Sankey fluxes.", ha="center", fontsize=7)
    save(fig, output, "Figure_14_NH_entry_to_NH3_turnover_cascade")
    reps.to_csv(output / "Figure_14_representative_cascade_source_data.csv", index=False)


def classify(turnover: pd.DataFrame) -> pd.DataFrame:
    out = turnover.copy()
    peak = out["nh3_productivity_cm-3s-1"].max()
    ratio, prod = out.NH3_loss_to_formation, out["nh3_productivity_cm-3s-1"]
    out["regime"] = np.select([ratio.ge(0.9), prod.ge(0.1 * peak) & ratio.lt(0.5), prod.lt(0.1 * peak) & ratio.lt(0.1)],
                              ["ion-loss cancellation", "productive hydrogenation", "formation-starved"], default="transition turnover")
    return out


def figure15(turnover: pd.DataFrame, p7: pd.DataFrame, output: Path) -> None:
    regimes = classify(turnover)
    local = p7.merge(regimes[["en_td", "n2_fraction", "NH3_loss_to_formation", "nh3_productivity_cm-3s-1", "regime"]],
                     on=["en_td", "n2_fraction"], validate="many_to_one")
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.1), constrained_layout=True)
    ax = axes[0, 0]
    regime_color = {"formation-starved": COLORS["starved"], "productive hydrogenation": COLORS["productive"],
                    "transition turnover": COLORS["transition"], "ion-loss cancellation": COLORS["cancellation"]}
    for regime, subset in regimes.groupby("regime", sort=False):
        ax.scatter(subset.NH3_loss_to_formation, subset["nh3_productivity_cm-3s-1"], s=22 + 35 * subset.n2_fraction,
                   color=regime_color[regime], edgecolor="black", linewidth=0.25, label=regime, alpha=0.9)
    ax.axvline(0.9, color="black", ls="--", lw=0.9, label=r"$D/P=0.9$")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(r"NH$_3$ destruction / formation, $D/P$")
    ax.set_ylabel(r"NH$_3$ outlet productivity (cm$^{-3}$ s$^{-1}$)")
    ax.set_title("Continuous product-turnover competition")
    ax.legend(frameon=False, fontsize=6.0, loc="lower left")
    panel(ax, "A")

    ax = axes[0, 1]
    for control, subset in local.groupby("control", sort=False):
        ax.scatter(subset.NH3_loss_to_formation, subset.sensitivity_productivity, s=20,
                   color=CONTROL_COLORS[control], edgecolor="black", linewidth=0.2, alpha=0.85,
                   label=CONTROL_SHORT[control])
    ax.axvline(0.9, color="black", ls="--", lw=0.9)
    ax.axhline(0, color="#555555", lw=0.7)
    ax.set_xscale("log"); ax.set_xlabel(r"$D/P$ at the unperturbed state")
    ax.set_ylabel(r"Local productivity sensitivity, $S$")
    ax.set_title("Selected kinetic levers across the turnover coordinate")
    ax.legend(frameon=False, fontsize=5.7, loc="upper left")
    panel(ax, "B")

    ax = axes[1, 0]
    for n2, subset in local.loc[local.control.eq("proton_NH3_ionization")].groupby("n2_fraction"):
        ax.plot(subset.en_td, subset.sensitivity_productivity, marker="o", ms=3.3, lw=1.1, label=f"$x_{{N_2}}$={n2:.1f}")
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xticks(sorted(local.en_td.unique())); ax.set_xlabel("E/N (Td)")
    ax.set_ylabel(r"$S$ of proton-driven NH$_3$ ionization")
    ax.set_title("Ionization becomes a high-field suppression lever")
    ax.legend(frameon=False, ncol=2, fontsize=5.8, loc="lower left")
    panel(ax, "C")

    ax = axes[1, 1]
    for n2, subset in local.loc[local.control.eq("NH_to_NH3_association")].groupby("n2_fraction"):
        ax.plot(subset.en_td, subset.sensitivity_productivity, marker="o", ms=3.3, lw=1.1, label=f"$x_{{N_2}}$={n2:.1f}")
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xticks(sorted(local.en_td.unique())); ax.set_xlabel("E/N (Td)")
    ax.set_ylabel(r"$S$ of NH+H$_2$+M $\rightarrow$ NH$_3$+M")
    ax.set_title("Hydrogenation control peaks near the productive window")
    panel(ax, "D")
    fig.text(0.5, -0.035, "Panels B--D use the explicit 30-point P7 atlas (five N$_2$ fractions × six E/N values); panel A uses all 108 accepted turnover points. No interpolation is used.", ha="center", fontsize=7)
    save(fig, output, "Figure_15_product_turnover_control_phase_space")
    local.to_csv(output / "Figure_15_turnover_control_phase_space_source_data.csv", index=False)


def figure16(dense: pd.DataFrame, turnover: pd.DataFrame, p7: pd.DataFrame, output: Path) -> None:
    table = add_margin(dense).merge(classify(turnover)[["en_td", "n2_fraction", "nh3_productivity_cm-3s-1", "NH3_loss_to_formation", "regime"]],
                                   on=["en_td", "n2_fraction"], validate="one_to_one")
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.0), constrained_layout=True)
    ax = axes[0, 0]
    residual = np.maximum(table.rel_state_max, table.rel_dnh3_cycle)
    pts = ax.scatter(residual, table.log10_named_to_secondary, c=np.log10(table["nh3_productivity_cm-3s-1"]),
                     cmap="cividis", s=25 + 42 * table.n2_fraction, edgecolor="black", linewidth=0.2)
    ax.axvline(1e-3, color="black", ls="--", lw=0.9, label="P0 threshold")
    ax.set_xscale("log"); ax.set_xlabel("max(P0 state, NH$_3$ residual)")
    ax.set_ylabel(r"log$_{10}(\Phi_{named}/\Phi_{largest\ secondary})$")
    ax.set_title("Numerical acceptance and mechanistic margin")
    ax.legend(frameon=False, loc="lower left")
    cb = fig.colorbar(pts, ax=ax, fraction=0.047, pad=0.04); cb.set_label(r"log$_{10}$ productivity")
    panel(ax, "A")

    ax = axes[0, 1]
    values, ens, n2s = matrix(table, "dominance_margin_pct_point")
    im = ax.imshow(values, origin="lower", aspect="auto", cmap="viridis")
    decorate_grid(ax, n2s, ens)
    ax.set_title("Absolute direct-NH dominance margin")
    cb = fig.colorbar(im, ax=ax, fraction=0.047, pad=0.04); cb.set_label("Named minus largest secondary source (percentage points)")
    panel(ax, "B")

    ax = axes[1, 0]
    local = p7.merge(classify(turnover)[["en_td", "n2_fraction", "regime"]], on=["en_td", "n2_fraction"], validate="many_to_one")
    order = ["formation-starved", "productive hydrogenation", "transition turnover", "ion-loss cancellation"]
    offsets = {"NH_to_NH3_association": -0.16, "NH2_hydrogenation": 0.0, "proton_NH3_ionization": 0.16}
    for i, regime in enumerate(order):
        for control, subset in local.loc[local.regime.eq(regime)].groupby("control"):
            jitter = np.linspace(-0.045, 0.045, len(subset)) if len(subset) > 1 else np.array([0.0])
            ax.scatter(np.full(len(subset), i + offsets[control]) + jitter, subset.sensitivity_productivity,
                       color=CONTROL_COLORS[control], s=20, edgecolor="black", linewidth=0.2)
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xticks(range(len(order)), ["formation\nstarved", "productive\nhydrogenation", "transition\nturnover", "ion-loss\ncancellation"])
    ax.set_ylabel(r"Local productivity sensitivity, $S$")
    ax.set_title("Control distribution by turnover regime (P7 points)")
    handles = [mpl.lines.Line2D([], [], marker="o", ls="", color=CONTROL_COLORS[key], label=CONTROL_SHORT[key]) for key in CONTROL_COLORS]
    ax.legend(handles=handles, frameon=False, fontsize=5.6, loc="upper left")
    panel(ax, "C")

    ax = axes[1, 1]; ax.axis("off")
    ax.set_title("Candidate-evidence summary", loc="left", pad=8)
    lines = [
        ("Numerical basis", "108/108 accepted P0 terminal windows; two documented retries."),
        ("Direct NH entry", f"Named H$_2$* exceeds the largest secondary family by {table.dominance_margin_pct_point.min():.1f}–{table.dominance_margin_pct_point.max():.1f} percentage points."),
        ("Product turnover", "Intermediate field maximizes product; high E/N can approach D/P = 1."),
        ("Kinetic leverage", "Explicit P7 perturbations distinguish NHx hydrogenation from ion-loss control."),
    ]
    y = 0.88
    for heading, body in lines:
        ax.add_patch(FancyBboxPatch((0.02, y - 0.11), 0.95, 0.13, boxstyle="round,pad=0.015", facecolor="#F7F7F7", edgecolor="#AAAAAA", lw=0.6))
        ax.text(0.05, y - 0.015, heading + ":", fontsize=7.2, weight="bold", va="top")
        ax.text(0.05, y - 0.065, body, fontsize=6.4, va="top", wrap=True)
        y -= 0.20
    ax.text(0.02, 0.035, "This overview is a candidate closing figure. It summarizes model evidence and does not add experimental validation or a self-consistent electron solution.", fontsize=6.4, wrap=True)
    panel(ax, "D")
    fig.text(0.5, -0.035, "Candidate Figure 16 combines accepted-grid diagnostics with the explicit P7 atlas. It is a transparency/robustness overview, not an independent data set.", ha="center", fontsize=7)
    save(fig, output, "Figure_16_numerical_mechanistic_robustness_overview")
    table.to_csv(output / "Figure_16_robustness_overview_source_data.csv", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate standalone Hong candidate Figures 13--16")
    parser.add_argument("--output", type=Path, default=FIG)
    args = parser.parse_args()
    dense, turnover = pd.read_csv(DENSE), pd.read_csv(TURNOVER)
    reactions, p7 = pd.read_csv(REACTIONS), pd.read_csv(P7)
    if len(dense) != 108 or len(turnover) != 108 or len(p7) != 90:
        raise ValueError("Expected 108 dense points, 108 turnover points, and 90 P7 sensitivity rows")
    if not dense.run_status.isin(["success", "retry_success"]).all():
        raise ValueError("Candidate figures refuse a non-accepted dense-scan point")
    configure()
    figure13(dense, turnover, args.output)
    figure14(dense, turnover, reactions, args.output)
    figure15(turnover, p7, args.output)
    figure16(dense, turnover, p7, args.output)
    print("Candidate Figures 13--16 written to", args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
