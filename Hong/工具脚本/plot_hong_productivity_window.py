#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot the NH3 field window from the validated 108-point CSTR data only."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


TOOL_DIR = Path(__file__).resolve().parent
PROJECT = TOOL_DIR.parent
DEFAULT_DATA = PROJECT / "analysis" / "figures_20260829" / "Figure_3_dense_EN_N2_source_data.csv"
DEFAULT_OUTPUT = PROJECT / "analysis" / "figures_20260829"


def configure_style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 9,
        "axes.titlesize": 9, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 6.8, "axes.linewidth": 0.7, "pdf.fonttype": 42,
        "ps.fonttype": 42, "savefig.transparent": False,
    })


def panel_label(ax, label: str) -> None:
    ax.text(-0.15, 1.08, label, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top")


def make_figure(data: pd.DataFrame, output: Path) -> pd.DataFrame:
    n2s = sorted(data["n2_fraction"].unique())
    ens = sorted(data["en_td"].unique())
    cmap = mpl.colormaps["viridis"]
    colors = {n2: cmap(i / (len(n2s) - 1)) for i, n2 in enumerate(n2s)}
    peaks = (data.loc[data.groupby("n2_fraction")["nh3_final_cm-3"].idxmax()]
             .sort_values("n2_fraction").copy())
    peaks["nh3_at_240_cm-3"] = [float(data.loc[(data["n2_fraction"] == n2) &
                                                (data["en_td"] == max(ens)),
                                                "nh3_final_cm-3"].iloc[0]) for n2 in peaks["n2_fraction"]]
    peaks["high_field_to_peak"] = peaks["nh3_at_240_cm-3"] / peaks["nh3_final_cm-3"]

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.55), constrained_layout=True)
    ax = axes[0]
    for n2 in n2s:
        sub = data.loc[data["n2_fraction"].eq(n2)].sort_values("en_td")
        ax.plot(sub["en_td"], sub["nh3_final_cm-3"], color=colors[n2], marker="o",
                ms=2.7, lw=1.25, label=f"{n2:.1f}")
    ax.set_yscale("log")
    ax.set_xticks(ens[::2])
    ax.set_xlabel("E/N (Td)")
    ax.set_ylabel(r"Steady [NH$_3$] (cm$^{-3}$)")
    ax.set_title("NH$_3$ response to reduced field")
    ax.legend(title="N$_2$ fraction", ncol=3, loc="lower center", frameon=False,
              handlelength=1.8, columnspacing=1.0)
    panel_label(ax, "A")
    optimum = data.loc[data["nh3_final_cm-3"].idxmax()]
    ax.plot(optimum["en_td"], optimum["nh3_final_cm-3"], marker="*", ms=10,
            mec="black", mfc="white", mew=0.8, zorder=5)
    ax.annotate(f"global maximum\n{int(optimum['en_td'])} Td, $x_{{N_2}}$={optimum['n2_fraction']:.1f}",
                (optimum["en_td"], optimum["nh3_final_cm-3"]), xytext=(8, -27), textcoords="offset points",
                arrowprops={"arrowstyle": "-", "lw": 0.6}, fontsize=6.2)

    ax = axes[1]
    points = ax.scatter(peaks["n2_fraction"], peaks["nh3_final_cm-3"],
                        c=peaks["en_td"], cmap="cividis", s=34, edgecolor="black", linewidth=0.35,
                        vmin=min(ens), vmax=max(ens), zorder=3)
    ax.plot(peaks["n2_fraction"], peaks["nh3_final_cm-3"], color="#777777", lw=0.8, zorder=1)
    for _, row in peaks.iterrows():
        ax.annotate(f"{int(row['en_td'])}", (row["n2_fraction"], row["nh3_final_cm-3"]),
                    xytext=(0, 5), textcoords="offset points", ha="center", fontsize=6)
    ax.set_yscale("log")
    ax.set_xlabel("N$_2$ inlet fraction")
    ax.set_ylabel(r"Maximum steady [NH$_3$] (cm$^{-3}$)")
    ax.set_title("Composition-resolved field of the NH$_3$ maximum")
    ax.set_xticks(n2s)
    cbar = fig.colorbar(points, ax=ax, fraction=0.047, pad=0.04)
    cbar.set_label("E/N at maximum (Td)")
    panel_label(ax, "B")

    ax = axes[2]
    points = ax.scatter(data["named_h2star_pct"], data["nh3_final_cm-3"],
                        c=data["en_td"], cmap="plasma", s=20 + 35 * data["n2_fraction"],
                        edgecolor="black", linewidth=0.25, alpha=0.9, vmin=min(ens), vmax=max(ens))
    ax.set_yscale("log")
    ax.set_xlabel("Named electronic H$_2$* NH-source share (%)")
    ax.set_ylabel(r"Steady [NH$_3$] (cm$^{-3}$)")
    ax.set_title("Persistent NH entry does not determine NH$_3$ inventory")
    ax.text(0.03, 0.05, "Colour: E/N\nMarker area: N$_2$ fraction", transform=ax.transAxes,
            va="bottom", fontsize=6.5)
    cbar = fig.colorbar(points, ax=ax, fraction=0.047, pad=0.04)
    cbar.set_label("E/N (Td)")
    panel_label(ax, "C")

    fig.text(0.5, 0.005,
             "CSTR/CW, T$_g$ = 300 K, fixed $n_e$ = 1.17×10$^8$ cm$^{-3}$, τ = 10 ms; "
             "108 explicit deterministic kinetic-model calculations (no replicate error bars). ★ = global inventory maximum.",
             ha="center", va="top", fontsize=7)
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / "Figure_4_productivity_window.pdf", bbox_inches="tight")
    fig.savefig(output / "Figure_4_productivity_window.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    return peaks


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Plot the validated Hong CSTR NH3 field window")
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    configure_style()
    data = pd.read_csv(args.data)
    required = {"en_td", "n2_fraction", "nh3_final_cm-3", "named_h2star_pct"}
    missing = required - set(data.columns)
    if missing or len(data) != 108:
        raise RuntimeError(f"Expected the validated 108-point source table; missing={sorted(missing)}, rows={len(data)}")
    if set(data["run_status"]) - {"success", "retry_success"}:
        raise RuntimeError("Source table contains a point without a valid steady-state status")
    args.output.mkdir(parents=True, exist_ok=True)
    data.to_csv(args.output / "Figure_4_productivity_window_source_data.csv", index=False, float_format="%.12g")
    peaks = make_figure(data, args.output)
    peaks.to_csv(args.output / "Figure_4_peak_conditions.csv", index=False, float_format="%.12g")
    print("Figure 4 written to", args.output.resolve())
    print(peaks[["n2_fraction", "en_td", "nh3_final_cm-3", "high_field_to_peak"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
