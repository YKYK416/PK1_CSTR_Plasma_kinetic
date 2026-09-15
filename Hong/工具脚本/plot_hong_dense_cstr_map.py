#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create dense, validated CSTR NH-production/pathway maps for the Hong mechanism.

Only status=success or retry_success cases are admitted.  Source fractions are
integrated from reaction-of-progress data over the final three period-steady
cycles, never over the transient start-up interval.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


TOOL_DIR = Path(__file__).resolve().parent
PROJECT = TOOL_DIR.parent
DEFAULT_ROOT = (PROJECT / "Reproduction" / "2017Hong" / "方案2_脉冲能效地图"
                / "build_pulse" / "runs" / "p4_dense_cstr_20260829")
sys.path.insert(0, str(TOOL_DIR))

from hong_nh_diagnostics import integrate_nh_channels, steady_window_start  # noqa: E402


def configure_style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 9,
        "axes.titlesize": 9, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "axes.linewidth": 0.7, "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.transparent": False,
    })


def sources_for_case(root: Path, tag: str) -> tuple[dict[str, float], int, int]:
    outdir = root / tag
    first, last = steady_window_start(outdir / f"pulse_cycles_{tag}.csv", keep_cycles=3)
    rows, _start, _end = integrate_nh_channels(outdir / f"pulse_rates_{tag}.csv", first_cycle=first)
    return {str(row["channel"]): float(row["share_pct"]) for row in rows}, first, last


def build_table(root: Path) -> pd.DataFrame:
    summary = pd.read_csv(root / "scan_summary.csv")
    allowed = {"success", "retry_success"}
    bad = summary.loc[~summary["status"].isin(allowed)]
    if not bad.empty:
        raise RuntimeError("Refusing to plot non-steady scan points:\n" +
                           bad[["en", "n2frac", "status"]].to_string(index=False))
    records: list[dict] = []
    for index, row in summary.sort_values(["en", "n2frac"]).iterrows():
        tag = str(row["tag"])
        shares, cycle_first, cycle_last = sources_for_case(root, tag)
        records.append({
            "en_td": int(row["en"]), "n2_fraction": float(row["n2frac"]), "tag": tag,
            "run_status": row["status"], "attempts": int(row["attempts"]),
            "cycles_used": int(row["cycles_used"]), "cycle_first": cycle_first,
            "cycle_last": cycle_last, "nh3_final_cm-3": float(row["NH3_final_cm-3"]),
            "nh3_per_j": float(row["NH3_per_J"]), "rel_state_max": float(row["rel_state_max"]),
            "rel_dnh3_cycle": float(row["rel_dNH3_cycle"]),
            "named_h2star_pct": shares["Named-H2*"], "rydberg_h2star_pct": shares["Rydberg-H2*"],
            "n2d_h2_pct": shares["N(2D)+H2"], "n2p_h2_pct": shares["N(2P)+H2"],
            "vib_h2_pct": shares["Vib-H2"], "association_pct": shares["Association H+N+M"],
        })
        print(f"[{index + 1:03d}/{len(summary):03d}] {tag}", flush=True)
    return pd.DataFrame(records)


def heatmap(ax, data: np.ndarray, ens: list[int], n2s: list[float], *, title: str,
            cbar_label: str, cmap: str, vmin: float | None = None, vmax: float | None = None,
            log10: bool = False):
    shown = np.log10(np.maximum(data, 1e-300)) if log10 else data
    im = ax.imshow(shown, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax, origin="upper")
    ax.set_xticks(range(len(n2s)), [f"{x:.1f}" for x in n2s])
    ax.set_yticks(range(len(ens)), [str(x) for x in ens])
    ax.set_xlabel("N$_2$ inlet fraction")
    ax.set_ylabel("E/N (Td)")
    ax.set_title(title)
    ax.set_xticks(np.arange(-0.5, len(n2s), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(ens), 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.32, alpha=0.72)
    ax.tick_params(which="minor", bottom=False, left=False)
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(cbar_label)
    return im


def label_panel(ax, label: str) -> None:
    ax.text(-0.15, 1.07, label, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top")


def mark_retries(axes, table: pd.DataFrame, ens: list[int], n2s: list[float]) -> None:
    retry = table.loc[table["run_status"].eq("retry_success")]
    for _, row in retry.iterrows():
        x, y = n2s.index(float(row["n2_fraction"])), ens.index(int(row["en_td"]))
        for ax in axes.flat:
            ax.plot(x, y, marker="x", ms=7, mew=1.25, color="white")
            ax.plot(x, y, marker="x", ms=5, mew=0.75, color="black")


def make_figure(table: pd.DataFrame, output: Path) -> None:
    ens = sorted(table["en_td"].unique().tolist())
    n2s = sorted(table["n2_fraction"].unique().tolist())

    def grid(column: str) -> np.ndarray:
        return (table.pivot(index="en_td", columns="n2_fraction", values=column)
                .reindex(index=ens, columns=n2s).to_numpy(float))

    nh3 = grid("nh3_final_cm-3")
    named = grid("named_h2star_pct")
    rydberg = grid("rydberg_h2star_pct")
    n2p = grid("n2p_h2_pct")

    fig, axes = plt.subplots(2, 2, figsize=(8.8, 6.5), constrained_layout=True)
    heatmap(axes[0, 0], nh3, ens, n2s, title="Steady NH$_3$ inventory",
            cbar_label=r"log$_{10}$[NH$_3$] (cm$^{-3}$)", cmap="cividis",
            vmin=np.floor(np.log10(max(nh3[nh3 > 0].min(), 1e-1))),
            vmax=np.ceil(np.log10(nh3.max())), log10=True)
    label_panel(axes[0, 0], "A")
    best = np.unravel_index(np.argmax(nh3), nh3.shape)
    axes[0, 0].plot(best[1], best[0], marker="*", ms=10, mec="black", mfc="white", mew=0.8)
    axes[0, 0].text(0.02, 0.03, "★ global maximum", transform=axes[0, 0].transAxes,
                    fontsize=6.5, bbox={"facecolor": "white", "edgecolor": "#AAAAAA", "alpha": 0.9, "pad": 1.2})

    heatmap(axes[0, 1], named, ens, n2s, title="Named electronic H$_2$* → NH",
            cbar_label="Integrated NH-source share (%)", cmap="viridis",
            vmin=min(90.0, float(named.min())), vmax=100.0)
    label_panel(axes[0, 1], "B")
    axes[0, 1].text(0.02, 0.03, f"range {named.min():.1f}–{named.max():.1f}%", transform=axes[0, 1].transAxes,
                    fontsize=6.5, bbox={"facecolor": "white", "edgecolor": "#AAAAAA", "alpha": 0.9, "pad": 1.2})

    heatmap(axes[1, 0], rydberg, ens, n2s, title="Rydberg-H$_2$* → NH",
            cbar_label="Integrated NH-source share (%)", cmap="magma",
            vmin=0.0, vmax=max(1.0, float(np.ceil(rydberg.max()))))
    label_panel(axes[1, 0], "C")
    axes[1, 0].text(0.02, 0.03, f"maximum {rydberg.max():.2f}%", transform=axes[1, 0].transAxes,
                    fontsize=6.5, bbox={"facecolor": "white", "edgecolor": "#AAAAAA", "alpha": 0.9, "pad": 1.2})

    heatmap(axes[1, 1], n2p, ens, n2s, title="N(2P) + H$_2$ → NH",
            cbar_label="Integrated NH-source share (%)", cmap="plasma",
            vmin=0.0, vmax=max(1.0, float(np.ceil(n2p.max()))))
    label_panel(axes[1, 1], "D")
    axes[1, 1].text(0.02, 0.03, f"maximum {n2p.max():.2f}%", transform=axes[1, 1].transAxes,
                    fontsize=6.5, bbox={"facecolor": "white", "edgecolor": "#AAAAAA", "alpha": 0.9, "pad": 1.2})

    mark_retries(axes, table, ens, n2s)
    fig.text(0.5, 0.005,
             "CSTR/CW, T$_g$ = 300 K, $n_e$ = 1.17×10$^8$ cm$^{-3}$, τ = 10 ms; "
             "each cell is one accepted calculation; ROP integrated over final 3 P0-steady cycles. × = successful numerical retry.",
             ha="center", va="top", fontsize=7)
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / "Figure_3_dense_EN_N2_CSTR_maps.pdf", bbox_inches="tight")
    fig.savefig(output / "Figure_3_dense_EN_N2_CSTR_maps.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Plot dense Hong CSTR NH mechanism maps")
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--output", type=Path, default=PROJECT / "analysis" / "figures_20260829")
    args = parser.parse_args(argv)
    configure_style()
    root = args.root.resolve()
    table = build_table(root)
    args.output.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output / "Figure_3_dense_EN_N2_source_data.csv", index=False, float_format="%.12g")
    make_figure(table, args.output)
    print("\nValidated points:", len(table))
    print("NH3 range:", f"{table['nh3_final_cm-3'].min():.6g}", "to", f"{table['nh3_final_cm-3'].max():.6g}", "cm^-3")
    print("Named-H2* range:", f"{table['named_h2star_pct'].min():.6g}", "to", f"{table['named_h2star_pct'].max():.6g}", "%")
    print("Output:", args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
