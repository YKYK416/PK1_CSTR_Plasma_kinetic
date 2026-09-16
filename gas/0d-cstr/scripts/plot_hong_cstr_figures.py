#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build publication-ready CSTR mechanism figures from validated Hong outputs.

Only the CSTR/CW P2-P3 data sets listed below are read.  Long-afterglow pulse
attempts are intentionally excluded because they did not pass the numerical
validity gate.  NH source fractions are right-rule ROP integrals over the
final three period-steady cycles.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from _paths import LEGACY_HONG_ROOT

TOOL_DIR = Path(__file__).resolve().parent
PROJECT = LEGACY_HONG_ROOT
BUILD_RUNS = PROJECT / "Reproduction" / "2017Hong" / "方案2_脉冲能效地图" / "build_pulse" / "runs"
P3_CAUSAL = PROJECT / "analysis" / "p3_cstr_20260828_v2"
P3_RYDBERG = PROJECT / "analysis" / "p3_rydberg_loss_cstr_20260828"

sys.path.insert(0, str(TOOL_DIR))
from hong_nh_diagnostics import integrate_nh_channels, steady_window_start  # noqa: E402
import zdp_gui as gui  # noqa: E402


OKABE_ITO = {
    "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73",
    "yellow": "#F0E442", "blue": "#0072B2", "vermillion": "#D55E00",
    "purple": "#CC79A7", "black": "#000000", "gray": "#777777",
}


def configure_style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 9,
        "axes.titlesize": 9.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 7, "axes.linewidth": 0.7, "pdf.fonttype": 42,
        "ps.fonttype": 42, "savefig.transparent": False,
    })


def source_fractions(outdir: Path, tag: str) -> dict[str, float]:
    """Calculate primary NH source fractions only in the terminal steady window."""
    first, _last = steady_window_start(outdir / f"pulse_cycles_{tag}.csv")
    rows, _start, _end = integrate_nh_channels(outdir / f"pulse_rates_{tag}.csv", first)
    return {row["channel"]: float(row["share_pct"]) for row in rows}


def metrics(outdir: Path, tag: str) -> dict:
    return gui.extract_metrics("pulse", outdir, tag)


def p2_point(en: int, n2: float) -> dict:
    code = {0.1: "01", 0.3333: "3333", 0.8: "80"}[n2]
    tag = f"p2_cstr_e{en}_n2{code}"
    outdir = BUILD_RUNS / tag
    return {"en": en, "n2": n2, "tag": tag, "metrics": metrics(outdir, tag),
            "shares": source_fractions(outdir, tag)}


def temperature_point(en: int) -> dict:
    tag = f"p2_cstr_t800_e{en}"
    outdir = BUILD_RUNS / tag
    return {"en": en, "tag": tag, "metrics": metrics(outdir, tag),
            "shares": source_fractions(outdir, tag)}


def electron_point(label: str, tag: str) -> dict:
    outdir = BUILD_RUNS / tag
    return {"factor": float(label), "tag": tag, "metrics": metrics(outdir, tag),
            "shares": source_fractions(outdir, tag)}


def summary_row(root: Path, branch: str) -> dict[str, str]:
    manifest = next((root / branch).rglob("nh_channel_summary.csv"), None)
    if manifest is None:
        raise FileNotFoundError(f"Missing summary for {branch}: {root}")
    with manifest.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"Empty summary: {manifest}")
    return rows[0]


def panel_label(ax, label: str) -> None:
    ax.text(-0.13, 1.15, label, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top", clip_on=False)


def save(fig: plt.Figure, output: Path, stem: str) -> None:
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(output / f"{stem}.png", dpi=600, bbox_inches="tight")


def make_figure1(output: Path) -> list[dict]:
    ens, n2s = [60, 120, 240], [0.1, 0.3333, 0.8]
    points = [p2_point(en, n2) for en in ens for n2 in n2s]
    by_key = {(p["en"], p["n2"]): p for p in points}
    nh3 = np.array([[by_key[(en, n2)]["metrics"]["NH3_final_cm-3"] for n2 in n2s] for en in ens])
    named = np.array([[by_key[(en, n2)]["shares"]["Named-H2*"] for n2 in n2s] for en in ens])

    causal_names = ["Base", "−Rydberg", "−named H$_2$*", "−excited N"]
    causal_keys = ["baseline", "no_rydberg_nh", "no_named_h2star_nh", "no_excited_n_nh"]
    causal = [float(summary_row(P3_CAUSAL, key)["NH3_final_cm-3"]) for key in causal_keys]
    causal_rel = 100 * np.asarray(causal) / causal[0]

    loss_mult = np.array([1, 10, 100, 1000], dtype=float)
    loss_rows = [summary_row(P3_RYDBERG, f"rydberg_loss_x{int(x)}") for x in loss_mult]
    ryd_share = []
    nh3_rel = []
    for row in loss_rows:
        branch = row["branch"]
        branch_rows_file = next((P3_RYDBERG / branch).rglob("nh_channel_summary.csv"), None)
        with branch_rows_file.open(encoding="utf-8", newline="") as fh:
            shares = {r["channel"]: float(r["share_pct"]) for r in csv.DictReader(fh)}
        ryd_share.append(shares["Rydberg-H2*"])
        nh3_rel.append(100 * float(row["NH3_final_cm-3"]) / causal[0])

    fig, axes = plt.subplots(2, 2, figsize=(8.35, 5.85), constrained_layout=True)
    ax = axes[0, 0]
    im = ax.imshow(np.log10(nh3), cmap="cividis", aspect="auto", vmin=9, vmax=15)
    ax.set_xticks(range(3), ["0.10", "0.3333", "0.80"])
    ax.set_yticks(range(3), ["60", "120", "240"])
    ax.set_xlabel("N$_2$ inlet fraction")
    ax.set_ylabel("E/N (Td)")
    ax.set_title("Steady NH$_3$ inventory (nine validation states)")
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{nh3[i, j]:.1e}", ha="center", va="center", fontsize=6,
                    color="white" if np.log10(nh3[i, j]) < 12.1 else "black")
    cbar = fig.colorbar(im, ax=ax, fraction=0.047, pad=0.04)
    cbar.set_label(r"log$_{10}$[NH$_3$] (cm$^{-3}$)")
    best = np.unravel_index(np.argmax(nh3), nh3.shape)
    ax.plot(best[1], best[0], marker="*", ms=10, mec="black", mfc="white", mew=0.8)
    panel_label(ax, "A")

    ax = axes[0, 1]
    im = ax.imshow(named, cmap="viridis", aspect="auto", vmin=90, vmax=100)
    ax.set_xticks(range(3), ["0.10", "0.3333", "0.80"])
    ax.set_yticks(range(3), ["60", "120", "240"])
    ax.set_xlabel("N$_2$ inlet fraction")
    ax.set_ylabel("E/N (Td)")
    ax.set_title("Direct NH source: named electronic H$_2$*")
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{named[i, j]:.1f}%", ha="center", va="center", fontsize=6,
                    color="black" if named[i, j] > 96 else "white")
    cbar = fig.colorbar(im, ax=ax, fraction=0.047, pad=0.04)
    cbar.set_label("Integrated NH source (%)")
    ax.text(0.03, 0.04, "All states: named H$_2$* is the largest\ndirect NH-entry family.",
            transform=ax.transAxes, fontsize=6.5, va="bottom",
            bbox={"facecolor": "white", "edgecolor": "#AAAAAA", "alpha": 0.9, "pad": 1.5})
    panel_label(ax, "B")

    ax = axes[1, 0]
    bars = ax.bar(range(len(causal_rel)), causal_rel,
                  color=[OKABE_ITO["blue"], OKABE_ITO["sky"], OKABE_ITO["vermillion"], OKABE_ITO["orange"]])
    ax.axhline(100, color=OKABE_ITO["black"], lw=0.7)
    ax.set_xticks(range(len(causal_rel)), causal_names, rotation=12, ha="right")
    ax.set_ylim(0, 115)
    ax.set_ylabel("Steady NH$_3$ relative to baseline (%)")
    ax.set_title("Causal pathway perturbations (120 Td, $x_{N_2}=0.3333$)")
    for bar, value in zip(bars, causal_rel):
        ax.text(bar.get_x() + bar.get_width()/2, value + 2, f"{value:.1f}", ha="center", va="bottom", fontsize=7)
    panel_label(ax, "C")

    ax = axes[1, 1]
    ax.plot(loss_mult, ryd_share, color=OKABE_ITO["purple"], marker="o", lw=1.5,
            label="Rydberg-H$_2$* NH share")
    ax.set_xscale("log")
    ax.set_xlabel("Added Rydberg wall-relaxation multiplier")
    ax.set_ylabel("Rydberg NH source (%)", color=OKABE_ITO["purple"])
    ax.tick_params(axis="y", labelcolor=OKABE_ITO["purple"])
    ax.set_title("Rydberg lifetime sensitivity (same reference state)")
    ax2 = ax.twinx()
    ax2.plot(loss_mult, nh3_rel, color=OKABE_ITO["black"], marker="s", ls="--", lw=1.2,
             label="Steady NH$_3$")
    ax2.set_ylabel("NH$_3$ relative to baseline (%)", color=OKABE_ITO["black"])
    ax2.tick_params(axis="y", labelcolor=OKABE_ITO["black"])
    ax.set_ylim(0, max(ryd_share) * 1.2)
    ax2.set_ylim(98.5, 100.2)
    panel_label(ax, "D")

    fig.text(0.5, -0.045,
             "CSTR/CW, $T_g$ = 300 K, $n_e$ = 1.17×10$^8$ cm$^{-3}$, $\\tau$ = 10 ms; "
             "direct-NH fluxes are integrated over the terminal three P0-steady cycles. ★ = largest inventory in panel A.",
             ha="center", va="bottom", fontsize=7)
    save(fig, output, "Figure_2_CSTR_mechanism_map")
    plt.close(fig)
    return points


def make_figure2(output: Path) -> None:
    ens = [60, 120, 240]
    temp = [temperature_point(en) for en in ens]
    electron = [electron_point("0.1", "p2_cstr_ne0p1"),
                electron_point("1", "p2_cstr_e120_n23333"),
                electron_point("10", "p2_cstr_ne10")]

    fig, axes = plt.subplots(1, 2, figsize=(8.35, 3.65), constrained_layout=True)
    ax = axes[0]
    for channel, color, marker in [("Named-H2*", OKABE_ITO["blue"], "o"),
                                   ("Rydberg-H2*", OKABE_ITO["purple"], "s"),
                                   ("N(2P)+H2", OKABE_ITO["orange"], "^"),
                                   ("Vib-H2", OKABE_ITO["green"], "D")]:
        values = [p["shares"][channel] for p in temp]
        ax.plot(ens, values, marker=marker, color=color, lw=1.5, label=channel)
    ax.set_yscale("log")
    ax.set_xlabel("E/N (Td)")
    ax.set_ylabel("Integrated NH source (%)")
    ax.set_title("800 K: no direct-NH pathway crossover")
    ax.set_xticks(ens)
    ax.set_ylim(1e-4, 150)
    ax.legend(frameon=False, loc="lower left")
    panel_label(ax, "A")

    ax = axes[1]
    x = np.array([p["factor"] for p in electron])
    y = np.array([p["metrics"]["NH3_final_cm-3"] for p in electron])
    ax.plot(x, y, marker="o", color=OKABE_ITO["blue"], lw=1.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(x, ["0.1", "1", "10"])
    ax.set_xlabel(r"Fixed $n_e$ / 1.17$\times$10$^8$")
    ax.set_ylabel(r"Steady [NH$_3$] (cm$^{-3}$)", color=OKABE_ITO["blue"])
    ax.tick_params(axis="y", labelcolor=OKABE_ITO["blue"])
    ax.set_title("Electron-density sensitivity (120 Td, $x_{N_2}=0.3333$)")
    ax2 = ax.twinx()
    named = [p["shares"]["Named-H2*"] for p in electron]
    ax2.plot(x, named, marker="s", ls="--", color=OKABE_ITO["vermillion"], lw=1.3)
    ax2.set_ylabel("Named H$_2$* NH source (%)", color=OKABE_ITO["vermillion"])
    ax2.tick_params(axis="y", labelcolor=OKABE_ITO["vermillion"])
    ax2.set_ylim(97.5, 100.0)
    panel_label(ax, "B")

    fig.text(0.5, -0.055,
             "Sensitivity checks retain the CSTR/P0 protocol; $n_e$ is an imposed model input, not a self-consistent plasma solution.",
             ha="center", va="bottom", fontsize=7)
    save(fig, output, "Figure_5_temperature_and_electron_robustness")
    plt.close(fig)


def write_data_table(output: Path, points: list[dict]) -> None:
    fields = ["en_td", "n2_fraction", "cycles_used", "nh3_final_cm-3", "named_h2star_pct",
              "rydberg_pct", "n2d_pct", "n2p_pct", "vib_pct"]
    with (output / "Figure_2_CSTR_mechanism_map_source_data.csv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for p in points:
            s, m = p["shares"], p["metrics"]
            writer.writerow({"en_td": p["en"], "n2_fraction": p["n2"],
                             "cycles_used": m["cycles_used"], "nh3_final_cm-3": m["NH3_final_cm-3"],
                             "named_h2star_pct": s["Named-H2*"], "rydberg_pct": s["Rydberg-H2*"],
                             "n2d_pct": s["N(2D)+H2"], "n2p_pct": s["N(2P)+H2"],
                             "vib_pct": s["Vib-H2"]})


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Generate validated Hong CSTR manuscript figures")
    parser.add_argument("--output", type=Path,
                        default=PROJECT / "analysis" / "figures_20260829")
    args = parser.parse_args(argv)
    configure_style()
    points = make_figure1(args.output)
    make_figure2(args.output)
    write_data_table(args.output, points)
    print(f"Figures written to: {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
