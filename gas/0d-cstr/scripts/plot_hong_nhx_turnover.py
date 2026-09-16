#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quantify terminal NHx source--sink turnover and CSTR outlet performance.

Only the 108 accepted CSTR/CW cases are processed.  For every elementary
reaction in the final three P0 cycles, signed stoichiometric contributions to
neutral NH, NH2, and NH3 are integrated.  This separates direct NH-source
topology from the subsequent source--sink balance that determines NH3 outlet
inventory.  The recorded electrical power is exported as a numerical audit
only: under fixed-ne operation it is explicitly not labelled an experimental
or physical energy efficiency.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd

from _paths import LEGACY_HONG_ROOT
from hong_reaction_hypergraph import PROJECT, RUN_ROOT, case_tag, parse_reaction


FIGURES = PROJECT / "analysis" / "figures_20260829"
SUMMARY = FIGURES / "Figure_3_dense_EN_N2_source_data.csv"
TARGETS = ("NH", "NH2", "NH3")
REPRESENTATIVES = [(20, 0.1, "Low field"), (140, 0.1, r"NH$_3$ maximum"),
                   (240, 0.9, r"High field, N$_2$-rich")]
COLORS = {"formation": "#009E73", "loss": "#D55E00", "net": "#0072B2", "gray": "#777777"}


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 9,
        "axes.titlesize": 9, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 6.8, "axes.linewidth": 0.7, "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def panel(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.14, 1.07, letter, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top")


def save(fig: plt.Figure, output: Path, name: str) -> None:
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(output / f"{name}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def species_coefficient(items: tuple[str, ...], target: str) -> int:
    """Coefficient of one neutral target species, excluding ionic relatives."""
    total = 0
    for item in items:
        match = re.match(r"^(\d+)?(.+)$", item.strip())
        coefficient, label = int(match.group(1) or "1"), match.group(2)
        if label == target:
            total += coefficient
    return total


def short_reaction(reaction: str) -> str:
    text = reaction.split(":", 1)[-1]
    text = text.replace("=>", " → ").replace("RYDBERG_SUM", "Ryd.")
    text = text.replace("H2(B3SIG)", "H2(B3)").replace("H2(B1SIG)", "H2(B1)")
    text = text.replace("H2(C3PI)", "H2(C3)").replace("H2(A3SIG)", "H2(A3)")
    return text


def terminal_power_and_inventory(case: Path, tag: str, n2_fraction: float) -> dict[str, float]:
    cycles_path = case / f"pulse_cycles_{tag}.csv"
    with cycles_path.open(encoding="utf-8", newline="") as handle:
        cycles = list(csv.DictReader(handle))
    last3 = cycles[-3:]
    time_values = [float(row["t_end_s"]) for row in last3]
    # CW period is also independently recorded as the last period increment.
    cycle_period = float(cycles[-1]["t_end_s"]) - float(cycles[-2]["t_end_s"])
    e_cycle = np.mean([float(row["E_on_Jcm3"]) + float(row["E_off_Jcm3"]) for row in last3])
    power = e_cycle / cycle_period if cycle_period > 0 else float("nan")

    series_path = case / f"pulse_series_{tag}.csv"
    with series_path.open(encoding="utf-8", newline="") as handle:
        series = list(csv.DictReader(handle))
    final = series[-1]
    ntot = 2.446e19
    n2_feed = n2_fraction * ntot
    n2_out = float(final["N2"])
    nh3_out = float(final["NH3"])
    tau = float(last3[-1]["tau_res_s"])
    productivity = nh3_out / tau if tau > 0 else float("nan")
    n2_conversion = max(0.0, (n2_feed - n2_out) / n2_feed) if n2_feed > 0 else float("nan")
    n_utilization = nh3_out / (2.0 * max(n2_feed - n2_out, 1e-300))
    apparent_per_j = productivity / power if power > 0 else float("nan")
    return {"cycle_period_s": cycle_period, "power_record_Wcm-3": power,
            "nh3_out_cm-3": nh3_out, "n2_out_cm-3": n2_out, "tau_res_s": tau,
            "nh3_productivity_cm-3s-1": productivity, "n2_conversion_frac": n2_conversion,
            "n_utilization_frac": n_utilization, "apparent_molecules_J_not_physical": apparent_per_j,
            "window_end_s": time_values[-1]}


def analyze_case(row: pd.Series) -> tuple[dict, list[dict]]:
    en_td, n2_fraction, tag = int(row.en_td), float(row.n2_fraction), str(row.tag)
    if str(row.run_status) not in {"success", "retry_success"}:
        raise ValueError(f"Non-accepted case requested: {tag}")
    case = RUN_ROOT / tag
    first, last = int(row.cycle_first), int(row.cycle_last)
    rates_path = case / f"pulse_rates_{tag}.csv"
    totals = {target: {"formation": 0.0, "loss": 0.0} for target in TARGETS}
    by_reaction: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    parsed = {}
    duration, time_token, parse_failures = 0.0, None, 0

    with rates_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            try:
                if int(raw["cycle"]) < first:
                    continue
                rate, dt = float(raw["rate_cm-3s-1"]), float(raw["dt_s"])
            except (KeyError, TypeError, ValueError):
                continue
            if not (math.isfinite(rate) and math.isfinite(dt) and dt > 0):
                continue
            token = (raw.get("cycle"), raw.get("time_s"))
            if token != time_token:
                duration += dt
                time_token = token
            reaction = (raw.get("reaction") or "").strip()
            if reaction not in parsed:
                parsed[reaction] = parse_reaction(reaction)
            equation = parsed[reaction]
            if equation is None:
                parse_failures += 1
                continue
            reactants, products = equation
            for target in TARGETS:
                delta = (species_coefficient(products, target) - species_coefficient(reactants, target)) * rate * dt
                if delta > 0:
                    totals[target]["formation"] += delta
                    by_reaction[reaction][f"{target}_formation"] += delta
                elif delta < 0:
                    totals[target]["loss"] += -delta
                    by_reaction[reaction][f"{target}_loss"] += -delta
    if duration <= 0:
        raise ValueError(f"No terminal duration in {rates_path}")
    out = {"en_td": en_td, "n2_fraction": n2_fraction, "tag": tag, "cycle_first": first,
           "cycle_last": last, "window_duration_s": duration, "parse_failures": parse_failures}
    for target in TARGETS:
        formation = totals[target]["formation"] / duration
        loss = totals[target]["loss"] / duration
        out[f"{target}_formation_cm-3s-1"] = formation
        out[f"{target}_loss_cm-3s-1"] = loss
        out[f"{target}_net_chem_cm-3s-1"] = formation - loss
        out[f"{target}_loss_to_formation"] = loss / formation if formation > 0 else float("nan")
    out.update(terminal_power_and_inventory(case, tag, n2_fraction))
    out["NH3_cstr_closure_rel"] = abs(out["NH3_net_chem_cm-3s-1"] - out["nh3_productivity_cm-3s-1"]) / max(
        abs(out["nh3_productivity_cm-3s-1"]), 1e4)

    reaction_rows: list[dict] = []
    for reaction, values in by_reaction.items():
        for target in TARGETS:
            for direction in ("formation", "loss"):
                integrated = values.get(f"{target}_{direction}", 0.0)
                if integrated > 0:
                    reaction_rows.append({"en_td": en_td, "n2_fraction": n2_fraction, "tag": tag,
                                          "target": target, "direction": direction, "reaction": reaction,
                                          "reaction_short": short_reaction(reaction),
                                          "flux_cm-3s-1": integrated / duration})
    return out, reaction_rows


def draw_turnover(ax: plt.Axes, reaction_table: pd.DataFrame, state: pd.Series, title: str) -> None:
    source = reaction_table.query("target == 'NH3' and direction == 'formation'").nlargest(3, "flux_cm-3s-1")
    sink = reaction_table.query("target == 'NH3' and direction == 'loss'").nlargest(3, "flux_cm-3s-1")
    entries = [(row["reaction_short"], row["flux_cm-3s-1"], "formation") for _, row in source.iterrows()] + \
              [(row["reaction_short"], -row["flux_cm-3s-1"], "loss") for _, row in sink.iterrows()]
    if not entries:
        ax.text(0.5, 0.5, "No NH$_3$ turnover records", ha="center", va="center")
        return
    labels, values, directions = zip(*entries)
    order = np.argsort(np.abs(values))
    labels = [labels[i] for i in order]; values = [values[i] for i in order]; directions = [directions[i] for i in order]
    colors = [COLORS[direction] for direction in directions]
    ax.barh(range(len(values)), values, color=colors)
    ax.axvline(0, color="black", lw=0.7)
    ax.set_yticks(range(len(values)), labels, fontsize=5.6)
    ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    ax.set_xlabel(r"Integrated reaction contribution (cm$^{-3}$ s$^{-1}$)")
    ax.set_title(title)
    ax.legend(handles=[mpl.patches.Patch(color=COLORS["formation"], label="formation"),
                       mpl.patches.Patch(color=COLORS["loss"], label="loss")], frameon=False, loc="lower right")
    ax.text(0.02, 0.03,
            rf"$P={state['NH3_formation_cm-3s-1']:.2e}$" + "\n" +
            rf"$D={state['NH3_loss_cm-3s-1']:.2e}$" + "\n" +
            rf"$D/P={state['NH3_loss_to_formation']:.3f}$",
            transform=ax.transAxes, va="bottom", fontsize=5.7,
            bbox={"facecolor": "white", "edgecolor": "#AAAAAA", "alpha": 0.92, "pad": 1.2})


def heatmap(ax: plt.Axes, data: pd.DataFrame, value: str, title: str, label: str,
            *, log: bool = True, cmap: str = "cividis") -> None:
    n2s, ens = sorted(data.n2_fraction.unique()), sorted(data.en_td.unique())
    matrix = data.pivot(index="en_td", columns="n2_fraction", values=value).loc[ens, n2s]
    payload = matrix.to_numpy(float)
    if log:
        payload = np.log10(np.clip(payload, 1e-300, None))
    im = ax.imshow(payload, origin="lower", aspect="auto", cmap=cmap)
    ax.set_xticks(range(len(n2s)), [f"{x:.1f}" for x in n2s])
    ax.set_yticks(range(0, len(ens), 2), [str(ens[i]) for i in range(0, len(ens), 2)])
    ax.set_xlabel(r"N$_2$ inlet fraction"); ax.set_ylabel("E/N (Td)"); ax.set_title(title)
    ax.set_xticks(np.arange(-0.5, len(n2s), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(ens), 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.32, alpha=0.72)
    ax.tick_params(which="minor", bottom=False, left=False)
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.047, pad=0.04); cb.set_label(label)


def make_figure8(summary: pd.DataFrame, reactions: pd.DataFrame, output: Path) -> None:
    fig = plt.figure(figsize=(11.65, 6.45), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, height_ratios=(1.0, 0.9))
    for ax, (en, n2, title), letter in zip([fig.add_subplot(grid[0, i]) for i in range(3)], REPRESENTATIVES, "ABC"):
        tag = case_tag(en, n2)
        state = summary.loc[summary.tag.eq(tag)].iloc[0]
        draw_turnover(ax, reactions.loc[reactions.tag.eq(tag)], state,
                      f"{title}: {en} Td, $x_{{N_2}}$={n2:.1f}")
        panel(ax, letter)
    ax = fig.add_subplot(grid[1, 0])
    heatmap(ax, summary, "NH3_formation_cm-3s-1", r"NH$_3$ formation flux", r"log$_{10}$ $P_{NH_3}$ (cm$^{-3}$ s$^{-1}$)")
    panel(ax, "D")
    ax = fig.add_subplot(grid[1, 1])
    heatmap(ax, summary, "NH3_loss_cm-3s-1", r"NH$_3$ destruction flux", r"log$_{10}$ $D_{NH_3}$ (cm$^{-3}$ s$^{-1}$)")
    panel(ax, "E")
    ax = fig.add_subplot(grid[1, 2])
    heatmap(ax, summary, "NH3_loss_to_formation", r"NH$_3$ loss / formation", r"log$_{10}$($D_{NH_3}/P_{NH_3}$)", cmap="magma")
    panel(ax, "F")
    fig.text(0.5, -0.045,
             "Top: three largest elementary sources and sinks at each stated operating point. Bottom: all 108 accepted CSTR/CW points; "
             "$P$ and $D$ are terminal-three-cycle reaction-of-progress fluxes.",
             ha="center", fontsize=7)
    save(fig, output, "Figure_8_NHx_source_sink_turnover")


def make_figure9(summary: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 6.75), constrained_layout=True)
    heatmap(axes[0, 0], summary, "nh3_productivity_cm-3s-1", r"CSTR NH$_3$ outlet productivity",
            r"log$_{10}$ $dot n_{NH_3,out}$ (cm$^{-3}$ s$^{-1}$)")
    panel(axes[0, 0], "A")
    best = summary.loc[summary["nh3_productivity_cm-3s-1"].idxmax()]
    n2s, ens = sorted(summary.n2_fraction.unique()), sorted(summary.en_td.unique())
    axes[0, 0].plot(n2s.index(best.n2_fraction), ens.index(best.en_td), marker="*", ms=10,
                     mec="black", mfc="white", mew=0.8)
    heatmap(axes[0, 1], summary, "n2_conversion_frac", r"CSTR N$_2$ conversion",
            r"log$_{10}$ conversion fraction")
    panel(axes[0, 1], "B")
    heatmap(axes[1, 0], summary, "n_utilization_frac", r"N-atom utilization to NH$_3$",
            r"log$_{10}$ $n_{NH_3}/[2(n_{N_2,in}-n_{N_2,out})]$")
    panel(axes[1, 0], "C")
    ax = axes[1, 1]
    points = ax.scatter(summary["n2_conversion_frac"], summary["nh3_productivity_cm-3s-1"],
                        c=summary["en_td"], s=22 + 30 * summary["n2_fraction"], cmap="cividis",
                        edgecolor="black", linewidth=0.25)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"N$_2$ conversion fraction"); ax.set_ylabel(r"NH$_3$ outlet productivity (cm$^{-3}$ s$^{-1}$)")
    ax.set_title("Performance frontier within the CSTR model")
    cb = fig.colorbar(points, ax=ax, fraction=0.047, pad=0.04); cb.set_label("E/N (Td)")
    panel(ax, "D")
    ax.plot(best["n2_conversion_frac"], best["nh3_productivity_cm-3s-1"], marker="*", ms=10,
            mec="black", mfc="white", mew=0.8, zorder=5)
    ax.text(0.03, 0.04, "★ global productivity maximum\nsize = N$_2$ inlet fraction",
            transform=ax.transAxes, fontsize=6.5, bbox={"facecolor": "white", "edgecolor": "#AAAAAA", "alpha": 0.9, "pad": 1.2})
    fig.text(0.5, -0.045,
             "All panels show 108 accepted CSTR/CW calculations. N-atom utilization is a material-balance descriptor; "
             "fixed-$n_e$ power is not used to claim physical energy efficiency.",
             ha="center", fontsize=7)
    save(fig, output, "Figure_9_CSTR_productivity_and_conversion")


def main() -> int:
    parser = argparse.ArgumentParser(description="Hong terminal NHx source-sink and CSTR performance analysis")
    parser.add_argument("--output", type=Path, default=FIGURES)
    parser.add_argument("--limit", type=int, default=0, help="process first N summary rows; 0 means all")
    parser.add_argument("--plot-only", action="store_true", help="render figures from existing 108-point CSV outputs")
    args = parser.parse_args()
    configure()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.plot_only:
        summary = pd.read_csv(args.output / "Figure_8_9_NHx_CSTR_summary.csv")
        reactions = pd.read_csv(args.output / "Figure_8_NHx_reaction_turnover.csv")
        if len(summary) != 108:
            raise ValueError("--plot-only requires the complete 108-point summary")
        make_figure8(summary, reactions, args.output)
        make_figure9(summary, args.output)
        return 0
    dense = pd.read_csv(SUMMARY)
    dense = dense.loc[dense.run_status.isin(["success", "retry_success"])].sort_values(["en_td", "n2_fraction"])
    if args.limit:
        dense = dense.head(args.limit)
    summaries, reaction_rows = [], []
    for number, (_, row) in enumerate(dense.iterrows(), start=1):
        result, reactions = analyze_case(row)
        summaries.append(result); reaction_rows.extend(reactions)
        print(f"[{number}/{len(dense)}] {result['tag']} closure={result['NH3_cstr_closure_rel']:.3e}")
    summary = pd.DataFrame(summaries)
    reactions = pd.DataFrame(reaction_rows)
    summary.to_csv(args.output / "Figure_8_9_NHx_CSTR_summary.csv", index=False)
    reactions.to_csv(args.output / "Figure_8_NHx_reaction_turnover.csv", index=False)
    if len(summary) == 108:
        make_figure8(summary, reactions, args.output)
        make_figure9(summary, args.output)
    else:
        print("Partial run: CSV outputs only; figures require all 108 accepted cases.")
    print(summary[["NH3_cstr_closure_rel", "apparent_molecules_J_not_physical"]].describe().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
