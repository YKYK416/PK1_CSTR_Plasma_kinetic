#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Summarise and plot the completed Hong closed-0D finite-exposure campaign.

The input is the active completed attempt for every case in ``cases_overview``.
Outputs deliberately describe 100 s *closed-0D finite exposure* results: they
are not CSTR outlet quantities, per-pass yields, or steady-state claims.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
SELECTED_X_N2 = (0.1, 0.5, 0.9)
REPRESENTATIVE_EN = (20.0, 80.0, 140.0, 240.0)
REPRESENTATIVE_T = 400.0
REPRESENTATIVE_X_N2 = 0.5


def read_series(path: Path) -> list[dict[str, float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [{name: float(value.replace("D", "E").replace("d", "e"))
                 for name, value in row.items()}
                for row in csv.DictReader(handle)]


def style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 9, "axes.titlesize": 9,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.linewidth": 0.8, "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })


def despine(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=350, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def nearest(rows: list[dict[str, float]], time_s: float) -> dict[str, float]:
    return min(rows, key=lambda row: abs(row["time_s"] - time_s))


def load_campaign(root: Path) -> tuple[list[dict[str, object]], dict[tuple[float, float, float], list[dict[str, float]]]]:
    overview = root / "cases_overview.csv"
    cases: list[dict[str, object]] = []
    trajectories: dict[tuple[float, float, float], list[dict[str, float]]] = {}
    with overview.open(encoding="utf-8", newline="") as handle:
        for ledger in csv.DictReader(handle):
            if ledger["status"] != "completed":
                raise RuntimeError(f"campaign contains non-completed case: {ledger['case_id']} ({ledger['status']})")
            output = Path(ledger["output_dir"])
            params_path = output / "params.json"
            if not params_path.is_file():
                raise RuntimeError(f"missing active params.json: {params_path}")
            params = json.loads(params_path.read_text(encoding="utf-8"))
            if params.get("returncode") != 0 or params.get("timed_out"):
                raise RuntimeError(f"completed ledger points to unsuccessful run: {ledger['case_id']}")
            tag = str(params["tag"])
            series_path = output / f"cw_surface_{tag}.csv"
            rows = read_series(series_path)
            final, t80, t10 = rows[-1], nearest(rows, 80.0), nearest(rows, 10.0)
            if abs(final["time_s"] - 100.0) > 1.0e-6:
                raise RuntimeError(f"{ledger['case_id']} does not end at 100 s")
            nh3_final = final["NH3_cm-3"]
            time_to_90pct = next((row["time_s"] for row in rows if row["NH3_cm-3"] >= 0.9 * nh3_final), math.nan)
            key = (float(ledger["temperature_K"]), float(ledger["en_Td"]), float(ledger["n2_fraction"]))
            trajectories[key] = rows
            cases.append({
                "case_id": ledger["case_id"], "temperature_K": key[0], "en_Td": key[1], "x_N2": key[2],
                "x_H2": float(ledger["h2_fraction"]), "phase": ledger["phase"],
                "attempt": int(ledger["active_attempt"]), "output_dir": str(output),
                "NH3_80s_cm-3": t80["NH3_cm-3"], "NH3_100s_cm-3": nh3_final,
                "NH3_mean_net_rate_80_100_cm-3_s-1": (nh3_final - t80["NH3_cm-3"]) / (final["time_s"] - t80["time_s"]),
                "NH3_relative_change_80_100": (nh3_final - t80["NH3_cm-3"]) / max(abs(nh3_final), 1.0),
                "NH3_time_to_90pct_final_s": time_to_90pct,
                "N_100s_cm-3": final["N_cm-3"], "H_100s_cm-3": final["H_cm-3"],
                "NH_100s_cm-3": final["NH_cm-3"], "NH2_100s_cm-3": final["NH2_cm-3"],
                "theta_Surf_100s": final["theta_Surf"], "theta_HSurf_100s": final["theta_HSurf"],
                "theta_NSurf_100s": final["theta_NSurf"], "theta_NHSurf_100s": final["theta_NHSurf"],
                "theta_NH2Surf_100s": final["theta_NH2Surf"],
                "theta_Surf_10s": t10["theta_Surf"], "theta_HSurf_10s": t10["theta_HSurf"],
                "theta_HSurf_change_10_100": final["theta_HSurf"] - t10["theta_HSurf"],
                "site_balance_error_max": max(abs(row["site_balance_error"]) for row in rows),
                "site_projection": bool(params.get("site_projection", False)),
                "preprojection_site_balance_error_max": (params.get("solver_completion") or {}).get("site_projection_max_rel"),
            })
    if len(cases) != 972 or len(trajectories) != 972:
        raise RuntimeError(f"expected 972 unique completed cases, got {len(cases)} records / {len(trajectories)} trajectories")
    return cases, trajectories


def matrix(cases: list[dict[str, object]], x_n2: float, field: str,
           temperatures: list[float], fields: list[float]) -> np.ndarray:
    lookup = {(float(row["temperature_K"]), float(row["en_Td"]), float(row["x_N2"])): float(row[field]) for row in cases}
    return np.array([[lookup[(temperature, en, x_n2)] for en in fields] for temperature in temperatures], dtype=float)


def heatmap_triptych(cases: list[dict[str, object]], temperatures: list[float], fields: list[float],
                     value_field: str, label: str, stem: Path, transform) -> None:
    values = [transform(matrix(cases, x_n2, value_field, temperatures, fields)) for x_n2 in SELECTED_X_N2]
    vmin, vmax = min(np.ma.min(item) for item in values), max(np.ma.max(item) for item in values)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#d9d9d9")
    fig, axes = plt.subplots(1, 3, figsize=(7.25, 2.5), sharey=True, constrained_layout=True)
    image = None
    for index, (ax, x_n2, value) in enumerate(zip(axes, SELECTED_X_N2, values)):
        ax.set_facecolor("#d9d9d9")
        image = ax.imshow(value, origin="lower", aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                          extent=(min(fields) - 10, max(fields) + 10, min(temperatures) - 12.5, max(temperatures) + 12.5))
        ax.set_title(rf"$x_{{N_2}}$ = {x_n2:.1f}")
        ax.set_xlabel(r"E/N (Td)")
        ax.set_xticks(fields[::2])
        if index == 0:
            ax.set_ylabel("Gas temperature (K)")
        despine(ax)
    assert image is not None
    colorbar = fig.colorbar(image, ax=axes, shrink=0.9, pad=0.02)
    colorbar.set_label(label)
    save(fig, stem)


def surface_heatmaps(cases: list[dict[str, object]], temperatures: list[float], fields: list[float], stem: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(7.25, 4.25), sharex=True, sharey=True, constrained_layout=True)
    configs = (("theta_Surf_100s", r"$\theta_{\mathrm{Surf}}$ at 100 s"),
               ("theta_HSurf_100s", r"$\theta_{\mathrm{HSurf}}$ at 100 s"))
    for row_index, (field, label) in enumerate(configs):
        values = [matrix(cases, x_n2, field, temperatures, fields) for x_n2 in SELECTED_X_N2]
        vmin, vmax = min(np.min(value) for value in values), max(np.max(value) for value in values)
        image = None
        for col_index, (ax, x_n2, value) in enumerate(zip(axes[row_index], SELECTED_X_N2, values)):
            image = ax.imshow(value, origin="lower", aspect="auto", cmap="cividis", vmin=vmin, vmax=vmax,
                              extent=(min(fields) - 10, max(fields) + 10, min(temperatures) - 12.5, max(temperatures) + 12.5))
            if row_index == 0:
                ax.set_title(rf"$x_{{N_2}}$ = {x_n2:.1f}")
            if row_index == 1:
                ax.set_xlabel(r"E/N (Td)")
                ax.set_xticks(fields[::2])
            if col_index == 0:
                ax.set_ylabel("Gas temperature (K)")
            despine(ax)
        assert image is not None
        colorbar = fig.colorbar(image, ax=axes[row_index], shrink=0.86, pad=0.02)
        colorbar.set_label(label)
    save(fig, stem)


def representative_trajectories(trajectories: dict[tuple[float, float, float], list[dict[str, float]]], stem: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 2.65), constrained_layout=True)
    for color, en in zip(OKABE_ITO, REPRESENTATIVE_EN):
        rows = trajectories[(REPRESENTATIVE_T, en, REPRESENTATIVE_X_N2)]
        axes[0].plot([row["time_s"] for row in rows], [row["NH3_cm-3"] for row in rows],
                     color=color, lw=1.8, marker="o", markersize=2.4, label=f"{en:.0f} Td")
        axes[1].plot([row["time_s"] for row in rows], [row["theta_HSurf"] for row in rows],
                     color=color, lw=1.8, marker="o", markersize=2.4, label=f"{en:.0f} Td")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Exposure time (s)")
    axes[0].set_ylabel(r"[NH$_3$] (cm$^{-3}$)")
    axes[0].set_title(rf"NH$_3$; {REPRESENTATIVE_T:.0f} K, $x_{{N_2}}$={REPRESENTATIVE_X_N2:.1f}")
    axes[1].set_xlabel("Exposure time (s)")
    axes[1].set_ylabel(r"$\theta_{\mathrm{{HSurf}}}$")
    axes[1].set_title(rf"H-covered sites; {REPRESENTATIVE_T:.0f} K, $x_{{N_2}}$={REPRESENTATIVE_X_N2:.1f}")
    for ax in axes:
        ax.legend(title="E/N", frameon=False, ncol=2, loc="best")
        ax.set_xlim(0, 100)
        despine(ax)
    save(fig, stem)


def cross_sections(cases: list[dict[str, object]], stem: Path) -> None:
    lookup = {(float(row["temperature_K"]), float(row["en_Td"]), float(row["x_N2"])): float(row["NH3_100s_cm-3"]) for row in cases}
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 2.65), constrained_layout=True)
    for color, x_n2 in zip(OKABE_ITO, [0.1, 0.3, 0.5, 0.7, 0.9]):
        fields = list(range(20, 241, 20))
        values = [lookup[(400.0, float(en), x_n2)] for en in fields]
        axes[0].plot(fields, values, color=color, lw=1.7, marker="o", markersize=2.5, label=rf"$x_{{N_2}}$={x_n2:.1f}")
    for color, en in zip(OKABE_ITO, [40.0, 80.0, 140.0, 200.0]):
        compositions = [round(index / 10, 1) for index in range(1, 10)]
        values = [lookup[(400.0, en, x_n2)] for x_n2 in compositions]
        axes[1].plot(compositions, values, color=color, lw=1.7, marker="o", markersize=2.5, label=f"{en:.0f} Td")
    axes[0].set_xlabel(r"E/N (Td)")
    axes[0].set_title(r"Composition sweep at 400 K")
    axes[1].set_xlabel(r"$x_{N_2}$ ($x_{H_2}=1-x_{N_2}$)")
    axes[1].set_title(r"Field sweep at 400 K")
    for ax in axes:
        ax.set_yscale("log")
        ax.set_ylabel(r"[NH$_3$] at 100 s (cm$^{-3}$)")
        ax.legend(frameon=False, ncol=2, loc="best")
        despine(ax)
    save(fig, stem)


FACTORS = ("temperature_K", "en_Td", "x_N2")
FACTOR_LABELS = {"temperature_K": "Gas temperature (K)", "en_Td": "E/N (Td)", "x_N2": r"$x_{N_2}$"}
VARIANCE_RESPONSES = {
    "log10_NH3_100s": r"log$_{10}$([NH$_3$] at 100 s / cm$^{-3}$)",
    "log10_N_100s": r"log$_{10}$([N] at 100 s / cm$^{-3}$)",
    "log10_H_100s": r"log$_{10}$([H] at 100 s / cm$^{-3}$)",
    "theta_Surf_100s": r"$\theta_{\mathrm{Surf}}$ at 100 s",
}


def response_value(row: dict[str, object], response: str) -> float:
    if response == "log10_NH3_100s":
        return math.log10(max(float(row["NH3_100s_cm-3"]), 1.0))
    if response == "log10_N_100s":
        return math.log10(max(float(row["N_100s_cm-3"]), 1.0))
    if response == "log10_H_100s":
        return math.log10(max(float(row["H_100s_cm-3"]), 1.0))
    return float(row[response])


def grouped_values(cases: list[dict[str, object]], factor: str, response: str) -> list[dict[str, object]]:
    groups: dict[float, list[float]] = defaultdict(list)
    for row in cases:
        groups[float(row[factor])].append(response_value(row, response))
    return [{"factor": factor, "factor_value": value, "response": response, "n": len(values),
             "mean": float(np.mean(values)), "min": float(np.min(values)), "max": float(np.max(values)),
             "range": float(np.max(values) - np.min(values))}
            for value, values in sorted(groups.items())]


def factorial_variance_partition(cases: list[dict[str, object]], response: str) -> list[dict[str, object]]:
    """Balanced 3-factor sum-of-squares decomposition; not an inferential ANOVA."""
    y = np.array([response_value(row, response) for row in cases], dtype=float)
    grand = float(np.mean(y))
    level = {factor: sorted({float(row[factor]) for row in cases}) for factor in FACTORS}
    main: dict[str, dict[float, float]] = {}
    for factor in FACTORS:
        main[factor] = {value: float(np.mean([response_value(row, response) for row in cases if float(row[factor]) == value]))
                        for value in level[factor]}
    pair: dict[tuple[str, str], dict[tuple[float, float], float]] = {}
    for first_index, first in enumerate(FACTORS):
        for second in FACTORS[first_index + 1:]:
            pair[(first, second)] = {}
            for a in level[first]:
                for b in level[second]:
                    pair[(first, second)][(a, b)] = float(np.mean([
                        response_value(row, response) for row in cases
                        if float(row[first]) == a and float(row[second]) == b
                    ]))
    main_prediction = np.array([
        grand + sum(main[factor][float(row[factor])] - grand for factor in FACTORS)
        for row in cases
    ])
    pair_prediction = np.zeros(len(cases), dtype=float)
    for (first, second), values in pair.items():
        pair_prediction += np.array([
            values[(float(row[first]), float(row[second]))] - main[first][float(row[first])]
            - main[second][float(row[second])] + grand for row in cases
        ])
    total_ss = float(np.sum((y - grand) ** 2))
    pieces: list[tuple[str, float]] = []
    for factor in FACTORS:
        pieces.append((FACTOR_LABELS[factor], float(np.sum(np.array([
            main[factor][float(row[factor])] - grand for row in cases
        ]) ** 2))))
    for (first, second), values in pair.items():
        pieces.append((f"{FACTOR_LABELS[first]} × {FACTOR_LABELS[second]}", float(np.sum(np.array([
            values[(float(row[first]), float(row[second]))] - main[first][float(row[first])]
            - main[second][float(row[second])] + grand for row in cases
        ]) ** 2))))
    pieces.append(("Three-factor residual", float(np.sum((y - main_prediction - pair_prediction) ** 2))))
    return [{"response": response, "response_label": VARIANCE_RESPONSES[response], "component": name,
             "sum_squares": value, "fraction_total_ss": value / total_ss if total_ss else math.nan,
             "total_sum_squares": total_ss}
            for name, value in pieces]


def main_effects_plot(effects: list[dict[str, object]], stem: Path) -> None:
    subset = [row for row in effects if row["response"] == "log10_NH3_100s"]
    fig, axes = plt.subplots(1, 3, figsize=(7.25, 2.5), constrained_layout=True)
    for panel, (ax, factor) in enumerate(zip(axes, FACTORS)):
        rows = [row for row in subset if row["factor"] == factor]
        x = [float(row["factor_value"]) for row in rows]
        mean = [float(row["mean"]) for row in rows]
        low = [float(row["min"]) for row in rows]
        high = [float(row["max"]) for row in rows]
        ax.fill_between(x, low, high, color="#56B4E9", alpha=0.30, label="min–max across other factors")
        ax.plot(x, mean, color="#0072B2", lw=1.8, marker="o", markersize=3, label="grid mean")
        ax.set_xlabel(FACTOR_LABELS[factor])
        ax.set_ylabel(r"log$_{10}$([NH$_3$] at 100 s / cm$^{-3}$)")
        ax.text(-0.16, 1.04, "ABC"[panel], transform=ax.transAxes, fontweight="bold", fontsize=10)
        if panel == 0:
            ax.legend(frameon=False, fontsize=6, loc="best")
        despine(ax)
    save(fig, stem)


def variance_partition_plot(records: list[dict[str, object]], stem: Path) -> None:
    responses = list(VARIANCE_RESPONSES)
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 4.75), constrained_layout=True)
    colors = ["#0072B2", "#56B4E9", "#009E73", "#E69F00", "#D55E00", "#CC79A7", "#777777"]
    for panel, (ax, response) in enumerate(zip(axes.flat, responses)):
        rows = [row for row in records if row["response"] == response]
        values = [100.0 * float(row["fraction_total_ss"]) for row in rows]
        labels = [str(row["component"]).replace("Gas temperature (K)", "T").replace("E/N (Td)", "E/N").replace(r"$x_{N_2}$", "xN2") for row in rows]
        ax.bar(range(len(rows)), values, color=colors, edgecolor="none")
        ax.set_xticks(range(len(rows)), labels, rotation=45, ha="right")
        ax.set_ylabel("Response SS contribution (%)")
        ax.set_title(VARIANCE_RESPONSES[response])
        ax.text(-0.02, 1.04, "ABCD"[panel], transform=ax.transAxes, fontweight="bold", fontsize=10)
        ax.set_ylim(0, max(values) * 1.18 if max(values) else 1)
        despine(ax)
    save(fig, stem)


def nh3_interaction_maps(cases: list[dict[str, object]], fields: list[float], compositions: list[float], stem: Path) -> None:
    temperatures = (300.0, 400.0, 500.0)
    lookup = {(float(row["temperature_K"]), float(row["en_Td"]), float(row["x_N2"])): response_value(row, "log10_NH3_100s") for row in cases}
    matrices = [np.array([[lookup[(temperature, en, x_n2)] for en in fields] for x_n2 in compositions]) for temperature in temperatures]
    vmin, vmax = min(np.min(value) for value in matrices), max(np.max(value) for value in matrices)
    fig, axes = plt.subplots(1, 3, figsize=(7.25, 2.7), sharey=True, constrained_layout=True)
    image = None
    for panel, (ax, temperature, value) in enumerate(zip(axes, temperatures, matrices)):
        image = ax.imshow(value, origin="lower", aspect="auto", cmap="viridis", vmin=vmin, vmax=vmax,
                          extent=(min(fields) - 10, max(fields) + 10, min(compositions) - .05, max(compositions) + .05))
        ax.set_title(f"{temperature:.0f} K")
        ax.set_xlabel("E/N (Td)")
        ax.set_xticks(fields[::2])
        if panel == 0:
            ax.set_ylabel(r"$x_{N_2}$")
        ax.text(-0.16, 1.04, "ABC"[panel], transform=ax.transAxes, fontweight="bold", fontsize=10)
        despine(ax)
    assert image is not None
    colorbar = fig.colorbar(image, ax=axes, shrink=.90, pad=.02)
    colorbar.set_label(r"log$_{10}$([NH$_3$] at 100 s / cm$^{-3}$)")
    save(fig, stem)


def gas_species_maps(cases: list[dict[str, object]], fields: list[float], compositions: list[float], stem: Path) -> None:
    temperature = 400.0
    metrics = (("N_100s_cm-3", "N"), ("H_100s_cm-3", "H"), ("NH2_100s_cm-3", "NH$_2$"), ("NH3_100s_cm-3", "NH$_3$"))
    lookup = {(float(row["temperature_K"]), float(row["en_Td"]), float(row["x_N2"])): row for row in cases}
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 4.7), sharex=True, sharey=True, constrained_layout=True)
    for panel, (ax, (field, label)) in enumerate(zip(axes.flat, metrics)):
        values = np.array([[math.log10(max(float(lookup[(temperature, en, x_n2)][field]), 1.0)) for en in fields] for x_n2 in compositions])
        image = ax.imshow(values, origin="lower", aspect="auto", cmap="plasma",
                          extent=(min(fields) - 10, max(fields) + 10, min(compositions) - .05, max(compositions) + .05))
        ax.set_title(f"{label}, {temperature:.0f} K")
        ax.set_xlabel("E/N (Td)")
        ax.set_ylabel(r"$x_{N_2}$")
        ax.text(-0.16, 1.04, "ABCD"[panel], transform=ax.transAxes, fontweight="bold", fontsize=10)
        colorbar = fig.colorbar(image, ax=ax, shrink=.82)
        colorbar.set_label(r"log$_{10}$ concentration (cm$^{-3}$)")
        despine(ax)
    save(fig, stem)


def temporal_feature_maps(cases: list[dict[str, object]], fields: list[float], compositions: list[float], stem: Path) -> None:
    temperatures = (300.0, 400.0, 500.0)
    lookup = {(float(row["temperature_K"]), float(row["en_Td"]), float(row["x_N2"])): row for row in cases}
    time_values = [np.array([[float(lookup[(temperature, en, x_n2)]["NH3_time_to_90pct_final_s"]) for en in fields] for x_n2 in compositions]) for temperature in temperatures]
    change_values = [np.ma.log10(np.ma.masked_less_equal(np.array([[float(lookup[(temperature, en, x_n2)]["NH3_relative_change_80_100"]) for en in fields] for x_n2 in compositions]), 0.0)) for temperature in temperatures]
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#d9d9d9")
    fig, axes = plt.subplots(2, 3, figsize=(7.25, 4.5), sharex=True, sharey=True, constrained_layout=True)
    time_image = change_image = None
    for panel, (temperature, value) in enumerate(zip(temperatures, time_values)):
        time_image = axes[0, panel].imshow(value, origin="lower", aspect="auto", cmap="cividis", vmin=0, vmax=100,
                                           extent=(min(fields) - 10, max(fields) + 10, min(compositions) - .05, max(compositions) + .05))
        axes[0, panel].set_title(f"{temperature:.0f} K")
        axes[0, panel].text(-0.16, 1.04, "ABC"[panel], transform=axes[0, panel].transAxes, fontweight="bold", fontsize=10)
    vmin = min(np.ma.min(value) for value in change_values)
    vmax = max(np.ma.max(value) for value in change_values)
    for panel, value in enumerate(change_values):
        axes[1, panel].set_facecolor("#d9d9d9")
        change_image = axes[1, panel].imshow(value, origin="lower", aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                                              extent=(min(fields) - 10, max(fields) + 10, min(compositions) - .05, max(compositions) + .05))
        axes[1, panel].set_xlabel("E/N (Td)")
        axes[1, panel].text(-0.16, 1.04, "DEF"[panel], transform=axes[1, panel].transAxes, fontweight="bold", fontsize=10)
    for row in axes:
        row[0].set_ylabel(r"$x_{N_2}$")
        for ax in row:
            despine(ax)
    assert time_image is not None and change_image is not None
    time_bar = fig.colorbar(time_image, ax=axes[0], shrink=.86, pad=.02)
    time_bar.set_label("Time to 90% final [NH$_3$] (s)")
    change_bar = fig.colorbar(change_image, ax=axes[1], shrink=.86, pad=.02)
    change_bar.set_label(r"log$_{10}$ relative [NH$_3$] change, 80–100 s")
    save(fig, stem)


def top_conditions(cases: list[dict[str, object]]) -> list[dict[str, object]]:
    metrics = ("NH3_100s_cm-3", "NH3_mean_net_rate_80_100_cm-3_s-1", "N_100s_cm-3", "theta_Surf_100s")
    records: list[dict[str, object]] = []
    for metric in metrics:
        for rank, row in enumerate(sorted(cases, key=lambda item: float(item[metric]), reverse=True)[:20], start=1):
            records.append({"metric": metric, "rank": rank, "case_id": row["case_id"], "temperature_K": row["temperature_K"],
                            "en_Td": row["en_Td"], "x_N2": row["x_N2"], "x_H2": row["x_H2"], "value": row[metric]})
    return records


def write_summary(path: Path, cases: list[dict[str, object]]) -> None:
    fields = list(cases[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(cases)


def write_records(path: Path, records: list[dict[str, object]]) -> None:
    fields = sorted({field for record in records for field in record})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot final results from completed Hong closed-0D campaign")
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root, output = args.campaign_root.resolve(), args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    figures = output / "figures"
    figures.mkdir(exist_ok=True)
    style()
    cases, trajectories = load_campaign(root)
    temperatures = sorted({float(row["temperature_K"]) for row in cases})
    fields = sorted({float(row["en_Td"]) for row in cases})
    compositions = sorted({float(row["x_N2"]) for row in cases})
    write_summary(output / "hong_closed0d_100s_summary.csv", cases)
    write_summary(output / "systematic_response_metrics.csv", cases)
    main_effect_records = [record for response in VARIANCE_RESPONSES for factor in FACTORS
                           for record in grouped_values(cases, factor, response)]
    variance_records = [record for response in VARIANCE_RESPONSES
                        for record in factorial_variance_partition(cases, response)]
    write_records(output / "factor_main_effects.csv", main_effect_records)
    write_records(output / "factorial_variance_partition.csv", variance_records)
    write_records(output / "top_conditions_by_response.csv", top_conditions(cases))
    heatmap_triptych(cases, temperatures, fields, "NH3_100s_cm-3", r"log$_{10}$([NH$_3$] at 100 s / cm$^{-3}$)",
                     figures / "figure_01_nh3_endpoint_maps", lambda value: np.log10(np.maximum(value, 1.0)))
    heatmap_triptych(cases, temperatures, fields, "NH3_mean_net_rate_80_100_cm-3_s-1",
                     r"log$_{10}$(mean d[NH$_3$]/dt, 80–100 s / cm$^{-3}$ s$^{-1}$)",
                     figures / "figure_02_nh3_terminal_rate_maps",
                     lambda value: np.ma.log10(np.ma.masked_less_equal(value, 0.0)))
    surface_heatmaps(cases, temperatures, fields, figures / "figure_03_surface_coverage_maps")
    representative_trajectories(trajectories, figures / "figure_04_representative_time_traces")
    cross_sections(cases, figures / "figure_05_nh3_cross_sections_400K")
    main_effects_plot(main_effect_records, figures / "figure_06_nh3_main_effects")
    variance_partition_plot(variance_records, figures / "figure_07_factorial_variance_partition")
    nh3_interaction_maps(cases, fields, compositions, figures / "figure_08_nh3_field_composition_interactions")
    gas_species_maps(cases, fields, compositions, figures / "figure_09_gas_species_maps_400K")
    temporal_feature_maps(cases, fields, compositions, figures / "figure_10_temporal_features")
    ordered = sorted(cases, key=lambda row: float(row["NH3_100s_cm-3"]), reverse=True)
    findings = {
        "model_scope": "Hong 2017/2018 corrected mechanism; closed 0D; continuous E/N; 100 s finite exposure; no CSTR; no pulse.",
        "case_count": len(cases), "temperatures_K": temperatures, "en_Td": fields, "x_N2": compositions,
        "maximum_NH3_100s": ordered[0], "minimum_NH3_100s": ordered[-1],
        "site_balance_error_max_over_campaign": max(float(row["site_balance_error_max"]) for row in cases),
        "site_projection_case_count": sum(bool(row["site_projection"]) for row in cases),
        "systematic_analysis": {
            "method": "balanced full-factorial response sum-of-squares partition; deterministic grid, not inferential ANOVA",
            "responses": VARIANCE_RESPONSES,
            "main_effect_records": len(main_effect_records), "variance_partition_records": len(variance_records),
        },
        "figures": [str(path.name) for path in sorted(figures.glob("*.pdf"))],
    }
    (output / "analysis_findings.json").write_text(json.dumps(findings, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(findings, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
