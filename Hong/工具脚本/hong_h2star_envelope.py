#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Mechanism-scenario envelope for the Hong CW-CSTR direct-NH hierarchy.

This is deliberately a transparent deterministic envelope, not a statistical
uncertainty quantification: no probability distributions or confidence bounds
are assigned to rate constants. Three mechanistically salient reaction groups
are independently scaled by 1/3 and 3 at two contrasting accepted CSTR states.
Every scenario must independently satisfy the terminal-window P0 gate.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

import zdp_gui as gui
from hong_nh_diagnostics import PRIMARY_PATTERNS, integrate_nh_channels, steady_window_start


TOOL_DIR = Path(__file__).resolve().parent
PROJECT = TOOL_DIR.parent
BASE = PROJECT / "Reproduction" / "2017Hong" / "方案2_脉冲能效地图" / "build_pulse"
DEFAULT_OUTPUT = PROJECT / "analysis" / "p9_h2star_scenario_envelope_20260831"

STATES = (
    ("productive", "140", "0.1", "Productivity maximum: 140 Td, x_N2 = 0.1"),
    ("loss_dominated", "240", "0.9", "Loss dominated: 240 Td, x_N2 = 0.9"),
)
CONTROLS = {
    "named_H2star_to_NH": {
        "label": "Named H2* to NH",
        "needles": tuple(f"N + H2({state}) => H + NH" for state in ("B3SIG", "B1SIG", "C3PI", "A3SIG")),
    },
    "named_H2star_wall_loss": {
        "label": "Named H2* wall relaxation",
        "needles": tuple(f"H2({state}) => H2" for state in ("B3SIG", "B1SIG", "C3PI", "A3SIG")),
    },
    "rydberg_to_NH": {
        "label": "Rydberg-H2* to NH",
        "needles": ("N + H2(RYDBERG_SUM) => H + NH",),
    },
}
FACTORS = (1.0 / 3.0, 3.0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_build_inputs(base: Path, target: Path) -> None:
    target.mkdir(parents=True, exist_ok=False)
    keep_names = {"preprocessor.exe", "bolsigdb.dat", "main_pulse.F90", "dvode_f90_m.F90"}
    keep_suffixes = {".dll", ".lib", ".dat", ".DAT"}
    for source in base.iterdir():
        if source.is_file() and (source.name in keep_names or source.suffix in keep_suffixes):
            shutil.copy2(source, target / source.name)


def scale_group(kinet_text: str, control: str, factor: float) -> tuple[str, list[str]]:
    spec = CONTROLS[control]
    counts = {needle: 0 for needle in spec["needles"]}
    output: list[str] = []
    for line in kinet_text.splitlines(keepends=True):
        reaction = line.split("!", 1)[0].strip()
        compact = re.sub(r"\s+", " ", reaction)
        matched = next((needle for needle in counts if compact == re.sub(r"\s+", " ", needle)), None)
        if matched is None:
            output.append(line)
            continue
        if "!" not in line:
            raise RuntimeError(f"Selected reaction has no rate expression: {reaction}")
        left, expression = line.split("!", 1)
        payload = expression.rstrip("\r\n").strip()
        if not payload:
            raise RuntimeError(f"Selected reaction has an empty rate expression: {reaction}")
        output.append(f"{left}! ({payload})*{factor:.16g}d0\n")
        counts[matched] += 1
    missing = [needle for needle, count in counts.items() if count != 1]
    if missing:
        raise RuntimeError(f"Invalid reaction selector for {control}: {missing}")
    return "".join(output), list(counts)


def prepare_branch(base: Path, root: Path, control: str | None, factor: float) -> tuple[Path, str, dict]:
    label = "baseline" if control is None else f"{control}_{'down' if factor < 1 else 'up'}"
    build_dir = root / "builds" / label
    copy_build_inputs(base, build_dir)
    source = base / "kinet.inp"
    if control is None:
        text, selectors = source.read_text(encoding="utf-8"), []
    else:
        text, selectors = scale_group(source.read_text(encoding="utf-8"), control, factor)
    kinet = build_dir / "kinet.inp"
    kinet.write_text(text, encoding="utf-8", newline="\n")
    manifest = {
        "scenario": label, "control": control or "baseline", "factor": factor,
        "selectors": selectors, "source_kinet_sha256": sha256(source),
        "branch_kinet_sha256": sha256(kinet), "created_at": datetime.now().isoformat(timespec="seconds"),
    }
    (build_dir / "scenario_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    # Build once per distinct mechanism scenario. The runtime CSTR parameters
    # are supplied afterwards and do not alter the compiled reaction mechanism.
    params = {
        "en": "140", "tg": "300", "n2frac": "0.1", "tend": "1", "en_off": "0.1",
        "freq": "1000", "duty": "0.2", "cycles": "120", "ne_mode": "fix", "tau_res": "1e-2",
        "atol": "", "rtol": "", "en_list": "", "recompile": False, "pulse_enabled": False, "tag": "p9_build",
    }
    def safe_log(message) -> None:
        print(str(message).encode("ascii", errors="backslashreplace").decode("ascii"), flush=True)

    rc = gui.execute_pipeline(build_dir, params, do_build=True, do_run=False, log_cb=safe_log)
    if rc != 0 or not (build_dir / "main_pulse.exe").is_file():
        raise RuntimeError(f"Could not build scenario {label}")
    return build_dir, label, manifest


def copy_runtime(build_dir: Path, case_dir: Path) -> None:
    case_dir.mkdir(parents=True, exist_ok=False)
    for source in build_dir.iterdir():
        if source.is_file() and (
            source.name in {"main_pulse.exe", "main_pulse.F90", "kinet.inp", "bolsigdb.dat"}
            or source.suffix.lower() in {".dll", ".lib", ".dat"}
        ):
            shutil.copy2(source, case_dir / source.name)


def source_fluxes(rates: Path, first_cycle: int) -> dict[str, float]:
    totals = defaultdict(float)
    with rates.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            try:
                if int(row["cycle"]) < first_cycle:
                    continue
                rate, dt = float(row["rate_cm-3s-1"]), float(row["dt_s"])
            except (KeyError, TypeError, ValueError):
                continue
            if not (math.isfinite(rate) and math.isfinite(dt) and dt > 0):
                continue
            reaction = (row.get("reaction") or "").replace(" ", "")
            for label, pattern in PRIMARY_PATTERNS:
                if re.search(pattern, reaction):
                    totals[label] += rate * dt
                    break
    return dict(totals)


def terminal_metrics(outdir: Path, tag: str) -> dict[str, float]:
    first, last = steady_window_start(outdir / f"pulse_cycles_{tag}.csv", keep_cycles=3)
    channels, start, end = integrate_nh_channels(outdir / f"pulse_rates_{tag}.csv", first_cycle=first)
    shares = {str(row["channel"]): float(row["share_pct"]) for row in channels}
    fluxes = source_fluxes(outdir / f"pulse_rates_{tag}.csv", first)
    named = fluxes.get("Named-H2*", 0.0)
    secondary = max((value for name, value in fluxes.items() if name != "Named-H2*"), default=0.0)
    series = list(csv.DictReader((outdir / f"pulse_series_{tag}.csv").open(encoding="utf-8", newline="")))
    return {
        "cycle_first": first, "cycle_last": last, "terminal_duration_s": end - start,
        "NH3_final_cm-3": float(series[-1]["NH3"]),
        "named_H2star_pct": shares.get("Named-H2*", float("nan")),
        "named_H2star_flux_cm-3": named, "largest_secondary_flux_cm-3": secondary,
        "named_H2star_dominance_margin": named / secondary if secondary > 0 else float("inf"),
    }


def run_case(build_dir: Path, root: Path, scenario: str, manifest: dict, state: tuple[str, str, str, str]) -> dict:
    code, en, n2frac, title = state
    case_id = f"{scenario}__{code}"
    case_dir = root / "cases" / case_id
    copy_runtime(build_dir, case_dir)
    tag = f"p9_{scenario[:8]}_{code[:4]}"
    params = {
        "en": en, "tg": "300", "n2frac": n2frac, "tend": "1", "en_off": "0.1",
        "freq": "1000", "duty": "0.2", "cycles": "120", "ne_mode": "fix", "tau_res": "1e-2",
        "atol": "", "rtol": "", "en_list": "", "recompile": False, "pulse_enabled": False, "tag": tag,
    }
    def safe_log(message) -> None:
        print(str(message).encode("ascii", errors="backslashreplace").decode("ascii"), flush=True)

    rc = gui.execute_pipeline(case_dir, params, do_build=False, do_run=True, log_cb=safe_log)
    outdir = case_dir / gui.out_rel_for(params, "pulse")
    record_path = outdir / "params.json"
    record = json.loads(record_path.read_text(encoding="utf-8")) if record_path.is_file() else {}
    row = {
        "case_id": case_id, "scenario": scenario, "control": manifest["control"], "factor": manifest["factor"],
        "state": code, "state_title": title, "en_td": int(en), "n2_fraction": float(n2frac),
        "run_rc": rc, "run_status": record.get("status", "missing_record"),
        "failure_reason": record.get("failure_reason") or "", "outdir": str(outdir),
        "mechanism_sha256": manifest["branch_kinet_sha256"],
    }
    if rc == 0 and row["run_status"] == "success":
        row.update(terminal_metrics(outdir, tag))
    return row


def plot_results(table: pd.DataFrame, output: Path) -> None:
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "pdf.fonttype": 42})
    accepted = table.loc[table.run_status.eq("success")].copy()
    if len(accepted) != len(table):
        raise RuntimeError("Refusing to plot a scenario envelope containing non-P0-accepted cases")
    controls = ["baseline", *CONTROLS]
    labels = {"baseline": "Baseline", **{key: value["label"] for key, value in CONTROLS.items()}}
    colors = {"baseline": "#222222", "named_H2star_to_NH": "#1f77b4",
              "named_H2star_wall_loss": "#ff7f0e", "rydberg_to_NH": "#2ca02c"}
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.2), constrained_layout=True)
    for ax, metric, ylabel in zip(
        axes,
        ("named_H2star_pct", "named_H2star_dominance_margin"),
        ("Named H$_2$* direct-NH share (%)", "Named H$_2$* dominance margin"),
    ):
        for state_index, (state, group) in enumerate(accepted.groupby("state", sort=False)):
            x0 = state_index * (len(controls) + 1)
            for offset, control in enumerate(controls):
                selected = group.loc[group.control.eq(control)].sort_values("factor")
                if control == "baseline":
                    value = float(selected[metric].iloc[0])
                    ax.scatter(x0 + offset, value, color=colors[control], marker="D", s=27,
                               label="Baseline" if state_index == 0 else None, zorder=3)
                else:
                    ax.plot([x0 + offset - 0.16, x0 + offset + 0.16], selected[metric].to_numpy(float),
                            color=colors[control], lw=1.1)
                    ax.scatter([x0 + offset - 0.16, x0 + offset + 0.16], selected[metric].to_numpy(float),
                               color=colors[control], s=18, label=labels[control] if state_index == 0 else None)
        ax.set_ylabel(ylabel)
        ax.set_xlim(-0.5, 8.5)
        ax.set_xticks([1.5, 6.5], [
            "Productivity maximum\n140 Td, x$_{N_2}$ = 0.1",
            "Loss dominated\n240 Td, x$_{N_2}$ = 0.9",
        ])
        ax.tick_params(axis="x", labelsize=7, pad=5)
        ax.grid(axis="y", alpha=0.25)
    axes[1].set_yscale("log")
    axes[1].axhline(1.0, color="black", lw=0.8, ls="--")
    axes[0].legend(fontsize=6.4, frameon=False, loc="best")
    fig.text(0.5, -0.07,
             "Each colored pair is an independent 1/3 and 3× rate-group scenario; "
             "all cases pass P0. This is a deterministic scenario envelope, not statistical UQ.",
             ha="center", va="top", fontsize=7)
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / "Figure_18_H2star_scenario_envelope.png", dpi=600, bbox_inches="tight")
    fig.savefig(output / "Figure_18_H2star_scenario_envelope.pdf", bbox_inches="tight")
    plt.close(fig)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Hong H2* deterministic scenario envelope")
    parser.add_argument("--base", type=Path, default=BASE)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plot-existing", action="store_true",
                        help="regenerate Figure 18 from an existing complete scenario table")
    args = parser.parse_args(argv)
    base, root = args.base.resolve(), args.output_root.resolve()
    if args.plot_existing:
        table_path = root / "P9_scenario_cases.csv"
        if not table_path.is_file():
            raise SystemExit(f"Existing scenario table not found: {table_path}")
        plot_results(pd.read_csv(table_path), root / "figures")
        return 0
    if root.exists():
        raise SystemExit(f"Output root already exists; refusing to overwrite: {root}")
    root.mkdir(parents=True)
    branch_specs = [(None, 1.0), *[(control, factor) for control in CONTROLS for factor in FACTORS]]
    rows = []
    for number, (control, factor) in enumerate(branch_specs, start=1):
        build_dir, scenario, manifest = prepare_branch(base, root, control, factor)
        for state in STATES:
            print(f"=== [{number}/{len(branch_specs)}] {scenario} | {state[0]} ===", flush=True)
            rows.append(run_case(build_dir, root, scenario, manifest, state))
            pd.DataFrame(rows).to_csv(root / "P9_scenario_cases.partial.csv", index=False)
    table = pd.DataFrame(rows)
    table.to_csv(root / "P9_scenario_cases.csv", index=False)
    accepted = table.loc[table.run_status.eq("success")]
    report = [
        "# P9 H2* deterministic scenario envelope", "",
        "This is not statistical UQ. Each named reaction group is independently scaled by 1/3 and 3.",
        f"Requested cases: {len(table)}. P0-accepted cases: {len(accepted)}.",
    ]
    if len(accepted) == len(table):
        report += [
            f"Minimum named-H2* share: {accepted.named_H2star_pct.min():.4f}%.",
            f"Minimum named-H2* dominance margin: {accepted.named_H2star_dominance_margin.min():.6g}.",
        ]
        plot_results(table, root / "figures")
    (root / "P9_H2star_scenario_envelope.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    if len(accepted) != len(table):
        print(table.loc[~table.run_status.eq('success'), ['case_id', 'run_status', 'failure_reason']].to_string(index=False))
        return 1
    print(f"P0-accepted cases: {len(table)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
