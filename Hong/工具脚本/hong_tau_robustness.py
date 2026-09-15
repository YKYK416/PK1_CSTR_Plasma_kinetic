#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Residence-time robustness screen for the accepted Hong CW-CSTR mechanism.

This script does not change the validated baseline input.  It evaluates a
small, mechanistically contrasting panel of CSTR states across residence time,
and admits a case only when the existing terminal-window P0 gate succeeds.
It reports direct-NH family shares and a named-H2* dominance margin alongside
NH3 productivity and source–sink turnover.  The screen is deliberately not a
full reactor-design response surface or kinetic-uncertainty quantification.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

import zdp_gui as gui
from hong_nh_diagnostics import PRIMARY_PATTERNS, integrate_nh_channels, steady_window_start
from hong_reaction_control import coefficient
from hong_reaction_hypergraph import parse_reaction


TOOL_DIR = Path(__file__).resolve().parent
PROJECT = TOOL_DIR.parent
BASE = PROJECT / "Reproduction" / "2017Hong" / "方案2_脉冲能效地图" / "build_pulse"
DEFAULT_OUTPUT = PROJECT / "analysis" / "p8_tau_robustness_20260831"

STATES = (
    ("low_field", "60", "0.1", "Low field: 60 Td, x_N2 = 0.1"),
    ("productive", "140", "0.1", "Productivity maximum: 140 Td, x_N2 = 0.1"),
    ("mid_composition", "140", "0.5", "Intermediate composition: 140 Td, x_N2 = 0.5"),
    ("loss_dominated", "240", "0.9", "Loss dominated: 240 Td, x_N2 = 0.9"),
)
TAU_MS = (0.5, 1.0, 3.0, 10.0, 30.0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_case_inputs(base: Path, target: Path) -> None:
    target.mkdir(parents=True, exist_ok=False)
    # Residence time is passed at run time; every case uses the unmodified
    # baseline mechanism.  Reuse the already validated baseline executable
    # rather than recompiling identical source for each tau value.
    keep_names = {"main_pulse.exe", "main_pulse.F90", "bolsigdb.dat"}
    keep_suffixes = {".dll", ".lib", ".dat", ".DAT"}
    for source in base.iterdir():
        if source.is_file() and (source.name in keep_names or source.suffix in keep_suffixes):
            shutil.copy2(source, target / source.name)


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
    first_cycle, last_cycle = steady_window_start(outdir / f"pulse_cycles_{tag}.csv", keep_cycles=3)
    channels, start_s, end_s = integrate_nh_channels(
        outdir / f"pulse_rates_{tag}.csv", first_cycle=first_cycle
    )
    shares = {str(row["channel"]): float(row["share_pct"]) for row in channels}
    fluxes = source_fluxes(outdir / f"pulse_rates_{tag}.csv", first_cycle)
    named = fluxes.get("Named-H2*", 0.0)
    secondary = max((value for key, value in fluxes.items() if key != "Named-H2*"), default=0.0)

    formation = 0.0
    destruction = 0.0
    with (outdir / f"pulse_rates_{tag}.csv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            try:
                if int(row["cycle"]) < first_cycle:
                    continue
                rate, dt = float(row["rate_cm-3s-1"]), float(row["dt_s"])
            except (KeyError, TypeError, ValueError):
                continue
            if not (math.isfinite(rate) and math.isfinite(dt) and dt > 0):
                continue
            equation = parse_reaction(row.get("reaction") or "")
            if equation is None:
                continue
            reactants, products = equation
            contribution = (coefficient(products, "NH3") - coefficient(reactants, "NH3")) * rate * dt
            if contribution > 0:
                formation += contribution
            elif contribution < 0:
                destruction -= contribution

    cycles = list(csv.DictReader((outdir / f"pulse_cycles_{tag}.csv").open(encoding="utf-8", newline="")))
    series = list(csv.DictReader((outdir / f"pulse_series_{tag}.csv").open(encoding="utf-8", newline="")))
    duration = end_s - start_s
    tau_s = float(cycles[-1]["tau_res_s"])
    return {
        "cycle_first": first_cycle,
        "cycle_last": last_cycle,
        "terminal_duration_s": duration,
        "NH3_final_cm-3": float(series[-1]["NH3"]),
        "NH3_productivity_cm-3s-1": float(series[-1]["NH3"]) / tau_s,
        "NH3_formation_cm-3s-1": formation / duration,
        "NH3_loss_cm-3s-1": destruction / duration,
        "NH3_loss_to_formation": destruction / formation if formation > 0 else float("nan"),
        "named_H2star_pct": shares.get("Named-H2*", float("nan")),
        "named_H2star_flux_cm-3": named,
        "largest_secondary_flux_cm-3": secondary,
        "named_H2star_dominance_margin": named / secondary if secondary > 0 else float("inf"),
    }


def run_case(base: Path, root: Path, state: tuple[str, str, str, str], tau_ms: float) -> dict:
    code, en, n2frac, title = state
    tau_code = f"{tau_ms:g}".replace(".", "p")
    case_id = f"{code}_tau{tau_code}ms"
    case_dir = root / "cases" / case_id
    copy_case_inputs(base, case_dir)
    source_kinet = base / "kinet.inp"
    kinet = case_dir / "kinet.inp"
    shutil.copy2(source_kinet, kinet)
    tag = f"p8_{code[:4]}_{tau_code}"
    # One CW analysis window is 1 ms at the retained 1 kHz bookkeeping
    # frequency.  Require at least 12 residence times while retaining the
    # established 90-window baseline for short residence times.  P0 remains
    # the acceptance gate; this only prevents long-tau cases from being
    # rejected before sufficient CSTR flushing is possible.
    requested_cycles = max(90, int(math.ceil(12.0 * tau_ms)))
    params = {
        "en": en, "tg": "300", "n2frac": n2frac, "tend": "1", "en_off": "0.1",
        "freq": "1000", "duty": "0.2", "cycles": str(requested_cycles), "ne_mode": "fix",
        "tau_res": f"{tau_ms * 1e-3:.8g}", "atol": "", "rtol": "", "en_list": "",
        "recompile": False, "pulse_enabled": False, "tag": tag,
    }
    manifest = {
        "case_id": case_id, "state": {"code": code, "title": title, "en_td": en, "n2_fraction": n2frac},
        "tau_ms": tau_ms, "requested_cycles": requested_cycles,
        "calculation_boundary": "CW, fixed electron density, 0D CSTR, 300 K",
        "source_kinet_sha256": sha256(source_kinet), "branch_kinet_sha256": sha256(kinet),
        "created_at": datetime.now().isoformat(timespec="seconds"),
    }
    (case_dir / "case_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    def log(message) -> None:
        print(str(message).encode("ascii", errors="backslashreplace").decode("ascii"), flush=True)

    rc = gui.execute_pipeline(case_dir, params, do_build=False, do_run=True, log_cb=log)
    outdir = case_dir / gui.out_rel_for(params, "pulse")
    record_path = outdir / "params.json"
    record = json.loads(record_path.read_text(encoding="utf-8")) if record_path.is_file() else {}
    row = {
        "case_id": case_id, "state": code, "state_title": title, "en_td": int(en),
        "n2_fraction": float(n2frac), "tau_ms": tau_ms, "run_rc": rc,
        "run_status": record.get("status", "missing_record"),
        "failure_reason": record.get("failure_reason") or "", "outdir": str(outdir),
        "mechanism_sha256": manifest["branch_kinet_sha256"],
    }
    if rc == 0 and row["run_status"] == "success":
        row.update(terminal_metrics(outdir, tag))
    return row


def plot_results(table: pd.DataFrame, output: Path) -> None:
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "pdf.fonttype": 42})
    accepted = table.loc[table.run_status.eq("success")].copy()
    if accepted.empty:
        raise RuntimeError("No P0-accepted residence-time cases are available for plotting")
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.1), constrained_layout=True)
    metrics = (
        ("NH3_productivity_cm-3s-1", r"NH$_3$ productivity (cm$^{-3}$ s$^{-1}$)", True),
        ("named_H2star_pct", "Named H$_2$* direct-NH share (%)", False),
        ("named_H2star_dominance_margin", "Named H$_2$* dominance margin", True),
    )
    for axis, (metric, ylabel, logarithmic) in zip(axes, metrics):
        for state, group in accepted.groupby("state", sort=False):
            group = group.sort_values("tau_ms")
            axis.plot(group.tau_ms, group[metric], marker="o", lw=1.35, ms=3.8,
                      label=group.state_title.iloc[0])
        axis.set_xscale("log")
        if logarithmic:
            axis.set_yscale("log")
        axis.set_xlabel(r"Residence time, $\tau$ (ms)")
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.25, which="both")
    axes[0].legend(fontsize=6.5, frameon=False, loc="best")
    axes[2].axhline(1.0, color="black", lw=0.8, ls="--")
    fig.text(0.5, -0.06,
             "CW CSTR, 300 K, imposed n$_e$ = 1.17×10$^8$ cm$^{-3}$; "
             "only cases passing the terminal-window P0 gate are shown.",
             ha="center", va="top", fontsize=7)
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / "Figure_17_residence_time_robustness.png", dpi=600, bbox_inches="tight")
    fig.savefig(output / "Figure_17_residence_time_robustness.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(root: Path, table: pd.DataFrame) -> None:
    accepted = table.loc[table.run_status.eq("success")].copy()
    lines = [
        "# P8 residence-time robustness screen", "",
        "## Scope", "",
        "This is a targeted continuous-wave CSTR transport-robustness screen, not a full response surface. ",
        "All reported values use terminal-window integrated ROP and only P0-accepted cases.", "",
        "## Acceptance", "",
        f"Requested cases: {len(table)}. P0-accepted cases: {len(accepted)}.", "",
    ]
    if not accepted.empty:
        lines += [
            "## Direct-NH hierarchy", "",
            f"Minimum named-H2* direct-NH share: {accepted.named_H2star_pct.min():.4f}%.",
            f"Minimum named-H2* dominance margin: {accepted.named_H2star_dominance_margin.min():.6g}.",
            "",
            "## Interpretation guardrail", "",
            "A positive margin in this targeted panel supports transport robustness inside the fixed mechanism. ",
            "It neither substitutes for a global kinetic-uncertainty analysis nor extends the conclusion to pulsed, ",
            "surface-coupled, or self-consistent plasma conditions.",
        ]
    (root / "P8_residence_time_robustness.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Hong CW-CSTR residence-time robustness screen")
    parser.add_argument("--base", type=Path, default=BASE)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--smoke", action="store_true", help="run the productive state at 1 and 10 ms only")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    base, root = args.base.resolve(), args.output_root.resolve()
    if not (base / "kinet.inp").is_file():
        raise SystemExit(f"Base mechanism not found: {base / 'kinet.inp'}")
    if root.exists():
        raise SystemExit(f"Output root already exists; refusing to overwrite: {root}")
    root.mkdir(parents=True)
    states = (STATES[1],) if args.smoke else STATES
    taus = (1.0, 10.0) if args.smoke else TAU_MS
    jobs = [(state, tau) for state in states for tau in taus]
    rows = []
    for index, (state, tau) in enumerate(jobs, start=1):
        print(f"=== [{index}/{len(jobs)}] {state[0]} | tau={tau:g} ms ===", flush=True)
        rows.append(run_case(base, root, state, tau))
        pd.DataFrame(rows).to_csv(root / "P8_tau_cases.partial.csv", index=False)
    table = pd.DataFrame(rows)
    table.to_csv(root / "P8_tau_cases.csv", index=False)
    write_report(root, table)
    failures = table.loc[~table.run_status.eq("success")]
    if not failures.empty:
        print(failures[["case_id", "run_rc", "run_status", "failure_reason"]].to_string(index=False))
        return 1
    plot_results(table, root / "figures")
    print(f"P0-accepted cases: {len(table)}; figure: {root / 'figures' / 'Figure_17_residence_time_robustness.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
