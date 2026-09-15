#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Targeted CSTR kinetic-control atlas for manuscript Figure 11.

The six mechanism branches (three controls x symmetric perturbations) are
compiled once each.  Their copied pulse driver accepts an extra ``metrics``
argument that preserves the governing equations, terminal species CSV and P0
cycle CSV, but suppresses the 469-reaction ROP rows.  Thus the 180-point
productivity sensitivity atlas is auditable without generating unused
multi-gigabyte ROP archives.  The original driver is never modified.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

import zdp_gui as gui
import zdp_runtime as runtime
from hong_reaction_control import CONTROLS, copy_case_inputs, scale_reaction_group


TOOL_DIR = Path(__file__).resolve().parent
PROJECT = TOOL_DIR.parent
BASE = PROJECT / "Reproduction" / "2017Hong" / "方案2_脉冲能效地图" / "build_pulse"
DEFAULT_OUTPUT = PROJECT / "analysis" / "p7_control_switching_atlas_20260830"
N2_VALUES = (0.1, 0.3, 0.5, 0.7, 0.9)
EN_VALUES = (20, 60, 100, 140, 180, 240)
CONTROLS_ATLAS = ("NH_to_NH3_association", "NH2_hydrogenation", "proton_NH3_ionization")
CONTROL_CODE = {"NH_to_NH3_association": "nh", "NH2_hydrogenation": "nh2",
                "proton_NH3_ionization": "hion"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_log(message) -> None:
    """Avoid a legacy-codepage glyph from aborting a long batch."""
    print(str(message).encode("ascii", errors="backslashreplace").decode("ascii"), flush=True)


def patch_metrics_driver(path: Path) -> None:
    """Add a copied-driver-only ROP suppression flag with exact source guards."""
    text = path.read_text(encoding="utf-8")
    old_decl = "integer :: converged, stable_streak, min_cycles, stable_required, flow_enabled"
    new_decl = old_decl + ", write_rates_enabled"
    old_default = "  stable_required = 3\n"
    new_default = old_default + "  write_rates_enabled = 1\n"
    old_args = "  if (iargc() .ge. 12) then; call getarg(12, arg); read(arg, *) tau_res; end if\n"
    new_args = old_args + (
        "  if (iargc() .ge. 13) then\n"
        "    call getarg(13, arg)\n"
        "    if (trim(adjustl(arg)) .eq. 'metrics') write_rates_enabled = 0\n"
        "  end if\n")
    old_rate_block = """      if (jj .eq. 1) then
        if (iphase .eq. 1) then
          call write_rates(ur, tcur, dt1, k, 'on')
        else
          call write_rates(ur, tcur, dt1, k, 'off')
        end if
      else
        if (iphase .eq. 1) then
          call write_rates(ur, tcur, dti, k, 'on')
        else
          call write_rates(ur, tcur, dti, k, 'off')
        end if
      end if
"""
    new_rate_block = """      if (write_rates_enabled .eq. 1) then
        if (jj .eq. 1) then
          if (iphase .eq. 1) then
            call write_rates(ur, tcur, dt1, k, 'on')
          else
            call write_rates(ur, tcur, dt1, k, 'off')
          end if
        else
          if (iphase .eq. 1) then
            call write_rates(ur, tcur, dti, k, 'on')
          else
            call write_rates(ur, tcur, dti, k, 'off')
          end if
        end if
      end if
"""
    for old, new, name in ((old_decl, new_decl, "integer declaration"),
                           (old_default, new_default, "default flag"),
                           (old_args, new_args, "argument parser"),
                           (old_rate_block, new_rate_block, "rate output block")):
        if text.count(old) != 1:
            raise RuntimeError(f"Metrics-driver patch expected one {name}, found {text.count(old)}")
        text = text.replace(old, new, 1)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def build_branch(base: Path, root: Path, control: str, factor: float) -> dict:
    sign = "p" if factor > 1 else "m"
    branch_id = f"{CONTROL_CODE[control]}_{sign}"
    branch = root / "branches" / branch_id
    copy_case_inputs(base, branch)
    source_kinet = base / "kinet.inp"
    text, selected = scale_reaction_group(source_kinet.read_text(encoding="utf-8"), control, factor)
    kinet = branch / "kinet.inp"
    with kinet.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
    main = branch / "main_pulse.F90"
    patch_metrics_driver(main)
    manifest = {
        "branch_id": branch_id, "control": control, "control_label": CONTROLS[control]["label"],
        "factor": factor, "selectors": selected, "metrics_only": True,
        "source_kinet_sha256": sha256(source_kinet), "branch_kinet_sha256": sha256(kinet),
        "source_main_sha256": sha256(base / "main_pulse.F90"), "branch_main_sha256": sha256(main),
        "created_at": datetime.now().isoformat(timespec="seconds"),
    }
    (branch / "atlas_branch_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    params = {"en": "120", "tg": "300", "n2frac": "0.5", "tend": "1", "en_off": "0.1",
              "freq": "1000", "duty": "0.2", "cycles": "90", "ne_mode": "fix", "tau_res": "1e-2",
              "atol": "", "rtol": "", "en_list": "", "recompile": False, "pulse_enabled": False,
              "tag": "p7_build"}
    rc = gui.execute_pipeline(branch, params, do_build=True, do_run=False, log_cb=safe_log)
    if rc != 0:
        raise RuntimeError(f"Compilation failed for {branch_id}: rc={rc}")
    executable = runtime.find_executable(branch)
    if executable is None:
        raise RuntimeError(f"Compiled executable not found for {branch_id}")
    return {"branch_id": branch_id, "branch": branch, "control": control, "factor": factor,
            "executable": executable, "manifest": manifest}


def prepare_outdir(branch: Path, tag: str) -> Path:
    outdir = branch / "runs" / tag
    if outdir.exists():
        raise RuntimeError(f"Refusing to overwrite existing atlas point: {outdir}")
    outdir.mkdir(parents=True)
    bolsig = branch / "bolsigdb.dat"
    if bolsig.is_file():
        shutil.copy2(bolsig, outdir / bolsig.name)
    for data in branch.glob("*.DAT"):
        shutil.copy2(data, outdir / data.name)
    return outdir


def read_terminal_productivity(outdir: Path, tag: str) -> tuple[float, int, int, float, float]:
    cycles_path = outdir / f"pulse_cycles_{tag}.csv"
    series_path = outdir / f"pulse_series_{tag}.csv"
    with cycles_path.open(encoding="utf-8", newline="") as handle:
        cycles = list(csv.DictReader(handle))
    with series_path.open(encoding="utf-8", newline="") as handle:
        series = list(csv.DictReader(handle))
    if not cycles or not series:
        raise RuntimeError("metrics-only run did not produce terminal CSV records")
    last = cycles[-1]
    stable = int(float(last["stable_streak"]))
    if stable < 3:
        raise RuntimeError("metrics-only run did not reach a three-cycle P0 steady streak")
    tau = float(last["tau_res_s"])
    nh3 = float(series[-1]["NH3"])
    return nh3 / tau, int(last["cycle"]) - 2, int(last["cycle"]), float(last["rel_state_max"]), float(last["rel_dNH3_cycle"])


def run_point(branch_info: dict, en: int, n2: float, point_index: int) -> dict:
    sign = "p" if branch_info["factor"] > 1 else "m"
    tag = f"p7_{CONTROL_CODE[branch_info['control']]}{sign}_{point_index:02d}"
    outdir = prepare_outdir(branch_info["branch"], tag)
    args = [f"{n2:.1f}", "300", str(en), "0.1", "1000", "1.0", "90", "fix", tag,
            "1.0", "1e-4", "1e-2", "metrics"]
    env = runtime._runtime_env(runtime.resolve_gfortran_dir())
    rc = runtime.run_process([branch_info["executable"], *args], outdir, log_cb=safe_log, env=env,
                             log_path=outdir / "console.log")
    row = {"branch_id": branch_info["branch_id"], "control": branch_info["control"],
           "control_label": CONTROLS[branch_info["control"]]["label"], "role": CONTROLS[branch_info["control"]]["role"],
           "factor": branch_info["factor"], "en_td": en, "n2_fraction": n2, "tag": tag,
           "run_rc": rc, "outdir": str(outdir), "mechanism_sha256": branch_info["manifest"]["branch_kinet_sha256"]}
    if rc == 0:
        try:
            product, first, last, rel_state, rel_nh3 = read_terminal_productivity(outdir, tag)
            row.update({"run_status": "success", "NH3_productivity_cm-3s-1": product,
                        "cycle_first": first, "cycle_last": last, "rel_state_max": rel_state,
                        "rel_dNH3_cycle": rel_nh3})
        except Exception as exc:  # noqa: BLE001
            row.update({"run_status": "invalid_output", "failure_reason": str(exc)})
    else:
        row.update({"run_status": "failed", "failure_reason": "non-zero executable exit"})
    return row


def sensitivities(cases: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for control in CONTROLS_ATLAS:
        for en in EN_VALUES:
            for n2 in N2_VALUES:
                sub = cases.loc[(cases.control.eq(control)) & cases.en_td.eq(en) & np.isclose(cases.n2_fraction, n2)]
                if len(sub) != 2 or not sub.run_status.eq("success").all():
                    raise RuntimeError(f"Cannot form central sensitivity for {control}, {en} Td, N2={n2}")
                minus = sub.loc[sub.factor.lt(1)].iloc[0]; plus = sub.loc[sub.factor.gt(1)].iloc[0]
                y_minus, y_plus = float(minus["NH3_productivity_cm-3s-1"]), float(plus["NH3_productivity_cm-3s-1"])
                if not (math.isfinite(y_minus) and math.isfinite(y_plus) and y_minus > 0 and y_plus > 0):
                    raise RuntimeError("Non-positive productivity in local sensitivity")
                rows.append({"control": control, "control_label": CONTROLS[control]["label"],
                             "role": CONTROLS[control]["role"], "en_td": en, "n2_fraction": n2,
                             "sensitivity_productivity": (math.log(y_plus)-math.log(y_minus)) /
                             (math.log(float(plus.factor))-math.log(float(minus.factor))),
                             "productivity_minus_cm-3s-1": y_minus, "productivity_plus_cm-3s-1": y_plus})
    return pd.DataFrame(rows)


def write_report(root: Path, cases: pd.DataFrame, sensitivity: pd.DataFrame) -> None:
    grouped = []
    for control in CONTROLS_ATLAS:
        sub = sensitivity.loc[sensitivity.control.eq(control)]
        grouped.append((CONTROLS[control]["label"], float(sub.sensitivity_productivity.min()),
                        float(sub.sensitivity_productivity.max())))
    lines = ["# P7 kinetic-control switching atlas", "",
             "## Scope", "",
             "This targeted 5 x 6 (N2 fraction x E/N) atlas contains central local sensitivities of CSTR NH3 outlet productivity. ",
             "Each sensitivity uses 1/1.10 and 1.10 rate multipliers. The copied metric-only driver writes no time-resolved ROP rows, ",
             "but retains the identical kinetics, integration, terminal species output and P0 convergence diagnostics.", "",
             "## Acceptance", "", f"All {len(cases)} requested perturbed cases have status `success`.", "",
             "## Sensitivity ranges", ""]
    for label, low, high in grouped:
        lines.append(f"- `{label}`: S(productivity) from {low:.4f} to {high:.4f}.")
    lines.extend(["", "## Interpretation boundary", "",
                  "The atlas identifies local productivity control only for the three pre-specified groups. It does not constitute a global ",
                  "all-reaction sensitivity analysis or a rate-coefficient uncertainty quantification."])
    (root / "P7_control_switching_atlas.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Hong targeted kinetic-control switching atlas")
    parser.add_argument("--base", type=Path, default=BASE)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--smoke", action="store_true", help="compile one branch and run one ± pair at 140 Td, N2=0.1")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    base, root = args.base.resolve(), args.output_root.resolve()
    if not (base / "kinet.inp").is_file() or not (base / "main_pulse.F90").is_file():
        raise SystemExit("Validated Hong pulse/CSTR source inputs were not found")
    if root.exists():
        raise SystemExit(f"Output root already exists; refusing to overwrite: {root}")
    root.mkdir(parents=True)
    controls = ("NH_to_NH3_association",) if args.smoke else CONTROLS_ATLAS
    points = ((140, 0.1),) if args.smoke else tuple((en, n2) for en in EN_VALUES for n2 in N2_VALUES)
    branches = []
    for control in controls:
        for factor in (1.0 / 1.10, 1.10):
            print(f"\n=== Compiling branch {control}, factor={factor:.8g} ===")
            branches.append(build_branch(base, root, control, factor))
    results = []
    total = len(branches) * len(points)
    sequence = 0
    for branch in branches:
        for index, (en, n2) in enumerate(points, start=1):
            sequence += 1
            print(f"\n=== [{sequence}/{total}] {branch['branch_id']} | {en} Td | N2={n2:.1f} ===")
            results.append(run_point(branch, en, n2, index))
            pd.DataFrame(results).to_csv(root / "P7_control_cases.partial.csv", index=False)
            if results[-1]["run_status"] != "success":
                raise RuntimeError(f"Atlas run failed: {results[-1]}")
    cases = pd.DataFrame(results)
    cases.to_csv(root / "P7_control_cases.csv", index=False)
    if args.smoke:
        return 0
    sensitivity = sensitivities(cases)
    sensitivity.to_csv(root / "P7_control_sensitivities.csv", index=False)
    write_report(root, cases, sensitivity)
    print(sensitivity.groupby("control").sensitivity_productivity.agg(["min", "max"]).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
