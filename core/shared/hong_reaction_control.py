#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Local kinetic-control analysis for the validated Hong CSTR calculations.

Every control group is perturbed in an isolated copy of the validated kinetic
case.  A symmetric logarithmic finite difference, using factors 1.10 and
1/1.10, gives the local kinetic sensitivity S=d ln(y)/d ln(k).  This is a
numerical-control diagnostic rather than an uncertainty analysis: it does not
assign an uncertainty range to any elementary rate coefficient.

The two representative states are deliberately complementary: the global
NH3-productivity maximum (140 Td, x_N2=0.1) and the NH3-destruction dominated
state (240 Td, x_N2=0.9).  All accepted cases use the same CW, fixed-electron-
density, 0D CSTR/P0 boundary as the 108-point primary map.
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

import numpy as np
import pandas as pd

from _paths import LEGACY_HONG_ROOT
import zdp_gui as gui
from hong_reaction_hypergraph import parse_reaction


TOOL_DIR = Path(__file__).resolve().parent
PROJECT = LEGACY_HONG_ROOT
BASE = PROJECT / "Reproduction" / "2017Hong" / "方案2_脉冲能效地图" / "build_pulse"
DEFAULT_OUTPUT = PROJECT / "analysis" / "p6_local_reaction_control_20260830"

# The selectors identify a deliberately small, mechanistically interpretable
# panel.  They are not a blind all-reaction screen.  Each member of a group is
# multiplied by the same scalar in a copied kinet.inp.
CONTROLS: dict[str, dict] = {
    "named_H2star_to_NH": {
        "label": "Named electronic H2* → NH",
        "role": "primary NH source",
        "needles": tuple(f"N + H2({state}) => H + NH" for state in ("B3SIG", "B1SIG", "C3PI", "A3SIG")),
    },
    "rydberg_H2star_to_NH": {
        "label": "Rydberg-H2* → NH",
        "role": "secondary NH source",
        "needles": ("N + H2(RYDBERG_SUM) => H + NH",),
    },
    "excited_N_to_NH": {
        "label": "N(2D,2P) + H2 → NH",
        "role": "secondary NH source",
        "needles": ("N(2D) + H2 => H + NH", "N(2P) + H2 => H + NH"),
    },
    "termolecular_NH": {
        "label": "H + N + M → NH + M",
        "role": "associative NH source",
        "needles": ("H + N + @M => NH + @M",),
    },
    "NH2_hydrogenation": {
        "label": "H + NH2 + M → NH3 + M",
        "role": "NH3 formation",
        "needles": ("H + NH2 + @M => NH3 + @M",),
    },
    "NH_to_NH3_association": {
        "label": "NH + H2 + M → NH3 + M",
        "role": "NH3 formation",
        "needles": ("NH + H2 + @M => NH3 + @M",),
    },
    "proton_NH3_ionization": {
        "label": "H+ + NH3 → NH3+ + H",
        "role": "NH3 ionization loss",
        "needles": ("H^+ + NH3 => NH3^+ + H",),
    },
    "NH3plus_NH3_conversion": {
        "label": "NH3+ + NH3 → NH4+ + NH2",
        "role": "ion-mediated NH3 loss",
        "needles": ("NH3^+ + NH3 => NH4^+ + NH2",),
    },
}

REPRESENTATIVES = {
    "productivity_maximum": {"en": "140", "n2frac": "0.1", "title": r"NH$_3$ productivity maximum (140 Td, $x_{N_2}$=0.1)"},
    "loss_dominated": {"en": "240", "n2frac": "0.9", "title": r"Loss-dominated high-field state (240 Td, $x_{N_2}$=0.9)"},
}
REPRESENTATIVE_CODE = {"productivity_maximum": "opt", "loss_dominated": "high"}
CONTROL_CODE = {
    "baseline": "base", "named_H2star_to_NH": "named", "rydberg_H2star_to_NH": "ryd",
    "excited_N_to_NH": "exn", "termolecular_NH": "term", "NH2_hydrogenation": "nh2",
    "NH_to_NH3_association": "nh", "proton_NH3_ionization": "hion",
    "NH3plus_NH3_conversion": "qconv",
}

NH_SOURCE_PATTERNS = {
    "named_H2star": re.compile(r"^N\+H2\((B3SIG|B1SIG|C3PI|A3SIG)\)=>H\+NH$"),
    "rydberg": re.compile(r"^N\+H2\(RYDBERG_SUM\)=>H\+NH$"),
    "excited_N": re.compile(r"^N\((2D|2P)\)\+H2=>H\+NH$"),
    "association": re.compile(r"^H\+N\+(N2|H2)=>NH\+\1$"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_case_inputs(base: Path, target: Path) -> None:
    """Copy only compiler/runtime inputs; never carry previous outputs forward."""
    target.mkdir(parents=True, exist_ok=False)
    keep_names = {"preprocessor.exe", "bolsigdb.dat", "main_pulse.F90", "dvode_f90_m.F90"}
    keep_suffixes = {".dll", ".lib", ".dat", ".DAT"}
    for source in base.iterdir():
        if source.is_file() and (source.name in keep_names or source.suffix in keep_suffixes):
            shutil.copy2(source, target / source.name)


def normalize_reaction(text: str) -> str:
    return re.sub(r"\s+", " ", text.strip())


def scale_reaction_group(kinet_text: str, control: str, factor: float) -> tuple[str, list[str]]:
    """Multiply the selected elementary-rate expressions in the input comments."""
    spec = CONTROLS[control]
    pending = {normalize_reaction(needle): 0 for needle in spec["needles"]}
    output: list[str] = []
    for line in kinet_text.splitlines(keepends=True):
        reaction = line.split("!", 1)[0].strip()
        key = normalize_reaction(reaction)
        if key not in pending:
            output.append(line)
            continue
        if "!" not in line:
            raise RuntimeError(f"Selected reaction has no rate expression: {reaction}")
        left, expression = line.split("!", 1)
        payload, newline = expression.rstrip("\r\n"), "\n" if line.endswith("\n") else ""
        if not payload.strip():
            raise RuntimeError(f"Selected reaction has an empty rate expression: {reaction}")
        output.append(f"{left}! ({payload.strip()})*{factor:.16g}d0{newline}")
        pending[key] += 1
    missing = [reaction for reaction, count in pending.items() if count != 1]
    if missing:
        raise RuntimeError(f"Expected exactly one occurrence for {control}; got invalid matches: {missing}")
    return "".join(output), list(pending)


def coefficient(items: tuple[str, ...], target: str) -> int:
    value = 0
    for item in items:
        matched = re.match(r"^(\d+)?(.+)$", item.strip())
        amount, species = int(matched.group(1) or "1"), matched.group(2)
        if species == target:
            value += amount
    return value


def terminal_window(cycles_path: Path, keep: int = 3) -> tuple[int, int, float]:
    with cycles_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) < keep:
        raise RuntimeError("P0 cycle record is shorter than the requested terminal window")
    last = rows[-1]
    if int(float(last["stable_streak"])) < keep:
        raise RuntimeError("Run did not supply three terminal P0-steady cycles")
    start = int(last["cycle"]) - keep + 1
    # The pulse driver records cycle-end times (not t_start_s).  The preceding
    # record is therefore the left boundary of the final `keep` cycles.
    duration = float(last["t_end_s"]) - float(rows[-keep - 1]["t_end_s"])
    if duration <= 0:
        raise RuntimeError("Invalid terminal P0-window duration")
    return start, int(last["cycle"]), duration


def terminal_metrics(outdir: Path, tag: str) -> dict[str, float]:
    """Integrate terminal ROP to obtain NH3 formation/destruction and NH sources."""
    cycle_first, cycle_last, duration = terminal_window(outdir / f"pulse_cycles_{tag}.csv")
    formation = loss = 0.0
    sources = defaultdict(float)
    with (outdir / f"pulse_rates_{tag}.csv").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                if int(row["cycle"]) < cycle_first:
                    continue
                rate, dt = float(row["rate_cm-3s-1"]), float(row["dt_s"])
            except (KeyError, TypeError, ValueError):
                continue
            if not (math.isfinite(rate) and math.isfinite(dt) and dt > 0):
                continue
            reaction = row.get("reaction", "")
            equation = parse_reaction(reaction)
            if equation is not None:
                reactants, products = equation
                contribution = (coefficient(products, "NH3") - coefficient(reactants, "NH3")) * rate * dt
                if contribution > 0:
                    formation += contribution
                elif contribution < 0:
                    loss -= contribution
            compact = re.sub(r"\s+", "", reaction.split(":", 1)[-1])
            for family, pattern in NH_SOURCE_PATTERNS.items():
                if pattern.search(compact):
                    sources[family] += rate * dt
                    break
    cycles = list(csv.DictReader((outdir / f"pulse_cycles_{tag}.csv").open(encoding="utf-8", newline="")))
    series = list(csv.DictReader((outdir / f"pulse_series_{tag}.csv").open(encoding="utf-8", newline="")))
    tau = float(cycles[-1]["tau_res_s"])
    product = float(series[-1]["NH3"]) / tau
    source_total = sum(sources.values())
    return {
        "cycle_first": cycle_first, "cycle_last": cycle_last, "terminal_duration_s": duration,
        "NH3_productivity_cm-3s-1": product,
        "NH3_formation_cm-3s-1": formation / duration,
        "NH3_loss_cm-3s-1": loss / duration,
        "NH3_loss_to_formation": loss / formation if formation > 0 else float("nan"),
        "named_H2star_source_share": sources["named_H2star"] / source_total if source_total > 0 else float("nan"),
    }


def run_case(base: Path, root: Path, representative: str, control: str, factor: float) -> dict:
    factor_label = "baseline" if control == "baseline" else ("plus" if factor > 1 else "minus")
    case_id = f"{representative}__{control}__{factor_label}"
    case_dir = root / "cases" / case_id
    copy_case_inputs(base, case_dir)
    source_kinet = base / "kinet.inp"
    source_text = source_kinet.read_text(encoding="utf-8")
    if control == "baseline":
        kinetic_text, matched = source_text, []
    else:
        kinetic_text, matched = scale_reaction_group(source_text, control, factor)
    kinet = case_dir / "kinet.inp"
    # Preserve LF line endings so an unmodified baseline also preserves the
    # source-file byte hash across Windows hosts.
    with kinet.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(kinetic_text)
    rep = REPRESENTATIVES[representative]
    # main_pulse.F90 stores output names in fixed-length character variables.
    # Keep tags short enough for every derived `pulse_*_<tag>.csv` filename.
    factor_code = "base" if control == "baseline" else ("p" if factor > 1 else "m")
    tag = f"p6_{REPRESENTATIVE_CODE[representative]}_{CONTROL_CODE[control]}_{factor_code}"
    params = {
        "en": rep["en"], "tg": "300", "n2frac": rep["n2frac"], "tend": "1",
        "en_off": "0.1", "freq": "1000", "duty": "0.2", "cycles": "90", "ne_mode": "fix",
        "tau_res": "1e-2", "atol": "", "rtol": "", "en_list": "", "recompile": False,
        "pulse_enabled": False, "tag": tag,
    }
    manifest = {
        "case_id": case_id, "representative": representative, "representative_title": rep["title"],
        "control": control, "factor": factor, "selectors": matched,
        "diagnostic_only": True, "source_kinet_sha256": sha256(source_kinet),
        "branch_kinet_sha256": sha256(kinet), "created_at": datetime.now().isoformat(timespec="seconds"),
    }
    (case_dir / "diagnostic_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    # Compiler/preprocessor output can contain legacy Windows encodings.  Keep
    # it observable without allowing an undecodable glyph to abort a solver run.
    def safe_log(message) -> None:
        print(str(message).encode("ascii", errors="backslashreplace").decode("ascii"), flush=True)

    rc = gui.execute_pipeline(case_dir, params, do_build=True, do_run=True, log_cb=safe_log)
    outdir = case_dir / gui.out_rel_for(params, "pulse")
    record_path = outdir / "params.json"
    record = json.loads(record_path.read_text(encoding="utf-8")) if record_path.is_file() else {}
    row = {"case_id": case_id, "representative": representative, "control": control, "factor": factor,
           "run_rc": rc, "run_status": record.get("status", "missing_record"),
           "failure_reason": record.get("failure_reason") or "", "outdir": str(outdir),
           "mechanism_sha256": manifest["branch_kinet_sha256"]}
    if rc == 0 and row["run_status"] == "success":
        row.update(terminal_metrics(outdir, tag))
    return row


def sensitivity_table(cases: pd.DataFrame) -> pd.DataFrame:
    outputs = ("NH3_productivity_cm-3s-1", "NH3_loss_to_formation", "named_H2star_source_share")
    records = []
    for representative in REPRESENTATIVES:
        subset = cases.loc[cases.representative.eq(representative)]
        baseline = subset.loc[subset.control.eq("baseline")].iloc[0]
        for control, spec in CONTROLS.items():
            lower = subset.loc[(subset.control.eq(control)) & (subset.factor < 1)].iloc[0]
            upper = subset.loc[(subset.control.eq(control)) & (subset.factor > 1)].iloc[0]
            for output in outputs:
                y_minus, y_plus = float(lower[output]), float(upper[output])
                value = (math.log(y_plus) - math.log(y_minus)) / (math.log(float(upper.factor)) - math.log(float(lower.factor)))
                records.append({
                    "representative": representative, "representative_title": REPRESENTATIVES[representative]["title"],
                    "control": control, "control_label": spec["label"], "role": spec["role"],
                    "output": output, "sensitivity": value,
                    "baseline_value": float(baseline[output]), "minus_value": y_minus, "plus_value": y_plus,
                })
    return pd.DataFrame(records)


def write_report(root: Path, cases: pd.DataFrame, sensitivity: pd.DataFrame) -> None:
    lines = [
        "# P6 local reaction-control analysis", "",
        "## Calculation boundary", "",
        "Symmetric local sensitivities were computed with factors 1/1.10 and 1.10 in isolated copies of `kinet.inp`. ",
        "The responses are terminal P0-window CSTR quantities in the established CW, fixed-electron-density boundary ",
        "(Tg = 300 K, tau = 10 ms).  This is a numerical control analysis, not a kinetic uncertainty quantification.", "",
        "## Accepted-run check", "",
        f"All {len(cases)} requested cases have status `success`.", "",
        "## Interpretation rule", "",
        "For a response y and a control-group multiplier k, `S = d ln(y) / d ln(k)`. Positive S means that increasing ",
        "the group increases y; negative S means suppression. The NH3 loss/formation response must be read as a turnover ",
        "balance, rather than as an energy-efficiency measure.", "",
        "## Largest local controls", "",
    ]
    for rep in REPRESENTATIVES:
        lines.extend([f"### {REPRESENTATIVES[rep]['title']}", ""])
        part = sensitivity.loc[(sensitivity.representative.eq(rep)) & sensitivity.output.eq("NH3_productivity_cm-3s-1")]
        for _, row in part.reindex(part.sensitivity.abs().sort_values(ascending=False).index).head(5).iterrows():
            lines.append(f"- `{row.control_label}`: S(productivity) = {row.sensitivity:.3f}.")
        lines.append("")
    (root / "P6_local_reaction_control.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Hong CSTR local reaction-control diagnostics")
    parser.add_argument("--base", type=Path, default=BASE)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--smoke", action="store_true", help="run baseline plus one named-H2* pair at the productivity maximum")
    parser.add_argument("--summarize-existing", action="store_true",
                        help="recompute the sensitivity CSV and report from an existing complete case table")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    base, root = args.base.resolve(), args.output_root.resolve()
    if args.summarize_existing:
        source = root / "P6_control_case_results.csv"
        if not source.is_file():
            raise SystemExit(f"Existing complete case table not found: {source}")
        cases = pd.read_csv(source)
        if len(cases) != 34 or not cases.run_status.eq("success").all():
            raise SystemExit("Existing table is not a complete 34-case accepted control panel")
        sensitivity = sensitivity_table(cases)
        sensitivity.to_csv(root / "P6_local_sensitivities.csv", index=False)
        write_report(root, cases, sensitivity)
        return 0
    if not (base / "kinet.inp").is_file():
        raise SystemExit(f"Base mechanism not found: {base / 'kinet.inp'}")
    if root.exists():
        raise SystemExit(f"Output root already exists; refusing to overwrite: {root}")
    root.mkdir(parents=True)
    jobs = [(representative, "baseline", 1.0) for representative in REPRESENTATIVES]
    controls = ("named_H2star_to_NH",) if args.smoke else tuple(CONTROLS)
    reps = ("productivity_maximum",) if args.smoke else tuple(REPRESENTATIVES)
    if args.smoke:
        jobs = [("productivity_maximum", "baseline", 1.0)]
    jobs += [(rep, control, factor) for rep in reps for control in controls for factor in (1.0 / 1.10, 1.10)]
    results = []
    for number, (representative, control, factor) in enumerate(jobs, start=1):
        print(f"\n=== [{number}/{len(jobs)}] {representative} | {control} | factor={factor:.8g} ===")
        results.append(run_case(base, root, representative, control, factor))
        pd.DataFrame(results).to_csv(root / "P6_control_case_results.partial.csv", index=False)
    cases = pd.DataFrame(results)
    cases.to_csv(root / "P6_control_case_results.csv", index=False)
    failures = cases.loc[~cases.run_status.eq("success")]
    if not failures.empty:
        print(failures[["case_id", "run_rc", "run_status", "failure_reason"]].to_string(index=False))
        return 1
    if args.smoke:
        return 0
    sensitivity = sensitivity_table(cases)
    sensitivity.to_csv(root / "P6_local_sensitivities.csv", index=False)
    write_report(root, cases, sensitivity)
    print(sensitivity.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
