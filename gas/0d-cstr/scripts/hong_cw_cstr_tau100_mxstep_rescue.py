#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rescue the tau=100 ms CW-CSTR cases that exhausted DVODE MXSTEP.

The ordinary final-uniform protocol advances from 10 to 20 s in one outer
10 s call.  For the affected low-N2 cases this lets DVODE consume MXSTEP in a
single call even though the NH3 trajectory is already nearly stationary.
This script leaves all prior results untouched and retries only those recorded
MXSTEP failures with 1 s outer calls after the existing 0.2--10 s transition.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Any


TOOL_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOL_DIR))
import hong_cw_longtime_scan as cw  # noqa: E402
import hong_cw_cstr_finalize_uniform_cases as uniform  # noqa: E402


SCAN_ROOT = Path(r"F:\Codex\PK1\GasReaction\CSTR\cw_gasphase_cstr_tau_1to100ms_20260907")
UNIFORM_ROOT = SCAN_ROOT / "final_uniform_validation_adaptive_r4_20260912"
DEFAULT_OUTPUT_ROOT = SCAN_ROOT / "tau100_mxstep_rescue_20260914"
MXSTEP = 500000
RTOL = 1.0e-4
OUTER_DT_S = 1.0
T_END_S = 20.0
REPORT_DT_S = 1.0
STARTUP_DT_S = 1.0e-3
STARTUP_DURATION_S = 2.0e-1
TRANSITION_DT_S = 1.0e-1
TRANSITION_DURATION_S = 10.0
ATOL_SEQUENCE = (1.0e2, 1.0e6, 1.0e8)


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    fields = sorted({key for record in records for key in record})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def load_mxstep_targets() -> list[dict[str, Any]]:
    summary_path = UNIFORM_ROOT / "uniform_validation_summary.csv"
    if not summary_path.is_file():
        raise RuntimeError(f"Missing uniform-validation summary: {summary_path}")
    targets = [
        dict(row) for row in read_csv(summary_path)
        if row.get("final_status") == "NUMERICAL_FAILURE_UNIFORM_VALIDATION"
        and row.get("failure_reason") == "process return code 2"
        and abs(float(row["tau_res_s"]) - 0.1) < 1.0e-12
    ]
    expected = {
        "T350_EN60_N20.1_TAU100ms",
        "T375_EN60_N20.1_TAU100ms",
        "T400_EN60_N20.1_TAU100ms",
    }
    found = {str(item["case_id"]) for item in targets}
    if found != expected:
        raise RuntimeError(f"Expected exactly the known three MXSTEP targets; found {sorted(found)}")
    return sorted(targets, key=lambda item: float(item["tg_K"]))


def policy_name(atol: float) -> str:
    return f"segmented_dt1s_rtol1e-4_atol{atol:.0e}".replace("+", "")


def run_target(work: Path, root: Path, target: dict[str, Any], index: int) -> dict[str, Any]:
    result: dict[str, Any] = {
        "case_id": target["case_id"],
        "original_case_id": target["original_case_id"],
        "tg_K": float(target["tg_K"]),
        "en_Td": float(target["en_Td"]),
        "n2_fraction": float(target["n2_fraction"]),
        "tau_res_s": float(target["tau_res_s"]),
        "reference_NH3_cm3": float(target["reference_NH3_cm3"]),
        "rescue_index": index,
        "outer_dt_s": OUTER_DT_S,
        "t_end_s": T_END_S,
        "rtol": RTOL,
        "started_at": now(),
    }
    for attempt, atol in enumerate(ATOL_SEQUENCE, 1):
        group = f"{policy_name(atol)}/attempt_{attempt:02d}"
        record = cw.run_case(
            work=work, root=root, tg=float(target["tg_K"]), en=float(target["en_Td"]),
            n2=float(target["n2_fraction"]), dt=OUTER_DT_S, t_end=T_END_S,
            report_dt=REPORT_DT_S, group=group, tau_res_s=float(target["tau_res_s"]),
            atol=atol, rtol=RTOL, startup_dt_s=STARTUP_DT_S,
            startup_duration_s=STARTUP_DURATION_S, transition_dt_s=TRANSITION_DT_S,
            transition_duration_s=TRANSITION_DURATION_S, timeout_s=21600.0,
            idle_timeout_s=300.0, watchdog_poll_s=5.0,
        )
        result.update(record)
        result["attempt"] = attempt
        result["atol"] = atol
        if record.get("returncode") != 0 or bool(record.get("timed_out")):
            result["last_failure_reason"] = str(record.get("failure_reason", "unknown"))
            continue
        result.update(uniform.terminal_block_test(root, record))
        final_nh3 = float(record["NH3_cm-3"])
        reference = float(target["reference_NH3_cm3"])
        relative_difference = abs(final_nh3 - reference) / max(abs(reference), 1.0)
        result["NH3_rel_diff_to_reference"] = relative_difference
        result["reference_rel_diff_limit"] = uniform.REFERENCE_REL_DIFF_LIMIT
        if bool(result["stationary_by_uniform_error_aware_test"]) and relative_difference <= uniform.REFERENCE_REL_DIFF_LIMIT:
            result["final_status"] = "ACCEPTED_MXSTEP_SEGMENTED_UNIFORM"
        elif not bool(result["stationary_by_uniform_error_aware_test"]):
            result["final_status"] = "SEGMENTED_RUN_REQUIRES_LONGER_TIME"
        else:
            result["final_status"] = "SEGMENTED_RUN_REFERENCE_DISAGREEMENT"
        result["ended_at"] = now()
        return result
    result["final_status"] = "MXSTEP_RESCUE_FAILED_ALL_ATOL"
    result["ended_at"] = now()
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Rescue the three tau=100 ms MXSTEP-limited CSTR cases")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.workers <= 0:
        parser.error("workers must be positive")
    targets = load_mxstep_targets()
    if args.dry_run:
        print(json.dumps({"targets": targets, "policies": list(ATOL_SEQUENCE)}, ensure_ascii=False, indent=2))
        return 0
    root = args.output_root.resolve()
    existing: list[dict[str, Any]] = []
    summary_path = root / "mxstep_rescue_summary.csv"
    if root.exists():
        if not args.resume:
            raise RuntimeError(f"Refusing to overwrite existing rescue root: {root}; use --resume")
        if summary_path.is_file():
            existing = [dict(row) for row in read_csv(summary_path)]
    else:
        root.mkdir(parents=True)
    accepted_ids = {
        str(row["case_id"]) for row in existing
        if row.get("final_status") == "ACCEPTED_MXSTEP_SEGMENTED_UNIFORM"
    }
    pending = [(index, target) for index, target in enumerate(targets, 1) if str(target["case_id"]) not in accepted_ids]
    (root / "PROTOCOL.json").write_text(json.dumps({
        "created_at": now(), "source_uniform_root": str(UNIFORM_ROOT),
        "purpose": "rescue low-N2 tau=100 ms MXSTEP failures by reducing outer integration segment to 1 s",
        "targets": [item["case_id"] for item in targets], "outer_dt_s": OUTER_DT_S,
        "t_end_s": T_END_S, "rtol": RTOL, "atol_sequence": ATOL_SEQUENCE,
        "mxstep": MXSTEP, "original_results_modified": False,
        "acceptance": "same terminal NH3 stationarity and reference-relative-difference criteria as final uniform validation",
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    work = cw.prepare_worktree(root / "build", MXSTEP)
    records = existing
    print(json.dumps({"accepted_saved": len(accepted_ids), "pending": len(pending), "total": len(targets)}, ensure_ascii=False), flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(run_target, work, root, target, index): target for index, target in pending}
        for completed, future in enumerate(as_completed(futures), 1):
            target = futures[future]
            try:
                result = future.result()
            except Exception as exc:  # Preserve the other rescue cases if one wrapper fails.
                result = {**target, "final_status": "MXSTEP_RESCUE_RUNTIME_ERROR", "runtime_error": repr(exc), "ended_at": now()}
            case_id = str(result["case_id"])
            records = [row for row in records if str(row.get("case_id")) != case_id] + [result]
            write_csv(summary_path, records)
            print(json.dumps({"completed": completed, "total_pending": len(pending), "case_id": case_id,
                              "status": result["final_status"]}, ensure_ascii=False), flush=True)
    records.sort(key=lambda item: str(item["case_id"]))
    write_csv(summary_path, records)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
