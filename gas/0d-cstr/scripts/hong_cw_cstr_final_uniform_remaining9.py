#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Complete the nine interrupted final-uniform CW-CSTR validation cases.

Targets are derived from the r4 manifest after excluding both r4 accepted
records and the three separately accepted 1 s-segment MXSTEP rescues.  Outputs
are isolated from every prior r4 attempt.  The tau=100 ms subset uses the
validated 1 s outer segment; the tau=10 ms case retains the original 10 s
outer segment and two-stage r4 horizon policy.
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
R4_ROOT = SCAN_ROOT / "final_uniform_validation_adaptive_r4_20260912"
MXSTEP_RESCUE_ROOT = SCAN_ROOT / "tau100_mxstep_rescue_20260914"
DEFAULT_OUTPUT_ROOT = SCAN_ROOT / "final_uniform_validation_remaining9_20260914"
MXSTEP = 500000
ATOL = 100.0
RTOL = 1.0e-4
REPORT_DT_S = 1.0
STARTUP_DT_S = 1.0e-3
STARTUP_DURATION_S = 2.0e-1
TRANSITION_DT_S = 1.0e-1
TRANSITION_DURATION_S = 10.0
CASE_WALL_TIMEOUT_S = 86400.0
EXPECTED_IDS = {
    "T400_EN60_N20.2_TAU100ms", "T425_EN60_N20.1_TAU100ms", "T425_EN60_N20.2_TAU100ms",
    "T450_EN60_N20.2_TAU100ms", "T475_EN60_N20.1_TAU100ms", "T475_EN60_N20.2_TAU100ms",
    "T500_EN60_N20.1_TAU100ms", "T500_EN60_N20.2_TAU100ms", "T500_EN100_N20.5_TAU10ms",
}


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


def load_targets() -> list[dict[str, Any]]:
    manifest = read_csv(R4_ROOT / "target_manifest.csv")
    r4 = read_csv(R4_ROOT / "uniform_validation_summary.csv")
    rescue = read_csv(MXSTEP_RESCUE_ROOT / "mxstep_rescue_summary.csv")
    accepted = {
        str(row["case_id"]) for row in r4
        if row.get("final_status") == "ACCEPTED_FINAL_UNIFORM"
    }
    accepted.update(
        str(row["case_id"]) for row in rescue
        if row.get("final_status") == "ACCEPTED_MXSTEP_SEGMENTED_UNIFORM"
    )
    targets = [dict(row) for row in manifest if str(row["case_id"]) not in accepted]
    found = {str(row["case_id"]) for row in targets}
    if found != EXPECTED_IDS:
        raise RuntimeError(f"Expected exactly nine interrupted targets, found {sorted(found)}")
    return sorted(targets, key=lambda row: (float(row["tau_res_s"]), float(row["tg_K"]), float(row["n2_fraction"])))


def retry_number(root: Path, case_id: str) -> int:
    return 1 + sum(1 for item in root.glob(f"{case_id}_attempt_*") if item.is_dir())


def outer_dt(tau_res_s: float) -> float:
    return 1.0 if tau_res_s >= 0.1 else 10.0


def run_target(work: Path, root: Path, target: dict[str, Any], index: int, retry: int) -> dict[str, Any]:
    tau = float(target["tau_res_s"])
    dt = outer_dt(tau)
    result: dict[str, Any] = {
        **target, "remaining_index": index, "retry": retry, "outer_dt_s": dt,
        "solver_policy": "uniform_NH3_rtol1e-4_atol100_segmented_tau100", "started_at": now(),
    }
    base = f"{target['case_id']}_attempt_{retry:02d}"
    for stage_s in uniform.validation_stages(tau):
        record = cw.run_case(
            work=work, root=root, tg=float(target["tg_K"]), en=float(target["en_Td"]),
            n2=float(target["n2_fraction"]), dt=dt, t_end=stage_s, report_dt=REPORT_DT_S,
            group=f"{base}/stage_{stage_s:g}s", tau_res_s=tau, atol=ATOL, rtol=RTOL,
            startup_dt_s=STARTUP_DT_S, startup_duration_s=STARTUP_DURATION_S,
            transition_dt_s=TRANSITION_DT_S, transition_duration_s=TRANSITION_DURATION_S,
            timeout_s=CASE_WALL_TIMEOUT_S, idle_timeout_s=300.0, watchdog_poll_s=5.0,
        )
        result.update(record)
        result["accepted_stage_s"] = stage_s
        if record.get("returncode") != 0 or bool(record.get("timed_out")):
            result["final_status"] = "NUMERICAL_FAILURE_REMAINING9"
            result["ended_at"] = now()
            return result
        result.update(uniform.terminal_block_test(root, record))
        final_nh3 = float(record["NH3_cm-3"])
        reference = float(target["reference_NH3_cm3"])
        rel_diff = abs(final_nh3 - reference) / max(abs(reference), 1.0)
        result["NH3_rel_diff_to_reference"] = rel_diff
        result["reference_rel_diff_limit"] = uniform.REFERENCE_REL_DIFF_LIMIT
        if bool(result["stationary_by_uniform_error_aware_test"]) and rel_diff <= uniform.REFERENCE_REL_DIFF_LIMIT:
            result["final_status"] = "ACCEPTED_REMAINING9_UNIFORM"
            result["ended_at"] = now()
            return result
    result["final_status"] = (
        "REQUIRES_LONGER_TIME_REMAINING9"
        if not bool(result.get("stationary_by_uniform_error_aware_test"))
        else "REQUIRES_REFERENCE_REVIEW_REMAINING9"
    )
    result["ended_at"] = now()
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Complete the nine interrupted r4 final-uniform cases")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.workers <= 0:
        parser.error("workers must be positive")
    targets = load_targets()
    if args.dry_run:
        print(json.dumps({"targets": targets, "outer_dt_s": {item["case_id"]: outer_dt(float(item["tau_res_s"])) for item in targets}}, ensure_ascii=False, indent=2))
        return 0
    root = args.output_root.resolve()
    summary_path = root / "remaining9_summary.csv"
    existing: list[dict[str, Any]] = []
    if root.exists():
        if not args.resume:
            raise RuntimeError(f"Refusing to overwrite existing output root: {root}; use --resume")
        if summary_path.is_file():
            existing = [dict(item) for item in read_csv(summary_path)]
    else:
        root.mkdir(parents=True)
    accepted_ids = {
        str(row["case_id"]) for row in existing
        if row.get("final_status") == "ACCEPTED_REMAINING9_UNIFORM"
    }
    pending = [(index, target) for index, target in enumerate(targets, 1) if str(target["case_id"]) not in accepted_ids]
    (root / "PROTOCOL.json").write_text(json.dumps({
        "created_at": now(), "purpose": "complete nine interrupted r4 uniform-validation targets",
        "targets": [item["case_id"] for item in targets], "rtol": RTOL, "atol": ATOL,
        "case_wall_timeout_s": CASE_WALL_TIMEOUT_S,
        "tau100_outer_dt_s": 1.0, "other_outer_dt_s": 10.0,
        "stages": "max(6 s, 200*tau), then 100 s only if required", "original_results_modified": False,
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    work = cw.prepare_worktree(root / "build", MXSTEP)
    records = existing
    print(json.dumps({"accepted_saved": len(accepted_ids), "pending": len(pending), "total": len(targets)}, ensure_ascii=False), flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(run_target, work, root, target, index, retry_number(root, str(target["case_id"]))): target
            for index, target in pending
        }
        for completed, future in enumerate(as_completed(futures), 1):
            target = futures[future]
            try:
                result = future.result()
            except Exception as exc:  # A single wrapper error must not abandon the other eight cases.
                result = {**target, "final_status": "RUNTIME_ERROR_REMAINING9", "runtime_error": repr(exc), "ended_at": now()}
            case_id = str(result["case_id"])
            records = [row for row in records if str(row.get("case_id")) != case_id] + [result]
            write_csv(summary_path, records)
            print(json.dumps({"completed": completed, "pending_at_start": len(pending), "case_id": case_id,
                              "status": result["final_status"]}, ensure_ascii=False), flush=True)
    records.sort(key=lambda item: str(item["case_id"]))
    write_csv(summary_path, records)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
