#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Uniform-NH3-precision validation for the remaining CW-CSTR cases.

This is the final numerical-quality layer for a completed 5832-case pure-gas
CW-CSTR grid.  It selects only cases whose final original calculation used
the relaxed third-attempt policy, excluding five cases already validated at
RTOL=1e-4.  The four original numerical failures are included through their
matching parameter tuples and prior rescue values.

All products are written to a new directory.  The original SQLite registry
and prior result directories are read-only inputs.  A manifest and a summary
are updated after every completed case, so completed work is retained if the
process is interrupted.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sqlite3
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Any


TOOL_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOL_DIR))
import hong_cw_longtime_scan as cw  # noqa: E402


SCAN_ROOT = Path(r"F:\Codex\PK1\GasReaction\CSTR\cw_gasphase_cstr_tau_1to100ms_20260907")
DEFAULT_OUTPUT_ROOT = SCAN_ROOT / "final_uniform_validation_adaptive_r4_20260912"
R3_RESCUE_SUMMARY = SCAN_ROOT / "numerical_rescue_20260912_r3" / "rescue_summary.csv"
UNIFORM5_SUMMARY = SCAN_ROOT / "uniform_relative_tolerance_rescue_20260912" / "uniform_tolerance_summary.csv"
R2_SUMMARY = SCAN_ROOT / "final_uniform_validation_adaptive_r2_20260912" / "uniform_validation_summary.csv"

# The absolute tolerance prevents unresolved trace radicals from stalling the
# integrator.  The high-concentration NH3 result is controlled by RTOL.
ATOL = 100.0
RTOL = 1.0e-4
MXSTEP = 500000
MIN_INITIAL_STAGE_S = 6.0
RESIDENCE_TIMES_AT_INITIAL_STAGE = 200.0
MAX_STAGE_S = 100.0
REPORT_DT_S = 1.0
STARTUP_DT_S = 1.0e-3
STARTUP_DURATION_S = 2.0e-1
TRANSITION_DT_S = 1.0e-1
TRANSITION_DURATION_S = 10.0
REFERENCE_REL_DIFF_LIMIT = 5.0e-4


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def canonical_id(tg: float, en: float, n2: float, tau: float) -> str:
    return f"T{tg:g}_EN{en:g}_N2{n2:g}_TAU{1000.0 * tau:g}ms"


def parameter_key(tg: float, en: float, n2: float, tau: float) -> tuple[float, float, float, float]:
    return (round(tg, 8), round(en, 8), round(n2, 8), round(tau, 8))


def validation_stages(tau_res_s: float) -> tuple[float, ...]:
    """Use at least 200 residence times, with 6 s minimum for chemistry."""
    initial = max(MIN_INITIAL_STAGE_S, RESIDENCE_TIMES_AT_INITIAL_STAGE * tau_res_s)
    if initial >= MAX_STAGE_S:
        return (MAX_STAGE_S,)
    return (initial, MAX_STAGE_S)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    fields = sorted({key for record in records for key in record})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def replace_case_record(records: list[dict[str, Any]], replacement: dict[str, Any]) -> list[dict[str, Any]]:
    """Keep one current validation record per case while preserving old attempts on disk."""
    case_id = str(replacement["case_id"])
    return [record for record in records if str(record.get("case_id")) != case_id] + [replacement]


def final_attempt_policies(connection: sqlite3.Connection) -> dict[tuple[str, int], str]:
    values: dict[tuple[str, int], str] = {}
    for event in connection.execute("SELECT case_id,attempt,details_json FROM events WHERE event='ATTEMPT_STARTED'"):
        details = json.loads(event["details_json"])
        values[(event["case_id"], int(event["attempt"]))] = str(details.get("solver_policy", "unknown"))
    return values


def load_targets(scan_root: Path) -> list[dict[str, Any]]:
    connection = sqlite3.connect(scan_root / "case_registry.sqlite")
    connection.row_factory = sqlite3.Row
    policies = final_attempt_policies(connection)
    uniform5_keys = {
        parameter_key(float(row["tg_K"]), float(row["en_Td"]), float(row["n2_fraction"]), float(row["tau_res_s"]))
        for row in read_csv(UNIFORM5_SUMMARY)
    }
    targets: dict[tuple[float, float, float, float], dict[str, Any]] = {}
    for row in connection.execute("SELECT * FROM cases WHERE status='ACCEPTED'"):
        policy = policies[(row["case_id"], int(row["total_attempt"]))]
        if policy != "relaxed_third_attempt":
            continue
        key = parameter_key(row["tg_K"], row["en_Td"], row["n2_fraction"], row["tau_res_s"])
        if key in uniform5_keys:
            raise RuntimeError(f"Already-uniform case unexpectedly appears in direct accepted selection: {row['case_id']}")
        targets[key] = {
            "case_id": canonical_id(*key),
            "original_case_id": row["case_id"],
            "source_class": "original_relaxed_accepted",
            "tg_K": row["tg_K"], "en_Td": row["en_Td"], "n2_fraction": row["n2_fraction"],
            "tau_res_s": row["tau_res_s"], "reference_NH3_cm3": row["NH3_cm3"],
            "reference_path": str(scan_root / row["output_group"]),
        }
    rescue_by_key = {
        parameter_key(float(row["tg_K"]), float(row["en_Td"]), float(row["n2_fraction"]), float(row["tau_res_s"])): row
        for row in read_csv(R3_RESCUE_SUMMARY)
    }
    for row in connection.execute("SELECT * FROM cases WHERE status='FAILED'"):
        key = parameter_key(row["tg_K"], row["en_Td"], row["n2_fraction"], row["tau_res_s"])
        rescue = rescue_by_key.get(key)
        if rescue is None or rescue.get("final_status") != "ACCEPTED":
            raise RuntimeError(f"No accepted r3 rescue for original failure {row['case_id']}")
        if key in targets:
            raise RuntimeError(f"Duplicate target key {key}")
        targets[key] = {
            "case_id": canonical_id(*key),
            "original_case_id": row["case_id"],
            "source_class": "r3_rescue_requires_uniform_validation",
            "tg_K": row["tg_K"], "en_Td": row["en_Td"], "n2_fraction": row["n2_fraction"],
            "tau_res_s": row["tau_res_s"], "reference_NH3_cm3": float(rescue["NH3_cm-3"]),
            "reference_path": str(R3_RESCUE_SUMMARY),
        }
    result = sorted(targets.values(), key=lambda item: (item["en_Td"], item["tau_res_s"], item["tg_K"], item["n2_fraction"]))
    expected = 78
    if len(result) != expected:
        raise RuntimeError(f"Expected {expected} remaining non-uniform cases, found {len(result)}")
    return result


def load_r2_seed_records(targets: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Reuse only r2 cases that already passed the identical solver criterion."""
    if not R2_SUMMARY.is_file():
        return []
    target_ids = {str(item["case_id"]) for item in targets}
    records = [dict(row) for row in read_csv(R2_SUMMARY)]
    if any(record.get("final_status") != "ACCEPTED_FINAL_UNIFORM" for record in records):
        raise RuntimeError("r2 summary contains a non-accepted record and cannot be used as a seed")
    if not all(str(record.get("case_id")) in target_ids for record in records):
        raise RuntimeError("r2 summary contains a case outside the current 78-case manifest")
    for record in records:
        record["seeded_from"] = str(R2_SUMMARY)
    return records


def terminal_block_test(root: Path, record: dict[str, Any]) -> dict[str, Any]:
    tag = str(record["tag"])
    series_path = root / str(record["group"]) / tag / f"cw_concentrations_{tag}.csv"
    series = cw.read_series(series_path)
    t_end = float(record["t_final_s"])
    window = min(cw.TERMINAL_WINDOW_S, 0.2 * t_end)
    split = t_end - window / 2.0
    first = [row["NH3_cm-3"] for row in series if t_end - window <= row["time_s"] < split]
    second = [row["NH3_cm-3"] for row in series if split <= row["time_s"] <= t_end]
    if first and second:
        mean_first = math.fsum(first) / len(first)
        mean_second = math.fsum(second) / len(second)
        mode = "terminal_half_window_means"
    else:
        # The driver outputs at every 1 s only during the 0.2--10 s
        # transition.  Later outer 10 s calls record their endpoint, so a
        # 20 s trajectory contains valid independent endpoints at 10 and 20 s
        # but no samples in the nominal 16--18 s subwindow.  Use those two
        # adjacent outer endpoints rather than treating the output cadence as
        # a physical convergence failure.
        if len(series) < 2:
            raise RuntimeError("Fewer than two saved concentration outputs")
        previous, latest = series[-2], series[-1]
        if latest["time_s"] <= previous["time_s"]:
            raise RuntimeError("Concentration output time is not increasing")
        mean_first = previous["NH3_cm-3"]
        mean_second = latest["NH3_cm-3"]
        mode = "last_two_saved_outer_endpoints"
    delta = abs(mean_second - mean_first)
    allowed = max(ATOL, RTOL * abs(mean_second))
    return {
        "NH3_first_block_mean_cm3": mean_first,
        "NH3_second_block_mean_cm3": mean_second,
        "NH3_block_mean_abs_delta_cm3": delta,
        "NH3_block_mean_rel_delta": delta / max(abs(mean_second), 1.0),
        "allowed_abs_delta_cm3": allowed,
        "stationary_by_uniform_error_aware_test": delta <= allowed,
        "stationarity_test_mode": mode,
    }


def retry_number(root: Path, target: dict[str, Any], target_index: int) -> int:
    """Choose a new group if an interrupted run left an incomplete directory."""
    tag = cw.case_name(float(target["tg_K"]), float(target["en_Td"]), float(target["n2_fraction"]), 10.0,
                       float(target["tau_res_s"]))
    prefix = f"attempt_{target_index:03d}"
    base = root / "uniform_rtol1e-4"
    attempts = 0
    if base.is_dir():
        for group in base.glob(f"{prefix}*"):
            if any((candidate / tag).exists() for candidate in group.glob("stage_*s")):
                attempts += 1
    return attempts + 1


def run_target(work: Path, root: Path, target: dict[str, Any], target_index: int, retry: int) -> dict[str, Any]:
    suffix = f"attempt_{target_index:03d}" if retry == 1 else f"attempt_{target_index:03d}_retry_{retry:03d}"
    base_group = f"uniform_rtol1e-4/{suffix}"
    result: dict[str, Any] = {
        **target,
        "target_index": target_index,
        "retry": retry,
        "solver_policy": "uniform_NH3_rtol1e-4_atol100_with_mxstep500000",
        "stages_attempted": 0,
    }
    stages = validation_stages(float(target["tau_res_s"]))
    for stage_s in stages:
        group = f"{base_group}/stage_{stage_s:g}s"
        record = cw.run_case(
            work=work, root=root, tg=float(target["tg_K"]), en=float(target["en_Td"]),
            n2=float(target["n2_fraction"]), dt=10.0, t_end=stage_s, report_dt=REPORT_DT_S,
            group=group, tau_res_s=float(target["tau_res_s"]), atol=ATOL, rtol=RTOL,
            startup_dt_s=STARTUP_DT_S, startup_duration_s=STARTUP_DURATION_S,
            transition_dt_s=TRANSITION_DT_S, transition_duration_s=TRANSITION_DURATION_S,
            timeout_s=21600.0, idle_timeout_s=300.0, watchdog_poll_s=5.0,
        )
        result.update(record)
        result["stages_attempted"] = int(result["stages_attempted"]) + 1
        result["accepted_stage_s"] = stage_s
        if record.get("returncode") != 0 or bool(record.get("timed_out")):
            result["final_status"] = "NUMERICAL_FAILURE_UNIFORM_VALIDATION"
            return result
        result.update(terminal_block_test(root, record))
        final_nh3 = float(record["NH3_cm-3"])
        reference = float(target["reference_NH3_cm3"])
        relative_difference = abs(final_nh3 - reference) / max(abs(reference), 1.0)
        result["NH3_rel_diff_to_reference"] = relative_difference
        result["reference_rel_diff_limit"] = REFERENCE_REL_DIFF_LIMIT
        if bool(result["stationary_by_uniform_error_aware_test"]) and relative_difference <= REFERENCE_REL_DIFF_LIMIT:
            result["final_status"] = "ACCEPTED_FINAL_UNIFORM"
            return result
    if not bool(result.get("stationary_by_uniform_error_aware_test")):
        result["final_status"] = "REQUIRES_LONGER_OR_MORE_STABLE_UNIFORM_RUN"
    else:
        result["final_status"] = "REQUIRES_REFERENCE_DISAGREEMENT_REVIEW"
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Finalize the remaining CW-CSTR relaxed-tolerance cases")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true", help="continue only targets absent from the saved summary")
    args = parser.parse_args()
    if args.workers <= 0:
        parser.error("workers must be positive")
    targets = load_targets(SCAN_ROOT)
    seed_records = load_r2_seed_records(targets)
    if args.dry_run:
        print(json.dumps({
            "targets": len(targets),
            "source_classes": {
                key: sum(item["source_class"] == key for item in targets)
                for key in sorted({item["source_class"] for item in targets})
            },
            "first_targets": targets[:5],
        }, ensure_ascii=False, indent=2))
        return 0
    root = args.output_root.resolve()
    existing_records: list[dict[str, Any]] = []
    if root.exists():
        if not args.resume:
            raise RuntimeError(f"Refusing to overwrite existing output root: {root}; use --resume to continue it")
        manifest_path = root / "target_manifest.csv"
        if not manifest_path.is_file():
            raise RuntimeError(f"Cannot resume without target manifest: {manifest_path}")
        expected_ids = {str(item["case_id"]) for item in targets}
        stored_ids = {str(item["case_id"]) for item in read_csv(manifest_path)}
        if expected_ids != stored_ids:
            raise RuntimeError("Current target selection does not match the saved manifest")
        summary_path = root / "uniform_validation_summary.csv"
        if summary_path.is_file():
            existing_records = [dict(item) for item in read_csv(summary_path)]
            interrupted = [
                item for item in existing_records
                if str(item.get("final_status")) != "ACCEPTED_FINAL_UNIFORM"
            ]
            if interrupted:
                audit_path = root / "uniform_validation_requeue_audit.csv"
                prior_audit = read_csv(audit_path) if audit_path.is_file() else []
                prior_signatures = {
                    json.dumps({key: value for key, value in item.items() if key != "requeued_at"},
                               ensure_ascii=False, sort_keys=True)
                    for item in prior_audit
                }
                additions = []
                for item in interrupted:
                    audit_record = {**item, "requeued_at": now()}
                    signature = json.dumps(item, ensure_ascii=False, sort_keys=True)
                    if signature not in prior_signatures:
                        additions.append(audit_record)
                if additions:
                    write_csv(audit_path, prior_audit + additions)
    else:
        root.mkdir(parents=True)
        write_csv(root / "target_manifest.csv", targets)
        existing_records = seed_records
        if seed_records:
            write_csv(root / "seeded_r2_accepted_records.csv", seed_records)
    protocol = {
        "created_at": now(), "source_scan_root": str(SCAN_ROOT),
        "purpose": "uniform NH3 relative-precision validation of the remaining 78 non-uniform cases",
        "targets": len(targets), "seeded_r2_accepted_records": len(seed_records), "workers": args.workers,
        "atol_cm3": ATOL, "rtol": RTOL, "mxstep": MXSTEP,
        "initial_stage_rule": {
            "minimum_s": MIN_INITIAL_STAGE_S,
            "minimum_residence_times": RESIDENCE_TIMES_AT_INITIAL_STAGE,
            "formula": "max(6 s, 200*tau_res_s)",
        },
        "fallback_stage_s": MAX_STAGE_S, "report_dt_s": REPORT_DT_S,
        "startup": {"dt_s": STARTUP_DT_S, "duration_s": STARTUP_DURATION_S},
        "transition": {"dt_s": TRANSITION_DT_S, "duration_s": TRANSITION_DURATION_S},
        "criterion": "at each stage, compare terminal half-window means when available; otherwise compare the final two saved outer-step endpoints; require difference <= max(ATOL, RTOL*abs(mean_last))",
        "reference_relative_difference_limit": REFERENCE_REL_DIFF_LIMIT,
        "original_registry_modified": False,
    }
    protocol_path = root / "PROTOCOL.json"
    if not protocol_path.exists():
        protocol_path.write_text(json.dumps(protocol, ensure_ascii=False, indent=2), encoding="utf-8")
    work = cw.prepare_worktree(root / "build", MXSTEP)
    records: list[dict[str, Any]] = existing_records
    # A manual pause terminates external solvers intentionally.  Their rows
    # must remain auditable, but cannot be treated as scientifically completed
    # cases on a later --resume invocation.
    completed_ids = {
        str(item["case_id"])
        for item in existing_records
        if str(item.get("final_status")) == "ACCEPTED_FINAL_UNIFORM"
    }
    pending = [(index, target) for index, target in enumerate(targets, 1) if str(target["case_id"]) not in completed_ids]
    print(json.dumps({"saved_accepted": len(completed_ids), "pending": len(pending), "total": len(targets)}, ensure_ascii=False), flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(run_target, work, root, target, index, retry_number(root, target, index)): (index, target)
            for index, target in pending
        }
        for completed, future in enumerate(as_completed(futures), 1):
            target_index, target = futures[future]
            try:
                result = future.result()
            except Exception as exc:  # Keep other concurrent cases running and requeue this one on --resume.
                result = {
                    **target,
                    "target_index": target_index,
                    "final_status": "RUNTIME_ERROR_REQUEUED_ON_RESUME",
                    "runtime_error": repr(exc),
                }
            records = replace_case_record(records, result)
            write_csv(root / "uniform_validation_summary.csv", records)
            print(json.dumps({
                "completed_this_run": completed, "pending_at_start": len(pending), "total": len(targets), "case_id": result["case_id"],
                "final_status": result["final_status"],
            }, ensure_ascii=False), flush=True)
    records.sort(key=lambda item: str(item["case_id"]))
    write_csv(root / "uniform_validation_summary.csv", records)
    counts: dict[str, int] = {}
    for record in records:
        status = str(record["final_status"])
        counts[status] = counts.get(status, 0) + 1
    (root / "COMPLETION.json").write_text(json.dumps({
        "finished_at": now(), "total": len(records), "status_counts": counts,
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
