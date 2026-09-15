#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recover four MXSTEP-limited pure-gas CW-CSTR cases without touching the grid.

The original full scan is preserved as evidence.  This utility creates a new
worktree and retries only the four failed points using the grid's documented
relaxed tolerance policy, the validated 1 ms startup through 0.2 s, then a
0.1 s transition through 10 s before returning to 10 s segments.  Each stage
starts from t=0, matching the adaptive grid policy.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

TOOL_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOL_DIR))
import hong_cw_longtime_scan as cw  # noqa: E402

DEFAULT_ROOT = Path(r"F:\Codex\PK1\GasReaction\CSTR\cw_gasphase_cstr_tau_1to100ms_20260907\numerical_rescue_20260912_r3")
TARGETS = (
    (450.0, 60.0, 0.2, 0.005),
    (500.0, 60.0, 0.2, 0.005),
    (475.0, 60.0, 0.2, 0.010),
    (500.0, 60.0, 0.2, 0.010),
)
STAGES = (100.0, 500.0, 1000.0)
STARTUP_DT_S = 1.0e-3
STARTUP_DURATION_S = 2.0e-1
TRANSITION_DT_S = 1.0e-1
TRANSITION_DURATION_S = 10.0
MXSTEP = 500000


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def successful(record: dict[str, object]) -> bool:
    return record.get("returncode") == 0 and not bool(record.get("timed_out"))


def run_stage(work: Path, root: Path, target: tuple[float, float, float, float], stage_s: float) -> dict[str, object]:
    tg, en, n2, tau = target
    common = dict(
        work=work, root=root, tg=tg, en=en, n2=n2, dt=10.0, t_end=stage_s,
        report_dt=1.0, tau_res_s=tau, startup_dt_s=STARTUP_DT_S,
        startup_duration_s=STARTUP_DURATION_S, transition_dt_s=TRANSITION_DT_S,
        transition_duration_s=TRANSITION_DURATION_S, timeout_s=21600.0,
        idle_timeout_s=300.0, watchdog_poll_s=5.0,
    )
    relaxed = cw.run_case(group=f"relaxed/stage_{int(stage_s)}s", atol=100.0, rtol=1.0e-3, **common)
    relaxed["solver_policy"] = "relaxed_tolerance_meso_transition"
    relaxed["rescue_method"] = "startup_1ms_0p2s_then_0p1s_to_10s_mxstep500000"
    return relaxed


def run_target(work: Path, root: Path, target: tuple[float, float, float, float]) -> dict[str, object]:
    history: list[dict[str, object]] = []
    final: dict[str, object] | None = None
    final_status = "NUMERICAL_FAILURE"
    for stage_s in STAGES:
        record = run_stage(work, root, target, stage_s)
        history.append(record)
        final = record
        if not successful(record):
            break
        if bool(record.get("NH3_tail_pass")) and bool(record.get("other_tail_pass")):
            final_status = "ACCEPTED"
            break
        if stage_s == STAGES[-1] and bool(record.get("other_tail_pass")):
            final_status = "ACCEPTED_AT_MAX_TIME_WITH_NH3_FLAG"
        elif stage_s == STAGES[-1]:
            final_status = "NUMERICAL_NONCONVERGED"
    assert final is not None
    summary = {
        "case_id": f"T{target[0]:g}_EN{target[1]:g}_N2{target[2]:g}_TAU{target[3] * 1e3:g}ms",
        "tg_K": target[0], "en_Td": target[1], "n2_fraction": target[2],
        "h2_fraction": 1.0 - target[2], "tau_res_s": target[3],
        "final_status": final_status, "stages_attempted": len(history),
        "stage_history_json": json.dumps([
            {key: value for key, value in item.items() if key not in {"stdout", "stderr"}}
            for item in history
        ], ensure_ascii=False),
    }
    summary.update(final)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Numerical rescue for four CW-CSTR MXSTEP failures")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    if args.workers <= 0:
        parser.error("workers must be positive")
    root = args.output_root.resolve()
    # Keep the compile worktree separate so scheduler logs and rescue products
    # can coexist safely at the result root without weakening prepare_worktree's
    # non-overwrite guard.
    work = cw.prepare_worktree(root / "build", MXSTEP)
    protocol = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "purpose": "dedicated numerical rescue for four original-grid MXSTEP failures",
        "targets": TARGETS,
        "stages_s": STAGES,
        "solver_tolerances": {"atol": 100.0, "rtol": 1.0e-3},
        "startup_dt_s": STARTUP_DT_S,
        "startup_duration_s": STARTUP_DURATION_S,
        "transition_dt_s": TRANSITION_DT_S,
        "transition_duration_s": TRANSITION_DURATION_S,
        "report_interval_s": 1.0,
        "mxstep": MXSTEP,
        "original_grid_is_unchanged": True,
    }
    (root / "RESCUE_PROTOCOL.json").write_text(json.dumps(protocol, ensure_ascii=False, indent=2), encoding="utf-8")
    records: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run_target, work, root, target) for target in TARGETS]
        for future in as_completed(futures):
            record = future.result()
            records.append(record)
            print(json.dumps({"case_id": record["case_id"], "final_status": record["final_status"]}, ensure_ascii=False), flush=True)
            write_csv(root / "rescue_summary.csv", records)
    write_csv(root / "rescue_summary.csv", sorted(records, key=lambda item: str(item["case_id"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
