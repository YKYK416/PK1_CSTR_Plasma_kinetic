#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recalculate the five non-baseline CW-CSTR cases with uniform NH3 accuracy.

These cases were numerically completed in the original scan only after its
relaxed fallback tolerance was used.  The script preserves that scan and
performs an isolated 100 s precision verification using the baseline NH3
relative tolerance and a larger DVODE MXSTEP budget.  A 100 cm-3 absolute
tolerance is retained for trace radicals because an all-species 1 cm-3
absolute tolerance stalls before the first output; it does not limit the
relative accuracy of the high-concentration NH3 target.  Existing
time-history evidence shows 100 s is sufficient for the physical steady
state; this is a numerical precision check, not a longer-time extension.
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


TOOL_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOL_DIR))
import hong_cw_longtime_scan as cw  # noqa: E402


DEFAULT_ROOT = Path(
    r"F:\Codex\PK1\GasReaction\CSTR\cw_gasphase_cstr_tau_1to100ms_20260907"
    r"\uniform_relative_tolerance_rescue_20260912"
)
TARGETS = (
    (375.0, 60.0, 0.2, 0.050),
    (500.0, 100.0, 0.5, 0.005),
    (500.0, 100.0, 0.5, 0.030),
    (500.0, 100.0, 0.5, 0.050),
    (500.0, 100.0, 0.5, 0.100),
)
T_END_S = 100.0
ATOL = 100.0
RTOL = 1.0e-4
MXSTEP = 500000
STARTUP_DT_S = 1.0e-3
STARTUP_DURATION_S = 2.0e-1
TRANSITION_DT_S = 1.0e-1
TRANSITION_DURATION_S = 10.0


def block_stationarity(root: Path, record: dict[str, object]) -> dict[str, float | bool]:
    tag = str(record["tag"])
    series_path = root / str(record["group"]) / tag / f"cw_concentrations_{tag}.csv"
    series = cw.read_series(series_path)
    t_end = float(record["t_final_s"])
    terminal_window = min(cw.TERMINAL_WINDOW_S, 0.2 * t_end)
    split = t_end - terminal_window / 2.0
    first = [row["NH3_cm-3"] for row in series if t_end - terminal_window <= row["time_s"] < split]
    second = [row["NH3_cm-3"] for row in series if split <= row["time_s"] <= t_end]
    if not first or not second:
        raise RuntimeError("No samples in one of the terminal convergence blocks")
    first_mean = math.fsum(first) / len(first)
    second_mean = math.fsum(second) / len(second)
    delta = abs(second_mean - first_mean)
    allowed = max(ATOL, RTOL * abs(second_mean))
    return {
        "NH3_first_block_mean_cm3": first_mean,
        "NH3_second_block_mean_cm3": second_mean,
        "NH3_block_mean_abs_delta_cm3": delta,
        "allowed_abs_delta_cm3": allowed,
        "stationary_by_uniform_error_aware_test": delta <= allowed,
    }


def write_csv(path: Path, records: list[dict[str, object]]) -> None:
    fields = sorted({key for record in records for key in record})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def run_target(work: Path, root: Path, target: tuple[float, float, float, float]) -> dict[str, object]:
    tg, en, n2, tau = target
    record = cw.run_case(
        work=work, root=root, tg=tg, en=en, n2=n2, dt=10.0, t_end=T_END_S,
        report_dt=1.0, group="uniform_tolerance/t100s", tau_res_s=tau,
        atol=ATOL, rtol=RTOL, startup_dt_s=STARTUP_DT_S,
        startup_duration_s=STARTUP_DURATION_S, transition_dt_s=TRANSITION_DT_S,
        transition_duration_s=TRANSITION_DURATION_S, timeout_s=21600.0,
        idle_timeout_s=300.0, watchdog_poll_s=5.0,
    )
    result: dict[str, object] = {
        "case_id": f"T{tg:g}_EN{en:g}_N2{n2:g}_TAU{tau * 1e3:g}ms",
        "tg_K": tg, "en_Td": en, "n2_fraction": n2, "tau_res_s": tau,
        "solver_policy": "uniform_NH3_rtol1e-4_atol100_with_mxstep500000",
        **record,
    }
    if record.get("returncode") == 0 and not bool(record.get("timed_out")):
        result.update(block_stationarity(root, record))
        result["final_status"] = (
            "ACCEPTED_UNIFORM_NH3_RELATIVE_TOLERANCE"
            if bool(result["stationary_by_uniform_error_aware_test"])
            else "REQUIRES_LONGER_OR_MORE_STABLE_UNIFORM_RUN"
        )
    else:
        result["final_status"] = "NUMERICAL_FAILURE_UNIFORM_NH3_RELATIVE_TOLERANCE"
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Uniform-NH3-relative-tolerance verification for five CW-CSTR cases")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers <= 0:
        parser.error("workers must be positive")
    root = args.output_root.resolve()
    if root.exists():
        raise RuntimeError(f"Refusing to overwrite existing output directory: {root}")
    root.mkdir(parents=True)
    work = cw.prepare_worktree(root / "build", MXSTEP)
    protocol = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "purpose": "uniform NH3-relative-tolerance numerical precision verification; original scan remains unchanged",
        "targets": TARGETS,
        "t_end_s": T_END_S,
        "atol_cm3": ATOL,
        "rtol": RTOL,
        "mxstep": MXSTEP,
        "startup": {"dt_s": STARTUP_DT_S, "duration_s": STARTUP_DURATION_S},
        "transition": {"dt_s": TRANSITION_DT_S, "duration_s": TRANSITION_DURATION_S},
    }
    (root / "PROTOCOL.json").write_text(json.dumps(protocol, ensure_ascii=False, indent=2), encoding="utf-8")
    records: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(run_target, work, root, target) for target in TARGETS]
        for future in as_completed(futures):
            record = future.result()
            records.append(record)
            write_csv(root / "uniform_tolerance_summary.csv", records)
            print(json.dumps({"case_id": record["case_id"], "final_status": record["final_status"]}, ensure_ascii=False), flush=True)
    write_csv(root / "uniform_tolerance_summary.csv", sorted(records, key=lambda item: str(item["case_id"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
