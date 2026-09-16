#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recheck CW-CSTR terminal stationarity without altering the scan registry.

The original adaptive scan used the range of every terminal sample as a
convergence gate.  This tool instead tests the change between the means of
two adjacent terminal blocks against the actual DVODE absolute/relative
tolerances recorded for each run.  It only reads existing scan products and
writes a separate reconciliation folder.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sqlite3
from collections import Counter
from datetime import datetime
from pathlib import Path


DEFAULT_SCAN_ROOT = Path(
    r"F:\Codex\PK1\GasReaction\CSTR\cw_gasphase_cstr_tau_1to100ms_20260907"
)
DEFAULT_OUTPUT_ROOT = DEFAULT_SCAN_ROOT / "convergence_recheck_20260912"
FLAGGED_STATUS = "ACCEPTED_AT_MAX_TIME_WITH_NH3_FLAG"


def token(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")


def case_tag(row: sqlite3.Row) -> str:
    return (
        f"T{token(row['tg_K'])}_EN{token(row['en_Td'])}_"
        f"N2{token(row['n2_fraction'])}_dt10s_tau{token(1.0e3 * row['tau_res_s'])}ms"
    )


def read_json(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def read_series(path: Path) -> list[dict[str, float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(handle)]


def assess_case(scan_root: Path, row: sqlite3.Row) -> dict[str, object]:
    tag = case_tag(row)
    run_root = scan_root / row["output_group"] / tag
    params = read_json(run_root / "params.json")
    series = read_series(run_root / f"cw_concentrations_{tag}.csv")
    t_end = float(params["t_final_s"])
    terminal_window = float(params["terminal_window_s"])
    split = t_end - terminal_window / 2.0
    first = [sample["NH3_cm-3"] for sample in series if t_end - terminal_window <= sample["time_s"] < split]
    second = [sample["NH3_cm-3"] for sample in series if split <= sample["time_s"] <= t_end]
    if not first or not second:
        raise RuntimeError(f"Insufficient terminal samples for {row['case_id']}")
    first_mean = math.fsum(first) / len(first)
    second_mean = math.fsum(second) / len(second)
    delta = abs(second_mean - first_mean)
    atol = float(params["atol"])
    rtol = float(params["rtol"])
    numerical_limit = max(atol, rtol * abs(second_mean))
    stationary = delta <= numerical_limit
    policy = "baseline" if (atol, rtol) == (1.0, 1.0e-4) else "fallback_or_nonbaseline"
    if policy == "baseline" and stationary:
        recheck_status = "ACCEPTED_RECHECKED"
    elif policy == "baseline":
        recheck_status = "REQUIRES_TARGETED_EXTENSION_OR_TIGHTER_SOLVER"
    elif stationary:
        recheck_status = "STATIONARY_BUT_REQUIRES_UNIFORM_TOLERANCE_RECALC"
    else:
        recheck_status = "REQUIRES_UNIFORM_TOLERANCE_RECALC"
    return {
        "case_id": row["case_id"],
        "tg_K": row["tg_K"],
        "en_Td": row["en_Td"],
        "n2_fraction": row["n2_fraction"],
        "tau_res_s": row["tau_res_s"],
        "original_status": row["status"],
        "output_group": row["output_group"],
        "atol_cm3": atol,
        "rtol": rtol,
        "terminal_window_s": terminal_window,
        "first_block_samples": len(first),
        "second_block_samples": len(second),
        "NH3_first_block_mean_cm3": first_mean,
        "NH3_second_block_mean_cm3": second_mean,
        "NH3_block_mean_abs_delta_cm3": delta,
        "NH3_block_mean_rel_delta": delta / max(abs(second_mean), 1.0),
        "allowed_abs_delta_cm3": numerical_limit,
        "stationary_by_error_aware_test": stationary,
        "solver_policy": policy,
        "recheck_status": recheck_status,
        "original_NH3_cm3": row["NH3_cm3"],
        "original_NH3_tail_rel_range": row["NH3_tail_rel_range"],
    }


def write_csv(path: Path, records: list[dict[str, object]]) -> None:
    fields = list(records[0]) if records else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def main() -> int:
    parser = argparse.ArgumentParser(description="Non-destructive stationarity recheck for CW-CSTR cases")
    parser.add_argument("--scan-root", type=Path, default=DEFAULT_SCAN_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args()
    scan_root = args.scan_root.resolve()
    output_root = args.output_root.resolve()
    if output_root.exists():
        raise RuntimeError(f"Refusing to overwrite existing output directory: {output_root}")
    output_root.mkdir(parents=True)
    connection = sqlite3.connect(scan_root / "case_registry.sqlite")
    connection.row_factory = sqlite3.Row
    rows = list(connection.execute("SELECT * FROM cases WHERE status=? ORDER BY case_id", (FLAGGED_STATUS,)))
    records = [assess_case(scan_root, row) for row in rows]
    write_csv(output_root / "convergence_recheck.csv", records)
    status_counts = Counter(str(record["recheck_status"]) for record in records)
    report = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "source_scan_root": str(scan_root),
        "source_registry": str(scan_root / "case_registry.sqlite"),
        "source_registry_modified": False,
        "original_status_rechecked": FLAGGED_STATUS,
        "criterion": {
            "description": "absolute difference between means of adjacent terminal half-windows",
            "pass_when": "abs(mean_last - mean_previous) <= max(ATOL, RTOL * abs(mean_last))",
        },
        "n_cases": len(records),
        "status_counts": dict(sorted(status_counts.items())),
    }
    (output_root / "recheck_protocol.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    lines = [
        "# CW-CSTR convergence recheck",
        "",
        "This folder is a read-only reassessment of completed trajectories; the source SQLite registry was not changed.",
        "",
        "## Criterion",
        "",
        "For NH3, compare the means of the two adjacent halves of the recorded terminal window. "
        "A case is stationary when `abs(mean_last - mean_previous) <= max(ATOL, RTOL * abs(mean_last))`, "
        "using the exact ATOL and RTOL recorded in that run's `params.json`.",
        "",
        "## Result counts",
        "",
    ]
    lines.extend(f"- `{key}`: {value}" for key, value in sorted(status_counts.items()))
    lines.extend(["", "Detailed per-case values: `convergence_recheck.csv`."])
    (output_root / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
