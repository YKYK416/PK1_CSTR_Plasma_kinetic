#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recoverable adaptive-time scan for the pure-gas Hong CW-CSTR model.

Every grid case starts at 100 s. Cases that do not satisfy the terminal
convergence gates are rerun from t=0 at 500 s, then at the user-selected
maximum of 1000 s. SQLite is the authoritative per-case registry; a readable
CSV snapshot and compact JSON progress file are exported periodically.

The solver does not provide state checkpoints. An interrupted attempt is
preserved and the same stage is restarted in a new attempt folder.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sqlite3
import sys
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime
from pathlib import Path

TOOL_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOL_DIR))
import hong_cw_longtime_scan as scan  # noqa: E402

DEFAULT_ROOT = Path(r"F:\Codex\PK1\GasReaction\CSTR\cw_gasphase_cstr_tau_1to100ms_20260907")
TEMPERATURES_K = tuple(float(value) for value in range(300, 501, 25))
EN_VALUES_TD = tuple(float(value) for value in range(20, 241, 20))
N2_FRACTIONS = tuple(value / 10.0 for value in range(1, 10))
TAU_VALUES_S = (0.001, 0.005, 0.010, 0.030, 0.050, 0.100)
STAGES_S = (100, 500, 1000)
EXPECTED_CASES = 9 * 12 * 9 * 6
BASE_ATOL = 1.0
BASE_RTOL = 1.0e-4
FALLBACK_ATOL = 100.0
FALLBACK_RTOL = 1.0e-3


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def case_id(tg_k: float, en_td: float, n2_fraction: float, tau_res_s: float) -> str:
    return (
        f"T{int(round(tg_k))}_EN{int(round(en_td))}_"
        f"N2{int(round(10.0 * n2_fraction)):02d}_TAU{scan.token(1000.0 * tau_res_s)}ms"
    )


def connect_registry(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path, timeout=60.0)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA synchronous=FULL")
    connection.executescript(
        """
        CREATE TABLE IF NOT EXISTS cases (
            case_id TEXT PRIMARY KEY,
            tg_K REAL NOT NULL,
            en_Td REAL NOT NULL,
            n2_fraction REAL NOT NULL,
            h2_fraction REAL NOT NULL,
            tau_res_s REAL NOT NULL,
            stage_s INTEGER NOT NULL,
            status TEXT NOT NULL,
            stage_attempt INTEGER NOT NULL DEFAULT 0,
            total_attempt INTEGER NOT NULL DEFAULT 0,
            last_outcome TEXT,
            output_group TEXT,
            started_at TEXT,
            ended_at TEXT,
            returncode INTEGER,
            timed_out INTEGER,
            t_final_s REAL,
            NH3_cm3 REAL,
            NH3_tail_rel_range REAL,
            NH3_tail_pass INTEGER,
            other_tail_pass INTEGER,
            failure_reason TEXT,
            updated_at TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS events (
            event_id INTEGER PRIMARY KEY AUTOINCREMENT,
            event_time TEXT NOT NULL,
            case_id TEXT NOT NULL,
            event TEXT NOT NULL,
            stage_s INTEGER,
            attempt INTEGER,
            details_json TEXT
        );
        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
        """
    )
    return connection


def add_event(
    connection: sqlite3.Connection,
    cid: str,
    event: str,
    stage_s: int | None,
    attempt: int | None,
    details: dict[str, object] | None = None,
) -> None:
    connection.execute(
        "INSERT INTO events(event_time,case_id,event,stage_s,attempt,details_json) "
        "VALUES(?,?,?,?,?,?)",
        (now(), cid, event, stage_s, attempt, json.dumps(details or {}, ensure_ascii=False)),
    )


def initialize_registry(connection: sqlite3.Connection) -> None:
    timestamp = now()
    for tau_res_s in TAU_VALUES_S:
        for tg_k in TEMPERATURES_K:
            for n2_fraction in N2_FRACTIONS:
                for en_td in EN_VALUES_TD:
                    cid = case_id(tg_k, en_td, n2_fraction, tau_res_s)
                    connection.execute(
                        """
                        INSERT OR IGNORE INTO cases(
                            case_id,tg_K,en_Td,n2_fraction,h2_fraction,tau_res_s,
                            stage_s,status,updated_at
                        ) VALUES(?,?,?,?,?,?,?,?,?)
                        """,
                        (
                            cid,
                            tg_k,
                            en_td,
                            n2_fraction,
                            1.0 - n2_fraction,
                            tau_res_s,
                            STAGES_S[0],
                            "PENDING",
                            timestamp,
                        ),
                    )
    count = connection.execute("SELECT COUNT(*) FROM cases").fetchone()[0]
    if count != EXPECTED_CASES:
        raise RuntimeError(f"案例台账数量错误：{count}，预期 {EXPECTED_CASES}")
    metadata = {
        "expected_cases": EXPECTED_CASES,
        "stages_s": STAGES_S,
        "temperature_K": TEMPERATURES_K,
        "en_Td": EN_VALUES_TD,
        "n2_fraction": N2_FRACTIONS,
        "tau_res_s": TAU_VALUES_S,
        "model": "pure-gas 0D CW-CSTR; surface reactions disabled",
        "created_or_checked_at": timestamp,
    }
    for key, value in metadata.items():
        connection.execute(
            "INSERT OR REPLACE INTO metadata(key,value) VALUES(?,?)",
            (key, json.dumps(value, ensure_ascii=False)),
        )
    connection.commit()


def recover_stale_running(connection: sqlite3.Connection) -> int:
    rows = connection.execute(
        "SELECT case_id,stage_s,total_attempt FROM cases WHERE status='RUNNING'"
    ).fetchall()
    timestamp = now()
    for row in rows:
        connection.execute(
            """
            UPDATE cases SET status='INTERRUPTED',last_outcome='INTERRUPTED',
                failure_reason='scheduler stopped before attempt completed',
                ended_at=?,updated_at=? WHERE case_id=?
            """,
            (timestamp, timestamp, row["case_id"]),
        )
        add_event(
            connection,
            row["case_id"],
            "INTERRUPTED_RECOVERED",
            row["stage_s"],
            row["total_attempt"],
        )
    connection.commit()
    return len(rows)


def import_completed_results(connection: sqlite3.Connection, source: Path) -> int:
    imported = 0
    for params_path in source.rglob("params.json"):
        try:
            record = json.loads(params_path.read_text(encoding="utf-8"))
            if record.get("returncode") != 0 or record.get("timed_out"):
                continue
            cid = case_id(
                float(record["tg_K"]),
                float(record["en_Td"]),
                float(record["n2_fraction"]),
                float(record["tau_res_s"]),
            )
            row = connection.execute(
                "SELECT status,total_attempt FROM cases WHERE case_id=?", (cid,)
            ).fetchone()
            if row is None or row["status"] in (
                "ACCEPTED",
                "ACCEPTED_AT_MAX_TIME_WITH_NH3_FLAG",
            ):
                continue
            t_final = float(record.get("t_final_s", 0.0))
            if t_final + 1.0e-8 < STAGES_S[-1]:
                continue
            nh3_pass = int(float(record.get("NH3_tail_pass", 0.0)) == 1.0)
            other_pass = int(float(record.get("other_tail_pass", 0.0)) == 1.0)
            status = (
                "ACCEPTED"
                if nh3_pass and other_pass
                else "ACCEPTED_AT_MAX_TIME_WITH_NH3_FLAG"
                if other_pass
                else "FAILED_AT_MAX_TIME"
            )
            timestamp = now()
            connection.execute(
                """
                UPDATE cases SET stage_s=?,status=?,last_outcome='IMPORTED',
                    output_group=?,ended_at=?,returncode=0,timed_out=0,t_final_s=?,
                    NH3_cm3=?,NH3_tail_rel_range=?,NH3_tail_pass=?,
                    other_tail_pass=?,failure_reason=NULL,updated_at=?
                WHERE case_id=?
                """,
                (
                    STAGES_S[-1],
                    status,
                    str(params_path.parent),
                    timestamp,
                    t_final,
                    float(record.get("NH3_cm-3", 0.0)),
                    float(record.get("NH3_cm-3_tail_rel_range", 0.0)),
                    nh3_pass,
                    other_pass,
                    timestamp,
                    cid,
                ),
            )
            add_event(
                connection,
                cid,
                "IMPORTED_COMPLETED_RESULT",
                STAGES_S[-1],
                row["total_attempt"],
                {"source": str(params_path), "status": status},
            )
            imported += 1
        except (KeyError, TypeError, ValueError, json.JSONDecodeError):
            continue
    connection.commit()
    return imported


def export_status(connection: sqlite3.Connection, root: Path) -> None:
    rows = connection.execute(
        """
        SELECT case_id,tg_K,en_Td,n2_fraction,h2_fraction,tau_res_s,stage_s,
               status,stage_attempt,total_attempt,last_outcome,output_group,
               started_at,ended_at,returncode,timed_out,t_final_s,NH3_cm3,
               NH3_tail_rel_range,NH3_tail_pass,other_tail_pass,failure_reason,
               updated_at
        FROM cases ORDER BY tau_res_s,tg_K,n2_fraction,en_Td
        """
    ).fetchall()
    csv_path = root / "case_status.csv"
    csv_temp = root / "case_status.csv.tmp"
    with csv_temp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(rows[0].keys())
        writer.writerows(tuple(row) for row in rows)
    os.replace(csv_temp, csv_path)

    counts = {
        row["status"]: row["count"]
        for row in connection.execute(
            "SELECT status,COUNT(*) AS count FROM cases GROUP BY status"
        ).fetchall()
    }
    stage_counts = {
        str(row["stage_s"]): row["count"]
        for row in connection.execute(
            "SELECT stage_s,COUNT(*) AS count FROM cases GROUP BY stage_s"
        ).fetchall()
    }
    progress = {
        "updated_at": now(),
        "expected_cases": EXPECTED_CASES,
        "status_counts": counts,
        "current_stage_counts": stage_counts,
        "terminal_cases": sum(
            counts.get(status, 0)
            for status in (
                "ACCEPTED",
                "ACCEPTED_AT_MAX_TIME_WITH_NH3_FLAG",
                "FAILED",
                "FAILED_AT_MAX_TIME",
            )
        ),
    }
    json_path = root / "adaptive_progress.json"
    json_temp = root / "adaptive_progress.json.tmp"
    json_temp.write_text(json.dumps(progress, ensure_ascii=False, indent=2), encoding="utf-8")
    os.replace(json_temp, json_path)


def acquire_lock(root: Path) -> tuple[Path, int]:
    lock_path = root / "adaptive_scan.lock"
    if lock_path.exists():
        try:
            old_pid = int(lock_path.read_text(encoding="ascii").strip())
            os.kill(old_pid, 0)
        except (OSError, ValueError):
            lock_path.unlink()
        else:
            raise RuntimeError(f"已有自适应扫描调度器在运行，PID={old_pid}")
    descriptor = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(descriptor, str(os.getpid()).encode("ascii"))
    os.close(descriptor)
    return lock_path, os.getpid()


def next_stage(stage_s: int) -> int | None:
    index = STAGES_S.index(stage_s)
    return STAGES_S[index + 1] if index + 1 < len(STAGES_S) else None


def select_batch(
    connection: sqlite3.Connection, workers: int, max_stage_attempts: int
) -> list[sqlite3.Row]:
    return connection.execute(
        """
        SELECT * FROM cases
        WHERE status IN ('PENDING','INTERRUPTED','RETRY')
          AND stage_attempt < ?
        ORDER BY stage_s ASC,tau_res_s DESC,en_Td DESC,tg_K DESC,n2_fraction ASC
        LIMIT ?
        """,
        (max_stage_attempts, workers),
    ).fetchall()


def mark_running(
    connection: sqlite3.Connection, row: sqlite3.Row
) -> tuple[dict[str, object], str]:
    stage_attempt = int(row["stage_attempt"]) + 1
    total_attempt = int(row["total_attempt"]) + 1
    group = f"adaptive_scan/stage_{int(row['stage_s'])}s/attempt_{total_attempt:03d}"
    timestamp = now()
    connection.execute(
        """
        UPDATE cases SET status='RUNNING',stage_attempt=?,total_attempt=?,
            output_group=?,started_at=?,ended_at=NULL,returncode=NULL,
            timed_out=NULL,failure_reason=NULL,updated_at=?
        WHERE case_id=?
        """,
        (
            stage_attempt,
            total_attempt,
            group,
            timestamp,
            timestamp,
            row["case_id"],
        ),
    )
    add_event(
        connection,
        row["case_id"],
        "ATTEMPT_STARTED",
        row["stage_s"],
        total_attempt,
        {"group": group},
    )
    task = dict(row)
    if stage_attempt >= 3:
        solver_policy = "relaxed_third_attempt"
        atol, rtol = FALLBACK_ATOL, FALLBACK_RTOL
    else:
        solver_policy = "baseline"
        atol, rtol = BASE_ATOL, BASE_RTOL
    task.update(
        {
            "stage_attempt": stage_attempt,
            "total_attempt": total_attempt,
            "output_group": group,
            "solver_policy": solver_policy,
            "atol": atol,
            "rtol": rtol,
        }
    )
    connection.execute(
        "UPDATE events SET details_json=? WHERE rowid=last_insert_rowid()",
        (json.dumps({"group": group, "solver_policy": solver_policy,
                     "atol": atol, "rtol": rtol}, ensure_ascii=False),),
    )
    return task, group


def run_task(
    work: Path,
    root: Path,
    task: dict[str, object],
    window_dt_s: float,
    report_dt_s: float,
    timeout_s: float,
    idle_timeout_s: float,
    watchdog_poll_s: float,
) -> dict[str, object]:
    return scan.run_case(
        work,
        root,
        float(task["tg_K"]),
        float(task["en_Td"]),
        float(task["n2_fraction"]),
        window_dt_s,
        float(task["stage_s"]),
        report_dt_s,
        str(task["output_group"]),
        float(task["tau_res_s"]),
        atol=float(task["atol"]),
        rtol=float(task["rtol"]),
        timeout_s=timeout_s,
        idle_timeout_s=idle_timeout_s,
        watchdog_poll_s=watchdog_poll_s,
    )


def apply_result(
    connection: sqlite3.Connection,
    task: dict[str, object],
    result: dict[str, object],
    max_stage_attempts: int,
) -> str:
    cid = str(task["case_id"])
    stage_s = int(task["stage_s"])
    total_attempt = int(task["total_attempt"])
    stage_attempt = int(task["stage_attempt"])
    success = result.get("returncode") == 0 and not bool(result.get("timed_out"))
    nh3_pass = int(float(result.get("NH3_tail_pass", 0.0)) == 1.0) if success else 0
    other_pass = int(float(result.get("other_tail_pass", 0.0)) == 1.0) if success else 0
    following_stage = next_stage(stage_s)

    if success and nh3_pass and other_pass:
        status = "ACCEPTED"
        outcome = "CONVERGED"
        new_stage = stage_s
        new_stage_attempt = stage_attempt
    elif success and following_stage is not None:
        status = "PENDING"
        outcome = "NEEDS_EXTENSION"
        new_stage = following_stage
        new_stage_attempt = 0
    elif success and other_pass:
        status = "ACCEPTED_AT_MAX_TIME_WITH_NH3_FLAG"
        outcome = "MAX_TIME_NH3_FLAG"
        new_stage = stage_s
        new_stage_attempt = stage_attempt
    elif success:
        status = "FAILED_AT_MAX_TIME"
        outcome = "MAX_TIME_OTHER_SPECIES_FAIL"
        new_stage = stage_s
        new_stage_attempt = stage_attempt
    elif stage_attempt < max_stage_attempts:
        status = "RETRY"
        outcome = "NUMERICAL_RETRY"
        new_stage = stage_s
        new_stage_attempt = stage_attempt
    else:
        status = "FAILED"
        outcome = "NUMERICAL_FAILURE"
        new_stage = stage_s
        new_stage_attempt = stage_attempt

    timestamp = now()
    connection.execute(
        """
        UPDATE cases SET stage_s=?,status=?,stage_attempt=?,last_outcome=?,
            ended_at=?,returncode=?,timed_out=?,t_final_s=?,NH3_cm3=?,
            NH3_tail_rel_range=?,NH3_tail_pass=?,other_tail_pass=?,
            failure_reason=?,updated_at=? WHERE case_id=?
        """,
        (
            new_stage,
            status,
            new_stage_attempt,
            outcome,
            timestamp,
            result.get("returncode"),
            int(bool(result.get("timed_out"))),
            result.get("t_final_s"),
            result.get("NH3_cm-3"),
            result.get("NH3_cm-3_tail_rel_range"),
            nh3_pass,
            other_pass,
            result.get("failure_reason"),
            timestamp,
            cid,
        ),
    )
    add_event(
        connection,
        cid,
        outcome,
        stage_s,
        total_attempt,
        {
            "status": status,
            "next_stage_s": new_stage,
            "returncode": result.get("returncode"),
            "timed_out": bool(result.get("timed_out")),
            "NH3_tail_rel_range": result.get("NH3_cm-3_tail_rel_range"),
            "other_tail_pass": other_pass,
        },
    )
    connection.commit()
    return outcome


def run_scheduler(
    connection: sqlite3.Connection,
    work: Path,
    root: Path,
    workers: int,
    window_dt_s: float,
    report_dt_s: float,
    timeout_s: float,
    idle_timeout_s: float,
    watchdog_poll_s: float,
    max_stage_attempts: int,
    snapshot_every: int,
    max_cases: int | None,
) -> None:
    attempts_finished = 0
    attempts_started = 0
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {}

        def submit_next() -> bool:
            nonlocal attempts_started
            if max_cases is not None and attempts_started >= max_cases:
                return False
            rows = select_batch(connection, 1, max_stage_attempts)
            if not rows:
                return False
            task, _ = mark_running(connection, rows[0])
            connection.commit()
            future = executor.submit(
                run_task,
                work,
                root,
                task,
                window_dt_s,
                report_dt_s,
                timeout_s,
                idle_timeout_s,
                watchdog_poll_s,
            )
            futures[future] = task
            attempts_started += 1
            return True

        while len(futures) < workers and submit_next():
            pass
        export_status(connection, root)

        while futures:
            completed, _ = wait(tuple(futures), return_when=FIRST_COMPLETED)
            for future in completed:
                task = futures.pop(future)
                try:
                    result = future.result()
                except Exception as exc:  # registry must survive unexpected task errors
                    result = {
                        "returncode": None,
                        "timed_out": False,
                        "failure_reason": f"{type(exc).__name__}: {exc}",
                    }
                outcome = apply_result(
                    connection, task, result, max_stage_attempts
                )
                attempts_finished += 1
                print(
                    f"[{attempts_finished}] {task['case_id']} "
                    f"stage={task['stage_s']}s outcome={outcome}",
                    flush=True,
                )
                submit_next()
                if attempts_finished % snapshot_every == 0:
                    export_status(connection, root)
        export_status(connection, root)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Recoverable adaptive-time pure-gas CW-CSTR scan"
    )
    parser.add_argument("--output-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--window-dt-s", type=float, default=10.0)
    parser.add_argument("--report-dt-s", type=float, default=10.0)
    parser.add_argument("--mxstep", type=int, default=50000)
    parser.add_argument("--timeout-s", type=float, default=21600.0)
    parser.add_argument("--idle-timeout-s", type=float, default=120.0)
    parser.add_argument("--watchdog-poll-s", type=float, default=5.0)
    parser.add_argument("--max-stage-attempts", type=int, default=3)
    parser.add_argument("--snapshot-every", type=int, default=32)
    parser.add_argument("--max-cases", type=int)
    parser.add_argument("--initialize-only", action="store_true")
    parser.add_argument("--import-results", type=Path, action="append", default=[])
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if (
        args.workers <= 0
        or args.window_dt_s <= 0
        or args.report_dt_s <= 0
        or args.timeout_s <= 0
        or args.idle_timeout_s <= 0
        or args.watchdog_poll_s <= 0
        or args.max_stage_attempts <= 0
        or args.snapshot_every <= 0
        or (args.max_cases is not None and args.max_cases <= 0)
    ):
        raise SystemExit("所有时间、计数和并行参数必须为正数")

    root = args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    lock_path, _ = acquire_lock(root)
    connection = connect_registry(root / "case_registry.sqlite")
    try:
        initialize_registry(connection)
        recovered = recover_stale_running(connection)
        imported = sum(
            import_completed_results(connection, source.resolve())
            for source in args.import_results
        )
        export_status(connection, root)
        print(
            json.dumps(
                {
                    "registry_cases": EXPECTED_CASES,
                    "recovered_stale_running": recovered,
                    "imported_results": imported,
                },
                ensure_ascii=False,
            ),
            flush=True,
        )
        if args.initialize_only:
            return 0
        work = scan.prepare_worktree(root, args.mxstep)
        run_scheduler(
            connection,
            work,
            root,
            args.workers,
            args.window_dt_s,
            args.report_dt_s,
            args.timeout_s,
            args.idle_timeout_s,
            args.watchdog_poll_s,
            args.max_stage_attempts,
            args.snapshot_every,
            args.max_cases,
        )
        export_status(connection, root)
        return 0
    finally:
        connection.close()
        try:
            lock_path.unlink()
        except FileNotFoundError:
            pass


if __name__ == "__main__":
    raise SystemExit(main())
