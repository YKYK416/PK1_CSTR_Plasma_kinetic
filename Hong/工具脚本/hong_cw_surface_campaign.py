#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recoverable, resource-controlled campaign runner for closed-0D Hong scans.

This module deliberately contains no plasma chemistry.  It prepares and
controls calls to :mod:`hong_cw_surface_longtime_scan`, which remains the
single source of truth for the Hong 2017/2018 corrected mechanism and the
closed-0D CW physical model.  The campaign ledger makes long scans restartable
after an interruption without overwriting partial attempts.

The default production plan contains 9 temperatures x 9 N2 fractions x 12
reduced fields.  The 60 Td column is scheduled in an isolated one-worker phase
because its DVODE trajectory is materially stiffer than the other columns.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import queue
import sqlite3
import sys
import time
from collections import Counter, deque
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Any, Iterator

import hong_cw_surface_longtime_scan as model


SCRIPT_PATH = Path(__file__).resolve()
TOOL_DIR = SCRIPT_PATH.parent
DEFAULT_ROOT = model.PROJECT / "analysis" / "campaign_hong_closed0d_100s"
DB_NAME = "campaign.sqlite"
MANIFEST_NAME = "campaign_manifest.json"
OVERVIEW_NAME = "cases_overview.csv"
EVENTS_NAME = "progress.jsonl"
LOCK_NAME = ".campaign.lock"
LIVE_LEDGER_NAME = "live_ledger.md"
LIVE_CSV_NAME = "live_ledger.csv"
CAMPAIGN_VERSION = 1
HEARTBEAT_SECONDS = 5.0
NORMAL_WORKERS = 2
HIGH_STIFFNESS_80_WORKERS = 2
HIGH_STIFFNESS_FIELDS_TD = (60.0, 80.0)

TERMINAL_STATES = {"completed", "failed", "invalid_output"}
RUNNABLE_STATES = {"planned", "interrupted", "failed", "invalid_output"}


def utc_now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    os.replace(temporary, path)


def safe_json(value: str | None) -> dict[str, Any]:
    if not value:
        return {}
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        return {}


def connect(root: Path) -> sqlite3.Connection:
    database = root / DB_NAME
    conn = sqlite3.connect(database, timeout=30.0)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA journal_mode = WAL")
    return conn


def initialise_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS campaign (
          campaign_id TEXT PRIMARY KEY,
          created_at TEXT NOT NULL,
          manifest_json TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS cases (
          case_id TEXT PRIMARY KEY,
          temperature_K REAL NOT NULL,
          en_Td REAL NOT NULL,
          n2_fraction REAL NOT NULL,
          h2_fraction REAL NOT NULL,
          phase TEXT NOT NULL,
          priority INTEGER NOT NULL,
          solver_json TEXT NOT NULL,
          status TEXT NOT NULL,
          attempt_count INTEGER NOT NULL DEFAULT 0,
          active_attempt INTEGER,
          last_error TEXT,
          created_at TEXT NOT NULL,
          updated_at TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS attempts (
          attempt_id INTEGER PRIMARY KEY AUTOINCREMENT,
          case_id TEXT NOT NULL REFERENCES cases(case_id),
          attempt_number INTEGER NOT NULL,
          status TEXT NOT NULL,
          output_dir TEXT NOT NULL,
          pid INTEGER,
          started_at TEXT,
          heartbeat_at TEXT,
          ended_at TEXT,
          last_reported_time_s REAL,
          return_code INTEGER,
          timed_out INTEGER,
          validation_json TEXT,
          error_text TEXT,
          UNIQUE(case_id, attempt_number)
        );
        CREATE TABLE IF NOT EXISTS events (
          event_id INTEGER PRIMARY KEY AUTOINCREMENT,
          created_at TEXT NOT NULL,
          case_id TEXT,
          attempt_number INTEGER,
          event_type TEXT NOT NULL,
          payload_json TEXT NOT NULL
        );
        CREATE INDEX IF NOT EXISTS idx_cases_queue ON cases(status, phase, priority, temperature_K, n2_fraction);
        CREATE INDEX IF NOT EXISTS idx_attempts_case ON attempts(case_id, attempt_number);
        """
    )
    conn.commit()


def case_id(temperature_K: float, en_Td: float, n2_fraction: float) -> str:
    return f"T{temperature_K:03.0f}_EN{en_Td:03.0f}_N2{n2_fraction:.1f}".replace(".", "p")


def solver_profile(en_Td: float) -> dict[str, Any]:
    """Numerical controls established by the prior convergence investigation."""
    # 80 Td was initially placed in the normal queue.  Production evidence at
    # 300--325 K showed repeated inability to leave t=0 under that branch.
    # A direct 100 s rescue validation for the hardest failed point
    # (300 K, x_N2=0.4) subsequently passed with ATOL=1e9 and a 0.1 s outer
    # advance.  The change is numerical only: it does not alter kinetics,
    # geometry, gas composition, or the closed-0D boundary condition.
    if en_Td == 80.0:
        return {
            "outer_target_dt_s": 0.1,
            "finite_exposure_time_s": 100.0,
            "report_interval_s": 10.0,
            "ATOL_cm-3": 1.0e9,
            "RTOL": 1.0e-4,
            "DVODE_HMAX_s": 5.0e-3,
            "MXSTEP": 100000,
            "early_diagnostics": False,
            "timeout_s": 900.0,
            "branch": "EoverN_80Td_rescue_ATOL1e9_dt0p1s",
        }

    # The pre-identified 60 Td column retains its established conservative
    # controls.  Other fields use the ordinary production profile.
    high_stiffness = en_Td in HIGH_STIFFNESS_FIELDS_TD
    return {
        "outer_target_dt_s": 1.0,
        "finite_exposure_time_s": 100.0,
        "report_interval_s": 10.0,
        "ATOL_cm-3": 1.0e8 if high_stiffness else 1.0,
        "RTOL": 1.0e-4,
        "DVODE_HMAX_s": 5.0e-3,
        "MXSTEP": 100000 if high_stiffness else 5000,
        "early_diagnostics": False,
        "timeout_s": 3600.0,
        "branch": f"EoverN_{en_Td:g}Td_high_stiffness" if high_stiffness else "default",
    }


def campaign_manifest() -> dict[str, Any]:
    mechanism = model.BASE / model.KINETIC_INPUT
    return {
        "campaign_version": CAMPAIGN_VERSION,
        "created_at": utc_now(),
        "physical_model": {
            "surface_mechanism": "Hong 2017/2018 corrected, Table-5 metal column",
            "reactor_boundary": "closed 0D; no CSTR/feed/outlet/residence-time term",
            "forcing": "continuous E/N; no pulse waveform or cycle count",
            "finite_exposure_interpretation": "NH3 inventory and terminal net rate at t_obs=100 s; not an all-species steady-state claim.",
        },
        "grid": {
            "temperature_K": list(range(300, 501, 25)),
            "n2_fraction": [index / 10 for index in range(1, 10)],
            "h2_fraction": [1.0 - index / 10 for index in range(1, 10)],
            "EoverN_Td": list(range(20, 241, 20)),
            "case_count": 972,
        },
        "scheduling": {
            "phase_order": ["high_stiffness_60Td", "high_stiffness_80Td", "normal"],
            "high_stiffness_workers": 1,
            "high_stiffness_80Td_workers": HIGH_STIFFNESS_80_WORKERS,
            "normal_workers": NORMAL_WORKERS,
            "high_stiffness_isolation": "60 Td is isolated to one worker; verified 80 Td reruns use two isolated workers and are never co-scheduled with normal cases.",
            "cpu_affinity": "Disjoint logical-core groups are requested per concurrently active solver; a warning is logged if unavailable.",
        },
        "provenance": {
            "campaign_script": str(SCRIPT_PATH),
            "campaign_script_sha256": sha256_file(SCRIPT_PATH),
            "driver_script": str(Path(model.__file__).resolve()),
            "driver_script_sha256": sha256_file(Path(model.__file__).resolve()),
            "kinetic_input": str(mechanism),
            "kinetic_input_sha256": sha256_file(mechanism),
        },
    }


def append_event(conn: sqlite3.Connection, root: Path, event_type: str, payload: dict[str, Any],
                 case: sqlite3.Row | dict[str, Any] | None = None, attempt_number: int | None = None) -> None:
    case_identifier = None if case is None else str(case["case_id"])
    created_at = utc_now()
    serialized = json.dumps(payload, ensure_ascii=False, sort_keys=True)
    conn.execute("INSERT INTO events(created_at, case_id, attempt_number, event_type, payload_json) VALUES(?,?,?,?,?)",
                 (created_at, case_identifier, attempt_number, event_type, serialized))
    with (root / EVENTS_NAME).open("a", encoding="utf-8") as handle:
        handle.write(json.dumps({"created_at": created_at, "case_id": case_identifier,
                                 "attempt_number": attempt_number, "event_type": event_type,
                                 "payload": payload}, ensure_ascii=False) + "\n")


def export_overview(conn: sqlite3.Connection, root: Path, filename: str = OVERVIEW_NAME) -> None:
    rows = conn.execute(
        """SELECT c.case_id, c.temperature_K, c.en_Td, c.n2_fraction, c.h2_fraction,
                  c.phase, c.status, c.attempt_count, c.active_attempt, c.last_error,
                  a.output_dir, a.pid, a.started_at, a.heartbeat_at, a.ended_at,
                  a.last_reported_time_s, a.return_code, a.timed_out
           FROM cases c LEFT JOIN attempts a
             ON a.case_id=c.case_id AND a.attempt_number=c.active_attempt
           ORDER BY c.priority, c.temperature_K, c.n2_fraction, c.en_Td"""
    ).fetchall()
    fields = [column[0] for column in conn.execute(
        """SELECT c.case_id, c.temperature_K, c.en_Td, c.n2_fraction, c.h2_fraction,
                  c.phase, c.status, c.attempt_count, c.active_attempt, c.last_error,
                  a.output_dir, a.pid, a.started_at, a.heartbeat_at, a.ended_at,
                  a.last_reported_time_s, a.return_code, a.timed_out
           FROM cases c LEFT JOIN attempts a
             ON a.case_id=c.case_id AND a.attempt_number=c.active_attempt LIMIT 0"""
    ).description]
    temporary = root / f".{filename}.{os.getpid()}.tmp"
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(dict(row) for row in rows)
    os.replace(temporary, root / filename)


def write_live_ledger(conn: sqlite3.Connection, root: Path) -> dict[str, Any]:
    """Write a human-readable live status file from the SQLite source of truth.

    This is intentionally read-only with respect to campaign execution.  Files
    are atomically replaced, so a viewer never sees a partially written ledger.
    """
    payload = status_payload(conn)
    phase_rows = conn.execute(
        "SELECT phase, status, COUNT(*) AS count FROM cases GROUP BY phase, status ORDER BY phase, status"
    ).fetchall()
    active_rows = conn.execute(
        """SELECT c.case_id, c.temperature_K, c.en_Td, c.n2_fraction, c.h2_fraction, c.phase,
                  c.attempt_count, a.pid, a.started_at, a.heartbeat_at, a.last_reported_time_s, a.output_dir
           FROM cases c JOIN attempts a ON a.case_id=c.case_id AND a.attempt_number=c.active_attempt
           WHERE c.status='running' ORDER BY c.updated_at"""
    ).fetchall()
    issue_rows = conn.execute(
        """SELECT c.case_id, c.temperature_K, c.en_Td, c.n2_fraction, c.phase, c.status,
                  c.attempt_count, c.last_error, a.output_dir
           FROM cases c LEFT JOIN attempts a ON a.case_id=c.case_id AND a.attempt_number=c.active_attempt
           WHERE c.status IN ('failed','invalid_output','interrupted')
           ORDER BY c.updated_at DESC LIMIT 20"""
    ).fetchall()
    event = conn.execute("SELECT created_at, event_type, payload_json FROM events ORDER BY event_id DESC LIMIT 1").fetchone()
    manifest = safe_json((root / MANIFEST_NAME).read_text(encoding="utf-8"))
    lines = [
        "# Hong closed-0D 972-case live ledger",
        "",
        f"- **Refreshed:** {utc_now()}",
        "- **Model:** Hong 2017/2018 corrected mechanism, Table-5 metal column; closed 0D; CW E/N; no CSTR/feed/outlet/residence-time; no pulse.",
        f"- **Campaign directory:** `{root}`",
        f"- **Grid:** {manifest.get('grid', {}).get('case_count', '?')} cases (9 T × 9 xN2 × 12 E/N).",
        "",
        "## Overall status",
        "",
        "| planned | running | completed | failed | invalid output | interrupted | total |",
        "|---:|---:|---:|---:|---:|---:|---:|",
        "| {planned} | {running} | {completed} | {failed} | {invalid_output} | {interrupted} | {total} |".format(
            planned=payload["counts"].get("planned", 0), running=payload["counts"].get("running", 0),
            completed=payload["counts"].get("completed", 0), failed=payload["counts"].get("failed", 0),
            invalid_output=payload["counts"].get("invalid_output", 0),
            interrupted=payload["counts"].get("interrupted", 0), total=payload["total"]),
        "",
        "## Queue breakdown",
        "",
        "| phase | status | cases |",
        "|---|---|---:|",
    ]
    lines.extend(f"| {row['phase']} | {row['status']} | {row['count']} |" for row in phase_rows)
    lines.extend(["", "## Active case(s)", ""])
    if active_rows:
        lines.extend(["| case | Tgas (K) | E/N (Td) | xN2 | xH2 | PID | latest physical time (s) | heartbeat | attempt |",
                      "|---|---:|---:|---:|---:|---:|---:|---|---:|"])
        lines.extend(
            f"| {row['case_id']} | {row['temperature_K']:g} | {row['en_Td']:g} | {row['n2_fraction']:.1f} | "
            f"{row['h2_fraction']:.1f} | {row['pid'] or ''} | {row['last_reported_time_s'] or 0:g} | "
            f"{row['heartbeat_at'] or ''} | {row['attempt_count']} |" for row in active_rows)
    else:
        lines.append("No active solver case is recorded.")
    lines.extend(["", "## Cases needing review or later resume", ""])
    if issue_rows:
        lines.extend(["| case | phase | status | attempt | reason | output directory |",
                      "|---|---|---|---:|---|---|"])
        for row in issue_rows:
            reason = (row["last_error"] or "").replace("|", "\\|").replace("\n", " ")
            lines.append(f"| {row['case_id']} | {row['phase']} | {row['status']} | {row['attempt_count']} | {reason} | `{row['output_dir'] or ''}` |")
    else:
        lines.append("None.")
    if event is not None:
        lines.extend(["", "## Latest ledger event", "", f"`{event['created_at']}` — `{event['event_type']}`", ""])
    atomic_json(root / "live_ledger.json", {"refreshed_at": utc_now(), "overall": payload,
                                               "active_cases": [dict(row) for row in active_rows],
                                               "issues": [dict(row) for row in issue_rows]})
    temporary = root / f".{LIVE_LEDGER_NAME}.{os.getpid()}.tmp"
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    os.replace(temporary, root / LIVE_LEDGER_NAME)
    export_overview(conn, root, LIVE_CSV_NAME)
    return payload


def watch_ledger(root: Path, interval_s: float) -> dict[str, Any]:
    """Refresh live files until the campaign controller releases its lock."""
    if interval_s < 1.0:
        raise ValueError("台账刷新间隔不得小于 1 s")
    conn = open_existing(root)
    try:
        while True:
            result = write_live_ledger(conn, root)
            lock_path = root / LOCK_NAME
            lock = safe_json(lock_path.read_text(encoding="utf-8")) if lock_path.is_file() else {}
            if not lock_path.is_file() or not pid_alive(lock.get("pid")):
                result["watcher_state"] = "campaign_not_active_final_snapshot_written"
                return result
            time.sleep(interval_s)
    finally:
        conn.close()


def initialise_campaign(root: Path) -> None:
    root = root.resolve()
    root.mkdir(parents=True, exist_ok=False)
    manifest = campaign_manifest()
    atomic_json(root / MANIFEST_NAME, manifest)
    (root / "cases").mkdir()
    conn = connect(root)
    try:
        initialise_schema(conn)
        campaign_id = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode("utf-8")).hexdigest()[:16]
        conn.execute("INSERT INTO campaign(campaign_id, created_at, manifest_json) VALUES(?,?,?)",
                     (campaign_id, utc_now(), json.dumps(manifest, ensure_ascii=False, sort_keys=True)))
        cases: list[tuple[Any, ...]] = []
        for temperature_K in range(300, 501, 25):
            for n2_index in range(1, 10):
                n2_fraction = n2_index / 10
                for en_Td in range(20, 241, 20):
                    profile = solver_profile(float(en_Td))
                    phase = ("high_stiffness_60Td" if en_Td == 60 else
                             "high_stiffness_80Td" if en_Td == 80 else "normal")
                    priority = 0 if phase == "high_stiffness_60Td" else 1 if phase == "high_stiffness_80Td" else 2
                    cases.append((case_id(temperature_K, en_Td, n2_fraction), temperature_K, en_Td,
                                  n2_fraction, 1.0 - n2_fraction, phase, priority,
                                  json.dumps(profile, ensure_ascii=False, sort_keys=True), "planned", utc_now(), utc_now()))
        conn.executemany(
            """INSERT INTO cases(case_id,temperature_K,en_Td,n2_fraction,h2_fraction,phase,priority,
                                  solver_json,status,created_at,updated_at)
               VALUES(?,?,?,?,?,?,?,?,?,?,?)""", cases)
        append_event(conn, root, "campaign_initialized", {"case_count": len(cases), "campaign_id": campaign_id})
        conn.commit()
        export_overview(conn, root)
    finally:
        conn.close()


def open_existing(root: Path) -> sqlite3.Connection:
    if not (root / MANIFEST_NAME).is_file() or not (root / DB_NAME).is_file():
        raise RuntimeError(f"不是有效的 campaign 目录：{root}")
    conn = connect(root)
    initialise_schema(conn)
    return conn


def status_payload(conn: sqlite3.Connection) -> dict[str, Any]:
    counts = Counter({row["status"]: row["count"] for row in conn.execute(
        "SELECT status, COUNT(*) AS count FROM cases GROUP BY status")})
    return {"total": sum(counts.values()), "counts": dict(sorted(counts.items()))}


def write_state(attempt_dir: Path, payload: dict[str, Any]) -> None:
    atomic_json(attempt_dir / "state.json", payload)


def pid_alive(pid: int | None) -> bool:
    if pid is None or pid <= 0:
        return False
    try:
        import psutil
        return psutil.pid_exists(pid)
    except Exception:  # noqa: BLE001
        if os.name == "nt":
            return False
        try:
            os.kill(pid, 0)
        except OSError:
            return False
        return True


def recover_orphaned(conn: sqlite3.Connection, root: Path) -> int:
    """Convert stale running rows to interrupted while preserving all artifacts."""
    rows = conn.execute(
        """SELECT c.*, a.output_dir, a.pid, a.attempt_number FROM cases c
           JOIN attempts a ON a.case_id=c.case_id AND a.attempt_number=c.active_attempt
           WHERE c.status='running'"""
    ).fetchall()
    recovered = 0
    for row in rows:
        if pid_alive(row["pid"]):
            raise RuntimeError(f"case {row['case_id']} 仍由 PID {row['pid']} 运行；拒绝启动第二个调度器")
        now = utc_now()
        conn.execute("UPDATE cases SET status='interrupted', last_error=?, updated_at=? WHERE case_id=?",
                     ("controller restart found no live recorded solver PID", now, row["case_id"]))
        conn.execute("UPDATE attempts SET status='interrupted', ended_at=?, error_text=? WHERE case_id=? AND attempt_number=?",
                     (now, "controller restart found no live recorded solver PID", row["case_id"], row["attempt_number"]))
        write_state(Path(row["output_dir"]), {"state": "interrupted", "case_id": row["case_id"],
                                               "attempt_number": row["attempt_number"], "ended_at": now,
                                               "reason": "controller restart found no live recorded solver PID"})
        append_event(conn, root, "attempt_interrupted", {"reason": "no live recorded solver PID"}, row,
                     int(row["attempt_number"]))
        recovered += 1
    conn.commit()
    return recovered


def reclassify_unfinished_80td(conn: sqlite3.Connection, root: Path) -> dict[str, Any]:
    """Move unresolved 80 Td cases to the conservative one-worker branch.

    Completed cases are deliberately left untouched: their result and exact
    numerical provenance remain immutable.  Planned, failed, invalid, and
    interrupted cases receive the high-stiffness solver profile for their next
    attempt.  The migration is recorded as a campaign event and a standalone
    JSON note so the original manifest remains historically accurate.
    """
    rows = conn.execute(
        """SELECT * FROM cases WHERE en_Td=80.0 AND status != 'completed'
           ORDER BY temperature_K, n2_fraction"""
    ).fetchall()
    now = utc_now()
    profile = json.dumps(solver_profile(80.0), ensure_ascii=False, sort_keys=True)
    for row in rows:
        conn.execute("UPDATE cases SET phase='high_stiffness_80Td', priority=1, solver_json=?, updated_at=? WHERE case_id=?",
                     (profile, now, row["case_id"]))
        append_event(conn, root, "numerical_branch_reclassified", {
            "from_phase": row["phase"], "to_phase": "high_stiffness_80Td",
            "reason": "80 Td normal branch repeatedly remained at t=0 or failed at 300--325 K",
            "solver": safe_json(profile),
        }, row, int(row["active_attempt"]) if row["active_attempt"] is not None else None)
    update = {
        "recorded_at": now,
        "scope": "unfinished 80 Td cases only; completed results unchanged",
        "case_count": len(rows),
        "new_phase": "high_stiffness_80Td",
        "solver": safe_json(profile),
        "reason": ("80 Td normal-branch stalls were followed by a passing direct 100 s "
                   "rescue validation at 300 K, x_N2=0.4 using ATOL=1e9 and a 0.1 s outer advance."),
        "campaign_script_sha256": sha256_file(SCRIPT_PATH),
        "driver_script_sha256": sha256_file(Path(model.__file__).resolve()),
        "kinetic_input_sha256": sha256_file(model.BASE / model.KINETIC_INPUT),
    }
    atomic_json(root / "numerical_strategy_update_80Td_rescue_ATOL1e9_dt0p1s.json", update)
    conn.commit()
    export_overview(conn, root)
    write_live_ledger(conn, root)
    return {"reclassified": len(rows), **status_payload(conn)}


def requeue_80td_timeouts(conn: sqlite3.Connection, root: Path) -> dict[str, Any]:
    """Requeue only 80 Td timeouts after the stdout-streaming repair.

    Completed and site-balance-review cases remain immutable.  The original
    80 Td production profile (ATOL=1e9, 0.1 s outer advance) is retained, so
    this is an orchestration repair rather than a physical-model change.
    """
    rows = conn.execute(
        """SELECT * FROM cases WHERE en_Td=80.0 AND status='failed'
           AND last_error='timeout' ORDER BY temperature_K, n2_fraction"""
    ).fetchall()
    now = utc_now()
    phase = "high_stiffness_80Td_timeout_rescue"
    for row in rows:
        conn.execute("UPDATE cases SET phase=?, priority=1, status='planned', last_error=NULL, updated_at=? WHERE case_id=?",
                     (phase, now, row["case_id"]))
        append_event(conn, root, "timeout_rescue_requeued", {
            "from_phase": row["phase"], "to_phase": phase,
            "reason": "prior timeout attributed to undrained stdout pipe; rerun uses streamed per-case console.log",
            "solver": safe_json(row["solver_json"]),
        }, row, int(row["active_attempt"]) if row["active_attempt"] is not None else None)
    update = {"recorded_at": now, "case_count": len(rows), "scope": "80 Td timeout cases only",
              "phase": phase, "numerical_profile": solver_profile(80.0),
              "orchestration_change": "stream solver stdout directly to each attempt console.log; avoid PIPE backpressure",
              "campaign_script_sha256": sha256_file(SCRIPT_PATH),
              "driver_script_sha256": sha256_file(Path(model.__file__).resolve())}
    atomic_json(root / "stdout_streaming_timeout_rescue_80Td.json", update)
    conn.commit()
    export_overview(conn, root)
    write_live_ledger(conn, root)
    return {"requeued": len(rows), **status_payload(conn)}


def normal_failure_rescue_profile(en_Td: float) -> dict[str, Any]:
    """Strict controls validated on the 400 K, 100 Td, x_N2=0.8 failure."""
    return {
        "outer_target_dt_s": 0.1,
        "finite_exposure_time_s": 100.0,
        "report_interval_s": 10.0,
        "ATOL_cm-3": 1.0e8,
        "RTOL": 1.0e-4,
        "DVODE_HMAX_s": 5.0e-3,
        "MXSTEP": 100000,
        "early_diagnostics": False,
        "timeout_s": 1800.0,
        "branch": f"EoverN_{en_Td:g}Td_normal_failure_rescue_ATOL1e8_dt0p1s",
    }


def requeue_normal_failures_strict(conn: sqlite3.Connection, root: Path) -> dict[str, Any]:
    """Requeue only the legacy normal-phase failures with validated strict controls."""
    rows = conn.execute(
        """SELECT * FROM cases WHERE phase='normal' AND status='failed'
           ORDER BY temperature_K, en_Td, n2_fraction"""
    ).fetchall()
    now, phase = utc_now(), "normal_failure_rescue_strict"
    for row in rows:
        profile = normal_failure_rescue_profile(float(row["en_Td"]))
        conn.execute("UPDATE cases SET phase=?, priority=2, status='planned', solver_json=?, last_error=NULL, updated_at=? WHERE case_id=?",
                     (phase, json.dumps(profile, ensure_ascii=False, sort_keys=True), now, row["case_id"]))
        append_event(conn, root, "normal_failure_rescue_requeued", {
            "from_phase": row["phase"], "to_phase": phase,
            "reason": "strict profile validated at 400 K, 100 Td, x_N2=0.8 after legacy DVODE failures",
            "solver": profile,
        }, row, int(row["active_attempt"]) if row["active_attempt"] is not None else None)
    update = {"recorded_at": now, "case_count": len(rows), "scope": "legacy normal-phase failed cases only",
              "phase": phase, "profile_template": normal_failure_rescue_profile(100.0),
              "validation_anchor": "400 K, 100 Td, x_N2=0.8; 100 s completed with site-balance error 1.39929104e-10",
              "campaign_script_sha256": sha256_file(SCRIPT_PATH),
              "driver_script_sha256": sha256_file(Path(model.__file__).resolve())}
    atomic_json(root / "normal_failure_rescue_strict.json", update)
    conn.commit()
    export_overview(conn, root)
    write_live_ledger(conn, root)
    return {"requeued": len(rows), **status_payload(conn)}


def final_recovery_profile(en_Td: float) -> dict[str, Any]:
    """Field-aware strict controls for the remaining non-accepted cases."""
    if en_Td == 60.0:
        return {
            "outer_target_dt_s": 0.1, "finite_exposure_time_s": 100.0, "report_interval_s": 10.0,
            "ATOL_cm-3": 1.0e7, "RTOL": 1.0e-4, "DVODE_HMAX_s": 2.5e-3, "MXSTEP": 200000,
            "early_diagnostics": False, "timeout_s": 1800.0,
            "branch": "EoverN_60Td_final_recovery_ATOL1e7_hmax0p0025s",
        }
    if en_Td == 80.0:
        return {
            "outer_target_dt_s": 0.1, "finite_exposure_time_s": 100.0, "report_interval_s": 10.0,
            "ATOL_cm-3": 1.0e8, "RTOL": 1.0e-4, "DVODE_HMAX_s": 5.0e-3, "MXSTEP": 100000,
            "early_diagnostics": False, "timeout_s": 1800.0,
            "branch": "EoverN_80Td_final_recovery_ATOL1e8_dt0p1s",
        }
    return {
        "outer_target_dt_s": 0.1, "finite_exposure_time_s": 100.0, "report_interval_s": 10.0,
        "ATOL_cm-3": 1.0e7, "RTOL": 1.0e-4, "DVODE_HMAX_s": 5.0e-3, "MXSTEP": 200000,
        "early_diagnostics": False, "timeout_s": 1800.0,
        "branch": f"EoverN_{en_Td:g}Td_final_recovery_ATOL1e7_dt0p1s",
    }


def requeue_remaining_final_recovery(conn: sqlite3.Connection, root: Path) -> dict[str, Any]:
    """Requeue every currently non-accepted case without changing prior attempts."""
    rows = conn.execute(
        """SELECT * FROM cases WHERE status IN ('failed', 'invalid_output', 'interrupted')
           ORDER BY CASE WHEN en_Td=80.0 THEN 0 WHEN en_Td=60.0 THEN 1 ELSE 2 END,
                    temperature_K, en_Td, n2_fraction"""
    ).fetchall()
    now, phase = utc_now(), "final_strict_recovery_2workers"
    for row in rows:
        profile = final_recovery_profile(float(row["en_Td"]))
        conn.execute("UPDATE cases SET phase=?, priority=0, status='planned', solver_json=?, last_error=NULL, updated_at=? WHERE case_id=?",
                     (phase, json.dumps(profile, ensure_ascii=False, sort_keys=True), now, row["case_id"]))
        append_event(conn, root, "final_recovery_requeued", {
            "from_phase": row["phase"], "to_phase": phase, "from_status": row["status"],
            "reason": "user-requested rerun of all non-accepted cases using field-aware strict controls and streamed stdout",
            "solver": profile,
        }, row, int(row["active_attempt"]) if row["active_attempt"] is not None else None)
    update = {
        "recorded_at": now, "case_count": len(rows), "scope": "all cases not currently accepted",
        "phase": phase, "workers": 2, "runtime": "runtime_build_streamed_stdout_v3",
        "profiles": {"60_Td": final_recovery_profile(60.0), "80_Td": final_recovery_profile(80.0),
                     "other_fields": final_recovery_profile(100.0)},
        "campaign_script_sha256": sha256_file(SCRIPT_PATH),
        "driver_script_sha256": sha256_file(Path(model.__file__).resolve()),
    }
    atomic_json(root / "final_strict_recovery_2workers.json", update)
    conn.commit()
    export_overview(conn, root)
    write_live_ledger(conn, root)
    return {"requeued": len(rows), **status_payload(conn)}


def requeue_80td_conservation_final(conn: sqlite3.Connection, root: Path) -> dict[str, Any]:
    """Requeue only final 80 Td site-balance review cases with the passing profile."""
    rows = conn.execute(
        """SELECT * FROM cases WHERE en_Td=80.0 AND status='invalid_output'
           ORDER BY temperature_K, n2_fraction"""
    ).fetchall()
    now, phase = utc_now(), "final_80Td_conservation_recovery"
    profile = {
        "outer_target_dt_s": 0.1, "finite_exposure_time_s": 100.0, "report_interval_s": 10.0,
        "ATOL_cm-3": 1.0e7, "RTOL": 1.0e-4, "DVODE_HMAX_s": 2.5e-3, "MXSTEP": 200000,
        "early_diagnostics": False, "timeout_s": 1800.0,
        "branch": "EoverN_80Td_final_conservation_ATOL1e7_hmax0p0025s",
    }
    for row in rows:
        conn.execute("UPDATE cases SET phase=?, priority=0, status='planned', solver_json=?, last_error=NULL, updated_at=? WHERE case_id=?",
                     (phase, json.dumps(profile, ensure_ascii=False, sort_keys=True), now, row["case_id"]))
        append_event(conn, root, "80td_conservation_final_requeued", {
            "from_phase": row["phase"], "to_phase": phase,
            "reason": "validated at 450 K, x_N2=0.2 with site-balance error 5.1744e-14",
            "solver": profile,
        }, row, int(row["active_attempt"]) if row["active_attempt"] is not None else None)
    atomic_json(root / "final_80Td_conservation_recovery.json", {
        "recorded_at": now, "case_count": len(rows), "phase": phase, "workers": 1,
        "validation_anchor": "450 K, 80 Td, x_N2=0.2; 100 s completed; site-balance error 5.1744e-14",
        "solver": profile, "campaign_script_sha256": sha256_file(SCRIPT_PATH),
        "driver_script_sha256": sha256_file(Path(model.__file__).resolve()),
    })
    conn.commit()
    export_overview(conn, root)
    write_live_ledger(conn, root)
    return {"requeued": len(rows), **status_payload(conn)}


def site_constrained_recovery_profile(en_Td: float) -> dict[str, Any]:
    """Use a report-boundary projection only for numerically audited fields."""
    profile = dict(solver_profile(en_Td) if en_Td == 80.0 else normal_failure_rescue_profile(en_Td))
    profile["site_projection"] = True
    profile["branch"] = f"{profile['branch']}_site_projection_report_boundary"
    return profile


def requeue_site_constrained_recovery(conn: sqlite3.Connection, root: Path) -> dict[str, Any]:
    """Recover the audited 80/100 Td tail cases without changing prior attempts."""
    recovered = recover_orphaned(conn, root)
    rows = conn.execute(
        """SELECT * FROM cases
           WHERE en_Td IN (80.0, 100.0) AND status IN ('failed', 'invalid_output', 'interrupted')
           ORDER BY CASE WHEN en_Td=80.0 THEN 0 ELSE 1 END, temperature_K, n2_fraction"""
    ).fetchall()
    now, phase = utc_now(), "site_constrained_recovery"
    for row in rows:
        profile = site_constrained_recovery_profile(float(row["en_Td"]))
        conn.execute("UPDATE cases SET phase=?, priority=0, status='planned', solver_json=?, last_error=NULL, updated_at=? WHERE case_id=?",
                     (phase, json.dumps(profile, ensure_ascii=False, sort_keys=True), now, row["case_id"]))
        append_event(conn, root, "site_constrained_recovery_requeued", {
            "from_phase": row["phase"], "to_phase": phase, "from_status": row["status"],
            "reason": ("report-boundary projection of the conserved five-species surface-site inventory; "
                       "validated at 475 K, 80 Td, x_N2=0.7 and 475 K, 100 Td, x_N2=0.8"),
            "solver": profile,
        }, row, int(row["active_attempt"]) if row["active_attempt"] is not None else None)
    update = {
        "recorded_at": now, "recovered_orphaned_attempts": recovered, "case_count": len(rows),
        "scope": "only unfinished E/N=80 or 100 Td cases; prior attempts remain immutable",
        "phase": phase, "workers": NORMAL_WORKERS, "runtime": "runtime_build_site_projection_v4",
        "profiles": {"80_Td": site_constrained_recovery_profile(80.0),
                     "100_Td": site_constrained_recovery_profile(100.0)},
        "validation_anchors": [
            "475 K, 80 Td, x_N2=0.7: completed 100 s; raw maximum site deviation 4.35985e-08; projected CSV maximum 2.24e-16",
            "475 K, 100 Td, x_N2=0.8: completed 100 s; raw maximum site deviation 1.11786e-08; projected CSV maximum 1.12e-16; all 12 reported states agree with the unprojected strict reference at CSV precision",
        ],
        "campaign_script_sha256": sha256_file(SCRIPT_PATH),
        "driver_script_sha256": sha256_file(Path(model.__file__).resolve()),
    }
    atomic_json(root / "site_constrained_recovery.json", update)
    conn.commit()
    export_overview(conn, root)
    write_live_ledger(conn, root)
    return {"requeued": len(rows), "recovered_orphaned": recovered, **status_payload(conn)}


@contextmanager
def campaign_lock(root: Path) -> Iterator[None]:
    lock_path = root / LOCK_NAME
    data = {"pid": os.getpid(), "created_at": utc_now(), "script": str(SCRIPT_PATH)}
    try:
        descriptor = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        existing = safe_json(lock_path.read_text(encoding="utf-8"))
        if pid_alive(existing.get("pid")):
            raise RuntimeError(f"已有活动 campaign 调度器：PID {existing.get('pid')}，锁文件 {lock_path}")
        stale_path = root / f"{LOCK_NAME}.stale_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        os.replace(lock_path, stale_path)
        descriptor = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(data, handle, ensure_ascii=False, indent=2)
        yield
    finally:
        try:
            lock_path.unlink()
        except FileNotFoundError:
            pass


def assign_affinity(pid: int, cpu_group: list[int]) -> dict[str, Any]:
    """Best-effort affinity/priority control; no solver result depends on it."""
    result: dict[str, Any] = {"requested_cpu_group": cpu_group, "applied": False}
    try:
        import psutil
        process = psutil.Process(pid)
        process.cpu_affinity(cpu_group)
        if hasattr(psutil, "BELOW_NORMAL_PRIORITY_CLASS"):
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
        result["applied"] = True
    except Exception as exc:  # noqa: BLE001
        result["warning"] = f"could not apply affinity/priority: {exc}"
    return result


def cpu_groups(worker_count: int) -> list[list[int]]:
    logical = max(os.cpu_count() or 1, 1)
    if worker_count == 1:
        return [list(range(min(4, logical)))]
    width = max(1, min(4, logical // worker_count))
    groups = [list(range(index * width, min((index + 1) * width, logical))) for index in range(worker_count)]
    return [group or [index % logical] for index, group in enumerate(groups)]


def most_recent_time(output_dir: Path) -> float | None:
    files = list(output_dir.glob("cw_surface_*.csv"))
    if not files:
        return None
    try:
        with files[0].open("rb") as handle:
            handle.seek(max(0, files[0].stat().st_size - 65536))
            lines = handle.read().decode("utf-8", errors="replace").splitlines()
        for line in reversed(lines):
            if line and not line.startswith("time_s"):
                return float(line.split(",", 1)[0])
    except (OSError, ValueError):
        return None
    return None


def create_attempt(conn: sqlite3.Connection, root: Path, case: sqlite3.Row) -> tuple[int, Path]:
    now = utc_now()
    number = int(case["attempt_count"]) + 1
    output_dir = root / "cases" / str(case["case_id"]) / f"attempt_{number:03d}"
    output_dir.mkdir(parents=True, exist_ok=False)
    conn.execute("UPDATE cases SET status='running', attempt_count=?, active_attempt=?, last_error=NULL, updated_at=? WHERE case_id=?",
                 (number, number, now, case["case_id"]))
    conn.execute(
        """INSERT INTO attempts(case_id,attempt_number,status,output_dir,started_at,heartbeat_at)
           VALUES(?,?,?,?,?,?)""", (case["case_id"], number, "running", str(output_dir), now, now))
    write_state(output_dir, {"state": "launching", "case_id": case["case_id"], "attempt_number": number,
                             "started_at": now, "physical_parameters": {"temperature_K": case["temperature_K"],
                             "en_Td": case["en_Td"], "n2_fraction": case["n2_fraction"], "h2_fraction": case["h2_fraction"]},
                             "solver": safe_json(case["solver_json"])})
    append_event(conn, root, "attempt_launched", {"output_dir": str(output_dir)}, case, number)
    conn.commit()
    return number, output_dir


def update_event(conn: sqlite3.Connection, root: Path, event: dict[str, Any]) -> None:
    case_identifier, number = str(event["case_id"]), int(event["attempt_number"])
    now = event.get("timestamp", utc_now())
    if event["event"] == "started":
        conn.execute("UPDATE attempts SET pid=?, heartbeat_at=? WHERE case_id=? AND attempt_number=?",
                     (int(event["pid"]), now, case_identifier, number))
    elif event["event"] == "heartbeat":
        conn.execute("UPDATE attempts SET heartbeat_at=?, last_reported_time_s=? WHERE case_id=? AND attempt_number=?",
                     (now, event.get("last_reported_time_s"), case_identifier, number))
    conn.execute("UPDATE cases SET updated_at=? WHERE case_id=?", (now, case_identifier))
    append_event(conn, root, f"solver_{event['event']}", event,
                 {"case_id": case_identifier}, number)
    conn.commit()


def validate_record(output_dir: Path, record: dict[str, Any]) -> tuple[bool, dict[str, Any], str | None]:
    if record.get("returncode") != 0 or record.get("timed_out"):
        return False, {}, str(record.get("failure_reason", "solver did not exit successfully"))
    # New attempts carry evidence parsed from the Fortran stdout. Historical
    # attempts predate this field, so their classification is reported by the
    # dedicated audit without silently changing their existing ledger status.
    completion = record.get("solver_completion")
    if completion is not None and not bool(completion.get("pass")):
        return False, {"solver_completion": completion}, "missing or inconsistent solver completion marker"
    tag = str(record["tag"])
    series = output_dir / f"cw_surface_{tag}.csv"
    if not series.is_file():
        return False, {}, "missing concentration CSV"
    try:
        rows = model.read_series(series)
        final = rows[-1]
        target = float(record["t_end_requested_s"])
        if abs(final["time_s"] - target) > max(1.0e-8, target * 1.0e-6):
            return False, {}, f"final physical time {final['time_s']:g} s does not reach {target:g} s"
        site_error = max(abs(row["site_balance_error"]) for row in rows)
    except Exception as exc:  # noqa: BLE001
        return False, {}, f"could not validate concentration CSV: {exc}"
    projected = bool(record.get("site_projection", False))
    raw_projection_error = completion.get("site_projection_max_rel") if isinstance(completion, dict) else None
    validation = {"final_time_s": final["time_s"], "site_balance_error_max": site_error,
                  "site_balance_limit": model.SITE_BALANCE_LIMIT,
                  "site_projection": projected,
                  "preprojection_site_balance_error_max": raw_projection_error,
                  "row_count": len(rows)}
    if site_error > model.SITE_BALANCE_LIMIT:
        return False, validation, f"site balance error {site_error:g} exceeds limit"
    return True, validation, None


def finish_attempt(conn: sqlite3.Connection, root: Path, case: sqlite3.Row, number: int,
                   output_dir: Path, record: dict[str, Any] | None, error: str | None) -> None:
    now = utc_now()
    if error is not None:
        valid, validation, reason = False, {}, error
    else:
        assert record is not None
        valid, validation, reason = validate_record(output_dir, record)
    status = "completed" if valid else ("invalid_output" if record and record.get("returncode") == 0 else "failed")
    write_state(output_dir, {"state": status, "case_id": case["case_id"], "attempt_number": number,
                             "ended_at": now, "record": record, "validation": validation, "reason": reason})
    conn.execute("UPDATE cases SET status=?, last_error=?, updated_at=? WHERE case_id=?",
                 (status, reason, now, case["case_id"]))
    conn.execute(
        """UPDATE attempts SET status=?, ended_at=?, return_code=?, timed_out=?, validation_json=?, error_text=?
           WHERE case_id=? AND attempt_number=?""",
        (status, now, None if record is None else record.get("returncode"),
         None if record is None else int(bool(record.get("timed_out"))),
         json.dumps(validation, ensure_ascii=False, sort_keys=True), reason, case["case_id"], number))
    append_event(conn, root, f"attempt_{status}", {"validation": validation, "reason": reason}, case, number)
    conn.commit()


def ordered_cases(conn: sqlite3.Connection, phase: str) -> deque[sqlite3.Row]:
    placeholders = ",".join("?" for _ in RUNNABLE_STATES)
    return deque(conn.execute(
        f"""SELECT * FROM cases WHERE phase=? AND status IN ({placeholders})
             ORDER BY temperature_K, n2_fraction, en_Td""",
        (phase, *sorted(RUNNABLE_STATES))).fetchall())


def run_phase(conn: sqlite3.Connection, root: Path, work: Path, phase: str, workers: int) -> None:
    pending = ordered_cases(conn, phase)
    if not pending:
        return
    event_queue: queue.Queue[dict[str, Any]] = queue.Queue()
    active: dict[Any, tuple[sqlite3.Row, int, Path]] = {}
    affinity_groups = cpu_groups(workers)

    def submit(executor: ThreadPoolExecutor, case: sqlite3.Row, lane: int) -> None:
        number, output_dir = create_attempt(conn, root, case)
        profile = safe_json(case["solver_json"])
        affinity = affinity_groups[lane]

        def execute() -> dict[str, Any]:
            def callback(event: dict[str, Any]) -> None:
                payload = dict(event)
                payload["case_id"] = case["case_id"]
                payload["attempt_number"] = number
                if payload["event"] == "started":
                    payload["cpu_control"] = assign_affinity(int(payload["pid"]), affinity)
                payload["last_reported_time_s"] = most_recent_time(output_dir)
                write_state(output_dir, {"state": payload["event"], **payload})
                event_queue.put(payload)

            return model.run_case(work, root, float(case["temperature_K"]), float(case["en_Td"]),
                                  float(case["n2_fraction"]), float(profile["outer_target_dt_s"]),
                                  float(profile["finite_exposure_time_s"]), float(profile["report_interval_s"]),
                                  "campaign", atol=float(profile["ATOL_cm-3"]),
                                  rtol=float(profile["RTOL"]), internal_hmax_s=float(profile["DVODE_HMAX_s"]),
                                  early_diagnostics=bool(profile["early_diagnostics"]),
                                  internal_mxstep=int(profile["MXSTEP"]), timeout_s=float(profile["timeout_s"]),
                                  site_projection=bool(profile.get("site_projection", False)),
                                  case_output_dir=output_dir, progress_callback=callback)

        future = executor.submit(execute)
        active[future] = (case, number, output_dir)

    with ThreadPoolExecutor(max_workers=workers, thread_name_prefix=f"hong-{phase}") as executor:
        for lane in range(min(workers, len(pending))):
            submit(executor, pending.popleft(), lane)
        while active:
            completed, _ = wait(active, timeout=HEARTBEAT_SECONDS, return_when=FIRST_COMPLETED)
            while not event_queue.empty():
                update_event(conn, root, event_queue.get())
            for future in completed:
                case, number, output_dir = active.pop(future)
                try:
                    finish_attempt(conn, root, case, number, output_dir, future.result(), None)
                except Exception as exc:  # noqa: BLE001
                    finish_attempt(conn, root, case, number, output_dir, None, f"orchestrator exception: {exc}")
                if pending:
                    submit(executor, pending.popleft(), len(active) % workers)
    while not event_queue.empty():
        update_event(conn, root, event_queue.get())


def run_campaign(root: Path, phase: str) -> None:
    conn = open_existing(root)
    try:
        with campaign_lock(root):
            recovered = recover_orphaned(conn, root)
            append_event(conn, root, "campaign_run_started", {"phase": phase, "recovered_orphaned": recovered})
            conn.commit()
            # Build only when the user explicitly invokes run/resume, never for init/status/self-test.
            # Completion-marker validation was added after the first campaign
            # runtime had been built.  A new, adjacent worktree preserves all
            # old attempts while ensuring future reruns compile the verified
            # driver rather than silently reusing the older executable.
            runtime_root = root / "runtime_build_site_projection_v4"
            work = model.prepare_worktree(runtime_root)
            phases = ["high_stiffness_60Td", "high_stiffness_80Td", "normal"] if phase == "all" else [phase]
            for current_phase in phases:
                workers = (1 if current_phase == "final_80Td_conservation_recovery"
                           else HIGH_STIFFNESS_80_WORKERS if current_phase.startswith("high_stiffness_80Td")
                           else 1 if current_phase.startswith("high_stiffness_") else NORMAL_WORKERS)
                run_phase(conn, root, work, current_phase, workers)
                export_overview(conn, root)
            append_event(conn, root, "campaign_run_finished", {"phase": phase, **status_payload(conn)})
            conn.commit()
            export_overview(conn, root)
    finally:
        conn.close()


def validate_campaign(root: Path) -> dict[str, Any]:
    conn = open_existing(root)
    try:
        invalidated = 0
        rows = conn.execute(
            """SELECT c.*, a.output_dir, a.attempt_number FROM cases c
               JOIN attempts a ON a.case_id=c.case_id AND a.attempt_number=c.active_attempt
               WHERE c.status='completed'"""
        ).fetchall()
        for row in rows:
            params_path = Path(row["output_dir"]) / "params.json"
            record = safe_json(params_path.read_text(encoding="utf-8")) if params_path.is_file() else None
            valid, validation, reason = validate_record(Path(row["output_dir"]), record or {})
            if not valid:
                conn.execute("UPDATE cases SET status='invalid_output', last_error=?, updated_at=? WHERE case_id=?",
                             (reason, utc_now(), row["case_id"]))
                conn.execute("UPDATE attempts SET status='invalid_output', validation_json=?, error_text=? WHERE case_id=? AND attempt_number=?",
                             (json.dumps(validation, ensure_ascii=False), reason, row["case_id"], row["attempt_number"]))
                append_event(conn, root, "attempt_invalidated", {"reason": reason, "validation": validation}, row,
                             int(row["attempt_number"]))
                invalidated += 1
        conn.commit()
        export_overview(conn, root)
        return {"invalidated": invalidated, **status_payload(conn)}
    finally:
        conn.close()


def audit_80td(root: Path) -> dict[str, Any]:
    """Write a non-mutating audit for the latest attempt of every 80 Td case."""
    conn = open_existing(root)
    try:
        rows = conn.execute(
            """SELECT c.case_id, c.status, c.temperature_K, c.n2_fraction,
                      a.attempt_number, a.output_dir
               FROM cases c JOIN attempts a
                 ON a.case_id=c.case_id AND a.attempt_number=c.active_attempt
               WHERE c.en_Td=80.0 ORDER BY c.temperature_K, c.n2_fraction"""
        ).fetchall()
        entries: list[dict[str, Any]] = []
        for row in rows:
            output_dir = Path(row["output_dir"])
            params_path, console_path = output_dir / "params.json", output_dir / "console.log"
            record = safe_json(params_path.read_text(encoding="utf-8")) if params_path.is_file() else {}
            tag = str(record.get("tag", ""))
            target = float(record.get("t_end_requested_s", 100.0))
            console = console_path.read_text(encoding="utf-8", errors="replace") if console_path.is_file() else ""
            if tag:
                evidence = model.solver_completion_evidence(console, tag, target)
            else:
                evidence = {"done_marker_present": False, "reported_t_end_s": None,
                            "requested_t_end_s": target, "endpoint_tolerance_s": None,
                            "pass": False, "dvode_step_size_warning_count": 0,
                            "dvode_error_marker_present": False}
            if bool(evidence["pass"]):
                classification = "solver_completion_verified"
            elif record.get("timed_out"):
                classification = "solver_timeout_no_completion_marker"
            elif record.get("returncode") not in (0, None):
                classification = "solver_nonzero_exit_no_completion_marker"
            else:
                classification = "unverified_no_completion_marker"
            entries.append({
                "case_id": row["case_id"], "ledger_status": row["status"],
                "temperature_K": row["temperature_K"], "n2_fraction": row["n2_fraction"],
                "attempt_number": row["attempt_number"], "output_dir": str(output_dir),
                "returncode": record.get("returncode"), "timed_out": record.get("timed_out"),
                "csv_final_time_s": record.get("t_final_s"), "classification": classification,
                **evidence,
            })
        counts = Counter(item["classification"] for item in entries)
        report = {
            "generated_at": utc_now(), "scope": "latest active attempt for all E/N=80 Td cases",
            "mutation": "none; SQLite statuses and existing outputs are not modified",
            "case_count": len(entries), "classification_counts": dict(sorted(counts.items())),
            "entries": entries, "campaign_script_sha256": sha256_file(SCRIPT_PATH),
            "driver_script_sha256": sha256_file(Path(model.__file__).resolve()),
        }
        atomic_json(root / "audit_80Td_solver_completion.json", report)
        fields = ["case_id", "ledger_status", "temperature_K", "n2_fraction", "attempt_number", "classification",
                  "returncode", "timed_out", "csv_final_time_s", "done_marker_present", "reported_t_end_s",
                  "requested_t_end_s", "endpoint_tolerance_s", "pass", "dvode_step_size_warning_count",
                  "dvode_error_marker_present", "site_projection_max_rel", "output_dir"]
        with (root / "audit_80Td_solver_completion.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(entries)
        lines = ["# 80 Td solver-completion audit", "",
                 "This report does not modify case states or previous outputs.", "",
                 "| classification | cases |", "| --- | ---: |"]
        lines.extend(f"| {name} | {count} |" for name, count in sorted(counts.items()))
        lines.extend(["", "A DVODE `T + H = T` warning is diagnostic metadata, not a failure by itself. A verified run must emit both `DONE tag=...` and a matching `t_end_s`."])
        (root / "audit_80Td_solver_completion.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
        return {"case_count": len(entries), "classification_counts": dict(sorted(counts.items())),
                "report": str(root / "audit_80Td_solver_completion.json")}
    finally:
        conn.close()


def run_self_test(root: Path) -> dict[str, Any]:
    """Exercise ledger creation and interruption recovery without a solver call."""
    initialise_campaign(root)
    conn = open_existing(root)
    try:
        total = conn.execute("SELECT COUNT(*) FROM cases").fetchone()[0]
        distinct = conn.execute("SELECT COUNT(DISTINCT case_id) FROM cases").fetchone()[0]
        high60 = conn.execute("SELECT COUNT(*) FROM cases WHERE phase='high_stiffness_60Td'").fetchone()[0]
        high80 = conn.execute("SELECT COUNT(*) FROM cases WHERE phase='high_stiffness_80Td'").fetchone()[0]
        normal = conn.execute("SELECT COUNT(*) FROM cases WHERE phase='normal'").fetchone()[0]
        if (total, distinct, high60, high80, normal) != (972, 972, 81, 81, 810):
            raise AssertionError(f"unexpected grid partition: {(total, distinct, high60, high80, normal)}")
        selected = conn.execute("SELECT * FROM cases ORDER BY priority, temperature_K, n2_fraction, en_Td LIMIT 1").fetchone()
        attempt_number, output_dir = create_attempt(conn, root, selected)
        # An impossible PID represents a controller that died after marking a
        # case running.  No executable is launched in this test.
        conn.execute("UPDATE attempts SET pid=? WHERE case_id=? AND attempt_number=?",
                     (999999999, selected["case_id"], attempt_number))
        conn.commit()
        recovered = recover_orphaned(conn, root)
        after = conn.execute("SELECT status, attempt_count FROM cases WHERE case_id=?", (selected["case_id"],)).fetchone()
        report = {"solver_called": False, "grid_cases": total, "unique_case_ids": distinct,
                  "high_stiffness_60Td_cases": high60, "high_stiffness_80Td_cases": high80,
                  "normal_cases": normal,
                  "recovered_interrupted_cases": recovered, "selected_case": selected["case_id"],
                  "selected_case_status": after["status"], "selected_case_attempt_count": after["attempt_count"],
                  "passed": recovered == 1 and after["status"] == "interrupted" and after["attempt_count"] == 1}
        atomic_json(root / "self_test_report.json", report)
        append_event(conn, root, "self_test_finished", report)
        conn.commit()
        export_overview(conn, root)
        if not report["passed"]:
            raise AssertionError(report)
        return report
    finally:
        conn.close()


def main() -> int:
    parser = argparse.ArgumentParser(description="Recoverable Hong closed-0D finite-exposure campaign controller")
    parser.add_argument("command", choices=("init", "status", "run", "resume", "validate", "export", "watch", "reclassify-80td", "requeue-80td-timeouts", "requeue-normal-failures-strict", "requeue-final-recovery", "requeue-80td-conservation-final", "requeue-site-constrained-recovery", "audit-80td", "self-test"))
    parser.add_argument("--output-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--phase", choices=("all", "high_stiffness_60Td", "high_stiffness_80Td", "high_stiffness_80Td_timeout_rescue", "normal", "normal_failure_rescue_strict", "final_strict_recovery_2workers", "final_80Td_conservation_recovery", "site_constrained_recovery"), default="all")
    parser.add_argument("--interval-s", type=float, default=5.0,
                        help="watch 模式刷新实时台账的间隔（秒，默认 5）")
    args = parser.parse_args()
    root = args.output_root.resolve()
    if args.command == "init":
        initialise_campaign(root)
        result: dict[str, Any] = {"created": str(root), "solver_called": False}
    elif args.command == "self-test":
        result = run_self_test(root)
    elif args.command == "status":
        conn = open_existing(root)
        try:
            result = status_payload(conn)
        finally:
            conn.close()
    elif args.command in {"run", "resume"}:
        run_campaign(root, args.phase)
        conn = open_existing(root)
        try:
            result = status_payload(conn)
        finally:
            conn.close()
    elif args.command == "validate":
        result = validate_campaign(root)
    elif args.command == "audit-80td":
        result = audit_80td(root)
    elif args.command == "watch":
        result = watch_ledger(root, args.interval_s)
    elif args.command == "reclassify-80td":
        conn = open_existing(root)
        try:
            result = reclassify_unfinished_80td(conn, root)
        finally:
            conn.close()
    elif args.command == "requeue-80td-timeouts":
        conn = open_existing(root)
        try:
            result = requeue_80td_timeouts(conn, root)
        finally:
            conn.close()
    elif args.command == "requeue-normal-failures-strict":
        conn = open_existing(root)
        try:
            result = requeue_normal_failures_strict(conn, root)
        finally:
            conn.close()
    elif args.command == "requeue-final-recovery":
        conn = open_existing(root)
        try:
            result = requeue_remaining_final_recovery(conn, root)
        finally:
            conn.close()
    elif args.command == "requeue-80td-conservation-final":
        conn = open_existing(root)
        try:
            result = requeue_80td_conservation_final(conn, root)
        finally:
            conn.close()
    elif args.command == "requeue-site-constrained-recovery":
        conn = open_existing(root)
        try:
            result = requeue_site_constrained_recovery(conn, root)
        finally:
            conn.close()
    else:  # export
        conn = open_existing(root)
        try:
            export_overview(conn, root)
            result = {"exported": str(root / OVERVIEW_NAME), **status_payload(conn)}
        finally:
            conn.close()
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
