#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recoverable 972-case Hong gas--surface CSTR campaign at tau=10 ms.

The ledger engine is reused from the closed-0D campaign, while the physical
driver is replaced by :mod:`hong_cw_surface_cstr`.  Existing campaigns are not
modified.  The default 100 s horizon is a stationarity test horizon, not a
closed-system finite-exposure interpretation.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import hong_cw_surface_campaign as ledger
import hong_cw_surface_cstr as model


SCRIPT_PATH = Path(__file__).resolve()
DEFAULT_ROOT = model.PROJECT / "analysis" / "campaign_hong_surface_cstr_tau10ms_20260915"
TAU_RES_S = model.TAU_RES_S
EXTENSION_HORIZON_S = 200.0
EXTENSION_PHASE = "stationarity_extension_200s"
STATIONARITY_PREFIX = "terminal CSTR stationarity gate failed"


def solver_profile(en_td: float) -> dict[str, Any]:
    """One conservative baseline used for all fields before field-specific evidence exists."""
    return {
        "outer_target_dt_s": 1.0e-2,
        "finite_exposure_time_s": 100.0,
        "report_interval_s": 1.0,
        "ATOL_cm-3": 1.0e8,
        "RTOL": 1.0e-4,
        "DVODE_HMAX_s": 1.0e-3,
        "MXSTEP": 200000,
        "early_diagnostics": True,
        "timeout_s": 3600.0,
        "residence_time_s": TAU_RES_S,
        "branch": f"surface_cstr_tau10ms_baseline_EN{en_td:g}Td",
    }


def campaign_manifest() -> dict[str, Any]:
    mechanism = model.BASE / model.KINETIC_INPUT
    return {
        "campaign_version": 1,
        "created_at": ledger.utc_now(),
        "physical_model": {
            "surface_mechanism": "Hong 2017/2018 corrected, Table-5 metal column; 37 explicit heterogeneous reactions",
            "reactor_boundary": "0D CSTR; gas heavy species use (C_in-C)/tau; electrons and SURF states are not flowed",
            "residence_time_s": TAU_RES_S,
            "feed": "N2/H2 only at the scanned inlet composition; reactive heavy species have zero inlet density",
            "forcing": "continuous E/N; no pulse waveform",
            "interpretation": "terminal stationarity is assessed over the last 20 s of a 100 s integration; not a closed-0D inventory.",
        },
        "grid": {"temperature_K": list(range(300, 501, 25)), "n2_fraction": [i / 10 for i in range(1, 10)],
                 "h2_fraction": [1.0 - i / 10 for i in range(1, 10)], "EoverN_Td": list(range(20, 241, 20)),
                 "case_count": 972},
        "scheduling": {"workers": 2, "cpu_policy": "two concurrent solver processes with disjoint affinity groups"},
        "provenance": {"campaign_script": str(SCRIPT_PATH), "driver_module": str(Path(model.__file__).resolve()),
                       "kinetic_input": str(mechanism), "kinetic_input_sha256": ledger.sha256_file(mechanism)},
    }


_base_validate_record = ledger.validate_record
_base_live_ledger = ledger.write_live_ledger


def initialise_campaign(root: Path) -> None:
    """Initialize a CSTR ledger while retaining the pre-written handoff note."""
    root = root.resolve()
    if (root / ledger.MANIFEST_NAME).is_file() and (root / ledger.DB_NAME).is_file():
        return
    handoff = "交接文档_Hong气表面CSTR_tau10ms.md"
    if root.exists():
        unexpected = [item.name for item in root.iterdir() if item.name != handoff]
        if unexpected:
            raise RuntimeError(f"CSTR campaign root contains unexpected pre-initialization files: {unexpected}")
    else:
        root.mkdir(parents=True)
    manifest = campaign_manifest()
    ledger.atomic_json(root / ledger.MANIFEST_NAME, manifest)
    (root / "cases").mkdir()
    conn = ledger.connect(root)
    try:
        ledger.initialise_schema(conn)
        campaign_id = ledger.hashlib.sha256(json.dumps(manifest, sort_keys=True).encode("utf-8")).hexdigest()[:16]
        conn.execute("INSERT INTO campaign(campaign_id, created_at, manifest_json) VALUES(?,?,?)",
                     (campaign_id, ledger.utc_now(), json.dumps(manifest, ensure_ascii=False, sort_keys=True)))
        cases: list[tuple[Any, ...]] = []
        for temperature_k in range(300, 501, 25):
            for n2_index in range(1, 10):
                n2_fraction = n2_index / 10
                for en_td in range(20, 241, 20):
                    profile = solver_profile(float(en_td))
                    phase = "high_stiffness_60Td" if en_td == 60 else "high_stiffness_80Td" if en_td == 80 else "normal"
                    priority = 0 if phase == "high_stiffness_60Td" else 1 if phase == "high_stiffness_80Td" else 2
                    cases.append((ledger.case_id(temperature_k, en_td, n2_fraction), temperature_k, en_td,
                                  n2_fraction, 1.0 - n2_fraction, phase, priority,
                                  json.dumps(profile, ensure_ascii=False, sort_keys=True), "planned",
                                  ledger.utc_now(), ledger.utc_now()))
        conn.executemany(
            """INSERT INTO cases(case_id,temperature_K,en_Td,n2_fraction,h2_fraction,phase,priority,
                                  solver_json,status,created_at,updated_at)
               VALUES(?,?,?,?,?,?,?,?,?,?,?)""", cases)
        ledger.append_event(conn, root, "campaign_initialized", {"case_count": len(cases), "campaign_id": campaign_id,
                                                                    "residence_time_s": TAU_RES_S})
        conn.commit()
        ledger.export_overview(conn, root)
    finally:
        conn.close()


def validate_record(output_dir: Path, record: dict[str, Any]) -> tuple[bool, dict[str, Any], str | None]:
    valid, validation, reason = _base_validate_record(output_dir, record)
    if not valid:
        return valid, validation, reason
    stationarity_keys = ("NH3_tail_pass", "gas_tail_pass", "coverage_tail_pass")
    stationarity = {key: bool(record.get(key, False)) for key in stationarity_keys}
    validation.update({"residence_time_s": record.get("residence_time_s"), "stationarity": stationarity,
                       "terminal_window_s": record.get("terminal_window_s")})
    if not all(stationarity.values()):
        failed = ", ".join(key for key, passed in stationarity.items() if not passed)
        return False, validation, f"terminal CSTR stationarity gate failed: {failed}"
    return True, validation, None


def write_live_ledger(conn, root: Path) -> dict[str, Any]:
    payload = _base_live_ledger(conn, root)
    path = root / ledger.LIVE_LEDGER_NAME
    text = path.read_text(encoding="utf-8")
    text = text.replace("# Hong closed-0D 972-case live ledger", "# Hong gas--surface CSTR (tau=10 ms) 972-case live ledger")
    text = text.replace("Hong 2017/2018 corrected mechanism, Table-5 metal column; closed 0D; CW E/N; no CSTR/feed/outlet/residence-time; no pulse.",
                        "Hong 2017/2018 corrected mechanism, Table-5 metal column; gas--surface 0D CSTR; tau=10 ms; CW E/N; N2/H2 feed; no pulse.")
    path.write_text(text, encoding="utf-8")
    return payload


def record_cases_needing_extension(root: Path, *, label: str, freeze: bool) -> list[dict[str, Any]]:
    """Export cases whose numerical solve completed but terminal stationarity failed.

    A prior 100 s attempt remains immutable.  At the final 200 s stage,
    ``freeze=True`` moves unresolved cases out of the runnable queue while
    retaining their final invalid attempt and its validation evidence.
    """
    conn = ledger.open_existing(root)
    try:
        rows = conn.execute(
            """SELECT c.case_id,c.temperature_K,c.en_Td,c.n2_fraction,c.h2_fraction,c.phase,
                      c.status,c.active_attempt,c.last_error,a.output_dir,a.validation_json
               FROM cases c LEFT JOIN attempts a ON a.case_id=c.case_id AND a.attempt_number=c.active_attempt
               WHERE c.status='invalid_output' AND c.last_error LIKE ?
               ORDER BY c.temperature_K,c.en_Td,c.n2_fraction""",
            (STATIONARITY_PREFIX + "%",),
        ).fetchall()
        records: list[dict[str, Any]] = []
        for row in rows:
            validation = ledger.safe_json(row["validation_json"])
            records.append({
                "case_id": row["case_id"], "temperature_K": row["temperature_K"], "en_Td": row["en_Td"],
                "x_N2": row["n2_fraction"], "x_H2": row["h2_fraction"], "source_phase": row["phase"],
                "source_attempt": row["active_attempt"], "source_status": row["status"],
                "reason": row["last_error"], "terminal_window_s": validation.get("terminal_window_s"),
                "stationarity": json.dumps(validation.get("stationarity", {}), ensure_ascii=False, sort_keys=True),
                "output_dir": row["output_dir"],
            })
        csv_path = root / f"needs_extension_after_{label}.csv"
        fields = list(records[0]) if records else ["case_id", "temperature_K", "en_Td", "x_N2", "x_H2", "source_phase", "source_attempt", "source_status", "reason", "terminal_window_s", "stationarity", "output_dir"]
        with csv_path.open("w", encoding="utf-8", newline="") as handle:
            writer = ledger.csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(records)
        md_path = root / f"needs_extension_after_{label}.md"
        lines = [f"# CSTR stationarity extension register after {label}", "",
                 f"- Generated: {ledger.utc_now()}", f"- Residence time: `{TAU_RES_S:g} s`", f"- Cases: {len(records)}", "",
                 "Each listed case reached its requested physical endpoint and preserved surface sites, but did not satisfy at least one terminal stationarity gate. Its source attempt is retained without overwrite."]
        if freeze:
            lines.extend(["", "These cases have been set to `needs_extension`; they are intentionally excluded from ordinary resume queues until a further horizon is selected."])
        md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        if freeze and rows:
            now = ledger.utc_now()
            for row in rows:
                conn.execute("UPDATE cases SET status='needs_extension', updated_at=? WHERE case_id=?", (now, row["case_id"]))
                ledger.append_event(conn, root, "stationarity_extension_required", {"after_horizon_s": label}, row,
                                    int(row["active_attempt"]) if row["active_attempt"] is not None else None)
            conn.commit()
            ledger.export_overview(conn, root)
        return records
    finally:
        conn.close()


def schedule_200s_extensions(root: Path) -> int:
    """Promote only 100 s stationarity failures to an isolated 200 s phase."""
    conn = ledger.open_existing(root)
    try:
        rows = conn.execute(
            """SELECT * FROM cases WHERE status='invalid_output' AND last_error LIKE ?
               ORDER BY temperature_K,en_Td,n2_fraction""", (STATIONARITY_PREFIX + "%",)
        ).fetchall()
        now = ledger.utc_now()
        for row in rows:
            # Use the validated uniform profile rather than copying a legacy
            # 100 s attempt that may predate the startup-tolerance refresh.
            profile = solver_profile(float(row["en_Td"]))
            profile.update({"finite_exposure_time_s": EXTENSION_HORIZON_S,
                            "branch": str(profile.get("branch", "surface_cstr")) + "_stationarity_extension_200s"})
            conn.execute("UPDATE cases SET phase=?,priority=?,status='planned',solver_json=?,last_error=NULL,updated_at=? WHERE case_id=?",
                         (EXTENSION_PHASE, 3, json.dumps(profile, ensure_ascii=False, sort_keys=True), now, row["case_id"]))
            ledger.append_event(conn, root, "stationarity_extension_scheduled", {
                "from_horizon_s": 100.0, "to_horizon_s": EXTENSION_HORIZON_S,
                "source_attempt": row["active_attempt"], "solver": profile,
            }, row, int(row["active_attempt"]) if row["active_attempt"] is not None else None)
        conn.commit()
        ledger.export_overview(conn, root)
        return len(rows)
    finally:
        conn.close()


def refresh_baseline_profiles(root: Path) -> dict[str, Any]:
    """Apply the validated CSTR startup profile to every unfinished stage-1 case.

    Previous attempts are intentionally retained.  Only their successor case
    state and solver JSON are changed, so the ledger preserves the numerical
    evidence that motivated this refresh.
    """
    conn = ledger.open_existing(root)
    try:
        recovered = ledger.recover_orphaned(conn, root)
        rows = conn.execute(
            """SELECT * FROM cases WHERE status IN ('planned','failed','interrupted')
               AND phase != ? ORDER BY temperature_K,en_Td,n2_fraction""", (EXTENSION_PHASE,)
        ).fetchall()
        now = ledger.utc_now()
        for row in rows:
            profile = solver_profile(float(row["en_Td"]))
            conn.execute("UPDATE cases SET status='planned',solver_json=?,last_error=NULL,updated_at=? WHERE case_id=?",
                         (json.dumps(profile, ensure_ascii=False, sort_keys=True), now, row["case_id"]))
            ledger.append_event(conn, root, "baseline_profile_refreshed", {
                "reason": "replace pre-validation ATOL=1 profile with validated ATOL=1e8, MXSTEP=200000 profile",
                "solver": profile,
            }, row, int(row["active_attempt"]) if row["active_attempt"] is not None else None)
        conn.commit()
        report = {"recorded_at": now, "recovered_orphaned": recovered, "case_count": len(rows),
                  "profile_example": solver_profile(80.0), "prior_attempts_preserved": True}
        ledger.atomic_json(root / "baseline_profile_refresh_20260915.json", report)
        ledger.export_overview(conn, root)
        write_live_ledger(conn, root)
        return report
    finally:
        conn.close()


def run_campaign(root: Path, phase: str) -> None:
    conn = ledger.open_existing(root)
    try:
        with ledger.campaign_lock(root):
            recovered = ledger.recover_orphaned(conn, root)
            ledger.append_event(conn, root, "campaign_run_started", {"phase": phase, "recovered_orphaned": recovered,
                                                                         "residence_time_s": TAU_RES_S})
            conn.commit()
            work = model.prepare_worktree(root / "runtime_build_surface_cstr_flow_v4")
            phases = ["high_stiffness_60Td", "high_stiffness_80Td", "normal"] if phase == "all" else [phase]
            for current_phase in phases:
                workers = 1 if current_phase.startswith("high_stiffness_") else 2
                ledger.run_phase(conn, root, work, current_phase, workers)
                ledger.export_overview(conn, root)
                write_live_ledger(conn, root)
            if phase == "all":
                # Stage 1 is complete only after every 100 s case has an
                # immutable attempt record.  Export the register before any
                # longer trajectory is launched.
                record_cases_needing_extension(root, label="100s", freeze=False)
                extension_count = schedule_200s_extensions(root)
                if extension_count:
                    ledger.run_phase(conn, root, work, EXTENSION_PHASE, 2)
                    ledger.export_overview(conn, root)
                    write_live_ledger(conn, root)
                    record_cases_needing_extension(root, label="200s", freeze=True)
            ledger.append_event(conn, root, "campaign_run_finished", {"phase": phase, **ledger.status_payload(conn)})
            conn.commit()
            ledger.export_overview(conn, root)
            write_live_ledger(conn, root)
    finally:
        conn.close()


def main() -> int:
    ledger.model = model
    ledger.SCRIPT_PATH = SCRIPT_PATH
    ledger.DEFAULT_ROOT = DEFAULT_ROOT
    ledger.solver_profile = solver_profile
    ledger.campaign_manifest = campaign_manifest
    ledger.initialise_campaign = initialise_campaign
    ledger.validate_record = validate_record
    ledger.write_live_ledger = write_live_ledger
    ledger.run_campaign = run_campaign
    return ledger.main()


if __name__ == "__main__":
    raise SystemExit(main())
