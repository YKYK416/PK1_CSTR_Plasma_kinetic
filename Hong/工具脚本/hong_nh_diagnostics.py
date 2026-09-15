#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Isolated Hong NH-source diagnostic branches and ROP integration.

Each branch receives a copied, rebuilt kinetic case.  The source ``kinet.inp``
is never modified.  Results are valid as diagnostic perturbations only: a
disabled pathway or parameterized Rydberg loss is not a replacement mechanism.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path

import zdp_gui as gui


TOOL_DIR = Path(__file__).resolve().parent
PROJECT_DIR = TOOL_DIR.parent
DEFAULT_BASE = PROJECT_DIR / "Reproduction" / "2017Hong" / "方案2_脉冲能效地图" / "build_pulse"

BRANCHES = {
    "baseline": "Unmodified kinetic input.",
    "no_rydberg_nh": "Disable N + H2(RYDBERG_SUM) -> H + NH.",
    "no_named_h2star_nh": "Disable the four named electronic H2* NH channels.",
    "no_excited_n_nh": "Disable N(2D)+H2 and N(2P)+H2 NH channels.",
    "rydberg_loss_x1": "Add named-H2* equivalent direct Rydberg wall relaxation.",
    "rydberg_loss_x10": "Add 10x named-H2* equivalent direct Rydberg wall relaxation.",
    "rydberg_loss_x100": "Add 100x named-H2* equivalent direct Rydberg wall relaxation.",
    "rydberg_loss_x1000": "Add 1000x named-H2* equivalent direct Rydberg wall relaxation.",
}

PRIMARY_PATTERNS = [
    ("Vib-H2", r"^N\+H2\(V[123]\)=>H\+NH$"),
    ("Named-H2*", r"^N\+H2\((B3SIG|B1SIG|C3PI|A3SIG)\)=>H\+NH$"),
    ("Rydberg-H2*", r"^N\+H2\(RYDBERG_SUM\)=>H\+NH$"),
    ("N(2D)+H2", r"^N\(2D\)\+H2=>H\+NH$"),
    ("N(2P)+H2", r"^N\(2P\)\+H2=>H\+NH$"),
    ("Association H+N+M", r"^H\+N\+(N2|H2)=>NH\+\1$"),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(65536), b""):
            h.update(block)
    return h.hexdigest()


def _disable_line(text: str, needle: str) -> tuple[str, int]:
    out, count = [], 0
    for line in text.splitlines(keepends=True):
        if needle in line and not line.lstrip().startswith("#DIAGNOSTIC_OFF#"):
            out.append("#DIAGNOSTIC_OFF# " + line)
            count += 1
        else:
            out.append(line)
    return "".join(out), count


def branch_kinet_text(text: str, branch: str) -> tuple[str, list[str]]:
    """Return a minimally perturbed input plus a human-readable change list."""
    if branch not in BRANCHES:
        raise ValueError(f"Unknown branch: {branch}")
    changes: list[str] = []
    if branch == "baseline":
        return text, changes
    if branch == "no_rydberg_nh":
        text, n = _disable_line(text, "N + H2(RYDBERG_SUM) => H + NH")
        if n != 1:
            raise RuntimeError(f"Expected one Rydberg NH reaction, found {n}")
        return text, ["Disabled N + H2(RYDBERG_SUM) => H + NH"]
    if branch == "no_named_h2star_nh":
        for state in ("B3SIG", "B1SIG", "C3PI", "A3SIG"):
            text, n = _disable_line(text, f"N + H2({state}) => H + NH")
            if n != 1:
                raise RuntimeError(f"Expected one {state} NH reaction, found {n}")
            changes.append(f"Disabled N + H2({state}) => H + NH")
        return text, changes
    if branch == "no_excited_n_nh":
        for state in ("N(2D)", "N(2P)"):
            text, n = _disable_line(text, f"{state} + H2 => H + NH")
            if n != 1:
                raise RuntimeError(f"Expected one {state}+H2 NH reaction, found {n}")
            changes.append(f"Disabled {state} + H2 => H + NH")
        return text, changes

    factor = int(branch.rsplit("x", 1)[1])
    anchor = "H2(A3SIG) => H2"
    lines = text.splitlines(keepends=True)
    index = next((i for i, line in enumerate(lines) if line.lstrip().startswith(anchor)), None)
    if index is None:
        raise RuntimeError("Cannot locate named-H2* wall-relaxation anchor")
    # The ZDPlasKin reactions section accepts #OFF# but not arbitrary # comment
    # records, so keep the diagnostic provenance exclusively in the manifest.
    addition = f"H2(RYDBERG_SUM) => H2 ! {factor}.0d0/(GAMMA_D+H2_WALL_SECOND_PART_E)\n"
    lines.insert(index + 1, addition)
    return "".join(lines), [f"Added diagnostic Rydberg direct wall relaxation at {factor}x named-H2* rate"]


def copy_case_inputs(base: Path, target: Path) -> None:
    """Copy only required source/runtime inputs, never old outputs or run folders."""
    target.mkdir(parents=True, exist_ok=False)
    keep_names = {"preprocessor.exe", "bolsigdb.dat", "main_pulse.F90", "dvode_f90_m.F90"}
    keep_suffixes = {".dll", ".lib", ".dat", ".DAT"}
    for src in base.iterdir():
        if not src.is_file():
            continue
        if src.name in keep_names or src.suffix in keep_suffixes:
            shutil.copy2(src, target / src.name)


def _finite(value) -> bool:
    try:
        return float(value) == float(value) and abs(float(value)) != float("inf")
    except (TypeError, ValueError):
        return False


def integrate_nh_channels(rates_csv: Path, first_cycle: int | None = None) -> tuple[list[dict], float, float]:
    """Right-rule integrate primary NH sources, optionally from a late-cycle window."""
    import re

    totals = {name: 0.0 for name, _ in PRIMARY_PATTERNS}
    starts, ends = [], []
    with rates_csv.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        required = {"time_s", "dt_s", "reaction", "rate_cm-3s-1"}
        if not required <= set(reader.fieldnames or []):
            raise RuntimeError("ROP CSV 缺少 time_s/dt_s/reaction/rate_cm-3s-1，无法积分")
        for row in reader:
            try:
                t, dt, rate = float(row["time_s"]), float(row["dt_s"]), float(row["rate_cm-3s-1"])
            except (KeyError, TypeError, ValueError):
                continue
            if first_cycle is not None:
                try:
                    if int(row["cycle"]) < first_cycle:
                        continue
                except (KeyError, TypeError, ValueError):
                    continue
            if not (_finite(t) and _finite(dt) and _finite(rate)) or dt <= 0:
                continue
            starts.append(t - dt)
            ends.append(t)
            for name, pattern in PRIMARY_PATTERNS:
                if re.search(pattern, row["reaction"] or ""):
                    totals[name] += rate * dt
                    break
    total = sum(totals.values())
    rows = [{"channel": name, "phi_cm-3": value,
             "share_pct": 100.0 * value / total if total else float("nan")}
            for name, value in totals.items()]
    if not starts:
        raise RuntimeError("ROP CSV 没有有效的带 dt 样本")
    return rows, min(starts), max(ends)


def steady_window_start(cycles_csv: Path, keep_cycles: int = 3) -> tuple[int, int]:
    """Return the final contiguous cycle window after the solver declares P0 steady."""
    with cycles_csv.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise RuntimeError("pulse_cycles CSV 为空，无法定位稳态 ROP 窗口")
    last = rows[-1]
    try:
        last_cycle = int(last["cycle"])
        stable = int(float(last["stable_streak"]))
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeError("pulse_cycles CSV 缺少周期稳态列") from exc
    if stable < keep_cycles:
        raise RuntimeError(f"末周期 stable_streak={stable}，不足 {keep_cycles} 个周期")
    return max(1, last_cycle - keep_cycles + 1), last_cycle


def write_summary(outdir: Path, branch: str, params: dict, manifest: dict) -> Path:
    rates = outdir / f"pulse_rates_{params['tag']}.csv"
    record = json.loads((outdir / "params.json").read_text(encoding="utf-8"))
    if record.get("status") != "success":
        raise RuntimeError(f"算例未通过周期稳态验证：{record.get('status')}；拒绝积分启动阶段 ROP")
    cycle_first, cycle_last = steady_window_start(outdir / f"pulse_cycles_{params['tag']}.csv")
    rows, start, end = integrate_nh_channels(rates, first_cycle=cycle_first)
    metrics = gui.extract_metrics("pulse", outdir, params["tag"])
    path = outdir / "nh_channel_summary.csv"
    fields = ["branch", "channel", "phi_cm-3", "share_pct", "window_start_s", "window_end_s",
              "run_status", "failure_reason", "mechanism_sha256", "en_Td", "Tgas_K", "n2_fraction",
              "frequency_Hz", "duty", "cycles_requested", "cycles_used", "cycle_first", "cycle_last",
              "tau_res_s", "ne_mode", "NH3_final_cm-3", "rel_state_max", "rel_dNH3_cycle"]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "branch": branch, **row, "window_start_s": start, "window_end_s": end,
                "run_status": record.get("status"), "failure_reason": record.get("failure_reason") or "",
                "mechanism_sha256": manifest["branch_kinet_sha256"], "en_Td": params["en"],
                "Tgas_K": params["tg"], "n2_fraction": params["n2frac"], "frequency_Hz": params["freq"],
                "duty": "1.0" if not params["pulse_enabled"] else params["duty"],
                "cycles_requested": params["cycles"], "cycles_used": metrics.get("cycles_used", ""),
                "cycle_first": cycle_first, "cycle_last": cycle_last, "tau_res_s": params["tau_res"],
                "ne_mode": params["ne_mode"], "NH3_final_cm-3": metrics.get("NH3_final_cm-3", ""),
                "rel_state_max": metrics.get("rel_state_max", ""),
                "rel_dNH3_cycle": metrics.get("rel_dNH3_cycle", ""),
            })
    return path


def run_branch(base: Path, root: Path, branch: str, params: dict) -> dict:
    branch_dir = root / branch
    copy_case_inputs(base, branch_dir)
    source_kinet = base / "kinet.inp"
    text, changes = branch_kinet_text(source_kinet.read_text(encoding="utf-8"), branch)
    kinet = branch_dir / "kinet.inp"
    kinet.write_text(text, encoding="utf-8")
    manifest = {
        "branch": branch,
        "description": BRANCHES[branch],
        "diagnostic_only": True,
        "source_kinet": str(source_kinet),
        "source_kinet_sha256": sha256(source_kinet),
        "branch_kinet_sha256": sha256(kinet),
        "changes": changes,
        "created_at": datetime.now().isoformat(timespec="seconds"),
    }
    (branch_dir / "diagnostic_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    tag = f"p3_{branch}"
    run_params = {**params, "tag": tag}
    rc = gui.execute_pipeline(branch_dir, run_params, do_build=True, do_run=True, log_cb=print)
    flavor = gui.detect_flavor(branch_dir) or "pulse"
    out_rel = gui.out_rel_for(run_params, flavor)
    outdir = branch_dir / out_rel
    record = json.loads((outdir / "params.json").read_text(encoding="utf-8")) if (outdir / "params.json").is_file() else {}
    summary = write_summary(outdir, branch, run_params, manifest) if (rc == 0 and record.get("status") == "success") else None
    return {"branch": branch, "rc": rc, "run_status": record.get("status", "missing_record"),
            "outdir": str(outdir), "summary": str(summary) if summary else "",
            "mechanism_sha256": manifest["branch_kinet_sha256"]}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Hong NH mechanism diagnostic branches")
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE, help="validated build_pulse directory")
    parser.add_argument("--output-root", type=Path, required=True, help="new empty diagnostic root")
    parser.add_argument("--branches", nargs="+", choices=tuple(BRANCHES), default=list(BRANCHES))
    parser.add_argument("--en", default="120")
    parser.add_argument("--tg", default="300")
    parser.add_argument("--n2-frac", dest="n2frac", default="0.3333")
    parser.add_argument("--freq", default="1000")
    parser.add_argument("--cycles", default="3")
    parser.add_argument("--ne-mode", default="fix", choices=("fix", "decay"))
    parser.add_argument("--tau-res", dest="tau_res", default="1e-2",
                        help="CSTR 停留时间 (s); 0 表示封闭批式，不能用于稳态机制结论")
    parser.add_argument("--pulse", action="store_true", help="use a pulse instead of continuous field")
    parser.add_argument("--duty", default="0.2")
    parser.add_argument("--en-off", dest="en_off", default="0.1")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    base, root = args.base.resolve(), args.output_root.resolve()
    if not (base / "kinet.inp").is_file():
        raise SystemExit(f"Base mechanism not found: {base / 'kinet.inp'}")
    if root.exists():
        raise SystemExit(f"Output root already exists; refusing to overwrite: {root}")
    root.mkdir(parents=True)
    params = {
        "en": args.en, "tg": args.tg, "n2frac": args.n2frac, "tend": "1",
        "en_off": args.en_off, "freq": args.freq, "duty": args.duty, "cycles": args.cycles,
        "ne_mode": args.ne_mode, "tau_res": args.tau_res, "atol": "", "rtol": "", "en_list": "", "recompile": False,
        "pulse_enabled": bool(args.pulse),
    }
    results = []
    for branch in args.branches:
        print(f"\n=== Diagnostic branch: {branch} ===")
        results.append(run_branch(base, root, branch, params))
    with (root / "diagnostic_runs.json").open("w", encoding="utf-8") as fh:
        json.dump(results, fh, ensure_ascii=False, indent=2)
    failed = [r for r in results if r["rc"] != 0 or r.get("run_status") != "success"]
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
