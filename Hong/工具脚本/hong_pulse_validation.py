"""Isolated numerical acceptance gate for the pulsed corrected-Hong CSTR model.

This program never edits the validated continuous-CSTR build.  It copies the
runtime inputs into a new empty directory, patches only that copied pulse
driver for prompt checkpointing, compiles it, and executes a small
reproducibility matrix.  A result is accepted only when the executable exits,
the periodic CSTR criterion is declared, and retained CSV values are finite
and non-negative within a documented numerical tolerance.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


TOOL_DIR = Path(__file__).resolve().parent
PROJECT = TOOL_DIR.parent
BASE = PROJECT / "Reproduction" / "2017Hong" / "方案2_脉冲能效地图" / "build_pulse"
DEFAULT_OUTPUT = PROJECT / "analysis" / "p0_pulse_validation_20260831"
KEEP_NAMES = {"preprocessor.exe", "bolsigdb.dat", "kinet.inp", "main_pulse.F90", "dvode_f90_m.F90"}
KEEP_SUFFIXES = {".dll", ".lib", ".dat", ".DAT"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_runtime_inputs(source: Path, target: Path) -> None:
    target.mkdir(parents=True, exist_ok=False)
    for item in source.iterdir():
        if item.is_file() and (item.name in KEEP_NAMES or item.suffix in KEEP_SUFFIXES):
            shutil.copy2(item, target / item.name)


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"Expected exactly one {label} anchor, found {count}")
    return text.replace(old, new, 1)


def patch_validation_driver(path: Path) -> None:
    """Patch a copied driver only; governing reaction equations are untouched."""
    text = path.read_text(encoding="utf-8")
    text = replace_once(
        text,
        "integer :: converged, stable_streak, min_cycles, stable_required, flow_enabled\n",
        "integer :: converged, stable_streak, min_cycles, stable_required, flow_enabled, write_rates_enabled\n",
        "integer declaration",
    )
    text = replace_once(
        text,
        "  stable_required = 3\n",
        "  stable_required = 3\n  write_rates_enabled = 1\n",
        "default configuration",
    )
    parser_anchor = "  if (iargc() .ge. 12) then; call getarg(12, arg); read(arg, *) tau_res; end if\n"
    parser_patch = parser_anchor + (
        "  if (iargc() .ge. 13) then\n"
        "    call getarg(13, arg)\n"
        "    if (trim(adjustl(arg)) .eq. 'metrics') write_rates_enabled = 0\n"
        "  end if\n"
    )
    text = replace_once(
        text,
        "  n_sub_off = 8\n\n  ! ---------- init ----------",
        "  n_sub_off = 8\n"
        "  if (iargc() .ge. 14) then; call getarg(14, arg); read(arg, *) n_sub_on; end if\n"
        "  if (iargc() .ge. 15) then; call getarg(15, arg); read(arg, *) n_sub_off; end if\n"
        "\n  ! ---------- init ----------",
        "substep override",
    )
    text = replace_once(text, parser_anchor, parser_patch, "argument parser")
    text = replace_once(
        text,
        "    call write_row(u, time, 1)\n\n    ! ---- afterglow phase ----",
        "    call write_row(u, time, 1)\n    call flush(u)\n\n    ! ---- afterglow phase ----",
        "on-phase checkpoint",
    )
    text = replace_once(
        text,
        "    call write_row(u, time, 2)\n\n    dnh3_on",
        "    call write_row(u, time, 2)\n    call flush(u)\n\n    dnh3_on",
        "off-phase checkpoint",
    )
    text = replace_once(
        text,
        "stable_streak, tau_res\n    cycles_used = k",
        "stable_streak, tau_res\n    call flush(uc)\n    cycles_used = k",
        "cycle checkpoint",
    )
    old_rate_block = """      if (jj .eq. 1) then
        if (iphase .eq. 1) then
          call write_rates(ur, tcur, dt1, k, 'on')
        else
          call write_rates(ur, tcur, dt1, k, 'off')
        end if
      else
        if (iphase .eq. 1) then
          call write_rates(ur, tcur, dti, k, 'on')
        else
          call write_rates(ur, tcur, dti, k, 'off')
        end if
      end if
"""
    new_rate_block = """      if (write_rates_enabled .eq. 1) then
        if (jj .eq. 1) then
          if (iphase .eq. 1) then
            call write_rates(ur, tcur, dt1, k, 'on')
          else
            call write_rates(ur, tcur, dt1, k, 'off')
          end if
        else
          if (iphase .eq. 1) then
            call write_rates(ur, tcur, dti, k, 'on')
          else
            call write_rates(ur, tcur, dti, k, 'off')
          end if
        end if
      end if
"""
    text = replace_once(text, old_rate_block, new_rate_block, "rate-output block")
    path.write_text(text, encoding="utf-8", newline="\n")


def build_case(build_dir: Path) -> None:
    sys.path.insert(0, str(TOOL_DIR))
    import zdp_gui as gui

    def safe_log(message: object) -> None:
        print(str(message).encode("ascii", errors="backslashreplace").decode("ascii"), flush=True)

    params = {
        "en": "120", "tg": "300", "n2frac": "0.3333", "tend": "1", "en_off": "0.1",
        "freq": "100", "duty": "0.2", "cycles": "90", "ne_mode": "fix", "tau_res": "1e-2",
        "atol": "", "rtol": "", "en_list": "", "recompile": False, "pulse_enabled": True,
        "tag": "p0_compile",
    }
    rc = gui.execute_pipeline(build_dir, params, do_build=True, do_run=False, log_cb=safe_log)
    if rc:
        raise RuntimeError(f"Validation driver build failed with rc={rc}")


def finite_nonnegative_csv(path: Path, negative_tolerance: float = 1.0e-9) -> dict:
    if not path.is_file() or path.stat().st_size == 0:
        return {"ok": False, "reason": "missing_or_empty", "rows": 0, "minimum": None}
    rows, minimum = 0, math.inf
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows += 1
            for key, value in row.items():
                if key in {"time_s", "phase", "Te_eV", "EN_Td", "Ptot"} or value in (None, ""):
                    continue
                try:
                    # ES13.5 is one character too narrow for a three-digit
                    # exponent, so legacy ZDPlasKin output may contain
                    # ``3.10-188`` rather than ``3.10E-188``.
                    normalized = re.sub(r"^([+-]?(?:\d+\.?\d*|\.\d+))([+-]\d{2,3})$", r"\1E\2", value.strip())
                    numeric = float(normalized)
                except ValueError:
                    return {"ok": False, "reason": f"non_numeric:{key}", "rows": rows, "minimum": None}
                if not math.isfinite(numeric):
                    return {"ok": False, "reason": f"non_finite:{key}", "rows": rows, "minimum": None}
                minimum = min(minimum, numeric)
    return {"ok": bool(rows and minimum >= -negative_tolerance), "reason": "ok" if rows and minimum >= -negative_tolerance else "negative_or_empty", "rows": rows, "minimum": minimum if rows else None}


def summarize_cycles(path: Path) -> dict:
    if not path.is_file() or path.stat().st_size == 0:
        return {"rows": 0, "periodic": False, "last_streak": None, "last_residuals": {}}
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return {"rows": 0, "periodic": False, "last_streak": None, "last_residuals": {}}
    last = rows[-1]
    residuals = {name: float(last[name]) for name in ("rel_state_max", "rel_dNH3_cycle", "rel_energy_cycle")}
    streak = int(last["stable_streak"])
    return {"rows": len(rows), "periodic": streak >= 3 and all(value <= 1.0e-3 for value in residuals.values()), "last_streak": streak, "last_residuals": residuals}


def run_variant(build_dir: Path, runs_root: Path, variant: dict, timeout_s: int) -> dict:
    sys.path.insert(0, str(TOOL_DIR))
    import zdp_runtime as runtime

    tag = variant["tag"]
    outdir = runs_root / tag
    outdir.mkdir(parents=True, exist_ok=False)
    for item in build_dir.glob("*.DAT"):
        shutil.copy2(item, outdir / item.name)
    shutil.copy2(build_dir / "bolsigdb.dat", outdir / "bolsigdb.dat")
    executable = build_dir / "main_pulse.exe"
    command = [
        str(executable), "0.3333", "300", "120", "0.1", "100", "0.2", "90", "fix", tag,
        variant["atol"], variant["rtol"], "1e-2", "metrics", str(variant["n_sub_on"]), str(variant["n_sub_off"]),
    ]
    started = datetime.now().isoformat(timespec="seconds")
    status = {"variant": variant, "command": command, "started_at": started, "timeout_s": timeout_s}
    try:
        proc = subprocess.run(command, cwd=outdir, env=runtime._runtime_env(runtime.resolve_gfortran_dir()),
                              text=True, capture_output=True, timeout=timeout_s, check=False)
        status.update({"returncode": proc.returncode, "timed_out": False, "stdout": proc.stdout, "stderr": proc.stderr})
    except subprocess.TimeoutExpired as exc:
        status.update({"returncode": None, "timed_out": True, "stdout": exc.stdout or "", "stderr": exc.stderr or ""})
    (outdir / "console.log").write_text((status["stdout"] or "") + (status["stderr"] or ""), encoding="utf-8", errors="replace")
    series = finite_nonnegative_csv(outdir / f"pulse_series_{tag}.csv")
    cycles = summarize_cycles(outdir / f"pulse_cycles_{tag}.csv")
    status["ended_at"] = datetime.now().isoformat(timespec="seconds")
    status["series_check"] = series
    status["cycle_check"] = cycles
    status["accepted"] = bool(status["returncode"] == 0 and not status["timed_out"] and series["ok"] and cycles["periodic"])
    (outdir / "run_status.json").write_text(json.dumps(status, ensure_ascii=False, indent=2), encoding="utf-8")
    return status


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="P0 pulse-baseline acceptance gate")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout-s", type=int, default=60)
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    root = args.output_root.resolve()
    if root.exists():
        raise SystemExit(f"Output root already exists; refusing to overwrite: {root}")
    root.mkdir(parents=True)
    build_dir = root / "build"
    copy_runtime_inputs(BASE, build_dir)
    patch_validation_driver(build_dir / "main_pulse.F90")
    manifest = {
        "purpose": "P0 pulsed-CSTR numerical acceptance gate; no kinetic-input modification",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "source_build": str(BASE),
        "source_kinet_sha256": sha256(BASE / "kinet.inp"),
        "source_driver_sha256": sha256(BASE / "main_pulse.F90"),
        "validation_driver_sha256": sha256(build_dir / "main_pulse.F90"),
    }
    (root / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    build_case(build_dir)
    variants = [
        {"tag": "nominal", "atol": "1", "rtol": "1e-4", "n_sub_on": 12, "n_sub_off": 8},
        {"tag": "refined", "atol": "0.1", "rtol": "1e-5", "n_sub_on": 48, "n_sub_off": 32},
    ]
    results = [run_variant(build_dir, root / "runs", variant, args.timeout_s) for variant in variants]
    payload = {"manifest": manifest, "results": results, "accepted": all(item["accepted"] for item in results)}
    (root / "validation_summary.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"output_root": str(root), "accepted": payload["accepted"]}, ensure_ascii=False))
    return 0 if payload["accepted"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
