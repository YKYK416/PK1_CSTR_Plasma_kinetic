#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pure-gas-phase CW-CSTR long-time convergence and parameter scans.

The legacy ``main_pulse`` driver writes reaction rates at every sub-step and
stops after a short periodic-window criterion.  That is appropriate for ROP
analysis but not for a 10^3 s concentration-convergence study.  This utility
creates an isolated worktree below a new analysis directory, compiles a
concentration-only CW driver there, and never modifies the validated source
tree or any existing results.

The CW driver has no off phase, no waveform, and no surface chemistry.  Its
``window_dt_s`` is a numerical segmentation interval passed to ZDPlasKin; the
underlying stiff ODE solver remains adaptive.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

TOOL_DIR = Path(__file__).resolve().parent
PROJECT = TOOL_DIR.parent
BASE = PROJECT / "Reproduction" / "2017Hong" / "方案2_脉冲能效地图" / "build_pulse"
DEFAULT_OUTPUT = PROJECT / "analysis" / "cw_gasphase_longtime_20260906"
DEFAULT_TAU_RES_S = 1.0e-2
NE_CM3 = 1.17e8
T_END_S = 2000.0
REPORT_DT_S = 10.0
NH3_DRIFT_LIMIT = 1.0e-4
OTHER_DRIFT_LIMIT = 1.0e-2
RESOLUTION_LIMIT = 1.0e-3
TERMINAL_WINDOW_S = 200.0

sys.path.insert(0, str(TOOL_DIR))
import zdp_runtime as runtime  # noqa: E402


LONG_CW_DRIVER = r'''! Concentration-only pure-gas CW-CSTR long-time driver.
! CLI: n2_frac Tgas_K EN_Td window_dt_s t_end_s report_dt_s tag atol rtol tau_res_s startup_dt_s startup_duration_s transition_dt_s transition_duration_s
program main_cw_longtime
  use ZDPlasKin
  implicit none

  double precision :: time, dt_step, next_report, window_dt, t_end, report_dt
  double precision :: startup_dt, startup_duration, transition_dt, transition_duration
  double precision :: ntot, n2_frac, Tgas, EN, ne0, tau_res
  double precision :: atol_in, rtol_in, feed_density(species_max)
  double precision :: nh3, nh2, nh, n_atom, h_atom, n2, h2, elec
  integer :: iargc, ios, u, samples
  character(len=128) :: arg, tag, series_csv

  n2_frac = 1.0d0/3.0d0
  Tgas = 300.0d0
  EN = 120.0d0
  window_dt = 1.0d0
  t_end = 2000.0d0
  report_dt = 10.0d0
  tag = 'cw_longtime'
  atol_in = 1.0d0
  rtol_in = 1.0d-4
  tau_res = 1.0d-2
  startup_dt = 1.0d-3
  startup_duration = 2.0d-1
  transition_dt = window_dt
  transition_duration = startup_duration

  if (iargc() .ge. 1) then; call getarg(1,arg); read(arg,*) n2_frac; end if
  if (iargc() .ge. 2) then; call getarg(2,arg); read(arg,*) Tgas; end if
  if (iargc() .ge. 3) then; call getarg(3,arg); read(arg,*) EN; end if
  if (iargc() .ge. 4) then; call getarg(4,arg); read(arg,*) window_dt; end if
  if (iargc() .ge. 5) then; call getarg(5,arg); read(arg,*) t_end; end if
  if (iargc() .ge. 6) then; call getarg(6,arg); read(arg,*) report_dt; end if
  if (iargc() .ge. 7) then; call getarg(7,tag); end if
  if (iargc() .ge. 8) then; call getarg(8,arg); read(arg,*) atol_in; end if
  if (iargc() .ge. 9) then; call getarg(9,arg); read(arg,*) rtol_in; end if
  if (iargc() .ge. 10) then; call getarg(10,arg); read(arg,*) tau_res; end if
  if (iargc() .ge. 11) then; call getarg(11,arg); read(arg,*) startup_dt; end if
  if (iargc() .ge. 12) then; call getarg(12,arg); read(arg,*) startup_duration; end if
  if (iargc() .ge. 13) then; call getarg(13,arg); read(arg,*) transition_dt; end if
  if (iargc() .ge. 14) then; call getarg(14,arg); read(arg,*) transition_duration; end if
  if (n2_frac .le. 0.0d0 .or. n2_frac .ge. 1.0d0) stop 'n2_frac must be in (0,1)'
  if (window_dt .le. 0.0d0 .or. t_end .le. 0.0d0 .or. report_dt .le. 0.0d0 .or. &
      startup_dt .le. 0.0d0 .or. startup_duration .lt. 0.0d0 .or. &
      transition_dt .le. 0.0d0 .or. transition_duration .lt. startup_duration) stop 'time inputs must be positive'

  series_csv = 'cw_concentrations_' // trim(tag) // '.csv'
  call ZDPlasKin_init()
  call ZDPlasKin_set_config(ATOL=atol_in, RTOL=rtol_in)
  ntot = 2.446d19 * (300.0d0 / Tgas)
  ne0 = 1.17d8
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=Tgas, REDUCED_FIELD=EN)
  call ZDPlasKin_set_density('N2', n2_frac * ntot)
  call ZDPlasKin_set_density('H2', (1.0d0 - n2_frac) * ntot)
  call ZDPlasKin_set_density('E', ne0, ldens_const=.true.)
  feed_density = density
  call ZDPlasKin_set_cstr_flow(tau_res, feed_density)

  u = 20
  open(unit=u, file=series_csv, status='replace', iostat=ios)
  if (ios .ne. 0) stop 'cannot open concentration output'
  write(u,'(A)') 'time_s,NH3_cm-3,NH2_cm-3,NH_cm-3,N_cm-3,H_cm-3,N2_cm-3,H2_cm-3,E_cm-3'
  time = 0.0d0
  samples = 0
  call write_sample(u, time)
  samples = samples + 1
  next_report = report_dt

  do while (time .lt. t_end - 1.0d-10 * max(1.0d0, t_end))
    dt_step = min(window_dt, t_end - time)
    ! The initial plasma-chemical transient is exceptionally stiff.  The
    ! default values preserve the validated 1 ms segmentation over 0.2 s.
    ! Dedicated numerical-rescue calculations may additionally use a finer
    ! first interval followed by a millisecond-scale transition interval.
    if (time .lt. startup_duration) dt_step = min(dt_step, startup_dt)
    if (time .ge. startup_duration .and. time .lt. transition_duration) dt_step = min(dt_step, transition_dt)
    call ZDPlasKin_timestep(time, dt_step)
    time = time + dt_step
    if (time + 1.0d-10 * max(1.0d0, t_end) .ge. next_report .or. &
        time + 1.0d-10 * max(1.0d0, t_end) .ge. t_end) then
      call write_sample(u, time)
      samples = samples + 1
      do while (next_report .le. time + 1.0d-10 * max(1.0d0, t_end))
        next_report = next_report + report_dt
      end do
    end if
  end do
  close(u)
  write(*,'(A,A)') 'DONE tag=', trim(tag)
  write(*,'(A,ES13.5)') 't_end_s=', time
  write(*,'(A,I0)') 'samples=', samples
  call ZDPlasKin_get_density('NH3', nh3)
  write(*,'(A,ES13.5)') 'NH3_final_cm-3=', nh3

contains
  subroutine write_sample(unit, t)
    integer, intent(in) :: unit
    double precision, intent(in) :: t
    call ZDPlasKin_get_density('NH3', nh3)
    call ZDPlasKin_get_density('NH2', nh2)
    call ZDPlasKin_get_density('NH', nh)
    call ZDPlasKin_get_density('N', n_atom)
    call ZDPlasKin_get_density('H', h_atom)
    call ZDPlasKin_get_density('N2', n2)
    call ZDPlasKin_get_density('H2', h2)
    call ZDPlasKin_get_density('E', elec)
    write(unit,'(9(ES18.10,:,","))') t, nh3, nh2, nh, n_atom, h_atom, n2, h2, elec
    call flush(unit)
  end subroutine write_sample
end program main_cw_longtime
'''


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def surface_reactions_are_disabled(kinet: Path) -> None:
    active = []
    for lineno, line in enumerate(kinet.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        stripped = line.lstrip()
        if "=>" in line and "surf" in line.lower() and not stripped.startswith("#"):
            active.append(lineno)
    if active:
        raise RuntimeError(f"发现未关闭的表面反应（行 {active[:8]}），拒绝开始纯气相计算")


def ignore_base(directory: str, names: list[str]) -> set[str]:
    ignored = {"runs", "__pycache__"}
    for name in names:
        if name.endswith((".o", ".mod")) or name in {"main_pulse.exe", "main_cw_longtime.exe", "zdplaskin_m.F90", "zdplaskin_m.f90"}:
            ignored.add(name)
    return ignored


def configure_mxstep(work: Path, mxstep: int) -> None:
    if mxstep <= 5000:
        raise RuntimeError("MXSTEP 必须大于默认值 5000")
    solver = work / "zdplaskin_m.F90"
    source = solver.read_text(encoding="utf-8")
    old = "dense_j=.true.,user_supplied_jacobian=.true., &\n                                          constrained="
    new = f"dense_j=.true.,user_supplied_jacobian=.true., mxstep={mxstep}, &\n                                          constrained="
    if new in source:
        return
    if old not in source:
        raise RuntimeError("未能定位 DVODE 配置位置，拒绝修改重算工作目录")
    solver.write_text(source.replace(old, new, 1), encoding="utf-8", newline="\n")
    rc = runtime.build(work, mode="main")
    if rc:
        raise RuntimeError(f"加固 DVODE 配置后重新编译失败，退出码 {rc}")


def prepare_worktree(root: Path, mxstep: int | None = None) -> Path:
    if not BASE.is_dir():
        raise RuntimeError(f"基础算例目录不存在：{BASE}")
    surface_reactions_are_disabled(BASE / "kinet.inp")
    work = root / "worktree_cw_longtime"
    if work.exists():
        driver = work / "main_cw_longtime.F90"
        exe = work / "main_cw_longtime.exe"
        if driver.is_file() and exe.is_file():
            if mxstep is not None:
                configure_mxstep(work, mxstep)
            return work
        raise RuntimeError(f"工作目录已存在但不完整：{work}；请使用新的输出目录")

    if root.exists():
        allowed = {".claude_resources.json"}
        unexpected = [item.name for item in root.iterdir() if item.name not in allowed]
        if unexpected:
            raise RuntimeError(f"输出目录已存在非资源检查文件：{unexpected[:8]}；请使用新的输出目录")
    else:
        root.mkdir(parents=True)
    shutil.copytree(BASE, work, ignore=ignore_base)
    old_main = work / "main_pulse.F90"
    if not old_main.is_file():
        raise RuntimeError("复制后未找到 main_pulse.F90")
    old_main.unlink()
    (work / "main_cw_longtime.F90").write_text(LONG_CW_DRIVER, encoding="utf-8", newline="\n")
    print("[build] 编译独立 CW 长时浓度驱动器")
    rc = runtime.build(work, mode="full")
    if rc:
        raise RuntimeError(f"CW 长时驱动器编译失败，退出码 {rc}")
    if mxstep is not None:
        configure_mxstep(work, mxstep)
    manifest = {
        "purpose": "pure-gas-phase CW-CSTR long-time NH3 convergence study",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "base_directory": str(BASE),
        "base_kinet_sha256": sha256(BASE / "kinet.inp"),
        "surface_reactions": "verified disabled in source kinet.inp",
        "electron_density_cm-3": NE_CM3,
        "residence_time_s": "provided per trajectory; see params.json",
        "waveform": "continuous wave; no off phase, duty cycle, or physical pulse frequency",
        "driver": "main_cw_longtime.F90",
    }
    (root / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return work


def token(value: float) -> str:
    return (f"{value:g}").replace(".", "p").replace("-", "m")


def case_name(tg: float, en: float, n2: float, dt: float, tau_res_s: float) -> str:
    return f"T{token(tg)}_EN{token(en)}_N2{token(n2)}_dt{token(dt)}s_tau{token(1.0e3 * tau_res_s)}ms"


def read_series(path: Path) -> list[dict[str, float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(handle)]


def relative_range(values: list[float], floor: float = 1.0) -> float:
    if not values:
        return math.nan
    return (max(values) - min(values)) / max(abs(values[-1]), floor)


def summarize_series(rows: list[dict[str, float]], t_end: float) -> dict[str, float]:
    if not rows:
        raise RuntimeError("浓度文件为空")
    last_t = rows[-1]["time_s"]
    if abs(last_t - t_end) > max(1.0e-6, 1.0e-8 * t_end):
        raise RuntimeError(f"轨迹未到达请求终止时间：{last_t} s / {t_end} s")
    terminal_window = min(TERMINAL_WINDOW_S, 0.2 * t_end)
    tail = [row for row in rows if row["time_s"] >= t_end - terminal_window - 1.0e-9]
    result = {"samples": len(rows), "t_final_s": last_t, "terminal_window_s": terminal_window}
    for field in ("NH3_cm-3", "NH2_cm-3", "NH_cm-3", "N_cm-3", "H_cm-3", "N2_cm-3", "H2_cm-3"):
        result[field] = rows[-1][field]
        result[f"{field}_tail_rel_range"] = relative_range([row[field] for row in tail])
    result["NH3_tail_pass"] = float(result["NH3_cm-3_tail_rel_range"] <= NH3_DRIFT_LIMIT)
    result["other_tail_pass"] = float(max(result[f"{field}_tail_rel_range"] for field in ("NH2_cm-3", "NH_cm-3", "N_cm-3", "H_cm-3")) <= OTHER_DRIFT_LIMIT)
    return result


def _process_cpu_seconds(process: subprocess.Popen[str]) -> float | None:
    """Return live child CPU time without adding a psutil dependency."""
    if os.name == "nt":
        import ctypes
        from ctypes import wintypes

        creation = wintypes.FILETIME()
        exit_time = wintypes.FILETIME()
        kernel = wintypes.FILETIME()
        user = wintypes.FILETIME()
        ok = ctypes.windll.kernel32.GetProcessTimes(
            wintypes.HANDLE(int(process._handle)),  # type: ignore[attr-defined]
            ctypes.byref(creation),
            ctypes.byref(exit_time),
            ctypes.byref(kernel),
            ctypes.byref(user),
        )
        if not ok:
            return None

        def seconds(value: wintypes.FILETIME) -> float:
            ticks = (value.dwHighDateTime << 32) | value.dwLowDateTime
            return ticks / 10_000_000.0

        return seconds(kernel) + seconds(user)

    stat_path = Path(f"/proc/{process.pid}/stat")
    try:
        fields = stat_path.read_text(encoding="ascii").split()
        clock_ticks = float(os.sysconf("SC_CLK_TCK"))
        return (float(fields[13]) + float(fields[14])) / clock_ticks
    except (OSError, ValueError, IndexError):
        return None


def _activity_signature(paths: tuple[Path, ...]) -> tuple[tuple[int, int], ...]:
    signature = []
    for path in paths:
        try:
            stat = path.stat()
            signature.append((stat.st_size, stat.st_mtime_ns))
        except FileNotFoundError:
            signature.append((-1, -1))
    return tuple(signature)


def _run_with_watchdog(
    command: list[str],
    cwd: Path,
    env: dict[str, str],
    timeout_s: float,
    idle_timeout_s: float | None,
    watchdog_poll_s: float,
    activity_paths: tuple[Path, ...],
) -> tuple[int | None, bool, str | None, str, str]:
    process = subprocess.Popen(
        command,
        cwd=cwd,
        env=env,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    wall_started = time.monotonic()
    last_activity = wall_started
    previous_cpu = _process_cpu_seconds(process)
    previous_signature = _activity_signature(activity_paths)
    timeout_reason: str | None = None

    while True:
        elapsed = time.monotonic() - wall_started
        remaining = timeout_s - elapsed
        if remaining <= 0.0:
            timeout_reason = "wall timeout"
            break
        try:
            stdout, stderr = process.communicate(
                timeout=min(watchdog_poll_s, remaining)
            )
            return process.returncode, False, None, stdout or "", stderr or ""
        except subprocess.TimeoutExpired:
            current_time = time.monotonic()
            current_cpu = _process_cpu_seconds(process)
            current_signature = _activity_signature(activity_paths)
            cpu_active = (
                current_cpu is not None
                and previous_cpu is not None
                and current_cpu > previous_cpu + 1.0e-3
            )
            output_active = current_signature != previous_signature
            if cpu_active or output_active:
                last_activity = current_time
            previous_cpu = current_cpu
            previous_signature = current_signature
            if (
                idle_timeout_s is not None
                and current_time - last_activity >= idle_timeout_s
            ):
                timeout_reason = "idle timeout"
                break

    process.kill()
    stdout, stderr = process.communicate()
    return None, True, timeout_reason, stdout or "", stderr or ""


def run_case(work: Path, root: Path, tg: float, en: float, n2: float, dt: float,
             t_end: float, report_dt: float, group: str, tau_res_s: float, *, atol: float = 1.0,
             rtol: float = 1.0e-4, startup_dt_s: float = 1.0e-3,
             startup_duration_s: float = 2.0e-1, transition_dt_s: float | None = None,
             transition_duration_s: float | None = None, timeout_s: float = 1800.0,
             idle_timeout_s: float | None = None,
             watchdog_poll_s: float = 5.0) -> dict[str, object]:
    if tau_res_s <= 0.0:
        raise ValueError("tau_res_s 必须为正数")
    tag = case_name(tg, en, n2, dt, tau_res_s)
    out = root / group / tag
    if out.exists():
        series = out / f"cw_concentrations_{tag}.csv"
        params = out / "params.json"
        if series.is_file() and params.is_file():
            record = json.loads(params.read_text(encoding="utf-8"))
            record.update(summarize_series(read_series(series), t_end))
            return record
        raise RuntimeError(f"目标目录已存在但缺少完整产物：{out}")
    out.mkdir(parents=True)
    bolsig = work / "bolsigdb.dat"
    if not bolsig.is_file():
        raise RuntimeError(f"工作目录缺少 BOLSIG 数据库：{bolsig}")
    shutil.copy2(bolsig, out / bolsig.name)
    for data_file in work.glob("*.DAT"):
        shutil.copy2(data_file, out / data_file.name)
    exe = work / "main_cw_longtime.exe"
    if transition_dt_s is None:
        transition_dt_s = dt
    if transition_duration_s is None:
        transition_duration_s = startup_duration_s
    if (atol <= 0.0 or rtol <= 0.0 or startup_dt_s <= 0.0 or startup_duration_s < 0.0
            or transition_dt_s <= 0.0 or transition_duration_s < startup_duration_s
            or timeout_s <= 0.0 or watchdog_poll_s <= 0.0
            or (idle_timeout_s is not None and idle_timeout_s <= 0.0)):
        raise ValueError("容差、超时和看门狗轮询参数必须为正数")
    command = [str(exe), f"{n2:.12g}", f"{tg:.12g}", f"{en:.12g}", f"{dt:.12g}",
               f"{t_end:.12g}", f"{report_dt:.12g}", tag, f"{atol:.12g}", f"{rtol:.12g}",
               f"{tau_res_s:.12g}", f"{startup_dt_s:.12g}", f"{startup_duration_s:.12g}",
               f"{transition_dt_s:.12g}", f"{transition_duration_s:.12g}"]
    env = runtime._runtime_env(runtime.resolve_gfortran_dir())
    started = datetime.now().isoformat(timespec="seconds")
    series = out / f"cw_concentrations_{tag}.csv"
    rc, timed_out, timeout_reason, stdout, stderr = _run_with_watchdog(
        command, out, env, timeout_s, idle_timeout_s, watchdog_poll_s,
        (series, out / "bolsigdb.log"),
    )
    (out / "console.log").write_text(stdout + stderr, encoding="utf-8")
    record: dict[str, object] = {
        "tag": tag, "group": group, "tg_K": tg, "en_Td": en, "n2_fraction": n2,
        "h2_fraction": 1.0 - n2, "window_dt_s": dt, "t_end_requested_s": t_end,
        "startup_dt_s": startup_dt_s, "startup_duration_s": startup_duration_s,
        "transition_dt_s": transition_dt_s, "transition_duration_s": transition_duration_s,
        "report_interval_s": report_dt, "tau_res_s": tau_res_s, "ne_cm-3": NE_CM3,
        "forcing": "CW", "pulse_enabled": False, "returncode": rc, "timed_out": timed_out,
        "atol": atol, "rtol": rtol, "timeout_s": timeout_s,
        "idle_timeout_s": idle_timeout_s, "watchdog_poll_s": watchdog_poll_s,
        "started_at": started, "ended_at": datetime.now().isoformat(timespec="seconds"),
    }
    if rc == 0 and not timed_out:
        if series.is_file():
            record.update(summarize_series(read_series(series), t_end))
        else:
            record["failure_reason"] = "missing concentration CSV"
    else:
        record["failure_reason"] = timeout_reason if timed_out else f"process return code {rc}"
    (out / "params.json").write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
    return record


def write_csv(path: Path, records: list[dict[str, object]]) -> None:
    fields = sorted({key for record in records for key in record})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def time_horizon_summary(root: Path, records: list[dict[str, object]], horizons: list[float]) -> None:
    """Extract several integration-time cutoffs from each single long trajectory."""
    summary = []
    for record in records:
        if record.get("returncode") != 0 or record.get("timed_out"):
            continue
        tag = str(record["tag"])
        series = root / str(record["group"]) / tag / f"cw_concentrations_{tag}.csv"
        rows = read_series(series)
        final_nh3 = rows[-1]["NH3_cm-3"]
        for horizon in horizons:
            if horizon > rows[-1]["time_s"] + 1.0e-8:
                continue
            row = min(rows, key=lambda item: abs(item["time_s"] - horizon))
            summary.append({
                "tag": tag,
                "tg_K": record["tg_K"],
                "en_Td": record["en_Td"],
                "n2_fraction": record["n2_fraction"],
                "tau_res_s": record["tau_res_s"],
                "window_dt_s": record["window_dt_s"],
                "requested_horizon_s": horizon,
                "recorded_time_s": row["time_s"],
                "NH3_cm-3": row["NH3_cm-3"],
                "NH3_rel_diff_to_final": abs(row["NH3_cm-3"] - final_nh3) / max(abs(final_nh3), 1.0),
            })
    write_csv(root / "time_horizon_summary.csv", summary)


def time_study(work: Path, root: Path, windows: list[float], t_end: float, report_dt: float,
               horizons: list[float], tau_values_s: list[float], workers: int) -> None:
    conditions = [(tg, en, n2) for tg in (300.0, 400.0, 500.0)
                  for en, n2 in ((20.0, 0.1), (20.0, 0.9), (140.0, 0.1),
                                 (140.0, 0.5), (240.0, 0.9))]
    records = []
    tasks = [(tg, en, n2, dt, tau_res_s)
             for tau_res_s in tau_values_s for tg, en, n2 in conditions for dt in windows]
    if workers == 1:
        for tg, en, n2, dt, tau_res_s in tasks:
            print(f"[time-study] tau={1.0e3 * tau_res_s:g} ms, Tg={tg:g} K, E/N={en:g} Td, xN2={n2:g}, window={dt:g} s")
            records.append(run_case(work, root, tg, en, n2, dt, t_end, report_dt, "time_resolution", tau_res_s))
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(run_case, work, root, tg, en, n2, dt, t_end, report_dt,
                                "time_resolution", tau_res_s): (tg, en, n2, dt, tau_res_s)
                for tg, en, n2, dt, tau_res_s in tasks
            }
            for index, future in enumerate(as_completed(futures), 1):
                tg, en, n2, dt, tau_res_s = futures[future]
                print(f"[time-study {index}/{len(tasks)}] tau={1.0e3 * tau_res_s:g} ms, Tg={tg:g} K, E/N={en:g} Td, xN2={n2:g}, window={dt:g} s")
                records.append(future.result())
    grouped: dict[tuple[object, object, object], list[dict[str, object]]] = defaultdict(list)
    for record in records:
        grouped[(record["tg_K"], record["en_Td"], record["n2_fraction"], record["tau_res_s"])].append(record)
    for group in grouped.values():
        good = [r for r in group if r.get("returncode") == 0 and not r.get("timed_out")]
        if not good:
            continue
        ref = min(good, key=lambda r: float(r["window_dt_s"]))
        ref_nh3 = float(ref["NH3_cm-3"])
        for record in good:
            record["NH3_resolution_rel_diff"] = abs(float(record["NH3_cm-3"]) - ref_nh3) / max(abs(ref_nh3), 1.0)
            record["resolution_pass"] = float(record["NH3_resolution_rel_diff"] <= RESOLUTION_LIMIT)
    write_csv(root / "time_resolution_summary.csv", records)
    time_horizon_summary(root, records, horizons)
    (root / "time_resolution_protocol.json").write_text(json.dumps({
        "temperature_K": [300, 400, 500],
        "composition_and_field": [{"E/N_Td": en, "x_N2": n2} for en, n2 in ((20, .1), (20, .9), (140, .1), (140, .5), (240, .9))],
        "residence_time_ms": [1.0e3 * item for item in tau_values_s],
        "window_dt_s": windows, "t_end_s": t_end, "report_interval_s": report_dt,
        "NH3_tail_rel_range_limit": NH3_DRIFT_LIMIT,
        "other_species_tail_rel_range_limit": OTHER_DRIFT_LIMIT,
        "NH3_resolution_rel_diff_limit": RESOLUTION_LIMIT,
        "terminal_window_s": TERMINAL_WINDOW_S,
        "integration_time_checkpoints_s": horizons,
    }, ensure_ascii=False, indent=2), encoding="utf-8")


def full_grid(work: Path, root: Path, dt: float, t_end: float, report_dt: float,
              workers: int, tau_values_s: list[float]) -> None:
    records = []
    cases = [(tau_res_s, float(tg), float(en), n2_int / 10)
             for tau_res_s in tau_values_s
             for tg in range(300, 501, 25)
             for n2_int in range(1, 10)
             for en in range(20, 241, 20)]
    if workers == 1:
        for tau_res_s, tg, en, n2 in cases:
            group = f"grid_5832/tau_{token(1.0e3 * tau_res_s)}ms"
            print(f"[grid] tau={1.0e3 * tau_res_s:g} ms, Tg={tg:g} K, E/N={en:g} Td, xN2={n2:g}, window={dt:g} s")
            records.append(run_case(work, root, tg, en, n2, dt, t_end, report_dt, group, tau_res_s))
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            future_map = {
                executor.submit(run_case, work, root, tg, en, n2, dt, t_end, report_dt,
                                f"grid_5832/tau_{token(1.0e3 * tau_res_s)}ms", tau_res_s):
                (tau_res_s, tg, en, n2)
                for tau_res_s, tg, en, n2 in cases
            }
            for index, future in enumerate(as_completed(future_map), 1):
                tau_res_s, tg, en, n2 = future_map[future]
                print(f"[grid {index}/{len(cases)}] tau={1.0e3 * tau_res_s:g} ms, Tg={tg:g} K, E/N={en:g} Td, xN2={n2:g}")
                records.append(future.result())
    for tau_res_s in tau_values_s:
        subset = [record for record in records if float(record["tau_res_s"]) == tau_res_s]
        write_csv(root / f"grid_972_tau_{token(1.0e3 * tau_res_s)}ms_summary.csv", subset)
    write_csv(root / "grid_5832_summary.csv", records)


def long_horizon_records(root: Path) -> list[dict[str, object]]:
    """Read the five prescribed 2000 s validation trajectories, if complete."""
    specs = [
        ("long_horizon_baseline", "T400_EN140_N20p5_dt10s", 400.0, 140.0, 0.5),
        ("long_horizon_boundaries", "T300_EN20_N20p1_dt10s", 300.0, 20.0, 0.1),
        ("long_horizon_boundaries", "T300_EN240_N20p9_dt10s", 300.0, 240.0, 0.9),
        ("long_horizon_boundaries", "T500_EN20_N20p9_dt10s", 500.0, 20.0, 0.9),
        ("long_horizon_boundaries", "T500_EN240_N20p1_dt10s", 500.0, 240.0, 0.1),
    ]
    records = []
    for group, tag, tg, en, n2 in specs:
        series = root / group / tag / f"cw_concentrations_{tag}.csv"
        record: dict[str, object] = {
            "group": group, "tag": tag, "tg_K": tg, "en_Td": en, "n2_fraction": n2,
            "h2_fraction": 1.0 - n2, "window_dt_s": 10.0, "t_end_requested_s": 2000.0,
        }
        if not series.is_file():
            record.update({"complete": False, "failure_reason": "missing concentration CSV"})
        else:
            try:
                record.update(summarize_series(read_series(series), 2000.0))
                record["complete"] = True
                record["acceptance_pass"] = bool(record["NH3_tail_pass"] and record["other_tail_pass"])
            except RuntimeError as exc:
                record.update({"complete": False, "failure_reason": str(exc)})
        records.append(record)
    return records


def gate_and_grid(work: Path, root: Path, dt: float, grid_t_end: float, report_dt: float,
                  workers: int, poll_s: float, timeout_s: float, tau_values_s: list[float]) -> None:
    """Wait for the long-horizon gate, then launch the grid only on a pass."""
    deadline = time.monotonic() + timeout_s
    while True:
        records = long_horizon_records(root)
        write_csv(root / "long_horizon_acceptance.csv", records)
        if all(bool(record.get("complete")) for record in records):
            break
        if time.monotonic() >= deadline:
            raise RuntimeError("等待长时轨迹完成超时；未启动主网格")
        print("[gate] 长时轨迹尚未全部完成，继续等待")
        time.sleep(poll_s)
    if not all(bool(record.get("acceptance_pass")) for record in records):
        failed = [str(record["tag"]) for record in records if not record.get("acceptance_pass")]
        raise RuntimeError("长时 NH3 验收失败；未启动主网格：" + ", ".join(failed))
    print("[gate] 5 条长时轨迹均通过；启动 972 点主网格")
    full_grid(work, root, dt, grid_t_end, report_dt, workers, tau_values_s)


def parse_windows(value: str) -> list[float]:
    windows = [float(item.strip()) for item in value.split(",") if item.strip()]
    if not windows or any(item <= 0 for item in windows):
        raise argparse.ArgumentTypeError("窗口长度应为逗号分隔的正数")
    return sorted(set(windows))


def parse_positive_list(value: str) -> list[float]:
    values = [float(item.strip()) for item in value.split(",") if item.strip()]
    if not values or any(item <= 0 for item in values):
        raise argparse.ArgumentTypeError("数值应为逗号分隔的正数")
    return sorted(set(values))


def main() -> int:
    parser = argparse.ArgumentParser(description="Pure-gas CW-CSTR long-time NH3 convergence study")
    parser.add_argument("mode", choices=("prepare", "time-study", "grid", "gate-grid"))
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--t-end-s", type=float, default=T_END_S)
    parser.add_argument("--report-dt-s", type=float, default=REPORT_DT_S)
    parser.add_argument("--windows", type=parse_windows, default=[0.05, 0.2, 1.0, 5.0])
    parser.add_argument("--horizons-s", type=parse_positive_list, default=[10.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0])
    parser.add_argument("--grid-window-dt-s", type=float)
    parser.add_argument("--tau-ms", type=parse_positive_list, default=[10.0],
                        help="停留时间列表，单位 ms；默认 10")
    parser.add_argument("--mxstep", type=int,
                        help="DVODE 单次调用允许的最大内部步数；仅用于数值加固重试")
    parser.add_argument("--workers", type=int, default=1, help="grid mode: concurrently running independent cases")
    parser.add_argument("--poll-s", type=float, default=60.0, help="gate-grid mode: long-horizon polling interval")
    parser.add_argument("--gate-timeout-s", type=float, default=7200.0, help="gate-grid mode: maximum wait before failing")
    args = parser.parse_args()
    if args.t_end_s <= 0 or args.report_dt_s <= 0:
        parser.error("t-end-s 和 report-dt-s 必须为正数")
    if args.workers <= 0:
        parser.error("workers 必须为正整数")
    tau_values_s = [item * 1.0e-3 for item in args.tau_ms]
    root = args.output_root.resolve()
    work = prepare_worktree(root, args.mxstep)
    if args.mode == "time-study":
        time_study(work, root, args.windows, args.t_end_s, args.report_dt_s, args.horizons_s,
                   tau_values_s, args.workers)
    elif args.mode in ("grid", "gate-grid"):
        if args.grid_window_dt_s is None or args.grid_window_dt_s <= 0:
            parser.error("grid/gate-grid 模式必须设置 --grid-window-dt-s")
        if args.mode == "grid":
            full_grid(work, root, args.grid_window_dt_s, args.t_end_s, args.report_dt_s,
                      args.workers, tau_values_s)
        else:
            if len(tau_values_s) != 1:
                parser.error("gate-grid 只支持单一停留时间；多停留时间请先运行 time-study 后使用 grid")
            if args.poll_s <= 0 or args.gate_timeout_s <= 0:
                parser.error("poll-s 和 gate-timeout-s 必须为正数")
            gate_and_grid(work, root, args.grid_window_dt_s, args.t_end_s, args.report_dt_s,
                          args.workers, args.poll_s, args.gate_timeout_s, tau_values_s)
    print(json.dumps({"output_root": str(root), "mode": args.mode}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
