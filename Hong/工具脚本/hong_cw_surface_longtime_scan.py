#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CW N2/H2 DBD gas--surface convergence and parameter scans.

This driver uses only the literature-traceable Hong 2017/2018 corrected
mechanism, including the Table-5 metal-column surface reaction set.  It does
not use the Shao--Mesbah Fe-DFT microkinetic model.  Each calculation is
executed in an isolated copied worktree; no validated input or prior result is
changed.

``window_dt_s`` is the maximum outer call interval to ZDPlasKin's adaptive
stiff integrator, not a physical pulse period.  ``internal_hmax_s`` is passed
to DVODE as its true maximum internal integration step.  The imposed field is
continuous and no pulse waveform is present.
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
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Callable

TOOL_DIR = Path(__file__).resolve().parent
PROJECT = TOOL_DIR.parent
BASE = PROJECT / "Reproduction" / "2017Hong" / "literature_extract_2017_2018_corrected"
DVODE_SOURCE = PROJECT / "Reproduction" / "2017Hong" / "build_full" / "dvode_f90_m.F90"
BOLSIG_LIBRARY_SOURCE = PROJECT / "Reproduction" / "2017Hong" / "build_full" / "bolsig_x86_64_g.lib"
# This full Hong-run database is the paired source for the 41 electron-impact
# processes in the base gas-phase mechanism.  The compact ``build`` database
# omits several of the vibrational/electronic channels (for example N2(v1)).
BOLSIG_DATABASE_SOURCE = PROJECT / "Reproduction" / "2017Hong" / "build_full" / "bolsigdb.dat"
KINETIC_INPUT = "kinet_hong_2017_2018_corrected.inp"
DEFAULT_OUTPUT = PROJECT / "analysis" / "cw_hong2017_2018_surface_convergence_20260908"

# The Hong input is treated as a closed 0D baseline: this driver deliberately
# has no feed, outlet, residence-time, or CSTR branch.  The imposed field is CW.
NE_CM3 = 1.17e8
# Hong 2017/2018: ST=1e15 cm^-2 and V/A=0.007 cm.
SURFACE_SITE_DENSITY_CM3 = 1.0e15 / 7.0e-3
T_END_S = 4000.0
REPORT_DT_S = 10.0
NH3_DRIFT_LIMIT = 1.0e-4
GAS_DRIFT_LIMIT = 1.0e-2
COVERAGE_DRIFT_LIMIT = 1.0e-4
SITE_BALANCE_LIMIT = 1.0e-8
RESOLUTION_LIMIT = 1.0e-3
TERMINAL_WINDOW_S = 400.0
ACTIVE_STATE_DRIFT_LIMIT = 1.0e-3
PRODUCT_RATE_STABILITY_LIMIT = 1.0e-2
SURFACE_COVERAGE_FLOOR = 1.0e-8

sys.path.insert(0, str(TOOL_DIR))
import zdp_runtime as runtime  # noqa: E402


CW_SURFACE_DRIVER = r'''! Closed-0D CW gas--surface driver.
! CLI: n2_frac Tgas_K EN_Td window_dt_s t_end_s report_dt_s tag atol rtol internal_hmax_s early_diagnostics_flag mxstep site_projection_flag
program main_cw_surface
  use ZDPlasKin
  implicit none
  double precision :: time, dt_step, next_report, window_dt, target_window, t_end, report_dt
  double precision :: ntot, n2_frac, Tgas, EN, ne0
  double precision :: atol_in, rtol_in, hmax_in
  double precision :: n2, h2, nh3, nh2, nh, n_atom, h_atom, elec
  double precision :: surf, hsurf, nsurf, nhsurf, nh2surf, sites, site_projection_max_rel
  integer :: iargc, ios, u, samples, early_diagnostics_flag, mxstep_in, site_projection_flag
  character(len=128) :: arg, tag, series_csv

  n2_frac = 0.5d0
  Tgas = 400.0d0
  EN = 140.0d0
  window_dt = 1.0d0
  t_end = 4000.0d0
  report_dt = 10.0d0
  tag = 'cw_surface'
  atol_in = 1.0d0
  rtol_in = 1.0d-4
  hmax_in = -1.0d0
  early_diagnostics_flag = 1
  mxstep_in = 5000
  site_projection_flag = 0
  site_projection_max_rel = 0.0d0

  if (iargc() .ge. 1) then; call getarg(1,arg); read(arg,*) n2_frac; end if
  if (iargc() .ge. 2) then; call getarg(2,arg); read(arg,*) Tgas; end if
  if (iargc() .ge. 3) then; call getarg(3,arg); read(arg,*) EN; end if
  if (iargc() .ge. 4) then; call getarg(4,arg); read(arg,*) window_dt; end if
  if (iargc() .ge. 5) then; call getarg(5,arg); read(arg,*) t_end; end if
  if (iargc() .ge. 6) then; call getarg(6,arg); read(arg,*) report_dt; end if
  if (iargc() .ge. 7) then; call getarg(7,tag); end if
  if (iargc() .ge. 8) then; call getarg(8,arg); read(arg,*) atol_in; end if
  if (iargc() .ge. 9) then; call getarg(9,arg); read(arg,*) rtol_in; end if
  if (iargc() .ge. 10) then; call getarg(10,arg); read(arg,*) hmax_in; end if
  if (iargc() .ge. 11) then; call getarg(11,arg); read(arg,*) early_diagnostics_flag; end if
  if (iargc() .ge. 12) then; call getarg(12,arg); read(arg,*) mxstep_in; end if
  if (iargc() .ge. 13) then; call getarg(13,arg); read(arg,*) site_projection_flag; end if
  if (n2_frac .le. 0.0d0 .or. n2_frac .ge. 1.0d0) stop 'n2_frac must be in (0,1)'
  if (window_dt .le. 0.0d0 .or. t_end .le. 0.0d0 .or. report_dt .le. 0.0d0) stop 'time inputs must be positive'
  if (early_diagnostics_flag .ne. 0 .and. early_diagnostics_flag .ne. 1) stop 'early_diagnostics_flag must be 0 or 1'
  if (mxstep_in .le. 0) stop 'mxstep must be positive'
  if (site_projection_flag .ne. 0 .and. site_projection_flag .ne. 1) stop 'site_projection_flag must be 0 or 1'

  series_csv = 'cw_surface_' // trim(tag) // '.csv'
  call ZDPlasKin_init()
  if (hmax_in .gt. 0.0d0) then
    call ZDPlasKin_set_config(ATOL=atol_in, RTOL=rtol_in, HMAX=hmax_in, MXSTEP=mxstep_in)
  else
    call ZDPlasKin_set_config(ATOL=atol_in, RTOL=rtol_in, MXSTEP=mxstep_in)
  endif
  ntot = 101325.0d0 / (1.38064852d-17 * Tgas)
  ne0 = 1.17d8
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=Tgas, REDUCED_FIELD=EN)
  call ZDPlasKin_set_density('N2', n2_frac * ntot)
  call ZDPlasKin_set_density('H2', (1.0d0 - n2_frac) * ntot)
  ! Hong 2017/2018 Table-5 geometry: 1e15 cm^-2 / (V/A = 0.007 cm).
  call ZDPlasKin_set_density('Surf', 1.4285714285714286d17)
  call ZDPlasKin_set_density('E', ne0, ldens_const=.true.)
  u = 20
  open(unit=u, file=series_csv, status='replace', iostat=ios)
  if (ios .ne. 0) stop 'cannot open concentration output'
  write(u,'(A)') 'time_s,N2_cm-3,H2_cm-3,NH3_cm-3,NH2_cm-3,NH_cm-3,N_cm-3,H_cm-3,E_cm-3,theta_Surf,theta_HSurf,theta_NSurf,theta_NHSurf,theta_NH2Surf,site_balance_error'
  time = 0.0d0
  samples = 0
  call write_sample(u, time)
  samples = samples + 1
  ! Long-horizon diagnostic studies retain logarithmic early outputs.  The
  ! fixed-exposure grid needs only its terminal 80--100 s window, so it avoids
  ! externally restarting DVODE throughout the stiffest initial transient.
  if (early_diagnostics_flag .eq. 1) then
    next_report = min(1.0d-4, report_dt)
  else
    next_report = report_dt
  endif

  do while (time .lt. t_end - 1.0d-10 * max(1.0d0, t_end))
    ! DVODE still chooses its own adaptive substeps.  Positive targets prevent
    ! the expensive one-internal-step bootstrap loop while retaining decade
    ! resolved early-time observations.
    target_window = window_dt
    ! At 60 Td the initial Hong surface/gas burst cannot be advanced directly
    ! to a 1 s outer target with DVODE's default MXSTEP.  A single 0.1 s
    ! startup target crosses that burst; later targets retain window_dt.
    if (early_diagnostics_flag .eq. 0 .and. time .lt. 1.0d-1) target_window = min(target_window, 1.0d-1)
    ! The high-MXSTEP branch is used only for the 60 Td stiff column.  Give
    ! its final 10 percent the same short external targets so no single
    ! 1-second request exhausts the internal step budget near t_end.
    if (early_diagnostics_flag .eq. 0 .and. mxstep_in .gt. 5000 .and. time .ge. 0.9d0*t_end) &
      target_window = min(target_window, 1.0d-1)
    dt_step = min(target_window, t_end - time, next_report - time)
    if (dt_step .le. 1.0d-14 * max(1.0d0, time)) stop 'non-positive outer target step'
    call ZDPlasKin_timestep(time, dt_step)
    time = time + dt_step
    if (time + 1.0d-10 * max(1.0d0, t_end) .ge. next_report .or. &
        time + 1.0d-10 * max(1.0d0, t_end) .ge. t_end) then
      if (site_projection_flag .eq. 1) call project_surface_sites()
      call write_sample(u, time)
      samples = samples + 1
      do while (next_report .le. time + 1.0d-10 * max(1.0d0, t_end))
        if (next_report .lt. report_dt) then
          next_report = min(report_dt, 1.0d1 * next_report)
        else
          next_report = next_report + report_dt
        endif
      end do
    end if
  end do
  close(u)
  write(*,'(A,A)') 'DONE tag=', trim(tag)
  write(*,'(A,ES13.5)') 't_end_s=', time
  write(*,'(A,I0)') 'samples=', samples
  write(*,'(A,ES13.5)') 'site_projection_max_rel=', site_projection_max_rel

contains
  subroutine write_sample(unit, t)
    integer, intent(in) :: unit
    double precision, intent(in) :: t
    call ZDPlasKin_get_density('N2', n2)
    call ZDPlasKin_get_density('H2', h2)
    call ZDPlasKin_get_density('NH3', nh3)
    call ZDPlasKin_get_density('NH2', nh2)
    call ZDPlasKin_get_density('NH', nh)
    call ZDPlasKin_get_density('N', n_atom)
    call ZDPlasKin_get_density('H', h_atom)
    call ZDPlasKin_get_density('E', elec)
    call ZDPlasKin_get_density('Surf', surf)
    call ZDPlasKin_get_density('HSurf', hsurf)
    call ZDPlasKin_get_density('NSurf', nsurf)
    call ZDPlasKin_get_density('NHSurf', nhsurf)
    call ZDPlasKin_get_density('NH2Surf', nh2surf)
    sites = surf + hsurf + nsurf + nhsurf + nh2surf
    if (sites .le. 0.0d0) stop 'surface site inventory is not positive'
    write(unit,'(15(ES18.10,:,","))') t, n2, h2, nh3, nh2, nh, n_atom, h_atom, elec, &
      surf/sites, hsurf/sites, nsurf/sites, nhsurf/sites, nh2surf/sites, &
      abs(sites-1.4285714285714286d17)/1.4285714285714286d17
    call flush(unit)
  end subroutine write_sample

  subroutine project_surface_sites()
    double precision :: raw_sites, scale, correction
    call ZDPlasKin_get_density('Surf', surf)
    call ZDPlasKin_get_density('HSurf', hsurf)
    call ZDPlasKin_get_density('NSurf', nsurf)
    call ZDPlasKin_get_density('NHSurf', nhsurf)
    call ZDPlasKin_get_density('NH2Surf', nh2surf)
    raw_sites = surf + hsurf + nsurf + nhsurf + nh2surf
    if (raw_sites .le. 0.0d0) stop 'surface site inventory is not positive for projection'
    correction = abs(raw_sites-1.4285714285714286d17)/1.4285714285714286d17
    site_projection_max_rel = max(site_projection_max_rel, correction)
    scale = 1.4285714285714286d17 / raw_sites
    call ZDPlasKin_set_density('Surf', surf * scale)
    call ZDPlasKin_set_density('HSurf', hsurf * scale)
    call ZDPlasKin_set_density('NSurf', nsurf * scale)
    call ZDPlasKin_set_density('NHSurf', nhsurf * scale)
    call ZDPlasKin_set_density('NH2Surf', nh2surf * scale)
  end subroutine project_surface_sites
end program main_cw_surface
'''


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def assert_active_surface_model(base: Path) -> None:
    """Verify the selected Hong 2017/2018 corrected surface mechanism."""
    source = base / KINETIC_INPUT
    if not source.is_file() or not (base / "preprocessor.exe").is_file():
        raise RuntimeError("缺少 Hong 2017/2018 校正机理输入或 ZDPlasKin 预处理器")
    text = source.read_text(encoding="utf-8", errors="replace")
    active_surface = [line for line in text.splitlines()
                      if "=>" in line and "surf" in line.lower() and not line.lstrip().startswith("#")]
    # The corrected input has S1--S18 expanded across vibrational/electronic
    # states: 37 active heterogeneous reaction lines in this representation.
    if len(active_surface) != 37:
        raise RuntimeError(f"Hong 校正机理的活性表面反应数应为 37，实际为 {len(active_surface)}")
    required = (
        "N + Surf => NSurf", "H + Surf => HSurf", "NSurf + HSurf => NHSurf + Surf",
        "NH2Surf + HSurf => NH3 + 2Surf", "N2 + Surf + Surf => NSurf + NSurf",
        "H2 + Surf + Surf => HSurf + HSurf",
    )
    missing = [reaction for reaction in required if reaction not in text]
    if missing:
        raise RuntimeError("Hong 校正机理缺少关键表面反应：" + "; ".join(missing))
    if not DVODE_SOURCE.is_file():
        raise RuntimeError(f"缺少 Hong 工程配套的 DVODE 源文件：{DVODE_SOURCE}")
    if not BOLSIG_LIBRARY_SOURCE.is_file():
        raise RuntimeError(f"缺少 Hong 工程配套的 BOLSIG 链接库：{BOLSIG_LIBRARY_SOURCE}")
    if (not BOLSIG_DATABASE_SOURCE.is_file()
            or "H2 -> H2(RYDBERG_SUM)" not in BOLSIG_DATABASE_SOURCE.read_text(encoding="utf-8", errors="replace")
            or "N2 -> N2(v1res)" not in BOLSIG_DATABASE_SOURCE.read_text(encoding="utf-8", errors="replace")):
        raise RuntimeError(f"缺少与 Hong 校正机理兼容的 BOLSIG 数据库：{BOLSIG_DATABASE_SOURCE}")


def ignore_base(directory: str, names: list[str]) -> set[str]:
    ignored = {"__pycache__", "runs"}
    for name in names:
        lowered = name.lower()
        if lowered.endswith((".o", ".mod", ".exe")) or name == "Const_E.F90":
            ignored.add(name)
    return ignored


HONG_THERMAL_VELOCITY_MARKER = "# HONG_CW_WORKTREE_THERMAL_VELOCITY_PATCH_V1"


def patch_hong_surface_thermal_velocities(kinet_path: Path) -> None:
    """Repair a declaration omission in the copied corrected Hong input.

    The corrected input invokes the atomic/radical thermal velocities in its
    Hong Table-5 adsorption and Eley--Rideal coefficients but declares only
    N2/H2 velocities.  The definitions below are the standard Maxwell mean
    speed expressions used in the same Hong project's ``build_full/kinet.inp``
    (declaration lines 1534 and formulas 1548--1551).
    """
    text = kinet_path.read_text(encoding="utf-8", errors="strict")
    if HONG_THERMAL_VELOCITY_MARKER in text:
        return
    declaration_anchor = "$ double precision :: STVOL, DSURF, EZ, KADS_N, KADS_H, KADS_NH, KADS_NH2"
    calculation_anchor = "$ STVOL = SITE_DENSITY / V_OVER_A"
    declarations = """# HONG_CW_WORKTREE_THERMAL_VELOCITY_PATCH_V1
# Missing atomic/radical Maxwell mean speeds used by Table-5 KADS/KER terms.
$ double precision, parameter :: H_MASS_KG = 1.6735575d-27, N_MASS_KG = 2.3258671d-26
$ double precision, parameter :: NH_MASS_KG = H_MASS_KG + N_MASS_KG, NH2_MASS_KG = 2.d0*H_MASS_KG + N_MASS_KG
$ double precision :: H_THERMAL_VEL, N_THERMAL_VEL, NH_THERMAL_VEL, NH2_THERMAL_VEL
"""
    calculations = """$ H_THERMAL_VEL = sqrt(8.d0*KB*Tgas/(H_MASS_KG*MYPI))*100.d0
$ N_THERMAL_VEL = sqrt(8.d0*KB*Tgas/(N_MASS_KG*MYPI))*100.d0
$ NH_THERMAL_VEL = sqrt(8.d0*KB*Tgas/(NH_MASS_KG*MYPI))*100.d0
$ NH2_THERMAL_VEL = sqrt(8.d0*KB*Tgas/(NH2_MASS_KG*MYPI))*100.d0
"""
    if text.count(declaration_anchor) != 1 or text.count(calculation_anchor) != 1:
        raise RuntimeError("无法定位 Hong 校正机理的热速度补丁锚点")
    text = text.replace(declaration_anchor, declarations + declaration_anchor, 1)
    text = text.replace(calculation_anchor, calculations + calculation_anchor, 1)
    kinet_path.write_text(text, encoding="utf-8", newline="\n")


def patch_hong_bolsig_mixture(kinet_path: Path) -> None:
    """Keep H2 Rydberg in reactions, but exclude it from the EEDF gas mixture."""
    marker = "# HONG_CW_WORKTREE_BOLSIG_MIXTURE_PATCH_V1"
    text = kinet_path.read_text(encoding="utf-8", errors="strict")
    if marker in text:
        return
    bolsig_start = text.find("BOLSIG\n")
    bolsig_end = text.find("\nEND", bolsig_start)
    if bolsig_start < 0 or bolsig_end < 0:
        raise RuntimeError("无法定位 Hong 输入文件的 BOLSIG 组分段")
    block = text[bolsig_start:bolsig_end]
    target = "H2(B3SIG) H2(B1SIG) H2(C3PI) H2(A3SIG) H2(RYDBERG_SUM)"
    if target not in block:
        raise RuntimeError("Hong 校正输入的 Rydberg BOLSIG 组分行与预期不符")
    block = block.replace(target, "H2(B3SIG) H2(B1SIG) H2(C3PI) H2(A3SIG)", 1)
    text = text[:bolsig_start] + block + text[bolsig_end:]
    text = text.replace("BOLSIG\n", f"BOLSIG\n{marker}\n", 1)
    kinet_path.write_text(text, encoding="utf-8", newline="\n")


def instrument_module_with_hmax(module_path: Path) -> None:
    """Expose DVODE HMAX and MXSTEP only in the copied simulation worktree.

    The reproduced surface module in ``BASE`` is retained byte-for-byte as a
    scientific reference.  This controlled patch extends its public config
    wrapper so the convergence study can constrain its internal step and
    internal-step budget without changing the reference mechanism.
    """
    source = module_path.read_text(encoding="utf-8", errors="strict")
    signature = (
        "subroutine ZDPlasKin_set_config(ATOL,RTOL,SILENCE_MODE,STAT_ACCUM,QTPLASKIN_SAVE,"
        "BOLSIG_EE_FRAC,BOLSIG_IGNORE_GAS_TEMPERATURE)"
    )
    declaration = "double precision, optional, intent(in) :: ATOL, RTOL, BOLSIG_EE_FRAC"
    locals_line = "double precision :: atol_loc, rtol_loc"
    saves_line = "double precision, save :: atol_save = -1.0d0, rtol_save = -1.0d0"
    configure_block = """  if(atol_loc/=atol_save .or. rtol_loc/=rtol_save) then
    atol_save = atol_loc
    rtol_save = rtol_loc
    if( lprint ) write(*,\"(2(A,1pd9.2),A)\") \"ZDPlasKin INFO: set accuracy\", atol_save, \" (absolute) &\", rtol_save, \" (relative)\"
    dens_loc(:,0) = 0.0d0
    dens_loc(:,1) = huge(dens_loc)
    vode_options  = set_intermediate_opts(abserr=atol_save,relerr=rtol_save, &
                                          dense_j=.true.,user_supplied_jacobian=.true., &
                                          constrained=bounded_components(:),clower=dens_loc(:,0),cupper=dens_loc(:,1))
    if(vode_istate /= 1) vode_istate = 3
  endif"""
    replacement_block = """  if( present(HMAX) ) then
    hmax_loc = HMAX
  else
    hmax_loc = hmax_save
  endif
  if( present(MXSTEP) ) then
    mxstep_loc = MXSTEP
  else
    mxstep_loc = mxstep_save
  endif
  if(hmax_loc == 0.0d0 .or. hmax_loc < -1.0d0) &
    call ZDPlasKin_stop(\"ZDPlasKin ERROR: HMAX must be positive when specified (ZDPlasKin_set_config)\")
  if(mxstep_loc <= 0) &
    call ZDPlasKin_stop(\"ZDPlasKin ERROR: MXSTEP must be positive (ZDPlasKin_set_config)\")
  if(atol_loc/=atol_save .or. rtol_loc/=rtol_save .or. hmax_loc/=hmax_save .or. mxstep_loc/=mxstep_save) then
    atol_save = atol_loc
    rtol_save = rtol_loc
    hmax_save = hmax_loc
    mxstep_save = mxstep_loc
    if( lprint ) write(*,\"(2(A,1pd9.2),A)\") \"ZDPlasKin INFO: set accuracy\", atol_save, \" (absolute) &\", rtol_save, \" (relative)\"
    dens_loc(:,0) = 0.0d0
    dens_loc(:,1) = huge(dens_loc)
    if(hmax_save > 0.0d0) then
      vode_options = set_intermediate_opts(abserr=atol_save,relerr=rtol_save,hmax=hmax_save,mxstep=mxstep_save, &
                                            dense_j=.true.,user_supplied_jacobian=.true., &
                                            constrained=bounded_components(:),clower=dens_loc(:,0),cupper=dens_loc(:,1))
    else
      vode_options = set_intermediate_opts(abserr=atol_save,relerr=rtol_save,mxstep=mxstep_save, &
                                            dense_j=.true.,user_supplied_jacobian=.true., &
                                            constrained=bounded_components(:),clower=dens_loc(:,0),cupper=dens_loc(:,1))
    endif
    if(vode_istate /= 1) vode_istate = 3
  endif"""
    anchors = (signature, declaration, locals_line, saves_line, configure_block)
    if any(source.count(anchor) != 1 for anchor in anchors):
        raise RuntimeError("无法对复制的 ZDPlasKin 模块施加 HMAX 控制：源锚点不唯一或缺失")
    source = source.replace(signature, signature[:-1] + ",HMAX,MXSTEP)")
    source = source.replace(declaration, declaration + ", HMAX")
    source = source.replace(declaration + ", HMAX", declaration + ", HMAX\n  integer, optional, intent(in) :: MXSTEP")
    source = source.replace(locals_line, "double precision :: atol_loc, rtol_loc, hmax_loc\n  integer :: mxstep_loc")
    source = source.replace(saves_line, saves_line + ", hmax_save = -1.0d0\n  integer, save :: mxstep_save = 5000")
    source = source.replace(configure_block, replacement_block)
    module_path.write_text(source, encoding="utf-8", newline="\n")


def prepare_worktree(root: Path) -> Path:
    assert_active_surface_model(BASE)
    work = root / "worktree_cw_surface"
    if work.is_dir():
        driver = work / "main_cw_surface.F90"
        exe = work / "main_cw_surface.exe"
        if driver.is_file() and exe.is_file():
            return work
        raise RuntimeError(f"工作目录已存在但不完整：{work}；请使用新的输出目录")
    if root.exists():
        permitted = {".claude_resources.json"}
        unexpected = [item.name for item in root.iterdir() if item.name not in permitted]
        if unexpected:
            raise RuntimeError(f"输出根目录已存在且包含非资源报告文件：{root}；请使用新的输出目录")
    else:
        root.mkdir(parents=True)
    shutil.copytree(BASE, work, ignore=ignore_base)
    # The corrected Hong archive is input-focused.  Use the DVODE source from
    # the same Hong reproduction tree, then regenerate the module from the
    # selected corrected input in this isolated copy.
    shutil.copy2(DVODE_SOURCE, work / "dvode_f90_m.F90")
    shutil.copy2(BASE / "preprocessor.exe", work / "preprocessor.exe")
    shutil.copy2(BOLSIG_LIBRARY_SOURCE, work / "bolsig_x86_64_g.lib")
    shutil.copy2(BOLSIG_DATABASE_SOURCE, work / "bolsigdb.dat")
    shutil.copy2(work / KINETIC_INPUT, work / "kinet.inp")
    patch_hong_surface_thermal_velocities(work / "kinet.inp")
    patch_hong_bolsig_mixture(work / "kinet.inp")
    hmax_config = """  if (hmax_in .gt. 0.0d0) then
    call ZDPlasKin_set_config(ATOL=atol_in, RTOL=rtol_in, HMAX=hmax_in, MXSTEP=mxstep_in)
  else
    call ZDPlasKin_set_config(ATOL=atol_in, RTOL=rtol_in, MXSTEP=mxstep_in)
  endif"""
    bootstrap_driver = CW_SURFACE_DRIVER.replace(hmax_config, "  call ZDPlasKin_set_config(ATOL=atol_in, RTOL=rtol_in)")
    if bootstrap_driver == CW_SURFACE_DRIVER:
        raise RuntimeError("无法生成 Hong 机制的 HMAX 两阶段编译驱动")
    # Stage 1: preprocess the selected input and build its unmodified module.
    (work / "main_cw_surface.F90").write_text(bootstrap_driver, encoding="utf-8", newline="\n")
    rc = runtime.build(work, mode="full")
    if rc:
        raise RuntimeError(f"Hong 校正机理的基础 CW 驱动器编译失败，退出码 {rc}")
    # Stage 2: add HMAX only in this copied module and rebuild the formal driver.
    instrument_module_with_hmax(work / "zdplaskin_m.F90")
    (work / "main_cw_surface.F90").write_text(CW_SURFACE_DRIVER, encoding="utf-8", newline="\n")
    rc = runtime.build(work, mode="main")
    if rc:
        raise RuntimeError(f"Hong 校正机理的 HMAX CW 驱动器编译失败，退出码 {rc}")
    manifest = {
        "purpose": "CW DBD gas--surface convergence and parameter scan using Hong 2017/2018 corrected mechanism",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "base_directory": str(BASE),
        "base_kinet_source": str(BASE / KINETIC_INPUT),
        "base_kinet_source_sha256": sha256(BASE / KINETIC_INPUT),
        "base_generated_module_sha256": sha256(BASE / "zdplaskin_m.F90"),
        "bolsig_database": str(BOLSIG_DATABASE_SOURCE),
        "bolsig_database_sha256": sha256(BOLSIG_DATABASE_SOURCE),
        "worktree_generated_module_sha256": sha256(work / "zdplaskin_m.F90"),
        "surface_model": "Hong 2017/2018 corrected surface mechanism, Table-5 metal column",
        "surface_site_density_cm-3": SURFACE_SITE_DENSITY_CM3,
        "surface_geometry": {"V_over_A_cm": 7.0e-3, "diffusion_length_cm": 1.0e-2,
                             "site_density_cm-2": 1.0e15},
        "electron_density_cm-3": NE_CM3,
        "reactor_boundary": "closed 0D; no feed, outlet, residence-time, or CSTR branch",
        "forcing": "continuous reduced field; no physical pulse waveform",
        "driver": "main_cw_surface.F90",
        "numerical_control": "DVODE relative/absolute tolerances plus worktree-only HMAX extension",
        "worktree_input_patches": [
            "added missing N/H/NH/NH2 thermal-velocity declarations from the same Hong build_full implementation",
            "excluded H2(RYDBERG_SUM) from BOLSIG EEDF mixture while retaining all Rydberg reactions",
        ],
    }
    (root / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return work


def token(value: float) -> str:
    return f"{value:g}".replace(".", "p").replace("-", "m")


def case_name(tg: float, en: float, n2: float, dt: float, *, hmax: float | None = None,
              rtol: float | None = None) -> str:
    tag = f"T{token(tg)}_EN{token(en)}_N2{token(n2)}_dt{token(dt)}s"
    if hmax is not None and hmax > 0.0:
        tag += f"_hmax{token(hmax)}s"
    if rtol is not None:
        tag += f"_rtol{token(rtol)}"
    return tag


def parse_fortran_float(value: str) -> float:
    """Parse standard and exponent-letter-omitted Fortran real literals."""
    text = value.strip()
    try:
        return float(text)
    except ValueError:
        match = re.fullmatch(r"([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d{2,3})", text)
        if not match:
            raise
        return float(f"{match.group(1)}e{match.group(2)}")


def solver_completion_evidence(stdout: str, tag: str, requested_t_end_s: float) -> dict[str, object]:
    """Extract the explicit normal-completion evidence emitted by the driver.

    The CSV time column follows the outer driver target. A usable solver
    result must also print its final ``DONE`` marker and matching ``t_end_s``.
    A DVODE ``T + H = T`` warning remains diagnostic metadata: it is not, by
    itself, a failed solve if the driver subsequently completes normally.
    """
    marker = f"DONE tag={tag}" in stdout
    reported: float | None = None
    matches = re.findall(r"(?im)^\s*t_end_s=\s*([^\s]+)", stdout)
    if matches:
        try:
            reported = parse_fortran_float(matches[-1])
        except ValueError:
            reported = None
    tolerance = max(1.0e-8, abs(requested_t_end_s) * 1.0e-6)
    endpoint_pass = reported is not None and abs(reported - requested_t_end_s) <= tolerance
    projection_max_rel: float | None = None
    projection_matches = re.findall(r"(?im)^\s*site_projection_max_rel=\s*([^\s]+)", stdout)
    if projection_matches:
        try:
            projection_max_rel = parse_fortran_float(projection_matches[-1])
        except ValueError:
            projection_max_rel = None
    return {
        "done_marker_present": marker,
        "reported_t_end_s": reported,
        "requested_t_end_s": requested_t_end_s,
        "endpoint_tolerance_s": tolerance,
        "pass": bool(marker and endpoint_pass),
        "dvode_step_size_warning_count": stdout.count("T + H = T"),
        "dvode_error_marker_present": "DVODE solver issued an error" in stdout,
        "site_projection_max_rel": projection_max_rel,
    }


def read_series(path: Path) -> list[dict[str, float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [{key: parse_fortran_float(value) for key, value in row.items()} for row in csv.DictReader(handle)]


def relative_range(values: list[float], floor: float) -> float:
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
    result: dict[str, float] = {"samples": float(len(rows)), "t_final_s": last_t, "terminal_window_s": terminal_window}
    gas = ("NH3_cm-3", "NH2_cm-3", "NH_cm-3", "N_cm-3", "H_cm-3")
    coverages = ("theta_Surf", "theta_HSurf", "theta_NSurf", "theta_NHSurf", "theta_NH2Surf")
    for field in gas:
        result[field] = rows[-1][field]
        result[f"{field}_tail_rel_range"] = relative_range([row[field] for row in tail], 1.0)
    for field in coverages:
        result[field] = rows[-1][field]
        result[f"{field}_tail_rel_range"] = relative_range([row[field] for row in tail], 1.0e-30)
    result["site_balance_error_max"] = max(abs(row["site_balance_error"]) for row in rows)
    result["NH3_tail_pass"] = float(result["NH3_cm-3_tail_rel_range"] <= NH3_DRIFT_LIMIT)
    result["gas_tail_pass"] = float(max(result[f"{field}_tail_rel_range"] for field in gas[1:]) <= GAS_DRIFT_LIMIT)
    result["coverage_tail_pass"] = float(max(result[f"{field}_tail_rel_range"] for field in coverages) <= COVERAGE_DRIFT_LIMIT)
    result["site_balance_pass"] = float(result["site_balance_error_max"] <= SITE_BALANCE_LIMIT)
    result["acceptance_pass"] = float(all(result[key] for key in ("NH3_tail_pass", "gas_tail_pass", "coverage_tail_pass", "site_balance_pass")))
    return result


def run_case(work: Path, root: Path, tg: float, en: float, n2: float, dt: float,
             t_end: float, report_dt: float, group: str, *, atol: float = 1.0,
             rtol: float = 1.0e-4, internal_hmax_s: float | None = None,
             early_diagnostics: bool = True,
             internal_mxstep: int = 5000,
             site_projection: bool = False,
             residence_time_s: float | None = None,
             timeout_s: float = 3600.0,
             case_output_dir: Path | None = None,
             progress_callback: Callable[[dict[str, object]], None] | None = None) -> dict[str, object]:
    """Run one isolated physical case.

    ``case_output_dir`` and ``progress_callback`` are optional campaign-layer
    hooks.  Existing studies keep their original directory layout and blocking
    behavior when neither is supplied.
    """
    tag = case_name(tg, en, n2, dt, hmax=internal_hmax_s, rtol=rtol)
    out = case_output_dir if case_output_dir is not None else root / group / tag
    series = out / f"cw_surface_{tag}.csv"
    params = out / "params.json"
    if out.exists():
        if series.is_file():
            if params.is_file():
                record = json.loads(params.read_text(encoding="utf-8"))
                if record.get("timed_out") or record.get("returncode") not in (0, None):
                    return record
            else:
                record = {
                    "tag": tag, "group": group, "tg_K": tg, "en_Td": en, "n2_fraction": n2,
                    "h2_fraction": 1.0 - n2, "window_dt_s": dt, "t_end_requested_s": t_end,
                    "report_interval_s": report_dt, "ne_cm-3": NE_CM3,
                    "surface_site_density_cm-3": SURFACE_SITE_DENSITY_CM3, "forcing": "CW",
                    "pulse_enabled": False, "returncode": 0, "timed_out": False,
                    "atol": atol, "rtol": rtol, "timeout_s": timeout_s,
                    "internal_hmax_s": internal_hmax_s,
                    "early_diagnostics": early_diagnostics,
                    "internal_mxstep": internal_mxstep,
                    "site_projection": site_projection,
                    "recovered_from_completed_csv": True,
                }
            record.update(summarize_series(read_series(series), t_end))
            params.write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
            return record
        # The campaign controller atomically creates ``state.json`` before it
        # launches a solver, so an interruption can always be recovered.  A
        # directory containing that marker and nothing else is a new attempt,
        # not a partial physical output.  Keep the stricter legacy behavior
        # for standalone calls and for any other pre-existing artifact.
        allowed_campaign_markers = {"state.json"}
        existing_names = {item.name for item in out.iterdir()}
        if case_output_dir is None or not existing_names.issubset(allowed_campaign_markers):
            raise RuntimeError(f"目标目录已存在但缺少完整产物：{out}")
    else:
        out.mkdir(parents=True)
    # Runtime inputs used by the Hong 2017/2018 corrected mechanism.  The
    # entropy/activation-energy side files belong to the unrelated Shao model.
    for name in ("bolsigdb.dat", "bolsig_x86_64_g.dll"):
        source = work / name
        if not source.is_file():
            raise RuntimeError(f"工作目录缺少所需输入：{source}")
        shutil.copy2(source, out / name)
    exe = work / "main_cw_surface.exe"
    command = [str(exe), f"{n2:.12g}", f"{tg:.12g}", f"{en:.12g}", f"{dt:.12g}",
               f"{t_end:.12g}", f"{report_dt:.12g}", tag, f"{atol:.12g}", f"{rtol:.12g}",
               f"{internal_hmax_s if internal_hmax_s is not None else -1.0:.12g}",
               "1" if early_diagnostics else "0", str(internal_mxstep),
               "1" if site_projection else "0"]
    if residence_time_s is not None:
        if residence_time_s <= 0.0:
            raise ValueError("residence_time_s 必须为正数")
        command.append(f"{residence_time_s:.12g}")
    started = datetime.now().isoformat(timespec="seconds")
    env = runtime._runtime_env(runtime.resolve_gfortran_dir())
    # ZDPlasKin asks for ENTER after a fatal DVODE return.  A detached campaign
    # has nobody to press it; DEVNULL supplies EOF so the solver exits promptly
    # and the ledger records the true numerical failure rather than a later
    # one-hour orchestration timeout.
    # Do not use a PIPE here.  High-stiffness DVODE trajectories may emit many
    # warnings; a parent that only consumes stdout after process exit can fill
    # the Windows pipe buffer and deadlock an otherwise finished solver.  A
    # per-case log file streams those diagnostics without sharing output files
    # between concurrent cases.
    console_path = out / "console.log"
    with console_path.open("w", encoding="utf-8", newline="") as console_handle:
        proc = subprocess.Popen(command, cwd=out, env=env, stdin=subprocess.DEVNULL,
                                stdout=console_handle, stderr=subprocess.STDOUT,
                                text=True, encoding="utf-8", errors="replace")
        if progress_callback is not None:
            progress_callback({"event": "started", "pid": proc.pid, "tag": tag,
                               "output_dir": str(out), "timestamp": started})
        try:
            deadline = time.monotonic() + timeout_s
            next_heartbeat = time.monotonic() + 5.0
            while proc.poll() is None:
                now = time.monotonic()
                if now >= deadline:
                    proc.kill()
                    proc.wait()
                    rc, timed_out = None, True
                    break
                if progress_callback is not None and now >= next_heartbeat:
                    progress_callback({"event": "heartbeat", "pid": proc.pid, "tag": tag,
                                       "output_dir": str(out),
                                       "timestamp": datetime.now().isoformat(timespec="seconds")})
                    next_heartbeat = now + 5.0
                time.sleep(min(0.5, max(0.01, deadline - now)))
            else:
                proc.wait()
                rc, timed_out = proc.returncode, False
        except BaseException:
            # Do not hide an orchestration interruption.  The campaign layer
            # has already recorded this attempt as running and can resume it.
            if proc.poll() is None:
                proc.kill()
                proc.wait()
            raise
    console_text = console_path.read_text(encoding="utf-8", errors="replace")
    completion = solver_completion_evidence(console_text, tag, t_end)
    record: dict[str, object] = {
        "tag": tag, "group": group, "tg_K": tg, "en_Td": en, "n2_fraction": n2,
        "h2_fraction": 1.0 - n2, "window_dt_s": dt, "t_end_requested_s": t_end,
        "report_interval_s": report_dt, "ne_cm-3": NE_CM3,
        "surface_site_density_cm-3": SURFACE_SITE_DENSITY_CM3, "forcing": "CW",
        "pulse_enabled": False, "returncode": rc, "timed_out": timed_out,
        "atol": atol, "rtol": rtol, "timeout_s": timeout_s,
        "internal_hmax_s": internal_hmax_s,
        "early_diagnostics": early_diagnostics,
        "internal_mxstep": internal_mxstep,
        "site_projection": site_projection,
        "residence_time_s": residence_time_s,
        "started_at": started, "ended_at": datetime.now().isoformat(timespec="seconds"),
        "solver_completion": completion,
    }
    if rc == 0 and not timed_out and series.is_file() and bool(completion["pass"]):
        record.update(summarize_series(read_series(series), t_end))
    else:
        if timed_out:
            record["failure_reason"] = "timeout"
        elif rc != 0:
            record["failure_reason"] = f"process return code {rc}"
        elif not bool(completion["pass"]):
            record["failure_reason"] = "missing or inconsistent solver completion marker"
        else:
            record["failure_reason"] = "missing concentration CSV"
        if rc == 0 and not series.is_file():
            record["failure_reason"] = "missing concentration CSV"
    params.write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
    return record


def write_csv(path: Path, records: list[dict[str, object]]) -> None:
    fields = sorted({key for record in records for key in record})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def row_at_time(rows: list[dict[str, float]], time_s: float) -> dict[str, float]:
    """Return an explicitly reported state at a requested physical exposure time."""
    row = min(rows, key=lambda candidate: abs(candidate["time_s"] - time_s))
    if abs(row["time_s"] - time_s) > max(1.0e-8, 1.0e-6 * time_s):
        raise RuntimeError(f"轨迹未包含所需暴露时间 {time_s:g} s（最近输出 {row['time_s']:g} s）")
    return row


def relative_drift(start: float, end: float, floor: float) -> float:
    return abs(end - start) / max(abs(end), floor)


def exposure_record(rows: list[dict[str, float]], base: dict[str, object], horizon_s: float) -> dict[str, object]:
    """Evaluate active-state quasi-steady behavior without requiring an NH3 plateau."""
    window_s = min(0.2 * horizon_s, 100.0)
    start = row_at_time(rows, horizon_s - window_s)
    middle = row_at_time(rows, horizon_s - 0.5 * window_s)
    end = row_at_time(rows, horizon_s)
    record: dict[str, object] = {
        "tag": base["tag"], "tg_K": base["tg_K"], "en_Td": base["en_Td"],
        "n2_fraction": base["n2_fraction"], "h2_fraction": base["h2_fraction"],
        "exposure_time_s": horizon_s, "evaluation_window_s": window_s,
        "nh3_cm-3": end["NH3_cm-3"],
        "nh3_net_rate_cm-3_s": (end["NH3_cm-3"] - start["NH3_cm-3"]) / window_s,
        "site_balance_error": end["site_balance_error"],
    }
    gas_drifts = []
    for field in ("N_cm-3", "H_cm-3"):
        drift = relative_drift(start[field], end[field], 1.0)
        record[f"{field}_window_rel_drift"] = drift
        gas_drifts.append(drift)
    coverage_drifts = []
    active_coverages = []
    for field in ("theta_Surf", "theta_HSurf", "theta_NSurf", "theta_NHSurf", "theta_NH2Surf"):
        relevant = max(abs(start[field]), abs(end[field])) >= SURFACE_COVERAGE_FLOOR
        record[f"{field}_relevant"] = float(relevant)
        if relevant:
            drift = relative_drift(start[field], end[field], SURFACE_COVERAGE_FLOOR)
            record[f"{field}_window_rel_drift"] = drift
            coverage_drifts.append(drift)
            active_coverages.append(field)
    rate_first_half = (middle["NH3_cm-3"] - start["NH3_cm-3"]) / (0.5 * window_s)
    rate_second_half = (end["NH3_cm-3"] - middle["NH3_cm-3"]) / (0.5 * window_s)
    rate_change = abs(rate_second_half - rate_first_half) / max(abs(rate_second_half), 1.0)
    record["nh3_rate_first_half_cm-3_s"] = rate_first_half
    record["nh3_rate_second_half_cm-3_s"] = rate_second_half
    record["nh3_rate_half_window_rel_change"] = rate_change
    record["reactive_gas_max_window_rel_drift"] = max(gas_drifts)
    record["major_coverage_max_window_rel_drift"] = max(coverage_drifts) if coverage_drifts else 0.0
    record["relevant_surface_coverages"] = ";".join(active_coverages)
    record["active_state_pass"] = float(
        record["reactive_gas_max_window_rel_drift"] <= ACTIVE_STATE_DRIFT_LIMIT
        and record["major_coverage_max_window_rel_drift"] <= ACTIVE_STATE_DRIFT_LIMIT
        and record["site_balance_error"] <= SITE_BALANCE_LIMIT
    )
    record["product_rate_pass"] = float(rate_change <= PRODUCT_RATE_STABILITY_LIMIT)
    record["hong_style_exposure_pass"] = float(record["active_state_pass"] and record["product_rate_pass"])
    return record


def hong_exposure_study(work: Path, root: Path, t_end: float, report_dt: float, workers: int) -> None:
    """Closed-0D Hong-style scan at six slow low-field composition-temperature corners."""
    horizons = [1.0, 10.0, 100.0, 1000.0]
    if t_end < max(horizons):
        raise ValueError("Hong 暴露时间扫描的 t_end 至少应为 1000 s")
    if t_end >= 9000.0:
        horizons.append(9000.0)
    conditions = [(float(tg), 20.0, n2) for tg in (300, 400, 500) for n2 in (0.1, 0.9)]

    def execute(case: tuple[float, float, float]) -> dict[str, object]:
        tg, en, n2 = case
        return run_case(work, root, tg, en, n2, 1.0, t_end, report_dt, f"hong_exposure_{token(t_end)}s",
                        atol=1.0, rtol=1.0e-4, internal_hmax_s=5.0e-3, timeout_s=7200.0)

    records: list[dict[str, object]] = []
    if workers == 1:
        for index, case in enumerate(conditions, 1):
            print(f"[hong-exposure {index}/{len(conditions)}] Tg={case[0]:g} K, E/N={case[1]:g} Td, xN2={case[2]:g}")
            records.append(execute(case))
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(execute, case): case for case in conditions}
            for index, future in enumerate(as_completed(futures), 1):
                tg, en, n2 = futures[future]
                print(f"[hong-exposure {index}/{len(conditions)}] Tg={tg:g} K, E/N={en:g} Td, xN2={n2:g}")
                records.append(future.result())
    write_csv(root / "hong_exposure_runs.csv", records)
    analysis: list[dict[str, object]] = []
    for base in records:
        if base.get("returncode") != 0 or base.get("timed_out"):
            continue
        tag = str(base["tag"])
        series = root / str(base["group"]) / tag / f"cw_surface_{tag}.csv"
        rows = read_series(series)
        analysis.extend(exposure_record(rows, base, horizon) for horizon in horizons)
    write_csv(root / "hong_exposure_time_summary.csv", analysis)
    (root / "hong_exposure_protocol.json").write_text(json.dumps({
        "surface_mechanism": "Hong 2017/2018 corrected, Table-5 metal column",
        "reactor_boundary": "closed 0D; no CSTR/feed/outlet/residence-time term",
        "forcing": "continuous E/N; no pulse waveform or cycle count",
        "representative_conditions": [{"Tgas_K": tg, "EoverN_Td": en, "x_N2": n2, "x_H2": 1.0 - n2}
                                      for tg, en, n2 in conditions],
        "integrated_to_s": t_end,
        "reported_exposure_times_s": horizons,
        "report_interval_s": report_dt,
        "solver": {"ATOL_cm-3": 1.0, "RTOL": 1.0e-4, "DVODE_HMAX_s": 5.0e-3},
        "acceptance": {
            "active_state_window_relative_drift": ACTIVE_STATE_DRIFT_LIMIT,
            "product_rate_half_window_relative_change": PRODUCT_RATE_STABILITY_LIMIT,
            "surface_coverage_relevance_floor": SURFACE_COVERAGE_FLOOR,
            "site_balance_error": SITE_BALANCE_LIMIT,
        },
        "product_interpretation": "NH3 is finite-exposure inventory plus terminal net rate, not a required steady concentration in closed 0D.",
    }, ensure_ascii=False, indent=2), encoding="utf-8")


def time_horizon_summary(root: Path, records: list[dict[str, object]], horizons: list[float]) -> None:
    """Compare stored time checkpoints to the final gas and surface state.

    This is a true integration-duration scan extracted from every completed
    long trajectory, so it does not repeat otherwise identical stiff solves.
    """
    fields = ("NH3_cm-3", "NH2_cm-3", "NH_cm-3", "N_cm-3", "H_cm-3",
              "theta_Surf", "theta_HSurf", "theta_NSurf", "theta_NHSurf", "theta_NH2Surf")
    summary: list[dict[str, object]] = []
    for record in records:
        if record.get("returncode") != 0 or record.get("timed_out"):
            continue
        tag = str(record["tag"])
        series = root / str(record["group"]) / tag / f"cw_surface_{tag}.csv"
        rows = read_series(series)
        final = rows[-1]
        for horizon in horizons:
            if horizon > final["time_s"] + 1.0e-8:
                continue
            row = min(rows, key=lambda candidate: abs(candidate["time_s"] - horizon))
            diffs = [abs(row[field] - final[field]) / max(abs(final[field]), 1.0e-30) for field in fields]
            summary.append({
                "tag": tag, "tg_K": record["tg_K"], "en_Td": record["en_Td"],
                "n2_fraction": record["n2_fraction"], "window_dt_s": record["window_dt_s"],
                "requested_horizon_s": horizon, "recorded_time_s": row["time_s"],
                "max_key_state_rel_diff_to_final": max(diffs),
            })
    write_csv(root / "time_horizon_summary.csv", summary)


def time_study(work: Path, root: Path, windows: list[float], t_end: float, report_dt: float,
               workers: int) -> None:
    conditions = [(tg, en, n2) for tg in (300.0, 400.0, 500.0)
                  for en, n2 in ((20.0, 0.1), (20.0, 0.9), (140.0, 0.5), (240.0, 0.1), (240.0, 0.9))]
    cases = [(tg, en, n2, dt) for tg, en, n2 in conditions for dt in windows]
    records = []
    if workers == 1:
        for tg, en, n2, dt in cases:
            print(f"[time-study] Tg={tg:g} K, E/N={en:g} Td, xN2={n2:g}, window={dt:g} s")
            records.append(run_case(work, root, tg, en, n2, dt, t_end, report_dt, "time_resolution"))
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(run_case, work, root, tg, en, n2, dt, t_end, report_dt, "time_resolution"):
                       (tg, en, n2, dt) for tg, en, n2, dt in cases}
            for index, future in enumerate(as_completed(futures), 1):
                tg, en, n2, dt = futures[future]
                print(f"[time-study {index}/{len(cases)}] Tg={tg:g} K, E/N={en:g} Td, xN2={n2:g}, window={dt:g} s")
                records.append(future.result())
    grouped: dict[tuple[object, object, object], list[dict[str, object]]] = defaultdict(list)
    for record in records:
        grouped[(record["tg_K"], record["en_Td"], record["n2_fraction"])].append(record)
    for group in grouped.values():
        good = [record for record in group if record.get("returncode") == 0 and not record.get("timed_out")]
        if not good:
            continue
        reference = min(good, key=lambda record: float(record["window_dt_s"]))
        compare_fields = ("NH3_cm-3", "theta_Surf", "theta_HSurf", "theta_NSurf", "theta_NHSurf", "theta_NH2Surf")
        for record in good:
            differences = [abs(float(record[field]) - float(reference[field])) / max(abs(float(reference[field])), 1.0e-30)
                           for field in compare_fields]
            record["max_resolution_rel_diff"] = max(differences)
            record["resolution_pass"] = float(record["max_resolution_rel_diff"] <= RESOLUTION_LIMIT)
    write_csv(root / "time_resolution_summary.csv", records)
    time_horizon_summary(root, records, [10.0, 25.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0, 4000.0])
    (root / "time_resolution_protocol.json").write_text(json.dumps({
        "temperature_K": [300, 400, 500],
        "composition_and_field": [{"E/N_Td": en, "x_N2": n2} for en, n2 in ((20, .1), (20, .9), (140, .5), (240, .1), (240, .9))],
        "window_dt_s": windows, "t_end_s": t_end, "report_interval_s": report_dt,
        "acceptance": {"NH3_tail_rel_range": NH3_DRIFT_LIMIT, "key_gas_tail_rel_range": GAS_DRIFT_LIMIT,
                       "surface_coverage_tail_rel_range": COVERAGE_DRIFT_LIMIT, "site_balance_error": SITE_BALANCE_LIMIT,
                       "resolution_rel_diff": RESOLUTION_LIMIT},
    }, ensure_ascii=False, indent=2), encoding="utf-8")


def low_field_study(work: Path, root: Path, windows: list[float], t_end: float, report_dt: float,
                    workers: int) -> None:
    """Resolution study dedicated to the 20 Td composition extremes."""
    conditions = [(float(tg), 20.0, n2) for tg in (300, 400, 500) for n2 in (0.1, 0.9)]
    cases = [(tg, en, n2, dt) for tg, en, n2 in conditions for dt in windows]
    records = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(run_case, work, root, tg, en, n2, dt, t_end, report_dt, "low_field_resolution"):
                   (tg, en, n2, dt) for tg, en, n2, dt in cases}
        for index, future in enumerate(as_completed(futures), 1):
            tg, en, n2, dt = futures[future]
            print(f"[low-field {index}/{len(cases)}] Tg={tg:g} K, xN2={n2:g}, window={dt:g} s")
            records.append(future.result())
    grouped: dict[tuple[object, object, object], list[dict[str, object]]] = defaultdict(list)
    for record in records:
        grouped[(record["tg_K"], record["en_Td"], record["n2_fraction"])].append(record)
    for group in grouped.values():
        good = [record for record in group if record.get("returncode") == 0 and not record.get("timed_out")]
        if not good:
            continue
        reference = min(good, key=lambda record: float(record["window_dt_s"]))
        fields = ("NH3_cm-3", "NH2_cm-3", "NH_cm-3", "N_cm-3", "H_cm-3",
                  "theta_Surf", "theta_HSurf", "theta_NSurf", "theta_NHSurf", "theta_NH2Surf")
        for record in good:
            record["max_resolution_rel_diff"] = max(
                abs(float(record[field]) - float(reference[field])) / max(abs(float(reference[field])), 1.0e-30)
                for field in fields
            )
            record["resolution_pass"] = float(record["max_resolution_rel_diff"] <= RESOLUTION_LIMIT)
    write_csv(root / "low_field_resolution_summary.csv", records)
    time_horizon_summary(root, records, [1.0, 10.0, 25.0, 50.0, 100.0])
    (root / "low_field_protocol.json").write_text(json.dumps({
        "field_Td": 20, "temperature_K": [300, 400, 500], "x_N2": [0.1, 0.9],
        "window_dt_s": windows, "t_end_s": t_end,
        "startup_method": "DVODE single adaptive steps through the first 0.2 s; every bootstrap state is written.",
        "next_stage": "For conditions failing the 100 s criterion, extend only those cases to 300 then 1000 s.",
    }, ensure_ascii=False, indent=2), encoding="utf-8")


def numerical_convergence_study(work: Path, root: Path, hmax_values: list[float], rtol_values: list[float],
                                t_end: float, report_dt: float, workers: int) -> None:
    """Scan true DVODE HMAX and tolerance controls at low-field composition extremes."""
    reference_hmax = min(hmax_values)
    reference_rtol = min(rtol_values)
    tolerance_hmax = hmax_values[len(hmax_values) // 2]
    settings = {(hmax, reference_rtol) for hmax in hmax_values}
    settings.update((tolerance_hmax, rtol) for rtol in rtol_values)
    conditions = [(float(tg), 20.0, n2) for tg in (300, 400, 500) for n2 in (0.1, 0.9)]
    cases = [(tg, en, n2, hmax, rtol) for tg, en, n2 in conditions for hmax, rtol in sorted(settings)]
    records: list[dict[str, object]] = []

    def execute(case: tuple[float, float, float, float, float]) -> dict[str, object]:
        tg, en, n2, hmax, rtol = case
        return run_case(work, root, tg, en, n2, 1.0, t_end, report_dt, "numerical_convergence",
                        atol=1.0, rtol=rtol, internal_hmax_s=hmax)

    if workers == 1:
        for index, case in enumerate(cases, 1):
            tg, _, n2, hmax, rtol = case
            print(f"[numerical {index}/{len(cases)}] Tg={tg:g} K, xN2={n2:g}, HMAX={hmax:g} s, RTOL={rtol:.0e}")
            records.append(execute(case))
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(execute, case): case for case in cases}
            for index, future in enumerate(as_completed(futures), 1):
                tg, _, n2, hmax, rtol = futures[future]
                print(f"[numerical {index}/{len(cases)}] Tg={tg:g} K, xN2={n2:g}, HMAX={hmax:g} s, RTOL={rtol:.0e}")
                records.append(future.result())

    fields = ("NH3_cm-3", "NH2_cm-3", "NH_cm-3", "N_cm-3", "H_cm-3",
              "theta_Surf", "theta_HSurf", "theta_NSurf", "theta_NHSurf", "theta_NH2Surf")
    grouped: dict[tuple[object, object, object], list[dict[str, object]]] = defaultdict(list)
    for record in records:
        grouped[(record["tg_K"], record["en_Td"], record["n2_fraction"])].append(record)
    for group in grouped.values():
        reference = next((record for record in group
                          if record.get("internal_hmax_s") == reference_hmax
                          and record.get("rtol") == reference_rtol
                          and record.get("returncode") == 0 and not record.get("timed_out")), None)
        if reference is None:
            continue
        for record in group:
            if record.get("returncode") != 0 or record.get("timed_out"):
                continue
            record["max_numerical_rel_diff"] = max(
                abs(float(record[field]) - float(reference[field])) / max(abs(float(reference[field])), 1.0e-30)
                for field in fields
            )
            record["numerical_reference"] = float(record is reference)
            record["numerical_convergence_pass"] = float(record["max_numerical_rel_diff"] <= RESOLUTION_LIMIT)
    write_csv(root / "numerical_convergence_summary.csv", records)
    (root / "numerical_convergence_protocol.json").write_text(json.dumps({
        "field_Td": 20, "temperature_K": [300, 400, 500], "x_N2": [0.1, 0.9],
        "outer_target_dt_s": 1.0, "t_end_s": t_end, "report_interval_s": report_dt,
        "internal_hmax_s": hmax_values, "rtol": rtol_values, "atol": 1.0,
        "reference": {"internal_hmax_s": reference_hmax, "rtol": reference_rtol},
        "comparison_fields": list(fields), "relative_difference_limit": RESOLUTION_LIMIT,
        "note": "HMAX constrains DVODE internal integration steps; outer targets are not treated as internal timesteps.",
    }, ensure_ascii=False, indent=2), encoding="utf-8")


def full_grid(work: Path, root: Path, dt: float, t_end: float, report_dt: float, workers: int) -> None:
    cases = [(float(tg), float(en), n2_int / 10)
             for tg in range(300, 501, 25)
             for n2_int in range(1, 10)
             for en in range(20, 241, 20)]
    records = []
    if workers == 1:
        for tg, en, n2 in cases:
            records.append(run_case(work, root, tg, en, n2, dt, t_end, report_dt, "grid_972"))
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(run_case, work, root, tg, en, n2, dt, t_end, report_dt, "grid_972"):
                       (tg, en, n2) for tg, en, n2 in cases}
            for index, future in enumerate(as_completed(futures), 1):
                tg, en, n2 = futures[future]
                print(f"[grid {index}/{len(cases)}] Tg={tg:g} K, E/N={en:g} Td, xN2={n2:g}")
                records.append(future.result())
    write_csv(root / "grid_972_summary.csv", records)


def hong_finite_exposure_grid(work: Path, root: Path, t_end: float, report_dt: float, workers: int) -> None:
    """Run the requested 9 x 9 x 12 Hong-mechanism grid at one fixed exposure time."""
    if abs(t_end - 100.0) > 1.0e-9:
        raise ValueError("本阶段固定使用 Hong 式 t_obs = 100 s；请勿混用不同有限暴露时间")
    if report_dt > 10.0:
        raise ValueError("为计算 80--100 s 的末窗生成速率，report_dt 不得大于 10 s")
    cases = [(float(tg), float(en), n2_int / 10)
             for tg in range(300, 501, 25)
             for n2_int in range(1, 10)
             for en in range(20, 241, 20)]

    def execute(case: tuple[float, float, float]) -> dict[str, object]:
        tg, en, n2 = case
        high_stiffness = en == 60.0
        return run_case(work, root, tg, en, n2, 1.0, t_end, report_dt, "hong_finite_exposure_grid_100s",
                        atol=1.0e8 if high_stiffness else 1.0, rtol=1.0e-4,
                        internal_hmax_s=5.0e-3, early_diagnostics=False,
                        internal_mxstep=100000 if high_stiffness else 5000, timeout_s=3600.0)

    records: list[dict[str, object]] = []
    if workers == 1:
        for index, case in enumerate(cases, 1):
            tg, en, n2 = case
            print(f"[hong-grid {index}/{len(cases)}] Tg={tg:g} K, E/N={en:g} Td, xN2={n2:g}")
            records.append(execute(case))
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(execute, case): case for case in cases}
            for index, future in enumerate(as_completed(futures), 1):
                tg, en, n2 = futures[future]
                print(f"[hong-grid {index}/{len(cases)}] Tg={tg:g} K, E/N={en:g} Td, xN2={n2:g}")
                records.append(future.result())
    write_csv(root / "hong_finite_exposure_grid_runs.csv", records)
    summary: list[dict[str, object]] = []
    for base in records:
        if base.get("returncode") != 0 or base.get("timed_out"):
            continue
        tag = str(base["tag"])
        series = root / str(base["group"]) / tag / f"cw_surface_{tag}.csv"
        summary.append(exposure_record(read_series(series), base, t_end))
    write_csv(root / "hong_finite_exposure_grid_972_summary.csv", summary)
    (root / "hong_finite_exposure_grid_protocol.json").write_text(json.dumps({
        "surface_mechanism": "Hong 2017/2018 corrected, Table-5 metal column",
        "reactor_boundary": "closed 0D; no CSTR/feed/outlet/residence-time term",
        "forcing": "continuous E/N; no pulse waveform or cycle count",
        "grid": {"Tgas_K": list(range(300, 501, 25)), "x_N2": [i / 10 for i in range(1, 10)],
                 "x_H2": [1.0 - i / 10 for i in range(1, 10)], "EoverN_Td": list(range(20, 241, 20)),
                 "case_count": len(cases)},
        "finite_exposure_time_s": t_end,
        "report_interval_s": report_dt,
        "early_diagnostics": False,
        "solver": {
            "default": {"ATOL_cm-3": 1.0, "RTOL": 1.0e-4, "DVODE_HMAX_s": 5.0e-3, "MXSTEP": 5000},
            "EoverN_60Td_high_stiffness": {"ATOL_cm-3": 1.0e8, "RTOL": 1.0e-4,
                                              "DVODE_HMAX_s": 5.0e-3, "MXSTEP": 100000},
        },
        "outputs": {
            "nh3_cm-3": "NH3 inventory at t_obs = 100 s",
            "nh3_net_rate_cm-3_s": "net NH3 rate averaged over 80--100 s",
            "active_state_metrics": "N/H and relevant surface-coverages drift over 80--100 s",
        },
        "interpretation": "Finite-exposure comparison, not a claim that all closed-0D species are at steady state.",
    }, ensure_ascii=False, indent=2), encoding="utf-8")


def hong_60td_stiffness_validation(work: Path, root: Path, t_end: float, report_dt: float,
                                   workers: int) -> None:
    """Validate the 60 Td solver controls at the temperature/composition bounds.

    The 60 Td column is the only one that requires relaxed absolute tolerances
    and a larger DVODE internal-step budget.  This bounded check establishes
    that those numerical controls complete cleanly at the parameter extremes
    before the same branch is used in the 972-case production scan.
    """
    if abs(t_end - 100.0) > 1.0e-9:
        raise ValueError("60 Td 高刚性验证固定使用 t_obs = 100 s")
    if report_dt > 10.0:
        raise ValueError("为计算 80--100 s 的末窗生成速率，report_dt 不得大于 10 s")
    conditions = [(tg, 60.0, n2) for tg in (300.0, 500.0) for n2 in (0.1, 0.9)]

    def execute(case: tuple[float, float, float]) -> dict[str, object]:
        tg, en, n2 = case
        return run_case(work, root, tg, en, n2, 1.0, t_end, report_dt, "hong_60td_stiffness_validation_100s",
                        atol=1.0e8, rtol=1.0e-4, internal_hmax_s=5.0e-3,
                        early_diagnostics=False, internal_mxstep=100000, timeout_s=3600.0)

    records: list[dict[str, object]] = []
    if workers == 1:
        for index, case in enumerate(conditions, 1):
            tg, en, n2 = case
            print(f"[60Td-validation {index}/{len(conditions)}] Tg={tg:g} K, E/N={en:g} Td, xN2={n2:g}")
            records.append(execute(case))
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(execute, case): case for case in conditions}
            for index, future in enumerate(as_completed(futures), 1):
                tg, en, n2 = futures[future]
                print(f"[60Td-validation {index}/{len(conditions)}] Tg={tg:g} K, E/N={en:g} Td, xN2={n2:g}")
                records.append(future.result())
    write_csv(root / "hong_60td_stiffness_validation_runs.csv", records)
    summary: list[dict[str, object]] = []
    for base in records:
        if base.get("returncode") != 0 or base.get("timed_out"):
            continue
        tag = str(base["tag"])
        series = root / str(base["group"]) / tag / f"cw_surface_{tag}.csv"
        summary.append(exposure_record(read_series(series), base, t_end))
    write_csv(root / "hong_60td_stiffness_validation_summary.csv", summary)
    (root / "hong_60td_stiffness_validation_protocol.json").write_text(json.dumps({
        "purpose": "Validate the 60 Td high-stiffness numerical branch before the production grid.",
        "surface_mechanism": "Hong 2017/2018 corrected, Table-5 metal column",
        "reactor_boundary": "closed 0D; no CSTR/feed/outlet/residence-time term",
        "forcing": "continuous E/N; no pulse waveform or cycle count",
        "conditions": [{"Tgas_K": tg, "EoverN_Td": en, "x_N2": n2, "x_H2": 1.0 - n2}
                       for tg, en, n2 in conditions],
        "finite_exposure_time_s": t_end,
        "report_interval_s": report_dt,
        "solver": {"ATOL_cm-3": 1.0e8, "RTOL": 1.0e-4,
                   "DVODE_HMAX_s": 5.0e-3, "MXSTEP": 100000},
        "acceptance": {"return_code": 0, "not_timed_out": True,
                       "site_balance_error": SITE_BALANCE_LIMIT},
    }, ensure_ascii=False, indent=2), encoding="utf-8")


def parse_positive_list(value: str) -> list[float]:
    values = sorted({float(item.strip()) for item in value.split(",") if item.strip()})
    if not values or any(item <= 0.0 for item in values):
        raise argparse.ArgumentTypeError("数值应为逗号分隔的正数")
    return values


def main() -> int:
    parser = argparse.ArgumentParser(description="CW gas--surface DBD convergence and parameter scans")
    parser.add_argument("mode", choices=("prepare", "hong-exposure-study", "hong-60td-validation", "hong-finite-grid", "time-study", "low-field-study", "numerical-study", "grid"))
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--t-end-s", type=float, default=T_END_S)
    parser.add_argument("--report-dt-s", type=float, default=REPORT_DT_S)
    parser.add_argument("--windows", type=parse_positive_list, default=[0.05, 0.2, 1.0, 5.0])
    parser.add_argument("--hmax-values", type=parse_positive_list, default=[0.01, 0.1, 1.0])
    parser.add_argument("--rtol-values", type=parse_positive_list, default=[1.0e-5, 1.0e-4, 1.0e-3])
    parser.add_argument("--grid-window-dt-s", type=float)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.t_end_s <= 0.0 or args.report_dt_s <= 0.0 or args.workers <= 0:
        parser.error("t-end-s、report-dt-s 与 workers 必须为正数")
    root = args.output_root.resolve()
    work = prepare_worktree(root)
    if args.mode == "hong-exposure-study":
        hong_exposure_study(work, root, args.t_end_s, args.report_dt_s, args.workers)
    elif args.mode == "hong-60td-validation":
        hong_60td_stiffness_validation(work, root, args.t_end_s, args.report_dt_s, args.workers)
    elif args.mode == "hong-finite-grid":
        hong_finite_exposure_grid(work, root, args.t_end_s, args.report_dt_s, args.workers)
    elif args.mode == "time-study":
        time_study(work, root, args.windows, args.t_end_s, args.report_dt_s, args.workers)
    elif args.mode == "low-field-study":
        low_field_study(work, root, args.windows, args.t_end_s, args.report_dt_s, args.workers)
    elif args.mode == "numerical-study":
        numerical_convergence_study(work, root, args.hmax_values, args.rtol_values,
                                    args.t_end_s, args.report_dt_s, args.workers)
    elif args.mode == "grid":
        if args.grid_window_dt_s is None or args.grid_window_dt_s <= 0.0:
            parser.error("grid 模式必须设置 --grid-window-dt-s")
        full_grid(work, root, args.grid_window_dt_s, args.t_end_s, args.report_dt_s, args.workers)
    print(json.dumps({"output_root": str(root), "mode": args.mode}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

