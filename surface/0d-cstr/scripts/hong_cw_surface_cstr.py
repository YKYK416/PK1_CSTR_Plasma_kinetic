#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Hong corrected gas--surface CW-CSTR physical driver.

This module is deliberately separate from the closed-0D driver.  It reuses
the isolated corrected-Hong worktree preparation, then adds a gas-heavy-
species CSTR source term.  Electrons, algebraic third bodies, and all SURF
states are excluded from the flow mask in the generated Fortran module.
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Callable

from _paths import LEGACY_HONG_ROOT
import hong_cw_surface_longtime_scan as closed


PROJECT = closed.PROJECT
BASE = closed.BASE
KINETIC_INPUT = closed.KINETIC_INPUT
NE_CM3 = closed.NE_CM3
SURFACE_SITE_DENSITY_CM3 = closed.SURFACE_SITE_DENSITY_CM3
SITE_BALANCE_LIMIT = closed.SITE_BALANCE_LIMIT
read_series = closed.read_series
solver_completion_evidence = closed.solver_completion_evidence
TAU_RES_S = 1.0e-2
MODEL_LABEL = "Hong 2017/2018 corrected gas--surface CW-CSTR, tau=10 ms"


def cstr_driver() -> str:
    """Derive the CSTR driver from the tested closed-0D surface driver."""
    source = closed.CW_SURFACE_DRIVER
    replacements = (
        ("! Closed-0D CW gas--surface driver.", "! CW gas--surface CSTR driver; gas residence time is supplied as CLI arg 14."),
        ("double precision :: ntot, n2_frac, Tgas, EN, ne0", "double precision :: ntot, n2_frac, Tgas, EN, ne0, tau_res, feed_density(species_max)"),
        ("site_projection_flag = 0\n  site_projection_max_rel", "site_projection_flag = 0\n  tau_res = 1.0d-2\n  site_projection_max_rel"),
        ("if (iargc() .ge. 13) then; call getarg(13,arg); read(arg,*) site_projection_flag; end if", "if (iargc() .ge. 13) then; call getarg(13,arg); read(arg,*) site_projection_flag; end if\n  if (iargc() .ge. 14) then; call getarg(14,arg); read(arg,*) tau_res; end if"),
        ("if (n2_frac .le. 0.0d0 .or. n2_frac .ge. 1.0d0) stop 'n2_frac must be in (0,1)'", "if (n2_frac .le. 0.0d0 .or. n2_frac .ge. 1.0d0) stop 'n2_frac must be in (0,1)'\n  if (tau_res .le. 0.0d0) stop 'tau_res must be positive'"),
        ("call ZDPlasKin_set_density('E', ne0, ldens_const=.true.)", "call ZDPlasKin_set_density('E', ne0, ldens_const=.true.)\n  ! Feed contains N2/H2 only; the generated flow mask excludes E and all SURF states.\n  feed_density = density\n  call ZDPlasKin_set_cstr_flow(tau_res, feed_density)"),
    )
    for old, new in replacements:
        if source.count(old) != 1:
            raise RuntimeError(f"CSTR driver transformation anchor is not unique: {old[:50]}")
        source = source.replace(old, new)
    source = source.replace("program main_cw_surface", "program main_cw_surface_cstr")
    source = source.replace("end program main_cw_surface", "end program main_cw_surface_cstr")
    return source


def instrument_module_with_cstr_flow(module_path: Path) -> None:
    """Add an explicit, gas-only CSTR source to an isolated generated module."""
    source = module_path.read_text(encoding="utf-8")
    marker = "! ZDP_CSTR_SURFACE_FLOW_PATCH_V2"
    if marker in source:
        return
    # The local runtime builder may already have injected its validated V1
    # source patch while compiling the corrected mechanism.  Reuse it only
    # after verifying that it has the essential immobile-surface exclusion.
    existing_marker = "! ZDP_CSTR_FLOW_PATCH_V1"
    if existing_marker in source:
        required = (
            "subroutine ZDPlasKin_set_cstr_flow(TAU_RES,FEED_DENSITY)",
            "if(index(trim(sname),'SURF') .gt. 0) cycle",
            "(cstr_flow_feed(:)-y(1:species_max))/cstr_flow_tau",
        )
        if not all(item in source for item in required):
            raise RuntimeError("existing CSTR flow patch does not meet the gas-only surface exclusion contract")
        return
    declaration_anchor = "! qtplaskin config"
    routine_anchor = "end subroutine ZDPlasKin_set_density"
    rhs_anchor = "  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0"
    if any(source.count(anchor) != 1 for anchor in (declaration_anchor, routine_anchor, rhs_anchor)):
        raise RuntimeError("CSTR flow patch anchors are missing or non-unique")
    declaration = f'''{marker}
! Gas-heavy-species CSTR boundary.  SURF states are immobile reactor states.
  logical, private :: cstr_flow_enabled = .false., cstr_flow_mask(species_max) = .false.
  double precision, private :: cstr_flow_tau = 0.0d0, cstr_flow_feed(species_max) = 0.0d0
!
'''
    routine = '''
! ZDP_CSTR_SURFACE_FLOW_PATCH_V2
subroutine ZDPlasKin_set_cstr_flow(TAU_RES,FEED_DENSITY)
  implicit none
  double precision, intent(in) :: TAU_RES, FEED_DENSITY(species_max)
  character(species_length) :: sname
  integer :: i
  if (TAU_RES .le. 0.0d0) stop 'CSTR residence time must be positive'
  cstr_flow_tau = TAU_RES
  cstr_flow_feed(:) = FEED_DENSITY(:)
  cstr_flow_enabled = .true.
  cstr_flow_mask(:) = .false.
  do i = 1, species_max
    sname = trim(adjustl(species_name(i)))
    if (i .eq. species_electrons) cycle
    if (trim(sname) .eq. 'M' .or. trim(sname) .eq. 'S') cycle
    if (index(trim(sname),'SURF') .gt. 0) cycle
    cstr_flow_mask(i) = .true.
  enddo
  vode_istate = 1
  return
end subroutine ZDPlasKin_set_cstr_flow
'''
    rhs = '''  ! ZDP_CSTR_SURFACE_FLOW_PATCH_V2: gas species only; E, M/S and SURF are masked.
  if(cstr_flow_enabled) then
    where(cstr_flow_mask(:))
      ydot(1:species_max) = ydot(1:species_max) + (cstr_flow_feed(:)-y(1:species_max))/cstr_flow_tau
    end where
  endif
'''
    source = source.replace(declaration_anchor, declaration + declaration_anchor)
    source = source.replace(routine_anchor, routine_anchor + routine)
    source = source.replace(rhs_anchor, rhs + rhs_anchor)
    module_path.write_text(source, encoding="utf-8", newline="\n")


def prepare_worktree(root: Path) -> Path:
    """Build a dedicated CSTR executable without modifying closed-0D runtime files."""
    work = closed.prepare_worktree(root)
    module_path = work / "zdplaskin_m.F90"
    driver_path = work / "main_cw_surface.F90"
    expected_driver = cstr_driver()
    if not driver_path.is_file() or driver_path.read_text(encoding="utf-8") != expected_driver:
        driver_path.write_text(expected_driver, encoding="utf-8", newline="\n")
        # ``runtime.build(..., mode='main')`` owns the validated V1 flow patch
        # injection.  Do not pre-inject a second module patch here.
        rc = closed.runtime.build(work, mode="main")
        if rc:
            raise RuntimeError(f"Hong gas--surface CSTR driver compilation failed, return code {rc}")
    module_source = module_path.read_text(encoding="utf-8")
    required_flow_contract = (
        "! ZDP_CSTR_FLOW_PATCH_V1",
        "subroutine ZDPlasKin_set_cstr_flow(TAU_RES,FEED_DENSITY)",
        "if(index(trim(sname),'SURF') .gt. 0) cycle",
        "(cstr_flow_feed(:)-y(1:species_max))/cstr_flow_tau",
    )
    if not all(item in module_source for item in required_flow_contract):
        raise RuntimeError("compiled CSTR worktree lacks the required gas-only V1 flow contract")
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.update({
        "purpose": "Hong corrected gas--surface CW-CSTR steady-state campaign",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "reactor_boundary": "0D CSTR; gas heavy species receive (C_in-C)/tau; electrons, M/S and SURF states are excluded",
        "residence_time_s": TAU_RES_S,
        "feed": "N2/H2 at the specified inlet composition; all reactive heavy species have zero inlet density",
        "driver": "main_cw_surface.F90",
        "worktree_input_patches": manifest.get("worktree_input_patches", []) + [
            "runtime CSTR flow patch V1: gas-heavy-species inflow/outflow only; SURF states excluded by generated uppercase species names"
        ],
    })
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return work


def run_case(work: Path, root: Path, tg: float, en: float, n2: float, dt: float,
             t_end: float, report_dt: float, group: str, *, atol: float = 1.0,
             rtol: float = 1.0e-4, internal_hmax_s: float | None = None,
             early_diagnostics: bool = True, internal_mxstep: int = 5000,
             site_projection: bool = False, residence_time_s: float = TAU_RES_S,
             timeout_s: float = 3600.0, case_output_dir: Path | None = None,
             progress_callback: Callable[[dict[str, object]], None] | None = None) -> dict[str, object]:
    """Run one CSTR condition using the shared, recoverable surface runner."""
    record = closed.run_case(
        work, root, tg, en, n2, dt, t_end, report_dt, group, atol=atol, rtol=rtol,
        internal_hmax_s=internal_hmax_s, early_diagnostics=early_diagnostics,
        internal_mxstep=internal_mxstep, site_projection=site_projection,
        residence_time_s=residence_time_s, timeout_s=timeout_s,
        case_output_dir=case_output_dir, progress_callback=progress_callback,
    )
    record.update({
        "reactor_boundary": "gas--surface 0D CSTR",
        "residence_time_s": residence_time_s,
        "feed": "N2/H2 only; reactive heavy species inlet density is zero",
        "surface_flow_policy": "all *SURF states excluded from CSTR flow source",
    })
    return record
