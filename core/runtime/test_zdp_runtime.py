#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python 原生运行时的无 Bash 单测。"""
import tempfile
from pathlib import Path

import zdp_runtime as r
import zdp_cli


ROOT = Path(__file__).resolve().parents[2] / "Hong"
BUILD = ROOT / "Reproduction" / "2017Hong" / "build"
TSCAN = ROOT / "Reproduction" / "2017Hong" / "方案1_振动开关温度" / "build_tscan"

assert r.resolve_gfortran_dir() and (r.resolve_gfortran_dir() / "gfortran.exe").is_file()
assert r.find_main_source(BUILD).name == "main_hong.F90"
assert r.find_main_source(TSCAN).name == "main_tscan.F90"
assert r.detect_flavor(BUILD)[0] == "hong"
assert r.detect_flavor(TSCAN)[0] == "tscan"

p = {"en": "60", "tg": "400", "n2frac": "0.5", "tend": "1e-6", "tag": "x",
     "en_off": "0.1", "freq": "5000", "duty": "0.2", "cycles": "2", "ne_mode": "fix"}
assert r._run_arguments("tscan", p) == (["0.5", "400", "60", "1e-6"], "")
assert r._run_arguments("pulse", p)[0][-1] == "fix"
assert r._run_arguments("auto", p) == (["400", "60", "1e-6"], "")
assert r._run_arguments("autopulse", p) == (["400", "60", "0.1", "5000", "0.2", "2"], "")
assert r._run_arguments("hong", p) == (["0.5", "1e-6"], "")
assert r._run_arguments("pulse", {**p, "tau_res": "1e-3"})[1] == "1e-3"
assert r._fortran_literal("45.1") == "45.1d0"
assert r._fortran_literal("1e-3") == "1d-3"
cli_args = zdp_cli._build_parser().parse_args(["run", "--n2-frac", "0.5"])
assert cli_args.n2frac == "0.5"
assert zdp_cli._build_parser().parse_args(["run", "--n2frac", "0.4"]).n2frac == "0.4"
assert zdp_cli._build_parser().parse_args(["run", "--tau-res", "1e-3"]).tau_res == "1e-3"

with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    try:
        r._safe_output_dir(td, "../outside", "x")
        raise AssertionError("输出路径越界应被拒绝")
    except RuntimeError:
        pass
    out = r._safe_output_dir(td, "runs/a", "x")
    assert out == (td / "runs" / "a").resolve()
    source = td / "main_hong.F90"
    source.write_text("program x\ncall z(GAS_TEMPERATURE=300.0d0, REDUCED_FIELD=45.1d0)\nend\n", encoding="utf-8")
    r._patch_hong_conditions(source, "60", "400", log_cb=lambda _x: None)
    changed = source.read_text(encoding="utf-8")
    assert "REDUCED_FIELD=60d0" in changed and "GAS_TEMPERATURE=400d0" in changed
    assert source.with_name("main_hong.F90.zdprunbak").is_file()
    r.restore_source(td, log_cb=lambda _x: None)
    assert "REDUCED_FIELD=45.1d0" in source.read_text(encoding="utf-8")

    # CSTR 补丁必须只在完整锚点存在时注入，且重复构建不重复修改模块。
    zd = td / "zdplaskin_m.F90"
    zd.write_text("""module ZDPlasKin
  integer, parameter :: species_max=3, species_length=16, species_electrons=1
  character(species_length) :: species_name(species_max)
  double precision :: y(species_max), ydot(species_max)
  logical :: ldensity_constant, density_constant(species_max)
  integer :: vode_istate
!
! qtplaskin config
!
subroutine ZDPlasKin_get_density(string,DENS,LDENS_CONST)
end subroutine ZDPlasKin_get_density
subroutine ZDPlasKin_fex()
  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0
end subroutine ZDPlasKin_fex
end module ZDPlasKin
""", encoding="utf-8")
    assert r._patch_cstr_flow(zd, True, log_cb=lambda _x: None)
    cstr = zd.read_text(encoding="utf-8")
    assert "subroutine ZDPlasKin_set_cstr_flow" in cstr
    assert "if(cstr_flow_enabled) then" in cstr and "end where" in cstr
    assert not r._patch_cstr_flow(zd, True, log_cb=lambda _x: None)

print("Python 原生运行时单测通过 ✓")
