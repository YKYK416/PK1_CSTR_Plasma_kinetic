#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ZDPlasKin 的 Windows 原生 Python 构建与运行运行时。

不依赖 Git Bash；直接使用 preprocessor.exe、gfortran.exe 和生成的
ZDPlasKin 主程序。GUI 与 zdp_cli 均通过本模块调用。
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
CONFIG_PATH = SCRIPT_DIR / "zdp_gui_config.json"
DEFAULT_GFORTRAN_DIR = Path(r"F:\Softwares\Ming64\mingw64\bin")
FFLAGS = ["-O2", "-ffree-line-length-none", "-static-libgfortran", "-static-libgcc"]


def load_config():
    try:
        import json
        return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
        return {}


def resolve_gfortran_dir(cfg=None):
    """定位 gfortran bin：config → GFORTRAN_DIR → PATH → 项目默认路径。"""
    cfg = cfg or load_config()
    candidates = []
    if cfg.get("gfortran_dir"):
        candidates.append(Path(cfg["gfortran_dir"]))
    if os.environ.get("GFORTRAN_DIR"):
        candidates.append(Path(os.environ["GFORTRAN_DIR"]))
    found = shutil.which("gfortran") or shutil.which("gfortran.exe")
    if found:
        candidates.append(Path(found).parent)
    candidates.append(DEFAULT_GFORTRAN_DIR)
    for directory in candidates:
        if (directory / "gfortran.exe").is_file():
            return directory
    return None


def _runtime_env(gfortran_dir=None):
    env = dict(os.environ)
    if gfortran_dir:
        env["PATH"] = str(gfortran_dir) + os.pathsep + env.get("PATH", "")
    # MinGW 在 TEMP/TMP 都不存在时可能错误回退到不可写系统目录。
    if not env.get("TEMP") and not env.get("TMP"):
        temp = tempfile.gettempdir()
        env["TEMP"] = temp
        env["TMP"] = temp
    return env


def _emit(log_cb, text):
    if log_cb:
        log_cb(str(text))


def run_process(args, cwd, log_cb=print, stop_holder=None, env=None, stdin_text=None,
                log_path=None):
    """运行外部程序并实时转发 stdout；停止时返回 -9。"""
    args = [str(x) for x in args]
    _emit(log_cb, "$ " + " ".join(args))
    try:
        proc = subprocess.Popen(args, cwd=str(cwd), stdin=subprocess.PIPE if stdin_text is not None else None,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                                encoding="utf-8", errors="replace", env=env)
    except OSError as exc:
        _emit(log_cb, f"[错误] 无法启动外部程序: {exc}")
        return 127
    if stop_holder is not None:
        lock = stop_holder.setdefault("lock", __import__("threading").Lock())
        with lock:
            stop_holder.setdefault("procs", set()).add(proc)
            stop_holder["proc"] = proc
    if stdin_text is not None and proc.stdin is not None:
        try:
            proc.stdin.write(stdin_text)
            proc.stdin.close()
        except OSError:
            pass
    handle = Path(log_path).open("w", encoding="utf-8", newline="") if log_path else None
    try:
        assert proc.stdout is not None
        for line in proc.stdout:
            line = line.rstrip("\r\n")
            if handle:
                handle.write(line + "\n")
                handle.flush()
            _emit(log_cb, line)
            if stop_holder and stop_holder.get("killed"):
                break
        rc = proc.wait()
    finally:
        if handle:
            handle.close()
        if stop_holder is not None:
            with stop_holder.setdefault("lock", __import__("threading").Lock()):
                stop_holder.setdefault("procs", set()).discard(proc)
                if stop_holder.get("proc") is proc:
                    stop_holder["proc"] = next(iter(stop_holder["procs"]), None)
    return -9 if stop_holder and stop_holder.get("killed") else rc


def find_main_source(directory):
    d = Path(directory)
    mains = sorted([p for p in d.glob("main*.[fF]90")
                    if re.search(r"(?im)^\s*program\s", p.read_text(encoding="utf-8", errors="ignore"))])
    if len(mains) == 1:
        return mains[0]
    if len(mains) > 1:
        raise RuntimeError("发现多个 main* 主程序：" + ", ".join(p.name for p in mains))
    candidates = []
    for p in sorted(d.glob("*.[fF]90")):
        if p.stem.lower() in ("dvode_f90_m", "zdplaskin_m"):
            continue
        if re.search(r"(?im)^\s*program\s", p.read_text(encoding="utf-8", errors="ignore")):
            candidates.append(p)
    if len(candidates) != 1:
        detail = "未找到" if not candidates else "发现多个：" + ", ".join(p.name for p in candidates)
        raise RuntimeError(detail + "含 PROGRAM 的主程序 .F90 文件")
    return candidates[0]


def find_executable(directory):
    d = Path(directory)
    exes = sorted(d.glob("main_*.exe"))
    if exes:
        return exes[0]
    excluded = {"preprocessor.exe", "bolsigplus.exe", "bolsigminus.exe", "t.exe"}
    exes = [p for p in sorted(d.glob("*.exe")) if p.name.lower() not in excluded]
    if not exes:
        raise RuntimeError("未找到主程序 exe（请先执行 Python 原生构建）")
    return exes[0]


def detect_flavor(directory, exe=None):
    d = Path(directory)
    src = None
    if exe:
        for suffix in (".F90", ".f90"):
            candidate = d / (Path(exe).stem + suffix)
            if candidate.is_file():
                src = candidate
                break
    try:
        src = src or find_main_source(d)
    except RuntimeError:
        src = None
    if src:
        text = src.read_text(encoding="utf-8", errors="ignore")
        if "MAIN_AUTO_PULSE" in text:
            return "autopulse", src
        if "MAIN_AUTO" in text:
            return "auto", src
        if "EN_on" in text:
            return "pulse", src
        if re.search(r"EN_set|Tgas_K EN_Td", text):
            return "tscan", src
        return "hong", src
    name = Path(exe).stem.lower() if exe else ""
    if "auto_pulse" in name or "autopulse" in name:
        return "autopulse", None
    if "pulse" in name:
        return "pulse", None
    if "tscan" in name:
        return "tscan", None
    if "auto" in name:
        return "auto", None
    return "hong", None


def _find_case_file(directory, names):
    d = Path(directory)
    for name in names:
        p = d / name
        if p.is_file():
            return p
    return None


def _link_library_arg(directory):
    d = Path(directory)
    preferred = d / "bolsig_x86_64_g.lib"
    if preferred.is_file():
        return "-lbolsig_x86_64_g"
    libs = sorted(d.glob("bolsig*.lib"))
    if libs:
        return "-l" + libs[0].stem
    archives = sorted(d.glob("bolsig*.a"))
    if archives:
        return "-l:" + archives[0].name
    raise RuntimeError("未找到 bolsig 库（bolsig_x86_64_g.lib / bolsig*.lib / bolsig*.a）")


def _patch_density(zd_path, log_cb):
    text = zd_path.read_text(encoding="utf-8", errors="replace")
    patched = re.sub(r"(?m)^(\s*)lreaction_block, rrt\s*$", r"\1lreaction_block, rrt, density", text)
    if patched != text:
        zd_path.write_text(patched, encoding="utf-8")
        _emit(log_cb, "已应用标准补丁: reac_rates use 列表加入 density")


_CSTR_MARKER = "! ZDP_CSTR_FLOW_PATCH_V1"


def _patch_cstr_flow(zd_path, enabled, log_cb):
    """向预处理生成的模块注入可选 CSTR 源项。

    该补丁只在主程序显式请求 ``ZDPlasKin_set_cstr_flow`` 时应用；流动
    项进入 DVODE 的 RHS，避免在主程序外部改写密度而反复重启刚性求解器。
    """
    if not enabled:
        return False
    text = zd_path.read_text(encoding="utf-8", errors="replace")
    if _CSTR_MARKER in text:
        return False
    decl_anchor = "!\n! qtplaskin config\n!"
    setter_anchor = "subroutine ZDPlasKin_get_density(string,DENS,LDENS_CONST)"
    rhs_anchor = "  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0"
    if decl_anchor not in text or setter_anchor not in text or rhs_anchor not in text:
        raise RuntimeError("无法定位 ZDPlasKin CSTR 补丁锚点；拒绝对预处理模块做不确定修改")
    declaration = """! ZDP_CSTR_FLOW_PATCH_V1
! optional well-mixed reactor boundary; configured by ZDPlasKin_set_cstr_flow
  logical, private :: cstr_flow_enabled = .false., cstr_flow_mask(species_max) = .false.
  double precision, private :: cstr_flow_tau = 0.0d0, cstr_flow_feed(species_max) = 0.0d0
!
"""
    setter = """! ZDP_CSTR_FLOW_PATCH_V1
subroutine ZDPlasKin_set_cstr_flow(TAU_RES,FEED_DENSITY)
  implicit none
  double precision, intent(in) :: TAU_RES, FEED_DENSITY(species_max)
  character(species_length) :: sname
  integer :: i
  cstr_flow_tau = TAU_RES
  cstr_flow_feed(:) = FEED_DENSITY(:)
  cstr_flow_enabled = (TAU_RES .gt. 0.0d0)
  cstr_flow_mask(:) = .false.
  if(cstr_flow_enabled) then
    do i = 1, species_max
      sname = trim(adjustl(species_name(i)))
      if(i .eq. species_electrons) cycle
      if(trim(sname) .eq. 'M' .or. trim(sname) .eq. 'S') cycle
      if(index(trim(sname),'SURF') .gt. 0) cycle
      cstr_flow_mask(i) = .true.
    enddo
  endif
  vode_istate = 1
  return
end subroutine ZDPlasKin_set_cstr_flow
!-----------------------------------------------------------------------------------------------------------------------------------
!
"""
    rhs = """  ! ZDP_CSTR_FLOW_PATCH_V1: gas heavy-species CSTR source; electrons,
  ! algebraic third bodies, and surface states are excluded by cstr_flow_mask.
  if(cstr_flow_enabled) then
    where(cstr_flow_mask(:))
    ydot(1:species_max) = ydot(1:species_max) + (cstr_flow_feed(:)-y(1:species_max))/cstr_flow_tau
    end where
  endif
"""
    text = text.replace(decl_anchor, declaration + decl_anchor, 1)
    text = text.replace(setter_anchor, setter + setter_anchor, 1)
    text = text.replace(rhs_anchor, rhs + rhs_anchor, 1)
    zd_path.write_text(text, encoding="utf-8")
    _emit(log_cb, "已注入 CSTR 流动 RHS 补丁（DVODE 耦合积分）")
    return True


def _smoke_gfortran(gfortran, env):
    with tempfile.TemporaryDirectory(prefix="zdp_gfortran_") as td:
        src = Path(td) / "t.f90"
        src.write_text("program t\nprint *,1\nend program t\n", encoding="ascii")
        proc = subprocess.run([str(gfortran), "-c", str(src)], cwd=td, env=env,
                              capture_output=True, text=True, encoding="utf-8", errors="replace")
        if proc.returncode:
            raise RuntimeError("gfortran 试编译失败，请检查 MinGW 工具链是否完整：" + proc.stderr.strip())


def build(directory, mode="full", log_cb=print, stop_holder=None, cfg=None):
    """Python 原生构建。mode 为 full / main / clean；返回退出码，失败抛 RuntimeError。"""
    d = Path(directory).resolve()
    if not d.is_dir():
        raise RuntimeError(f"构建目录不存在: {d}")
    if mode == "clean":
        removed = 0
        for pattern in ("*.o", "*.mod"):
            for p in d.glob(pattern):
                p.unlink()
                removed += 1
        _emit(log_cb, f"已删除 {removed} 个 .o/.mod 中间产物")
        return 0
    gf_dir = resolve_gfortran_dir(cfg)
    if not gf_dir:
        raise RuntimeError("找不到 gfortran。请在环境自检/配置中指定 gfortran bin 目录")
    gfortran = gf_dir / "gfortran.exe"
    env = _runtime_env(gf_dir)
    _emit(log_cb, f"[信息] gfortran: {gfortran}")
    _smoke_gfortran(gfortran, env)
    dvode = _find_case_file(d, ("dvode_f90_m.F90", "dvode_f90_m.f90"))
    if not dvode:
        raise RuntimeError("未找到 dvode_f90_m.F90")
    main = find_main_source(d)
    exe = d / (main.stem + ".exe")
    link_arg = _link_library_arg(d)
    if mode == "full":
        kinet = _find_case_file(d, ("kinet.inp", "kinetics.inp"))
        preprocessor = d / "preprocessor.exe"
        if kinet and preprocessor.is_file():
            _emit(log_cb, "[信息] ===== 1/2：预处理 =====")
            rc = run_process([preprocessor, kinet.name], d, log_cb, stop_holder, env=env, stdin_text=".\n\n")
            if rc:
                return rc
        elif not _find_case_file(d, ("zdplaskin_m.F90", "zdplaskin_m.f90")):
            raise RuntimeError("既无 kinet.inp+preprocessor.exe，也无已有 zdplaskin_m.F90")
    zd = _find_case_file(d, ("zdplaskin_m.F90", "zdplaskin_m.f90"))
    if not zd:
        raise RuntimeError("未找到 zdplaskin_m.F90（预处理产物）")
    _patch_density(zd, log_cb)
    needs_cstr_flow = "ZDPlasKin_set_cstr_flow" in main.read_text(encoding="utf-8", errors="ignore")
    cstr_patched = _patch_cstr_flow(zd, needs_cstr_flow, log_cb)
    if mode == "main":
        if not (d / "zdplaskin_m.o").is_file() or not (d / "dvode_f90_m.o").is_file():
            raise RuntimeError("main 模式需要已有 zdplaskin_m.o / dvode_f90_m.o；请先全量构建")
        zd_source = zd.read_text(encoding="utf-8", errors="ignore")
        # A caller may apply a controlled interface extension (currently DVODE
        # HMAX) to the generated module after the bootstrap full build.  In
        # that case its .mod interface must be regenerated before compiling
        # the driver, just as for the optional CSTR extension.
        needs_module_rebuild = cstr_patched or _CSTR_MARKER in zd_source or "HMAX" in zd_source
        compile_files = [zd.name, main.name] if needs_module_rebuild else [main.name]
        _emit(log_cb, f"[信息] 重编译{' zdplaskin_m + ' if len(compile_files) > 1 else ' '}主程序: {main.name}")
    else:
        compile_files = [dvode.name, zd.name, main.name]
        _emit(log_cb, "[信息] ===== 2/2：编译 dvode + zdplaskin_m + 主程序 =====")
    rc = run_process([gfortran, *FFLAGS, "-c", *compile_files], d, log_cb, stop_holder, env=env)
    if rc:
        return rc
    rc = run_process([gfortran, "-O2", "-o", exe.name, main.with_suffix(".o").name,
                      "zdplaskin_m.o", "dvode_f90_m.o", "-L.", link_arg,
                      "-static-libgfortran", "-static-libgcc"], d, log_cb, stop_holder, env=env)
    if rc:
        return rc
    if not exe.is_file():
        raise RuntimeError(f"链接未生成预期产物: {exe.name}")
    _emit(log_cb, f"[成功] 构建完成: {exe.name}（{exe.stat().st_size} B）")
    return 0


def _fortran_literal(value):
    s = str(value).strip()
    if not re.match(r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?$", s):
        raise RuntimeError(f"不是合法 Fortran 数值: {value!r}")
    return re.sub(r"[eE]", "d", s) if re.search(r"[eEdD]", s) else s + "d0"


def restore_source(directory, log_cb=print):
    """恢复 Python 原生运行时为 Hong 硬编码修改创建的源码备份。"""
    d = Path(directory)
    backups = sorted(d.glob("*.zdprunbak"))
    if not backups:
        _emit(log_cb, "没有找到 .zdprunbak 备份文件")
        return 0
    for backup in backups:
        source = backup.with_suffix("")
        shutil.copy2(backup, source)
        _emit(log_cb, f"已恢复: {backup.name} → {source.name}（请重新构建使其生效）")
    return 0


def _patch_hong_conditions(source, en, tg, log_cb):
    backup = source.with_name(source.name + ".zdprunbak")
    if not backup.is_file():
        shutil.copy2(source, backup)
        _emit(log_cb, f"已备份源码: {backup.name}")
    text = source.read_text(encoding="utf-8", errors="replace")
    changes = 0
    for key, value in (("REDUCED_FIELD", en), ("GAS_TEMPERATURE", tg)):
        pattern = rf"{key}=([0-9][0-9.eEdD+-]*)"
        text, n = re.subn(pattern, f"{key}={_fortran_literal(value)}", text)
        changes += n
    if not changes:
        raise RuntimeError("未匹配到 REDUCED_FIELD=/GAS_TEMPERATURE= 数值赋值，拒绝静默修改源码")
    source.write_text(text, encoding="utf-8")
    _emit(log_cb, f"已补丁 Hong 硬编码条件：E/N={en} Td，Tgas={tg} K")


def _run_arguments(flavor, params):
    """按五种主程序的位置参数签名构造运行参数。"""
    en = str(params.get("en", "45.1"))
    tg = str(params.get("tg", "300"))
    n2frac = str(params.get("n2frac", "0.3333"))
    tend = str(params.get("tend", "1"))
    en_off = str(params.get("en_off", "0.1"))
    freq = str(params.get("freq", "1000"))
    duty = str(params.get("duty", "0.5"))
    cycles = str(params.get("cycles", "3"))
    ne_mode = str(params.get("ne_mode", "fix"))
    tau_res = str(params.get("tau_res", "")).strip()
    if params.get("cw") and flavor in ("pulse", "autopulse"):
        duty = "1.0"
    if flavor == "tscan":
        args = [n2frac, tg, en, tend]
    elif flavor == "pulse":
        args = [n2frac, tg, en, en_off, freq, duty, cycles, ne_mode]
    elif flavor == "autopulse":
        args = [tg, en, en_off, freq, duty, cycles]
    elif flavor == "auto":
        args = [tg, en, tend]
    else:
        args = [n2frac, tend]
    return args, tau_res


def _safe_output_dir(directory, out_rel, tag):
    base = Path(directory).resolve()
    candidate = base / (out_rel or str(Path("runs") / tag))
    resolved = candidate.resolve()
    try:
        resolved.relative_to(base)
    except ValueError as exc:
        raise RuntimeError("输出目录必须位于算例目录内") from exc
    return resolved


def run(directory, params, log_cb=print, stop_holder=None, cfg=None):
    """Python 原生运行一个算例；params 使用 GUI 的同名键。"""
    d = Path(directory).resolve()
    if not d.is_dir():
        raise RuntimeError(f"运行目录不存在: {d}")
    if params.get("restore"):
        return restore_source(d, log_cb)
    exe = find_executable(d)
    flavor, source = detect_flavor(d, exe)
    en, tg = str(params.get("en", "45.1")), str(params.get("tg", "300"))
    if flavor == "hong" and (en != "45.1" or tg != "300"):
        if not params.get("recompile"):
            raise RuntimeError("hong 型主程序的 E/N 与温度为源码硬编码；请启用 recompile 才能修改")
        if not source:
            raise RuntimeError("Hong 参数重编译需要对应的 .F90 主程序源码")
        _patch_hong_conditions(source, en, tg, log_cb)
        rc = build(d, mode="main", log_cb=log_cb, stop_holder=stop_holder, cfg=cfg)
        if rc:
            return rc
        exe = find_executable(d)
    atol, rtol = str(params.get("atol", "")).strip(), str(params.get("rtol", "")).strip()
    if bool(atol) != bool(rtol):
        raise RuntimeError("atol 与 rtol 必须同时提供")
    if flavor == "hong" and atol:
        raise RuntimeError("hong 型主程序不读取 atol/rtol 位置参数")
    tag = str(params.get("tag") or datetime.now().strftime("r%m%d_%H%M%S"))
    outdir = _safe_output_dir(d, params.get("out_rel"), tag)
    outdir.mkdir(parents=True, exist_ok=True)
    bolsig = d / "bolsigdb.dat"
    if bolsig.is_file():
        shutil.copy2(bolsig, outdir / bolsig.name)
    for dat in d.glob("*.DAT"):
        shutil.copy2(dat, outdir / dat.name)
    gf_dir = resolve_gfortran_dir(cfg)
    env = _runtime_env(gf_dir)
    run_args, tau_res = _run_arguments(flavor, params)
    if tau_res:
        if flavor not in ("pulse", "autopulse"):
            raise RuntimeError("tau_res 仅适用于脉冲主程序")
        _fortran_literal(tau_res)
        if float(tau_res.replace("d", "e").replace("D", "E")) < 0:
            raise RuntimeError("tau_res 必须 ≥ 0")
    args = run_args + [tag]
    if atol or tau_res:
        args += [atol or "1.0", rtol or "1e-4"]
    if tau_res:
        args += [tau_res]
    _emit(log_cb, f"[信息] 运行 {exe.name}（{flavor} 型）→ {outdir}")
    log_path = None if params.get("no_log") else outdir / "console.log"
    rc = run_process([exe, *args], outdir, log_cb, stop_holder, env=env, log_path=log_path)
    if rc:
        _emit(log_cb, f"[失败] 运行退出码 {rc}；日志: {log_path or '未保存'}")
    else:
        _emit(log_cb, "[成功] 运行完成")
    return rc
