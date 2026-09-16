#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zdp_gui.py — ZDPlasKin 计算控制台（Tkinter GUI）

功能：选目录 → 自动检测输入文件 →
      Python 原生预处理 / 编译 / 运行 → 实时日志与进度 → 输出文件清单与快速绘图。

依赖：本目录的 zdp_cli.py 与 zdp_runtime.py；子进程直接调用
preprocessor.exe、gfortran.exe 与主程序，不依赖 Git Bash 或 shell 脚本。
"""
import os
import re
import csv
import sys
import argparse
import time
import math
import queue
import shutil
import threading
import subprocess
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PYTHON_CLI = SCRIPT_DIR / "zdp_cli.py"

FLAVOR_NAMES = {"hong": "恒定场·参数固化型（E/N、温度硬编码，改动需重编译）",
                "tscan": "恒定场·命令行型（全参数命令行直传）",
                "pulse": "脉冲型（方波脉冲，全参数命令行）",
                "auto": "向导生成型（组成内置，E/N/温度走命令行）",
                "autopulse": "脉冲·向导生成型（组成内置，脉冲参数走命令行）"}

# ============================================================
# 非界面逻辑（可独立测试）
# ============================================================

def _find_main_src(d: Path):
    """主程序探测：main* 且含 PROGRAM 优先。"""
    mains = sorted([p for p in d.glob("main*.[fF]90") if p.suffix.lower() == ".f90"])
    prog_mains = [p for p in mains
                  if re.search(r"(?im)^\s*program\s", p.read_text(encoding="utf-8", errors="ignore"))]
    if len(prog_mains) >= 1:
        return prog_mains[0]
    for p in sorted(d.glob("*.[fF]90")):
        if p.stem.lower() in ("dvode_f90_m", "zdplaskin_m"):
            continue
        if re.search(r"(?im)^\s*program\s", p.read_text(encoding="utf-8", errors="ignore")):
            return p
    return None


def detect_inputs(directory):
    """检测目录关键文件，返回 [(标签, 是否满足, 说明), ...]。"""
    d = Path(directory)
    checks = []
    inp = None
    for name in ("kinet.inp", "kinetics.inp"):
        if (d / name).is_file():
            inp = name
            break
    checks.append(("机理文件 kinet.inp", inp is not None, inp or "未找到"))
    checks.append(("截面库 bolsigdb.dat", (d / "bolsigdb.dat").is_file(), ""))
    main = _find_main_src(d)
    checks.append(("主程序 .f90/.F90", main is not None, main.name if main else "未找到含 PROGRAM 的源文件"))
    zd = any((d / n).is_file() for n in ("zdplaskin_m.F90", "zdplaskin_m.f90"))
    checks.append(("zdplaskin_m（预处理产物）", zd, "缺失时 Python 原生构建会先运行 preprocessor 生成" if not zd else ""))
    checks.append(("preprocessor.exe", (d / "preprocessor.exe").is_file(),
                   "无 kinet.inp 变更时非必需" if not (d / "preprocessor.exe").is_file() else ""))
    dv = any((d / n).is_file() for n in ("dvode_f90_m.F90", "dvode_f90_m.f90"))
    checks.append(("dvode_f90_m.F90", dv, ""))
    lib = any(d.glob("bolsig*.lib")) or any(d.glob("bolsig*.a"))
    dll = any(d.glob("bolsig*.dll"))
    checks.append(("bolsig 库 (.lib/.a + .dll)", lib and dll,
                   "优先 bolsig_x86_64_g 变体" if (d / "bolsig_x86_64_g.lib").is_file() else "未找到 x86_64_g 变体，原生构建将回退其他变体"))
    # 可编译的最低条件：机理或已生成模块 + 截面库 + 主程序 + dvode + bolsig 库
    can_build = ((inp is not None and (d / "preprocessor.exe").is_file()) or zd) \
        and (d / "bolsigdb.dat").is_file() and main is not None and dv and lib
    return checks, can_build


def detect_flavor(directory):
    """识别主程序类型 hong/tscan/pulse。"""
    d = Path(directory)
    src = _find_main_src(d)
    if src is None:
        return None
    text = src.read_text(encoding="utf-8", errors="ignore")
    # 顺序敏感：MAIN_AUTO_PULSE 含子串 MAIN_AUTO 且源码含 EN_on，必须最先判
    if "MAIN_AUTO_PULSE" in text:
        return "autopulse"
    if "MAIN_AUTO" in text:
        return "auto"
    if "EN_on" in text:
        return "pulse"
    if re.search(r"EN_set|Tgas_K EN_Td", text):
        return "tscan"
    return "hong"


_NUM = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?$")

def validate_params(p: dict) -> dict:
    """前端数值校验。不合法抛 ValueError（中文消息），合法返回清洗后的 dict。"""
    out = dict(p)
    def num(key, label, lo=None, hi=None, integer=False):
        v = str(out.get(key, "")).strip()
        if not _NUM.match(v):
            raise ValueError(f"{label} 不是合法数值: {v!r}")
        x = float(v.replace("d", "e").replace("D", "E"))
        if integer and abs(x - round(x)) > 1e-9:
            raise ValueError(f"{label} 必须是整数: {v!r}")
        if lo is not None and not (x > lo if label != "占空比" else x >= lo):
            raise ValueError(f"{label} 必须 > {lo}: {v!r}")
        if hi is not None and x > hi:
            raise ValueError(f"{label} 必须 ≤ {hi}: {v!r}")
        out[key] = str(int(round(x))) if integer else v
    num("en", "E/N (Td)", lo=0)
    num("tg", "气体温度 (K)", lo=0)
    num("n2frac", "N2 摩尔分数", lo=0, hi=1)
    num("tend", "模拟时长 (s)", lo=0)
    if out.get("en_off"):
        num("en_off", "EN_off (Td)", lo=0)
    if out.get("freq"):
        num("freq", "脉冲频率 (Hz)", lo=0)
    if out.get("duty"):
        num("duty", "占空比", lo=0, hi=1)
    if out.get("cycles"):
        num("cycles", "周期数", lo=0, integer=True)
    if out.get("tau_res"):
        num("tau_res", "停留时间 tau_res (s)", lo=-1)
        if _to_float(out["tau_res"]) < 0:
            raise ValueError("停留时间 tau_res (s) 必须 ≥ 0")
    for k, lab in (("atol", "atol"), ("rtol", "rtol")):
        v = str(out.get(k, "")).strip()
        if v:
            if not _NUM.match(v):
                raise ValueError(f"{lab} 不是合法数值: {v!r}")
            out[k] = v
    if bool(out.get("atol")) != bool(out.get("rtol")):
        raise ValueError("atol 与 rtol 必须同时填写（原生运行参数不能错位）")
    if out.get("retry_on_failure"):
        for k, lab in (("retry_atol", "二次 atol"), ("retry_rtol", "二次 rtol")):
            v = str(out.get(k, "")).strip()
            if not _NUM.match(v):
                raise ValueError(f"{lab} 不是合法数值: {v!r}")
            if _to_float(v) <= 0:
                raise ValueError(f"{lab} 必须 > 0: {v!r}")
            out[k] = v
    if out.get("scan_workers") not in (None, ""):
        num("scan_workers", "并行任务数", lo=0, integer=True)
    tag = str(out.get("tag", "")).strip()
    if not re.match(r"^[A-Za-z0-9_\-]+$", tag):
        raise ValueError(f"tag 只允许字母/数字/下划线/连字符: {tag!r}")
    out["tag"] = tag
    en_list = str(out.get("en_list", "")).strip()
    if en_list:
        items = re.split(r"[\s,;]+", en_list)
        for it in items:
            if not _NUM.match(it):
                raise ValueError(f"E/N 扫描列表含非法数值: {it!r}")
        out["en_list"] = " ".join(items)
    return out


def build_run_args(_unused: str, params: dict, flavor: str, out_rel: str = None) -> list:
    """按 Python 原生运行时的命名参数拼命令（列表形式，不拼 shell 字符串）。
    out_rel：自定义输出子目录（runs/<批次>/<参数组合> 或扁平 <参数组合>），空则默认 runs/<tag>。
    auto（向导生成型）组分固化在源码中：永不下发 --n2-frac（防止"改了以为生效"的静默失效）。"""
    args = [sys.executable, str(PYTHON_CLI), "run", "--directory", ".",
            "--en", str(params["en"]),
            "--tg", str(params["tg"])]
    if flavor not in ("auto", "autopulse"):
        args += ["--n2-frac", str(params["n2frac"])]
    args += ["--tend", str(params["tend"]),
             "--tag", str(params["tag"])]
    if flavor in ("pulse", "autopulse"):
        # CW is a pulse-program mode with duty forced to one, but its
        # frequency/cycle count still control the period-steady calculation.
        args += ["--en-off", str(params.get("en_off") or 0.1),
                 "--freq", str(params.get("freq") or 1000),
                 "--duty", str(params.get("duty") or 0.5),
                 "--cycles", str(params.get("cycles") or 3)]
        if flavor == "pulse":
            args += ["--ne-mode", str(params.get("ne_mode") or "fix")]
        if not params.get("pulse_enabled"):
            args += ["--cw"]  # 不勾选 = 连续恒定场（运行时强制 duty=1.0）
        if str(params.get("tau_res", "")).strip():
            args += ["--tau-res", str(params["tau_res"])]
    if flavor != "hong" and params.get("atol") and params.get("rtol"):
        args += ["--atol", str(params["atol"]), "--rtol", str(params["rtol"])]
    if params.get("en_list"):
        args += ["--en-list", str(params["en_list"])]
    if params.get("recompile"):
        args += ["--recompile"]
    if out_rel:
        args += ["--out-rel", str(out_rel)]
    if params.get("out_console") is False:
        args += ["--no-log"]  # 输出控制：不保存控制台日志
    return args


def build_build_args() -> list:
    return [sys.executable, str(PYTHON_CLI), "build", "--directory", "."]


_TIME_RE = re.compile(r"(?:^|\s)(?:time|t)\s*=\s*([0-9][0-9.eE+-]*)")

def parse_progress(line: str, t_end: float):
    """从输出行解析时间推进，返回 0-1 进度；解析不到返回 None。"""
    m = _TIME_RE.search(line)
    if not m:
        return None
    try:
        t = float(m.group(1))
    except ValueError:
        return None
    if t_end <= 0:
        return None
    return max(0.0, min(1.0, t / t_end))


def list_outputs(directory, tag: str):
    """列出 runs/<tag>/ 批次目录的全部产物（含 <参数组合>/ 子目录，兼容旧版平铺结构）。"""
    d = Path(directory) / "runs"
    results = []
    if not d.is_dir():
        return results
    for sub in sorted(d.iterdir()):
        if sub.is_dir() and (sub.name == tag or sub.name.startswith(tag + "_")):
            for f in sorted(sub.rglob("*")):
                if f.is_file():
                    results.append((f, f.stat().st_size))
    return results


# ============================================================
# 多维参数扫描（1D/2D/3D）：网格生成、指标注册表、逐点执行
# ============================================================

SCAN_AXES_BASE = [("en", "E/N (Td)"), ("tg", "气体温度 (K)"),
                  ("n2frac", "N2 摩尔分数"), ("tend", "模拟时长 (s)")]
SCAN_AXES_PULSE = [("freq", "脉冲频率 (Hz)"), ("duty", "占空比"), ("en_off", "EN_off (Td)"),
                   ("cycles", "最大周期数"), ("tau_res", "停留时间 (s；0=封闭)")]
SCAN_MAX_DIM = 3
SCAN_CONFIRM_POINTS = 100  # 超过此点数弹确认
_NAN = float("nan")
DEFAULT_ATOL = "1.0"
DEFAULT_RTOL = "1e-4"
DEFAULT_RETRY_ATOL = "10.0"
DEFAULT_RETRY_RTOL = "1e-3"


def default_scan_workers():
    """保留一个逻辑核心给桌面与 GUI，默认最多启动 4 个外部求解器。"""
    cpu = os.cpu_count() or 1
    return max(1, min(4, cpu - 1))


def scan_axes_for_flavor(flavor):
    """该主程序类型可选的扫描轴 [(key, 标签), ...]（面板动态过滤的唯一依据）。
    hong（恒定场·参数固化型）：E/N、Tgas 硬编码，仅 N2 摩尔分数 / 模拟时长可扫
    （hong 签名第一参数即 n2_frac、第二参数 t_end，命令行直传，不受 --en-list 限制影响）；
    auto（向导生成型）：组分固化在源码中，n2frac 不可扫；
    pulse：模拟时长由周期数决定（tend 不可扫），ne_mode 是枚举字符串、不做扫描轴。"""
    if flavor is None:
        return []
    if flavor == "hong":
        return [("n2frac", "N2 摩尔分数"), ("tend", "模拟时长 (s)")]
    if flavor == "auto":
        return [("en", "E/N (Td)"), ("tg", "气体温度 (K)"), ("tend", "模拟时长 (s)")]
    if flavor == "autopulse":
        # 组成内置（无 n2frac）；时长由 cycles/freq 决定（无 tend）；ne 恒定（无 ne_mode）
        return [("en", "EN_on (Td)"), ("tg", "气体温度 (K)")] + SCAN_AXES_PULSE
    axes = list(SCAN_AXES_BASE)
    if flavor == "pulse":
        axes = [(k, lab) for k, lab in axes if k != "tend"] + SCAN_AXES_PULSE
    return axes


def _to_float(s: str) -> float:
    return float(str(s).replace("d", "e").replace("D", "E"))


def parse_axis_values(spec: str):
    """解析一维扫描取值：显式列表（空格/逗号分隔，如 "30 60 90"）或
    区间写法 起:止:步数（线性均分，如 "30:120:5"）。返回 float 列表；非法抛 ValueError。"""
    s = str(spec).strip()
    if not s:
        raise ValueError("扫描轴取值不能为空")
    if ":" in s:
        parts = [x.strip() for x in s.split(":")]
        if len(parts) != 3:
            raise ValueError(f"区间写法应为 起:止:步数（如 30:120:5）: {s!r}")
        a_s, b_s, n_s = parts
        for v, lab in ((a_s, "区间起点"), (b_s, "区间终点")):
            if not _NUM.match(v):
                raise ValueError(f"{lab}不是合法数值: {v!r}")
        if not n_s.isdigit() or int(n_s) < 1:
            raise ValueError(f"区间步数必须是正整数: {n_s!r}")
        a, b, n = _to_float(a_s), _to_float(b_s), int(n_s)
        if n == 1:
            return [a]
        step = (b - a) / (n - 1)
        return [a + i * step for i in range(n)]
    vals = []
    for it in re.split(r"[\s,;]+", s):
        if not _NUM.match(it):
            raise ValueError(f"取值列表含非法数值: {it!r}")
        vals.append(_to_float(it))
    return vals


def fmt_axis(v: float) -> str:
    """扫描轴取值 → 命令行字符串（去尾零）。"""
    return f"{v:g}"


def build_scan_grid(axes):
    """axes: [(key, [v...]), ...]（1–3 维）。返回 [{key: v, ...}, ...] 笛卡尔积。"""
    if not 1 <= len(axes) <= SCAN_MAX_DIM:
        raise ValueError(f"扫描维度须为 1–{SCAN_MAX_DIM}")
    keys = [k for k, _ in axes]
    if len(set(keys)) != len(keys):
        raise ValueError("扫描轴参数不能重复")
    grid = [{}]
    for k, vals in axes:
        if not vals:
            raise ValueError(f"扫描轴 {k} 没有取值")
        grid = [{**g, k: v} for g in grid for v in vals]
    return grid


# ---------- 产物解析与指标注册表 ----------

def _read_csv_rows(path):
    """读 CSV → (表头 list, 数据 [[float,...], ...])。解析失败返回 (None, [])。"""
    try:
        with open(path, encoding="utf-8", errors="replace", newline="") as fh:
            rows = list(csv.reader(fh))
    except OSError:
        return None, []
    if len(rows) < 2:
        return None, []
    hdr = [h.strip() for h in rows[0]]
    data = []
    for r in rows[1:]:
        if len(r) < len(hdr):
            continue
        try:
            data.append([float(x) for x in r[:len(hdr)]])
        except ValueError:
            continue
    return hdr, data


def _col(hdr, data, name):
    if not hdr or name not in hdr:
        return None
    i = hdr.index(name)
    return [r[i] for r in data]


def _m_nh3_final(hdr, data):
    c = _col(hdr, data, "NH3")
    return c[-1] if c else _NAN


def _m_nh3_peak(hdr, data):
    c = _col(hdr, data, "NH3")
    return max(c) if c else _NAN


def _m_te_final(hdr, data):
    c = _col(hdr, data, "Te_eV")
    return c[-1] if c else _NAN


def _m_te_peak(hdr, data):
    c = _col(hdr, data, "Te_eV")
    return max(c) if c else _NAN


def _m_n2_conv(hdr, data):
    c = _col(hdr, data, "N2")
    if c and c[0] > 0:
        return (c[0] - c[-1]) / c[0] * 100.0
    return _NAN


# 时间序列指标注册表（三种主程序的序列 CSV 结构一致，同一组提取函数）
SERIES_METRICS = [
    ("NH3_final_cm-3", "末态 NH3 密度", _m_nh3_final),
    ("NH3_peak_cm-3", "峰值 NH3 密度", _m_nh3_peak),
    ("Te_final_eV", "末态电子温度", _m_te_final),
    ("Te_peak_eV", "峰值电子温度", _m_te_peak),
    ("N2_conv_pct", "N2 转化率 %", _m_n2_conv),
]

_SERIES_FILE = {"pulse": "pulse_series_{tag}.csv",
                "autopulse": "pulse_series_{tag}.csv",
                "tscan": "tscan_output_{tag}.csv",
                "hong": "hong_output_{tag}.csv",
                "auto": "auto_output_{tag}.csv"}


def _series_file(flavor, outdir: Path, tag: str):
    pat = _SERIES_FILE.get(flavor)
    if pat:
        f = outdir / pat.format(tag=tag)
        if f.is_file():
            return f
    cand = sorted(outdir.glob("*_output_*.csv")) or sorted(outdir.glob("pulse_series_*.csv"))
    return cand[0] if cand else None


def metric_names_for_flavor(flavor):
    """该类型汇总表的指标列（固定顺序）。"""
    names = [n for n, _d, _f in SERIES_METRICS]
    if flavor in ("pulse", "autopulse"):
        names += ["E_tot_Jcm3", "cycles_used", "NH3_per_J", "ne0_cm-3",
                  "rel_state_max", "rel_dNH3_cycle", "rel_energy_cycle", "stable_streak"]
    return names


def extract_metrics(flavor, outdir, tag: str) -> dict:
    """从 runs/<tag>/ 产物提取指标。解析不到的指标为 NaN，绝不抛异常。"""
    outdir = Path(outdir)
    out = {}
    hdr = data = None
    sf = _series_file(flavor, outdir, tag)
    if sf is not None:
        hdr, data = _read_csv_rows(sf)
    for name, _desc, fn in SERIES_METRICS:
        try:
            out[name] = fn(hdr, data) if data else _NAN
        except Exception:  # noqa: BLE001
            out[name] = _NAN
    if flavor in ("pulse", "autopulse"):
        chdr, cdata = (None, [])
        cf = outdir / f"pulse_cycles_{tag}.csv"
        if cf.is_file():
            chdr, cdata = _read_csv_rows(cf)
        e_on = _col(chdr, cdata, "E_on_Jcm3") or []
        e_off = _col(chdr, cdata, "E_off_Jcm3") or []
        e_tot = (sum(e_on) + sum(e_off)) if cdata else _NAN
        out["E_tot_Jcm3"] = e_tot
        out["cycles_used"] = float(len(cdata)) if cdata else _NAN
        nf = out.get("NH3_final_cm-3", _NAN)
        out["NH3_per_J"] = nf / e_tot if (e_tot == e_tot and e_tot > 0 and nf == nf) else _NAN
        for key in ("rel_state_max", "rel_dNH3_cycle", "rel_energy_cycle", "stable_streak"):
            values = _col(chdr, cdata, key) or []
            out[key] = values[-1] if values else _NAN
        try:
            console = (outdir / "console.log").read_text(encoding="utf-8", errors="replace")
        except OSError:
            console = ""
        out["ne0_cm-3"] = _console_float(_PULSE_NE0_RE.search(console))
    return out


def validate_scan_metrics(flavor, metrics):
    """校验扫描点的核心输出；返回 None 或可写入汇总表的失效原因。

    ZDPlasKin/DVODE 偶尔会以退出码 0 结束、却留下空文件或 NaN/Inf 序列。
    为避免把这类点当成成功结果，NH3 末态必须有限，且核心指标中至少一半
    非有限时判为无效。单个非关键指标缺失不会误伤正常算例。
    """
    core = ["NH3_final_cm-3", "NH3_peak_cm-3", "Te_final_eV", "N2_conv_pct"]
    if flavor in ("pulse", "autopulse"):
        core += ["E_tot_Jcm3", "cycles_used", "NH3_per_J"]
    bad = []
    for name in core:
        try:
            valid = math.isfinite(float(metrics.get(name, _NAN)))
        except (TypeError, ValueError):
            valid = False
        if not valid:
            bad.append(name)
    if "NH3_final_cm-3" in bad:
        return "NH3_final_cm-3 为 NaN/Inf 或主输出缺失"
    if len(bad) * 2 >= len(core):
        return f"核心指标 {len(bad)}/{len(core)} 项为 NaN/Inf：" + "、".join(bad)
    return None


_SERIES_NON_DENSITY = {"time_s", "phase", "te_ev", "en_td", "ptot", "pelast", "pinel"}
_PULSE_CONVERGED_RE = re.compile(r"cycles_used\s*=\s*(\d+)\s+converged\s*=\s*([01])", re.I)
_PULSE_TEND_RE = re.compile(r"\bt_end\s*=\s*([+-]?[\d.]+(?:[EeDd][+-]?\d+)?)", re.I)
_PULSE_NH3_RE = re.compile(r"NH3\s+final\s*\(cm-3\)\s*:\s*([+-]?[\d.]+(?:[EeDd][+-]?\d+)?)", re.I)
_PULSE_NE0_RE = re.compile(r"ne0_cm-3\s*=\s*([+-]?[\d.]+(?:[EeDd][+-]?\d+)?)", re.I)
_PULSE_STEADY_FIELDS = ("rel_state_max", "rel_dNH3_cycle", "rel_energy_cycle", "stable_streak")
_PULSE_STEADY_RTOL = 1.0e-3


def _console_float(match):
    """读取 Fortran 控制台科学计数法；解析失败时返回 NaN。"""
    try:
        return float(match.group(1).replace("D", "E").replace("d", "e"))
    except (AttributeError, TypeError, ValueError):
        return _NAN


def _number_param(params, name):
    """读取用户数值参数；解析失败返回 None（调用方转成明确的产物错误）。"""
    try:
        return float(str(params.get(name, "")).replace("d", "e").replace("D", "E"))
    except (AttributeError, TypeError, ValueError):
        return None


def _series_validity(hdr, data):
    """检查时间序列本身的基本物理约束，返回 (末时间, 错误原因)。

    密度中的极小负值可能来自数值舍入，因此按每一列最大绝对值的 1e-12
    给出容差；明显负密度、任意 NaN/Inf、时间倒退均视为不可用输出。
    """
    times = _col(hdr, data, "time_s")
    if not times or len(times) < 2:
        return None, "time_s 有效数据不足（至少需要 2 行）"
    if any(not math.isfinite(t) for t in times):
        return None, "time_s 含 NaN/Inf"
    if times[0] < 0 or any(b < a for a, b in zip(times, times[1:])):
        return None, "time_s 非单调或出现负时间"

    invalid, negative = [], []
    for i, name in enumerate(hdr or []):
        if name.strip().lower() in _SERIES_NON_DENSITY:
            continue
        values = [r[i] for r in data]
        if any(not math.isfinite(v) for v in values):
            invalid.append(name)
            continue
        scale = max((abs(v) for v in values), default=0.0)
        tol = max(1e-30, 1e-12 * scale)
        if any(v < -tol for v in values):
            negative.append(name)
    if invalid:
        return None, "序列含 NaN/Inf：" + "、".join(invalid[:4])
    if negative:
        return None, "物种密度显著为负：" + "、".join(negative[:4])
    return times[-1], None


def validate_run_output(flavor, outdir, tag, params, metrics=None):
    """将退出码 0 的产物分为成功、输出无效与周期未收敛。

    连续场程序只验证确实推进到 ``tend``；这不是、也不能代替稳态判定。
    脉冲程序则读取其原生 ``converged`` 标志，以区分周期稳态尚未达到的
    正常结束与数值损坏。返回 ``(status, failure_reason)``。
    """
    outdir = Path(outdir)
    metrics = metrics if metrics is not None else extract_metrics(flavor, outdir, tag)
    reason = validate_scan_metrics(flavor, metrics)
    if reason:
        return "invalid_output", reason

    sf = _series_file(flavor, outdir, tag)
    hdr, data = _read_csv_rows(sf) if sf is not None else (None, [])
    t_last, reason = _series_validity(hdr, data)
    if reason:
        return "invalid_output", reason

    if flavor not in ("pulse", "autopulse"):
        t_target = _number_param(params, "tend")
        if t_target is None or t_target <= 0:
            return "invalid_output", "无法核验目标终止时间 tend"
        tol = max(1e-15, abs(t_target) * 1e-6)
        if abs(t_last - t_target) > tol:
            return "invalid_output", f"序列末时间 {t_last:.6g} s 未达到目标 {t_target:.6g} s"
        return "success", ""

    cf = outdir / f"pulse_cycles_{tag}.csv"
    chdr, cdata = _read_csv_rows(cf) if cf.is_file() else (None, [])
    ctimes = _col(chdr, cdata, "t_end_s")
    if not ctimes or len(ctimes) < 1 or any(not math.isfinite(t) for t in ctimes):
        return "invalid_output", "pulse_cycles 缺少有效 t_end_s"
    if ctimes[0] <= 0 or any(b <= a for a, b in zip(ctimes, ctimes[1:])):
        return "invalid_output", "pulse_cycles 的 t_end_s 非严格递增"
    freq = _number_param(params, "freq")
    if freq is None or freq <= 0:
        return "invalid_output", "无法核验脉冲频率 freq"
    expected_end = len(ctimes) / freq
    tol = max(1e-15, abs(expected_end) * 1e-6)
    if abs(ctimes[-1] - expected_end) > tol or abs(t_last - ctimes[-1]) > tol:
        return "invalid_output", (f"脉冲时间不一致（series={t_last:.6g} s，"
                                  f"cycles={ctimes[-1]:.6g} s，期望={expected_end:.6g} s）")

    try:
        console = (outdir / "console.log").read_text(encoding="utf-8", errors="replace")
    except OSError:
        console = ""
    match = _PULSE_CONVERGED_RE.search(console)
    if not match:
        return "convergence_unknown", "未保留 console.log，无法确认周期稳态（请启用控制台日志）"
    console_cycles = int(match.group(1))
    if console_cycles != len(ctimes):
        return "invalid_output", (f"脉冲周期数不一致（console={console_cycles}，"
                                  f"cycles CSV={len(ctimes)}）")
    console_tend = _console_float(_PULSE_TEND_RE.search(console))
    if math.isfinite(console_tend) and abs(console_tend - ctimes[-1]) > tol:
        return "invalid_output", (f"脉冲终止时间不一致（console={console_tend:.6g} s，"
                                  f"cycles CSV={ctimes[-1]:.6g} s）")
    nh3 = _col(hdr, data, "NH3") or []
    console_nh3 = _console_float(_PULSE_NH3_RE.search(console))
    if nh3 and math.isfinite(console_nh3):
        # Fortran ES13.5 output keeps about six significant digits, so allow
        # its expected rounding error while still catching a whole extra cycle.
        nh3_tol = max(1e-30, abs(nh3[-1]) * 1e-5)
        if abs(console_nh3 - nh3[-1]) > nh3_tol:
            return "invalid_output", (f"NH3 终值不一致（console={console_nh3:.6g}，"
                                  f"series CSV={nh3[-1]:.6g} cm-3）")
    steady = {key: (_col(chdr, cdata, key) or []) for key in _PULSE_STEADY_FIELDS}
    if any(not values for values in steady.values()):
        return "convergence_unknown", "pulse_cycles 缺少 P0 周期稳态诊断列；请重编译最新脉冲主程序"
    last_steady = {key: values[-1] for key, values in steady.items()}
    if any(not math.isfinite(value) for value in last_steady.values()):
        return "invalid_output", "末周期稳态诊断含 NaN/Inf"
    if match.group(2) == "0":
        return "not_converged", f"脉冲程序在 {len(ctimes)} 个周期后未达到周期稳态"
    if last_steady["stable_streak"] < 3:
        return "invalid_output", "程序报告已收敛，但连续稳态周期数不足 3"
    if any(last_steady[key] > _PULSE_STEADY_RTOL * 1.01
           for key in ("rel_state_max", "rel_dNH3_cycle", "rel_energy_cycle")):
        return "invalid_output", "程序报告已收敛，但末周期稳态误差超过阈值"
    return "success", ""


def _scan_success(status):
    """扫描汇总中应视为可用数值结果的状态。"""
    return status in ("success", "retry_success")


def _fmt_cell(v):
    if isinstance(v, float):
        return "NaN" if not math.isfinite(v) else f"{v:.6g}"
    return str(v)


def write_scan_summary(path, axis_keys, metric_names, rows):
    """写 scan_summary.csv：各扫描轴、指标、状态、尝试次数、退出码和标签。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = list(axis_keys) + list(metric_names) + ["status", "failure_reason", "attempts", "elapsed_s", "rc", "tag"]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for row in rows:
            w.writerow([_fmt_cell(row.get(c, _NAN)) for c in cols])
    return path


def execute_scan(directory, base_params, axes, flavor, log_cb=print,
                 stop_holder=None, progress_cb=None, runner=None, flat_out=False,
                 max_workers=1, status_cb=None):
    """执行 1–3 维参数扫描，支持受控并行、失败跳过和一次宽松重试。

    每个点使用独立输出目录，故可由多个外部 ZDPlasKin 进程安全并行执行。
    ``progress_cb(done, total, row)`` 在每个已完成点调用；``status_cb`` 用于 GUI
    展示正在运行点。返回的 rows 始终按扫描网格顺序排列。
    """
    if flavor is None:
        raise ValueError("未识别主程序类型，无法扫描")
    axis_keys = [k for k, _ in axes]
    # 二次拒绝（面板剔除之外的硬守卫）
    if flavor == "hong":
        bad = sorted(set(axis_keys) - {"n2frac", "tend"})
        if bad:
            raise ValueError(f"恒定场·参数固化型（hong）主程序的 E/N、温度为硬编码，"
                             f"不能作为扫描轴: {bad}"
                             "（可扫 N2 摩尔分数 / 模拟时长；改 E/N、温度单跑可用 --recompile）")
    if flavor in ("auto", "autopulse") and "n2frac" in axis_keys:
        raise ValueError(f"向导生成型（{flavor}）主程序的气体组分已固化在源码中，"
                         "N2 摩尔分数不能作为扫描轴（修改组分请回第 4 步重新生成主程序）")
    if flavor in ("pulse", "autopulse") and not base_params.get("pulse_enabled"):
        pulse_only = {k for k, _ in SCAN_AXES_PULSE}
        bad = pulse_only & set(axis_keys)
        if bad:
            raise ValueError(f"未勾选脉冲模式（--cw 连续场）时，脉冲参数不能作为扫描轴: {sorted(bad)}")
    runner = runner or run_subprocess
    d = Path(directory)
    if not PYTHON_CLI.is_file():
        raise RuntimeError(f"Python 原生运行入口缺失: {PYTHON_CLI}")
    grid = build_scan_grid(axes)
    base_tag = str(base_params["tag"])
    metric_names = metric_names_for_flavor(flavor)
    rows_by_index, stopped = {}, False
    n = len(grid)
    try:
        max_workers = max(1, int(max_workers))
    except (TypeError, ValueError):
        max_workers = 1
    max_workers = min(max_workers, n)

    def run_point(i, point):
        """运行一个点；任何启动/运行异常均变为 failed 行，避免中断整个批次。"""
        if stop_holder and stop_holder.get("killed"):
            return None
        tag = combo_name(base_params, flavor, axis_keys, point)
        params = dict(base_params)
        for k in axis_keys:
            params[k] = fmt_axis(point[k])
        params["tag"] = tag
        # 批次目录用 base_tag（用户填的批次名），子目录用参数组合名（兼作 exe tag）
        out_rel = tag if flat_out else f"runs/{base_tag}/{tag}"
        args = build_run_args("", params, flavor, out_rel=out_rel)
        log_cb(f"[扫描] 第 {i}/{n} 点: " + ", ".join(f"{k}={fmt_axis(point[k])}" for k in axis_keys)
               + f"  →  {out_rel}")
        if status_cb:
            status_cb("started", i, n, tag, params)
        t0 = time.time()
        attempts, status = 1, "success"
        failure_reason = ""
        metrics = {}
        try:
            rc = runner(args, d, log_cb, stop_holder)
        except Exception as exc:  # noqa: BLE001
            log_cb(f"[失败] 第 {i}/{n} 点启动异常: {exc}")
            rc = 127
        if rc == -9:
            status = "stopped"
        elif rc != 0:
            status = "failed"
            failure_reason = f"求解器退出码 {rc}"
        else:
            metrics = extract_metrics(flavor, d / out_rel, tag)
            status, failure_reason = validate_run_output(flavor, d / out_rel, tag, params, metrics)
            if status == "invalid_output":
                log_cb(f"[失败] 第 {i}/{n} 点退出码 0 但数值输出无效：{failure_reason}")
            elif status == "not_converged":
                log_cb(f"[未收敛] 第 {i}/{n} 点已结束但未达到周期稳态：{failure_reason}")
            elif status == "convergence_unknown":
                log_cb(f"[警告] 第 {i}/{n} 点周期稳态无法核验：{failure_reason}")

        can_retry = (status in ("failed", "invalid_output")
                     and bool(base_params.get("retry_on_failure")) and flavor != "hong"
                     and not (stop_holder and stop_holder.get("killed")))
        if can_retry:
            outdir = d / out_rel
            first_log = outdir / "console.log"
            if first_log.is_file():
                shutil.copy2(first_log, outdir / "console_first_attempt.log")
            retry_params = dict(params)
            retry_params["atol"] = str(base_params["retry_atol"])
            retry_params["rtol"] = str(base_params["retry_rtol"])
            retry_params["retry_attempt"] = True
            retry_args = build_run_args("", retry_params, flavor, out_rel=out_rel)
            attempts = 2
            if status_cb:
                status_cb("retry", i, n, tag, retry_params)
            log_cb(f"[重试] 第 {i}/{n} 点首次{('数值无效' if status == 'invalid_output' else '失败')}"
                   f"（{failure_reason}），使用宽松容差 atol={retry_params['atol']}、"
                   f"rtol={retry_params['rtol']} 再试一次")
            try:
                rc = runner(retry_args, d, log_cb, stop_holder)
            except Exception as exc:  # noqa: BLE001
                log_cb(f"[失败] 第 {i}/{n} 点二次启动异常: {exc}")
                rc = 127
            args, params = retry_args, retry_params
            metrics = extract_metrics(flavor, d / out_rel, tag) if rc == 0 else {}
            if rc == 0:
                status, failure_reason = validate_run_output(flavor, d / out_rel, tag, params, metrics)
            else:
                status, failure_reason = "failed", f"二次求解器退出码 {rc}"
            if rc == 0 and status == "success":
                status = "retry_success"
            elif rc == -9:
                status, failure_reason = "stopped", "用户中止"
            elif rc == 0 and status == "invalid_output":
                status = "invalid_output"
                failure_reason = "二次收敛后仍数值无效：" + failure_reason
                log_cb(f"[失败] 第 {i}/{n} 点二次收敛后仍数值输出无效；标记并跳过：{failure_reason}")
            elif rc == 0:
                log_cb(f"[警告] 第 {i}/{n} 点二次运行结束但结果不可用于汇总：{failure_reason}")
            else:
                status = "failed"
                log_cb(f"[失败] 第 {i}/{n} 点二次收敛失败；标记并跳过：{failure_reason}")
        elapsed_s = time.time() - t0
        finalize_run_record(d, out_rel, flavor, params, args, t0, rc, log_cb,
                            scan_axes=axis_keys, point=point,
                            status=status, attempts=attempts, elapsed_s=elapsed_s,
                            failure_reason=failure_reason or None)
        if not _scan_success(status):
            metrics = {}
        row = dict(point)
        for name in metric_names:
            row[name] = metrics.get(name, _NAN)
        row["rc"] = rc
        row["tag"] = tag
        row["status"] = status
        row["failure_reason"] = failure_reason
        row["attempts"] = attempts
        row["elapsed_s"] = elapsed_s
        return i, row

    completed = 0
    if max_workers == 1:
        for i, point in enumerate(grid, 1):
            if stop_holder and stop_holder.get("killed"):
                stopped = True
                break
            result = run_point(i, point)
            if result is None:
                stopped = True
                break
            idx, row = result
            rows_by_index[idx] = row
            completed += 1
            if progress_cb:
                progress_cb(completed, n, row)
            if row["rc"] == -9:
                stopped = True
                break
    else:
        log_cb(f"[扫描] 启用并行：最多同时运行 {max_workers} 个独立任务")
        next_index = 1
        active = {}
        with ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix="zdp-scan") as pool:
            while next_index <= n and len(active) < max_workers:
                active[pool.submit(run_point, next_index, grid[next_index - 1])] = next_index
                next_index += 1
            while active:
                done, _ = wait(active, return_when=FIRST_COMPLETED)
                for future in done:
                    active.pop(future, None)
                    result = future.result()
                    if result is not None:
                        idx, row = result
                        rows_by_index[idx] = row
                        completed += 1
                        if progress_cb:
                            progress_cb(completed, n, row)
                        if row["rc"] == -9:
                            stopped = True
                if stop_holder and stop_holder.get("killed"):
                    stopped = True
                while not stopped and next_index <= n and len(active) < max_workers:
                    active[pool.submit(run_point, next_index, grid[next_index - 1])] = next_index
                    next_index += 1
    rows = [rows_by_index[i] for i in sorted(rows_by_index)]
    summary = d / "scan_summary.csv" if flat_out else d / "runs" / base_tag / "scan_summary.csv"
    if base_params.get("out_summary", True):
        summary = write_scan_summary(summary, axis_keys, metric_names, rows)
        n_failed = sum(1 for row in rows if row["status"] in ("failed", "invalid_output"))
        n_invalid = sum(1 for row in rows if row["status"] == "invalid_output")
        n_unconverged = sum(1 for row in rows if row["status"] == "not_converged")
        n_unknown = sum(1 for row in rows if row["status"] == "convergence_unknown")
        n_retried = sum(1 for row in rows if row["status"] == "retry_success")
        log_cb(f"[扫描] 完成 {len(rows)}/{n} 点（失败 {n_failed}，其中数值无效 {n_invalid}，"
               f"未周期收敛 {n_unconverged}，收敛未知 {n_unknown}，重试成功 {n_retried}），汇总: {summary}")
    else:
        log_cb(f"[扫描] 完成 {len(rows)}/{n} 点（输出控制：汇总 CSV 已关闭）")
    return rows, summary, stopped


# ============================================================
# 向导（wizard）：从零构建算例
# ============================================================

TEMPLATE_DIR = SCRIPT_DIR.parent / "templates"
TEMPLATE_MAIN = TEMPLATE_DIR / "main_auto.F90.tmpl"
TEMPLATE_MAIN_PULSE = TEMPLATE_DIR / "main_auto_pulse.F90.tmpl"
TEMPLATE_KINET = TEMPLATE_DIR / "kinet_min_N2H2.inp"
CONFIG_PATH = SCRIPT_DIR / "zdp_gui_config.json"

# 发行包最小程序文件集（必需 + 可选随行运行库）
DISTRO_REQUIRED = ["preprocessor.exe", "dvode_f90_m.F90",
                   "bolsig_x86_64_g.lib", "bolsig_x86_64_g.dll"]
DISTRO_OPTIONAL = ["libquadmath-0.dll", "libgcc_s_seh-1.dll",
                   "libgfortran-5.dll", "libwinpthread-1.dll"]
# 发行包里没有、但 gfortran 编译产物运行期可能需要的库：回退到 gfortran bin 目录取
RUNTIME_FALLBACK = ["libquadmath-0.dll"]


def find_gfortran_bin():
    """定位 gfortran 的 bin 目录。回退链：config['gfortran_dir'] → GFORTRAN_DIR 环境变量
    → PATH → 常见安装位置。返回 Path 或 None。"""
    r = resolve_gfortran_dir(load_config())
    return Path(r) if r else None


def deploy_program_files(distro_dir, case_dir):
    """把发行包最小程序文件集复制进算例目录（已存在则跳过）。
    返回 (copied, skipped, missing_required)。"""
    src, dst = Path(distro_dir), Path(case_dir)
    if not src.is_dir():
        raise FileNotFoundError(f"发行包目录不存在: {src}")
    dst.mkdir(parents=True, exist_ok=True)
    copied, skipped, missing = [], [], []
    for name in DISTRO_REQUIRED + DISTRO_OPTIONAL:
        s = src / name
        if not s.is_file() and name in RUNTIME_FALLBACK:
            gb = find_gfortran_bin()
            s = (gb / name) if gb else s
        if not s.is_file():
            if name in DISTRO_REQUIRED:
                missing.append(name)
            continue
        t = dst / name
        if t.is_file():
            skipped.append(name)
        else:
            shutil.copy2(s, t)
            copied.append(name)
    return copied, skipped, missing


_BOLSIG_KIND = re.compile(r"^(ELASTIC|EFFECTIVE|EXCITATION|IONIZATION|ATTACHMENT|TOTAL|MOMENTUM|ROTATION|VIBRATION)\s*$", re.I)


def merge_bolsigdb(txt_paths, out_path):
    """把若干 BOLSIG+ 格式（LXCat 下载）txt 顺序拼接成 bolsigdb.dat。
    返回 (out_path, 过程块数估计)。发行包没有官方转换工具，本项目做法是按块拼接。"""
    out_path = Path(out_path)
    n_blocks = 0
    with open(out_path, "w", encoding="utf-8") as out:
        for i, p in enumerate(txt_paths):
            p = Path(p)
            if not p.is_file():
                raise FileNotFoundError(f"截面文件不存在: {p}")
            text = p.read_text(encoding="utf-8", errors="replace")
            n_blocks += sum(1 for ln in text.splitlines() if _BOLSIG_KIND.match(ln.strip()))
            if i:
                out.write("\n")
            out.write(text.rstrip("\n") + "\n")
    return out_path, n_blocks


def parse_kinet_species(kinet_path):
    """解析 kinet.inp 的 SPECIES 段，返回物种名列表（保序去重）。解析失败返回 []。"""
    try:
        lines = Path(kinet_path).read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return []
    in_block, species = False, []
    for ln in lines:
        s = ln.strip()
        if not in_block:
            if re.match(r"(?i)^SPECIES\b", s):
                in_block = True
            continue
        if re.match(r"(?i)^END\b", s):
            break
        if not s or s.startswith(("#", "!")):
            continue
        for tok in s.split():
            if tok not in species:
                species.append(tok)
    return species


def parse_kinet_reaction_stats(kinet_path):
    """统计 ``REACTIONS`` 段的有效反应及其中 BOLSIG+ 计算的反应。

    ZDPlasKin 机理中的反应以 ``=>`` 表示；以 ``!`` 开始的行或以
    ``#`` 开始的注释不计入。BOLSIG+ 反应以同一反应行的 ``!`` 注释中
    出现 ``BOLSIG`` 为准，避免把固定速率的电子反应误算进去。
    解析失败时返回零计数和 ``parsed=False``，调用方应显示警告而不是猜测。
    """
    try:
        lines = Path(kinet_path).read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return {"total": 0, "bolsig": 0, "fixed": 0, "parsed": False}
    in_reactions = False
    total = bolsig = 0
    for raw in lines:
        stripped = raw.strip()
        if not in_reactions:
            if re.match(r"(?i)^REACTIONS\b", stripped):
                in_reactions = True
            continue
        if re.match(r"(?i)^END\b", stripped):
            break
        if not stripped or stripped.startswith(("!", "#", "$")):
            continue
        chemistry, sep, comment = raw.partition("!")
        if "=>" not in chemistry:
            continue
        total += 1
        if sep and re.search(r"\bBOLSIG(?:\+)?\b", comment, flags=re.I):
            bolsig += 1
    return {"total": total, "bolsig": bolsig, "fixed": total - bolsig,
            "parsed": in_reactions}


def validate_composition(rows):
    """校验初始气体组成 [(species, frac_str), ...]：
    物种非空且不重复、分数为 (0,1] 数值、总和 ≈1（容差 1e-3）。
    合法返回 [(species, float), ...]，非法抛 ValueError（中文）。"""
    if not rows:
        raise ValueError("初始气体组成不能为空")
    out = []
    seen = set()
    total = 0.0
    for sp, fs in rows:
        sp = str(sp).strip()
        fs = str(fs).strip()
        if not sp:
            raise ValueError("存在未选择物种的组分行")
        if sp in seen:
            raise ValueError(f"物种重复: {sp}")
        seen.add(sp)
        if not _NUM.match(fs):
            raise ValueError(f"摩尔分数不是合法数值（{sp}）: {fs!r}")
        x = _to_float(fs)
        if not (0.0 < x <= 1.0):
            raise ValueError(f"摩尔分数须在 (0,1] 区间（{sp}）: {fs!r}")
        total += x
        out.append((sp, x))
    if abs(total - 1.0) > 1e-3:
        raise ValueError(f"摩尔分数总和须等于 1（当前 {total:.6g}）")
    return out


def default_composition_for_species(species):
    """返回向导默认气体组成；N2/H2 同时存在时使用可编辑的 0.5/0.5。"""
    names = set(species or [])
    return [("N2", "0.5"), ("H2", "0.5")] if {"N2", "H2"} <= names else []


def _f90_num(x) -> str:
    """Python 数 → Fortran 双精度字面量（如 0.3333d0、1.17d8）。"""
    if isinstance(x, str):
        x = _to_float(x)
    s = f"{x:.6g}"
    if "e" in s or "E" in s:
        mant, ex = re.split(r"[eE]", s)
        return f"{mant}d{int(ex)}"
    return s + "d0"


DEFAULT_SITE_DENSITY = "1e16"  # 自由表面位点密度 [cm-3]，机理含表面物种 S 时写入，可生成后手改

# 输出控制（模块 E）→ 模板 @BLOCK 名
OUTPUT_BLOCK_MAP = {"species": "species_cols", "te": "te_col",
                    "rates_t": "rates_t", "rates": "rates"}


def _filter_blocks(src, disabled):
    """移除模板 @BLOCK 标记行；disabled 中的块连同内容一起删除。"""
    out, skip = [], False
    for ln in src.splitlines():
        s = ln.strip()
        m = re.match(r"!\s*@BLOCK:(\w+)", s)
        if m:
            skip = m.group(1) in disabled
            continue
        if re.match(r"!\s*@ENDBLOCK:\w+", s):
            skip = False
            continue
        if not skip:
            out.append(ln)
    return "\n".join(out) + "\n"


def render_main_auto(composition, en="45.1", tg="300", tend="1", ne="1.17e8",
                     species=None, template_path=None, output_opts=None):
    """用模板渲染通用主程序源码。composition: [(species, frac), ...]（已经 validate_composition）。
    species 为机理全物种列表（parse_kinet_species）时，自动补特殊物种初始密度：
    M（第三体）= ntot 恒定、S（自由表面位点）= DEFAULT_SITE_DENSITY——与 65 版参考主程序一致。
    output_opts（模块 E 输出控制）: {"species": bool, "te": bool, "rates_t": bool, "rates": bool}，
    为 False 的项对应 write 语句块被真实移除；None = 全开。
    返回源码字符串；保证无占位符与块标记残留。"""
    tpl = Path(template_path or TEMPLATE_MAIN).read_text(encoding="utf-8")
    lines = [f"  call ZDPlasKin_set_density('{sp}', {_f90_num(x)}*ntot)" for sp, x in composition]
    if species:
        comp_names = {sp for sp, _x in composition}
        if "M" in species and "M" not in comp_names:
            lines.append("  call ZDPlasKin_set_density('M', ntot, ldens_const=.true.)  ! 第三体（自动）")
        if "S" in species and "S" not in comp_names:
            lines.append(f"  call ZDPlasKin_set_density('S', {_f90_num(DEFAULT_SITE_DENSITY)})"
                         "  ! 自由表面位点（自动，可改）")
    init_lines = "\n".join(lines)
    src = (tpl.replace("{{INIT_LINES}}", init_lines)
              .replace("{{EN}}", _f90_num(en))
              .replace("{{TG}}", _f90_num(tg))
              .replace("{{TEND}}", _f90_num(tend))
              .replace("{{NE}}", _f90_num(ne)))
    disabled = set()
    if output_opts:
        disabled = {OUTPUT_BLOCK_MAP[k] for k, on in output_opts.items()
                    if k in OUTPUT_BLOCK_MAP and not on}
    src = _filter_blocks(src, disabled)
    if "{{" in src or "}}" in src:
        raise RuntimeError("模板占位符未全部替换")
    if "@BLOCK" in src or "@ENDBLOCK" in src:
        raise RuntimeError("模板块标记未全部处理")
    return src


def render_main_auto_pulse(composition, en_on="45.1", en_off="0.1", freq="1000",
                           duty="0.5", cycles="3", tg="300", ne="1.17e8",
                           species=None, template_path=None, output_opts=None):
    """渲染脉冲版通用主程序（main_auto_pulse）。CLI: Tgas EN_on EN_off freq duty cycles tag [atol rtol]。
    组成内置（无 n2 参数）、电子全程恒定（无 ne_mode）、模拟时长 = cycles/freq（无 tend）。
    composition/species/output_opts 语义同 render_main_auto；返回源码字符串，保证无占位符与块标记残留。"""
    tpl = Path(template_path or TEMPLATE_MAIN_PULSE).read_text(encoding="utf-8")
    lines = [f"  call ZDPlasKin_set_density('{sp}', {_f90_num(x)}*ntot)" for sp, x in composition]
    if species:
        comp_names = {sp for sp, _x in composition}
        if "M" in species and "M" not in comp_names:
            lines.append("  call ZDPlasKin_set_density('M', ntot, ldens_const=.true.)  ! 第三体（自动）")
        if "S" in species and "S" not in comp_names:
            lines.append(f"  call ZDPlasKin_set_density('S', {_f90_num(DEFAULT_SITE_DENSITY)})"
                         "  ! 自由表面位点（自动，可改）")
    init_lines = "\n".join(lines)
    src = (tpl.replace("{{INIT_LINES}}", init_lines)
              .replace("{{EN_ON}}", _f90_num(en_on))
              .replace("{{EN_OFF}}", _f90_num(en_off))
              .replace("{{FREQ}}", _f90_num(freq))
              .replace("{{DUTY}}", _f90_num(duty))
              .replace("{{CYCLES}}", str(int(_to_float(cycles))))
              .replace("{{TG}}", _f90_num(tg))
              .replace("{{NE}}", _f90_num(ne)))
    disabled = set()
    if output_opts:
        disabled = {OUTPUT_BLOCK_MAP[k] for k, on in output_opts.items()
                    if k in OUTPUT_BLOCK_MAP and not on}
    src = _filter_blocks(src, disabled)
    if "{{" in src or "}}" in src:
        raise RuntimeError("模板占位符未全部替换")
    if "@BLOCK" in src or "@ENDBLOCK" in src:
        raise RuntimeError("模板块标记未全部处理")
    return src


def eval_step_status(directory):
    """按目录内容评估向导前 4 步状态：{1: 程序文件齐, 2: 截面库, 3: 机理, 4: 主程序}。"""
    d = Path(directory)
    ok1 = d.is_dir() and (d / "preprocessor.exe").is_file() \
        and any((d / n).is_file() for n in ("dvode_f90_m.F90", "dvode_f90_m.f90")) \
        and any(d.glob("bolsig*.lib")) and any(d.glob("bolsig*.dll"))
    ok2 = d.is_dir() and (d / "bolsigdb.dat").is_file()
    inp = next(((d / n) for n in ("kinet.inp", "kinetics.inp") if (d / n).is_file()), None)
    ok3 = inp is not None and len(parse_kinet_species(inp)) >= 2
    ok4 = d.is_dir() and _find_main_src(d) is not None
    return {1: bool(ok1), 2: bool(ok2), 3: bool(ok3), 4: bool(ok4)}


def wizard_gate(status: dict) -> bool:
    """编译门禁：前 4 步全绿才允许。"""
    return all(status.get(s) for s in (1, 2, 3, 4))


def load_config():
    try:
        import json
        return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
        return {}


def save_config(cfg: dict):
    try:
        import json
        CONFIG_PATH.write_text(json.dumps(cfg, ensure_ascii=False, indent=2), encoding="utf-8")
    except Exception:  # noqa: BLE001
        pass


# ============================================================
# 输出目录命名（需求 1）：runs/<批次名>/<参数组合>/
# ============================================================

AXIS_PREFIX = {"en": "EN", "tg": "Tg", "n2frac": "n2", "tend": "t",
               "freq": "f", "duty": "d", "en_off": "ENoff", "cycles": "c"}


def combo_name(params, flavor, axis_keys=None, point=None):
    """参数组合子目录名。
    扫描点：按轴顺序拼接轴值（EN30_Tg400）；无扫描：EN<值>_Tg<值>_n2-<值>
    （pulse 追加 _f<频率>_d<占空比>，连续场追加 _cw）。"""
    if axis_keys and point is not None:
        return "_".join(f"{AXIS_PREFIX[k]}{fmt_axis(point[k])}" for k in axis_keys)
    parts = [f"EN{params['en']}", f"Tg{params['tg']}", f"n2-{params['n2frac']}"]
    if flavor in ("pulse", "autopulse"):
        if params.get("pulse_enabled"):
            parts.append(f"f{params.get('freq') or 1000}")
            parts.append(f"d{params.get('duty') or 0.5}")
        else:
            parts.append("cw")
    return "_".join(parts)


def out_rel_for(params, flavor, flat=False, axis_keys=None, point=None):
    """运行输出子目录（相对构建目录）。
    flat=False（原地目录）：runs/<批次名>/<参数组合>；
    flat=True（向导新建算例，构建目录本身即 <算例目录>/<批次名>）：<参数组合>。"""
    combo = combo_name(params, flavor, axis_keys, point)
    return combo if flat else f"runs/{params['tag']}/{combo}"


# ============================================================
# 参数档案（需求 1）：每个运行点子文件夹一份 计算参数.md + params.json
# ============================================================

def sha256_8(path):
    """文件 sha256 前 8 位（输入文件指纹）。"""
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()[:8]


def collect_input_files(directory):
    """算例输入清单 + sha256 前 8 位：kinet.inp、bolsigdb.dat、主程序源码、*.DAT。"""
    d = Path(directory)
    out = []
    if not d.is_dir():
        return out
    for name in ("kinet.inp", "kinetics.inp", "bolsigdb.dat"):
        p = d / name
        if p.is_file():
            out.append({"file": name, "sha256_8": sha256_8(p)})
    main = _find_main_src(d)
    if main is not None:
        out.append({"file": main.name, "sha256_8": sha256_8(main)})
    for p in sorted(d.glob("*.DAT")):
        out.append({"file": p.name, "sha256_8": sha256_8(p)})
    return out


def write_param_record(outdir, record):
    """把运行参数档案写进输出子文件夹：计算参数.md（人读）+ params.json（机读）。
    返回 (md_path, json_path)。"""
    import json
    from datetime import datetime
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rec = dict(record)
    rec.setdefault("written_at", datetime.now().isoformat(timespec="seconds"))
    js = outdir / "params.json"
    js.write_text(json.dumps(rec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    L = []
    L.append("# 计算参数档案")
    L.append("")
    L.append(f"- 主程序类型：{rec.get('flavor', '?')}")
    L.append(f"- 输出目录：{outdir}")
    L.append(f"- 开始时间：{rec.get('t_start', '?')}")
    L.append(f"- 结束时间：{rec.get('t_end', '?')}")
    rc = rec.get("rc")
    L.append(f"- 退出状态：{rc}（{'成功' if rc == 0 else ('中止' if rc == -9 else '失败')}）")
    L.append(f"- 批次状态：{rec.get('status', '?')}；尝试次数：{rec.get('attempts', 1)}")
    if rec.get("failure_reason"):
        L.append(f"- 失败原因：{rec['failure_reason']}")
    if rec.get("elapsed_s") is not None:
        L.append(f"- 本点耗时：{float(rec['elapsed_s']):.1f} s")
    if rec.get("scan_axes"):
        L.append(f"- 扫描轴：{', '.join(rec['scan_axes'])}")
    if rec.get("point"):
        L.append(f"- 本扫描点：{rec['point']}")
    L.append("")
    L.append("## 参数（含默认值）")
    L.append("")
    for k, v in (rec.get("params") or {}).items():
        L.append(f"- {k} = {v}")
    L.append("")
    L.append("## 实际命令行")
    L.append("")
    L.append("```")
    L.append(" ".join(str(a) for a in (rec.get("cmdline") or [])))
    L.append("```")
    L.append("")
    L.append("## 输入文件（sha256 前 8 位）")
    L.append("")
    for it in (rec.get("input_files") or []):
        L.append(f"- {it['file']}  `{it['sha256_8']}`")
    L.append("")
    L.append("## 产物清单")
    L.append("")
    for p in (rec.get("products") or []):
        L.append(f"- {p}")
    L.append("")
    md = outdir / "计算参数.md"
    md.write_text("\n".join(L), encoding="utf-8")
    return md, js


def finalize_run_record(directory, out_rel, flavor, params, cmdline, t0, rc, log_cb,
                        scan_axes=None, point=None, status=None, attempts=1, elapsed_s=None,
                        failure_reason=None):
    """运行（或扫描点）结束后写参数档案；任何异常只记日志，不影响主流程。"""
    try:
        from datetime import datetime
        outdir = Path(directory) / out_rel
        products = sorted(f.name for f in outdir.iterdir()
                          if f.is_file() and f.name not in ("计算参数.md", "params.json")) \
            if outdir.is_dir() else []
        record = {
            "flavor": flavor,
            "params": {k: v for k, v in params.items()},
            "scan_axes": list(scan_axes) if scan_axes else None,
            "point": ({k: fmt_axis(point[k]) for k in scan_axes}
                      if (scan_axes and point is not None) else None),
            "cmdline": [str(a) for a in cmdline],
            "input_files": collect_input_files(directory),
            "t_start": datetime.fromtimestamp(t0).isoformat(timespec="seconds"),
            "t_end": datetime.now().isoformat(timespec="seconds"),
            "rc": rc,
            "status": status or ("success" if rc == 0 else ("stopped" if rc == -9 else "failed")),
            "failure_reason": failure_reason,
            "attempts": attempts,
            "elapsed_s": elapsed_s,
            "products": products,
        }
        write_param_record(outdir, record)
    except Exception as exc:  # noqa: BLE001
        log_cb(f"[警告] 写参数档案失败（不影响计算结果）: {exc}")


# ============================================================
# 模块 A 警告（需求 3）：解析已编译主程序的硬编码 E/N、Tgas
# ============================================================

def parse_compiled_conditions(main_src):
    """解析主程序源码 REDUCED_FIELD=/GAS_TEMPERATURE= 字面值 → (en_str, tg_str)；无匹配为 None。"""
    try:
        text = Path(main_src).read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None, None
    en = re.search(r"REDUCED_FIELD=([0-9][0-9.eEdD+-]*)", text)
    tg = re.search(r"GAS_TEMPERATURE=([0-9][0-9.eEdD+-]*)", text)
    return (en.group(1) if en else None, tg.group(1) if tg else None)


def module_a_warning(flavor, main_src, en_val, tg_val):
    """模块 A 不一致警告（仅 hong 型硬编码主程序；命令行直传型返回 None）。
    返回警告字符串或 None。"""
    if flavor != "hong" or not main_src:
        return None
    en_c, tg_c = parse_compiled_conditions(main_src)
    msgs = []
    try:
        if en_c and abs(_to_float(en_c) - _to_float(en_val)) > 1e-9:
            msgs.append(f"E/N 编译值 {en_c}")
        if tg_c and abs(_to_float(tg_c) - _to_float(tg_val)) > 1e-9:
            msgs.append(f"Tgas 编译值 {tg_c}")
    except (ValueError, TypeError):
        return None
    if msgs:
        return "⚠ 与已编译主程序不一致（" + "、".join(msgs) + "）——需勾选 recompile 重编译后生效"
    return None


# ============================================================
# v3 需求 3：auto 型固化组分解析（第 6 步只读摘要的数据源）
# ============================================================

def parse_auto_composition(main_src):
    """从 main_auto.F90 解析生成时固化的气体组分：[(species, frac_str), ...]
    （只取 "<分数>*ntot" 形式；M/S 自动初始化行不算用户组分）。解析失败返回 []。"""
    try:
        text = Path(main_src).read_text(encoding="utf-8", errors="ignore")
    except (OSError, TypeError):
        return []
    return [(m.group(1), m.group(2)) for m in re.finditer(
        r"set_density\('([^']+)',\s*([0-9][0-9.eEdD+-]*)\s*\*\s*ntot\s*\)", text)]


def fmt_composition(comp):
    """[(species, frac_str), ...] → 'N2 33.33% / H2 66.67%'；空返回 ''。"""
    parts = []
    for sp, fr in comp or []:
        try:
            parts.append(f"{sp} {_to_float(fr) * 100:g}%")
        except (ValueError, TypeError):
            parts.append(sp)
    return " / ".join(parts)


# ============================================================
# v3 需求 1：第 1 步检测 → 第 2/3 步联动（手动优先）
# ============================================================

def auto_link_updates(slots, xs_manual=False, kinet_manual=False):
    """按槽位检测结果计算第 2/3 步应自动带入的文件路径。
    返回 {"bolsigdb": path|None, "kinet": path|None}：
    仅当槽位为 auto（第 1 步自动识别）且该步未被手动指定（★）过时给出路径——
    手动指定过的绝不覆盖。"""
    out = {"bolsigdb": None, "kinet": None}
    for key, manual in (("bolsigdb", xs_manual), ("kinet", kinet_manual)):
        st = (slots or {}).get(key, {})
        if not manual and st.get("status") == "auto" and st.get("path"):
            out[key] = st["path"]
    return out


# ============================================================
# 输入文件中心（需求 2）：槽位模型（自动识别 + 逐项覆盖/重置）
# ============================================================

# (key, 显示名, 是否必需)
SLOT_DEFS = [
    ("bolsigdb",    "截面库 bolsigdb.dat（或 txt 组拼接）", True),
    ("kinet",       "动力学机理 kinet.inp", True),
    ("program",     "程序文件 preprocessor/dvode/bolsig 库", True),
    ("main",        "主程序 .F90（可在第 4 步生成）", True),
    ("zdplaskin_m", "zdplaskin_m.F90（可选加速：跳过预处理）", False),
    ("entropy",     "熵表 *.DAT（按机理自动判断：含熵机理时必需，其余可缺）", False),
    ("libquadmath", "libquadmath-0.dll（缺失可从 gfortran 回退）", False),
]
SLOT_KEYS = [k for k, _n, _r in SLOT_DEFS]


def zdplaskin_m_species_max(zdm_path):
    """从 zdplaskin_m.F90 解析 species_max 参数值；解析失败返回 None。"""
    try:
        text = Path(zdm_path).read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None
    m = re.search(r"species_max\s*=\s*(\d+)", text)
    return int(m.group(1)) if m else None


def detect_slots(directory, overrides=None):
    """检测输入文件槽位。overrides: {slot_key: 手动指定路径}（"重置为自动"即删除对应键）。
    返回 {slot_key: {"status": auto|manual|missing|optional, "path": str|None,
    "detail": str, "ok": bool, "required": bool, "name": str}}。"""
    d = Path(directory) if directory else None
    dd = d if (d and d.is_dir()) else None
    overrides = dict(overrides or {})
    slots = {}

    def file_ok(p):
        return p and Path(p).is_file()

    for key, name, required in SLOT_DEFS:
        ov = overrides.get(key)
        st = {"name": name, "required": required, "status": "missing",
              "path": None, "detail": "", "ok": False}
        if key == "program":
            # 多文件槽：算例目录内齐备，或覆盖目录内齐备
            src_dir = Path(ov) if ov else dd
            if src_dir and src_dir.is_dir():
                have = [n for n in DISTRO_REQUIRED if (src_dir / n).is_file()]
                lack = [n for n in DISTRO_REQUIRED if not (src_dir / n).is_file()]
                if not lack:
                    st.update(status="manual" if ov else "auto", path=str(src_dir),
                              ok=True, detail="必需 4 件齐备")
                else:
                    st.update(status="manual" if ov else "missing", path=str(src_dir),
                              detail="缺: " + ", ".join(lack))
            elif ov:
                st.update(status="manual", path=str(ov), detail="覆盖目录不存在")
        elif key == "entropy":
            n_dat = len(list(dd.glob("*.DAT"))) if dd else 0
            if ov and file_ok(ov):
                st.update(status="manual", path=str(ov), ok=True, detail="手动指定")
            elif n_dat:
                st.update(status="auto", path=str(dd), ok=True,
                          detail=f"检测到 {n_dat} 个 .DAT")
            else:
                st.update(status="optional", detail="未检测到 .DAT（按机理自动判断：含熵机理时必需，如 515 版需 4 个）")
                st["ok"] = True  # 可选槽：缺失不卡门禁
        elif key == "libquadmath":
            if ov and file_ok(ov):
                st.update(status="manual", path=str(ov), ok=True, detail="手动指定")
            elif dd and (dd / "libquadmath-0.dll").is_file():
                st.update(status="auto", path=str(dd / "libquadmath-0.dll"), ok=True)
            else:
                gb = find_gfortran_bin()
                if gb and (gb / "libquadmath-0.dll").is_file():
                    st.update(status="optional", path=str(gb / "libquadmath-0.dll"), ok=True,
                              detail="算例目录缺失，可自 gfortran bin 回退")
                else:
                    st.update(status="optional", ok=True, detail="缺失且 gfortran 回退不可用")
        elif key == "zdplaskin_m":
            zm = None
            if ov and file_ok(ov):
                zm = Path(ov)
            elif dd:
                for n in ("zdplaskin_m.F90", "zdplaskin_m.f90"):
                    if (dd / n).is_file():
                        zm = dd / n
                        break
            if zm:
                st.update(status="manual" if ov else "auto", path=str(zm), ok=True)
                # 与 kinet.inp 物种数匹配检查
                kin = None
                for n in ("kinet.inp", "kinetics.inp"):
                    if dd and (dd / n).is_file():
                        kin = dd / n
                        break
                if kin:
                    n_zm = zdplaskin_m_species_max(zm)
                    n_kin = len(parse_kinet_species(kin))
                    if n_zm and n_kin and n_zm != n_kin:
                        st["detail"] = (f"⚠ species_max={n_zm} 与 kinet.inp 的 {n_kin} 个物种不匹配"
                                        "——建议重新执行完整预处理与构建")
                    else:
                        st["detail"] = f"species_max={n_zm}，与机理匹配" if n_zm else ""
            else:
                st.update(status="optional", ok=True,
                          detail="缺失时 Python 原生构建会先运行 preprocessor 生成")
        else:
            # 单文件必需槽：bolsigdb / kinet / main
            names = {"bolsigdb": ("bolsigdb.dat",), "kinet": ("kinet.inp", "kinetics.inp")}[key] \
                if key in ("bolsigdb", "kinet") else ()
            p = None
            if ov and file_ok(ov):
                p = Path(ov)
            elif dd:
                if key == "main":
                    p = _find_main_src(dd)
                else:
                    for n in names:
                        if (dd / n).is_file():
                            p = dd / n
                            break
            if p:
                st.update(status="manual" if ov else "auto", path=str(p), ok=True)
                if key == "kinet":
                    ns = len(parse_kinet_species(p))
                    rs = parse_kinet_reaction_stats(p)
                    if ns >= 2 and rs["parsed"]:
                        st["detail"] = (f"{ns} 个物种 · {rs['total']} 个反应"
                                        f"（BOLSIG+ {rs['bolsig']}，固定/其他 {rs['fixed']}）")
                    elif ns >= 2:
                        st["detail"] = f"{ns} 个物种 · ⚠ REACTIONS 段解析异常"
                    else:
                        st["detail"] = "⚠ SPECIES 段解析异常"
                    if ns < 2:
                        st["ok"] = False
            elif ov:
                st.update(status="manual", path=str(ov), detail="覆盖文件不存在")
        slots[key] = st
    return slots


def materialize_slots(directory, overrides):
    """把手动覆盖（★）槽位的文件复制进算例目录（原始文件只读复制）。
    program 槽位从覆盖目录复制必需+可选程序文件。返回复制文件名列表。"""
    d = Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    copied = []
    for key, ov in (overrides or {}).items():
        if not ov:
            continue
        src = Path(ov)
        if key == "program":
            if src.is_dir():
                for n in DISTRO_REQUIRED + DISTRO_OPTIONAL:
                    if (src / n).is_file() and not (d / n).is_file():
                        shutil.copy2(src / n, d / n)
                        copied.append(n)
        elif key == "entropy":
            if src.is_file() and not (d / src.name).is_file():
                shutil.copy2(src, d / src.name)
                copied.append(src.name)
        elif key in ("bolsigdb", "kinet", "main", "zdplaskin_m", "libquadmath"):
            if src.is_file():
                dst = d / ("kinet.inp" if key == "kinet" else
                           "bolsigdb.dat" if key == "bolsigdb" else src.name)
                if not dst.is_file():
                    shutil.copy2(src, dst)
                    copied.append(dst.name)
    return copied


# ============================================================
# 依赖解析（需求 6 便携性）：config → PATH → 常见安装位置 回退链
# ============================================================

COMMON_GFORTRAN_DIRS = [
    r"F:\Softwares\Ming64\mingw64\bin",
    r"C:\mingw64\bin",
    r"C:\msys64\mingw64\bin",
    r"C:\Program Files\mingw-w64\bin",
]


def _resolve_chain(config_val, which_name, commons, which_fn=None, file_pred=None):
    """回退链：config 值 → PATH(where) → 常见安装位置。全失败返回 None。
    which_fn/file_pred 可注入以便测试 mock。"""
    which_fn = which_fn or shutil.which
    file_pred = file_pred or (lambda p: Path(p).is_file())
    if config_val and file_pred(config_val):
        return config_val
    w = which_fn(which_name)
    if w:
        return w
    for c in commons:
        if file_pred(c):
            return c
    return None


def resolve_gfortran_dir(cfg=None, which_fn=None, file_pred=None, env=None):
    """定位 gfortran bin 目录：config['gfortran_dir'] → GFORTRAN_DIR 环境变量
    → PATH → 常见安装位置。返回目录路径字符串或 None。"""
    cfg = cfg or {}
    env = os.environ if env is None else env
    file_pred = file_pred or (lambda p: Path(p).is_file())
    exe_pred = lambda d: file_pred(str(Path(d) / "gfortran.exe"))
    if cfg.get("gfortran_dir") and exe_pred(cfg["gfortran_dir"]):
        return cfg["gfortran_dir"]
    gd = env.get("GFORTRAN_DIR")
    if gd and exe_pred(gd):
        return gd
    w = (which_fn or shutil.which)("gfortran")
    if w:
        return str(Path(w).parent)
    for c in COMMON_GFORTRAN_DIRS:
        if exe_pred(c):
            return c
    return None


def self_check(cfg=None):
    """启动环境自检 → [(名称, ok, 说明)]。"""
    cfg = cfg or {}
    checks = []
    gf = resolve_gfortran_dir(cfg)
    checks.append(("gfortran", gf is not None,
                   gf or "未找到——编译不可用（config/GFORTRAN_DIR/PATH/常见位置均无）"))
    tpl_ok = (TEMPLATE_MAIN.is_file() and TEMPLATE_MAIN_PULSE.is_file()
              and TEMPLATE_KINET.is_file())
    checks.append(("模板文件 templates", tpl_ok, str(TEMPLATE_DIR)))
    runtime_ok = (PYTHON_CLI.is_file() and (SCRIPT_DIR / "zdp_runtime.py").is_file())
    checks.append(("Python 原生运行时", runtime_ok,
                   str(PYTHON_CLI) if runtime_ok else "缺 zdp_cli.py 或 zdp_runtime.py"))
    checks.append(("CPython 解释器", sys.implementation.name == "cpython", sys.executable))
    distro = cfg.get("distro_dir", "")
    if distro:
        ok = Path(distro).is_dir()
        checks.append(("发行包目录（config 记录）", True if ok else None,
                       distro if ok else f"{distro} —— 路径已失效，请重新选择"))
    return checks


def wizard_workdir(case_dir, batch):
    """需求 4：批次名非空 → <算例目录>/<批次名>/（新建算例，构建产物归此）；
    空 → 算例目录本身（已成型目录原地模式）。"""
    case_dir = str(case_dir).strip()
    batch = str(batch).strip()
    if case_dir and batch:
        return str(Path(case_dir) / batch)
    return case_dir


# ============================================================
# 子进程执行引擎（GUI 线程与无界面测试共用）
# ============================================================

def run_subprocess(args, cwd, log_cb=print, stop_holder=None):
    """运行子进程并逐行回传输出。返回退出码；被中止返回 -9。"""
    log_cb("$ " + " ".join(str(a) for a in args))
    env = dict(os.environ)
    try:
        proc = subprocess.Popen(
            [str(a) for a in args], cwd=str(cwd),
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace",
            env=env,
        )
    except Exception as exc:  # noqa: BLE001
        log_cb(f"[错误] 无法启动子进程: {exc}")
        return 127
    if stop_holder is not None:
        with stop_holder.setdefault("lock", threading.Lock()):
            stop_holder.setdefault("procs", set()).add(proc)
            stop_holder["proc"] = proc  # 保留单任务旧接口兼容性
    assert proc.stdout is not None
    for line in proc.stdout:
        if stop_holder is not None and stop_holder.get("killed"):
            break
        log_cb(line.rstrip("\n"))
    rc = proc.wait()
    if stop_holder is not None:
        with stop_holder.setdefault("lock", threading.Lock()):
            stop_holder.setdefault("procs", set()).discard(proc)
            if stop_holder.get("proc") is proc:
                stop_holder["proc"] = next(iter(stop_holder["procs"]), None)
        if stop_holder.get("killed"):
            return -9
    return rc


def kill_process_tree(stop_holder):
    """中止当前所有子进程及各自进程树（并行扫描时一次停止整批活动点）。"""
    if not stop_holder:
        return
    stop_holder["killed"] = True
    with stop_holder.setdefault("lock", threading.Lock()):
        procs = list(stop_holder.get("procs", set()))
        if stop_holder.get("proc") is not None and stop_holder["proc"] not in procs:
            procs.append(stop_holder["proc"])
    for proc in procs:
        try:
            subprocess.run(["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                           capture_output=True, timeout=10)
        except Exception:  # noqa: BLE001
            try:
                proc.kill()
            except Exception:  # noqa: BLE001
                pass


def pulse_model_warnings(directory, flavor, params):
    """返回脉冲结果解释前必须显示的、由输入可直接判定的模型边界。"""
    if flavor not in ("pulse", "autopulse"):
        return []
    warnings = []
    ne_mode = str(params.get("ne_mode", "fix")).split("@", 1)[0].strip().lower()
    if flavor == "autopulse" or ne_mode == "fix":
        warnings.append("当前采用固定电子密度；该设置不强制电荷自洽，不能单独作为自洽放电解解释。")
    try:
        kinet = next((Path(directory) / n for n in ("kinet.inp", "kinetics.inp")
                      if (Path(directory) / n).is_file()), None)
        text = kinet.read_text(encoding="utf-8", errors="ignore") if kinet else ""
    except (OSError, StopIteration):
        text = ""
    if "H2(RYDBERG_SUM)" in text and not re.search(
            r"^\s*H2\(RYDBERG_SUM\)\s*=>\s*H2\b", text, re.M | re.I):
        warnings.append("H2(RYDBERG_SUM) 未检出直接壁弛豫到 H2 的反应；NH 来源应做 Rydberg 寿命/关断敏感性分析。")
    return warnings


def execute_pipeline(directory, params, do_build=True, do_run=True,
                      log_cb=print, stop_holder=None, flat_out=False):
    """Python 原生构建 → （可选）运行 → 写参数档案。返回最终退出码。
    flat_out=True 时输出到 <构建目录>/<参数组合>/（向导新建算例：构建目录本身即批次）；
    默认输出到 runs/<批次名>/<参数组合>/。
    GUI 工作线程与无界面链路测试都走这个入口。"""
    d = Path(directory)
    if not PYTHON_CLI.is_file():
        log_cb(f"[错误] Python 原生运行入口缺失: {PYTHON_CLI}")
        return 126
    log_cb(f"[信息] 使用 Python 原生运行时: {PYTHON_CLI.name}")
    rc = 0
    if do_build:
        log_cb("[信息] ===== 阶段 1/2：编译 =====")
        rc = run_subprocess(build_build_args(), d, log_cb, stop_holder)
        if rc != 0:
            log_cb(f"[失败] 编译退出码 {rc}")
            return rc
        log_cb("[信息] 编译成功")
    if do_run:
        flavor = detect_flavor(d) or "hong"
        out_rel = out_rel_for(params, flavor, flat=flat_out)
        log_cb(f"[信息] ===== 阶段 {'2/2' if do_build else '1/1'}：运行（{flavor} 型）→ {out_rel} =====")
        for warning in pulse_model_warnings(d, flavor, params):
            log_cb(f"[模型边界] {warning}")
        args = build_run_args("", params, flavor, out_rel=out_rel)
        t0 = time.time()
        rc = run_subprocess(args, d, log_cb, stop_holder)
        status, failure_reason = None, None
        if rc == 0:
            metrics = extract_metrics(flavor, d / out_rel, params["tag"])
            status, failure_reason = validate_run_output(
                flavor, d / out_rel, params["tag"], params, metrics)
            if status != "success":
                log_cb(f"[输出验证] {status}: {failure_reason}")
        finalize_run_record(d, out_rel, flavor, params, args, t0, rc, log_cb,
                            status=status, failure_reason=failure_reason)
        if rc == 0:
            log_cb(f"[成功] 运行完成，产物在 {out_rel}/")
        else:
            log_cb(f"[失败] 运行退出码 {rc}")
    return rc


# ============================================================
# 快速绘图（增强项：任何失败只弹提示，不影响主流程）
# ============================================================

def _load_density_series(file: Path):
    """读取密度-时间序列。支持 qt_species_density.txt 与本项目 *_output_*.csv。
    返回 (t, {物种名: 序列})。"""
    if file.suffix.lower() == ".csv":
        with open(file, encoding="utf-8", errors="replace", newline="") as fh:
            rows = list(csv.reader(fh))
        if not rows:
            raise ValueError("CSV 文件为空")
        hdr = [h.strip() for h in rows[0]]
        if not hdr:
            raise ValueError("CSV 缺少表头")
        tcol = hdr.index("time_s") if "time_s" in hdr else 0
        start = 6 if len(hdr) > 6 and hdr[1].startswith("Te") else 1
        indices = list(range(start, len(hdr)))
        t, series = [], {hdr[i]: [] for i in indices}
        for row in rows[1:]:
            if len(row) <= tcol:
                continue
            try:
                time_value = float(row[tcol])
            except ValueError:
                continue
            t.append(time_value)
            for i in indices:
                try:
                    series[hdr[i]].append(float(row[i]))
                except (IndexError, ValueError):
                    series[hdr[i]].append(_NAN)
        return t, series
    # qt_species_density.txt：首行表头（可带 #），空白分隔，首列时间
    with open(file, encoding="utf-8", errors="replace") as fh:
        lines = [ln for ln in fh if ln.strip()]
    if len(lines) < 2:
        raise ValueError("物种密度文件缺少数据行")
    hdr = lines[0].lstrip("#").split()
    if len(hdr) < 2:
        raise ValueError("物种密度文件表头无物种列")
    data = []
    for line in lines[1:]:
        try:
            data.append([float(x) for x in line.split()])
        except ValueError:
            continue
    if not data:
        raise ValueError("物种密度文件没有可解析数值")
    ncol = min(len(hdr), min(len(row) for row in data))
    t = [r[0] for r in data]
    series = {hdr[i]: [r[i] for r in data] for i in range(1, ncol)}
    return t, series


def find_plot_targets(outdir: Path):
    """返回可绘制序列候选，兼容单算例和 runs/<批次>/<点>/ 嵌套扫描目录。

    按最近修改时间降序排列；quick_plot 会再过滤掉无法解析或没有正数点的文件。
    """
    outdir = Path(outdir)
    if not outdir.is_dir():
        return []
    files = []
    patterns = ("qt_species_density.txt", "*_output_*.csv", "*output*.csv", "pulse_series_*.csv")
    for pattern in patterns:
        files.extend(p for p in outdir.rglob(pattern) if p.is_file())
    unique = {p.resolve(): p for p in files}
    return sorted(unique.values(), key=lambda p: p.stat().st_mtime, reverse=True)


def _positive_finite_points(t, ys):
    """为双对数坐标保留同一行中均为有限正数的 (time, density) 点。"""
    points = []
    for x, y in zip(t, ys):
        try:
            x, y = float(x), float(y)
        except (TypeError, ValueError):
            continue
        if math.isfinite(x) and math.isfinite(y) and x > 0 and y > 0:
            points.append((x, y))
    return points


def quick_plot(parent, outdir: Path):
    """在弹窗中画主要物种密度-时间对数曲线。返回 True/False。"""
    candidates = find_plot_targets(outdir)
    if not candidates:
        from tkinter import messagebox
        messagebox.showinfo("快速绘图", f"{outdir} 及其扫描子目录内未找到物种序列文件", parent=parent)
        return False
    try:
        sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
        import matplotlib
        matplotlib.use("TkAgg")
        from matplotlib import pyplot as plt
        try:
            from daimon_runtime import setup_plot  # 中文字体
            setup_plot()
        except Exception:  # noqa: BLE001
            plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
            plt.rcParams["axes.unicode_minus"] = False
        target = None
        top = []
        parse_errors = []
        for candidate in candidates:
            try:
                t, series = _load_density_series(candidate)
                valid = [(name, _positive_finite_points(t, ys)) for name, ys in series.items()]
                valid = [(name, points) for name, points in valid if points]
                if valid:
                    target = candidate
                    top = sorted(valid, key=lambda item: max(y for _x, y in item[1]), reverse=True)[:6]
                    break
            except Exception as exc:  # noqa: BLE001
                parse_errors.append(f"{candidate.name}: {exc}")
        if target is None:
            detail = "；".join(parse_errors[:2]) or "全部序列仅含 0、NaN 或 Inf"
            raise ValueError("未找到可用于对数坐标的有限正值数据（" + detail + "）")
        fig, ax = plt.subplots(figsize=(8, 5))
        for name, points in top:
            ax.plot([x for x, _y in points], [y for _x, y in points], label=name)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("时间 (s)"); ax.set_ylabel("密度 (cm$^{-3}$)")
        try:
            rel_target = target.relative_to(Path(outdir))
        except ValueError:
            rel_target = target.name
        ax.set_title(f"主要物种演化 — {rel_target}")
        ax.legend(); fig.tight_layout()
        import tkinter as tk
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        win = tk.Toplevel(parent)
        win.title(f"快速绘图 — {outdir.name}")
        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        return True
    except Exception as exc:  # noqa: BLE001
        from tkinter import messagebox
        messagebox.showerror("快速绘图失败", str(exc), parent=parent)
        return False


def scan_plot(parent, rows, axis_keys, metric_names):
    """扫描结果绘图：1D → 各指标-参数曲线；2D → 各指标热力图；3D → 只提示（看表）。
    返回 True/False。NaN 点在图中跳过。"""
    from tkinter import messagebox
    if not rows:
        messagebox.showinfo("扫描绘图", "没有扫描结果可绘制", parent=parent)
        return False
    ndim = len(axis_keys)
    if ndim >= 3:
        messagebox.showinfo("扫描绘图", "3D 扫描不出图，请直接查看汇总表 / scan_summary.csv", parent=parent)
        return False
    try:
        sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
        import matplotlib
        matplotlib.use("TkAgg")
        from matplotlib import pyplot as plt
        try:
            from daimon_runtime import setup_plot  # 中文字体
            setup_plot()
        except Exception:  # noqa: BLE001
            plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
            plt.rcParams["axes.unicode_minus"] = False
        import math
        metrics = [m for m in metric_names
                   if any(isinstance(r.get(m), float) and r.get(m) == r.get(m) for r in rows)]
        if not metrics:
            messagebox.showinfo("扫描绘图", "所有指标均为 NaN，无可绘制数据", parent=parent)
            return False
        ncol = min(3, len(metrics))
        nrow = math.ceil(len(metrics) / ncol)
        fig, axes_arr = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.4 * nrow), squeeze=False)
        k1 = axis_keys[0]
        if ndim == 1:
            xs_all = sorted({r[k1] for r in rows})
            for j, m in enumerate(metrics):
                ax = axes_arr[j // ncol][j % ncol]
                ys = {r[k1]: r.get(m, _NAN) for r in rows}
                pts = [(x, ys[x]) for x in xs_all if ys.get(x) == ys.get(x)]
                if pts:
                    ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", ms=4)
                ax.set_xlabel(k1); ax.set_title(m); ax.grid(alpha=0.3)
        else:  # 2D 热力图
            k2 = axis_keys[1]
            xs = sorted({r[k1] for r in rows})
            ys_ax = sorted({r[k2] for r in rows})
            for j, m in enumerate(metrics):
                ax = axes_arr[j // ncol][j % ncol]
                grid = [[_NAN] * len(xs) for _ in ys_ax]
                for r in rows:
                    v = r.get(m, _NAN)
                    if v == v:
                        grid[ys_ax.index(r[k2])][xs.index(r[k1])] = v
                masked = [[float("nan") if v != v else v for v in row] for row in grid]
                im = ax.imshow(masked, origin="lower", aspect="auto",
                               extent=[min(xs), max(xs), min(ys_ax), max(ys_ax)])
                ax.set_xlabel(k1); ax.set_ylabel(k2); ax.set_title(m)
                fig.colorbar(im, ax=ax)
        for j in range(len(metrics), nrow * ncol):
            axes_arr[j // ncol][j % ncol].axis("off")
        fig.suptitle("参数扫描 — " + " × ".join(axis_keys))
        fig.tight_layout()
        import tkinter as tk
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        win = tk.Toplevel(parent)
        win.title("扫描结果绘图")
        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        return True
    except Exception as exc:  # noqa: BLE001
        from tkinter import messagebox
        messagebox.showerror("扫描绘图失败", str(exc), parent=parent)
        return False


# ============================================================
# GUI（向导式：从零构建算例并计算）
# ============================================================

# 需求 5 · 样式常量：红色只用于关键输入；警告用橙色
FONT_KEY = ("Segoe UI", 11, "bold")   # 关键输入：11pt 加粗
COLOR_KEY = "#8B1A1A"                 # 关键输入深红
COLOR_WARN = "#B26A00"                # 警告橙
COLOR_OK = "#2E7D32"
COLOR_MANUAL = "#1565C0"
COLOR_MISS = "#C62828"
COLOR_GRAY = "#9E9E9E"


def make_scrollable_frame(tab, tk, ttk):
    """v3 需求 2：把 Notebook 页变成可垂直滚动容器（Canvas + 垂直 Scrollbar + 鼠标滚轮），
    窗口不最大化也能滚到底。返回 (内部内容 Frame, 外层 Canvas)。
    Text/Listbox/Treeview/Canvas 控件保留自身滚轮行为，不抢占。"""
    canvas = tk.Canvas(tab, highlightthickness=0)
    vsb = ttk.Scrollbar(tab, orient="vertical", command=canvas.yview)
    inner = ttk.Frame(canvas)
    win_id = canvas.create_window((0, 0), window=inner, anchor="nw")
    inner.bind("<Configure>",
               lambda _e: canvas.configure(scrollregion=canvas.bbox("all")))
    canvas.bind("<Configure>",
                lambda e: canvas.itemconfigure(win_id, width=e.width))
    canvas.configure(yscrollcommand=vsb.set)
    canvas.pack(side="left", fill="both", expand=True)
    vsb.pack(side="right", fill="y")

    def _on_wheel(e):
        canvas.yview_scroll(int(-e.delta / 120), "units")

    def _bind_wheel(widget):
        if not isinstance(widget, (tk.Text, tk.Canvas, tk.Listbox, ttk.Treeview)):
            widget.bind("<MouseWheel>", _on_wheel)
        for ch in widget.winfo_children():
            _bind_wheel(ch)

    inner.bind_wheel = lambda: _bind_wheel(inner)  # 内容全部建完后调用一次
    return inner, canvas


def make_gui(debug=False):
    """创建 GUI。

    ``debug=True`` 时仍完整初始化 GUI（便于复用同一条运行链路），但隐藏
    主窗口，并把编译/运行日志同步写到启动该进程的终端。
    """
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, scrolledtext

    app = tk.Tk()
    app.title("ZDPlasKin 通用计算控制台")
    app.geometry("1080x1180")
    if debug:
        app.withdraw()

    q = queue.Queue()
    stop_holder = {"proc": None, "procs": set(), "lock": threading.Lock(), "killed": False}
    state = {"running": False, "t0": None, "flavor": None, "last_outdir": None,
             "scan_rows": None, "scan_axes": None, "scan_metrics": None,
             "scan_summary": None, "scan_progress": None, "composition_dirty": False}

    def emit_log(kind, line):
        """向 GUI 日志队列投递；调试模式额外实时输出至终端。"""
        line = str(line)
        if debug:
            channel = "编译" if kind == "log5" else "运行"
            print(f"[GUI 调试·{channel}] {line}", flush=True)
        q.put((kind, line))

    cfg = load_config()

    # ---------- 启动环境自检（需求 6）：记录解析结果供 bat/后续运行复用 ----------
    try:
        (SCRIPT_DIR / "python_path.txt").write_text(sys.executable, encoding="utf-8")
    except OSError:
        pass
    _gf0 = resolve_gfortran_dir(cfg)
    _cfg_dirty = False
    if _gf0 and cfg.get("gfortran_dir") != _gf0:
        cfg["gfortran_dir"] = _gf0
        _cfg_dirty = True
    if _cfg_dirty:
        save_config(cfg)

    # ---------- 模块 E 输出控制变量（提前定义：第 4 步生成主程序要读取） ----------
    out_vars = {}
    for _k, _dv in [("out_species", True), ("out_te", True), ("out_rates", True),
                    ("out_rates_t", True), ("out_energy", True), ("out_console", True),
                    ("out_summary", True)]:
        out_vars[_k] = tk.BooleanVar(value=_dv)

    # ---------- 顶部 Notebook 六步向导 ----------
    nb = ttk.Notebook(app)
    nb.pack(fill="both", expand=True, padx=4, pady=4)
    tab1 = ttk.Frame(nb); tab2 = ttk.Frame(nb); tab3 = ttk.Frame(nb)
    tab4 = ttk.Frame(nb); tab5 = ttk.Frame(nb); tab6 = ttk.Frame(nb)
    nb.add(tab1, text=" 1 输入文件中心 ")
    nb.add(tab2, text=" 2 截面库来源 ")
    nb.add(tab3, text=" 3 机理 kinet.inp ")
    nb.add(tab4, text=" 4 主程序生成 ")
    nb.add(tab5, text=" 5 编译 ")
    nb.add(tab6, text=" 6 运行与扫描 ")
    # v3 需求 2：第 6 步内容超高，包一层可滚动容器（其余步高度均在窗口内，无需处理）
    tab6, tab6_canvas = make_scrollable_frame(tab6, tk, ttk)

    lights = {}

    def make_light(parent, step):
        lbl = ttk.Label(parent, text="● 未开始", foreground=COLOR_GRAY)
        lbl.pack(side="left", padx=8, pady=4)
        lights[step] = lbl
        return lbl

    def set_light(step, st, extra=""):
        txt = {None: "● 未开始", True: "● 通过", False: "● 未通过"}[st] + extra
        col = {None: COLOR_GRAY, True: COLOR_OK, False: COLOR_MISS}[st]
        lights[step].config(text=txt, foreground=col)

    # ========================================================
    # 第 1 步：输入文件中心（需求 2：槽位表 + 逐项覆盖/重置；需求 4：批次名）
    # ========================================================
    frm1_dir = ttk.LabelFrame(tab1, text="算例目录（关键输入；新建或选择。已有完整构建目录可直接跳到第 5/6 步）")
    frm1_dir.pack(fill="x", padx=8, pady=6)
    var_dir = tk.StringVar()
    tk.Entry(frm1_dir, textvariable=var_dir, font=FONT_KEY, fg=COLOR_KEY).pack(
        side="left", fill="x", expand=True, padx=4, pady=4)

    def browse_case():
        d = filedialog.askdirectory(title="选择算例目录")
        if d:
            var_dir.set(d)

    def new_case():
        from tkinter import simpledialog
        d = filedialog.askdirectory(title="选择父目录（在其下新建算例文件夹）")
        if not d:
            return
        name = simpledialog.askstring("新建算例", "文件夹名：", parent=app)
        if name:
            p = Path(d) / name
            p.mkdir(parents=True, exist_ok=True)
            var_dir.set(str(p))
            if not var_batch.get().strip():
                var_batch.set("case_" + time.strftime("%H%M%S"))

    ttk.Button(frm1_dir, text="浏览…", command=browse_case).pack(side="left", padx=4)
    ttk.Button(frm1_dir, text="新建文件夹…", command=new_case).pack(side="left", padx=4)

    frm1_batch = ttk.LabelFrame(tab1, text="批次名（需求 4：新建算例的构建子目录；留空 = 已成型目录原地编译运行）")
    frm1_batch.pack(fill="x", padx=8, pady=4)
    var_batch = tk.StringVar(value="")
    ttk.Entry(frm1_batch, textvariable=var_batch, width=28).pack(side="left", padx=4, pady=4)
    ttk.Label(frm1_batch,
              text="填写后：程序文件/机理/截面/主程序/exe 全部构建到 <算例目录>/<批次名>/，运行输出到其下 <参数组合>/；原始输入只读复制不覆盖",
              foreground="#37474F").pack(side="left", padx=8)

    def eff_dir():
        return wizard_workdir(var_dir.get(), var_batch.get())

    frm1_distro = ttk.LabelFrame(tab1, text="ZDPlasKin 发行包目录（提供 preprocessor / dvode / bolsig 库；config 记录，启动时重新校验）")
    frm1_distro.pack(fill="x", padx=8, pady=6)
    var_distro = tk.StringVar(value=cfg.get("distro_dir", ""))
    ttk.Entry(frm1_distro, textvariable=var_distro).pack(side="left", fill="x", expand=True, padx=4, pady=4)

    def browse_distro():
        d = filedialog.askdirectory(title="选择 ZDPlasKin 发行包目录")
        if d:
            var_distro.set(d)

    ttk.Button(frm1_distro, text="浏览…", command=browse_distro).pack(side="left", padx=4)

    # ---------- 输入文件中心槽位表 ----------
    frm_slots = ttk.LabelFrame(tab1, text="输入文件中心 —— 严格最小 = 截面库 + 机理 + 程序文件（✓ 自动 / ★ 手动 / ✗ 缺失 / ○ 可选缺）")
    frm_slots.pack(fill="x", padx=8, pady=6)
    slot_overrides = {}   # key -> 手动指定路径（"重置为自动"即删除）
    slot_rows = {}        # key -> (status_lbl, path_lbl)

    def slot_browse(key):
        if key == "program":
            p = filedialog.askdirectory(title="选择程序文件所在目录（发行包或既有构建目录）")
        elif key == "entropy":
            p = filedialog.askopenfilename(title="选择熵表 .DAT", filetypes=[("熵表", "*.DAT"), ("全部", "*.*")])
        elif key == "bolsigdb":
            p = filedialog.askopenfilename(title="选择 bolsigdb.dat", filetypes=[("截面库", "*.dat"), ("全部", "*.*")])
        elif key == "kinet":
            p = filedialog.askopenfilename(title="选择 kinet.inp", filetypes=[("机理文件", "*.inp"), ("全部", "*.*")])
        elif key == "main":
            p = filedialog.askopenfilename(title="选择主程序", filetypes=[("Fortran", "*.F90 *.f90"), ("全部", "*.*")])
        else:
            p = filedialog.askopenfilename(title=f"选择 {key}")
        if p:
            slot_overrides[key] = p
            refresh_slots()

    def slot_reset(key):
        slot_overrides.pop(key, None)
        refresh_slots()

    for r, (key, name, required) in enumerate(SLOT_DEFS):
        st_lbl = tk.Label(frm_slots, text="✗", font=("Segoe UI", 10, "bold"), fg=COLOR_MISS, width=2)
        st_lbl.grid(row=r, column=0, padx=(6, 2), pady=1, sticky="w")
        is_key_input = key in ("bolsigdb", "kinet")
        name_lbl = tk.Label(frm_slots, text=name + ("［必需］" if required else "［可选］"),
                            font=FONT_KEY if is_key_input else ("Segoe UI", 9),
                            fg=COLOR_KEY if is_key_input else "#212121", anchor="w")
        name_lbl.grid(row=r, column=1, padx=4, sticky="w")
        path_lbl = tk.Label(frm_slots, text="—", fg="#37474F", anchor="w")
        path_lbl.grid(row=r, column=2, padx=4, sticky="w")
        ttk.Button(frm_slots, text="浏览…", width=7,
                   command=lambda k=key: slot_browse(k)).grid(row=r, column=3, padx=2)
        ttk.Button(frm_slots, text="重置为自动", width=10,
                   command=lambda k=key: slot_reset(k)).grid(row=r, column=4, padx=(2, 6))
        slot_rows[key] = (st_lbl, path_lbl)

    def _trunc(s, n=72):
        s = str(s)
        return s if len(s) <= n else "…" + s[-(n - 1):]

    def refresh_slots():
        slots = detect_slots(eff_dir(), slot_overrides)
        for key, _n, _req in SLOT_DEFS:
            st = slots[key]
            icon = {"auto": "✓", "manual": "★", "missing": "✗", "optional": "○"}[st["status"]]
            color = {"auto": COLOR_OK, "manual": COLOR_MANUAL,
                     "missing": COLOR_MISS, "optional": COLOR_GRAY}[st["status"]]
            slot_rows[key][0].config(text=icon, fg=color)
            txt = st["path"] or "未找到"
            if st["detail"]:
                txt += "  — " + st["detail"] if st["path"] else st["detail"]
            slot_rows[key][1].config(text=_trunc(txt))
        sync_step23(slots)

    # ---------- v3 需求 1：第 1 步检测 → 第 2/3 步联动（手动优先，三处共享槽位状态） ----------
    link_state = {"xs": None, "kinet": None, "kinet_path": None}  # None|"auto"|"manual"

    def sync_step23(slots):
        upd = auto_link_updates(slots,
                                xs_manual=(link_state["xs"] == "manual"),
                                kinet_manual=(link_state["kinet"] == "manual"))
        # 截面库 → 第 2 步
        if upd["bolsigdb"]:
            if var_dbfile.get().strip() != upd["bolsigdb"]:
                var_dbfile.set(upd["bolsigdb"])
            link_state["xs"] = "auto"
            lbl_src2.config(text="来源：来自第 1 步自动识别（可在本步【浏览…】改手动指定）",
                            foreground=COLOR_OK)
        elif link_state["xs"] == "manual":
            lbl_src2.config(text="来源：手动指定（★ 优先，自动识别不覆盖）",
                            foreground=COLOR_MANUAL)
        # 机理 → 第 3 步
        if upd["kinet"]:
            p = upd["kinet"]
            if link_state["kinet_path"] != p:
                try:
                    content = Path(p).read_text(encoding="utf-8", errors="replace")
                except OSError:
                    content = None
                if content is not None:
                    kinet_text.delete("1.0", "end")
                    kinet_text.insert("1.0", content)
                    link_state["kinet_path"] = p
            link_state["kinet"] = "auto"
            ns = len(parse_kinet_species(p))
            lbl_src3.config(text=f"来源：来自第 1 步自动识别（{Path(p).name}"
                                 + (f"，{ns} 个物种" if ns else "") + "）",
                            foreground=COLOR_OK)
            # 第 4 步尚未人工添加组分时，N2/H2 机理自动给出可编辑的 0.5/0.5 默认值。
            if not comp_rows:
                refresh_comp_species()
        elif link_state["kinet"] == "manual":
            lbl_src3.config(text="来源：手动指定（★ 优先，自动识别不覆盖）",
                            foreground=COLOR_MANUAL)

    def deploy_program():
        base = var_dir.get().strip()
        if not base:
            messagebox.showerror("错误", "请先选择算例目录", parent=app)
            return
        wd = eff_dir()
        Path(wd).mkdir(parents=True, exist_ok=True)
        notes = []
        # 需求 4：批次模式下把容器根的原始输入只读复制进工作目录（绝不覆盖原文件）
        if Path(wd) != Path(base):
            for name in ("kinet.inp", "kinetics.inp", "bolsigdb.dat"):
                src = Path(base) / name
                dst = Path(wd) / ("kinet.inp" if name == "kinetics.inp" else name)
                if src.is_file() and not dst.is_file():
                    shutil.copy2(src, dst)
                    notes.append(f"原始输入 {name} → 工作目录")
            for p in Path(base).glob("*.DAT"):
                if not (Path(wd) / p.name).is_file():
                    shutil.copy2(p, Path(wd) / p.name)
                    notes.append(f"原始输入 {p.name} → 工作目录")
        try:
            copied, skipped, missing = deploy_program_files(var_distro.get().strip(), wd)
        except Exception as exc:  # noqa: BLE001
            messagebox.showerror("复制失败", str(exc), parent=app)
            set_light(1, False)
            return
        mat = materialize_slots(wd, slot_overrides)
        if mat:
            notes.append("槽位覆盖已复制: " + ", ".join(mat))
        if missing:
            messagebox.showwarning("文件缺失", "发行包内未找到必需文件：\n" + "\n".join(missing), parent=app)
            set_light(1, False)
        else:
            msg = f"已复制 {len(copied)} 个程序文件，跳过已有 {len(skipped)} 个（目标: {wd}）"
            if notes:
                msg += "\n" + "\n".join(notes)
            messagebox.showinfo("完成", msg, parent=app)
        c = dict(load_config())
        c["distro_dir"] = var_distro.get().strip()
        save_config(c)
        refresh_wizard_status()

    frm1_deploy = ttk.Frame(tab1)
    frm1_deploy.pack(fill="x", padx=8, pady=2)
    ttk.Button(frm1_deploy, text="部署程序文件 + 应用槽位覆盖到工作目录",
               command=deploy_program).pack(side="left", padx=4)
    make_light(frm1_deploy, 1)

    # ========================================================
    # 第 2 步：截面库来源（供"输入文件中心"的截面槽位）
    # ========================================================
    frm2 = ttk.LabelFrame(tab2, text="截面库来源（二选一；bolsigdb.dat 由 BOLSIG+ 格式数据块拼接而成，发行包无官方转换工具）")
    frm2.pack(fill="x", padx=8, pady=6)
    var_xs = tk.StringVar(value="db")
    ttk.Radiobutton(frm2, text="直接提供 bolsigdb.dat", variable=var_xs, value="db").grid(row=0, column=0, sticky="w", padx=4)
    ttk.Radiobutton(frm2, text="由 BOLSIG+ 格式 .txt（LXCat 下载）拼接生成", variable=var_xs, value="merge").grid(row=1, column=0, sticky="w", padx=4)

    var_dbfile = tk.StringVar()
    ttk.Entry(frm2, textvariable=var_dbfile, width=60).grid(row=0, column=1, padx=4, pady=2)

    def xs_browse_db():
        p = filedialog.askopenfilename(title="选择 bolsigdb.dat")
        if p:
            var_dbfile.set(p)
            link_state["xs"] = "manual"  # v3 需求 1：手动指定（★）后自动识别不再覆盖
            lbl_src2.config(text="来源：手动指定（★ 优先，自动识别不覆盖）", foreground=COLOR_MANUAL)

    ttk.Button(frm2, text="浏览…", command=xs_browse_db).grid(row=0, column=2, padx=2)
    lbl_src2 = ttk.Label(tab2, text="来源：—（等待第 1 步识别或本步手动指定）", foreground=COLOR_GRAY)
    lbl_src2.pack(anchor="w", padx=8, pady=(0, 2))

    def xs_copy_db():
        src, d = var_dbfile.get().strip(), eff_dir()
        if not d:
            messagebox.showerror("错误", "请先在第 1 步选择算例目录", parent=app)
            return
        if not Path(src).is_file():
            messagebox.showerror("错误", "bolsigdb.dat 文件不存在", parent=app)
            return
        Path(d).mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, Path(d) / "bolsigdb.dat")
        messagebox.showinfo("完成", f"已复制 bolsigdb.dat 到 {d}", parent=app)
        refresh_wizard_status()

    ttk.Button(frm2, text="复制到工作目录", command=xs_copy_db).grid(row=0, column=3, padx=4)

    lst_txt = tk.Listbox(frm2, height=5, width=58)
    lst_txt.grid(row=1, column=1, rowspan=3, padx=4, pady=2, sticky="w")

    def xs_add_txt():
        fs = filedialog.askopenfilenames(title="选择 BOLSIG+ 格式 .txt（可多选）")
        for f in fs:
            lst_txt.insert("end", f)

    def xs_del_txt():
        for i in reversed(lst_txt.curselection()):
            lst_txt.delete(i)

    def xs_merge():
        files = list(lst_txt.get(0, "end"))
        d = eff_dir()
        if not d:
            messagebox.showerror("错误", "请先在第 1 步选择算例目录", parent=app)
            return
        if not files:
            messagebox.showerror("错误", "请先添加至少一个 .txt 截面文件", parent=app)
            return
        try:
            Path(d).mkdir(parents=True, exist_ok=True)
            out, n = merge_bolsigdb(files, Path(d) / "bolsigdb.dat")
        except Exception as exc:  # noqa: BLE001
            messagebox.showerror("拼接失败", str(exc), parent=app)
            set_light(2, False)
            return
        messagebox.showinfo("完成", f"已生成 {out.name}（约 {n} 个过程块）", parent=app)
        link_state["xs"] = "manual"  # v3 需求 1：用户选择自行拼接 = 手动来源
        lbl_src2.config(text="来源：手动拼接生成（★ 优先，自动识别不覆盖）", foreground=COLOR_MANUAL)
        refresh_wizard_status()

    frm2_btn = ttk.Frame(frm2); frm2_btn.grid(row=1, column=2, columnspan=2, rowspan=3, sticky="n")
    ttk.Button(frm2_btn, text="添加…", command=xs_add_txt).pack(pady=2)
    ttk.Button(frm2_btn, text="移除选中", command=xs_del_txt).pack(pady=2)
    ttk.Button(frm2_btn, text="拼接生成 bolsigdb.dat", command=xs_merge).pack(pady=2)
    frm2_light = ttk.Frame(tab2); frm2_light.pack(fill="x")
    make_light(frm2_light, 2)
    ttk.Label(tab2, text="说明：发行包自带一个标准 bolsigdb.dat（Ar 等示例）；本项目的 N2/H2 库即由 LXCat 下载的\n"
                         "BOLSIG+ 格式 txt 直接拼接而成（merge_bolsigdb.py 同款做法）。",
              foreground="#37474F", justify="left").pack(anchor="w", padx=8, pady=6)

    # ========================================================
    # 第 3 步：动力学机理 kinet.inp
    # ========================================================
    frm3_btn = ttk.Frame(tab3)
    frm3_btn.pack(fill="x", padx=8, pady=4)
    kinet_text = scrolledtext.ScrolledText(tab3, height=22, font=("Consolas", 9))
    lbl_species = ttk.Label(tab3, text="物种校验：未保存", foreground="#37474F")
    lbl_src3 = ttk.Label(tab3, text="来源：—（等待第 1 步识别或本步手动载入）", foreground=COLOR_GRAY)

    def k_load_template():
        if not TEMPLATE_KINET.is_file():
            messagebox.showerror("错误", f"模板缺失: {TEMPLATE_KINET}", parent=app)
            return
        kinet_text.delete("1.0", "end")
        kinet_text.insert("1.0", TEMPLATE_KINET.read_text(encoding="utf-8"))
        lbl_species.config(text="物种校验：已载入模板（未保存）")
        link_state["kinet"] = "manual"  # v3 需求 1：手动载入后自动识别不再覆盖
        lbl_src3.config(text="来源：手动指定（★ 优先，自动识别不覆盖）", foreground=COLOR_MANUAL)

    def k_load_existing():
        f = filedialog.askopenfilename(title="选择现有 kinet.inp",
                                       filetypes=[("机理文件", "*.inp"), ("全部", "*.*")])
        if f:
            kinet_text.delete("1.0", "end")
            kinet_text.insert("1.0", Path(f).read_text(encoding="utf-8", errors="replace"))
            lbl_species.config(text=f"物种校验：已载入 {Path(f).name}（未保存）")
            link_state["kinet"] = "manual"
            lbl_src3.config(text="来源：手动指定（★ 优先，自动识别不覆盖）", foreground=COLOR_MANUAL)

    def k_save():
        d = eff_dir()
        if not d:
            messagebox.showerror("错误", "请先在第 1 步选择算例目录", parent=app)
            return None
        Path(d).mkdir(parents=True, exist_ok=True)
        out = Path(d) / "kinet.inp"
        out.write_text(kinet_text.get("1.0", "end-1c"), encoding="utf-8")
        species = parse_kinet_species(out)
        if len(species) < 2:
            lbl_species.config(text="物种校验：SPECIES 段解析少于 2 个物种，请检查格式")
            set_light(3, False)
        else:
            lbl_species.config(text=f"物种校验：{len(species)} 个物种 — " + " ".join(species[:8]) + (" …" if len(species) > 8 else ""))
            refresh_comp_species()
        refresh_wizard_status()
        return out

    ttk.Button(frm3_btn, text="从模板新建", command=k_load_template).pack(side="left", padx=4)
    ttk.Button(frm3_btn, text="打开现有 inp…", command=k_load_existing).pack(side="left", padx=4)
    ttk.Button(frm3_btn, text="保存到工作目录（kinet.inp）", command=k_save).pack(side="left", padx=4)
    kinet_text.pack(fill="both", expand=True, padx=8, pady=4)
    lbl_src3.pack(anchor="w", padx=8)
    lbl_species.pack(anchor="w", padx=8)
    frm3_light = ttk.Frame(tab3); frm3_light.pack(fill="x")
    make_light(frm3_light, 3)

    # ========================================================
    # 第 4 步：主程序生成
    # ========================================================
    frm4_type = ttk.LabelFrame(tab4, text="主程序类型")
    frm4_type.pack(fill="x", padx=8, pady=4)
    var_main_kind = tk.StringVar(value="cw")
    ttk.Radiobutton(frm4_type, text="恒定 E/N 单点版（main_auto，CLI: Tgas EN t_end tag）",
                    variable=var_main_kind, value="cw").pack(side="left", padx=6)
    ttk.Radiobutton(frm4_type, text="脉冲版（main_auto_pulse，CLI: Tgas EN_on EN_off freq duty cycles tag）",
                    variable=var_main_kind, value="pulse").pack(side="left", padx=6)

    frm_comp = ttk.LabelFrame(tab4, text="初始气体组成（物种取自第 3 步机理的 SPECIES 段；电子由模板单独设置，勿选 E）")
    frm_comp.pack(fill="x", padx=8, pady=4)
    comp_rows = []  # (frame, combo, entry)
    lbl_sum = ttk.Label(frm_comp, text="当前总和: 0", foreground="#37474F")

    def species_choices():
        d = eff_dir()
        for name in ("kinet.inp", "kinetics.inp"):
            p = Path(d) / name if d else None
            if p and p.is_file():
                return parse_kinet_species(p)
        return []

    def comp_update_sum(_evt=None):
        total = 0.0
        for _f, _cb, e in comp_rows:
            try:
                total += float(e.get().strip() or 0)
            except ValueError:
                pass
        lbl_sum.config(text=f"当前总和: {total:.6g}")

    def add_comp_row(sp="", frac=""):
        f = ttk.Frame(frm_comp)
        cb = ttk.Combobox(f, width=16, state="readonly", values=species_choices())
        e = ttk.Entry(f, width=10)
        cb.pack(side="left", padx=2, pady=1)
        e.pack(side="left", padx=2)
        if sp:
            cb.set(sp)
        if frac:
            e.insert(0, frac)
        e.bind("<KeyRelease>", comp_update_sum)
        cb.bind("<<ComboboxSelected>>", comp_update_sum)

        def remove():
            comp_rows.remove((f, cb, e))
            f.destroy()
            comp_update_sum()

        ttk.Button(f, text="✕", width=3, command=remove).pack(side="left", padx=2)
        f.pack(anchor="w", pady=1)
        comp_rows.append((f, cb, e))
        comp_update_sum()
        return cb, e

    def refresh_comp_species():
        vals = species_choices()
        for _f, cb, _e in comp_rows:
            cb.config(values=vals)
        if not comp_rows and vals:
            defaults = default_composition_for_species(vals)
            if defaults:
                for sp, frac in defaults:
                    add_comp_row(sp, frac)
            else:
                add_comp_row()

    def set_composition_rows(rows):
        """把已校验的组成写回第 4 步表格（第 6 步修改后复用）。"""
        for f, _cb, _e in list(comp_rows):
            f.destroy()
        comp_rows.clear()
        for sp, frac in rows:
            add_comp_row(sp, f"{float(frac):g}")
        state["composition_dirty"] = True

    frm4_btns = ttk.Frame(frm_comp); frm4_btns.pack(anchor="w", pady=2)
    ttk.Button(frm4_btns, text="添加组分", command=lambda: add_comp_row()).pack(side="left", padx=4)
    ttk.Button(frm4_btns, text="从机理刷新物种列表", command=refresh_comp_species).pack(side="left", padx=4)
    lbl_sum.pack(anchor="w", padx=4)

    frm4_def = ttk.LabelFrame(tab4, text="默认物理条件（运行时可经命令行覆盖）")
    frm4_def.pack(fill="x", padx=8, pady=4)
    def_entries = {}
    for i, (key, label, dv) in enumerate([("en", "E/N (Td)", "45.1"), ("tg", "气体温度 (K)", "300"),
                                          ("tend", "模拟时长 (s)", "1"), ("ne", "电子密度 (cm-3)", "1.17e8")]):
        ttk.Label(frm4_def, text=label).grid(row=0, column=i * 2, sticky="w", padx=4)
        e = ttk.Entry(frm4_def, width=10)
        e.insert(0, dv)
        e.grid(row=0, column=i * 2 + 1, sticky="w", padx=4)
        def_entries[key] = e
    for i, (key, label, dv) in enumerate([("en_off", "EN_off (Td)", "0.1"), ("freq", "频率 (Hz)", "1000"),
                                          ("duty", "占空比", "0.5"), ("cycles", "最大周期数", "60")]):
        ttk.Label(frm4_def, text=label).grid(row=1, column=i * 2, sticky="w", padx=4)
        e = ttk.Entry(frm4_def, width=10)
        e.insert(0, dv)
        e.grid(row=1, column=i * 2 + 1, sticky="w", padx=4)
        def_entries[key] = e
    ttk.Label(frm4_def, text="第 2 行仅脉冲版使用（模拟时长上限=最大周期数/频率；连续 3 周期稳态后自动结束）；恒定版只用第 1 行",
              foreground="#37474F").grid(row=2, column=0, columnspan=8, sticky="w", padx=4)

    def sync_main_kind(*_):
        """脉冲参数表单随主程序类型单选启用/禁用。"""
        st = "normal" if var_main_kind.get() == "pulse" else "disabled"
        for k in ("en_off", "freq", "duty", "cycles"):
            def_entries[k].config(state=st)
    var_main_kind.trace_add("write", sync_main_kind)
    sync_main_kind()

    preview = scrolledtext.ScrolledText(tab4, height=14, state="disabled", font=("Consolas", 9))
    preview.pack(fill="both", expand=True, padx=8, pady=4)

    def gen_main(silent=False):
        d = eff_dir()
        if not d:
            if not silent:
                messagebox.showerror("错误", "请先在第 1 步选择算例目录", parent=app)
            return None
        rows = [(cb.get(), e.get()) for _f, cb, e in comp_rows]
        try:
            comp = validate_composition(rows)
            kin = Path(d) / "kinet.inp"
            species = parse_kinet_species(kin) if kin.is_file() else None
            output_opts = {"species": out_vars["out_species"].get(),
                           "te": out_vars["out_te"].get(),
                           "rates_t": out_vars["out_rates_t"].get(),
                           "rates": out_vars["out_rates"].get()}
            if var_main_kind.get() == "pulse":
                src = render_main_auto_pulse(comp, en_on=def_entries["en"].get(),
                                             en_off=def_entries["en_off"].get(),
                                             freq=def_entries["freq"].get(),
                                             duty=def_entries["duty"].get(),
                                             cycles=def_entries["cycles"].get(),
                                             tg=def_entries["tg"].get(),
                                             ne=def_entries["ne"].get(),
                                             species=species or None, output_opts=output_opts)
                main_name = "main_auto_pulse.F90"
            else:
                src = render_main_auto(comp, def_entries["en"].get(), def_entries["tg"].get(),
                                       def_entries["tend"].get(), def_entries["ne"].get(),
                                       species=species or None, output_opts=output_opts)
                main_name = "main_auto.F90"
        except (ValueError, RuntimeError) as exc:
            if not silent:
                messagebox.showerror("生成失败", str(exc), parent=app)
            set_light(4, False)
            return None
        Path(d).mkdir(parents=True, exist_ok=True)
        (Path(d) / main_name).write_text(src, encoding="utf-8")
        preview.config(state="normal")
        preview.delete("1.0", "end")
        preview.insert("1.0", "\n".join(src.splitlines()[:40]))
        preview.config(state="disabled")
        state["composition_dirty"] = False
        refresh_wizard_status()
        if not silent:
            messagebox.showinfo("完成", f"{main_name} 已生成（预览为前 40 行）", parent=app)
        return src

    frm4_gen = ttk.Frame(tab4); frm4_gen.pack(fill="x")
    ttk.Button(frm4_gen, text="生成主程序（按上方所选类型）", command=gen_main).pack(side="left", padx=8, pady=4)
    make_light(frm4_gen, 4)
    ttk.Label(tab4, text="说明：生成内容受第 6 步「E · 输出控制」中 物种浓度/电子温度/反应速率 开关影响（真实取舍 write 语句）。",
              foreground="#37474F", justify="left").pack(anchor="w", padx=8, pady=2)

    # ========================================================
    # 第 5 步：编译（门禁：前 4 步全绿）
    # ========================================================
    frm5 = ttk.Frame(tab5)
    frm5.pack(fill="both", expand=True)
    lbl_gate = ttk.Label(frm5, text="", foreground=COLOR_WARN)
    lbl_gate.pack(anchor="w", padx=8, pady=4)
    txt_log5 = scrolledtext.ScrolledText(frm5, state="disabled", font=("Consolas", 9))
    txt_log5.pack(fill="both", expand=True, padx=8, pady=4)

    def log5(line):
        emit_log("log5", line)

    def start_build5():
        if state["running"]:
            return
        d = eff_dir()
        st = eval_step_status(d)
        if not wizard_gate(st):
            missing = [f"第 {s} 步" for s in (1, 2, 3, 4) if not st[s]]
            messagebox.showerror("编译门禁", "请先完成：" + "、".join(missing), parent=app)
            return
        state["running"] = True
        state["t0"] = time.time()
        state["sim_progress"] = None
        stop_holder["killed"] = False
        btn5_build.config(state="disabled")
        btn5_stop.config(state="normal")
        params = {"en": "45.1", "tg": "300", "n2frac": "0.3333", "tend": "1",
                  "tag": "build", "pulse_enabled": False, "recompile": False}

        def worker():
            rc = execute_pipeline(d, params, do_build=True, do_run=False,
                                  log_cb=log5, stop_holder=stop_holder)
            q.put(("done5", rc))

        threading.Thread(target=worker, daemon=True).start()

    frm5_btn = ttk.Frame(frm5); frm5_btn.pack(fill="x", side="bottom")
    btn5_build = ttk.Button(frm5_btn, text="部署脚本并编译", command=start_build5)
    btn5_build.pack(side="left", padx=8, pady=4)
    btn5_stop = ttk.Button(frm5_btn, text="中止", command=lambda: kill_process_tree(stop_holder), state="disabled")
    btn5_stop.pack(side="left", padx=4)
    make_light(frm5_btn, 5)

    # ========================================================
    # 第 6 步：运行与扫描（需求 3：参数五模块；需求 5：关键输入样式）
    # ========================================================
    frm_dir = ttk.LabelFrame(tab6, text="目录检测（工作目录由第 1 步 算例目录+批次名 决定）")
    frm_dir.pack(fill="x", padx=8, pady=4)
    txt_check = scrolledtext.ScrolledText(frm_dir, height=5, state="disabled", font=("Consolas", 9))
    txt_check.pack(fill="x", padx=4, pady=2)
    lbl_flavor = ttk.Label(frm_dir, text="主程序类型：未检测", foreground="#37474F")
    lbl_flavor.pack(anchor="w", padx=4)

    btn_state = {}
    entries = {}

    def refresh_detection(_evt=None):
        d = eff_dir()
        txt_check.config(state="normal"); txt_check.delete("1.0", "end")
        if not d or not Path(d).is_dir():
            txt_check.insert("end", "请先在第 1 步选择有效目录\n")
            txt_check.config(state="disabled")
            set_run_buttons(False)
            if "scan" in btn_state:
                btn_state["scan"].config(state="disabled")
            return
        checks, can_build = detect_inputs(d)
        for label, ok, detail in checks:
            mark = "✓" if ok else "✗"
            line = f" {mark}  {label}" + (f"    — {detail}" if detail else "") + "\n"
            txt_check.insert("end", line)
        flavor = detect_flavor(d)
        state["flavor"] = flavor
        lbl_flavor.config(text=f"主程序类型：{FLAVOR_NAMES.get(flavor, '未识别')}（工作目录: {d}）")
        txt_check.insert("end", "\n" + ("可以编译 ✓" if can_build else "缺关键文件，编译按钮已禁用（补齐上方 ✗ 项）") + "\n")
        txt_check.config(state="disabled")
        set_run_buttons(bool(can_build))
        if "scan" in btn_state:
            btn_state["scan"].config(state="normal" if flavor else "disabled")
        update_param_state()

    # ---------- 模块 A · 物理场与气体（关键输入，置顶最醒目） ----------
    frm_a = ttk.LabelFrame(tab6, text="A · 物理场与气体（关键输入）")
    frm_a.pack(fill="x", padx=8, pady=4)
    tk.Label(frm_a, text="E/N (Td)", font=FONT_KEY, fg=COLOR_KEY).grid(row=0, column=0, sticky="w", padx=4, pady=2)
    e_en = tk.Entry(frm_a, width=12, font=FONT_KEY, fg=COLOR_KEY)
    e_en.insert(0, "45.1")
    e_en.grid(row=0, column=1, sticky="w", padx=4, pady=2)
    entries["en"] = e_en
    tk.Label(frm_a, text="气体温度 Tgas (K)", font=FONT_KEY, fg=COLOR_KEY).grid(row=0, column=2, sticky="w", padx=4, pady=2)
    e_tg = tk.Entry(frm_a, width=12, font=FONT_KEY, fg=COLOR_KEY)
    e_tg.insert(0, "300")
    e_tg.grid(row=0, column=3, sticky="w", padx=4, pady=2)
    entries["tg"] = e_tg
    lbl_warn_a = tk.Label(frm_a, text="", fg=COLOR_WARN, justify="left")
    lbl_warn_a.grid(row=1, column=0, columnspan=4, sticky="w", padx=4)
    ttk.Label(frm_a, text="命令行型 / 脉冲型 / 向导生成型走命令行直传；参数固化型见上方警告。"
                          "默认值 45.1 Td / 300 K 为常用默认，可按算例修改。",
              foreground="#9E9E9E").grid(row=2, column=0, columnspan=4, sticky="w", padx=4)

    # ---------- 模块 B · DBD 放电参数 ----------
    frm_b = ttk.LabelFrame(tab6, text="B · DBD 放电参数")
    frm_b.pack(fill="x", padx=8, pady=4)
    for i, (key, label, dv) in enumerate([("n2frac", "N2 摩尔分数（常用默认 0.3333）", "0.3333"),
                                          ("tend", "模拟时长 (s)", "1"),
                                          ("tag", "批次名 tag（输出 runs/<tag>/）", "gui_run")]):
        ttk.Label(frm_b, text=label).grid(row=0, column=i * 2, sticky="w", padx=4, pady=2)
        e = ttk.Entry(frm_b, width=16)
        e.insert(0, dv)
        e.grid(row=0, column=i * 2 + 1, sticky="w", padx=4, pady=2)
        entries[key] = e
    entries["en_list"] = ttk.Entry(frm_b)  # 隐藏占位：兼容 validate_params（扫描走下面板）
    entries["atol"] = ttk.Entry(frm_b)     # 占位，模块 D 重新布局真正的输入框
    entries["rtol"] = ttk.Entry(frm_b)
    # v3 需求 3：auto（向导生成型）组分固化在源码中 → 只读摘要 + 回第 4 步跳转
    frm_b_auto = ttk.Frame(frm_b)
    lbl_b_auto = tk.Label(frm_b_auto, text="", fg=COLOR_MANUAL, justify="left")
    lbl_b_auto.pack(side="left", padx=4)
    def edit_runtime_composition():
        """第 6 步编辑向导生成型的初始组成，并标记为需自动重走 4→5→6。"""
        current = [(cb.get().strip(), e.get().strip()) for _f, cb, e in comp_rows
                   if cb.get().strip() or e.get().strip()]
        if not current:
            d = eff_dir()
            src = _find_main_src(Path(d)) if d and Path(d).is_dir() else None
            current = parse_auto_composition(src) if src else []
        if not current:
            current = default_composition_for_species(species_choices())
        win = tk.Toplevel(app)
        win.title("修改初始气体组成 — 修改后自动重新生成并编译")
        win.transient(app); win.grab_set()
        ttk.Label(win, text="每行：物种  摩尔分数；总和必须为 1。\n保存后，下一次运行/扫描将自动执行第 4 步生成 → 第 5 步编译 → 第 6 步计算。",
                  justify="left").pack(anchor="w", padx=12, pady=(12, 4))
        editor = scrolledtext.ScrolledText(win, width=42, height=8, font=("Consolas", 10))
        editor.pack(fill="both", expand=True, padx=12, pady=4)
        editor.insert("1.0", "\n".join(f"{sp} {frac}" for sp, frac in current))

        def apply_composition():
            rows = []
            for no, line in enumerate(editor.get("1.0", "end-1c").splitlines(), 1):
                line = line.strip()
                if not line:
                    continue
                parts = re.split(r"[\s,;]+", line)
                if len(parts) != 2:
                    messagebox.showerror("格式错误", f"第 {no} 行应为：物种 空格 摩尔分数", parent=win)
                    return
                rows.append((parts[0], parts[1]))
            try:
                comp = validate_composition(rows)
            except ValueError as exc:
                messagebox.showerror("组成不合法", str(exc), parent=win)
                return
            known = set(species_choices())
            unknown = [sp for sp, _frac in comp if known and sp not in known]
            if unknown:
                messagebox.showerror("物种不在机理中", "未在第 3 步机理识别到：" + ", ".join(unknown), parent=win)
                return
            set_composition_rows(comp)
            update_module_b_state()
            win.destroy()

        btns = ttk.Frame(win); btns.pack(fill="x", padx=12, pady=(4, 12))
        ttk.Button(btns, text="保存组成", command=apply_composition).pack(side="left")
        ttk.Button(btns, text="取消", command=win.destroy).pack(side="left", padx=6)

    btn_back4 = ttk.Button(frm_b_auto, text="修改初始组成…", command=edit_runtime_composition)
    btn_back4.pack(side="left", padx=8)
    frm_b_auto.grid(row=1, column=0, columnspan=6, sticky="w", padx=4, pady=2)
    frm_b_auto.grid_remove()  # 默认隐藏，检测到 auto 型才显示

    # ---------- 模块 C · 脉冲参数（勾选展开；不勾选 = --cw 连续场） ----------
    var_pulse = tk.BooleanVar(value=False)
    frm_c = ttk.LabelFrame(tab6, text="C · 脉冲参数（仅 pulse 型主程序使用）")
    frm_c.pack(fill="x", padx=8, pady=4)
    frm_pulse = ttk.Frame(frm_c)
    pulse_grid = [("en_off", "EN_off (Td)", "0.1"), ("freq", "频率 (Hz)", "1000"),
                  ("duty", "占空比", "0.5"), ("cycles", "最大周期数", "60"),
                  ("tau_res", "停留时间 (s；0=封闭)", "0"), ("ne_mode", "ne_mode", "fix")]
    for i, (key, label, default) in enumerate(pulse_grid):
        ttk.Label(frm_pulse, text=label).grid(row=0, column=i * 2, sticky="w", padx=4, pady=2)
        e = ttk.Entry(frm_pulse, width=10)
        e.insert(0, default)
        e.grid(row=0, column=i * 2 + 1, sticky="w", padx=4, pady=2)
        entries[key] = e
    lbl_ne_mode_note = ttk.Label(frm_pulse, text="", foreground="#37474F")
    lbl_ne_mode_note.grid(row=1, column=10, columnspan=2, sticky="w", padx=4)
    ttk.Label(frm_pulse, text="收敛：显著物种、NH3 周期增量和周期能量均 ≤0.1%，连续 3 周期。",
              foreground="#37474F").grid(row=1, column=0, columnspan=10, sticky="w", padx=4)

    def toggle_pulse():
        if var_pulse.get():
            frm_pulse.grid()
        else:
            frm_pulse.grid_remove()

    chk_pulse = ttk.Checkbutton(frm_c, text="启用脉冲模式（不勾选 = 连续恒定场 --cw）",
                                variable=var_pulse, command=toggle_pulse)
    chk_pulse.grid(row=0, column=0, sticky="w", padx=4, pady=2)
    frm_pulse.grid(row=1, column=0, sticky="ew", padx=4, pady=2)
    frm_pulse.grid_remove()

    # ---------- 模块 D · 收敛控制 ----------
    frm_d = ttk.LabelFrame(tab6, text="D · 收敛控制")
    frm_d.pack(fill="x", padx=8, pady=4)
    ttk.Label(frm_d, text="atol").grid(row=0, column=0, sticky="w", padx=4, pady=2)
    e_atol = ttk.Entry(frm_d, width=12)
    e_atol.insert(0, DEFAULT_ATOL)
    e_atol.grid(row=0, column=1, sticky="w", padx=4, pady=2)
    ttk.Label(frm_d, text="rtol").grid(row=0, column=2, sticky="w", padx=4, pady=2)
    e_rtol = ttk.Entry(frm_d, width=12)
    e_rtol.insert(0, DEFAULT_RTOL)
    e_rtol.grid(row=0, column=3, sticky="w", padx=4, pady=2)
    entries["atol"] = e_atol
    entries["rtol"] = e_rtol
    ttk.Label(frm_d, text="默认 atol=1.0、rtol=1e-4；需主程序支持位置参数（命令行型/脉冲型/向导生成型），参数固化型不可用。",
               foreground="#9E9E9E").grid(row=1, column=0, columnspan=4, sticky="w", padx=4)
    var_retry = tk.BooleanVar(value=False)
    chk_retry = ttk.Checkbutton(frm_d, text="失败后自动使用宽松容差重试一次", variable=var_retry)
    chk_retry.grid(row=2, column=0, columnspan=2, sticky="w", padx=4, pady=(4, 2))
    ttk.Label(frm_d, text="二次 atol").grid(row=2, column=2, sticky="w", padx=4, pady=2)
    e_retry_atol = ttk.Entry(frm_d, width=12)
    e_retry_atol.insert(0, DEFAULT_RETRY_ATOL)
    e_retry_atol.grid(row=2, column=3, sticky="w", padx=4, pady=2)
    ttk.Label(frm_d, text="二次 rtol").grid(row=2, column=4, sticky="w", padx=4, pady=2)
    e_retry_rtol = ttk.Entry(frm_d, width=12)
    e_retry_rtol.insert(0, DEFAULT_RETRY_RTOL)
    e_retry_rtol.grid(row=2, column=5, sticky="w", padx=4, pady=2)
    entries["retry_atol"] = e_retry_atol
    entries["retry_rtol"] = e_retry_rtol

    def toggle_retry():
        st = "normal" if var_retry.get() else "disabled"
        e_retry_atol.config(state=st)
        e_retry_rtol.config(state=st)
    chk_retry.config(command=toggle_retry)
    toggle_retry()

    # ---------- 模块 E · 输出控制 ----------
    frm_e = ttk.LabelFrame(tab6, text="E · 输出控制")
    frm_e.pack(fill="x", padx=8, pady=4)
    out_defs = [("out_species", "物种浓度时程"), ("out_te", "电子温度时程"),
                ("out_rates", "反应速率（末态）"), ("out_rates_t", "反应速率时程"),
                ("out_energy", "每周期能耗/能效（pulse）"), ("out_console", "保存控制台日志"),
                ("out_summary", "批次汇总 CSV")]
    for i, (key, label) in enumerate(out_defs):
        ttk.Checkbutton(frm_e, text=label, variable=out_vars[key]).grid(
            row=i // 4, column=i % 4, sticky="w", padx=6, pady=1)
    ttk.Label(frm_e, text="对向导生成的 main_auto：第 4 步生成时真实取舍 write 语句；对现成主程序：仅控制下方产物清单的展示筛选，不改主程序源码。",
              foreground=COLOR_WARN, justify="left").grid(row=2, column=0, columnspan=4, sticky="w", padx=4)

    var_recompile = tk.BooleanVar(value=False)
    chk_recompile = ttk.Checkbutton(tab6, text="--recompile（参数固化型修改硬编码 E/N/温度时勾选，将备份并补丁源码）",
                                    variable=var_recompile)
    chk_recompile.pack(anchor="w", padx=8)

    # ---------- 扫描面板 ----------
    frm_scan = ttk.LabelFrame(tab6, text="多维参数扫描（可选；启用后用“开始扫描”运行，单次运行不受影响）")
    frm_scan.pack(fill="x", padx=8, pady=4)
    var_scan = tk.BooleanVar(value=False)
    var_ndim = tk.IntVar(value=1)
    scan_axis_combos, scan_axis_entries = [], []

    chk_scan = ttk.Checkbutton(frm_scan, text="启用参数扫描", variable=var_scan)
    chk_scan.grid(row=0, column=0, sticky="w", padx=4, pady=2)
    ttk.Label(frm_scan, text="维度:").grid(row=0, column=1, sticky="e")
    cmb_ndim = ttk.Combobox(frm_scan, width=3, state="readonly", values=(1, 2, 3),
                            textvariable=var_ndim)
    cmb_ndim.grid(row=0, column=2, sticky="w", padx=2)
    lbl_points = ttk.Label(frm_scan, text="预计点数: —", foreground="#37474F")
    lbl_points.grid(row=0, column=3, sticky="w", padx=8)
    lbl_scan_hint = ttk.Label(frm_scan, text="", foreground=COLOR_WARN)
    lbl_scan_hint.grid(row=0, column=4, sticky="w", padx=8)
    ttk.Label(frm_scan, text="并行任务数:").grid(row=0, column=5, sticky="e", padx=(12, 2))
    e_workers = ttk.Entry(frm_scan, width=4)
    e_workers.insert(0, str(default_scan_workers()))
    e_workers.grid(row=0, column=6, sticky="w")
    entries["scan_workers"] = e_workers
    ttk.Label(frm_scan, text="默认保留 1 个 CPU 核，最多 4 个；设为 1 即串行", foreground="#9E9E9E").grid(
        row=0, column=7, sticky="w", padx=4)

    for r in range(SCAN_MAX_DIM):
        ttk.Label(frm_scan, text=f"轴{r + 1}:").grid(row=r + 1, column=0, sticky="e", padx=4)
        cb = ttk.Combobox(frm_scan, width=14, state="readonly")
        cb.grid(row=r + 1, column=1, sticky="w", padx=2, pady=2)
        e = tk.Entry(frm_scan, width=26, font=FONT_KEY, fg=COLOR_KEY)  # 关键输入：扫描轴取值
        e.grid(row=r + 1, column=2, columnspan=2, sticky="w", padx=2, pady=2)
        ttk.Label(frm_scan, text="取值: 列表（30 60 90）或 起:止:步数（30:120:5）").grid(
            row=r + 1, column=4, sticky="w", padx=4)
        scan_axis_combos.append(cb)
        scan_axis_entries.append(e)

    def _scan_axis_keys():
        nd = var_ndim.get()
        keys = []
        for r in range(nd):
            label = scan_axis_combos[r].get()
            key = next((k for k, lab in _scan_axis_map if lab == label), None)
            keys.append(key)
        return keys

    def read_scan_axes():
        nd = var_ndim.get()
        keys = _scan_axis_keys()
        axes = []
        for r in range(nd):
            if keys[r] is None:
                raise ValueError(f"轴{r + 1} 未选择参数")
            axes.append((keys[r], parse_axis_values(scan_axis_entries[r].get())))
        return axes

    def update_points_label(_evt=None):
        if not var_scan.get():
            lbl_points.config(text="预计点数: —")
            return
        try:
            axes = read_scan_axes()
            n = 1
            for _k, vals in axes:
                n *= len(vals)
            lbl_points.config(text=f"预计点数: {n}")
        except ValueError as exc:
            lbl_points.config(text=f"预计点数: ?（{exc}）")

    def refresh_scan_panel(_evt=None):
        flavor = state.get("flavor")
        allowed = scan_axes_for_flavor(flavor)
        if flavor in ("pulse", "autopulse") and not var_pulse.get():
            pulse_only = {k for k, _ in SCAN_AXES_PULSE}
            allowed = [(k, lab) for k, lab in allowed if k not in pulse_only]
        labels = [lab for _k, lab in allowed]
        nonlocal _scan_axis_map
        _scan_axis_map = allowed
        for r, cb in enumerate(scan_axis_combos):
            cb.config(values=labels)
            if labels and cb.get() not in labels:
                cb.current(min(r, len(labels) - 1))
        if flavor == "hong":
            # v3 需求 4：固化型可扫 N2 摩尔分数 / 模拟时长（命令行位置参数直传）
            chk_scan.config(state="normal")
            lbl_scan_hint.config(text="参数固化型：E/N、温度硬编码不可扫；可扫 N2 摩尔分数 / 模拟时长")
        elif flavor is None:
            chk_scan.config(state="disabled")
            lbl_scan_hint.config(text="")
        elif flavor == "auto":
            chk_scan.config(state="normal")
            lbl_scan_hint.config(text="向导生成型：组分已固化，N2 比例不可扫（回第 4 步重新生成可改）")
        elif flavor == "autopulse":
            chk_scan.config(state="normal")
            lbl_scan_hint.config(text="脉冲·向导生成型：组分已固化、电子恒定、时长=周期数/频率"
                                 if var_pulse.get() else "未勾选脉冲：按 --cw 连续场逐点运行")
        else:
            chk_scan.config(state="normal")
            lbl_scan_hint.config(text="" if flavor != "pulse" or var_pulse.get()
                                 else "未勾选脉冲：按 --cw 连续场逐点运行")
        nd = var_ndim.get()
        for r in range(SCAN_MAX_DIM):
            st = "normal" if r < nd else "disabled"
            scan_axis_combos[r].config(state="readonly" if r < nd else "disabled")
            scan_axis_entries[r].config(state=st)
        update_points_label()

    _scan_axis_map = scan_axes_for_flavor(None)
    chk_scan.config(command=lambda: (refresh_scan_panel(), update_points_label()))
    var_ndim.trace_add("write", lambda *_: refresh_scan_panel())
    for cb in scan_axis_combos:
        cb.bind("<<ComboboxSelected>>", update_points_label)
    for e in scan_axis_entries:
        e.bind("<KeyRelease>", update_points_label)
    var_pulse.trace_add("write", lambda *_: refresh_scan_panel())

    def update_module_a_warn():
        d = eff_dir()
        src = _find_main_src(Path(d)) if d and Path(d).is_dir() else None
        msg = module_a_warning(state.get("flavor"), src,
                               entries["en"].get(), entries["tg"].get())
        lbl_warn_a.config(text=msg or "")

    def update_module_b_state():
        """向导生成型在第 6 步显示可编辑组成；改动标记为自动重走 4→5→6。"""
        flavor = state.get("flavor")
        if flavor in ("auto", "autopulse"):
            src = None
            d = eff_dir()
            if d and Path(d).is_dir():
                src = _find_main_src(Path(d))
            comp = [(cb.get(), e.get()) for _f, cb, e in comp_rows] if state["composition_dirty"] else []
            comp = comp or (parse_auto_composition(src) if src else [])
            comp_txt = fmt_composition(comp) or "见主程序源码"
            entries["n2frac"].config(state="disabled")
            suffix = "已在第 6 步修改；下次运行将自动重走第 4→5→6 步" if state["composition_dirty"] else "可在此修改"
            lbl_b_auto.config(text=f"气体组分：{comp_txt} —— {suffix}")
            frm_b_auto.grid()
        else:
            entries["n2frac"].config(state="normal")
            frm_b_auto.grid_remove()

    def update_param_state():
        flavor = state.get("flavor")
        var_pulse.set(flavor in ("pulse", "autopulse"))
        # autopulse（脉冲·向导生成型）电子密度固定恒定，ne_mode 无意义：禁用并注明
        if flavor == "autopulse":
            entries["ne_mode"].config(state="disabled")
            lbl_ne_mode_note.config(text="生成版固定恒定电子密度")
        else:
            entries["ne_mode"].config(state="normal")
            lbl_ne_mode_note.config(text="")
        toggle_pulse()
        update_module_a_warn()
        update_module_b_state()
        refresh_scan_panel()

    entries["en"].bind("<KeyRelease>", lambda *_: update_module_a_warn())
    entries["tg"].bind("<KeyRelease>", lambda *_: update_module_a_warn())

    # ---------- 执行按钮区 ----------
    frm_btn = ttk.Frame(tab6)
    frm_btn.pack(fill="x", padx=8, pady=4)
    prog = ttk.Progressbar(frm_btn, mode="determinate", maximum=100)
    lbl_time = ttk.Label(frm_btn, text="就绪", width=30)

    def log(line):
        emit_log("log", line)

    def collect_params():
        if var_batch.get().strip() and not re.match(r"^[A-Za-z0-9_\-]+$", var_batch.get().strip()):
            raise ValueError(f"批次名只允许字母/数字/下划线/连字符: {var_batch.get()!r}")
        p = {k: e.get().strip() for k, e in entries.items()}
        p["recompile"] = var_recompile.get()
        p["pulse_enabled"] = var_pulse.get()
        p["retry_on_failure"] = var_retry.get()
        p["ne_mode"] = entries["ne_mode"].get().strip() or "fix"
        for k, v in out_vars.items():
            p[k] = v.get()
        if var_batch.get().strip():
            p["tag"] = var_batch.get().strip()  # 需求 4：批次模式下批次名即 tag
        return validate_params(p)

    def set_run_buttons(enabled):
        for name in ("build", "run", "both"):
            if name in btn_state:
                btn_state[name].config(state=("normal" if enabled else "disabled"))
        if "stop" in btn_state:
            btn_state["stop"].config(state=("normal" if not enabled and state["running"] else "disabled"))

    def start_task(do_build, do_run):
        if state["running"]:
            return
        d = eff_dir()
        if not d or not Path(d).is_dir():
            messagebox.showerror("错误", "请先在第 1 步选择有效的算例目录", parent=app)
            return
        try:
            params = collect_params()
        except ValueError as exc:
            messagebox.showerror("参数不合法", str(exc), parent=app)
            return
        # 第 6 步修改向导生成型组成后，不让旧 exe 静默继续运行：自动执行 4→5→6。
        regenerated = state.get("composition_dirty", False)
        if regenerated:
            if gen_main(silent=True) is None:
                messagebox.showerror("生成失败", "第 6 步组分修改后无法重新生成主程序", parent=app)
                return
            do_build = True
            log("[信息] 检测到第 6 步组分修改：自动重新生成主程序并重新编译")
        state["running"] = True
        state["t0"] = time.time()
        state["sim_progress"] = None
        stop_holder["killed"] = False
        try:
            t_end_hint["value"] = float(str(params["tend"]).replace("d", "e").replace("D", "E"))
        except (ValueError, KeyError):
            t_end_hint["value"] = None
        for name in ("build", "run", "both", "scan"):
            btn_state[name].config(state="disabled")
        btn_state["stop"].config(state="normal")
        prog.config(mode="indeterminate"); prog.start(12)
        flat = bool(var_batch.get().strip())

        def worker():
            rc = execute_pipeline(d, params, do_build=do_build, do_run=do_run,
                                  log_cb=log, stop_holder=stop_holder, flat_out=flat)
            q.put(("done", rc, params["tag"]))

        threading.Thread(target=worker, daemon=True).start()

    def start_scan():
        if state["running"]:
            return
        d = eff_dir()
        if not d or not Path(d).is_dir():
            messagebox.showerror("错误", "请先在第 1 步选择有效的算例目录", parent=app)
            return
        flavor = state.get("flavor")
        if not flavor:
            messagebox.showerror("错误", "未识别主程序类型，无法扫描", parent=app)
            return
        # v3 需求 4：hong 型允许 n2frac/tend 轴（execute_scan 内置二次守卫，非法轴报错弹窗）
        try:
            params = collect_params()
            axes = read_scan_axes()
        except ValueError as exc:
            messagebox.showerror("扫描参数不合法", str(exc), parent=app)
            return
        n = 1
        for _k, vals in axes:
            n *= len(vals)
        if n > SCAN_CONFIRM_POINTS and not messagebox.askyesno(
                "点数较多", f"扫描共 {n} 个点，可能耗时很长。确认继续？", parent=app):
            return
        regenerated = state.get("composition_dirty", False)
        if regenerated:
            if gen_main(silent=True) is None:
                messagebox.showerror("生成失败", "第 6 步组分修改后无法重新生成主程序", parent=app)
                return
        state["running"] = True
        state["t0"] = time.time()
        state["sim_progress"] = None
        stop_holder["killed"] = False
        for name in ("build", "run", "both", "scan"):
            btn_state[name].config(state="disabled")
        btn_state["stop"].config(state="normal")
        prog.stop()
        prog.config(mode="determinate", maximum=100, value=0)
        flat = bool(var_batch.get().strip())

        def worker():
            def progress(i, total, _row):
                q.put(("scan_progress", i, total))
            def scan_status(kind, i, total, tag, _params):
                q.put(("scan_status", kind, i, total, tag))
            try:
                if regenerated:
                    emit_log("log", "[信息] 第 6 步组分已变更：扫描前自动重新编译")
                    rc = execute_pipeline(d, params, do_build=True, do_run=False,
                                          log_cb=log, stop_holder=stop_holder, flat_out=flat)
                    if rc != 0:
                        q.put(("scan_error", f"自动重新编译失败（退出码 {rc}），扫描未启动"))
                        return
                rows, summary, stopped = execute_scan(
                    d, params, axes, flavor, log_cb=log,
                    stop_holder=stop_holder, progress_cb=progress, flat_out=flat,
                    max_workers=int(params["scan_workers"]), status_cb=scan_status)
                q.put(("scan_done", rows, str(summary), stopped,
                       [k for k, _ in axes], metric_names_for_flavor(flavor)))
            except Exception as exc:  # noqa: BLE001
                q.put(("scan_error", str(exc)))

        threading.Thread(target=worker, daemon=True).start()

    btn_state["build"] = ttk.Button(frm_btn, text="编译", command=lambda: start_task(True, False), state="disabled")
    btn_state["run"] = ttk.Button(frm_btn, text="运行", command=lambda: start_task(False, True), state="disabled")
    btn_state["both"] = tk.Button(frm_btn, text="一键编译并运行", font=FONT_KEY, fg=COLOR_KEY,
                                  command=lambda: start_task(True, True), state="disabled")
    btn_state["scan"] = tk.Button(frm_btn, text="开始扫描", font=FONT_KEY, fg=COLOR_KEY,
                                  command=start_scan, state="disabled")
    btn_state["stop"] = ttk.Button(frm_btn, text="中止", command=lambda: kill_process_tree(stop_holder), state="disabled")
    for name in ("build", "run", "both", "scan", "stop"):
        btn_state[name].pack(side="left", padx=4, pady=4)
    prog.pack(side="left", fill="x", expand=True, padx=8)
    lbl_time.pack(side="right", padx=4)

    # ---------- 日志区 ----------
    frm_log = ttk.LabelFrame(tab6, text="实时日志")
    frm_log.pack(fill="both", expand=True, padx=8, pady=4)
    txt_log = scrolledtext.ScrolledText(frm_log, state="disabled", font=("Consolas", 9))
    txt_log.pack(fill="both", expand=True)

    # ---------- 结果区 ----------
    frm_res = ttk.LabelFrame(tab6, text="输出文件（受模块 E 展示筛选）")
    frm_res.pack(fill="x", padx=8, pady=4)
    tree = ttk.Treeview(frm_res, columns=("size",), show="tree headings", height=5)
    tree.heading("#0", text="文件"); tree.heading("size", text="大小 (B)")
    tree.column("size", width=120, anchor="e")
    tree.pack(side="left", fill="x", expand=True, padx=4, pady=4)
    frm_res_btn = ttk.Frame(frm_res)
    frm_res_btn.pack(side="right", padx=4)

    def _product_visible(name):
        """模块 E 后处理筛选（仅影响展示，不影响已生成的文件）。"""
        n = name.lower()
        if n == "console.log":
            return out_vars["out_console"].get()
        if n == "scan_summary.csv":
            return out_vars["out_summary"].get()
        if "pulse_cycles" in n:
            return out_vars["out_energy"].get()
        if "_rates" in n:
            return out_vars["out_rates"].get() or out_vars["out_rates_t"].get()
        if ("_output_" in n) or n.startswith("pulse_series") or n == "qt_species_density.txt":
            return out_vars["out_species"].get() or out_vars["out_te"].get()
        return True  # bolsigdb.dat / 计算参数.md / params.json / 其他始终显示

    def refresh_results(tag):
        for item in tree.get_children():
            tree.delete(item)
        d = eff_dir()
        files = [(f, sz) for f, sz in list_outputs(d, tag) if _product_visible(f.name)]
        outdir = Path(d) / "runs" / tag
        if not outdir.is_dir() and (Path(d) / "scan_summary.csv").is_file():
            outdir = Path(d)  # 批次扁平模式
        state["last_outdir"] = outdir if outdir.is_dir() else (Path(d) / "runs" if (Path(d) / "runs").is_dir() else Path(d))
        if not files:
            tree.insert("", "end", text="（未发现 runs/<tag>/ 产物）", values=("",))
            return
        for f, sz in files:
            tree.insert("", "end", text=str(f.relative_to(Path(d))), values=(sz,))

    def open_outdir():
        outdir = state.get("last_outdir")
        if outdir and Path(outdir).is_dir():
            os.startfile(str(outdir))  # noqa: S606
        else:
            messagebox.showinfo("提示", "还没有可打开的输出目录", parent=app)

    def do_plot():
        outdir = state.get("last_outdir")
        if not outdir or not Path(outdir).is_dir():
            messagebox.showinfo("提示", "还没有可绘图的输出目录", parent=app)
            return
        quick_plot(app, Path(outdir))

    ttk.Button(frm_res_btn, text="打开输出文件夹", command=open_outdir).pack(pady=4)
    ttk.Button(frm_res_btn, text="快速绘图", command=do_plot).pack(pady=4)

    # ---------- 扫描汇总区 ----------
    frm_sum = ttk.LabelFrame(tab6, text="扫描汇总（每次扫描完成自动刷新）")
    frm_sum.pack(fill="x", padx=8, pady=4)
    sum_tree = ttk.Treeview(frm_sum, show="headings", height=5)
    sum_xscroll = ttk.Scrollbar(frm_sum, orient="horizontal", command=sum_tree.xview)
    sum_tree.configure(xscrollcommand=sum_xscroll.set)
    sum_tree.pack(side="top", fill="x", expand=True, padx=4, pady=(4, 0))
    sum_xscroll.pack(side="top", fill="x", padx=4)
    frm_sum_btn = ttk.Frame(frm_sum)
    frm_sum_btn.pack(side="bottom", anchor="e", padx=4, pady=2)

    def update_scan_tree(rows, axis_keys, metric_names):
        cols = list(axis_keys) + list(metric_names) + ["status", "failure_reason", "attempts", "elapsed_s", "rc", "tag"]
        sum_tree["columns"] = cols
        for c in cols:
            sum_tree.heading(c, text=c)
            sum_tree.column(c, width=(250 if c == "failure_reason" else 92), anchor="e", stretch=False)
        sum_tree.tag_configure("scan_failed", foreground=COLOR_MISS)
        sum_tree.tag_configure("scan_retry", foreground=COLOR_WARN)
        sum_tree.tag_configure("scan_unconverged", foreground="#8A6D00")
        for item in sum_tree.get_children():
            sum_tree.delete(item)
        for row in rows:
            tags = ()
            if row.get("status") in ("failed", "invalid_output"):
                tags = ("scan_failed",)
            elif row.get("status") in ("not_converged", "convergence_unknown"):
                tags = ("scan_unconverged",)
            elif row.get("status") == "retry_success":
                tags = ("scan_retry",)
            sum_tree.insert("", "end", values=[_fmt_cell(row.get(c, _NAN)) for c in cols], tags=tags)
        return cols

    def do_scan_plot():
        rows = state.get("scan_rows")
        if not rows:
            messagebox.showinfo("提示", "还没有扫描结果", parent=app)
            return
        scan_plot(app, rows, state.get("scan_axes") or [], state.get("scan_metrics") or [])

    def open_summary_dir():
        sp = state.get("scan_summary")
        if sp and Path(sp).is_file():
            os.startfile(str(Path(sp).parent))  # noqa: S606
        else:
            messagebox.showinfo("提示", "还没有扫描汇总文件", parent=app)

    ttk.Button(frm_sum_btn, text="扫描绘图", command=do_scan_plot).pack(side="left", padx=4)
    ttk.Button(frm_sum_btn, text="导出汇总CSV（打开所在文件夹）", command=open_summary_dir).pack(side="left", padx=4)

    # ---------- 向导状态总刷新 ----------
    def refresh_wizard_status(_evt=None):
        refresh_slots()
        st = eval_step_status(eff_dir())
        for s in (1, 2, 3, 4):
            set_light(s, True if st[s] else False)
        d = Path(eff_dir()) if eff_dir() else None
        has_exe = bool(d and d.is_dir() and any(d.glob("main_*.exe")))
        if not state["running"]:
            set_light(5, True if has_exe else None, "（已有 exe）" if has_exe else "")
        missing = [f"第 {s} 步" for s in (1, 2, 3, 4) if not st[s]]
        lbl_gate.config(text="" if not missing else "编译前需完成：" + "、".join(missing))
        refresh_detection()

    # ---------- 底部环境自检区（需求 6） ----------
    frm_env = ttk.Frame(app)
    frm_env.pack(fill="x", padx=8, pady=2)
    lbl_env = ttk.Label(frm_env, text="环境: 自检中…", foreground="#37474F")
    lbl_env.pack(side="left", padx=4)

    def run_env_check(popup=False):
        checks = self_check(load_config())
        parts = []
        for name, ok, _detail in checks:
            mark = "✓" if ok is True else ("⚠" if ok is None else "✗")
            parts.append(f"{name.split('（')[0]} {mark}")
        lbl_env.config(text="环境: " + " · ".join(parts))
        if popup:
            win = tk.Toplevel(app)
            win.title("环境自检 — 依赖解析结果")
            for name, ok, detail in checks:
                mark = "✓" if ok is True else ("⚠" if ok is None else "✗")
                col = COLOR_OK if ok is True else (COLOR_WARN if ok is None else COLOR_MISS)
                tk.Label(win, text=f"{mark} {name}", fg=col, font=("Segoe UI", 10, "bold"),
                         anchor="w").pack(fill="x", padx=10, pady=(6, 0))
                tk.Label(win, text=str(detail), fg="#37474F", wraplength=620,
                         justify="left", anchor="w").pack(fill="x", padx=28)
            ttk.Button(win, text="关闭", command=win.destroy).pack(pady=8)

    ttk.Button(frm_env, text="环境自检", command=lambda: run_env_check(popup=True)).pack(side="right", padx=4)

    # ---------- 底部导航 ----------
    frm_nav = ttk.Frame(app)
    frm_nav.pack(fill="x", padx=8, pady=2)

    def nav(delta):
        i = (nb.index("current") + delta) % 6
        nb.select(i)

    ttk.Button(frm_nav, text="← 上一步", command=lambda: nav(-1)).pack(side="left", padx=4)
    ttk.Button(frm_nav, text="下一步 →", command=lambda: nav(1)).pack(side="left", padx=4)
    ttk.Label(frm_nav, text="可跳步查看；编译前第 1–4 步须全部通过",
              foreground=COLOR_GRAY).pack(side="right", padx=8)

    # ---------- 队列轮询 ----------
    t_end_hint = {"value": None}

    def poll():
        try:
            while True:
                kind, *payload = q.get_nowait()
                if kind in ("log", "log5"):
                    target = txt_log if kind == "log" else txt_log5
                    line = payload[0]
                    target.config(state="normal")
                    target.insert("end", line + "\n")
                    target.see("end")
                    target.config(state="disabled")
                    if kind == "log":
                        te = t_end_hint.get("value")
                        if te:
                            frac = parse_progress(line, te)
                            if frac is not None:
                                if prog["mode"] == "indeterminate":
                                    prog.stop()
                                    prog.config(mode="determinate")
                                prog["value"] = frac * 100
                                state["sim_progress"] = (frac * te, te)
                elif kind == "done5":
                    rc = payload[0]
                    state["running"] = False
                    btn5_build.config(state="normal")
                    btn5_stop.config(state="disabled")
                    set_light(5, rc == 0)
                    refresh_wizard_status()
                    if rc == 0:
                        messagebox.showinfo("编译完成", "编译成功，可进入第 6 步运行", parent=app)
                    elif rc == -9:
                        messagebox.showwarning("已中止", "编译已被用户中止", parent=app)
                        set_light(5, False)
                    else:
                        messagebox.showerror("编译失败", f"退出码 {rc}，详见本页日志", parent=app)
                        set_light(5, False)
                elif kind == "done":
                    rc, tag = payload
                    state["running"] = False
                    prog.stop()
                    prog.config(mode="determinate", value=100 if rc == 0 else 0)
                    btn_state["stop"].config(state="disabled")
                    refresh_detection()
                    refresh_results(tag)
                    elapsed = time.time() - (state["t0"] or time.time())
                    if rc == 0:
                        messagebox.showinfo("完成", f"任务成功结束（退出码 0，耗时 {elapsed:.1f} s）", parent=app)
                    elif rc == -9:
                        messagebox.showwarning("已中止", "任务已被用户中止", parent=app)
                    else:
                        messagebox.showerror("失败", f"任务失败，退出码 {rc}。详见日志。", parent=app)
                elif kind == "scan_progress":
                    i, total = payload
                    state["scan_progress"] = (i, total)
                    prog.config(mode="determinate", maximum=100, value=(100 * i / total if total else 0))
                elif kind == "scan_status":
                    status, i, total, tag = payload
                    state["scan_active"] = (status, i, total, tag)
                elif kind == "scan_done":
                    rows, summary, stopped, axis_keys, metric_names = payload
                    state["running"] = False
                    state["scan_progress"] = None
                    state["scan_active"] = None
                    prog.stop()
                    prog.config(mode="determinate", maximum=100, value=100)
                    btn_state["stop"].config(state="disabled")
                    state["scan_rows"] = rows
                    state["scan_axes"] = axis_keys
                    state["scan_metrics"] = metric_names
                    state["scan_summary"] = summary
                    update_scan_tree(rows, axis_keys, metric_names)
                    refresh_detection()
                    refresh_results(Path(summary).parent.name)
                    elapsed = time.time() - (state["t0"] or time.time())
                    n_ok = sum(1 for r in rows if _scan_success(r.get("status")))
                    n_failed = sum(1 for r in rows if r.get("status") in ("failed", "invalid_output"))
                    n_invalid = sum(1 for r in rows if r.get("status") == "invalid_output")
                    n_unconverged = sum(1 for r in rows if r.get("status") == "not_converged")
                    n_unknown = sum(1 for r in rows if r.get("status") == "convergence_unknown")
                    n_retry = sum(1 for r in rows if r.get("status") == "retry_success")
                    if stopped:
                        messagebox.showwarning("扫描已中止",
                                               f"已完成 {len(rows)} 点（成功 {n_ok} 点、失败 {n_failed} 点，数值无效 {n_invalid} 点、未周期收敛 {n_unconverged} 点、收敛未知 {n_unknown} 点、重试成功 {n_retry} 点），其余已跳过。\n汇总: {summary}",
                                               parent=app)
                    else:
                        messagebox.showinfo("扫描完成",
                                            f"共 {len(rows)} 点（成功 {n_ok} 点、失败 {n_failed} 点，数值无效 {n_invalid} 点、未周期收敛 {n_unconverged} 点、收敛未知 {n_unknown} 点、重试成功 {n_retry} 点），耗时 {elapsed:.1f} s\n汇总: {summary}",
                                            parent=app)
                elif kind == "scan_error":
                    state["running"] = False
                    state["scan_progress"] = None
                    prog.stop()
                    prog.config(mode="determinate", maximum=100, value=0)
                    btn_state["stop"].config(state="disabled")
                    refresh_detection()
                    messagebox.showerror("扫描失败", payload[0], parent=app)
        except queue.Empty:
            pass
        if state["running"] and state["t0"]:
            sp = state.get("scan_progress")
            if sp:
                done, total = sp
                elapsed = time.time() - state["t0"]
                eta = (elapsed / done * (total - done)) if done else None
                active = state.get("scan_active")
                active_txt = ""
                if active:
                    phase, _i, _total, tag = active
                    active_txt = f" · {'重试' if phase == 'retry' else '运行'} {tag}"
                eta_txt = f" · 预计剩余 {eta:.0f} s" if eta is not None else " · 正在估算剩余时间"
                lbl_time.config(text=f"扫描 {done}/{total}（{100 * done / total:.0f}%）{active_txt} · 已耗时 {elapsed:.0f} s{eta_txt}")
            else:
                sim = state.get("sim_progress")
                if sim:
                    now, total = sim
                    lbl_time.config(text=f"计算 {now:g}/{total:g} s（{100 * now / total:.0f}%） · 已耗时 {time.time() - state['t0']:.0f} s")
                else:
                    lbl_time.config(text=f"运行中… 已耗时 {time.time() - state['t0']:.0f} s")
        app.after(120, poll)

    var_dir.trace_add("write", lambda *_: refresh_wizard_status())
    var_batch.trace_add("write", lambda *_: refresh_wizard_status())
    tab6.bind_wheel()  # v3 需求 2：内容建完后统一绑定鼠标滚轮
    app.after(120, poll)
    refresh_scan_panel()
    refresh_wizard_status()
    run_env_check(popup=False)

    # 供测试/冒烟访问
    app._state = state
    app._nb = nb
    app._lights = lights
    app._var_dir = var_dir
    app._var_batch = var_batch
    app._var_distro = var_distro
    app._eff_dir = eff_dir
    app._slot_overrides = slot_overrides
    app._refresh_slots = refresh_slots
    app._slot_rows = slot_rows
    app._link_state = link_state
    app._lbl_src2 = lbl_src2
    app._lbl_src3 = lbl_src3
    app._kinet_text = kinet_text
    app._k_save = k_save
    app._k_load_template = k_load_template
    app._preview = preview
    app._comp_rows = comp_rows
    app._add_comp_row = add_comp_row
    app._set_composition_rows = set_composition_rows
    app._gen_main = gen_main
    app._refresh_wizard_status = refresh_wizard_status
    app._var_pulse = var_pulse
    app._toggle_pulse = toggle_pulse
    app._frm_pulse = frm_pulse
    app._var_scan = var_scan
    app._var_ndim = var_ndim
    app._scan_axis_combos = scan_axis_combos
    app._scan_axis_entries = scan_axis_entries
    app._read_scan_axes = read_scan_axes
    app._lbl_points = lbl_points
    app._update_points = update_points_label
    app._refresh_scan_panel = refresh_scan_panel
    app._update_scan_tree = update_scan_tree
    app._out_vars = out_vars
    app._lbl_warn_a = lbl_warn_a
    app._update_module_a_warn = update_module_a_warn
    app._update_param_state = update_param_state
    app._frm_b_auto = frm_b_auto
    app._lbl_b_auto = lbl_b_auto
    app._btn_back4 = btn_back4
    app._tab6_canvas = tab6_canvas
    app._entries = entries
    app._var_main_kind = var_main_kind
    app._def_entries = def_entries
    app._lbl_ne_mode_note = lbl_ne_mode_note
    app._run_env_check = run_env_check
    app._lbl_env = lbl_env
    app._debug = bool(debug)
    if debug:
        print("[GUI 调试] 主窗口已隐藏；环境自检已完成。"
              " 可通过 zdp_cli.py 执行无界面构建/运行，或从自动化脚本驱动此 GUI 实例。", flush=True)
    return app


def parse_gui_args(argv=None):
    """解析 GUI 启动参数；单独保留便于无界面单测。"""
    parser = argparse.ArgumentParser(description="ZDPlasKin 图形控制台")
    parser.add_argument("--debug", action="store_true",
                        help="隐藏主窗口，并将 GUI 构建/运行日志同步输出到终端")
    return parser.parse_args(argv)


if __name__ == "__main__":
    gui_args = parse_gui_args()
    make_gui(debug=gui_args.debug).mainloop()
