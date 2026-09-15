#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zdp_gui 无界面逻辑单测（临时测试脚本）"""
import os
import sys
import tempfile
import threading
import time
from pathlib import Path

TOOL_DIR = Path(__file__).resolve().parent
PROJECT_DIR = TOOL_DIR.parent
sys.path.insert(0, str(TOOL_DIR))
import zdp_gui as g

TSCAN = PROJECT_DIR / "Reproduction" / "2017Hong" / "方案1_振动开关温度" / "build_tscan"
BUILD = PROJECT_DIR / "Reproduction" / "2017Hong" / "build"
PULSE = PROJECT_DIR / "Reproduction" / "2017Hong" / "方案2_脉冲能效地图" / "build_pulse"


def arg_value(args, flag):
    return args[args.index(flag) + 1]


def write_valid_scan_output(args, cwd, flavor="tscan", nh3="1e13"):
    """为 execute_scan 的假 runner 写最小但数值有效的运行产物。"""
    outdir = Path(cwd) / arg_value(args, "--out-rel")
    tag = arg_value(args, "--tag")
    outdir.mkdir(parents=True, exist_ok=True)
    if flavor in ("pulse", "autopulse"):
        series = outdir / f"pulse_series_{tag}.csv"
        cycles = outdir / f"pulse_cycles_{tag}.csv"
        freq = float(arg_value(args, "--freq"))
        ncycles = int(float(arg_value(args, "--cycles")))
        cycle_rows = ["cycle,t_end_s,E_on_Jcm3,E_off_Jcm3,rel_state_max,rel_dNH3_cycle,rel_energy_cycle,stable_streak"]
        for i in range(1, ncycles + 1):
            cycle_rows.append(f"{i},{i / freq:.12g},1.0,0.5,0,0,0,3")
        cycles.write_text("\n".join(cycle_rows) + "\n", encoding="utf-8")
        (outdir / "console.log").write_text(
            f"cycles_used = {ncycles} converged = 1\n", encoding="utf-8")
        tend = ncycles / freq
    else:
        name = {"hong": "hong_output", "auto": "auto_output"}.get(flavor, "tscan_output")
        series = outdir / f"{name}_{tag}.csv"
        tend = float(arg_value(args, "--tend"))
    series.write_text(
        f"time_s,N2,NH3,Te_eV\n0,1e19,0,1.0\n{tend:.12g},9.9e18,{nh3},1.1\n",
        encoding="utf-8")

# 1. Python 原生运行时（Git Bash 不是 GUI 依赖）
assert g.PYTHON_CLI.is_file() and (TOOL_DIR / "zdp_runtime.py").is_file()
assert g.parse_gui_args([]).debug is False
assert g.parse_gui_args(["--debug"]).debug is True
print("✓ Python 原生运行时入口存在；Git Bash 非必需")

# 2. detect_inputs（build_tscan 应全部满足）
checks, can_build = g.detect_inputs(TSCAN)
for label, ok, detail in checks:
    print(f"  {'✓' if ok else '✗'} {label} {detail}")
assert can_build, "build_tscan 应满足可编译条件"
print("✓ detect_inputs: build_tscan 可编译")

# 3. detect_flavor 三种类型
assert g.detect_flavor(TSCAN) == "tscan", g.detect_flavor(TSCAN)
assert g.detect_flavor(BUILD) == "hong", g.detect_flavor(BUILD)
assert g.detect_flavor(PULSE) == "pulse", g.detect_flavor(PULSE)
print("✓ detect_flavor: tscan/hong/pulse 识别正确")

# 4. validate_params
p = g.validate_params({"en": "60", "tg": "500", "n2frac": "0.3333", "tend": "1e-6",
                       "tag": "guitest", "en_off": "", "freq": "", "duty": "", "cycles": "",
                       "atol": "", "rtol": "", "en_list": ""})
assert p["en"] == "60" and p["tag"] == "guitest"
try:
    g.validate_params({"en": "abc", "tg": "300", "n2frac": "0.3", "tend": "1", "tag": "x",
                       "en_off": "", "freq": "", "duty": "", "cycles": "", "atol": "", "rtol": "", "en_list": ""})
    raise AssertionError("非法 EN 应报错")
except ValueError:
    pass
try:
    g.validate_params({"en": "60", "tg": "300", "n2frac": "0.3", "tend": "1", "tag": "x",
                       "en_off": "", "freq": "", "duty": "", "cycles": "",
                       "atol": "1e3", "rtol": "", "en_list": ""})
    raise AssertionError("只给 atol 应报错")
except ValueError:
    pass
try:
    g.validate_params({"en": "60", "tg": "300", "n2frac": "0.3", "tend": "1", "tag": "非法 tag",
                       "en_off": "", "freq": "", "duty": "", "cycles": "", "atol": "", "rtol": "", "en_list": ""})
    raise AssertionError("非法 tag 应报错")
except ValueError:
    pass
p2 = g.validate_params({"en": "45.1", "tg": "300", "n2frac": "0.3333", "tend": "1", "tag": "s",
                        "en_off": "", "freq": "", "duty": "", "cycles": "", "atol": "", "rtol": "",
                        "en_list": "30 60, 90"})
assert p2["en_list"] == "30 60 90"
print("✓ validate_params: 合法/非法分支正确")

# 6. build_run_args
args = g.build_run_args("", p, "tscan")
assert args[:3] == [sys.executable, str(g.PYTHON_CLI), "run"]
assert "--en" in args and "60" in args and "--tend" in args and "1e-6" in args
assert "--freq" not in args, "tscan 型不应带脉冲参数"
args_p_off = g.build_run_args("", {**p, "freq": "5000", "duty": "0.2", "cycles": "2", "en_off": "0.1"}, "pulse")
assert "--freq" in args_p_off and arg_value(args_p_off, "--cycles") == "2", "CW 仍须传递周期数"
assert "--cw" in args_p_off, "pulse 型未勾选脉冲模式时应带 --cw（连续恒定场）"
args_p = g.build_run_args("", {**p, "freq": "5000", "duty": "0.2", "cycles": "2",
                                        "en_off": "0.1", "pulse_enabled": True}, "pulse")
assert "--freq" in args_p and "5000" in args_p and "--ne-mode" in args_p
assert "--cw" not in args_p, "勾选脉冲模式时不应带 --cw"
args_p_flow = g.build_run_args("", {**p, "freq": "5000", "duty": "0.2", "cycles": "2",
                                             "en_off": "0.1", "tau_res": "1e-3", "pulse_enabled": True}, "pulse")
assert arg_value(args_p_flow, "--tau-res") == "1e-3"
args_r = g.build_run_args("", {**p, "recompile": True}, "hong")
assert "--recompile" in args_r
print("✓ build_run_args:", " ".join(args))

# 7. parse_progress
assert g.parse_progress(" time =   5.0E-04", 1e-3) == 0.5
assert g.parse_progress("NH3 final (cm-3): 1e14", 1e-3) is None
print("✓ parse_progress 正确")

# 8. list_outputs（空 tag 不炸）
assert isinstance(g.list_outputs(TSCAN, "no_such_tag_xyz"), list)
print("✓ list_outputs 正常")

# 8b. 快速绘图数据：递归发现扫描点产物，过滤 t=0 / NaN / Inf 后保留有效正值
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    old = td / "tscan_output_old.csv"
    old.write_text("time_s,Te_eV,EN_Td,Ptot,Pelast,Pinel,N2,NH3\n1,1,60,0,0,0,1e19,1e12\n", encoding="utf-8")
    nested = td / "runs" / "batch" / "EN60" / "tscan_output_EN60.csv"
    nested.parent.mkdir(parents=True)
    nested.write_text(
        "time_s,Te_eV,EN_Td,Ptot,Pelast,Pinel,N2,NH3\n"
        "0,1,60,0,0,0,1e19,0\n"
        "1,1,60,0,0,0,inf,nan\n"
        "2,1,60,0,0,0,9e18,1e13\n", encoding="utf-8")
    os.utime(old, (1, 1)); os.utime(nested, (2, 2))
    targets = g.find_plot_targets(td)
    assert targets[0] == nested and old in targets, targets
    t8, s8 = g._load_density_series(nested)
    assert len(t8) == 3 and "NH3" in s8
    assert g._positive_finite_points(t8, s8["NH3"]) == [(2.0, 1e13)]
    assert g._positive_finite_points(t8, s8["N2"]) == [(2.0, 9e18)]
print("✓ 快速绘图：扫描子目录发现与 0/NaN/Inf 过滤正确")

# 9. GUI 冒烟：实例化 + 脉冲勾选切换可见性
try:
    app = g.make_gui()
except Exception as e:
    print(f"- GUI 冒烟跳过（无显示环境）：{e}")
    app = None
if app is not None:
    app.update()
    assert app.title() == "ZDPlasKin 通用计算控制台"
    assert app._frm_pulse.winfo_manager() == "", "脉冲参数区应默认隐藏"
    app._var_pulse.set(True); app._toggle_pulse(); app.update()
    assert app._frm_pulse.winfo_manager() == "grid", "勾选后脉冲参数区应显示"
    app._var_pulse.set(False); app._toggle_pulse(); app.update()
    assert app._frm_pulse.winfo_manager() == "", "取消勾选后脉冲参数区应隐藏"
    print("✓ GUI 冒烟：默认隐藏 / 勾选显示 / 取消隐藏 均正确")

    # 10. 扫描面板冒烟：tscan 型 2D 轴选择与点数估算
    app._state["flavor"] = "tscan"
    app._refresh_scan_panel(); app.update()
    app._var_scan.set(True)
    app._var_ndim.set(2); app.update()
    app._scan_axis_combos[0].set("E/N (Td)")
    app._scan_axis_combos[1].set("气体温度 (K)")
    app._scan_axis_entries[0].delete(0, "end"); app._scan_axis_entries[0].insert(0, "30:90:3")
    app._scan_axis_entries[1].delete(0, "end"); app._scan_axis_entries[1].insert(0, "300 400")
    app._update_points(); app.update()
    axes = app._read_scan_axes()
    assert [k for k, _ in axes] == ["en", "tg"]
    assert len(axes[0][1]) == 3 and len(axes[1][1]) == 2
    assert "6" in app._lbl_points.cget("text"), app._lbl_points.cget("text")
    # pulse 型未勾选时脉冲轴不可选
    app._state["flavor"] = "pulse"
    app._var_pulse.set(False); app._refresh_scan_panel(); app.update()
    vals = list(app._scan_axis_combos[0].cget("values"))
    assert "占空比" not in vals and "脉冲频率 (Hz)" not in vals, vals
    # 汇总 Treeview 结构
    rows = [{"en": 30.0, "tg": 300.0, "NH3_final_cm-3": 1.2e13, "rc": 0, "tag": "t_i1"},
            {"en": 60.0, "tg": 300.0, "NH3_final_cm-3": float("nan"), "rc": 1, "tag": "t_i2"}]
    cols = app._update_scan_tree(rows, ["en", "tg"], ["NH3_final_cm-3"])
    assert cols == ["en", "tg", "NH3_final_cm-3", "status", "failure_reason", "attempts", "elapsed_s", "rc", "tag"]
    print("✓ 扫描面板冒烟：轴选择 / 点数估算 / 脉冲轴联动 / 汇总表结构 均正确")

    # 10b. 向导弹烟：逐页切换 + 模板载入保存 + 主程序生成预览（屏蔽弹窗）
    import tkinter.messagebox as _mb
    _mb.showinfo = lambda *a, **k: None
    _mb.showerror = lambda *a, **k: None
    _mb.showwarning = lambda *a, **k: None
    _mb.askyesno = lambda *a, **k: True
    for i in range(6):
        app._nb.select(i); app.update()
    with tempfile.TemporaryDirectory() as tdw:
        app._var_dir.set(tdw); app.update()
        app._nb.select(2)
        app._k_load_template(); app.update()
        app._k_save(); app.update()
        st = g.eval_step_status(tdw)
        assert st[3] and not st[1] and not st[4], st
        app._nb.select(3)
        app._set_composition_rows([("N2", 0.3333), ("H2", 0.6667)])
        src = app._gen_main()
        assert src and "{{" not in src and "program main_auto" in src
        assert (Path(tdw) / "main_auto.F90").is_file()
        pv = app._preview.get("1.0", "end")
        assert "MAIN_AUTO" in pv and pv.count("\n") <= 41
        st2 = g.eval_step_status(tdw)
        assert st2[4] and g.wizard_gate(st2) is False  # 缺第 1、2 步，门禁应拒绝
    print("✓ 向导弹烟：逐页切换 / 模板机理 / 主程序生成预览 / 门禁 均正确")
    app.withdraw()  # 保留 after 轮询，避免销毁 Tk 根窗口后留下的回调噪声

# 11. parse_axis_values：列表与区间
assert g.parse_axis_values("30 60 90") == [30.0, 60.0, 90.0]
assert g.parse_axis_values("30,60; 90") == [30.0, 60.0, 90.0]
v = g.parse_axis_values("30:120:5")
assert len(v) == 5 and v[0] == 30.0 and abs(v[-1] - 120.0) < 1e-9 and abs(v[1] - 52.5) < 1e-9
assert g.parse_axis_values("45:45:1") == [45.0]
for bad in ("", "30:120", "abc", "30:120:x", "30:120:0"):
    try:
        g.parse_axis_values(bad)
        raise AssertionError(f"应拒绝: {bad!r}")
    except ValueError:
        pass
print("✓ parse_axis_values: 列表/区间/非法分支正确")

# 12. build_scan_grid：1/2/3 维笛卡尔积与去重校验
g1 = g.build_scan_grid([("en", [30.0, 60.0, 90.0])])
assert len(g1) == 3 and g1[0] == {"en": 30.0}
g2 = g.build_scan_grid([("en", [60.0, 120.0]), ("duty", [0.2, 0.5])])
assert len(g2) == 4 and g2[0] == {"en": 60.0, "duty": 0.2} and g2[-1] == {"en": 120.0, "duty": 0.5}
g3 = g.build_scan_grid([("en", [1.0, 2.0]), ("tg", [3.0, 4.0]), ("n2frac", [5.0, 6.0])])
assert len(g3) == 8 and g3[-1] == {"en": 2.0, "tg": 4.0, "n2frac": 6.0}
for bad_axes in ([("en", [1.0]), ("en", [2.0])], [], [("a", [1]), ("b", [2]), ("c", [3]), ("d", [4])], [("en", [])]):
    try:
        g.build_scan_grid(bad_axes)
        raise AssertionError(f"应拒绝: {bad_axes}")
    except ValueError:
        pass
print("✓ build_scan_grid: 1/2/3 维笛卡尔积正确")

# 13. scan_axes_for_flavor（v3 需求 4：轴扩充 + 按型过滤）
assert [k for k, _ in g.scan_axes_for_flavor("tscan")] == ["en", "tg", "n2frac", "tend"]
assert [k for k, _ in g.scan_axes_for_flavor("pulse")] == \
    ["en", "tg", "n2frac", "freq", "duty", "en_off", "cycles", "tau_res"], "pulse 应含 cycles/tau_res 轴且不含 tend"
assert [k for k, _ in g.scan_axes_for_flavor("hong")] == ["n2frac", "tend"], "hong 仅 n2frac/tend 可扫"
assert [k for k, _ in g.scan_axes_for_flavor("auto")] == ["en", "tg", "tend"], "auto 组分固化，n2frac 不可扫"
assert g.scan_axes_for_flavor(None) == []
print("✓ scan_axes_for_flavor 正确")

# 14. extract_metrics：真实产物断言
m_pulse = g.extract_metrics("pulse", PULSE / "runs" / "guipulsetest", "guipulsetest")
assert abs(m_pulse["NH3_final_cm-3"] - 1.04845e13) / 1.04845e13 < 1e-3, m_pulse
assert m_pulse["NH3_peak_cm-3"] >= m_pulse["NH3_final_cm-3"]
assert abs(m_pulse["Te_peak_eV"] - 2.79723) < 1e-3, m_pulse
assert m_pulse["E_tot_Jcm3"] > 0 and m_pulse["cycles_used"] == 2.0
assert m_pulse["NH3_per_J"] > 0
m_tscan = g.extract_metrics("tscan", TSCAN / "runs" / "guitest", "guitest")
assert m_tscan["NH3_final_cm-3"] > 0 and m_tscan["NH3_peak_cm-3"] >= m_tscan["NH3_final_cm-3"]
assert abs(m_tscan["Te_final_eV"] - 1.19674) < 1e-3, m_tscan
assert 0 <= m_tscan["N2_conv_pct"] < 100
assert g.validate_run_output("tscan", TSCAN / "runs" / "guitest", "guitest",
                              {"tend": "1e-6"}, m_tscan) == ("success", "")
pulse_status, pulse_reason = g.validate_run_output(
    "pulse", PULSE / "runs" / "guipulsetest", "guipulsetest",
    {"freq": "5000"}, m_pulse)
assert pulse_status == "invalid_output" and "周期数不一致" in pulse_reason, (pulse_status, pulse_reason)
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    (td / "pulse_series_consistent.csv").write_text(
        "time_s,phase,Te_eV,EN_Td,Ptot,N2,NH3\n"
        "0,0,2.0,100,1e-8,1e19,0\n"
        "0.001,1,2.0,100,1e-8,1e19,1e13\n",
        encoding="utf-8")
    (td / "pulse_cycles_consistent.csv").write_text(
        "cycle,t_end_s,E_on_Jcm3,E_off_Jcm3,rel_state_max,rel_dNH3_cycle,rel_energy_cycle,stable_streak\n"
        "1,0.001,1e-11,0,0.1,0.1,0.1,0\n",
        encoding="utf-8")
    (td / "console.log").write_text(
        "cycles_used = 1 converged = 0\nt_end = 1.0E-3\nNH3 final (cm-3): 1.0E13\n",
        encoding="utf-8")
    m_consistent = g.extract_metrics("pulse", td, "consistent")
    s_consistent, r_consistent = g.validate_run_output(
        "pulse", td, "consistent", {"freq": "1000"}, m_consistent)
    assert s_consistent == "not_converged" and "周期稳态" in r_consistent
m_none = g.extract_metrics("pulse", PULSE / "runs" / "no_such_dir_xyz", "xx")
assert all(v != v for v in m_none.values()), "缺失产物应全部 NaN"
assert list(g.extract_metrics("tscan", TSCAN, "x").keys()) == g.metric_names_for_flavor("tscan")
assert "E_tot_Jcm3" in g.metric_names_for_flavor("pulse")
print("✓ extract_metrics / validate_run_output: 真实产物、时间完整性与脉冲未收敛标记正确")

# 15. write_scan_summary + execute_scan 失败点容错（假 runner，不跑真程序）
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    rows_in = [{"en": 30.0, "NH3_final_cm-3": 1.2e13, "rc": 0, "tag": "s_i1"},
               {"en": 60.0, "NH3_final_cm-3": float("nan"), "rc": 1, "tag": "s_i2"}]
    sp = g.write_scan_summary(td / "runs" / "s" / "scan_summary.csv", ["en"],
                              ["NH3_final_cm-3"], rows_in)
    lines = sp.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0] == "en,NH3_final_cm-3,status,failure_reason,attempts,elapsed_s,rc,tag"
    assert "NaN" in lines[2] and lines[1].startswith("30,")

    calls = []
    def fake_runner(args, cwd, log_cb, stop_holder):
        calls.append(list(args))
        if "EN60" in args:
            return 1
        write_valid_scan_output(args, cwd)
        return 0
    base = {"en": "30", "tg": "300", "n2frac": "0.3333", "tend": "1e-6",
            "tag": "s", "pulse_enabled": False, "recompile": False}
    rows, summary, stopped = g.execute_scan(
        td, base, [("en", [30.0, 60.0, 90.0])], "tscan",
        log_cb=lambda _l: None, runner=fake_runner)
    assert not stopped and len(rows) == 3
    assert rows[0]["rc"] == 0 and rows[1]["rc"] == 1 and rows[2]["rc"] == 0
    assert rows[1]["NH3_final_cm-3"] != rows[1]["NH3_final_cm-3"], "失败点指标应为 NaN"
    assert all("--en" in c for c in calls) and "60" in calls[1] and "90" in calls[2]
    assert [r["tag"] for r in rows] == ["EN30", "EN60", "EN90"], "扫描点 tag 应为轴参数组合名"
    assert all(any("--out-rel" == c[i] and c[i + 1] == f"runs/s/{t}"
                   for i in range(len(c) - 1)) for c, t in zip(calls, ["EN30", "EN60", "EN90"])), \
        "扫描点应带 --out-rel runs/<批次>/<参数组合>"
    # 每个扫描点都应写入参数档案（即使失败点也记录 rc）
    for t in ("EN30", "EN60", "EN90"):
        assert (td / "runs" / "s" / t / "计算参数.md").is_file(), f"{t} 缺 计算参数.md"
        assert (td / "runs" / "s" / t / "params.json").is_file(), f"{t} 缺 params.json"
    txt = summary.read_text(encoding="utf-8")
    assert txt.splitlines()[0].startswith("en,NH3_final_cm-3")
    assert "NaN" in txt
    # 中止：第二点被 kill
    holder = {"killed": False}
    def kill_runner(args, cwd, log_cb, stop_holder):
        stop_holder["killed"] = True
        return -9
    rows2, _s2, stopped2 = g.execute_scan(
        td, base, [("en", [30.0, 60.0, 90.0])], "tscan",
        log_cb=lambda _l: None, stop_holder=holder, runner=kill_runner)
    assert stopped2 and len(rows2) == 1 and rows2[0]["rc"] == -9
    # pulse 未勾选却选脉冲轴 → 拒绝
    try:
        g.execute_scan(td, base, [("duty", [0.2, 0.5])], "pulse",
                       log_cb=lambda _l: None, runner=fake_runner)
        raise AssertionError("应拒绝脉冲轴+--cw 组合")
    except ValueError:
        pass
    # hong 型拒绝扫描
    try:
        g.execute_scan(td, base, [("en", [30.0])], "hong", log_cb=lambda _l: None)
        raise AssertionError("应拒绝 hong 扫描")
    except ValueError:
        pass
print("✓ write_scan_summary / execute_scan: 失败容错、中止、守卫分支正确")

# 15b. 并行扫描、百分比回调、失败后二次宽松收敛
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    active = {"now": 0, "peak": 0}
    lock = threading.Lock()
    progress = []
    def parallel_runner(args, cwd, log_cb, stop_holder):
        with lock:
            active["now"] += 1
            active["peak"] = max(active["peak"], active["now"])
        time.sleep(0.04)
        write_valid_scan_output(args, cwd)
        with lock:
            active["now"] -= 1
        return 0
    base_parallel = {"en": "30", "tg": "300", "n2frac": "0.5", "tend": "1e-6",
                     "tag": "parallel", "pulse_enabled": False, "recompile": False,
                     "atol": "1.0", "rtol": "1e-4"}
    rows_p, _sp, stop_p = g.execute_scan(
        td, base_parallel, [("en", [30.0, 60.0, 90.0, 120.0])], "tscan",
        log_cb=lambda _l: None, runner=parallel_runner, max_workers=2,
        progress_cb=lambda done, total, row: progress.append((done, total, row["tag"])))
    assert not stop_p and len(rows_p) == 4 and active["peak"] == 2, active
    assert [p[0] for p in progress] == [1, 2, 3, 4] and all(p[1] == 4 for p in progress)

    retry_calls = []
    def retry_runner(args, cwd, log_cb, stop_holder):
        retry_calls.append(list(args))
        if len(retry_calls) == 1:
            return 1
        write_valid_scan_output(args, cwd)
        return 0
    rows_r, _sr, stop_r = g.execute_scan(
        td, {**base_parallel, "tag": "retry", "retry_on_failure": True,
             "retry_atol": "10", "retry_rtol": "1e-3"}, [("en", [60.0])], "tscan",
        log_cb=lambda _l: None, runner=retry_runner)
    assert not stop_r and rows_r[0]["status"] == "retry_success" and rows_r[0]["attempts"] == 2
    assert "--atol" in retry_calls[1] and "10" in retry_calls[1] and "1e-3" in retry_calls[1]

    # rc=0 但 NaN/Inf：也应进入二次收敛；二次仍无效不得阻断后续点。
    assert "NaN/Inf" in g.validate_scan_metrics("tscan", {
        "NH3_final_cm-3": float("inf"), "NH3_peak_cm-3": 1.0,
        "Te_final_eV": 1.0, "N2_conv_pct": 1.0})
    invalid_calls = []
    def invalid_runner(args, cwd, log_cb, stop_holder):
        invalid_calls.append(arg_value(args, "--tag"))
        if arg_value(args, "--tag") == "EN60":
            # 首次和二次都返回 0，但故意不写主输出，模拟大量 NaN。
            return 0
        write_valid_scan_output(args, cwd)
        return 0
    rows_i, summary_i, stop_i = g.execute_scan(
        td, {**base_parallel, "tag": "invalid", "retry_on_failure": True,
             "retry_atol": "10", "retry_rtol": "1e-3"},
        [("en", [60.0, 90.0])], "tscan", log_cb=lambda _l: None,
        runner=invalid_runner)
    assert not stop_i and [r["status"] for r in rows_i] == ["invalid_output", "success"]
    assert rows_i[0]["attempts"] == 2 and "二次收敛后仍数值无效" in rows_i[0]["failure_reason"]
    assert rows_i[0]["NH3_final_cm-3"] != rows_i[0]["NH3_final_cm-3"]
    assert invalid_calls == ["EN60", "EN60", "EN90"], invalid_calls
    assert "invalid_output" in summary_i.read_text(encoding="utf-8")

    # 正常退出但提前终止、或出现显著负密度，均不得作为可用扫描点。
    short_dir = td / "short"
    short_dir.mkdir()
    (short_dir / "tscan_output_short.csv").write_text(
        "time_s,N2,NH3,Te_eV\n0,1e19,0,1\n0.2,9e18,1e13,1.1\n", encoding="utf-8")
    short_metrics = g.extract_metrics("tscan", short_dir, "short")
    short_status, short_reason = g.validate_run_output(
        "tscan", short_dir, "short", {"tend": "1"}, short_metrics)
    assert short_status == "invalid_output" and "未达到目标" in short_reason
    (short_dir / "tscan_output_short.csv").write_text(
        "time_s,N2,NH3,Te_eV\n0,1e19,0,1\n1,-1e10,1e13,1.1\n", encoding="utf-8")
    neg_status, neg_reason = g.validate_run_output(
        "tscan", short_dir, "short", {"tend": "1"},
        g.extract_metrics("tscan", short_dir, "short"))
    assert neg_status == "invalid_output" and "显著为负" in neg_reason
print("✓ execute_scan: 受控并行、完成百分比回调、失败后二次宽松收敛正确")

# 16. parse_kinet_species（真实 65 版机理 + 模板机理）
K65 = PROJECT_DIR / "Reproduction" / "2017Hong" / "build" / "kinet.inp"
sp65 = g.parse_kinet_species(K65)
assert "N2" in sp65 and "H2" in sp65 and "E" in sp65 and len(sp65) >= 30, sp65[:8]
sp_tpl = g.parse_kinet_species(g.TEMPLATE_KINET)
assert sp_tpl == ["E", "N2", "H2", "N2^+", "H2^+"], sp_tpl
assert g.parse_kinet_species(K65.parent / "no_such.inp") == []
print(f"✓ parse_kinet_species: 65 版 {len(sp65)} 物种，模板 {len(sp_tpl)} 物种")

# 16b. 机理反应统计（总反应 / BOLSIG+ 反应）与 N2/H2 默认组成
rs65 = g.parse_kinet_reaction_stats(K65)
assert rs65 == {"total": 65, "bolsig": 10, "fixed": 55, "parsed": True}, rs65
assert g.default_composition_for_species(["E", "N2", "H2"]) == [("N2", "0.5"), ("H2", "0.5")]
assert g.default_composition_for_species(["E", "Ar"]) == []
print("✓ parse_kinet_reaction_stats / default_composition: 总反应、BOLSIG+ 反应与 0.5/0.5 默认值正确")

# 17. validate_composition
comp = g.validate_composition([("N2", "0.3333"), ("H2", "0.6667")])
assert comp == [("N2", 0.3333), ("H2", 0.6667)]
for bad in ([("N2", "0.5"), ("H2", "0.4")], [("N2", "0.5"), ("N2", "0.5")],
            [("N2", "abc"), ("H2", "1")], [], [("N2", "1.5")]):
    try:
        g.validate_composition(bad)
        raise AssertionError(f"应拒绝: {bad}")
    except ValueError:
        pass
print("✓ validate_composition: 合法/非法分支正确")

# 18. _f90_num 与 render_main_auto
assert g._f90_num("45.1") == "45.1d0"
assert g._f90_num("1.17e8") == "1.17d8"
assert g._f90_num("1e-6") == "1d-6"
src = g.render_main_auto(comp, "60", "400", "1e-3", "1.17e8")
assert "{{" not in src and "}}" not in src, "占位符残留"
assert "program main_auto" in src and "MAIN_AUTO" in src
assert "call ZDPlasKin_set_density('N2', 0.3333d0*ntot)" in src
assert "call ZDPlasKin_set_density('H2', 0.6667d0*ntot)" in src
assert "Tgas     = 400d0" in src and "EN_set   = 60d0" in src and "time_end = 0.001d0" in src
assert "species_name(species_electrons)" in src  # 电子物种名机制无关写法
assert "PROGRESS time =" in src and "next_progress" in src, "恒定场模板应每 1 s 输出模拟时间进度"
# 特殊物种自动初始化：机理含 M/S 时自动补密度；不含则不写
src_ms = g.render_main_auto(comp, "60", "400", "1e-3", "1.17e8",
                            species=["E", "N2", "H2", "M", "S", "NS"])
assert "call ZDPlasKin_set_density('M', ntot, ldens_const=.true.)" in src_ms
assert "call ZDPlasKin_set_density('S', 1d16)" in src_ms
src_noms = g.render_main_auto(comp, "60", "400", "1e-3", "1.17e8", species=["E", "N2", "H2"])
assert "'M', ntot" not in src_noms and "'S'," not in src_noms
# 用户已把 M 列入组分时，不重复自动补
src_dup = g.render_main_auto([("N2", 0.5), ("M", 0.5)], "60", "400", "1e-3", "1.17e8",
                             species=["E", "N2", "M", "S"])
assert src_dup.count("set_density('M'") == 1 and "set_density('S'" in src_dup
print("✓ render_main_auto: 无占位符残留，关键行与特殊物种自动初始化齐全")

assert "--no-log" in g.build_run_args("", {**p, "out_console": False}, "tscan")
print("✓ Python 原生运行命令支持控制台日志开关")

# 19. merge_bolsigdb + deploy_program_files + eval_step_status + wizard_gate
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    t1 = td / "a.txt"; t1.write_text("N2\nELASTIC\nN2\n 1.0 2.0\n-----\n", encoding="utf-8")
    t2 = td / "b.txt"; t2.write_text("H2\nIONIZATION\nH2\n 3.0 4.0\n-----\n", encoding="utf-8")
    out, nb_ = g.merge_bolsigdb([t1, t2], td / "bolsigdb.dat")
    merged = out.read_text(encoding="utf-8")
    assert "ELASTIC" in merged and "IONIZATION" in merged and nb_ == 2
    try:
        g.merge_bolsigdb([td / "no.txt"], td / "x.dat")
        raise AssertionError("应拒绝缺失文件")
    except FileNotFoundError:
        pass

    DISTRO = PROJECT_DIR / "1.ZDPlasKin" / "ZDPlasKin_2.0a_Windows"
    case = td / "case"
    copied, skipped, missing = g.deploy_program_files(DISTRO, case)
    assert not missing and "preprocessor.exe" in copied and "dvode_f90_m.F90" in copied
    assert (case / "bolsig_x86_64_g.lib").is_file() and (case / "bolsig_x86_64_g.dll").is_file()
    copied2, skipped2, _ = g.deploy_program_files(DISTRO, case)
    assert copied2 == [] and "preprocessor.exe" in skipped2
    try:
        g.deploy_program_files(td / "no_such_dir", td / "c2")
        raise AssertionError("应拒绝缺失发行包")
    except FileNotFoundError:
        pass

    st = g.eval_step_status(TSCAN)
    assert st == {1: True, 2: True, 3: True, 4: True} and g.wizard_gate(st)
    st_empty = g.eval_step_status(td / "nothing_here")
    assert not any(st_empty.values()) and not g.wizard_gate(st_empty)
print("✓ merge_bolsigdb / deploy_program_files / eval_step_status / wizard_gate 正确")

# 20. combo_name / out_rel_for（需求 1 命名规则）
_p = {"en": "60", "tg": "300", "n2frac": "0.5", "tend": "1e-6", "tag": "b1",
      "freq": "5000", "duty": "0.2", "pulse_enabled": True}
assert g.combo_name(_p, "tscan") == "EN60_Tg300_n2-0.5"
assert g.combo_name(_p, "auto") == "EN60_Tg300_n2-0.5"
assert g.combo_name(_p, "pulse") == "EN60_Tg300_n2-0.5_f5000_d0.2"
assert g.combo_name({**_p, "pulse_enabled": False}, "pulse") == "EN60_Tg300_n2-0.5_cw"
assert g.combo_name(_p, "tscan", ["en", "tg"], {"en": 30.0, "tg": 400.0}) == "EN30_Tg400"
assert g.combo_name(_p, "pulse", ["en", "duty"], {"en": 60.0, "duty": 0.5}) == "EN60_d0.5"
assert g.out_rel_for(_p, "tscan") == "runs/b1/EN60_Tg300_n2-0.5"
assert g.out_rel_for(_p, "tscan", flat=True) == "EN60_Tg300_n2-0.5"
assert g.out_rel_for(_p, "tscan", axis_keys=["en"], point={"en": 30.0}) == "runs/b1/EN30"
print("✓ combo_name / out_rel_for: 无扫描/脉冲/连续/扫描点/扁平 命名正确")

# 21. write_param_record + collect_input_files + sha256_8（需求 1 参数档案）
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    (td / "kinet.inp").write_text("SPECIES\nE N2 H2\nEND\n", encoding="utf-8")
    (td / "bolsigdb.dat").write_text("N2\nELASTIC\n", encoding="utf-8")
    outdir = td / "runs" / "b1" / "EN60_Tg300_n2-0.5"
    outdir.mkdir(parents=True)
    (outdir / "auto_output_x.csv").write_text("time_s,Te_eV\n0,1\n", encoding="utf-8")
    inps = g.collect_input_files(td)
    names = [i["file"] for i in inps]
    assert "kinet.inp" in names and "bolsigdb.dat" in names
    assert all(len(i["sha256_8"]) == 8 for i in inps)
    rec = {"flavor": "auto", "params": {"en": "60", "tg": "300"},
           "scan_axes": None, "point": None,
           "cmdline": [sys.executable, "zdp_cli.py", "run", "--en", "60"],
           "input_files": inps, "t_start": "2026-08-23T10:00:00",
           "t_end": "2026-08-23T10:00:05", "rc": 0, "products": ["auto_output_x.csv"]}
    md, js = g.write_param_record(outdir, rec)
    mdt = md.read_text(encoding="utf-8")
    assert "计算参数档案" in mdt and "E/N" not in mdt  # 键名原样
    assert "- en = 60" in mdt and "zdp_cli.py run --en 60" in mdt
    assert "kinet.inp" in mdt and "auto_output_x.csv" in mdt and "成功" in mdt
    import json as _json
    jr = _json.loads(js.read_text(encoding="utf-8"))
    assert jr["rc"] == 0 and jr["flavor"] == "auto" and jr["params"]["en"] == "60"
print("✓ write_param_record / collect_input_files / sha256_8: md/json 内容齐全")

# 22. detect_slots：自动/覆盖/重置 + zdplaskin_m 不匹配警告 + materialize（需求 2）
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    slots0 = g.detect_slots(td)
    assert slots0["bolsigdb"]["status"] == "missing" and not slots0["bolsigdb"]["ok"]
    assert slots0["kinet"]["status"] == "missing" and slots0["program"]["status"] == "missing"
    assert slots0["zdplaskin_m"]["ok"] and slots0["entropy"]["ok"]  # 可选槽不卡门禁
    (td / "kinet.inp").write_text("SPECIES\nE N2 H2\nEND\n", encoding="utf-8")
    (td / "bolsigdb.dat").write_text("x", encoding="utf-8")
    slots1 = g.detect_slots(td)
    assert slots1["kinet"]["status"] == "auto" and slots1["kinet"]["ok"]
    assert "3 个物种" in slots1["kinet"]["detail"]
    assert slots1["bolsigdb"]["status"] == "auto"
    # zdplaskin_m 与 kinet 物种数不匹配 → ⚠ 警告
    (td / "zdplaskin_m.F90").write_text(
        "  integer, parameter :: species_max = 5, species_electrons = 1\n", encoding="utf-8")
    slots2 = g.detect_slots(td)
    assert "⚠" in slots2["zdplaskin_m"]["detail"] and "不匹配" in slots2["zdplaskin_m"]["detail"]
    assert g.zdplaskin_m_species_max(td / "zdplaskin_m.F90") == 5
    # 手动覆盖（★）与重置
    ext = td / "ext.inp"
    ext.write_text("SPECIES\nE Ar\nEND\n", encoding="utf-8")
    slots3 = g.detect_slots(td, {"kinet": str(ext)})
    assert slots3["kinet"]["status"] == "manual" and slots3["kinet"]["path"] == str(ext)
    slots4 = g.detect_slots(td, {})  # 重置为自动
    assert slots4["kinet"]["status"] == "auto"
    # materialize：覆盖文件复制进算例目录（不覆盖已存在）
    dst = td / "case2"
    copied = g.materialize_slots(dst, {"kinet": str(ext)})
    assert copied == ["kinet.inp"] and (dst / "kinet.inp").read_text(encoding="utf-8") == ext.read_text(encoding="utf-8")
print("✓ detect_slots / materialize_slots: 自动/覆盖/重置/不匹配警告 正确")

# 23. parse_compiled_conditions + module_a_warning（需求 3 模块 A 警告分支）
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    src = td / "main_hong.F90"
    src.write_text("program main_hong\n call ZDPlasKin_set_conditions(GAS_TEMPERATURE=300.0d0, REDUCED_FIELD=45.1d0)\nend\n",
                   encoding="utf-8")
    assert g.parse_compiled_conditions(src) == ("45.1d0", "300.0d0")
    w = g.module_a_warning("hong", src, "60", "300")
    assert w and "⚠" in w and "E/N 编译值" in w
    w2 = g.module_a_warning("hong", src, "45.1", "300")
    assert w2 is None, "与编译值一致时不应警告"
    assert g.module_a_warning("tscan", src, "60", "500") is None, "命令行直传型不警告"
    assert g.module_a_warning("hong", None, "60", "300") is None
print("✓ module_a_warning: hong 不一致警告 / 一致不警告 / 直传型不警告 正确")

# 24. resolve_gfortran_dir 回退链（需求 6，mock 注入）
_gf_cfg = {"gfortran_dir": r"C:\cfg\mingw\bin"}
assert g.resolve_gfortran_dir(_gf_cfg, which_fn=lambda n: None, env={},
                              file_pred=lambda p: p == r"C:\cfg\mingw\bin\gfortran.exe") == r"C:\cfg\mingw\bin"
assert g.resolve_gfortran_dir({}, which_fn=lambda n: None, env={"GFORTRAN_DIR": r"E:\gf\bin"},
                              file_pred=lambda p: p == r"E:\gf\bin\gfortran.exe") == r"E:\gf\bin", "环境变量次之"
assert g.resolve_gfortran_dir({}, which_fn=lambda n: r"D:\x\bin\gfortran.exe", env={},
                              file_pred=lambda p: False) == r"D:\x\bin", "PATH 再次之"
assert g.resolve_gfortran_dir({}, which_fn=lambda n: None, env={},
                              file_pred=lambda p: p == r"C:\mingw64\bin\gfortran.exe") == r"C:\mingw64\bin", "常见位置兜底"
assert g.resolve_gfortran_dir({}, which_fn=lambda n: None, env={}, file_pred=lambda p: False) is None
_sc = g.self_check({})
assert any(n.startswith("模板") and ok for n, ok, _d in _sc), "templates 自检应通过"
assert any(n.startswith("Python 原生运行时") and ok for n, ok, _d in _sc), "Python 原生运行时自检应通过"
print("✓ resolve_gfortran_dir / self_check: Python 原生链路自检正确")

# 25. render_main_auto 输出开关（需求 3 模块 E 对 main_auto 真实接线）
src_off = g.render_main_auto(comp, "60", "300", "1e-3", "1.17e8",
                             output_opts={"species": False, "te": False, "rates_t": False, "rates": False})
assert "@BLOCK" not in src_off and "@ENDBLOCK" not in src_off, "块标记应全部清除"
assert "Te_eV" not in src_off and "species_name(i)" not in src_off
assert "open(unit=urt" not in src_off and "close(urt)" not in src_off, "rates_t 文件块应移除"
assert "open(unit=ur," not in src_off and "close(ur)" not in src_off, "rates 文件块应移除"
src_on = g.render_main_auto(comp, "60", "300", "1e-3", "1.17e8",
                            output_opts={"species": True, "te": True, "rates_t": True, "rates": True})
assert "@BLOCK" not in src_on and "Te_eV" in src_on and "open(unit=urt" in src_on
src_mid = g.render_main_auto(comp, "60", "300", "1e-3", "1.17e8",
                             output_opts={"rates_t": False})
assert "open(unit=urt" not in src_mid and "open(unit=ur," in src_mid and "Te_eV" in src_mid
print("✓ render_main_auto output_opts: 输出开关取舍 write 语句正确")

# 26. v3 需求 3：auto 型 build_run_args 永不含 n2 参数（组分固化在源码，防静默失效）
args_auto = g.build_run_args("", p, "auto")
assert "--n2-frac" not in args_auto, "auto 型不应下发 --n2-frac"
assert "--n2-frac" in g.build_run_args("", p, "tscan")
assert "--n2-frac" in g.build_run_args("", p, "hong")
print("✓ build_run_args: auto 型不含 --n2-frac，其余型含")

# 27. v3 需求 3/4：扫描轴守卫与 n2/cycles 轴网格
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    calls27 = []
    def fake_runner27(args, cwd, log_cb, stop_holder):
        calls27.append(list(args))
        write_valid_scan_output(args, cwd, "pulse" if "--freq" in args else "hong")
        return 0
    base27 = {"en": "45.1", "tg": "300", "n2frac": "0.3333", "tend": "1e-6",
              "tag": "s", "pulse_enabled": False, "recompile": False}
    # auto + n2frac 轴 → 二次拒绝
    try:
        g.execute_scan(td, base27, [("n2frac", [0.3, 0.5])], "auto",
                       log_cb=lambda _l: None, runner=fake_runner27)
        raise AssertionError("应拒绝 auto + n2frac 轴")
    except ValueError as e27:
        assert "固化" in str(e27)
    # hong + en 轴 → 拒绝；hong + n2frac 轴 → 允许（第一位置参数命令行直传）
    try:
        g.execute_scan(td, base27, [("en", [30.0])], "hong",
                       log_cb=lambda _l: None, runner=fake_runner27)
        raise AssertionError("应拒绝 hong + en 轴")
    except ValueError:
        pass
    calls27.clear()
    rows27, _s27, stop27 = g.execute_scan(
        td, base27, [("n2frac", [0.3333, 0.5])], "hong",
        log_cb=lambda _l: None, runner=fake_runner27)
    assert not stop27 and len(rows27) == 2 and all(r["rc"] == 0 for r in rows27)
    assert [r["tag"] for r in rows27] == ["n20.3333", "n20.5"]
    def _argval(call, flag):
        return call[call.index(flag) + 1]
    assert _argval(calls27[0], "--n2-frac") == "0.3333" and _argval(calls27[1], "--n2-frac") == "0.5", \
        "hong 扫描点 --n2-frac 应逐点变化"
    # pulse + cycles 轴：勾选脉冲允许；未勾选（--cw）拒绝
    calls27.clear()
    base27p = {**base27, "pulse_enabled": True, "freq": "5000", "duty": "0.2",
               "cycles": "2", "en_off": "0.1", "ne_mode": "fix"}
    rows27p, _s27p, _st27p = g.execute_scan(
        td, base27p, [("cycles", [2.0, 4.0])], "pulse",
        log_cb=lambda _l: None, runner=fake_runner27)
    assert [r["tag"] for r in rows27p] == ["c2", "c4"]
    assert _argval(calls27[0], "--cycles") == "2" and _argval(calls27[1], "--cycles") == "4"
    try:
        g.execute_scan(td, base27, [("cycles", [2.0])], "pulse",
                       log_cb=lambda _l: None, runner=fake_runner27)
        raise AssertionError("未勾选脉冲应拒绝 cycles 轴")
    except ValueError:
        pass
print("✓ execute_scan 守卫：auto n2frac 拒绝 / hong en 拒绝 / hong n2 与 pulse cycles 网格正确")

# 28. v3 需求 3：parse_auto_composition + fmt_composition
src28 = g.render_main_auto([("N2", 0.3333), ("H2", 0.6667)], "60", "300", "1e-3", "1.17e8",
                           species=["E", "N2", "H2", "M", "S"])
with tempfile.TemporaryDirectory() as td:
    f28 = Path(td) / "main_auto.F90"
    f28.write_text(src28, encoding="utf-8")
    comp28 = g.parse_auto_composition(f28)
    assert comp28 == [("N2", "0.3333d0"), ("H2", "0.6667d0")], comp28  # M/S 自动行不计入
    assert g.fmt_composition(comp28) == "N2 33.33% / H2 66.67%"
assert g.parse_auto_composition(Path(td) / "no_such.F90") == []
assert g.fmt_composition([]) == ""
print("✓ parse_auto_composition / fmt_composition 正确")

# 29. v3 需求 1：auto_link_updates 联动（自动带入 / 手动优先不覆盖 / 缺失不带入）
slots29 = {"bolsigdb": {"status": "auto", "path": "/d/bolsigdb.dat"},
           "kinet": {"status": "auto", "path": "/d/kinet.inp"}}
assert g.auto_link_updates(slots29) == {"bolsigdb": "/d/bolsigdb.dat", "kinet": "/d/kinet.inp"}
assert g.auto_link_updates(slots29, xs_manual=True)["bolsigdb"] is None, "手动指定过绝不覆盖"
assert g.auto_link_updates(slots29, kinet_manual=True)["kinet"] is None
assert g.auto_link_updates(slots29, xs_manual=True)["kinet"] == "/d/kinet.inp", "另一项不受影响"
slots29m = {"bolsigdb": {"status": "manual", "path": "/x.dat"},
            "kinet": {"status": "missing", "path": None}}
assert g.auto_link_updates(slots29m) == {"bolsigdb": None, "kinet": None}, \
    "manual(★)/missing 槽位不自动带入"
print("✓ auto_link_updates: 自动带入 / 手动优先 / 缺失不带入 正确")

# 30. GUI 冒烟：第 6 步滚动容器、auto 组成编辑入口与联动标签
try:
    app3 = g.make_gui()
except Exception as e:
    print(f"- v3 GUI 冒烟跳过（无显示环境）：{e}")
    app3 = None
if app3 is not None:
    import tkinter.messagebox as _mb3
    _mb3.showinfo = lambda *a, **k: None
    _mb3.showerror = lambda *a, **k: None
    _mb3.showwarning = lambda *a, **k: None
    app3.update()
    assert app3._tab6_canvas.winfo_exists(), "第 6 步应有滚动 Canvas"
    assert str(app3._entries["n2frac"].cget("state")) == "normal"
    with tempfile.TemporaryDirectory() as tda:
        src30 = g.render_main_auto([("N2", 0.3333), ("H2", 0.6667)], "60", "300", "1e-3", "1.17e8")
        (Path(tda) / "main_auto.F90").write_text(src30, encoding="utf-8")
        app3._var_dir.set(tda); app3.update()
        app3._state["flavor"] = "auto"
        app3._update_param_state(); app3.update()
        assert str(app3._entries["n2frac"].cget("state")) == "disabled", "auto 型单一 N2 参数应禁用，改用组成编辑入口"
        txt30 = app3._lbl_b_auto.cget("text")
        assert "可在此修改" in txt30 and "33.33%" in txt30, txt30
        assert app3._frm_b_auto.winfo_manager() == "grid", "auto 型应显示组成编辑条"
        # auto 型扫描轴面板剔除 N2 比例
        app3._refresh_scan_panel(); app3.update()
        vals30 = list(app3._scan_axis_combos[0].cget("values"))
        assert not any("N2" in v for v in vals30), vals30
        # 切回 tscan：恢复可编辑、摘要条隐藏
        app3._state["flavor"] = "tscan"
        app3._update_param_state(); app3.update()
        assert str(app3._entries["n2frac"].cget("state")) == "normal"
        assert app3._frm_b_auto.winfo_manager() == ""
        # 联动：kinet.inp / bolsigdb.dat 自动带入第 2/3 步
        (Path(tda) / "kinet.inp").write_text("SPECIES\nE N2 H2\nEND\n", encoding="utf-8")
        (Path(tda) / "bolsigdb.dat").write_text("N2\nELASTIC\n", encoding="utf-8")
        app3._refresh_slots(); app3.update()
        assert app3._link_state["xs"] == "auto" and "自动识别" in app3._lbl_src2.cget("text")
        assert app3._link_state["kinet"] == "auto" and "自动识别" in app3._lbl_src3.cget("text")
        assert "SPECIES" in app3._kinet_text.get("1.0", "end"), "第 3 步应自动载入机理内容"
        # 手动优先：标记手动后自动识别不覆盖
        app3._link_state["xs"] = "manual"
        app3._refresh_slots(); app3.update()
        assert "手动指定" in app3._lbl_src2.cget("text")
    app3.withdraw()  # 同上
    print("✓ GUI 冒烟：滚动容器 / auto 组成编辑入口 / 轴剔除 / 联动与手动优先 均正确")

# 31. render_main_auto_pulse：无占位符/块标记残留、关键结构在位
src31 = g.render_main_auto_pulse([("N2", 0.3333), ("H2", 0.6667)], en_on="60", en_off="0.1",
                                 freq="5000", duty="0.2", cycles="2",
                                 species=["E", "N2", "H2", "M", "S"])
assert "{{" not in src31 and "}}" not in src31, "占位符残留"
assert "@BLOCK" not in src31 and "@ENDBLOCK" not in src31, "块标记残留"
assert "MAIN_AUTO_PULSE" in src31, "探测标记缺失"
assert "program main_auto_pulse" in src31
assert "EN_on  = 60d0" in src31 and "freq   = 5000d0" in src31 and "duty   = 0.2d0" in src31, \
    [ln for ln in src31.splitlines() if "EN_on" in ln or "freq" in ln or "duty" in ln][:6]
assert "n_max  = 2" in src31, "cycles 应渲染为整数字面量"
assert "set_density('N2', 0.3333d0*ntot)" in src31
assert "set_density('M', ntot" in src31 and "set_density('S', 1d16)" in src31
assert "species_electrons" in src31, "电子应为机制无关写法"
assert "pulse_series_" in src31 and "pulse_cycles_" in src31 and "pulse_rates_" in src31
assert "time_s,dt_s,cycle,phase,reaction,rate_cm-3s-1" in src31
assert "call write_rates(ur, tcur, dt1, k, 'on')" in src31
assert "SOFT_RESET=.true." in src31, "脉冲 E/N 跳变必须重启 DVODE 历史"
assert "cycles_used =', cycles_used" in src31
assert "rel_state_max" in src31 and "stable_streak" in src31
assert "ZDPlasKin_set_cstr_flow" in src31 and "tau_res" in src31
src_pulse = (PULSE / "main_pulse.F90").read_text(encoding="utf-8")
assert "ZDPlasKin_set_cstr_flow" in src_pulse and "tau_res_s" in src_pulse
assert "feed_density" in src_pulse and "flow_enabled" in src_pulse
assert "nh3_idx" in src31, "NH3 应机制无关定位"
assert "PROGRESS time =" in src31 and "next_progress" in src31, "脉冲模板应每 1 s 输出模拟时间进度"
# output_opts 裁剪：species_cols/rates 可移除，rates_t=False 无害（无此块）
src31b = g.render_main_auto_pulse([("N2", 1.0)], output_opts={"species": False, "rates": False,
                                                              "te": False, "rates_t": False})
assert "time_s,dt_s,cycle,phase,reaction,rate_cm-3s-1" not in src31b, "rates 块内容应被裁掉"
assert "',' // trim(adjustl(species_name" not in src31b and "density(kk)" not in src31b, \
    "species_cols 块应被裁掉"
print("✓ render_main_auto_pulse: 渲染无残留、关键结构在位、块裁剪正确")

# 32. detect_flavor 五型：autopulse 不被 MAIN_AUTO/EN_on 子串误判
with tempfile.TemporaryDirectory() as td:
    f32 = Path(td) / "main_auto_pulse.F90"
    f32.write_text(src31, encoding="utf-8")
    assert g.detect_flavor(td) == "autopulse", g.detect_flavor(td)
with tempfile.TemporaryDirectory() as td:
    (Path(td) / "main_auto.F90").write_text(
        g.render_main_auto([("N2", 1.0)]), encoding="utf-8")
    assert g.detect_flavor(td) == "auto", "auto 不应被 autopulse 分支抢走"
print("✓ detect_flavor: autopulse/auto/pulse/tscan/hong 五型区分正确")

# 33. build_run_args autopulse 映射：无 --n2-frac、无 --ne-mode；未勾选带 --cw
p33 = {"en": "60", "tg": "300", "n2frac": "0.3333", "tend": "1", "tag": "t33",
       "en_off": "0.1", "freq": "5000", "duty": "0.2", "cycles": "2", "pulse_enabled": True}
a33 = g.build_run_args("", p33, "autopulse")
assert "--n2-frac" not in a33, "autopulse 组成内置，不应下发 --n2-frac"
assert "--ne-mode" not in a33, "autopulse 电子恒定，不应下发 --ne-mode"
for k in ("--en-off", "--freq", "--duty", "--cycles"):
    assert k in a33, f"{k} 缺失"
a33cw = g.build_run_args("", {**p33, "pulse_enabled": False}, "autopulse")
assert "--cw" in a33cw and "--freq" in a33cw and arg_value(a33cw, "--cycles") == "2"
# pulse 型仍带 --ne-mode（不回归）
a33p = g.build_run_args("", {**p33, "ne_mode": "fix"}, "pulse")
assert "--ne-mode" in a33p and "--n2-frac" in a33p
print("✓ build_run_args: autopulse 映射正确（无 n2/ne-mode，--cw 互斥）")

# 34. scan_axes_for_flavor autopulse + execute_scan 守卫
keys34 = [k for k, _ in g.scan_axes_for_flavor("autopulse")]
assert "n2frac" not in keys34 and "tend" not in keys34, keys34
for k in ("en", "tg", "freq", "duty", "en_off", "cycles"):
    assert k in keys34, f"轴 {k} 缺失"
try:
    g.execute_scan(Path(td), {"tag": "b"}, [("n2frac", [0.5])], "autopulse", runner=lambda *a, **k: 0)
    raise AssertionError("autopulse n2frac 轴应被拒绝")
except ValueError as e:
    assert "固化" in str(e)
try:
    g.execute_scan(Path(td), {"tag": "b", "pulse_enabled": False}, [("duty", [0.2])], "autopulse",
                   runner=lambda *a, **k: 0)
    raise AssertionError("未勾选脉冲时 duty 轴应被拒绝")
except ValueError as e:
    assert "脉冲" in str(e)
# 指标注册：autopulse 与 pulse 同套
assert g.metric_names_for_flavor("autopulse") == g.metric_names_for_flavor("pulse")
assert g._SERIES_FILE["autopulse"] == "pulse_series_{tag}.csv"
print("✓ scan_axes/execute_scan 守卫/指标注册: autopulse 正确")

# 35. GUI 联动：第 4 步脉冲表单随单选启停；autopulse 时 ne_mode 禁用并注明
try:
    app4 = g.make_gui()
except Exception as e:
    print(f"- GUI 联动冒烟跳过（无显示环境）：{e}")
    app4 = None
if app4 is not None:
    app4.update()
    # 默认恒定场：脉冲参数表单禁用
    assert app4._var_main_kind.get() == "cw"
    for k in ("en_off", "freq", "duty", "cycles"):
        assert str(app4._def_entries[k].cget("state")) == "disabled", k
    # 切脉冲：启用
    app4._var_main_kind.set("pulse"); app4.update()
    for k in ("en_off", "freq", "duty", "cycles"):
        assert str(app4._def_entries[k].cget("state")) == "normal", k
    # 切回恒定场：恢复禁用
    app4._var_main_kind.set("cw"); app4.update()
    assert str(app4._def_entries["freq"].cget("state")) == "disabled"
    # 第 6 步：pulse 型 ne_mode 可编辑、无注明
    app4._state["flavor"] = "pulse"
    app4._update_param_state(); app4.update()
    assert str(app4._entries["ne_mode"].cget("state")) == "normal"
    assert app4._lbl_ne_mode_note.cget("text") == ""
    # autopulse 型：ne_mode 禁用 + 注明
    app4._state["flavor"] = "autopulse"
    app4._update_param_state(); app4.update()
    assert str(app4._entries["ne_mode"].cget("state")) == "disabled", "autopulse ne_mode 应禁用"
    assert "恒定电子密度" in app4._lbl_ne_mode_note.cget("text")
    # 切回 pulse：恢复
    app4._state["flavor"] = "pulse"
    app4._update_param_state(); app4.update()
    assert str(app4._entries["ne_mode"].cget("state")) == "normal"
    app4.destroy()
    print("✓ GUI 联动：脉冲表单随单选启停 / autopulse ne_mode 禁用注明 均正确")

print("\n全部单测通过 ✓")
