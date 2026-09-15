#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""方案1：振动催化开关温度 —— T×E/N 扫描通道分解与交叉分析
输入: build_tscan/tscan_rates_t_T{T}_EN{EN}[_aXX].csv + tscan_output_*.csv
输出: data/results/tscan_results.json
      data/results/fig_i_NH_channel_share_vs_T.png
      data/results/fig_j_vdf_vs_T.png
      data/results/fig_k_nh3_yield_energy_vs_T.png
所有数字均读自 CSV，不虚构。"""
import sys, csv, json, re
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from daimon_runtime import setup_plot
setup_plot()

HERE = Path(__file__).parent.resolve()
BD = HERE / "build_tscan"
RES = HERE / "data" / "results"
TS = [300, 400, 500, 600, 700, 800]
ENS = [30, 45.1, 80, 120]
ALPHAS = [("a03", 0.3), ("a05", 0.5), ("a10", 1.0), ("a15", 1.5)]

# ---------------------------------------------------------------- 反应式解析（同 pathway_analysis.py）
def norm_name(s):
    return re.sub(r",V\d[^)]*\)", ")", s.strip())

def split_side(side):
    raw, toks, i = side.split("+"), [], 0
    while i < len(raw):
        t = raw[i]
        if t.endswith("^") and i + 1 < len(raw):
            t = t + "+" + raw[i + 1]; i += 1
        toks.append(t); i += 1
    out = {}
    for t in toks:
        t = norm_name(t)
        if not t: continue
        m = re.match(r"^(\d+)([A-Za-z].*)$", t)
        n, name = (int(m.group(1)), m.group(2)) if m else (1, t)
        out[name] = out.get(name, 0) + n
    return out

def parse_sign(sign):
    if sign.startswith("bolsig:"):
        body = sign[7:]
        if body == "N2->N2(SUM)":
            return {"N2": 1}, {"N": 2}
        a, b = body.split("->")
        return {norm_name(a): 1}, {norm_name(b): 1}
    lhs, rhs = sign.split("=>")
    return split_side(lhs), split_side(rhs)

# ---------------------------------------------------------------- 通道分类
H2V  = {"H2(V1)", "H2(V2)", "H2(V3)"}                                  # 振动
H2EL = {"H2(B3SIG)", "H2(B1SIG)", "H2(C3PI)", "H2(A3SIG)", "H2(RYDBERG_SUM)"}  # 电子态
NATM = {"N", "N(2D)", "N(2P)"}

def classify_nh_prod(R):
    if any(s in NATM for s in R):
        if any(s in H2V for s in R):  return "N+H2(v) 振动"
        if any(s in H2EL for s in R): return "N+H2* 电子态"
        if "H2" in R:                 return "N+H2 基态"
    if "NH3" in R:  return "NH3 解离"
    if "NH2" in R:  return "NH2 转化"
    return "其他"

NH_CLASSES = ["N+H2(v) 振动", "N+H2* 电子态", "N+H2 基态", "NH3 解离", "NH2 转化", "其他"]

def classify_nh3_prod(R):
    return "气相三体" if sum(R.values()) >= 3 else "气相两体"

# ---------------------------------------------------------------- 数据加载
def ffloat(s):
    """解析 Fortran ES 格式，容忍三位指数丢 'E'（如 3.27010-164）。"""
    s = s.strip()
    try:
        return float(s)
    except ValueError:
        m = re.match(r"^([+-]?\d*\.\d+)([+-]\d+)$", s)
        if m:
            return float(m.group(1) + "E" + m.group(2))
        m = re.match(r"^([+-]?\d+)([+-]\d{3,})$", s)
        if m:
            return float(m.group(1) + "E" + m.group(2))
        raise

def load_rates_at(path, t_target=100.0):
    rows = list(csv.reader(open(path, encoding="utf-8")))[1:]
    agg = defaultdict(float)
    for t, sign, rate in rows:
        if abs(ffloat(t) - t_target) / t_target < 1e-6:
            agg[sign] += ffloat(rate)
    out = []
    for sign, rate in agg.items():
        R, P = parse_sign(sign)
        out.append((sign, R, P, rate))
    return out

def load_output(path):
    rows = list(csv.reader(open(path, encoding="utf-8")))
    hdr = [h.strip() for h in rows[0]]
    data = {n: np.array([ffloat(r[i]) for r in rows[1:]]) for i, n in enumerate(hdr)}
    return data

def nh_shares(rxns):
    sh = defaultdict(float); tot = 0.0; detail = defaultdict(list)
    for sign, R, P, rate in rxns:
        net = P.get("NH", 0) - R.get("NH", 0)
        if net > 0:
            c = classify_nh_prod(R)
            sh[c] += rate * net; tot += rate * net
            detail[c].append((sign, rate * net))
    return ({c: 100.0 * sh[c] / tot if tot else 0.0 for c in NH_CLASSES},
            tot, detail)

def nh3_shares(rxns):
    sh = defaultdict(float); tot = 0.0; detail = defaultdict(list)
    for sign, R, P, rate in rxns:
        net = P.get("NH3", 0) - R.get("NH3", 0)
        if net > 0:
            c = classify_nh3_prod(R)
            sh[c] += rate * net; tot += rate * net
            detail[c].append((sign, rate * net))
    return ({c: 100.0 * sh[c] / tot if tot else 0.0 for c in ["气相三体", "气相两体"]},
            tot, detail)

# ---------------------------------------------------------------- 主分析
results = {"baseline": {}, "alpha": {}, "meta": {
    "note": "纯气相 515-46=469 反应；n_e=1.17e8 固定；1 atm；ntot=2.446e19*(300/T)",
    "tol_overrides": {"T600_EN80": "ATOL=1d2 RTOL=1d-3", "T700_EN80": "ATOL=1d2 RTOL=1d-3",
                      "T800_EN80": "ATOL=1d4 RTOL=1d-2"}}}

def tag_of(T, EN, a=None):
    return f"T{T}_EN{EN}" + (f"_{a}" if a else "")

def analyze_point(T, EN, a=None):
    tag = tag_of(T, EN, a)
    rxns = load_rates_at(BD / f"tscan_rates_t_{tag}.csv", 100.0)
    out = load_output(BD / f"tscan_output_{tag}.csv")
    nh_sh, nh_tot, nh_det = nh_shares(rxns)
    nh3_sh, nh3_tot, nh3_det = nh3_shares(rxns)
    # VDF @100s
    n2v = {f"v{i}": float(out[f"N2(V{i})"][-1]) for i in range(1, 9)}
    h2v = {f"v{i}": float(out[f"H2(V{i})"][-1]) for i in range(1, 4)}
    # NH 时序份额（用于追踪）
    nh_sh_t = {}
    for t in [1e-6, 1e-4, 1e-2, 0.1, 1.0, 10.0, 100.0]:
        r_t = load_rates_at(BD / f"tscan_rates_t_{tag}.csv", t)
        s_t, _, _ = nh_shares(r_t)
        nh_sh_t[f"{t:g}"] = s_t
    ptot = float(out["Ptot"][-1])
    nh3_fin = float(out["NH3"][-1])
    return {
        "Te_eV": float(out["Te_eV"][-1]),
        "NH3_100s": nh3_fin, "NH_100s": float(out["NH"][-1]),
        "N_100s": float(out["N"][-1]), "H_100s": float(out["H"][-1]),
        "nh_prod_share_pct": nh_sh, "nh_prod_total": nh_tot,
        "nh_prod_top": {c: sorted(detail, key=lambda x: -x[1])[:3] for c, detail in nh_det.items()},
        "nh_share_vs_time": nh_sh_t,
        "nh3_prod_share_pct": nh3_sh, "nh3_prod_total": nh3_tot,
        "nh3_prod_top": {c: sorted(detail, key=lambda x: -x[1])[:4] for c, detail in nh3_det.items()},
        "N2_VDF": n2v, "H2_VDF": h2v,
        "Ptot_Wcm3": ptot,
        "NH3_per_eV": nh3_fin * 1.602e-19 / (ptot * 100.0) if ptot > 0 else None,
    }

for T in TS:
    for EN in ENS:
        results["baseline"][tag_of(T, EN)] = analyze_point(T, EN)
for a, _ in ALPHAS[1:]:
    for T in TS:
        results["alpha"][tag_of(T, 45.1, a)] = analyze_point(T, 45.1, a)

# ---------------------------------------------------------------- 交叉温度 T*
def vib_share(res):  return res["nh_prod_share_pct"]["N+H2(v) 振动"]

def crossover(get_share, Ts):
    sh = [get_share(T) for T in Ts]
    if sh[0] >= 50.0:
        return float(Ts[0])  # 起点已越过 50%（α=1.5 时的非单调情形）
    for i in range(len(Ts) - 1):
        if sh[i] < 50.0 <= sh[i + 1]:
            return Ts[i] + (50.0 - sh[i]) / (sh[i + 1] - sh[i]) * (Ts[i + 1] - Ts[i])
    return None

tstar = {}
for EN in ENS:
    g = lambda T: vib_share(results["baseline"][tag_of(T, EN)])
    tstar[str(EN)] = {"T_star_K": crossover(g, TS),
                      "vib_share_at_Ts": [vib_share(results["baseline"][tag_of(T, EN)]) for T in TS]}
for a, aval in ALPHAS:
    if a == "a03":
        g = lambda T: vib_share(results["baseline"][tag_of(T, 45.1)])
    else:
        g = lambda T: vib_share(results["alpha"][tag_of(T, 45.1, a)])
    tstar[f"alpha_{aval}"] = {"T_star_K": crossover(g, TS),
                              "vib_share_at_Ts": [g(T) for T in TS]}
results["tstar"] = tstar

with open(RES / "tscan_results.json", "w", encoding="utf-8") as f:
    json.dump(results, f, ensure_ascii=False, indent=2)
print("saved:", RES / "tscan_results.json")

# ---- 摘要打印 ----
print("\n=== NH 振动通道份额 (%) @100s, baseline α=0.3 ===")
print("T\\EN  " + "".join(f"{en:>10}" for en in ENS))
for T in TS:
    print(f"{T:>4}  " + "".join(f"{vib_share(results['baseline'][tag_of(T, en)]):>10.3f}" for en in ENS))
print("\n=== T* (振动份额=50%) ===")
for k, v in tstar.items():
    ts = v["T_star_K"]
    print(f"  {k:>10}: " + (f"{ts:.0f} K" if ts else f"800 K 内未达（800K 时 {v['vib_share_at_Ts'][-1]:.2f}%）"))
print("\n=== Te (eV) @100s ===")
print("T\\EN  " + "".join(f"{en:>10}" for en in ENS))
for T in TS:
    print(f"{T:>4}  " + "".join(f"{results['baseline'][tag_of(T, en)]['Te_eV']:>10.3f}" for en in ENS))
print("\n=== NH3 @100s (cm-3) ===")
for T in TS:
    print(f"{T:>4}  " + "".join(f"{results['baseline'][tag_of(T, en)]['NH3_100s']:>10.2e}" for en in ENS))
print("\n=== α 变体 NH 振动份额 @45.1 Td ===")
print("T    " + "".join(f"α={a:>6}" for _, a in ALPHAS))
for T in TS:
    row = [vib_share(results["baseline"][tag_of(T, 45.1)])] + \
          [vib_share(results["alpha"][tag_of(T, 45.1, a)]) for a, _ in ALPHAS[1:]]
    print(f"{T:>4}  " + "".join(f"{v:>9.3f}" for v in row))

# ---------------------------------------------------------------- fig_i: NH 通道份额 vs T
COLORS = {"N+H2(v) 振动": "#C23B22", "N+H2* 电子态": "#4C9BD6", "N+H2 基态": "#9ECAE1",
          "NH3 解离": "#E8853D", "NH2 转化": "#77B255", "其他": "#999999"}
fig = plt.figure(figsize=(13, 8.5))
gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1], hspace=0.32, wspace=0.28)
# 左上大面板：45.1 Td 堆叠面积
ax0 = fig.add_subplot(gs[0, :2])
stacks = np.zeros(len(TS))
for c in NH_CLASSES:
    vals = np.array([results["baseline"][tag_of(T, 45.1)]["nh_prod_share_pct"][c] for T in TS])
    ax0.fill_between(TS, stacks, stacks + vals, color=COLORS[c], alpha=0.85, label=c)
    stacks += vals
ax0.axhline(50, color="k", lw=0.8, ls=":")
ts45 = tstar["45.1"]["T_star_K"]
if ts45:
    ax0.axvline(ts45, color="#C23B22", lw=1.4, ls="--")
    ax0.annotate(f"T* = {ts45:.0f} K", xy=(ts45, 52), fontsize=10, color="#C23B22")
ax0.set_xlim(300, 800); ax0.set_ylim(0, 100)
ax0.set_xlabel("气体温度 T (K)"); ax0.set_ylabel("NH 生成通道份额 (%)")
ax0.set_title("E/N = 45.1 Td（论文基线）· α=0.3", fontsize=11)
ax0.legend(fontsize=8, loc="center left", bbox_to_anchor=(0.01, 0.55))
# 右上：α 敏感性（振动份额曲线）
ax1 = fig.add_subplot(gs[0, 2])
for a, aval in ALPHAS:
    if a == "a03":
        vals = [vib_share(results["baseline"][tag_of(T, 45.1)]) for T in TS]
    else:
        vals = [vib_share(results["alpha"][tag_of(T, 45.1, a)]) for T in TS]
    ax1.plot(TS, vals, "o-", ms=4, lw=1.6, label=f"α={aval}")
ax1.axhline(50, color="k", lw=0.8, ls=":")
ax1.set_xlabel("T (K)"); ax1.set_ylabel("振动通道份额 (%)")
ax1.set_title("α 敏感性 @45.1 Td", fontsize=11)
ax1.legend(fontsize=8); ax1.grid(alpha=0.3)
ax1.set_yscale("log")
# 下行：其余 E/N
for j, EN in enumerate([30, 80, 120]):
    ax = fig.add_subplot(gs[1, j])
    stacks = np.zeros(len(TS))
    for c in NH_CLASSES:
        vals = np.array([results["baseline"][tag_of(T, EN)]["nh_prod_share_pct"][c] for T in TS])
        ax.fill_between(TS, stacks, stacks + vals, color=COLORS[c], alpha=0.85)
        stacks += vals
    tsE = tstar[str(EN)]["T_star_K"]
    if tsE:
        ax.axvline(tsE, color="#C23B22", lw=1.2, ls="--")
        ax.annotate(f"T*={tsE:.0f} K", xy=(tsE, 90), fontsize=9, color="#C23B22", ha="right")
    ax.axhline(50, color="k", lw=0.8, ls=":")
    ax.set_xlim(300, 800); ax.set_ylim(0, 100)
    ax.set_xlabel("T (K)"); ax.set_ylabel("份额 (%)" if j == 0 else "")
    ax.set_title(f"E/N = {EN} Td", fontsize=10)
fig.suptitle("NH 生成通道份额随气体温度的演化（纯气相完整机制，1 atm，N$_2$:H$_2$=1:2，n$_e$ 固定）", fontsize=12)
fig.savefig(RES / "fig_i_NH_channel_share_vs_T.png", bbox_inches="tight")
plt.close(fig)
print("saved:", RES / "fig_i_NH_channel_share_vs_T.png")

# ---------------------------------------------------------------- fig_j: VDF
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
TC = {300: "#4C9BD6", 500: "#E8853D", 700: "#C23B22"}
for T in [300, 500, 700]:
    r = results["baseline"][tag_of(T, 45.1)]
    n2tot = float(load_output(BD / f"tscan_output_{tag_of(T,45.1)}.csv")["N2"][-1])
    h2tot = float(load_output(BD / f"tscan_output_{tag_of(T,45.1)}.csv")["H2"][-1])
    axes[0].semilogy(range(1, 9), [max(r["N2_VDF"][f"v{i}"] / n2tot, 1e-30) for i in range(1, 9)],
                     "o-", color=TC[T], label=f"{T} K")
    axes[1].semilogy(range(1, 4), [max(r["H2_VDF"][f"v{i}"] / h2tot, 1e-30) for i in range(1, 4)],
                     "s-", color=TC[T], label=f"{T} K")
axes[0].set_xlabel("N$_2$ 振动态 v"); axes[0].set_ylabel("N$_2$(v)/N$_2$ 占比")
axes[0].set_title("N$_2$ VDF @100 s, 45.1 Td", fontsize=11); axes[0].legend(); axes[0].grid(alpha=0.3)
axes[1].set_xlabel("H$_2$ 振动态 v"); axes[1].set_ylabel("H$_2$(v)/H$_2$ 占比")
axes[1].set_title("H$_2$ VDF @100 s, 45.1 Td", fontsize=11); axes[1].legend(); axes[1].grid(alpha=0.3)
fig.suptitle("振动分布函数随气体温度的变化（Boltzmann 虚线供对照）", fontsize=12)
# Boltzmann 对照（H2: ΔE=5800 K；N2: ~3393 K 等间距近似）
for T in [300, 500, 700]:
    axes[0].semilogy(range(1, 9), [np.exp(-3393 * i / T) for i in range(1, 9)],
                     "--", color=TC[T], alpha=0.35, lw=1)
    axes[1].semilogy(range(1, 4), [np.exp(-5800 * i / T) for i in range(1, 4)],
                     "--", color=TC[T], alpha=0.35, lw=1)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(RES / "fig_j_vdf_vs_T.png", bbox_inches="tight")
plt.close(fig)
print("saved:", RES / "fig_j_vdf_vs_T.png")

# ---------------------------------------------------------------- fig_k: NH3 产率与能耗
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
ENC = {30: "#77B255", 45.1: "#4C9BD6", 80: "#E8853D", 120: "#C23B22"}
for EN in ENS:
    axes[0].semilogy(TS, [results["baseline"][tag_of(T, EN)]["NH3_100s"] for T in TS],
                     "o-", color=ENC[EN], label=f"{EN} Td")
    axes[1].semilogy(TS, [results["baseline"][tag_of(T, EN)]["NH3_per_eV"] for T in TS],
                     "o-", color=ENC[EN], label=f"{EN} Td")
axes[0].set_xlabel("T (K)"); axes[0].set_ylabel("NH$_3$ @100 s (cm$^{-3}$)")
axes[0].set_title("气相 NH$_3$ 产率 vs T", fontsize=11); axes[0].legend(); axes[0].grid(alpha=0.3)
axes[1].set_xlabel("T (K)"); axes[1].set_ylabel("NH$_3$ 分子 / eV")
axes[1].set_title("单位电子能耗产率（恒定功率近似）", fontsize=11); axes[1].legend(); axes[1].grid(alpha=0.3)
fig.suptitle("NH$_3$ 产率与能效随温度（1 atm, N$_2$:H$_2$=1:2, n$_e$=1.17e8 固定）", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(RES / "fig_k_nh3_yield_energy_vs_T.png", bbox_inches="tight")
plt.close(fig)
print("saved:", RES / "fig_k_nh3_yield_energy_vs_T.png")
