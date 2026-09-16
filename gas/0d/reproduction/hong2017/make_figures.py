#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Hong 2017 复现后处理：读取 build/ 下 ZDPlasKin 输出，生成对比图到 data/results/"""
import sys, csv, shutil
from pathlib import Path

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from daimon_runtime import setup_plot
setup_plot()

HERE = Path(__file__).parent.resolve()
BUILD = HERE / "build"
RES = HERE / "data" / "results"
RES.mkdir(parents=True, exist_ok=True)


def load_csv(path):
    rows = list(csv.reader(open(path, encoding="utf-8")))
    hdr = [h.strip() for h in rows[0]]
    data = {n: np.array([float(r[i]) for r in rows[1:]]) for i, n in enumerate(hdr)}
    return data


def load_rates(path):
    rows = list(csv.reader(open(path, encoding="utf-8")))[1:]
    return [(r[0].strip('"'), float(r[1])) for r in rows]


# ---------------------------------------------------------------- fig (a)
base = load_csv(BUILD / "hong_output_N2H2_1_2.csv")
t = base["time_s"]
fig, ax = plt.subplots(figsize=(8, 5.5))
show = ["E", "N2(V1)", "N2(A3)", "H2(V1)", "H", "N", "NH", "NH2", "NH3",
        "N2^+", "H3^+", "N2H^+"]
for s in show:
    y = np.maximum(base[s], 1e-1)
    ax.loglog(t, y, label=s, lw=1.6)
ax.set_xlabel("时间 t (s)")
ax.set_ylabel("数密度 (cm$^{-3}$)")
ax.set_title("Hong 2017 基线复现：主要物种密度随时间演化\n"
             "(1 atm, 300 K, N$_2$:H$_2$=1:2, E/N=45.1 Td, $n_e$=1.17e8 cm$^{-3}$)")
ax.legend(ncol=3, fontsize=9)
ax.set_xlim(1e-9, 1e2)
ax.set_ylim(1e2, 1e20)
fig.savefig(RES / "fig_a_species_vs_time.png", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- fig (b)
runs = [("N2H2_1_1", "1:1"), ("N2H2_1_2", "1:2 (基线)"), ("N2H2_1_3", "1:3"),
        ("N2H2_1_4", "1:4"), ("N2H2_3_1", "3:1")]
fig, ax = plt.subplots(figsize=(8, 5.5))
nh3_at_10s = {}
for tag, lab in runs:
    d = load_csv(BUILD / f"hong_output_{tag}.csv")
    ax.loglog(d["time_s"], np.maximum(d["NH3"], 1e-1), lw=1.8, label=f"N$_2$:H$_2$={lab}")
    i10 = np.argmin(np.abs(d["time_s"] - 10.0))
    nh3_at_10s[lab] = d["NH3"][i10]
ax.axhspan(1e14, 1e15, color="gray", alpha=0.15, label="验证目标 1e14–1e15")
ax.set_xlabel("时间 t (s)")
ax.set_ylabel("NH$_3$ 数密度 (cm$^{-3}$)")
ax.set_title("NH$_3$ 密度累积曲线（不同进气比例）")
ax.legend(fontsize=9)
ax.set_xlim(1e-9, 1e2)
fig.savefig(RES / "fig_b_nh3_buildup.png", bbox_inches="tight")
plt.close(fig)

# composition trend at t = 10 s (paper Fig.3 analog)
comp = [("1:1", 0.5), ("1:2", 1/3), ("1:3", 0.25), ("1:4", 0.2), ("3:1", 0.75)]
fig, ax = plt.subplots(figsize=(6.5, 4.5))
xs = [c[1] for c in comp]
labs = ["1:1", "1:2", "1:3", "1:4", "3:1"]
ys = [nh3_at_10s.get("1:1"), nh3_at_10s.get("1:2 (基线)"), nh3_at_10s.get("1:3"),
      nh3_at_10s.get("1:4"), nh3_at_10s.get("3:1")]
ax.plot(range(5), ys, "o-", lw=2, ms=8, color="tab:red")
ax.set_xticks(range(5)); ax.set_xticklabels(labs)
ax.set_yscale("log")
ax.set_xlabel("进气比例 N$_2$:H$_2$")
ax.set_ylabel("NH$_3$ 密度 @ t=10 s (cm$^{-3}$)")
ax.set_title("NH$_3$ 随气体组成变化（对应论文 Fig.3 趋势）")
for i, y in enumerate(ys):
    ax.annotate(f"{y:.1e}", (i, y), textcoords="offset points",
                xytext=(0, 8), ha="center", fontsize=9)
fig.savefig(RES / "fig_b2_nh3_vs_composition.png", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- fig (c)
rates = dict(load_rates(BUILD / "hong_rates_N2H2_1_2.csv"))
# rate [cm-3 s-1] x threshold [eV] -> energy channeling [eV cm-3 s-1]
chan = {}
def add(group, sign, thr):
    for k, v in rates.items():
        if sign.lower() in k.lower():
            chan[group] = chan.get(group, 0.0) + abs(v) * thr
            return
add("振动激发", "N2->N2(V1)", 0.29)
add("振动激发", "H2->H2(V1)", 0.516)
add("电子激发", "N2->N2(A3)", 6.17)
add("电子激发", "N2->N2(B3)", 7.35)
add("电子激发", "N2->N2(A1)", 8.55)
add("电子激发", "H2->H2(B3)", 8.9)
add("电子激发", "H2->H2(B1)", 11.3)
add("解离", "E+N2=>E+N+N", 9.75)
add("电离", "N2->N2^+", 15.6)
add("电离", "H2->H2^+", 15.4)
add("电离", "N->N^+", 14.5)
add("电离", "E+H=>E+E+H^+", 13.6)
add("电离", "E+NH3=>E+E+NH3^+", 10.2)

# elastic power from Ptot/Pelast columns of output (reduced powers, eV cm3 s-1)
pel = abs(base["Pelast"][-1]); pin = abs(base["Pinel"][-1])
inel_sum = sum(chan.values())
scale = pin / inel_sum if inel_sum > 0 else 1.0  # normalize to BOLSIG inelastic power

labels = ["弹性(动量转移)"] + list(chan.keys())
vals = [pel] + [chan[k] * scale for k in list(chan.keys())]
tot = sum(vals)
fig, ax = plt.subplots(figsize=(7, 4.5))
colors = ["#888888", "#2ca02c", "#1f77b4", "#ff7f0e", "#d62728"]
bars = ax.bar(labels, [v / tot * 100 for v in vals], color=colors)
for b, v in zip(bars, vals):
    ax.annotate(f"{v/tot*100:.2f}%", (b.get_x() + b.get_width()/2, b.get_height()),
                ha="center", va="bottom", fontsize=10)
ax.set_ylabel("电子能量沉积份额 (%)")
ax.set_title("电子能量通道分配 @ E/N=45.1 Td, t=100 s\n（振动激发主导——与论文“低电子能量”定性结论一致）")
fig.savefig(RES / "fig_c_energy_branching.png", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- summary
summary = {
    "baseline_N2H2_1_2": {
        "t_end_s": float(t[-1]),
        "Te_eV": float(base["Te_eV"][-1]),
        "NH3_final": float(base["NH3"][-1]),
        "NH3_at_1s": float(base["NH3"][np.argmin(np.abs(t-1.0))]),
        "NH3_at_10s": float(base["NH3"][np.argmin(np.abs(t-10.0))]),
        "N2v1_final": float(base["N2(V1)"][-1]),
        "N2A3_final": float(base["N2(A3)"][-1]),
        "H_final": float(base["H"][-1]),
        "N_final": float(base["N"][-1]),
    },
    "NH3_at_10s_vs_composition": {k: float(v) for k, v in
        zip(["1:1", "1:2", "1:3", "1:4", "3:1"], ys)},
    "energy_branching_percent": {l: float(v/tot*100) for l, v in zip(labels, vals)},
}
import json
with open(RES / "summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, ensure_ascii=False, indent=2)

# copy raw CSVs
for tag, _ in runs:
    shutil.copy(BUILD / f"hong_output_{tag}.csv", RES / f"hong_output_{tag}.csv")
    shutil.copy(BUILD / f"hong_rates_{tag}.csv", RES / f"hong_rates_{tag}.csv")

print(json.dumps(summary, ensure_ascii=False, indent=2))
print("Figures + CSVs saved to:", RES)
