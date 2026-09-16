#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第三阶段：Hong 完整反应集（48 物种/515 反应）vs 65 反应重构版对比。
读取 build_full/full_output_N2H2_1_2.csv（完整版）
与 data/results/hong_output_N2H2_1_2.csv（65 反应重构版，N2:H2=1:2, 1 atm, 300 K, 45.1 Td），
输出 full_vs_recon_comparison.json 与 fig_e_full_vs_recon.png。"""
import sys, csv, json
from pathlib import Path

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from daimon_runtime import setup_plot
setup_plot()

HERE = Path(__file__).parent.resolve()
RES = HERE / "data" / "results"


def load_csv(path):
    rows = list(csv.reader(open(path, encoding="utf-8")))
    hdr = [h.strip() for h in rows[0]]
    return {n: np.array([float(r[i]) for r in rows[1:]]) for i, n in enumerate(hdr)}


full = load_csv(HERE / "build_full" / "full_output_N2H2_1_2.csv")
recon = load_csv(RES / "hong_output_N2H2_1_2.csv")


def get(d, *names):
    """按候选列名取列（大小写已按文件实际，找不到返回 None）。"""
    for n in names:
        if n in d:
            return d[n]
    return None


def at(d, col, t):
    """最接近 t 时刻的取值（log 时间轴用最近邻即可）。"""
    if col is None:
        return None
    i = int(np.argmin(np.abs(d["time_s"] - t)))
    return float(col[i])


# (标签, 完整版列, 65版列)
SPECIES = [
    ("NH3", ["NH3"], ["NH3"]),
    ("NH", ["NH"], ["NH"]),
    ("NH2", ["NH2"], ["NH2"]),
    ("H", ["H"], ["H"]),
    ("N", ["N"], ["N"]),
    ("N2(v1)", ["N2(V1)"], ["N2(V1)"]),
    ("H2(v1)", ["H2(V1)"], ["H2(V1)"]),
    ("N2(A3)", ["N2(A3)"], ["N2(A3)"]),
    ("N(2D)", ["N(2D)"], ["N(2D)"]),
    ("Te_eV", ["Te_eV"], ["Te_eV"]),
]

snap = {}
for label, cf, cr in SPECIES:
    colf, colr = get(full, *cf), get(recon, *cr)
    snap[label] = {
        f"full_{t}s": at(full, colf, t) for t in (1, 10, 100)
    } | {f"recon_{t}s": at(recon, colr, t) for t in (1, 10, 100)}

# 能量分支（inelastic 功率占比，取 100 s 时刻）
ener = {
    "full": {"Ptot": at(full, get(full, "Ptot"), 100),
             "Pinel": at(full, get(full, "Pinel"), 100),
             "Pelast": at(full, get(full, "Pelast"), 100)},
    "recon": {"Ptot": at(recon, get(recon, "Ptot"), 100),
              "Pinel": at(recon, get(recon, "Pinel"), 100),
              "Pelast": at(recon, get(recon, "Pelast"), 100)},
    "recon_branching_from_make_figures": {"vib_pct": 66.1, "elec_pct": 33.7, "ioniz_pct": 0.07},
}
for k in ("full", "recon"):
    e = ener[k]
    e["Pinel_frac_pct"] = 100.0 * e["Pinel"] / e["Ptot"] if e["Ptot"] else None

out = {
    "condition": "N2:H2=1:2, 1 atm, 300 K, E/N=45.1 Td, ne=1.17e8 cm-3 固定",
    "mechanisms": {
        "full": "515 反应 / 48 物种 (Shao & Mesbah JACS Au 2024 GitHub, 气相=Hong 原始, 表面=Fe+DFT熵势垒)",
        "recon": "65 反应重构版 (Gap1+Gap2 增补后)",
    },
    "species_snapshots_cm-3": snap,
    "energy": ener,
}
with open(RES / "full_vs_recon_comparison.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("saved:", RES / "full_vs_recon_comparison.json")

# ---- 图：2x3 面板 loglog ----
PANELS = [
    ("NH3", "NH3", "NH3"),
    ("NH", "NH", "NH"),
    ("NH2", "NH2", "NH2"),
    ("H", "H", "H"),
    ("N", "N", "N"),
    ("Te_eV", "Te_eV", "Te_eV"),
]
fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=True)
for ax, (title, cf, cr) in zip(axes.ravel(), PANELS):
    colf, colr = get(full, cf), get(recon, cr)
    tf, tr = full["time_s"], recon["time_s"]
    if title == "Te_eV":
        ax.semilogx(tr, colr, lw=1.8, color="tab:blue", label="65 反应重构")
        ax.semilogx(tf, colf, lw=1.8, color="tab:red", ls="--", label="完整 515 反应")
        ax.set_ylabel("T$_e$ (eV)")
    else:
        ax.loglog(tr, np.maximum(colr, 1e-1), lw=1.8, color="tab:blue", label="65 反应重构")
        ax.loglog(tf, np.maximum(colf, 1e-1), lw=1.8, color="tab:red", ls="--",
                  label="完整 515 反应")
        if ax in axes[:, 0]:
            ax.set_ylabel("数密度 (cm$^{-3}$)")
    ax.set_title(title)
    ax.set_xlim(1e-9, 1.5e2)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
for ax in axes[-1]:
    ax.set_xlabel("时间 t (s)")
fig.suptitle("Hong 完整反应集 vs 65 反应重构（N$_2$:H$_2$=1:2, 1 atm, 300 K, 45.1 Td）\n"
             "完整版气相=Hong 原始机制；表面为 Shao & Mesbah 的 Fe 表面+DFT 熵势垒（非论文氧化铝）",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig(RES / "fig_e_full_vs_recon.png", bbox_inches="tight")
plt.close(fig)
print("saved:", RES / "fig_e_full_vs_recon.png")
