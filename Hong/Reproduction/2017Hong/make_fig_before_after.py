#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""2026-08-09 增补：机理增补（Gap1+Gap2）前后对比图。
读取 data/results/hong_output_N2H2_1_2_before_gap12.csv（增补前快照）
与 build/hong_output_N2H2_1_2.csv（增补后），输出 fig_d_before_after_gap12.png。"""
import sys, csv
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


before = load_csv(RES / "hong_output_N2H2_1_2_before_gap12.csv")
after = load_csv(HERE / "build" / "hong_output_N2H2_1_2.csv")

species = ["NH3", "NH", "H", "N"]
fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
for ax, sp in zip(axes.ravel(), species):
    tb, ta = before["time_s"], after["time_s"]
    ax.loglog(tb, np.maximum(before[sp], 1e-1), lw=1.8, color="tab:gray",
              label="增补前 (61 反应)")
    ax.loglog(ta, np.maximum(after[sp], 1e-1), lw=1.8, color="tab:red",
              ls="--", label="增补后 (65 反应)")
    ax.set_title(sp)
    ax.set_xlim(1e-9, 1e2)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
for ax in axes[-1]:
    ax.set_xlabel("时间 t (s)")
for ax in axes[:, 0]:
    ax.set_ylabel("数密度 (cm$^{-3}$)")
fig.suptitle("机理增补前后对比（N$_2$:H$_2$=1:2, 1 atm, 300 K, 45.1 Td）\n"
             "新增：N+H$_2$(v1)/H$_2$(b3)→NH+H；N$_2$(A3)、H$_2$(v1) 增强解离吸附",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(RES / "fig_d_before_after_gap12.png", bbox_inches="tight")
plt.close(fig)
print("saved:", RES / "fig_d_before_after_gap12.png")
