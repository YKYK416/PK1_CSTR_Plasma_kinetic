#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_fig_pulse.py - 方案2 三张图
fig_l_efficiency_map.png   duty×f 能效热力图（4 个 E/N 面板）
fig_m_afterglow_share.png  duty×f 余辉份额热力图（4 个 E/N 面板）
fig_n_table3_validation.png 模型 vs 论文 Table 3（6.4×/2.4×）+ 能效对比
"""
import json, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()
import matplotlib.pyplot as plt

BASE = Path(__file__).parent
FIG = BASE / "figures"
FIG.mkdir(exist_ok=True)
d = json.load(open(BASE / "data" / "results" / "pulse_scan_results.json"))
scan = d['scan']
t3 = d['table3']

freqs = [500, 1000, 5000, 20000, 100000]
dutys = [0.1, 0.2, 0.4, 0.7, 1.0]
ens = [30, 50, 80, 120]
flab = {500: '0.5', 1000: '1', 5000: '5', 20000: '20', 100000: '100'}

def grid(ne_mode, en, key):
    g = np.full((len(dutys), len(freqs)), np.nan)
    for s in scan.values():
        if s['ne_mode'] == ne_mode and s['EN'] == en:
            i = dutys.index(s['duty']); j = freqs.index(s['freq'])
            g[i, j] = s[key]
    return g

def heatmap_fig(key, cmap, title, cbar_label, fname, lognorm=False, vmin=None, vmax=None):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True, sharey=True)
    import matplotlib.colors as mcolors
    for ax, en in zip(axes.flat, ens):
        g = grid('fix', en, key)
        norm = None
        if lognorm:
            gg = np.where(g > 0, g, np.nan)
            norm = mcolors.LogNorm(vmin=np.nanmin(gg), vmax=np.nanmax(gg))
        im = ax.pcolormesh(range(len(freqs) + 1), range(len(dutys) + 1),
                           np.ma.masked_invalid(g), cmap=cmap, norm=norm,
                           vmin=None if lognorm else vmin, vmax=None if lognorm else vmax)
        ax.set_title(f'E/N$_{{on}}$ = {en} Td', fontsize=12)
        # 单元格数值
        for i in range(len(dutys)):
            for j in range(len(freqs)):
                v = g[i, j]
                if not np.isnan(v):
                    if lognorm or vmax is None:
                        txt = f'{v:.0e}' if v >= 100 else f'{v:.1e}'
                    else:
                        txt = f'{v:.2f}'
                    ax.text(j + 0.5, i + 0.5, txt, ha='center', va='center',
                            fontsize=7.5, color='black')
        cb = fig.colorbar(im, ax=ax, shrink=0.85)
        cb.set_label(cbar_label, fontsize=9)
    for ax in axes[1, :]:
        ax.set_xticks(np.arange(len(freqs)) + 0.5)
        ax.set_xticklabels([flab[f] for f in freqs])
        ax.set_xlabel('脉冲频率 f (kHz)', fontsize=11)
    for ax in axes[:, 0]:
        ax.set_yticks(np.arange(len(dutys)) + 0.5)
        ax.set_yticklabels([f'{dd:.1f}' for dd in dutys])
        ax.set_ylabel('占空比 duty', fontsize=11)
    fig.suptitle(title, fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG / fname, bbox_inches='tight', dpi=160)
    print('saved', FIG / fname)

# fig_l: 能效地图 (molec/eV, log)
heatmap_fig('eff_molec_per_eV', 'viridis',
            '方波脉冲 NH$_3$ 合成能效地图（纯气相，ne 固定 1.17×10$^8$ cm$^{-3}$，300 K）',
            '能效 (molec/eV)', 'fig_l_efficiency_map.png', lognorm=True)

# fig_m: 余辉份额地图
heatmap_fig('afterglow_share', 'RdYlBu_r',
            '余辉相 NH$_3$ 生成份额地图（dNH$_3$$_{off}$ / dNH$_3$$_{cycle}$）',
            '余辉份额', 'fig_m_afterglow_share.png', vmin=0, vmax=1)

# fig_n: Table 3 验证
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
ax = axes[0]
x = np.arange(2)
model = [t3['ratio_r2_r1'], t3['ratio_r3_r1']]
paper = [6.4, 2.4]
w = 0.36
ax.bar(x - w/2, paper, w, label='论文 Table 3', color='#888888')
ax.bar(x + w/2, model, w, label='本模型（纯气相）', color='#2b7bbb')
for xi, v in zip(x - w/2, paper):
    ax.text(xi, v + 0.1, f'{v:.1f}×', ha='center', fontsize=11)
for xi, v in zip(x + w/2, model):
    ax.text(xi, v + 0.1, f'{v:.2f}×', ha='center', fontsize=11)
ax.set_xticks(x)
ax.set_xticklabels(['条件2 vs 条件1\n(E/N+10%, n$_e$×2, 连续)',
                    '条件3 vs 条件1\n(E/N+10%, n$_e$×2, 脉冲 duty=0.4)'])
ax.set_ylabel('NH$_3$ 产率比值（准稳态速率）')
ax.set_title('Table 3 验证：产率增强倍数')
ax.legend()
ax.set_ylim(0, 7.6)

ax = axes[1]
labels = ['条件1\n500Hz 连续\n45.1Td, n$_e$',
          '条件2\n1kHz 连续\n49.6Td, 2n$_e$',
          '条件3\n1kHz duty0.4\n49.6Td, 2n$_e$']
rates = [t3['t3_r1']['rate_cm3s'], t3['t3_r2']['rate_cm3s'], t3['t3_r3']['rate_cm3s']]
pows = [t3['t3_r1']['power_Wcm3'], t3['t3_r2']['power_Wcm3'], t3['t3_r3']['power_Wcm3']]
xx = np.arange(3)
ax.bar(xx - w/2, np.array(rates) / rates[0], w, label='NH$_3$ 产率（归一）', color='#2b7bbb')
ax.bar(xx + w/2, np.array(pows) / pows[0], w, label='能耗功率（归一）', color='#e08e45')
for xi, v in zip(xx - w/2, np.array(rates) / rates[0]):
    ax.text(xi, v + 0.08, f'{v:.2f}', ha='center', fontsize=10)
for xi, v in zip(xx + w/2, np.array(pows) / pows[0]):
    ax.text(xi, v + 0.08, f'{v:.2f}', ha='center', fontsize=10)
ax.set_xticks(xx)
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel('相对条件1 归一化')
ax.set_title('脉冲条件3：产率↑ 且 能耗↓')
ax.legend()
fig.tight_layout()
fig.savefig(FIG / 'fig_n_table3_validation.png', bbox_inches='tight', dpi=160)
print('saved', FIG / 'fig_n_table3_validation.png')

# decay 敏感性摘要（报告用）
print('\n=== ne_mode decay vs fix（子网格 27 工况）===')
pairs = []
for s in scan.values():
    if s['ne_mode'] == 'decay':
        tag_fix = s['tag'].replace('scan_decay_', 'scan_')
        sf = scan.get(tag_fix)
        if sf and sf['dNH3_cycle'] > 0:
            pairs.append((s['tag'], s['dNH3_cycle'] / sf['dNH3_cycle'],
                          s['afterglow_share'], sf['afterglow_share']))
ratios = [p[1] for p in pairs]
print(f'配对 {len(pairs)} 组; dNH3(decay)/dNH3(fix): '
      f'min={min(ratios):.3f} med={sorted(ratios)[len(ratios)//2]:.3f} max={max(ratios):.3f}')
dsh = [p[2] - p[3] for p in pairs]
print(f'余辉份额变化(decay-fix): min={min(dsh):.3f} max={max(dsh):.3f}')
for p in sorted(pairs, key=lambda x: x[1])[:3]:
    print('  最低:', p[0], f'ratio={p[1]:.3f}')
for p in sorted(pairs, key=lambda x: -x[1])[:3]:
    print('  最高:', p[0], f'ratio={p[1]:.3f}')
