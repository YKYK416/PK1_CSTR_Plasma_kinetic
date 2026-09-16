#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_fig_dbuq.py - 方案3 三张图 + swarm 汇总表
fig_o_swarm_validation.png / fig_p_rate_divergence.png / fig_q_chemistry_bands.png
"""
import json, math, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()
import matplotlib.pyplot as plt

BASE = Path(__file__).parent            # Reproduction/2017Hong
RES = BASE / "data" / "results"
sw = json.load(open(RES / "swarm_scan.json"))
ch = json.load(open(RES / "dbuq_chem.json"))

C = dict(A='#d62728', B='#1f77b4', Bhe='#2ca02c', Ahe='#ff7f0e', fit='#555555')

# ============================================================ fig_o
fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6))

ax = axes[0, 0]
for tag, c, lb in [('A_siglo', C['A'], '库A (SIGLO 合并, 65反应版)'),
                   ('B_full', C['B'], '库B (完整机制自带)'),
                   ('B_h2elastic', C['Bhe'], '库B + H2 EFFECTIVE→ELASTIC'),
                   ('A_h2effective', C['Ahe'], '库A + 补 H2 动量转移')]:
    s = sw[tag]
    ax.plot(s['EN'], s['mean_energy_eV'], '-o', ms=4, color=c, label=lb)
ax.set_xscale('log'); ax.set_xlabel('E/N (Td)'); ax.set_ylabel('平均电子能量 (eV)')
ax.set_title('(a) N$_2$:H$_2$=1:2 混合气：平均能量', fontsize=11)
ax.legend(fontsize=8.5); ax.grid(alpha=0.3)

ax = axes[0, 1]
for tag, c, lb in [('A_siglo', C['A'], '库A'), ('B_full', C['B'], '库B'),
                   ('B_h2elastic', C['Bhe'], '库B+ELASTIC'), ('A_h2effective', C['Ahe'], '库A+H2 mt')]:
    s = sw[tag]
    ax.plot(s['EN'], s['drift_m_s'], '-s', ms=4, color=c, label=lb)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('E/N (Td)'); ax.set_ylabel('漂移速度 (m/s)')
ax.set_title('(b) 混合气：漂移速度', fontsize=11)
ax.legend(fontsize=8.5); ax.grid(alpha=0.3, which='both')

ax = axes[1, 0]
en_fit = np.logspace(math.log10(10), math.log10(200), 50)
vd_fit = 10**(5.5236702 + 0.7822439*np.log10(en_fit)) * 1e-2  # cm/s -> m/s
for tag, c, lb in [('A_siglo_N2pure', C['A'], '库A 纯N$_2$'), ('B_full_N2pure', C['B'], '库B 纯N$_2$')]:
    s = sw[tag]
    ax.plot(s['EN'], s['drift_m_s'], '-o', ms=5, color=c, label=lb)
ax.plot(en_fit, vd_fit, '--', color=C['fit'],
        label='Hake–Phelps 实验拟合 (arXiv:2010.07570 引[4])')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('E/N (Td)'); ax.set_ylabel('漂移速度 (m/s)')
ax.set_title('(c) 纯 N$_2$ 漂移速度 vs 经典实验锚点', fontsize=11)
ax.legend(fontsize=8.5); ax.grid(alpha=0.3, which='both')

ax = axes[1, 1]
for tag, c, lb in [('A_siglo_N2pure', C['A'], '库A 纯N$_2$'), ('B_full_N2pure', C['B'], '库B 纯N$_2$'),
                   ('B_full_H2pure', '#9467bd', '库B 纯H$_2$ (EFFECTIVE)'),
                   ('B_h2elastic_H2pure', C['Bhe'], '库B 纯H$_2$ (ELASTIC)')]:
    s = sw[tag]
    ax.plot(s['EN'], s['mean_energy_eV'], '-^', ms=4, color=c, label=lb)
ax.set_xscale('log'); ax.set_xlabel('E/N (Td)'); ax.set_ylabel('平均电子能量 (eV)')
ax.set_title('(d) 纯气体平均能量：H2 动量转移类型效应', fontsize=11)
ax.legend(fontsize=8.5); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(RES / 'fig_o_swarm_validation.png', bbox_inches='tight', dpi=160)
print('saved fig_o')

# ============================================================ fig_p
# 关键过程速率系数比值 A/B（按名匹配）
pairs = [
    ('N2:N2 -> N2(v1)',          'N2:N2 <-> N2(v1)',          'N2 v1 振动'),
    ('N2:N2 -> N2(v4)',          'N2:N2 <-> N2(v4)',          'N2 v4 振动'),
    ('N2:N2 -> N2(A3)',          'N2:N2 -> N2(A3,v0-4)',      'N2(A3) 电子态'),
    ('N2:N2 -> N2(B3)',          'N2:N2 -> N2(B3)',           'N2(B3)'),
    ('N2:N2 -> N2(C3)',          'N2:N2 -> N2(C3)',           'N2(C3)'),
    ('N2:N2 -> N2(SUM)',         'N2:N2 -> N2(SUM)',          'N2 13eV 解离和'),
    ('N2:N2 -> N2^+',            'N2:N2 -> N2^+',             'N2 电离'),
    ('H2:H2 -> H2(v1)',          'H2:H2 <-> H2(v1)',          'H2 v1 振动'),
    ('H2:H2 -> H2(b3)',          'H2:H2 <-> H2(B3SIG)',       'H2(b3) 三重态'),
    ('H2:H2 -> H2^+',            'H2:H2 -> H2^+',             'H2 电离'),
]
fig, ax = plt.subplots(figsize=(9.5, 6))
EN = sw['A_siglo']['EN']
for ka, kb, lb in pairs:
    ra = sw['A_siglo']['rates_m3_s'].get(ka)
    rb = sw['B_full']['rates_m3_s'].get(kb)
    if ra is None or rb is None:
        print('缺过程:', ka, kb); continue
    ratio = [x/y if y > 0 and x > 0 else np.nan for x, y in zip(ra, rb)]
    ax.plot(EN, ratio, '-o', ms=4, label=lb)
ax.axhline(1.0, color='gray', ls='--', lw=1)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('E/N (Td)'); ax.set_ylabel('速率系数比值  库A / 库B')
ax.set_title('关键电子碰撞过程速率系数的库间发散（N$_2$:H$_2$=1:2, 300 K）')
ax.legend(fontsize=9, ncol=2); ax.grid(alpha=0.3, which='both')
fig.tight_layout()
fig.savefig(RES / 'fig_p_rate_divergence.png', bbox_inches='tight', dpi=160)
print('saved fig_p')

# ============================================================ fig_q
fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8))

# (a) NH3 发散带
ax = axes[0]
conds = sorted(ch['bands'].keys(), key=lambda k: (float(k.split('_')[0][1:]), float(k.split('EN')[1])))
labels = {'B_full': '库B 基线', 'B_h2elastic': 'H2 mt 换 ELASTIC',
          'B_vibswapA': 'N2 vib 换 SIGLO', 'B_nov1res': '去 v1res'}
cols = {'B_full': C['B'], 'B_h2elastic': C['Bhe'], 'B_vibswapA': '#9467bd', 'B_nov1res': '#8c564b'}
runs = ch['runs']
xpos = np.arange(len(conds))

def run_tag(var, cond):
    # cond 形如 'T300_EN120.0'；run 键形如 'dbuq_B_full_T300_EN120'（EN 去掉多余的 .0）
    T, ENv = cond.split('_')
    ENnum = ENv[2:]
    if '.' in ENnum:
        ENnum = ENnum.rstrip('0').rstrip('.')
    return f"dbuq_{var}_{T}_EN{ENnum}"

for var in ['B_full', 'B_h2elastic', 'B_vibswapA', 'B_nov1res']:
    ys = []
    for c in conds:
        ys.append(runs.get(run_tag(var, c), {}).get('NH3', np.nan))
    ax.plot(xpos, ys, 'o-', ms=5, color=cols[var], label=labels[var])
ax.set_yscale('log')
ax.set_xticks(xpos)
ax.set_xticklabels([c.replace('_', '\n').replace('EN', '') for c in conds], fontsize=7.5)
ax.set_ylabel('NH$_3$ @100 s (cm$^{-3}$)')
ax.set_title('(a) NH$_3$ 产率的库结构发散带', fontsize=11)
ax.legend(fontsize=8); ax.grid(alpha=0.3, which='both')

# (b) Te 比较
ax = axes[1]
for var in ['B_full', 'B_h2elastic']:
    ys = []
    for c in conds:
        ys.append(runs.get(run_tag(var, c), {}).get('Te_eV', np.nan))
    ax.plot(xpos, ys, 'o-', ms=5, color=cols[var], label=labels[var])
ax.set_xticks(xpos)
ax.set_xticklabels([c.replace('_', '\n').replace('EN', '') for c in conds], fontsize=7.5)
ax.set_ylabel('T$_e$ @100 s (eV)')
ax.set_title('(b) 0D 运行中的 Te：H2 动量转移类型效应', fontsize=11)
ax.legend(fontsize=8.5); ax.grid(alpha=0.3)

# (c) Te 间隙归因瀑布 @43.1 Td（swarm 平均能量）
ax = axes[2]
me = {t: sw[t]['mean_energy_eV'][3] for t in
      ['A_siglo', 'A_h2elastic', 'A_h2effective', 'B_h2elastic', 'B_full', 'B_norot', 'B_nov1res', 'B_a3merged']}
steps = [
    ('库A\n(SIGLO 合并)', me['A_siglo'], C['A']),
    ('+H2 ELASTIC\n(Phelps)', me['A_h2elastic'], C['Bhe']),
    ('ELASTIC→EFFECTIVE\n(类型差异)', me['A_h2effective'], C['Ahe']),
    ('+v1res/转动/A3分支\n(结构差异)', me['B_full'], C['B']),
]
xs = np.arange(len(steps))
bottom = 0
for i, (lb, v, c) in enumerate(steps):
    ax.bar(i, v, 0.62, color=c)
    ax.text(i, v + 0.03, f'{v:.2f}', ha='center', fontsize=10)
ax.set_xticks(xs)
ax.set_xticklabels([s[0] for s in steps], fontsize=7.5)
ax.set_ylabel('平均电子能量 (eV)')
ax.set_title('(c) Te 间隙归因链 @E/N≈43 Td（N2:H2=1:2）', fontsize=11)
ax.grid(alpha=0.3, axis='y')
fig.tight_layout()
fig.savefig(RES / 'fig_q_chemistry_bands.png', bbox_inches='tight', dpi=160)
print('saved fig_q')

# ============================================================ swarm 汇总表 CSV
import csv
with open(RES / 'swarm_table.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['db', 'EN_Td', 'mean_energy_eV', 'Te_proxy_eV(2/3)', 'drift_m_s',
                'mobilityN_1/m/V/s', 'share_elastic', 'share_rot', 'share_vib',
                'share_electronic', 'share_ionization'])
    for tag, s in sorted(sw.items()):
        for j, en in enumerate(s['EN']):
            sh = s['energy_share']
            row = [tag, en, s['mean_energy_eV'][j], s['Te_proxy_eV'][j],
                   s['drift_m_s'][j], s['mobilityN'][j],
                   sh['elastic'][j], sh['rotational'][j], sh['vibrational'][j],
                   sh['electronic'][j], sh['ionization'][j]]
            w.writerow(row)
print('saved swarm_table.csv')

# 归因数字摘要
print('\n=== Te 间隙归因（平均能量 @43.1 Td, 混合气）===')
gap = me['A_siglo'] - me['B_full']
print(f"总间隙: {gap:.3f} eV (A={me['A_siglo']:.3f}, B={me['B_full']:.3f})")
print(f"H2 动量转移缺失(补 Phelps ELASTIC): {me['A_siglo']-me['A_h2elastic']:.3f} eV ({(me['A_siglo']-me['A_h2elastic'])/gap*100:.0f}%)")
print(f"H2 动量转移缺失(补 EFFECTIVE):    {me['A_siglo']-me['A_h2effective']:.3f} eV ({(me['A_siglo']-me['A_h2effective'])/gap*100:.0f}%)")
print(f"EFFECTIVE vs ELASTIC 类型差:       {me['B_h2elastic']-me['B_full']:.3f} eV")
print(f"v1res: {me['B_full']-me['B_nov1res']:.4f} eV | 转动: {me['B_norot']-me['B_full']:.4f} eV | A3分支: {me['B_full']-me['B_a3merged']:.4f} eV")
