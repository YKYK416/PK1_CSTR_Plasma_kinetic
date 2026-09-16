#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pulse_analysis.py - 方案2 方波脉冲扫描分析
读取 build_pulse/pulse_cycles_scan_*.csv 与 t3_r*.csv，
计算每工况：收敛周期数、准稳态 dNH3/cycle、产率速率、能耗、能效、余辉份额
输出 data/results/pulse_scan_results.json
"""
import json, re, sys
from pathlib import Path

BASE = Path(__file__).parent
BUILD = BASE / "build_pulse"
OUT = BASE / "data" / "results" / "pulse_scan_results.json"

EV_J = 1.602e-19
M_NH3 = 17.031 / 6.022e23  # g per molecule

def ffloat(s):
    """容错解析 Fortran ES13.5（三位指数丢 E，如 3.27010-164）"""
    s = s.strip()
    try:
        return float(s)
    except ValueError:
        m = re.match(r'^([+-]?\d*\.?\d+)([+-]\d{2,3})$', s)
        if m:
            return float(m.group(1) + 'e' + m.group(2))
        return float('nan')

def parse_tag(tag):
    m = re.match(r'scan_(decay_)?f(\d+)_d(\d+)_en(\d+)', tag)
    if not m:
        return None
    decay = bool(m.group(1))
    duty = int(m.group(3)) / 10.0   # d01->0.1, d10->1.0 (d04->0.4 etc.)
    return dict(ne_mode='decay' if decay else 'fix',
                freq=float(m.group(2)), duty=duty, EN=float(m.group(4)))

def read_cycles(path):
    rows = []
    lines = path.read_text().splitlines()
    for ln in lines[1:]:
        if not ln.strip():
            continue
        p = ln.split(',')
        if len(p) < 9:
            continue
        rows.append(dict(cycle=int(p[0]), t_end=ffloat(p[1]),
                         nh3_start=ffloat(p[2]), nh3_mid=ffloat(p[3]),
                         nh3_end=ffloat(p[4]), e_on=ffloat(p[5]),
                         e_off=ffloat(p[6]), d_on=ffloat(p[7]), d_off=ffloat(p[8])))
    return rows

def log_converged(tag):
    log = BUILD / "logs" / f"scan_{tag}.log"
    if not log.exists():
        log = BUILD / "logs" / f"{tag}.log"
    if not log.exists():
        return None, None
    txt = log.read_text(errors='replace')
    m = re.search(r'cycles_used =\s*(\d+)\s+converged =\s*(\d+)', txt)
    if m:
        return int(m.group(1)), int(m.group(2))
    return None, None

def summarize(tag):
    meta = parse_tag(tag)
    if meta is None:
        return None
    csv = BUILD / f"pulse_cycles_{tag}.csv"
    if not csv.exists():
        return None
    rows = read_cycles(csv)
    if not rows:
        return None
    last = rows[-1]
    prev = rows[-2] if len(rows) > 1 else None
    cyc_used, conv = log_converged(tag)
    f = meta['freq']
    d_tot = last['d_on'] + last['d_off']
    e_tot = last['e_on'] + last['e_off']
    rate = d_tot * f                       # NH3 molec cm-3 s-1
    power = e_tot * f                      # J cm-3 s-1
    eff_molec_per_J = d_tot / e_tot if e_tot > 0 else float('nan')
    eff_molec_per_eV = eff_molec_per_J * EV_J
    eff_g_per_kWh = eff_molec_per_J * M_NH3 * 3.6e6  # g/kWh
    share_off = last['d_off'] / d_tot if d_tot != 0 else float('nan')
    if prev is not None:
        d_prev = prev['d_on'] + prev['d_off']
        rel_change = abs(d_tot - d_prev) / abs(d_tot) if d_tot != 0 else float('nan')
    else:
        rel_change = float('nan')
    return dict(meta, tag=tag, cycles_used=cyc_used, converged=conv,
                dNH3_cycle=d_tot, dNH3_on=last['d_on'], dNH3_off=last['d_off'],
                dNH3_rel_change_last=rel_change,
                E_on_Jcm3=last['e_on'], E_off_Jcm3=last['e_off'],
                rate_cm3s=rate, power_Wcm3=power,
                eff_molec_per_J=eff_molec_per_J,
                eff_molec_per_eV=eff_molec_per_eV,
                eff_g_per_kWh=eff_g_per_kWh,
                afterglow_share=share_off,
                t_end=last['t_end'], NH3_final=last['nh3_end'])

def main():
    results = {}
    for csv in sorted(BUILD.glob("pulse_cycles_scan_*.csv")):
        tag = csv.stem.replace("pulse_cycles_", "")
        s = summarize(tag)
        if s:
            results[tag] = s
    # Table 3 validation rows
    t3 = {}
    for tag, desc in [("t3_r1", "500Hz duty1.0 EN45.1 ne1x"),
                      ("t3_r2", "1kHz duty1.0 EN49.6 ne2x"),
                      ("t3_r3", "1kHz duty0.4 EN49.6 ne2x")]:
        csv = BUILD / f"pulse_cycles_{tag}.csv"
        rows = read_cycles(csv)
        last = rows[-1]
        cyc, conv = log_converged(tag)
        freq = 500.0 if tag == "t3_r1" else 1000.0
        d_tot = last['d_on'] + last['d_off']
        e_tot = last['e_on'] + last['e_off']
        t3[tag] = dict(desc=desc, freq=freq, cycles_used=cyc, converged=conv,
                       dNH3_cycle=d_tot, rate_cm3s=d_tot * freq,
                       E_cycle_Jcm3=e_tot, power_Wcm3=e_tot * freq,
                       eff_molec_per_J=d_tot / e_tot,
                       afterglow_share=(last['d_off'] / d_tot if d_tot else 0.0))
    t3['ratio_r2_r1'] = t3['t3_r2']['rate_cm3s'] / t3['t3_r1']['rate_cm3s']
    t3['ratio_r3_r1'] = t3['t3_r3']['rate_cm3s'] / t3['t3_r1']['rate_cm3s']
    t3['paper_targets'] = {'r2_r1': 6.4, 'r3_r1': 2.4}

    out = dict(scan=results, table3=t3,
               notes=dict(n_max="f500:20, f1000:24, other:30",
                          tolerances="ATOL=1e2 RTOL=1e-3",
                          EN_off=0.1, Tgas=300, n2_frac=1/3,
                          ne0="1.17e8 cm-3 (fix) or decay tau=3us"))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(f"写入 {OUT}: scan {len(results)} 工况")
    print(f"Table3: r2/r1 = {t3['ratio_r2_r1']:.2f} (论文 6.4) | "
          f"r3/r1 = {t3['ratio_r3_r1']:.2f} (论文 2.4)")
    # 头条统计
    fix = [s for s in results.values() if s['ne_mode'] == 'fix']
    best = max(fix, key=lambda s: s['eff_molec_per_J'])
    print(f"能效最优: {best['tag']} eff={best['eff_molec_per_J']:.3e} molec/J "
          f"({best['eff_g_per_kWh']:.1f} g/kWh), 余辉份额={best['afterglow_share']:.2f}")
    shares = [s['afterglow_share'] for s in fix if s['duty'] < 1.0]
    print(f"余辉份额范围: {min(shares):.3f} - {max(shares):.3f}")
    nc = sum(1 for s in results.values() if s['converged'] == 0)
    print(f"未达收敛判据(n_max截断): {nc}/{len(results)}")
    rc = [s['dNH3_rel_change_last'] for s in results.values()]
    import math
    rc = [x for x in rc if not math.isnan(x)]
    print(f"末周期相对变化: 中位数 {sorted(rc)[len(rc)//2]:.4f}, 最大 {max(rc):.4f}")

if __name__ == "__main__":
    main()
