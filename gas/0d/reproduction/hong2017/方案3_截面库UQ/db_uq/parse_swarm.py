#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parse_swarm.py - 解析 bolsigminus output.dat → swarm_scan.json
每变体: E/N 网格, 平均能量, 迁移率, 功率损失, 逐过程速率系数与能量损失系数
过程名映射: C 索引按物种内文件顺序对应 db 块顺序
"""
import json, re, sys
from pathlib import Path

BASE = Path(__file__).parent
WS = BASE.parents[2]
sys.path.insert(0, str(BASE))
from db_tools import parse, filt, species_of, name_of, DB_A, DB_B, H2_PHELPS

TD = 1e-21  # V m2

def cindex_map(tag):
    """C 索引 -> (species, name, type)；变体构造与 db_tools 完全一致"""
    from db_tools import build_variants
    gas = None
    base = tag
    if tag.endswith('_N2pure'): base, gas = tag[:-7], 'N2'
    elif tag.endswith('_H2pure'): base, gas = tag[:-7], 'H2'
    blocks = build_variants()[base]
    if gas:
        blocks = [b for b in blocks if species_of(b) == gas]
    else:
        # BOLSIG 按脚本物种列表 "N2 H2" 顺序连续编号（已核对 output.dat 头），非文件顺序
        order = {'N2': 0, 'H2': 1}
        blocks = sorted(blocks, key=lambda b: order.get(species_of(b), 9))
    # BOLSIG 编号: 按文件顺序连续编号（已核对 output.dat 过程数 = 块数）
    m = {}
    for i, b in enumerate(blocks, start=1):
        m[f'C{i}'] = dict(species=species_of(b), name=name_of(b), type=b['type'])
    return m

def parse_output(path):
    lines = path.read_text(errors='replace').splitlines()
    n = len(lines)
    i = 0
    out = dict(series={}, rates={}, eloss={})
    cur_c = None
    while i < n:
        ln = lines[i]
        m = re.match(r'^(C\d+)\s+(\S+)\s+(.*)', ln)
        if m and 'Input cross section' not in ln:
            cur_c = m.group(1)
            i += 1
            continue
        if ln.startswith('E/N (Td)'):
            qty = ln.split('\t', 1)[1].strip()
            xs, ys = [], []
            i += 1
            while i < n and re.match(r'^\s*[\d.]+\s', lines[i] or ''):
                p = lines[i].split()
                if len(p) >= 2:
                    try:
                        xs.append(float(p[0])); ys.append(float(p[1]))
                    except ValueError:
                        break
                i += 1
            if qty.startswith('Rate coefficient') and cur_c:
                out['rates'][cur_c] = ys
            elif qty.startswith('Energy loss coefficient') and cur_c:
                out['eloss'][cur_c] = ys
            else:
                out['series'][qty] = ys
                if 'EN' not in out: out['EN'] = xs
            continue
        i += 1
    return out

def group_of(name, typ):
    nl = name.lower()
    if typ in ('ELASTIC', 'EFFECTIVE'): return 'elastic'
    if '(rot)' in nl or '(j0-2)' in nl or '(j1-3)' in nl: return 'rotational'
    if '(v' in nl: return 'vibrational'
    if typ == 'IONIZATION': return 'ionization'
    return 'electronic'

def main():
    result = {}
    for d in sorted((BASE / 'outputs').iterdir()):
        f = d / 'output.dat'
        if not f.exists(): continue
        tag = d.name
        o = parse_output(f)
        if 'EN' not in o: print('EMPTY', tag); continue
        cmap = cindex_map(tag)
        EN = o['EN']
        mean_e = o['series'].get('Mean energy (eV)')
        muN = o['series'].get('Mobility *N (1/m/V/s)')
        pel = o['series'].get('Elastic power loss /N (eV m3/s)')
        pin = o['series'].get('Inelastic power loss /N (eV m3/s)')
        drift = [m * e * TD for m, e in zip(muN, EN)] if muN else None
        # 能量分支: 逐过程能量损失求和分组
        branch = {g: [0.0]*len(EN) for g in ('elastic','rotational','vibrational','electronic','ionization')}
        for c, ys in o['eloss'].items():
            info = cmap.get(c)
            if not info: continue
            g = group_of(info['name'], info['type'])
            for j, y in enumerate(ys):
                if j < len(EN): branch[g][j] += y
        tot = [sum(branch[g][j] for g in branch) for j in range(len(EN))]
        share = {g: [branch[g][j]/tot[j] if tot[j] > 0 else None for j in range(len(EN))] for g in branch}
        rates_named = {}
        for c, ys in o['rates'].items():
            info = cmap.get(c)
            if info: rates_named[f"{info['species']}:{info['name']}"] = ys
        result[tag] = dict(EN=EN, mean_energy_eV=mean_e,
                           Te_proxy_eV=[2/3*x for x in mean_e] if mean_e else None,
                           mobilityN=muN, drift_m_s=drift,
                           power_elastic=pel, power_inelastic=pin,
                           energy_share=share, rates_m3_s=rates_named)
        print(f'{tag:16s} pts={len(EN)} rates={len(rates_named)} '
              f'meanE@43Td={mean_e[3]:.3f} eV' if mean_e else f'{tag} pts={len(EN)}')
    out = BASE.parent / 'data' / 'results' / 'swarm_scan.json'
    out.write_text(json.dumps(result, indent=1))
    print('写入', out)

if __name__ == '__main__':
    main()
