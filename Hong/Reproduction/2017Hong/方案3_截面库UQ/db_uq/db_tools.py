#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
db_tools.py - bolsigdb 解析与变体构建（方案3）
块结构: KEYWORD / 名称行 / 参数+注释行... / 虚线分隔 / 数据行...
"""
import re, sys
from pathlib import Path

BASE = Path(__file__).parent
WS = BASE.parents[2]                      # G:\Kimi_project\ZDPlaskin模拟
DB_A = WS / "Reproduction/2017Hong/data/bolsigdb.dat"       # SIGLO 合并集
DB_B = WS / "Reproduction/2017Hong/build_full/bolsigdb.dat" # 完整机制自带库
H2_PHELPS = WS / "bolsigdb_H2.dat"        # Phelps H2（含 ELASTIC H2）

KEYWORDS = {"ELASTIC","EFFECTIVE","EXCITATION","IONIZATION","ATTACHMENT",
            "ROTATION","VIBRATION","TOTAL","BREMSSTRAHLUNG"}
NUM = re.compile(r'^\s*[+-]?[\d.]+(?:[eE][+-]?\d+)?\s+[+-]?[\d.]+(?:[eE][+-]?\d+)?\s*$')

def parse(path):
    blocks, cur = [], None
    for line in path.read_text(errors='replace').splitlines():
        s = line.strip()
        if s in KEYWORDS:
            if cur: blocks.append(cur)
            cur = dict(type=s, header=[], data=[])
        elif cur is not None:
            if not cur['data'] and set(s) <= set('-') and s:
                continue                      # 分隔虚线
            if not cur['data'] and not NUM.match(line):
                cur['header'].append(line)    # 名称/参数/注释
            elif NUM.match(line):
                cur['data'].append(line)
            else:
                # 数据区后又出现非数值行 -> 块结束异常，忽略
                pass
    if cur: blocks.append(cur)
    return blocks

def species_of(b):
    name = b['header'][0].strip() if b['header'] else ''
    return re.split(r'\s*(?:->|<->)\s*|\s+', name)[0]

def name_of(b):
    return b['header'][0].strip() if b['header'] else ''

def emit(blocks):
    out = []
    for b in blocks:
        out.append(b['type'])
        out.extend(b['header'])
        out.append('-'*60)
        out.extend(b['data'])
        out.append('-'*30)          # 数据结束虚线（原格式必需）
    return '\n'.join(out) + '\n'

def filt(blocks, species=('N2','H2')):
    return [b for b in blocks if species_of(b) in species]

def drop(blocks, pred):
    return [b for b in blocks if not pred(b)]

def build_variants():
    """返回 {tag: blocks} 全部变体（main 与 parse_swarm 共用，保证 C 索引映射一致）"""
    import copy as _copy
    A = filt(parse(DB_A))
    B = filt(parse(DB_B))
    h2p = parse(H2_PHELPS)
    h2_elastic = [b for b in h2p if b['type']=='ELASTIC' and species_of(b)=='H2'][0]
    h2_eff_B   = [b for b in B if b['type']=='EFFECTIVE' and species_of(b)=='H2'][0]
    a3_A       = [b for b in parse(DB_A) if name_of(b)=='N2 -> N2(A3)'][0]

    is_n2vib = lambda b: species_of(b)=='N2' and ('(v' in name_of(b).lower())
    is_h2vib = lambda b: species_of(b)=='H2' and '(v' in name_of(b).lower()
    is_rot   = lambda b: '(rot)' in name_of(b).lower() or '(j0-2)' in name_of(b) or '(j1-3)' in name_of(b)
    is_v1res = lambda b: 'v1res' in name_of(b)
    is_h2mt  = lambda b: species_of(b)=='H2' and b['type'] in ('ELASTIC','EFFECTIVE')

    # B_vibswapA: B 的 N2 v1-v8 数据换成 A(SIGLO) 的同阈振动截面（保留 <-> 与块名，0D 可链接）
    a_vib = {}
    for b in A:
        m = re.match(r'N2 -> N2\(v(\d)\)', name_of(b))
        if m: a_vib[int(m.group(1))] = b['data']
    B_vibswapA = []
    for b in B:
        m = re.match(r'N2 <-> N2\(v(\d)\)', name_of(b))
        if m and int(m.group(1)) in a_vib:
            nb = _copy.deepcopy(b)
            nb['data'] = a_vib[int(m.group(1))]
            B_vibswapA.append(nb)
        else:
            B_vibswapA.append(b)

    return {
        # 基准两套
        'A_siglo':        A,
        'B_full':         B,
        # B 消融/替换（归因：哪些结构特征造成 Te 间隙）
        'B_h2elastic':    [h2_elastic if is_h2mt(b) else b for b in B],      # EFFECTIVE->Phelps ELASTIC
        'B_nov1res':      drop(B, is_v1res),
        'B_norot':        drop(B, is_rot),
        'B_novib':        drop(B, lambda b: is_n2vib(b) or is_h2vib(b)),
        'B_a3merged':     [a3_A if name_of(b)=='N2 -> N2(A3,v0-4)' else b
                           for b in drop(B, lambda b: name_of(b) in
                           ('N2 -> N2(A3,v5-9)','N2 -> N2(A3,v10-)'))],
        # A 反向增补（验证缺 H2 动量转移是否是 A 侧主因）
        'A_h2effective':  A + [h2_eff_B],
        'A_h2elastic':    A + [h2_elastic],
        # 振动截面库间替换（方案1 稳健性检验用；块名/双向性保留，0D 可链接）
        'B_vibswapA':     B_vibswapA,
    }

def main():
    variants = build_variants()
    outdir = BASE / 'dbs'
    for tag, bl in variants.items():
        (outdir / f'{tag}.dat').write_text(emit(bl))
        n2 = sum(1 for b in bl if species_of(b)=='N2')
        h2 = sum(1 for b in bl if species_of(b)=='H2')
        print(f'{tag:16s} N2={n2:2d} H2={h2:2d} 总计={len(bl):2d}')

if __name__ == '__main__':
    main()
