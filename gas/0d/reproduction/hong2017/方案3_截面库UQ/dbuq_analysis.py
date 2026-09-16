#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dbuq_analysis.py - 方案3 0D 化学层发散分析
读取 build_dbuq/tscan_output_dbuq_*.csv 与 tscan_rates_dbuq_*.csv
输出 data/results/dbuq_chem.json
"""
import csv, json, re, sys
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).parent                    # Reproduction/2017Hong
BD = BASE / "build_dbuq"
sys.path.insert(0, str(BASE))
from tscan_analysis import ffloat, parse_sign, classify_nh_prod, NH_CLASSES

VARIANTS = ["B_full", "B_h2elastic", "B_vibswapA", "B_nov1res"]
KEEP = ["NH3", "NH", "N", "H", "N2(V1)", "N2(A3)", "H2(V1)"]

def load_output(path):
    rows = list(csv.reader(open(path, encoding="utf-8")))
    hdr = [h.strip() for h in rows[0]]
    last = rows[-1]
    d = {n: ffloat(last[i]) for i, n in enumerate(hdr)}
    return d

def load_rates_final(path):
    rows = list(csv.reader(open(path, encoding="utf-8")))[1:]
    out = []
    for r in rows:
        if len(r) >= 2:
            sign = r[0].strip().strip('"')
            out.append((sign, ffloat(r[1])))
    return out

def vib_share(rates):
    sh = defaultdict(float); tot = 0.0
    for sign, rate in rates:
        R, P = parse_sign(sign)
        net = P.get("NH", 0) - R.get("NH", 0)
        if net > 0:
            c = classify_nh_prod(R)
            sh[c] += rate * net; tot += rate * net
    return {c: 100.0 * sh[c] / tot if tot else 0.0 for c in NH_CLASSES}

def main():
    runs = {}
    for f in sorted(BD.glob("tscan_output_dbuq_*.csv")):
        tag = f.stem.replace("tscan_output_", "")
        m = re.match(r'dbuq_([A-Za-z0-9_]+)_T(\d+)_EN([\d.]+)', tag)
        if not m: continue
        var, T, EN = m.group(1), float(m.group(2)), float(m.group(3))
        d = load_output(f)
        rec = dict(variant=var, T=T, EN=EN,
                   Te_eV=d["Te_eV"],
                   **{k: d.get(k) for k in KEEP})
        rf = BD / f"tscan_rates_{tag}.csv"
        if rf.exists():
            rec["nh_channel_share_pct"] = vib_share(load_rates_final(rf))
        runs[tag] = rec

    # 发散带: 每 (T,EN) 跨变体 min/max
    bands = {}
    for tag, r in runs.items():
        key = f"T{r['T']:.0f}_EN{r['EN']}"
        bands.setdefault(key, []).append(r)
    band_out = {}
    for key, rs in sorted(bands.items()):
        b = dict(n_variants=len(rs))
        for q in ["Te_eV"] + KEEP:
            vs = [r[q] for r in rs if r.get(q) is not None and r[q] > 0]
            if vs:
                b[q] = dict(min=min(vs), max=max(vs), spread_ratio=max(vs)/min(vs))
        band_out[key] = b

    # 方案1 稳健性: T=800 振动通道份额
    plan1 = {}
    for tag, r in runs.items():
        if r["T"] == 800 and "nh_channel_share_pct" in r:
            plan1[tag] = dict(Te_eV=r["Te_eV"],
                              vib_share_pct=r["nh_channel_share_pct"]["N+H2(v) 振动"],
                              NH3=r["NH3"], NH=r["NH"])

    out = dict(runs=runs, bands=band_out, plan1_robustness=plan1,
               variants_note=dict(B_full="完整机制库原样(基线)",
                                  B_h2elastic="H2 EFFECTIVE→Phelps ELASTIC",
                                  B_vibswapA="N2 v1-v8 截面换 SIGLO 值",
                                  B_nov1res="去掉 N2->N2(v1res)"))
    op = BASE / "data" / "results" / "dbuq_chem.json"
    op.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print("写入", op, f"({len(runs)} runs)")

    print("\n=== 发散带 (spread ratio = max/min 跨4变体) ===")
    print(f"{'cond':16s} {'Te':>6s} {'NH3':>8s} {'NH':>8s} {'N':>8s} {'N2(V1)':>8s}")
    for key, b in sorted(band_out.items()):
        def g(q): return b.get(q, {}).get('spread_ratio', float('nan'))
        print(f"{key:16s} {g('Te_eV'):6.3f} {g('NH3'):8.2f} {g('NH'):8.2f} {g('N'):8.2f} {g('N2(V1)'):8.2f}")

    print("\n=== 方案1 稳健性 (T=800 K 振动通道份额 %) ===")
    for tag, r in sorted(plan1.items()):
        print(f"  {tag:32s} Te={r['Te_eV']:.3f}  vib_share={r['vib_share_pct']:.3f}%  NH3={r['NH3']:.2e}")

if __name__ == "__main__":
    main()
