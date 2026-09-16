#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Hong 2017 复现 · 反应路径通量分析
输入: build/hong_rates_t_N2H2_1_2.csv        (65 反应重构版, 多时刻速率)
      build_full/full_rates_t_N2H2_1_2.csv   (完整 515 反应版, 多时刻速率)
      以及两套 *_output_N2H2_1_2.csv (取 t=100s 表面位占据)
输出: data/results/pathway_analysis.json
      data/results/fig_f_nh3_channels_time.png
      data/results/fig_g_pathway_flux.png
      data/results/fig_h_channel_contrib.png
所有数字均来自 CSV，不虚构。"""
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
RES = HERE / "data" / "results"
TIMES = [1e-6, 1e-4, 1e-2, 1e-1, 1.0, 10.0, 100.0]

# ---------------------------------------------------------------- 反应式解析
def norm_name(s):
    s = s.strip()
    # bolsig 标签中的振动细分: N2(A3,V0-4) -> N2(A3)
    s = re.sub(r",V\d[^)]*\)", ")", s)
    return s

def split_side(side):
    """按 '+' 分割但保护 '^+' 离子后缀；剥离前置化学计量系数。"""
    raw = side.split("+")
    toks, i = [], 0
    while i < len(raw):
        t = raw[i]
        if t.endswith("^") and i + 1 < len(raw):
            t = t + "+" + raw[i + 1]
            i += 1
        toks.append(t)
        i += 1
    out = {}
    for t in toks:
        t = norm_name(t)
        if not t:
            continue
        m = re.match(r"^(\d+)([A-Za-z].*)$", t)
        if m:
            n, name = int(m.group(1)), m.group(2)
        else:
            n, name = 1, t
        out[name] = out.get(name, 0) + n
    return out

def parse_sign(sign):
    if sign.startswith("bolsig:"):
        body = sign[len("bolsig:"):]
        # 特殊集总标签：N2(SUM) 是 LXCat 集总解离截面 e+N2=>e+2N
        if body == "N2->N2(SUM)":
            return {"N2": 1}, {"N": 2}
        a, b = body.split("->")
        return {norm_name(a): 1}, {norm_name(b): 1}
    lhs, rhs = sign.split("=>")
    return split_side(lhs), split_side(rhs)

# ---------------------------------------------------------------- 数据加载
def load_rates_t(path):
    """-> dict[float, list[(sign, R, P, rate)]]；重复 sign 按时间合并求和。"""
    rows = list(csv.reader(open(path, encoding="utf-8")))[1:]
    agg = defaultdict(lambda: defaultdict(float))
    for t, sign, rate in rows:
        agg[float(t)][sign] += float(rate)
    out = {}
    for t, d in agg.items():
        lst = []
        for sign, rate in d.items():
            R, P = parse_sign(sign)
            lst.append((sign, R, P, rate))
        out[t] = lst
    return out

MECHS = {
    "recon": {
        "label": "65 反应重构",
        "rates": HERE / "build" / "hong_rates_t_N2H2_1_2.csv",
        "output": HERE / "build" / "hong_output_N2H2_1_2.csv",
        "surface": {"S", "HS", "NS", "NHS", "NH2S"},
        "S": "S", "NS": "NS", "HS": "HS", "NHS": "NHS", "NH2S": "NH2S",
        "h2v": {"H2(V1)"},
        "h2el": {"H2(B3)", "H2(B1)"},
        "n_atm": {"N"},
        "n2_states": lambda s: s == "N2" or s.startswith("N2("),
        "h2_states": lambda s: s == "H2" or s.startswith("H2("),
    },
    "full": {
        "label": "完整 515 反应",
        "rates": HERE / "build_full" / "full_rates_t_N2H2_1_2.csv",
        "output": HERE / "build_full" / "full_output_N2H2_1_2.csv",
        "surface": {"SURF", "HSURF", "NSURF", "NHSURF", "NH2SURF"},
        "S": "SURF", "NS": "NSURF", "HS": "HSURF", "NHS": "NHSURF", "NH2S": "NH2SURF",
        "h2v": {"H2(V1)", "H2(V2)", "H2(V3)", "H2(RYDBERG_SUM)"},
        "h2el": {"H2(B3SIG)", "H2(B1SIG)", "H2(C3PI)", "H2(A3SIG)"},
        "n_atm": {"N", "N(2D)", "N(2P)"},
        "n2_states": lambda s: s == "N2" or s.startswith("N2("),
        "h2_states": lambda s: s == "H2" or s.startswith("H2("),
    },
}

# ---------------------------------------------------------------- 通道分类
def classify_nh3_prod(R, P, mech):
    surf = mech["surface"]
    if any(s in surf for s in list(R) + list(P)):
        if set(R) <= surf:
            return "表面 LH: NHx(s)+H(s)"
        return "表面 ER: 气相种+吸附种"
    if sum(R.values()) >= 3:
        return "气相三体"
    return "气相两体"

def classify_ns_prod(R, mech):
    if any(mech["n2_states"](s) for s in R):
        return "N2/N2* 解离吸附"
    if any(s in mech["n_atm"] for s in R):
        return "原子 N 直接吸附"
    return "其他"

def classify_hs_prod(R, mech):
    if any(mech["h2_states"](s) for s in R):
        return "H2/H2* 解离吸附"
    if "H" in R:
        return "原子 H 直接吸附"
    return "其他"

def classify_nh_prod(R, mech):
    if any(s in mech["n_atm"] for s in R):
        if any(s in mech["h2v"] for s in R):
            return "N + H2(v)"
        if any(s in mech["h2el"] for s in R):
            return "N + H2(电子态)"
        if "H2" in R:
            return "N + H2 基态"
    if "NH3" in R:
        return "NH3 解离"
    if "NH2" in R:
        return "NH2 转化"
    return "其他"

# ---------------------------------------------------------------- 分析核心
def prod_cons(species, rxns):
    prod, cons = [], []
    for sign, R, P, rate in rxns:
        net = P.get(species, 0) - R.get(species, 0)
        if net > 0:
            prod.append((sign, rate * net))
        elif net < 0:
            cons.append((sign, -rate * net))
    prod.sort(key=lambda x: -x[1])
    cons.sort(key=lambda x: -x[1])
    return prod, cons

def top5(pairs):
    tot = sum(v for _, v in pairs)
    return [[s, v, (100.0 * v / tot if tot > 0 else 0.0)] for s, v in pairs[:5]], tot

def class_shares(species, rxns, classifier, mech):
    shares = defaultdict(float)
    tot = 0.0
    for sign, R, P, rate in rxns:
        net = P.get(species, 0) - R.get(species, 0)
        if net > 0:
            c = classifier(R, mech) if classifier != classify_nh3_prod else classifier(R, P, mech)
            shares[c] += rate * net
            tot += rate * net
    return {k: (100.0 * v / tot if tot > 0 else 0.0) for k, v in shares.items()}, shares, tot

def integrate_powerlaw(ts, rs):
    """分段幂律精确积分（对数-对数线性插值的解析积分）。"""
    tot = 0.0
    for i in range(len(ts) - 1):
        t1, t2, r1, r2 = ts[i], ts[i + 1], rs[i], rs[i + 1]
        r1c, r2c = max(r1, 1e-300), max(r2, 1e-300)
        p = np.log(r2c / r1c) / np.log(t2 / t1)
        if abs(p + 1.0) < 1e-8:
            tot += r1c * t1 * np.log(t2 / t1)
        else:
            tot += r1c * t1 * ((t2 / t1) ** (p + 1.0) - 1.0) / (p + 1.0)
    return tot

def read_densities_at(path, t_target):
    rows = list(csv.reader(open(path, encoding="utf-8")))
    hdr = [h.strip() for h in rows[0]]
    ti = hdr.index("time_s")
    best = min(rows[1:], key=lambda r: abs(float(r[ti]) - t_target))
    return {hdr[i]: float(best[i]) for i in range(len(hdr))}

# ---------------------------------------------------------------- 主分析
result = {}
rates_data, dens_data = {}, {}
for key, cfg in MECHS.items():
    rates_data[key] = load_rates_t(cfg["rates"])
    dens_data[key] = read_densities_at(cfg["output"], 100.0)

SPECIES_OF = lambda cfg: [
    "NH3", "NH", "NH2", "N", "H", "N2(V1)", "H2(V1)",
    cfg["NS"], cfg["HS"], cfg["NHS"], cfg["NH2S"], cfg["S"],
]

for key, cfg in MECHS.items():
    surf = cfg["surface"]
    me = {"times": TIMES, "species": {}, "nh3_class_shares": {},
          "ns_source_split": {}, "hs_source_split": {}, "nh_source_split": {}}
    for t in TIMES:
        rxns = rates_data[key][t]
        tk = f"{t:g}"
        # 各物种 top5 产/耗
        for sp in SPECIES_OF(cfg):
            if sp not in me["species"]:
                me["species"][sp] = {"prod": {}, "cons": {}}
            prod, cons = prod_cons(sp, rxns)
            p5, ptot = top5(prod)
            c5, ctot = top5(cons)
            me["species"][sp]["prod"][tk] = {"top5": p5, "total": ptot}
            me["species"][sp]["cons"][tk] = {"top5": c5, "total": ctot}
        # NH3 通道分类占比
        sh, _, _ = class_shares("NH3", rxns, classify_nh3_prod, cfg)
        me["nh3_class_shares"][tk] = sh
        # N(s)/H(s)/NH 来源分解
        sh, _, _ = class_shares(cfg["NS"], rxns, classify_ns_prod, cfg)
        me["ns_source_split"][tk] = sh
        sh, _, _ = class_shares(cfg["HS"], rxns, classify_hs_prod, cfg)
        me["hs_source_split"][tk] = sh
        sh, _, _ = class_shares("NH", rxns, classify_nh_prod, cfg)
        me["nh_source_split"][tk] = sh
    # ---- 时间积分: NH3 各生成通道累计产量 (1e-6 s -> 100 s) ----
    nh3_signs = set()
    for t in TIMES:
        for sign, R, P, rate in rates_data[key][t]:
            if P.get("NH3", 0) - R.get("NH3", 0) > 0:
                nh3_signs.add(sign)
    integ = []
    for sign in nh3_signs:
        R, P = parse_sign(sign)
        net = P["NH3"] - R.get("NH3", 0)
        rs = []
        for t in TIMES:
            r = next((rate for s2, _, _, rate in rates_data[key][t] if s2 == sign), 0.0)
            rs.append(r * net)
        integ.append((sign, integrate_powerlaw(TIMES, rs), classify_nh3_prod(R, P, cfg)))
    integ.sort(key=lambda x: -x[1])
    tot_int = sum(v for _, v, _ in integ)
    me["nh3_integrated_1e-6_100s"] = {
        "total": tot_int,
        "top_channels": [[s, v, 100.0 * v / tot_int if tot_int else 0.0, c]
                         for s, v, c in integ[:10]],
        "by_class": {c: 100.0 * sum(v for _, v, cc in integ if cc == c) / tot_int
                     for c in {cc for _, _, cc in integ}} if tot_int else {},
    }
    # 表面位占据 @100s
    d = dens_data[key]
    me["site_densities_100s"] = {sp: d.get(sp) for sp in
                                 [cfg["S"], cfg["HS"], cfg["NS"], cfg["NHS"], cfg["NH2S"]]}
    result[key] = me

with open(RES / "pathway_analysis.json", "w", encoding="utf-8") as f:
    json.dump(result, f, ensure_ascii=False, indent=2)
print("saved:", RES / "pathway_analysis.json")

# ---------------------------------------------------------------- 打印摘要（供报告）
def fmt(v):
    return f"{v:.2e}" if v is not None else "-"

for key, cfg in MECHS.items():
    me = result[key]
    print("\n" + "=" * 78)
    print(f"【{cfg['label']}】")
    for tk in ["1e-06", "0.01", "100"]:
        print(f"\n-- t = {tk} s  NH3 生成 top5 --")
        for s, v, pct in me["species"]["NH3"]["prod"][tk]["top5"]:
            print(f"   {pct:6.2f}%  {fmt(v):>10s}  {s}")
        print(f"   NH3 生成通道分类占比: "
              + ", ".join(f"{k}={v:.1f}%" for k, v in sorted(me["nh3_class_shares"][tk].items(), key=lambda x: -x[1])))
        print(f"   NH 来源: " + ", ".join(f"{k}={v:.1f}%" for k, v in sorted(me["nh_source_split"][tk].items(), key=lambda x: -x[1])))
        print(f"   N(s) 来源: " + ", ".join(f"{k}={v:.1f}%" for k, v in sorted(me["ns_source_split"][tk].items(), key=lambda x: -x[1])))
        print(f"   H(s) 来源: " + ", ".join(f"{k}={v:.1f}%" for k, v in sorted(me["hs_source_split"][tk].items(), key=lambda x: -x[1])))
    print("\n-- N 原子生成 top5 @100s --")
    for s, v, pct in me["species"]["N"]["prod"]["100"]["top5"]:
        print(f"   {pct:6.2f}%  {fmt(v):>10s}  {s}")
    print("-- NH 生成 top5 @100s --")
    for s, v, pct in me["species"]["NH"]["prod"]["100"]["top5"]:
        print(f"   {pct:6.2f}%  {fmt(v):>10s}  {s}")
    print("\n-- NH3 累计生成 (1e-6→100 s, 分段幂律积分) top5 --")
    for s, v, pct, c in me["nh3_integrated_1e-6_100s"]["top_channels"][:5]:
        print(f"   {pct:6.2f}%  {fmt(v):>10s}  [{c}]  {s}")
    print("   按类别: " + ", ".join(f"{k}={v:.1f}%" for k, v in sorted(me["nh3_integrated_1e-6_100s"]["by_class"].items(), key=lambda x: -x[1])))
    print("   表面位密度 @100s: " + ", ".join(f"{k}={fmt(v)}" for k, v in me["site_densities_100s"].items()))

# ---------------------------------------------------------------- fig_f: NH3 生成通道速率随时间
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=False)
for ax, (key, cfg) in zip(axes, MECHS.items()):
    me = result[key]
    # 早期(1e-4)与晚期(100s) top5 的并集
    signs = []
    for tk in ["0.0001", "100"]:
        for s, _, _ in me["species"]["NH3"]["prod"][tk]["top5"]:
            if s not in signs:
                signs.append(s)
    for s in signs:
        R, P = parse_sign(s)
        net = P["NH3"] - R.get("NH3", 0)
        rs = []
        for t in TIMES:
            r = next((rate for s2, _, _, rate in rates_data[key][t] if s2 == s), 0.0)
            rs.append(max(r * net, 1e-30))
        ax.loglog(TIMES, rs, "o-", lw=1.6, ms=4, label=s)
    ax.set_title(f"{cfg['label']}", fontsize=11)
    ax.set_xlabel("时间 t (s)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7, loc="best")
axes[0].set_ylabel("NH$_3$ 生成速率 (cm$^{-3}$ s$^{-1}$)")
fig.suptitle("NH$_3$ 主要生成通道速率随时间演化（N$_2$:H$_2$=1:2, 1 atm, 300 K, 45.1 Td）", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(RES / "fig_f_nh3_channels_time.png", bbox_inches="tight")
plt.close(fig)
print("saved:", RES / "fig_f_nh3_channels_time.png")

# ---------------------------------------------------------------- fig_h: NH3 生成通道分类占比堆叠条
CLASSES = ["气相三体", "气相两体", "表面 ER: 气相种+吸附种", "表面 LH: NHx(s)+H(s)"]
COLORS = {"气相三体": "#4C9BD6", "气相两体": "#9ECAE1",
          "表面 ER: 气相种+吸附种": "#E8853D", "表面 LH: NHx(s)+H(s)": "#C23B22"}
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
for ax, (key, cfg) in zip(axes, MECHS.items()):
    me = result[key]
    bottoms = np.zeros(len(TIMES))
    for c in CLASSES:
        vals = np.array([me["nh3_class_shares"][f"{t:g}"].get(c, 0.0) for t in TIMES])
        ax.bar(range(len(TIMES)), vals, bottom=bottoms, color=COLORS[c], label=c, width=0.65)
        bottoms += vals
    ax.set_xticks(range(len(TIMES)))
    ax.set_xticklabels([f"{t:g}" for t in TIMES], fontsize=9)
    ax.set_xlabel("时间 t (s)")
    ax.set_title(cfg["label"], fontsize=11)
    ax.set_ylim(0, 100)
axes[0].set_ylabel("NH$_3$ 生成占比 (%)")
axes[1].legend(fontsize=8, loc="lower right")
fig.suptitle("NH$_3$ 生成通道分类占比随时间的切换（早期气相三体 → 晚期表面通道）", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(RES / "fig_h_channel_contrib.png", bbox_inches="tight")
plt.close(fig)
print("saved:", RES / "fig_h_channel_contrib.png")

# ---------------------------------------------------------------- fig_g: 晚期通量图（t=100 s）
def class_rate(key, species, classifier, t, classes):
    """t 时刻 species 生成中属于 classes 集合的速率总和。"""
    cfg, rxns = MECHS[key], rates_data[key][t]
    tot = 0.0
    for sign, R, P, rate in rxns:
        net = P.get(species, 0) - R.get(species, 0)
        if net > 0:
            c = classifier(R, P, cfg) if classifier == classify_nh3_prod else classifier(R, cfg)
            if c in classes:
                tot += rate * net
    return tot

def sign_rate(key, signs_want, t):
    rxns = rates_data[key][t]
    return sum(rate for sign, _, _, rate in rxns if sign in signs_want)

def draw_panel(ax, key, t=100.0):
    cfg = MECHS[key]
    ax.set_xlim(0, 10); ax.set_ylim(0, 7.2); ax.axis("off")
    def box(x, y, w, h, text, fc="#EAF2FB"):
        ax.add_patch(plt.Rectangle((x, y), w, h, fc=fc, ec="#335", lw=1.2, zorder=2))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8.5, zorder=3)
    def arrow(x1, y1, x2, y2):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="-|>", color="#444", lw=1.6), zorder=1)
    def lab(x, y, text, ha="center", color="#222", fs=7):
        ax.text(x, y, text, ha=ha, va="center", fontsize=fs, color=color, zorder=4)
    # ---- 速率 ----
    r_en2 = sign_rate(key, {"E+N2=>E+N+N", "bolsig:N2->N2(SUM)"}, t)
    r_nex = class_rate(key, "NH", classify_nh_prod, t, {"N + H2(v)", "N + H2(电子态)"})
    r_hnh2 = sign_rate(key, {"H+NH2=>H2+NH"}, t)
    r_ns_diss = class_rate(key, cfg["NS"], classify_ns_prod, t, {"N2/N2* 解离吸附"})
    r_ns_atom = class_rate(key, cfg["NS"], classify_ns_prod, t, {"原子 N 直接吸附"})
    r_hs_diss = class_rate(key, cfg["HS"], classify_hs_prod, t, {"H2/H2* 解离吸附"})
    r_hs_atom = class_rate(key, cfg["HS"], classify_hs_prod, t, {"原子 H 直接吸附"})
    r_lh = class_rate(key, "NH3", classify_nh3_prod, t, {"表面 LH: NHx(s)+H(s)"})
    r_er = class_rate(key, "NH3", classify_nh3_prod, t, {"表面 ER: 气相种+吸附种"})
    r_g3 = class_rate(key, "NH3", classify_nh3_prod, t, {"气相三体"})
    r_g2 = class_rate(key, "NH3", classify_nh3_prod, t, {"气相两体"})
    r_nh3d = sign_rate(key, {"E+NH3=>E+NH2+H", "E+NH3=>E+NH+H2", "H+NH3=>NH2+H2"}, t)
    # ---- boxes ----
    box(0.1, 5.8, 1.6, 0.9, "N$_2$ + H$_2$\n1 atm, 300 K")
    box(2.4, 5.8, 1.7, 0.9, "e$^-$ 碰撞\n解离·激发")
    box(4.8, 5.8, 2.0, 0.9, "N, H 原子\nN$_2$(v), H$_2$(v)")
    box(7.7, 5.9, 1.3, 0.8, "NH")
    box(7.7, 4.2, 1.3, 0.8, "NH$_2$")
    box(2.4, 0.7, 1.4, 0.8, "H(s)", fc="#FDECEA")
    box(5.0, 0.7, 1.4, 0.8, "N(s)", fc="#FDECEA")
    box(7.7, 0.7, 1.5, 0.8, "NH$_2$(s)", fc="#FDECEA")
    box(8.6, 2.6, 1.3, 0.9, "NH$_3$", fc="#E8F6E8")
    # ---- arrows ----
    arrow(1.7, 6.25, 2.4, 6.25)
    arrow(4.1, 6.25, 4.8, 6.25); lab(4.45, 6.6, f"e+N$_2$→2N\n{fmt(r_en2)}")
    arrow(6.8, 6.3, 7.7, 6.3); lab(7.25, 6.65, f"N+H$_2$*→NH\n{fmt(r_nex)}", color="#C23B22")
    arrow(8.0, 5.0, 8.0, 5.9); lab(7.9, 5.45, f"H+NH$_2$→NH\n{fmt(r_hnh2)}", ha="right")
    arrow(3.1, 5.8, 3.1, 1.5); lab(3.25, 3.5, f"H$_2$*解离吸附 {fmt(r_hs_diss)}\n原子H吸附 {fmt(r_hs_atom)}", ha="left")
    arrow(5.7, 5.8, 5.7, 1.5); lab(5.85, 3.5, f"N$_2$*解离吸附 {fmt(r_ns_diss)}\n原子N吸附 {fmt(r_ns_atom)}", ha="left")
    arrow(6.4, 1.1, 7.7, 1.1); lab(7.05, 0.78, "+H(s) 逐级加氢", fs=6.5)
    arrow(8.5, 1.5, 9.1, 2.6); lab(8.42, 2.2, f"LH {fmt(r_lh)}", ha="right", color="#C23B22")
    arrow(8.6, 4.6, 9.0, 4.1); lab(9.1, 4.5, f"气相 {fmt(r_g3 + r_g2)}", ha="left", color="#4C9BD6")
    # ---- 底部清单 ----
    lab(0.15, 0.32, f"NH$_3$ 生成: 气相三体 {fmt(r_g3)} | 气相两体 {fmt(r_g2)} | "
                    f"表面ER {fmt(r_er)} | 表面LH {fmt(r_lh)}", ha="left", fs=8)
    lab(0.15, 0.02, f"NH$_3$ 解离(→NH$_2$/NH, 含 H+NH$_3$): {fmt(r_nh3d)}", ha="left", fs=8)
    ax.set_title(f"{cfg['label']} · t=100 s 主要通量 (cm$^{{-3}}$ s$^{{-1}}$)", fontsize=10)

fig, axes = plt.subplots(1, 2, figsize=(15, 5.6))
draw_panel(axes[0], "full")
draw_panel(axes[1], "recon")
fig.tight_layout()
fig.savefig(RES / "fig_g_pathway_flux.png", bbox_inches="tight")
plt.close(fig)
print("saved:", RES / "fig_g_pathway_flux.png")
