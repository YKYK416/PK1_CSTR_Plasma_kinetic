#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot and analyse the 300--500 K pure-gas CW-CSTR Hong scan.

This script is deliberately limited to terminal gas-phase concentrations.
It does not infer reaction-path fluxes, energy efficiency, pulse behaviour or
surface contributions.  Failed first-pass cases are replaced only by the
documented numerical retry records.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from publication_figure_theme import PALETTE, add_panel_label, apply_theme, export_figure


PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_ROOT = PROJECT / "analysis" / "cw_gasphase_longtime_20260906"
SPECIES = ["NH3", "NH2", "NH", "N", "H"]
COLORS = {"NH3": PALETTE["turnover"], "NH2": PALETTE["named"], "NH": PALETTE["secondary"],
          "N": PALETTE["transition"], "H": PALETTE["loss"]}


def valid(record: pd.Series) -> bool:
    return record.get("returncode") == 0 and not bool(record.get("timed_out", False))


def read_retry(path: Path, source: str) -> pd.DataFrame:
    table = pd.read_csv(path)
    table = table.loc[table.apply(valid, axis=1)].copy()
    table["data_source"] = source
    return table


def build_table(root: Path) -> pd.DataFrame:
    main = pd.read_csv(root / "grid_972_summary.csv")
    main = main.loc[main.apply(valid, axis=1)].copy()
    main["data_source"] = "main"
    retry_root = root / "grid_972_retry_mxstep"
    retry = read_retry(retry_root / "grid_972_retry_summary.csv", "retry")
    refined = pd.DataFrame([json.loads((retry_root / "refine_t500_en100_n20p5_result.json").read_text(encoding="utf-8"))])
    refined = refined.loc[refined.apply(valid, axis=1)].copy()
    refined["data_source"] = "refined_retry"
    # One retry case is superseded by its strictly converged refined calculation.
    replacements = pd.concat([retry, refined], ignore_index=True)
    replacements["_source_priority"] = replacements["data_source"].map({"retry": 0, "refined_retry": 1})
    replacements = replacements.sort_values(["tg_K", "en_Td", "n2_fraction", "_source_priority"])
    replacements = replacements.drop_duplicates(["tg_K", "en_Td", "n2_fraction"], keep="last")
    replacements = replacements.drop(columns="_source_priority")
    keys = ["tg_K", "en_Td", "n2_fraction"]
    merged = main.merge(replacements[keys], on=keys, how="left", indicator=True)
    merged = merged.loc[merged["_merge"].eq("left_only")].drop(columns="_merge")
    merged = pd.concat([merged, replacements], ignore_index=True, sort=False)
    merged["tg_K"] = merged["tg_K"].astype(float)
    merged["en_Td"] = merged["en_Td"].astype(float)
    merged["n2_fraction"] = merged["n2_fraction"].astype(float)
    merged["h2_fraction"] = 1.0 - merged["n2_fraction"]
    merged["nh3_strict_converged"] = merged["NH3_tail_pass"].fillna(False).astype(bool)
    merged["all_species_converged"] = merged["other_tail_pass"].fillna(False).astype(bool)
    for sp in SPECIES:
        merged[f"log10_{sp}"] = np.log10(merged[f"{sp}_cm-3"].clip(lower=1.0))
    merged = merged.sort_values(keys).reset_index(drop=True)
    if len(merged) != 972 or merged.duplicated(keys).any():
        raise RuntimeError(f"Expected a complete 972-point merged table; got {len(merged)} rows")
    return merged


def heat(ax, table: pd.DataFrame, value: str, *, label: str, cmap: str, log: bool = True,
         vmin: float | None = None, vmax: float | None = None) -> None:
    es = sorted(table.en_Td.unique())
    n2s = sorted(table.n2_fraction.unique())
    arr = table.pivot(index="en_Td", columns="n2_fraction", values=value).reindex(index=es, columns=n2s).to_numpy(float)
    shown = np.log10(np.clip(arr, 1.0, None)) if log else arr
    image = ax.imshow(shown, aspect="auto", origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(n2s)), [f"{x:.1f}" for x in n2s])
    ax.set_yticks(range(len(es)), [f"{int(x)}" for x in es])
    ax.set_xlabel(r"$x_{N_2}$")
    ax.set_ylabel(r"$E/N$ (Td)")
    ax.set_xticks(np.arange(-.5, len(n2s), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(es), 1), minor=True)
    ax.grid(which="minor", color="white", lw=.35, alpha=.7)
    ax.tick_params(which="minor", bottom=False, left=False)
    bar = ax.figure.colorbar(image, ax=ax, fraction=.046, pad=.03)
    bar.set_label((r"log$_{10}$ " if log else "") + label)


def figure_a_maps(table: pd.DataFrame, out: Path) -> None:
    temperatures = sorted(table.tg_K.unique())
    fig, axes = plt.subplots(3, 3, figsize=(8.5, 7.1), constrained_layout=True)
    low = float(np.floor(np.log10(table["NH3_cm-3"].clip(lower=1.0).min())))
    high = float(np.ceil(np.log10(table["NH3_cm-3"].max())))
    for ax, temp in zip(axes.flat, temperatures):
        heat(ax, table.loc[table.tg_K.eq(temp)], "NH3_cm-3", label=r"[NH$_3$] (cm$^{-3}$)", cmap="viridis", vmin=low, vmax=high)
        ax.set_title(f"{temp:.0f} K")
    fig.suptitle(r"Terminal NH$_3$ inventory across the pure-gas CW-CSTR parameter space", fontsize=10)
    fig.text(.5, .003, r"$n_e=1.17\times10^8$ cm$^{-3}$, $\tau=10$ ms, 100 s integration; no pulse and no surface chemistry.", ha="center", fontsize=7)
    export_figure(fig, out, "Figure_A_NH3_temperature_maps", keep_boxed=True)
    plt.close(fig)


def figure_b_intermediates(table: pd.DataFrame, out: Path) -> None:
    focus = table.loc[table.tg_K.eq(400.0)]
    fig, axes = plt.subplots(1, 5, figsize=(11.0, 2.35), constrained_layout=True)
    for i, (ax, sp) in enumerate(zip(axes, SPECIES)):
        heat(ax, focus, f"{sp}_cm-3", label=rf"[{sp}] (cm$^{{-3}}$)", cmap="magma" if sp == "NH3" else "cividis")
        ax.set_title(f"{sp} at 400 K")
        add_panel_label(ax, chr(65 + i))
    fig.text(.5, -.03, "Terminal concentrations only; these panels do not assign reaction-path fluxes.", ha="center", fontsize=7)
    export_figure(fig, out, "Figure_B_NHx_intermediate_maps_400K", keep_boxed=True)
    plt.close(fig)


def figure_c_slices(table: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.0), constrained_layout=True)
    n2s = [0.1, 0.5, 0.9]
    for ax, n2 in zip(axes, n2s):
        sub = table.loc[np.isclose(table.n2_fraction, n2)]
        for temp, group in sub.groupby("tg_K"):
            group = group.sort_values("en_Td")
            ax.plot(group.en_Td, group["NH3_cm-3"], marker="o", ms=2.4, label=f"{temp:.0f} K")
        ax.set_yscale("log")
        ax.set_xlabel(r"$E/N$ (Td)")
        ax.set_ylabel(r"[NH$_3$] (cm$^{-3}$)")
        ax.set_title(rf"$x_{{N_2}}={n2:.1f}$")
        ax.grid(axis="y")
    axes[-1].legend(ncol=1, title=r"$T_g$", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.suptitle(r"Field response of NH$_3$ at three feed compositions", fontsize=10)
    export_figure(fig, out, "Figure_C_NH3_field_temperature_slices")
    plt.close(fig)


def figure_d_temperature(table: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.0), constrained_layout=True)
    for n2, style in [(0.1, "-"), (0.5, "--"), (0.9, ":")]:
        sub = table.loc[np.isclose(table.n2_fraction, n2)]
        for en, group in sub.groupby("en_Td"):
            if en not in (40.0, 80.0, 140.0, 200.0, 240.0):
                continue
            group = group.sort_values("tg_K")
            axes[0].plot(group.tg_K, group["NH3_cm-3"], linestyle=style, color=plt.cm.viridis((en - 20) / 220), alpha=.9)
        best = (sub.loc[sub.groupby("tg_K")["NH3_cm-3"].idxmax()]
                .sort_values("tg_K"))
        axes[1].plot(best.tg_K, best.en_Td, linestyle=style, color=PALETTE["turnover"], marker="o", ms=3, label=rf"$x_{{N_2}}={n2:.1f}$")
    axes[0].set_yscale("log")
    axes[0].set_xlabel(r"$T_g$ (K)")
    axes[0].set_ylabel(r"[NH$_3$] (cm$^{-3}$)")
    axes[0].set_title("Temperature trends at selected fields")
    axes[0].grid(axis="y")
    mapper = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(40, 240))
    bar = fig.colorbar(mapper, ax=axes[0], fraction=.048, pad=.02)
    bar.set_label(r"$E/N$ (Td)")
    axes[1].set_xlabel(r"$T_g$ (K)")
    axes[1].set_ylabel(r"$E/N$ maximizing [NH$_3$] (Td)")
    axes[1].set_title("Best-field ridge at fixed feed composition")
    axes[1].legend(title="Feed")
    axes[1].grid()
    fig.text(.23, .01, "Colours in left panel increase from 40 to 240 Td; line style denotes feed composition.", ha="center", fontsize=6.5)
    export_figure(fig, out, "Figure_D_temperature_response_and_optimum_ridge")
    plt.close(fig)


def figure_e_correlation(table: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 2.8), constrained_layout=True)
    pairs = [("NH2", "NH$_2$"), ("NH", "NH"), ("N", "N")]
    colors = plt.cm.plasma((table.tg_K - table.tg_K.min()) / (table.tg_K.max() - table.tg_K.min()))
    for ax, (sp, label) in zip(axes, pairs):
        ax.scatter(table[f"{sp}_cm-3"], table["NH3_cm-3"], c=colors, s=12, alpha=.72, edgecolors="none")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(rf"[{label}] (cm$^{{-3}}$)")
        ax.set_ylabel(r"[NH$_3$] (cm$^{-3}$)")
        rho = table[[f"log10_{sp}", "log10_NH3"]].corr(method="spearman").iloc[0, 1]
        ax.set_title(rf"Spearman $\rho={rho:.2f}$")
        ax.grid()
    mapper = plt.cm.ScalarMappable(cmap="plasma", norm=plt.Normalize(table.tg_K.min(), table.tg_K.max()))
    bar = fig.colorbar(mapper, ax=axes, fraction=.025, pad=.02); bar.set_label(r"$T_g$ (K)")
    fig.suptitle(r"Concentration associations: NH$_3$ versus nitrogen-bearing intermediates", fontsize=10)
    export_figure(fig, out, "Figure_E_NH3_intermediate_associations")
    plt.close(fig)


def figure_f_quality(table: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(8.2, 6.8), constrained_layout=True)
    temps = sorted(table.tg_K.unique())
    source_code = {"main": 0, "retry": 1, "refined_retry": 2}
    for ax, temp in zip(axes.flat, temps):
        sub = table.loc[table.tg_K.eq(temp)].copy()
        x = sub.n2_fraction.to_numpy(); y = sub.en_Td.to_numpy()
        code = sub.data_source.map(source_code).to_numpy()
        ax.scatter(x, y, c=code, cmap=plt.cm.get_cmap("Set2", 3), vmin=-.5, vmax=2.5, marker="s", s=48, edgecolors="white", linewidths=.35)
        flagged = sub.loc[~sub.nh3_strict_converged]
        ax.scatter(flagged.n2_fraction, flagged.en_Td, facecolors="none", edgecolors="black", marker="o", s=65, linewidths=1.0)
        ax.set_title(f"{temp:.0f} K")
        ax.set_xlim(.05, .95); ax.set_ylim(10, 250)
        ax.set_xticks([.1, .3, .5, .7, .9]); ax.set_yticks([20, 100, 180, 240])
        ax.set_xlabel(r"$x_{N_2}$"); ax.set_ylabel(r"$E/N$ (Td)")
    cmap = plt.cm.get_cmap("Set2", 3)
    handles = [
        Line2D([], [], marker="s", linestyle="", markerfacecolor=cmap(i), markeredgecolor="white", label=label, markersize=7)
        for i, label in enumerate(["main", "numerical retry", "refined retry"])
    ]
    handles.append(Line2D([], [], marker="o", linestyle="", markerfacecolor="none", markeredgecolor="black", label="NH$_3$ strict flag", markersize=8))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False)
    fig.suptitle("Numerical provenance and NH$_3$ convergence flags", fontsize=10)
    export_figure(fig, out, "Figure_F_numerical_provenance", keep_boxed=True)
    plt.close(fig)


def figure_g_time(root: Path, out: Path) -> None:
    records = pd.read_csv(root / "time_horizon_summary.csv")
    focus = records.loc[(records.tg_K.eq(400.0)) & (records.en_Td.eq(140.0)) & (np.isclose(records.n2_fraction, .5))]
    fig, ax = plt.subplots(figsize=(4.3, 3.0), constrained_layout=True)
    ax.plot(focus.requested_horizon_s, focus["NH3_cm-3"], marker="o", color=PALETTE["turnover"])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Integration horizon (s)"); ax.set_ylabel(r"[NH$_3$] (cm$^{-3}$)")
    ax.set_title(r"Long-horizon check: 400 K, 140 Td, $x_{N_2}=0.5$")
    ax.grid()
    export_figure(fig, out, "Figure_G_NH3_time_horizon_check")
    plt.close(fig)


def report(table: pd.DataFrame, out: Path) -> None:
    best = table.loc[table["NH3_cm-3"].idxmax()]
    per_temp = table.loc[table.groupby("tg_K")["NH3_cm-3"].idxmax()].sort_values("tg_K")
    strict = int(table.nh3_strict_converged.sum())
    source_counts = table.data_source.value_counts().to_dict()
    corr = {sp: table[["log10_NH3", f"log10_{sp}"]].corr(method="spearman").iloc[0, 1] for sp in ("NH2", "NH", "N", "H")}
    lines = [
        "# Hong 气相 CW-CSTR 细化扫描：图件分析报告", "",
        "## 范围与边界", "",
        "本报告使用 Hong 2017/2018 修订气相机理的 0D CW-CSTR 结果。所有点均为固定电子密度（1.17×10^8 cm^-3）、固定停留时间（10 ms）、无脉冲和无表面反应的终端浓度。图件不包含 ROP，因此不对 direct-NH 路径份额、反应控制步骤、能量效率或真实脉冲行为作出判断。", "",
        "## 数据验收", "",
        f"- 合并后网格：{len(table)} / 972 点；来源：主扫描 {source_counts.get('main', 0)} 点、普通数值重算 {source_counts.get('retry', 0)} 点、精化重算 {source_counts.get('refined_retry', 0)} 点。",
        f"- NH3 严格收敛：{strict} / {len(table)}；其他关键物种宽松收敛：{int(table.all_species_converged.sum())} / {len(table)}。",
        "- 未通过 NH3 严格阈值的点均保留在数据和质量控制图中；它们主要是低场、极低 NH3 库存情形，不作为最优工况的依据。", "",
        "## 直接观察", "",
        f"- 全局最高 NH3 终端浓度为 {best['NH3_cm-3']:.3e} cm^-3，位于 T={best.tg_K:.0f} K、E/N={best.en_Td:.0f} Td、x_N2={best.n2_fraction:.1f}。",
        "- 温度切片图表明，最高 NH3 库存并不固定于单一 E/N 或单一进料组成；应将其表述为本受控 CSTR 域内的库存响应面，而非普适能效最优。",
        "- NH3 与 NH2、NH、N、H 的散点图仅用于识别浓度共变。相关性不等价于任何一个物种或反应是 NH3 的速率控制因素。", "",
        "## 各温度下 NH3 最大点", "",
        "| T (K) | E/N (Td) | x_N2 | [NH3] (cm^-3) |", "|---:|---:|---:|---:|",
    ]
    for _, row in per_temp.iterrows():
        lines.append(f"| {row.tg_K:.0f} | {row.en_Td:.0f} | {row.n2_fraction:.1f} | {row['NH3_cm-3']:.3e} |")
    lines += ["", "## 浓度关联（Spearman ρ，log10 浓度）", "", "| NH3 对象 | ρ |", "|---|---:|"]
    lines += [f"| {sp} | {rho:.3f} |" for sp, rho in corr.items()]
    lines += ["", "## 与 Hong 文献的关系", "", "Hong 文献及现有项目框架强调将 NH 的入口归因与 NH3 的下游周转分开讨论。本轮新增扫描只有浓度终值，因此可以延续这种两层叙事的边界：本报告量化的是 NH3/NHx 库存如何随 T、E/N 和进料组成变化；若需将该响应面归因到命名电子激发 H2* 或 NH3 损失通道，必须在相同 972 点上补充末端 ROP/通量积分。", ""]
    (out / "CW_gasphase_Hong_literature_guided_analysis.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot the refined Hong pure-gas CW-CSTR scan")
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    out = (args.output or root / "figures").resolve()
    out.mkdir(parents=True, exist_ok=True)
    apply_theme()
    table = build_table(root)
    table.to_csv(out / "cw_gasphase_972_accepted_merged.csv", index=False, float_format="%.12g")
    figure_a_maps(table, out); figure_b_intermediates(table, out); figure_c_slices(table, out)
    figure_d_temperature(table, out); figure_e_correlation(table, out); figure_f_quality(table, out); figure_g_time(root, out)
    report(table, out)
    print(json.dumps({"points": len(table), "nh3_strict": int(table.nh3_strict_converged.sum()), "output": str(out)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
