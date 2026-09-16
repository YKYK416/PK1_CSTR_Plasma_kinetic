#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create main-text graph-theoretical reaction-path figures for Hong kinetics.

Figure 6 draws reaction hypergraphs from terminal steady-state ROP fluxes.
Figure 7 reports graph-derived pathway hierarchy together with the validated
108-point direct-NH source field.  No rejected afterglow cases are read.
"""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import networkx as nx
import numpy as np
import pandas as pd

from _paths import LEGACY_HONG_ROOT
from hong_reaction_hypergraph import (
    PROJECT, SUMMARY, ReactionFlux, bipartite_hypergraph, display_reactions,
    integrate_case, molecular_formula, nitrogen_lineage_graph, validate_records,
)


FIGURES = PROJECT / "analysis" / "figures_20260829"
OUT = FIGURES
CASES = [(20, 0.1, "Low field"), (140, 0.1, "NH$_3$ maximum"),
         (240, 0.9, "High field, N$_2$-rich")]
COLORS = {
    "named electronic H2*": "#0072B2", "Rydberg-H2*": "#CC79A7",
    "N(2D)+H2": "#56B4E9", "N(2P)+H2": "#E69F00",
    "vibrational H2": "#009E73", "H+N+M": "#D55E00", "other": "#777777",
    "downstream": "#009E73", "loss": "#D55E00",
}


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7.5, "axes.labelsize": 8.5,
        "axes.titlesize": 8.7, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 6.3, "axes.linewidth": 0.7, "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def panel(ax: plt.Axes, label: str) -> None:
    ax.text(-0.07, 1.07, label, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top")


def save(fig: plt.Figure, output: Path, stem: str) -> None:
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(output / f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def reaction_short(record: ReactionFlux) -> str:
    left = "+".join(record.reactants).replace("RYDBERG_SUM", "Ryd.")
    right = "+".join(record.products)
    return f"{left} → {right}"


def is_nhx_production(record: ReactionFlux, product: str) -> bool:
    products = [molecular_formula(x) for x in record.products]
    reactants = [molecular_formula(x) for x in record.reactants]
    return product in products and product not in reactants


def is_nh_loss(record: ReactionFlux) -> bool:
    reactants = [molecular_formula(x) for x in record.reactants]
    products = [molecular_formula(x) for x in record.products]
    return "NH" in reactants and "NH" not in products


def pathway_subset(records: list[ReactionFlux]) -> list[ReactionFlux]:
    """Keep readable source-family hyperedges plus NHx formation/loss steps.

    The four named electronic H2* reactions (and the two third-body variants)
    are merged only for drawing.  Their individual ROP values remain present
    in the reaction table and Figure 7D.
    """
    direct = [row for row in records if row.route != "other"]
    selected: list[ReactionFlux] = []
    family_reactants = {
        "named electronic H2*": (("N", "H2(named H2*)"), ("H", "NH")),
        "Rydberg-H2*": (("N", "H2(Ryd.)"), ("H", "NH")),
        "N(2D)+H2": (("N(2D)", "H2"), ("H", "NH")),
        "N(2P)+H2": (("N(2P)", "H2"), ("H", "NH")),
        "vibrational H2": (("N", "H2(v=1–3)"), ("H", "NH")),
        "H+N+M": (("H", "N", "M"), ("NH", "M")),
    }
    for route, (reactants, products) in family_reactants.items():
        rows = [row for row in direct if row.route == route]
        phi = sum(row.phi_cm3 for row in rows)
        if phi <= 0:
            continue
        selected.append(ReactionFlux(
            reaction=f"{route} family", reactants=tuple(reactants), products=tuple(products),
            phi_cm3=phi, route=route, n_balanced=True, h_balanced=True))
    for product, count in (("NH2", 1), ("NH3", 2)):
        candidates = [row for row in records if is_nhx_production(row, product)]
        selected.extend(candidates[:count])
    selected.extend([row for row in records if is_nh_loss(row)][:1])
    unique = {row.reaction: row for row in selected}
    return sorted(unique.values(), key=lambda row: row.phi_cm3, reverse=True)


def route_share(rows: list[ReactionFlux]) -> dict[str, float]:
    direct = [row for row in rows if row.route != "other"]
    total = sum(row.phi_cm3 for row in direct)
    result: dict[str, float] = defaultdict(float)
    for row in direct:
        result[row.route] += 100.0 * row.phi_cm3 / total if total else float("nan")
    return dict(result)


def reactions_for_node(records: list[ReactionFlux], formula: str, product: bool = True) -> list[ReactionFlux]:
    side = (lambda row: row.products) if product else (lambda row: row.reactants)
    return [row for row in records if any(molecular_formula(s) == formula for s in side(row))]


def representative_graph_metrics(records: list[ReactionFlux]) -> dict[str, float]:
    """Compute a reproducible N-lineage centrality proxy from the full hypergraph.

    The incidence graph is built and checked to retain reaction arity.  The
    N-lineage projection then allocates only conserved N flux among N-bearing
    participants; weighted betweenness is therefore a topology metric, not a
    substitute for the flux shares used to identify NH sources.
    """
    hyper = bipartite_hypergraph(records)
    lineage = nitrogen_lineage_graph(records)
    if lineage.number_of_nodes() == 0:
        return {"hyper_nodes": hyper.number_of_nodes(), "hyper_edges": hyper.number_of_edges()}
    max_flux = max(float(data["flux_cm3"]) for _, _, data in lineage.edges(data=True))
    for _, _, data in lineage.edges(data=True):
        data["cost"] = math.log(max_flux / max(float(data["flux_cm3"]), 1e-300)) + 1e-12
    centrality = nx.betweenness_centrality(lineage, weight="cost", normalized=True)
    metrics = {"hyper_nodes": hyper.number_of_nodes(), "hyper_edges": hyper.number_of_edges(),
               "lineage_nodes": lineage.number_of_nodes(), "lineage_edges": lineage.number_of_edges()}
    for species in ("N2", "N", "N(2D)", "N(2P)", "NH", "NH2", "NH3"):
        metrics[f"bc_{species}"] = centrality.get(species, 0.0)
    return metrics


def node_layer(species: str) -> tuple[float, str]:
    formula = molecular_formula(species)
    if species.startswith("N("):
        return 0.31, species
    if species == "M":
        return 0.17, "M"
    if formula == "N2":
        return 0.08, "N$_2$ reservoir"
    if formula in {"N", "NH"}:
        return (0.40 if formula == "N" else 0.66), formula
    if formula == "NH2":
        return 0.81, formula
    if formula == "NH3":
        return 0.94, formula
    if formula == "H2":
        return 0.20, species.replace("RYDBERG_SUM", "Ryd.")
    if formula == "H":
        return 0.53, "H"
    return 0.24, species.replace("RYDBERG_SUM", "Ryd.")


def draw_hypergraph(ax: plt.Axes, records: list[ReactionFlux], title: str) -> None:
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off"); ax.set_title(title, pad=2)
    selected = pathway_subset(records)
    direct = [row for row in selected if row.route != "other"]
    direct_phi = sum(row.phi_cm3 for row in records if row.route != "other")

    def box(x: float, y: float, label: str, *, face: str = "#F5F5F5", width: float = 0.12) -> None:
        ax.add_patch(FancyBboxPatch((x - width / 2, y - 0.028), width, 0.056,
                                    boxstyle="round,pad=0.008", facecolor=face,
                                    edgecolor="#555555", lw=0.5, zorder=6))
        ax.text(x, y, label, ha="center", va="center", fontsize=5.3, zorder=7)

    def arrow(src: tuple[float, float], dst: tuple[float, float], color: str, width: float,
              *, dashed: bool = False) -> None:
        ax.add_patch(FancyArrowPatch(src, dst, arrowstyle="-|>", mutation_scale=7.2, lw=width,
                                     color=color, alpha=0.82, linestyle="--" if dashed else "-",
                                     connectionstyle="arc3,rad=0.0", zorder=3))

    # Shared species nodes.  Reaction squares are explicit hypernodes.
    pos = {"N": (0.08, 0.68), "H": (0.08, 0.48), "M": (0.08, 0.28),
           "NH": (0.69, 0.68), "Hout": (0.69, 0.48), "Mout": (0.69, 0.28),
           "NH2": (0.84, 0.42), "NH3": (0.94, 0.24)}
    box(*pos["N"], "N", face="#EAF3F8", width=0.09)
    box(*pos["H"], "H", width=0.07); box(*pos["M"], "M", width=0.07)
    box(*pos["NH"], "NH", face="#DDEEDB", width=0.09)
    box(*pos["Hout"], "H", width=0.07); box(*pos["Mout"], "M", width=0.07)
    box(*pos["NH2"], r"NH$_2$", face="#DDEEDB", width=0.10)
    box(*pos["NH3"], r"NH$_3$", face="#DDEEDB", width=0.10)

    ys = np.linspace(0.90, 0.20, len(direct))
    labels: list[str] = []
    short_route = {"named electronic H2*": "named H$_2$*", "Rydberg-H2*": "Rydberg-H$_2$*",
                   "N(2D)+H2": "N($^2$D)+H$_2$", "N(2P)+H2": "N($^2$P)+H$_2$",
                   "vibrational H2": "vibrational H$_2$", "H+N+M": "H+N+M"}
    for number, (row, y) in enumerate(zip(direct, ys), start=1):
        route = row.route
        color = COLORS[route]
        rel = row.phi_cm3 / max(direct_phi, 1e-300)
        width = 0.45 + 2.45 * max(0.0, min(1.0, (math.log10(max(rel, 1e-12)) + 12) / 12))
        # State/co-reactant node is separate from N, so every required reactant is visible.
        label = next((item for item in row.reactants if item not in {"N", "H", "M"}), "")
        state_pos = (0.28, float(y))
        if route in {"N(2D)+H2", "N(2P)+H2"}:
            atom_pos, h2_pos = (0.21, float(y)), (0.34, float(y))
            box(*atom_pos, route.split("+")[0], width=0.11)
            box(*h2_pos, r"H$_2$", width=0.07)
            arrow(atom_pos, (0.46, y), color, width)
            arrow(h2_pos, (0.46, y), color, width * 0.55)
        elif route == "H+N+M":
            box(*state_pos, "H + N + M", width=0.14)
            for source in (pos["N"], pos["H"], pos["M"]):
                arrow(source, (0.46, y), color, width * 0.5)
        else:
            box(*state_pos, label, width=0.15)
            arrow(pos["N"], (0.46, y), color, width * 0.65)
            arrow(state_pos, (0.46, y), color, width)
        ax.add_patch(FancyBboxPatch((0.46, y - 0.021), 0.042, 0.042, boxstyle="round,pad=0.004",
                                    facecolor=color, edgecolor="white", lw=0.5, zorder=6))
        ax.text(0.481, y, str(number), ha="center", va="center", fontsize=5.3, color="white", zorder=7)
        arrow((0.505, y), pos["NH"], color, width)
        if route == "H+N+M":
            arrow((0.505, y), pos["Mout"], color, width * 0.5)
        else:
            arrow((0.505, y), pos["Hout"], color, width * 0.5)
        labels.append(f"{number}. {short_route[route]}")

    # The highest-flux downstream NHx conversion reactions are shown as a separate layer.
    downstream = [row for row in selected if row.route == "other" and
                  any(molecular_formula(item) in {"NH2", "NH3"} for item in row.products)]
    for number, row in enumerate(downstream[:3], start=len(labels) + 1):
        product = "NH3" if is_nhx_production(row, "NH3") else "NH2"
        y = 0.10 if product == "NH3" else 0.34
        rpos = (0.78, y)
        ax.add_patch(FancyBboxPatch((rpos[0] - 0.021, rpos[1] - 0.021), 0.042, 0.042,
                                    boxstyle="round,pad=0.004", facecolor=COLORS["downstream"],
                                    edgecolor="white", lw=0.5, zorder=6))
        ax.text(*rpos, str(number), ha="center", va="center", fontsize=5.3, color="white", zorder=7)
        source = pos["NH"] if any(molecular_formula(item) == "NH" for item in row.reactants) else pos["N"]
        target = pos[product]
        arrow(source, rpos, COLORS["downstream"], 1.1)
        arrow(rpos, target, COLORS["downstream"], 1.1)
        labels.append(f"{number}. NHx formation")

    ax.text(0.50, 0.985, r"species nodes → reaction hypernodes → products; edge width: terminal ROP",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.1)
    ax.text(0.01, 0.01, "\n".join(labels), transform=ax.transAxes, fontsize=4.8, va="bottom",
            bbox={"facecolor": "white", "edgecolor": "#BBBBBB", "alpha": 0.93, "pad": 1.4})


def make_figure6(analyses: list[tuple[list[ReactionFlux], dict, str]], output: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12.55, 4.45), constrained_layout=True)
    for ax, (records, meta, condition), letter in zip(axes, analyses, "ABC"):
        title = f"{condition}: {meta['en_td']} Td, $x_{{N_2}}$={meta['n2_fraction']:.1f}"
        draw_hypergraph(ax, records, title)
        panel(ax, letter)
    fig.text(0.5, -0.045,
             "Each coloured square is an explicit reaction hypernode; edge width is proportional to terminal-three-cycle ROP. "
             "The drawings retain co-reactants and separate direct NH-entry families from downstream NH$_x$ formation.",
             ha="center", fontsize=7)
    save(fig, output, "Figure_6_flux_weighted_reaction_hypergraphs")


def make_figure7(analyses: list[tuple[list[ReactionFlux], dict, str]], output: Path) -> pd.DataFrame:
    dense = pd.read_csv(SUMMARY)
    channels = ["rydberg_h2star_pct", "n2d_h2_pct", "n2p_h2_pct", "vib_h2_pct", "association_pct"]
    secondary = dense[channels].max(axis=1)
    dense["named_to_secondary"] = dense["named_h2star_pct"] / secondary.clip(lower=1e-15)
    n2s, ens = sorted(dense.n2_fraction.unique()), sorted(dense.en_td.unique())
    metrics_rows = []
    for records, meta, condition in analyses:
        validation = validate_records(records)
        metrics = representative_graph_metrics(records)
        shares = route_share(records)
        metrics_rows.append({**meta, "condition": condition, **validation, **metrics, **shares})
    metrics_df = pd.DataFrame(metrics_rows)
    metrics_df.to_csv(output / "Figure_7_hypergraph_metrics.csv", index=False)

    fig = plt.figure(figsize=(13.0, 6.35), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, width_ratios=(1.0, 1.0, 0.86), height_ratios=(1, 1.1))

    ax = fig.add_subplot(grid[0, 0])
    field = dense.pivot(index="en_td", columns="n2_fraction", values="named_to_secondary").loc[ens, n2s]
    im = ax.imshow(np.log10(field), origin="lower", aspect="auto", cmap="cividis")
    ax.set_xticks(range(len(n2s)), [f"{v:.1f}" for v in n2s])
    ax.set_yticks(range(0, len(ens), 2), [str(ens[i]) for i in range(0, len(ens), 2)])
    ax.set_xlabel(r"N$_2$ inlet fraction"); ax.set_ylabel("E/N (Td)")
    ax.set_title("Graph field: direct-NH dominance ratio")
    cb = fig.colorbar(im, ax=ax, fraction=0.047, pad=0.04)
    cb.set_label(r"log$_{10}$($\Phi_{\rm named}/\Phi_{\rm largest\ secondary}$)")
    panel(ax, "A")

    ax = fig.add_subplot(grid[0, 1])
    draw_channels = ["named electronic H2*", "Rydberg-H2*", "N(2P)+H2", "H+N+M"]
    x = np.arange(len(metrics_df))
    bottom = np.zeros(len(metrics_df))
    for channel in draw_channels:
        vals = metrics_df.get(channel, pd.Series(0.0, index=metrics_df.index)).fillna(0.0).to_numpy(float)
        ax.bar(x, vals, bottom=bottom, color=COLORS[channel], width=0.65, label=channel)
        bottom += vals
    ax.set_xticks(x, [f"{int(row.en_td)} Td\n{row.n2_fraction:.1f}" for row in metrics_df.itertuples()])
    ax.set_ylim(0, 101); ax.set_ylabel("Direct NH-source flux (%)")
    ax.set_title("Hypergraph source-family decomposition")
    ax.legend(frameon=False, ncol=2, loc="upper right")
    panel(ax, "B")

    ax = fig.add_subplot(grid[1, 0])
    species = ["N2", "N", "N(2D)", "N(2P)", "NH", "NH2", "NH3"]
    mat = np.array([[float(row.get(f"bc_{species_name}", 0.0)) for species_name in species]
                    for _, row in metrics_df.iterrows()])
    im = ax.imshow(mat, aspect="auto", cmap="magma")
    ax.set_xticks(range(len(species)), [r"N$_2$", "N", r"N($^2$D)", r"N($^2$P)", "NH", r"NH$_2$", r"NH$_3$"], rotation=30, ha="right")
    ax.set_yticks(range(len(metrics_df)), [f"{int(row.en_td)} Td, {row.n2_fraction:.1f}" for row in metrics_df.itertuples()])
    ax.set_title("N-lineage betweenness in the full hypergraph")
    cb = fig.colorbar(im, ax=ax, fraction=0.047, pad=0.04); cb.set_label("weighted betweenness")
    panel(ax, "C")

    ax = fig.add_subplot(grid[1, 1])
    # Each point is one raw direct NH reaction, so the pathway hierarchy stays reaction-resolved.
    for idx, (records, meta, condition) in enumerate(analyses):
        direct = [row for row in records if row.route != "other"]
        total = sum(row.phi_cm3 for row in direct)
        ordered = sorted(direct, key=lambda row: row.phi_cm3, reverse=True)
        y = [100 * row.phi_cm3 / total for row in ordered]
        colors = [COLORS[row.route] for row in ordered]
        ax.scatter(np.full(len(y), idx), y, s=18, c=colors, edgecolor="black", linewidth=0.25, zorder=3)
        # Ties point labels to the four named electronic components without cluttering minor routes.
        for rank, (row, value) in enumerate(zip(ordered, y), start=1):
            if rank <= 4:
                ax.annotate(str(rank), (idx, value), xytext=(3, 1), textcoords="offset points", fontsize=5.5)
    ax.set_yscale("log"); ax.set_ylim(1e-13, 110)
    ax.set_xticks(range(len(analyses)), [f"{m['en_td']} Td\n$m_{{N_2}}$={m['n2_fraction']:.1f}" for _, m, _ in analyses])
    ax.set_ylabel("Individual direct NH-reaction flux (%)")
    ax.set_title("Reaction-level pathway hierarchy (rank labels 1–4)")
    handles = [mpl.lines.Line2D([], [], marker="o", linestyle="", color=COLORS[name], label=name)
               for name in draw_channels]
    ax.legend(handles=handles, frameon=False, loc="lower right", ncol=1)
    panel(ax, "D")

    # The Hong labels are incomplete term symbols.  This is a classification
    # ledger, deliberately not a selection-rule or rate-prediction claim.
    ax = fig.add_subplot(grid[0, 2]); ax.axis("off")
    ax.set_title("State-symmetry evidence boundary", loc="left", pad=8)
    ax.text(0.00, 0.93, r"Isolated H$_2$ belongs to $D_{\infty h}$; collision kinetics", fontsize=6.7)
    ax.text(0.00, 0.87, "also require dynamical coupling and state-resolved cross sections.", fontsize=6.7)
    rows = [
        [r"H$_2$(v=0–3)", r"X $^1\Sigma_g^+$", "complete term"],
        ["B3SIG / A3SIG", r"$^3\Sigma$", "g/u, +/- omitted"],
        ["B1SIG", r"$^1\Sigma$", "g/u, +/- omitted"],
        ["C3PI", r"$^3\Pi$", "g/u, +/- omitted"],
        [r"Rydberg-H$_2$*", "unspecified", "term unresolved"],
        [r"N($^2$P), N($^2$D)", "atomic", r"not $D_{\infty h}$"],
    ]
    table = ax.table(cellText=rows, colLabels=["family", r"$S,\Lambda$", "evidence status"],
                     cellLoc="left", colLoc="left", colWidths=[0.37, 0.23, 0.40], bbox=[0.0, 0.23, 1.0, 0.56])
    table.auto_set_font_size(False); table.set_fontsize(5.8)
    for (row, _column), cell in table.get_celld().items():
        cell.set_linewidth(0.25)
        if row == 0:
            cell.set_facecolor("#DDEBF7"); cell.set_text_props(weight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#F7F7F7")
    ax.text(0.00, 0.13, "Use: classify spin/orbital labels. Do not infer", fontsize=6.5, weight="bold")
    ax.text(0.00, 0.07, "a selection rule or rate ordering from abbreviation alone.", fontsize=6.5)
    ax.text(0.00, 0.01, "The ROP field, not term-label completeness, ranks routes.", fontsize=6.5)
    panel(ax, "E")

    ax = fig.add_subplot(grid[1, 2]); ax.axis("off")
    ax.set_title("Graph-theoretical interpretation", loc="left", pad=8)
    n_balanced = bool((metrics_df["n_unbalanced"] == 0).all())
    h_balanced = bool((metrics_df["h_unbalanced"] == 0).all())
    source_range = (dense.named_h2star_pct.min(), dense.named_h2star_pct.max())
    ax.text(0.00, 0.86, "1. Hypergraph: every reaction is a bipartite", fontsize=7.0, weight="bold")
    ax.text(0.04, 0.78, "species–reaction hyperedge; stoichiometric arity is retained.", fontsize=6.8)
    ax.text(0.00, 0.63, "2. N-lineage: conserved N flux produces a weighted", fontsize=7.0, weight="bold")
    ax.text(0.04, 0.55, "projection used only for topological centrality (panel C).", fontsize=6.8)
    ax.text(0.00, 0.40, "3. Flux hierarchy: terminal ROP determines direct", fontsize=7.0, weight="bold")
    ax.text(0.04, 0.32, f"NH entry; named H$_2$* spans {source_range[0]:.1f}–{source_range[1]:.1f}% (panel A).", fontsize=6.8)
    ax.text(0.00, 0.15, f"Atom-balance check: N = {n_balanced}; H = {h_balanced}.", fontsize=6.7,
            bbox={"facecolor": "#F5F5F5", "edgecolor": "#AAAAAA", "pad": 1.5})
    panel(ax, "F")
    fig.text(0.5, -0.055,
             "Graph metrics describe the reaction-network topology; integrated reaction-of-progress fluxes provide the pathway ranking. "
             "No unrecorded molecular parity or collision-selection rule is assumed.", ha="center", fontsize=7)
    save(fig, output, "Figure_7_graph_theoretical_path_hierarchy")
    return metrics_df


def write_reaction_table(analyses: list[tuple[list[ReactionFlux], dict, str]], output: Path) -> None:
    fields = ["condition", "en_td", "n2_fraction", "reaction", "route", "phi_cm-3", "direct_nh_share_pct",
              "reactants", "products", "N_atom_balanced", "H_atom_balanced"]
    with (output / "Figure_6_reaction_flux_table.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader()
        for records, meta, condition in analyses:
            direct = [row for row in records if row.route != "other"]
            total = sum(row.phi_cm3 for row in direct)
            for row in records:
                writer.writerow({"condition": condition, "en_td": meta["en_td"], "n2_fraction": meta["n2_fraction"],
                                 "reaction": row.reaction, "route": row.route, "phi_cm-3": row.phi_cm3,
                                 "direct_nh_share_pct": 100 * row.phi_cm3 / total if row.route != "other" and total else "",
                                 "reactants": " + ".join(row.reactants), "products": " + ".join(row.products),
                                 "N_atom_balanced": row.n_balanced, "H_atom_balanced": row.h_balanced})


def main() -> int:
    parser = argparse.ArgumentParser(description="Hong ROP hypergraph figure generator")
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    configure()
    analyses = []
    for en_td, n2_fraction, name in CASES:
        records, meta = integrate_case(en_td, n2_fraction)
        analyses.append((records, meta, name))
        print(f"{meta['tag']}: {validate_records(records)}")
    write_reaction_table(analyses, args.output)
    make_figure6(analyses, args.output)
    make_figure7(analyses, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
