#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Topology-only audit of the explicit Fe(110) surface subnetwork.

The Shao--Mesbah mechanism is parsed as supplied.  This script intentionally
does not load rate coefficients, simulate a reactor, or infer a surface flux.
Its claims are restricted to the declared reaction topology and to the first
surface N--H bond-forming reactions encoded in the frozen input.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from collections import Counter
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import networkx as nx

from _paths import LEGACY_HONG_ROOT

PROJECT = LEGACY_HONG_ROOT
SOURCE = (PROJECT / "Reproduction" / "external_models" / "DFT-microkinetic" /
          "Model_SA_Const_Entropy_base" /
          "kinet_varyT_metal_auto_sens_DFT_in_entropy_verying_basis.txt")
DEFAULT_OUTPUT = PROJECT / "analysis" / "p14_fe110_surface_graph_audit_20260901_r3"

from hong_reaction_hypergraph import parse_reaction  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def formula_token(token: str) -> str:
    """Remove a leading stoichiometric coefficient only."""
    return token.strip().lstrip("0123456789")


def is_surface(token: str) -> bool:
    return "Surf" in formula_token(token)


def load_surface_reactions(path: Path) -> list[dict]:
    rows = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#") or "=>" not in line:
            continue
        reaction = line.split("!", 1)[0].strip()
        parsed = parse_reaction(reaction)
        if parsed is None:
            continue
        reactants, products = parsed
        if not any(is_surface(item) for item in reactants + products):
            continue
        rows.append({"line_number": line_number, "reaction": reaction,
                     "reactants": tuple(formula_token(item) for item in reactants),
                     "products": tuple(formula_token(item) for item in products)})
    if len(rows) != 46:
        raise RuntimeError(f"Expected 46 explicit surface reactions, found {len(rows)}. Source format changed.")
    return rows


def first_nh_family(row: dict) -> str | None:
    left, right = set(row["reactants"]), set(row["products"])
    if "NHSurf" not in right:
        return None
    if {"NSurf", "HSurf"}.issubset(left):
        return "Langmuir–Hinshelwood: NSurf + HSurf"
    if "HSurf" in left and any(item in left for item in ("N", "N(2D)", "N(2P)")):
        return "Eley–Rideal: gas N(x) + HSurf"
    if "NSurf" in left and "H" in left:
        return "Eley–Rideal: H + NSurf"
    return None


def reaction_family(row: dict) -> str:
    left, right = set(row["reactants"]), set(row["products"])
    if first_nh_family(row):
        return "first surface N–H formation"
    if "NSurf" in right and any(item.startswith("N2") for item in left):
        return "N2 dissociative adsorption"
    if "HSurf" in right and any(item.startswith("H2") for item in left):
        return "H2 dissociative adsorption"
    if "NSurf" in right and any(item in left for item in ("N", "N(2D)", "N(2P)")):
        return "N adsorption"
    if "HSurf" in right and "H" in left:
        return "H adsorption"
    if "NHSurf" in right and "NH" in left:
        return "NH adsorption / hydrogenation"
    if "NH2Surf" in right and "NH2" in left:
        return "NH2 adsorption / hydrogenation"
    if "NH2Surf" in right and "NHSurf" in left:
        return "NHSurf → NH2Surf"
    if "NH3" in right:
        return "NHxSurf → NH3 release"
    if "N2" in right and "NSurf" in left:
        return "N removal / N2 formation"
    if "H2" in right and "HSurf" in left:
        return "H recombination / H2 formation"
    return "other explicit surface process"


def build_bipartite(rows: list[dict]) -> nx.DiGraph:
    graph = nx.DiGraph()
    for index, row in enumerate(rows, start=1):
        reaction_node = f"r{index:02d}"
        graph.add_node(reaction_node, kind="reaction", reaction=row["reaction"], line=row["line_number"])
        for species in row["reactants"]:
            graph.add_node(species, kind="species")
            graph.add_edge(species, reaction_node)
        for species in row["products"]:
            graph.add_node(species, kind="species")
            graph.add_edge(reaction_node, species)
    return graph


def write_ledgers(output: Path, rows: list[dict], graph: nx.DiGraph) -> tuple[list[dict], list[dict]]:
    ledger = []
    first = []
    for row in rows:
        family = reaction_family(row)
        record = {"line_number": row["line_number"], "reaction": row["reaction"], "family": family,
                  "reactants": " + ".join(row["reactants"]), "products": " + ".join(row["products"])}
        ledger.append(record)
        entry = first_nh_family(row)
        if entry:
            first.append({**record, "first_nh_class": entry})
    with (output / "P14_Fe110_surface_reaction_ledger.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=ledger[0].keys())
        writer.writeheader(); writer.writerows(ledger)
    with (output / "P14_first_surface_NH_entries.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=first[0].keys())
        writer.writeheader(); writer.writerows(first)
    species = [node for node, data in graph.nodes(data=True) if data["kind"] == "species"]
    nontrivial_scc = [nodes for nodes in nx.strongly_connected_components(graph) if len(nodes) > 1]
    stats = [{"metric": "explicit_surface_reactions", "value": len(rows)},
             {"metric": "species_nodes_in_bipartite_projection", "value": len(species)},
             {"metric": "reaction_nodes_in_bipartite_projection", "value": len(rows)},
             {"metric": "first_surface_NH_entries", "value": len(first)},
             {"metric": "nontrivial_strongly_connected_components", "value": len(nontrivial_scc)}]
    with (output / "P14_surface_graph_metrics.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=("metric", "value")); writer.writeheader(); writer.writerows(stats)
    return ledger, first


def node(ax, xy, text, color, width=0.16, height=0.085, fontsize=7.0):
    x, y = xy
    ax.add_patch(FancyBboxPatch((x - width / 2, y - height / 2), width, height, boxstyle="round,pad=0.008",
                                facecolor=color, edgecolor="#555555", linewidth=0.55, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, zorder=4)


def arrow(ax, start, end, text="", color="#6D6D6D", bend=0.0, fontsize=5.8):
    patch = FancyArrowPatch(start, end, connectionstyle=f"arc3,rad={bend}", arrowstyle="-|>",
                            mutation_scale=8.5, linewidth=0.85, color=color, zorder=2)
    ax.add_patch(patch)
    if text:
        ax.text((start[0] + end[0]) / 2, (start[1] + end[1]) / 2 + (0.025 if bend >= 0 else -0.025), text,
                ha="center", va="center", fontsize=fontsize, color="#303030")


def draw_figure(output: Path, ledger: list[dict], first: list[dict], graph: nx.DiGraph) -> None:
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.3, "axes.titlesize": 9.2,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(12.2, 8.4))
    grid = fig.add_gridspec(2, 2, left=0.07, right=0.98, bottom=0.085, top=0.93, hspace=0.47, wspace=0.28,
                            width_ratios=(1.10, 0.90))

    ax = fig.add_subplot(grid[0, 0]); ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title("Explicit Fe(110) surface branch: structural N–H-entry routes", loc="left", pad=8)
    node(ax, (0.09, 0.74), "N$_2$\n+N$_2$(x)", "#DDEBF7")
    node(ax, (0.09, 0.27), "H$_2$\n+H$_2$(x)", "#FCE4D6")
    node(ax, (0.36, 0.74), "NSurf", "#E2F0D9")
    node(ax, (0.36, 0.27), "HSurf", "#FFF2CC")
    node(ax, (0.60, 0.52), "NHSurf\n(first N–H)", "#F4CCCC", width=0.20)
    node(ax, (0.79, 0.52), "NH$_2$Surf", "#F4CCCC")
    node(ax, (0.94, 0.52), "NH$_3$\n(gas)", "#EADCF8", width=0.11)
    arrow(ax, (0.17, 0.74), (0.28, 0.74), "13 N$_2$(x) dissociative\nadsorption reactions")
    arrow(ax, (0.17, 0.27), (0.28, 0.27), "8 H$_2$(x) dissociative\nadsorption reactions")
    arrow(ax, (0.43, 0.27), (0.49, 0.50), "LH", color="#C55A11", bend=0.10)
    arrow(ax, (0.43, 0.74), (0.49, 0.54), "ER: H + NSurf", color="#C55A11", bend=-0.10)
    arrow(ax, (0.23, 0.50), (0.49, 0.54), "ER: N(x) + HSurf\n(3 N states)", color="#C55A11", bend=-0.20)
    ax.text(0.10, 0.51, "gas N, N(2D), N(2P)", fontsize=6.1, ha="left", va="center")
    arrow(ax, (0.70, 0.52), (0.705, 0.52))
    arrow(ax, (0.87, 0.52), (0.885, 0.52))
    ax.text(0.78, 0.62, "H / HSurf", fontsize=5.9, ha="center")
    ax.text(0.94, 0.62, "H / HSurf / H$_2$(x)", fontsize=5.65, ha="center")
    ax.text(0.55, 0.09, "ER = Eley–Rideal; LH = Langmuir–Hinshelwood. Arrows show declared reaction classes,\nnot their rates, probabilities, or net fluxes.", fontsize=6.25, ha="center", style="italic")
    ax.text(-0.07, 1.06, "A", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    ax = fig.add_subplot(grid[0, 1]); ax.axis("off")
    ax.set_title("Five explicit routes create surface NH (first N–H bond)", loc="left", pad=8)
    rows = [[entry["line_number"], entry["reaction"].replace("Surf", "*").replace("  ", " "),
             entry["first_nh_class"].replace("Langmuir–Hinshelwood: ", "LH: ").replace("Eley–Rideal: ", "ER: ")]
            for entry in first]
    table = ax.table(cellText=rows, colLabels=["Source\nline", "Declared reaction", "Class"], cellLoc="left", colLoc="left",
                     colWidths=[0.13, 0.51, 0.36], bbox=[0.0, 0.15, 1.0, 0.70])
    table.auto_set_font_size(False); table.set_fontsize(6.2)
    for (row, _col), cell in table.get_celld().items():
        cell.set_linewidth(0.35)
        if row == 0: cell.set_facecolor("#DDEBF7"); cell.set_text_props(weight="bold")
        elif row % 2 == 0: cell.set_facecolor("#F7F7F7")
    ax.text(0.0, 0.055, "These are structural competitors for the first surface N–H bond. Their ranking\nrequires a field-transferred, boundary-harmonized rate-of-progress calculation.", fontsize=6.35, style="italic")
    ax.text(-0.13, 1.06, "B", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    ax = fig.add_subplot(grid[1, 0])
    family_counts = Counter(row["family"] for row in ledger)
    order = ["N2 dissociative adsorption", "H2 dissociative adsorption", "first surface N–H formation",
             "NHxSurf → NH3 release", "N adsorption", "N removal / N2 formation", "NHSurf → NH2Surf",
             "NH adsorption / hydrogenation", "NH2 adsorption / hydrogenation", "H adsorption", "H recombination / H2 formation"]
    labels = [name for name in order if name in family_counts]
    display = {"N2 dissociative adsorption": "N$_2$(x) dissociative adsorption", "H2 dissociative adsorption": "H$_2$(x) dissociative adsorption",
               "first surface N–H formation": "first surface N–H formation", "NHxSurf → NH3 release": "NH$_x$Surf → NH$_3$ release",
               "NHSurf → NH2Surf": "NHSurf → NH$_2$Surf", "NH adsorption / hydrogenation": "NH adsorption / hydrogenation",
               "NH2 adsorption / hydrogenation": "NH$_2$ adsorption / hydrogenation", "H recombination / H2 formation": "H recombination / H$_2$ formation"}
    colors = ["#9DC3E6" if "adsorption" in item else "#ED7D31" if "N–H" in item else "#A9D18E" if "NH" in item else "#B7B7B7" for item in labels]
    ax.barh([display.get(item, item) for item in labels[::-1]], [family_counts[item] for item in labels[::-1]],
            color=colors[::-1], edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Declared reaction records (not flux)")
    ax.set_title("Surface-branch reaction-family census")
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(-0.13, 1.07, "C", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    ax = fig.add_subplot(grid[1, 1]); ax.axis("off")
    species = [node for node, data in graph.nodes(data=True) if data["kind"] == "species"]
    nontrivial_scc = [members for members in nx.strongly_connected_components(graph) if len(members) > 1]
    items = [
        ("Frozen input", "46 explicit Fe(110) surface reactions;\nSHA-256 retained"),
        ("Surface N–H entry", "5 declared routes:\n4 ER and 1 LH"),
        ("Topology", f"{len(species)} species nodes; {len(nontrivial_scc)} nontrivial SCCs\nin the directed reaction bipartite graph"),
        ("What topology\ncan say", "which molecular routes are encoded\nand structurally available"),
        ("What it cannot\nsay", "which route dominates, or whether surface\novertakes gas-phase NH entry"),
        ("Required next\ngate", "calibrated F$_s$=βE$_{bulk}$ + shared CSTR,\nsite balance and ROP fractions"),
    ]
    ax.set_title("Interpretation boundary: a topology result, not a flux result", loc="left", pad=8)
    y = 0.87
    for label, body in items:
        ax.add_patch(FancyBboxPatch((0.02, y - 0.10), 0.96, 0.105, boxstyle="round,pad=0.009", facecolor="#F5F6F7", edgecolor="#9E9E9E", lw=0.45))
        ax.text(0.05, y - 0.045, label, fontweight="bold", fontsize=6.35, va="center")
        ax.text(0.35, y - 0.045, body, fontsize=6.15, va="center")
        y -= 0.135
    ax.text(-0.13, 1.06, "D", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    fig.text(0.5, 0.018, "Figure 23. Graph-theoretic audit of the explicit Fe(110) surface branch in the Shao–Mesbah DFT–microkinetic input.\nIt records topology and first-N–H-entry alternatives only; it does not report surface fluxes or a gas–surface pathway ranking.", ha="center", fontsize=7.5)
    fig.savefig(output / "Figure_23_Fe110_surface_NH_entry_graph.png", dpi=600, bbox_inches="tight")
    fig.savefig(output / "Figure_23_Fe110_surface_NH_entry_graph.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(output: Path, rows: list[dict], first: list[dict], graph: nx.DiGraph) -> None:
    species = [node for node, data in graph.nodes(data=True) if data["kind"] == "species"]
    nontrivial_scc = [members for members in nx.strongly_connected_components(graph) if len(members) > 1]
    lines = [
        "# Fe(110) surface-N–H-entry graph audit",
        "",
        "**Status:** source-locked, topology-only audit. No surface rate-of-progress, hybrid CSTR, or gas–surface flux fraction is reported.",
        "",
        "## Frozen source",
        "",
        f"- Input: `{SOURCE}`",
        f"- SHA-256: `{sha256(SOURCE)}`",
        f"- Parsed explicit surface reactions: {len(rows)}",
        "",
        "## Structural result",
        "",
        "Five reactions in the declared Fe(110) subnetwork form `NHSurf` while joining an N-bearing and an H-bearing precursor. Four are Eley–Rideal alternatives: `N + HSurf`, `N(2D) + HSurf`, `N(2P) + HSurf`, and `H + NSurf`. One is Langmuir–Hinshelwood: `NSurf + HSurf -> NHSurf + Surf`.",
        "",
        "The network also encodes 13 N2-state dissociative-adsorption entries and 8 H2-state dissociative-adsorption entries. These multiplicities are reaction-record counts, not probabilities or flux shares.",
        "",
        "## Graph checks",
        "",
        f"- Bipartite graph species nodes: {len(species)}",
        f"- Reaction nodes: {len(rows)}",
        f"- Nontrivial strongly connected components: {len(nontrivial_scc)}",
        "- The absence of a nontrivial directed SCC merely reflects the forward reaction records selected in this input; it must not be interpreted as a proof that physical reverse chemistry is absent.",
        "",
        "## Interpretation guardrail",
        "",
        "This figure expands the manuscript’s mechanism discussion by distinguishing Eley–Rideal and Langmuir–Hinshelwood structural alternatives at the surface first-N–H step. It cannot establish which alternative is kinetically dominant, whether the surface contribution exceeds the Hong gas-phase direct-NH contribution, or whether a real device realizes the Shao local field. Those claims remain prohibited until the Figure 22 field-transfer, shared-boundary and numerical-acceptance gates are passed.",
        "",
        "## Deliverables",
        "",
        "- `Figure_23_Fe110_surface_NH_entry_graph.png/.pdf`",
        "- `P14_Fe110_surface_reaction_ledger.csv`",
        "- `P14_first_surface_NH_entries.csv`",
        "- `P14_surface_graph_metrics.csv`",
    ]
    (output / "P14_Fe110_surface_graph_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Topology-only Fe(110) surface graph audit")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not SOURCE.is_file():
        raise FileNotFoundError(SOURCE)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    rows = load_surface_reactions(SOURCE)
    graph = build_bipartite(rows)
    ledger, first = write_ledgers(output, rows, graph)
    draw_figure(output, ledger, first, graph)
    write_report(output, rows, first, graph)
    print(f"surface_reactions={len(rows)}")
    print(f"first_surface_NH_entries={len(first)}")
    print(f"output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
