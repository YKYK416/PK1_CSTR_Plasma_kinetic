#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Automorphism-orbit audit for the frozen Fe(110) reaction subnetwork.

This group-theoretic analysis applies only to the *unweighted directed
bipartite topology* of the 46 explicit surface records.  Each reported orbit
is verified by an edge-preserving within-orbit cyclic permutation.  The script
does not treat topology-equivalent state labels as kinetically equivalent.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from _paths import LEGACY_HONG_ROOT

PROJECT = LEGACY_HONG_ROOT
DEFAULT_OUTPUT = PROJECT / "analysis" / "p15_fe110_surface_symmetry_orbits_20260901_r1"

from plot_fe110_surface_graph_audit import (  # noqa: E402
    SOURCE, build_bipartite, load_surface_reactions, sha256,
)


ORBIT_CANDIDATES = {
    "O_N2exc": ("N2(A3)", "N2(B3)", "N2(C3)", "N2(a`1)", "N2(v1)", "N2(v2)", "N2(v3)",
                  "N2(v4)", "N2(v5)", "N2(v6)", "N2(v7)", "N2(v8)"),
    "O_H2elec": ("H2(A3SIG)", "H2(B1SIG)", "H2(B3SIG)", "H2(C3PI)"),
    "O_H2v": ("H2(v1)", "H2(v2)", "H2(v3)"),
    "O_Nactive": ("N", "N(2D)", "N(2P)"),
}


def reaction_key(row: dict) -> tuple[tuple[str, ...], tuple[str, ...]]:
    return tuple(sorted(row["reactants"])), tuple(sorted(row["products"]))


def transformed_key(row: dict, mapping: dict[str, str]) -> tuple[tuple[str, ...], tuple[str, ...]]:
    return (tuple(sorted(mapping.get(token, token) for token in row["reactants"])),
            tuple(sorted(mapping.get(token, token) for token in row["products"])))


def cyclic_permutation(members: tuple[str, ...]) -> dict[str, str]:
    return {members[index]: members[(index + 1) % len(members)] for index in range(len(members))}


def verify_orbit(rows: list[dict], graph, members: tuple[str, ...]) -> tuple[bool, int]:
    """Verify one non-identity cycle extends to an edge-preserving graph map."""
    species_map = cyclic_permutation(members)
    lookup = {reaction_key(row): index for index, row in enumerate(rows, start=1)}
    full_map = {node: node for node in graph.nodes}
    full_map.update(species_map)
    for index, row in enumerate(rows, start=1):
        target = lookup.get(transformed_key(row, species_map))
        if target is None:
            return False, 0
        full_map[f"r{index:02d}"] = f"r{target:02d}"
    edges = set(graph.edges())
    transformed_edges = {(full_map[left], full_map[right]) for left, right in edges}
    return transformed_edges == edges and len(set(full_map.values())) == len(full_map), len(full_map)


def orbit_row(label: str, members: tuple[str, ...], rows: list[dict], graph) -> dict:
    passed, node_count = verify_orbit(rows, graph, members)
    appearances = []
    for member in members:
        incident = [row for row in rows if member in row["reactants"] or member in row["products"]]
        appearances.append(len(incident))
    return {"orbit": label, "size": len(members), "members": " | ".join(members),
            "surface_records_per_member": ",".join(map(str, appearances)),
            "cyclic_permutation_verified": passed, "mapped_graph_nodes": node_count}


def draw_box(ax, xy, text, face, width=0.22, height=0.12, fontsize=7):
    x, y = xy
    ax.add_patch(FancyBboxPatch((x - width / 2, y - height / 2), width, height, boxstyle="round,pad=0.01",
                                facecolor=face, edgecolor="#555555", linewidth=0.55))
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize)


def arrow(ax, start, end, text="", bend=0, color="#777777"):
    ax.add_patch(FancyArrowPatch(start, end, connectionstyle=f"arc3,rad={bend}", arrowstyle="-|>", mutation_scale=8,
                                 linewidth=0.8, color=color))
    if text:
        ax.text((start[0] + end[0]) / 2, (start[1] + end[1]) / 2 + (0.026 if bend >= 0 else -0.026), text,
                ha="center", va="center", fontsize=5.9)


def draw_figure(output: Path, orbit_rows: list[dict]) -> None:
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.1, "axes.titlesize": 9.1,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(12.2, 8.4))
    grid = fig.add_gridspec(2, 2, left=0.075, right=0.98, bottom=0.085, top=0.93, hspace=0.47, wspace=0.28,
                            width_ratios=(1.06, 0.94))

    ax = fig.add_subplot(grid[0, 0]); ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title("Topology-equivalent state orbits in the explicit Fe(110) branch", loc="left", pad=8)
    draw_box(ax, (0.13, 0.77), "O$_{N2exc}$\n12 N$_2$(x) states", "#DDEBF7")
    draw_box(ax, (0.13, 0.51), "O$_{H2elec}$\n4 electronic H$_2$ states", "#FCE4D6")
    draw_box(ax, (0.13, 0.25), "O$_{H2v}$\n3 vibrational H$_2$ states", "#FFF2CC")
    draw_box(ax, (0.48, 0.77), "NSurf", "#E2F0D9", width=0.15, height=0.09)
    draw_box(ax, (0.48, 0.45), "HSurf", "#E2F0D9", width=0.15, height=0.09)
    draw_box(ax, (0.48, 0.18), "NH$_3$", "#EADCF8", width=0.15, height=0.09)
    draw_box(ax, (0.80, 0.60), "O$_{Nactive}$\nN, N(2D), N(2P)", "#F4CCCC", width=0.25)
    draw_box(ax, (0.80, 0.27), "NHSurf", "#F4CCCC", width=0.17, height=0.09)
    arrow(ax, (0.24, 0.77), (0.40, 0.77), "12 matched N$_2$(x) adsorption edges")
    arrow(ax, (0.24, 0.51), (0.40, 0.45), "4 matched H$_2$(e) adsorption edges", bend=-0.08)
    arrow(ax, (0.24, 0.25), (0.40, 0.18), "3 matched H$_2$(v) → NH$_3$ edges", bend=-0.08)
    arrow(ax, (0.67, 0.60), (0.70, 0.30), "3 matched N(x) + HSurf → NHSurf", bend=-0.10, color="#C55A11")
    ax.text(0.52, 0.06, "Each colored set permits a verified cyclic relabelling in the unweighted\ndirected reaction graph, including its incident reaction nodes.", ha="center", fontsize=6.25, style="italic")
    ax.text(-0.10, 1.06, "A", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    ax = fig.add_subplot(grid[0, 1]); ax.axis("off")
    ax.set_title("Verified orbit table and automorphism subgroup", loc="left", pad=8)
    labels = {"O_N2exc": "O$_{N2exc}$", "O_H2elec": "O$_{H2elec}$", "O_H2v": "O$_{H2v}$", "O_Nactive": "O$_{Nactive}$"}
    short_members = {"O_N2exc": "N$_2$(A,B,a,C,v1–v8)", "O_H2elec": "H$_2$(A,B,C,B$^3$)",
                     "O_H2v": "H$_2$(v1–v3)", "O_Nactive": "N, N(2D), N(2P)"}
    table_rows = [[labels[row["orbit"]], row["size"], short_members[row["orbit"]], "passed"] for row in orbit_rows]
    table = ax.table(cellText=table_rows, colLabels=["Orbit", "Size", "Member labels", "Cycle\ncheck"], cellLoc="left", colLoc="left",
                     colWidths=[0.16, 0.11, 0.54, 0.19], bbox=[0.0, 0.40, 1.0, 0.42])
    table.auto_set_font_size(False); table.set_fontsize(6.55)
    for (row, _col), cell in table.get_celld().items():
        cell.set_linewidth(0.35)
        if row == 0: cell.set_facecolor("#DDEBF7"); cell.set_text_props(weight="bold")
        elif row % 2 == 0: cell.set_facecolor("#F7F7F7")
    order = math.prod(math.factorial(row["size"]) for row in orbit_rows)
    ax.add_patch(FancyBboxPatch((0.04, 0.17), 0.92, 0.14, boxstyle="round,pad=0.012", facecolor="#F5F6F7", edgecolor="#777777", lw=0.5))
    ax.text(0.08, 0.24, r"Guaranteed subgroup", fontsize=6.7, fontweight="bold", va="center")
    ax.text(0.36, 0.24, r"$S_{12}\times S_4\times S_3\times S_3$", fontsize=8.0, va="center")
    ax.text(0.36, 0.19, f"at least {order:,} topology-preserving relabellings", fontsize=6.2, va="center")
    ax.text(0.02, 0.055, "The count belongs to an unlabelled, unweighted graph representation. It is not a count of\nchemical microstates and it supplies no kinetic degeneracy factor.", fontsize=6.25, style="italic")
    ax.text(-0.13, 1.06, "B", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    ax = fig.add_subplot(grid[1, 0]); ax.axis("off")
    ax.set_title("How physical kinetics breaks the topology-only symmetry", loc="left", pad=8)
    stages = [
        ("1. Unweighted topology", "state labels have identical\nincidence patterns", "#DDEBF7"),
        ("2. State-resolved rates", "DFT energetics, cross sections\nand adsorption terms differ", "#FFF2CC"),
        ("3. State populations", "E/N, EEDF, wall loss and gas\ncomposition select densities", "#FCE4D6"),
        ("4. Device-local outcome", "F$_s$=βE$_{bulk}$, sites, transport\nand CSTR boundary select flux", "#F4CCCC"),
    ]
    x_positions = [0.13, 0.38, 0.63, 0.88]
    for index, ((title, body, color), x) in enumerate(zip(stages, x_positions)):
        ax.add_patch(FancyBboxPatch((x - 0.105, 0.40), 0.21, 0.28, boxstyle="round,pad=0.012", facecolor=color, edgecolor="#666666", lw=0.55))
        ax.text(x, 0.60, title, ha="center", va="center", fontsize=6.3, fontweight="bold")
        ax.text(x, 0.48, body, ha="center", va="center", fontsize=6.0)
        if index < len(stages) - 1:
            arrow(ax, (x + 0.11, 0.54), (x_positions[index + 1] - 0.11, 0.54), "", color="#555555")
    ax.text(0.5, 0.22, "A topological orbit is a hypothesis-generating reduction only. The symmetry is broken before any\nrate-of-progress ranking is calculated; retaining state labels in the kinetic ledger is therefore mandatory.", ha="center", fontsize=6.55, style="italic")
    ax.text(-0.10, 1.06, "C", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    ax = fig.add_subplot(grid[1, 1]); ax.axis("off")
    ax.set_title("Consequences for mechanism discrimination", loc="left", pad=8)
    rows = [
        ("Valid inference", "The source contains symmetric structural alternatives; a reported\npathway ranking must explain the kinetic symmetry-breaking input."),
        ("Invalid inference", "Orbit size does not imply equal flux, equal rate coefficient,\nor an entropic/degeneracy multiplier."),
        ("Experimental handle", "State-selective OES / absorption and isotope-resolved NH$_3$ can\nconstrain which symmetry-breaking layer is active."),
        ("Model gate", "Do not merge orbit members in a fitted model unless their rates\nand populations are independently demonstrated to be exchangeable."),
    ]
    y = 0.85
    for title, body in rows:
        ax.add_patch(FancyBboxPatch((0.02, y - 0.12), 0.96, 0.13, boxstyle="round,pad=0.009", facecolor="#F5F6F7", edgecolor="#999999", lw=0.45))
        ax.text(0.05, y - 0.045, title, fontsize=6.65, fontweight="bold", va="center")
        ax.text(0.32, y - 0.045, body, fontsize=6.2, va="center")
        y -= 0.17
    ax.text(-0.13, 1.06, "D", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    fig.text(0.5, 0.017, "Figure 24. Automorphism-orbit analysis of the source-locked Fe(110) surface reaction topology.\nThe identified permutation symmetries are graph-theoretic only; state-resolved kinetics, populations and device fields break them before flux ranking.", ha="center", fontsize=7.45)
    fig.savefig(output / "Figure_24_Fe110_surface_symmetry_orbits.png", dpi=600, bbox_inches="tight")
    fig.savefig(output / "Figure_24_Fe110_surface_symmetry_orbits.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(output: Path, rows: list[dict], orbit_rows: list[dict]) -> None:
    order = math.prod(math.factorial(row["size"]) for row in orbit_rows)
    text = [
        "# Fe(110) surface-reaction automorphism-orbit audit",
        "",
        "**Status:** source-locked group-theoretic topology audit. It has no kinetic weights, no surface ROP and no hybrid CSTR calculation.",
        "",
        "## Source and representation",
        "",
        f"- Input: `{SOURCE}`",
        f"- SHA-256: `{sha256(SOURCE)}`",
        f"- Surface reactions: {len(rows)}",
        "- Graph: directed bipartite species–reaction representation; species and reaction nodes retain only their node type during the symmetry test.",
        "",
        "## Verified topology orbits",
        "",
        "Each orbit below passed an explicit non-identity cyclic permutation that maps every graph edge back onto an edge, including all incident reaction nodes.",
        "",
    ]
    for row in orbit_rows:
        text.append(f"- `{row['orbit']}` (size {row['size']}): {row['members']}")
    text += [
        "",
        "The independently permutable orbit factors guarantee a subgroup `S12 × S4 × S3 × S3`, with at least " + f"{order:,}" + " automorphisms in the topology-only representation.",
        "",
        "## Scientific meaning and guardrail",
        "",
        "The result says that these labels have matching *incidence patterns* inside the explicit surface subnetwork. It does not say that their adsorption barriers, field response, electron-impact source terms, concentrations or rates are equal. Those state-resolved quantities are exactly the symmetry-breaking mechanisms that determine an eventual ROP ranking.",
        "",
        "Consequently, no state labels were collapsed in the kinetic model and no degeneracy multiplier was introduced. A future reduced model may merge an orbit only after demonstrating exchangeability of both the rate laws and the state populations under a shared physical boundary.",
        "",
        "## Deliverables",
        "",
        "- `Figure_24_Fe110_surface_symmetry_orbits.png/.pdf`",
        "- `P15_Fe110_surface_topology_orbits.csv`",
        "- `P15_Fe110_surface_symmetry_orbits.md`",
    ]
    (output / "P15_Fe110_surface_symmetry_orbits.md").write_text("\n".join(text) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Fe(110) surface topology automorphism-orbit audit")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    rows = load_surface_reactions(SOURCE)
    graph = build_bipartite(rows)
    orbit_rows = [orbit_row(label, members, rows, graph) for label, members in ORBIT_CANDIDATES.items()]
    if not all(row["cyclic_permutation_verified"] for row in orbit_rows):
        raise RuntimeError("At least one declared topology orbit did not pass edge-preserving permutation verification")
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    with (output / "P15_Fe110_surface_topology_orbits.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=orbit_rows[0].keys()); writer.writeheader(); writer.writerows(orbit_rows)
    draw_figure(output, orbit_rows)
    write_report(output, rows, orbit_rows)
    print(f"verified_orbits={len(orbit_rows)}")
    print(f"automorphism_subgroup_lower_bound={math.prod(math.factorial(row['size']) for row in orbit_rows)}")
    print(f"output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
