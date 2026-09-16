#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stoichiometry-preserving ROP reaction-hypergraph utilities.

The raw ZDPlasKin ROP output is converted to a directed species--reaction
bipartite graph.  Reaction nodes retain every reactant and product, avoiding
the chemically invalid shortcut of a simple species-only graph.  A separate
nitrogen-lineage projection is supplied solely for graph metrics; it is
explicitly derived from N-atom stoichiometry and does not replace the
hypergraph in reaction-path visualizations.
"""
from __future__ import annotations

import csv
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import networkx as nx

from _paths import LEGACY_HONG_ROOT

PROJECT = LEGACY_HONG_ROOT
RUN_ROOT = (PROJECT / "Reproduction" / "2017Hong" / "方案2_脉冲能效地图" /
            "build_pulse" / "runs" / "p4_dense_cstr_20260829")
SUMMARY = PROJECT / "analysis" / "figures_20260829" / "Figure_3_dense_EN_N2_source_data.csv"

DIRECT_NH = {
    "named electronic H2*": re.compile(r"^N\+H2\((B3SIG|B1SIG|C3PI|A3SIG)\)=>H\+NH$"),
    "Rydberg-H2*": re.compile(r"^N\+H2\(RYDBERG_SUM\)=>H\+NH$"),
    "N(2D)+H2": re.compile(r"^N\(2D\)\+H2=>H\+NH$"),
    "N(2P)+H2": re.compile(r"^N\(2P\)\+H2=>H\+NH$"),
    "vibrational H2": re.compile(r"^N\+H2\(V[123]\)=>H\+NH$"),
    "H+N+M": re.compile(r"^H\+N\+(N2|H2)=>NH\+\1$"),
}


@dataclass(frozen=True)
class ReactionFlux:
    reaction: str
    reactants: tuple[str, ...]
    products: tuple[str, ...]
    phi_cm3: float
    route: str
    n_balanced: bool
    h_balanced: bool


def case_tag(en_td: int, n2_fraction: float) -> str:
    return f"EN{int(en_td)}_n2{n2_fraction:.1f}"


def split_side(side: str) -> tuple[str, ...]:
    """Split on + outside state-label parentheses."""
    parts, buf, depth = [], [], 0
    for char in side.strip():
        if char == "(":
            depth += 1
        elif char == ")":
            depth = max(0, depth - 1)
        if char == "+" and depth == 0:
            token = "".join(buf).strip()
            if token:
                parts.append(token)
            buf = []
        else:
            buf.append(char)
    token = "".join(buf).strip()
    if token:
        parts.append(token)
    return tuple(parts)


def parse_reaction(text: str) -> tuple[tuple[str, ...], tuple[str, ...]] | None:
    """Parse a ZDPlasKin ROP label without flattening its stoichiometry."""
    text = (text or "").strip()
    if ":" in text:
        text = text.split(":", 1)[1].strip()
    separator = "=>" if "=>" in text else ("->" if "->" in text else None)
    if separator is None:
        return None
    left, right = text.split(separator, 1)
    reactants, products = split_side(left), split_side(right)
    return (reactants, products) if reactants and products else None


def molecular_formula(species: str) -> str:
    """Return the neutral formula stem, keeping only conventional atoms."""
    value = species.strip().replace("*", "")
    value = re.sub(r"\([^)]*\)", "", value)
    value = value.replace("+", "").replace("-", "")
    value = re.sub(r"^\d+", "", value)
    return value


def atom_count(species: str, atom: str) -> int:
    coefficient_match = re.match(r"^\s*(\d+)", species)
    coefficient = int(coefficient_match.group(1)) if coefficient_match else 1
    formula = molecular_formula(species)
    total = 0
    for symbol, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if symbol == atom:
            total += int(count or "1")
    return coefficient * total


def atom_total(species: Iterable[str], atom: str) -> int:
    return sum(atom_count(item, atom) for item in species)


def classify_direct_nh(reaction: str) -> str:
    return next((label for label, pattern in DIRECT_NH.items() if pattern.search(reaction)), "other")


def terminal_window(cycles_csv: Path, keep_cycles: int = 3) -> tuple[int, int]:
    with cycles_csv.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Empty cycle file: {cycles_csv}")
    last = rows[-1]
    last_cycle = int(last["cycle"])
    if int(float(last["stable_streak"])) < keep_cycles:
        raise ValueError(f"P0 not established in {cycles_csv}")
    return max(1, last_cycle - keep_cycles + 1), last_cycle


def integrate_case(en_td: int, n2_fraction: float, run_root: Path = RUN_ROOT) -> tuple[list[ReactionFlux], dict]:
    """Integrate positive ROP over exactly the final three validated P0 cycles."""
    tag = case_tag(en_td, n2_fraction)
    case = run_root / tag
    cycles = case / f"pulse_cycles_{tag}.csv"
    rates = case / f"pulse_rates_{tag}.csv"
    first, last = terminal_window(cycles)
    total: Counter[str] = Counter()
    parsed: dict[str, tuple[tuple[str, ...], tuple[str, ...]]] = {}
    samples = 0
    with rates.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        expected = {"cycle", "reaction", "rate_cm-3s-1", "dt_s"}
        if not expected.issubset(reader.fieldnames or set()):
            raise ValueError(f"Missing ROP columns in {rates}")
        for row in reader:
            try:
                if int(row["cycle"]) < first:
                    continue
                rate, dt = float(row["rate_cm-3s-1"]), float(row["dt_s"])
            except (TypeError, ValueError):
                continue
            if not (rate > 0.0 and dt > 0.0):
                continue
            reaction = (row["reaction"] or "").strip()
            if reaction not in parsed:
                parsed_value = parse_reaction(reaction)
                if parsed_value is None:
                    continue
                parsed[reaction] = parsed_value
            total[reaction] += rate * dt
            samples += 1
    if not total:
        raise ValueError(f"No positive terminal ROP data in {rates}")
    records = []
    for reaction, phi in total.items():
        reactants, products = parsed[reaction]
        records.append(ReactionFlux(
            reaction=reaction, reactants=reactants, products=products, phi_cm3=float(phi),
            route=classify_direct_nh(reaction),
            n_balanced=atom_total(reactants, "N") == atom_total(products, "N"),
            h_balanced=atom_total(reactants, "H") == atom_total(products, "H"),
        ))
    records.sort(key=lambda item: item.phi_cm3, reverse=True)
    meta = {"en_td": int(en_td), "n2_fraction": float(n2_fraction), "tag": tag,
            "cycle_first": first, "cycle_last": last, "positive_samples": samples,
            "rates_csv": str(rates)}
    return records, meta


def bipartite_hypergraph(records: Iterable[ReactionFlux]) -> nx.DiGraph:
    """Species--reaction incidence graph that retains full reaction arity."""
    graph = nx.DiGraph()
    for index, record in enumerate(records):
        rnode = f"R:{index}"
        graph.add_node(rnode, kind="reaction", reaction=record.reaction,
                       phi_cm3=record.phi_cm3, route=record.route)
        for species in record.reactants:
            snode = f"S:{species}"
            graph.add_node(snode, kind="species", label=species)
            graph.add_edge(snode, rnode, role="reactant", capacity=record.phi_cm3)
        for species in record.products:
            snode = f"S:{species}"
            graph.add_node(snode, kind="species", label=species)
            graph.add_edge(rnode, snode, role="product", capacity=record.phi_cm3)
    return graph


def nitrogen_lineage_graph(records: Iterable[ReactionFlux]) -> nx.DiGraph:
    """N-atom-flow projection used only for centrality metrics.

    Each hyperreaction flux is split among N-bearing reactants and products in
    proportion to their N-atom stoichiometric coefficients.  Co-reactant
    requirements remain visible in the parent hypergraph.
    """
    graph = nx.DiGraph()
    for record in records:
        n_r = npairs(record.reactants, "N")
        n_p = npairs(record.products, "N")
        if not n_r or not n_p:
            continue
        total_r, total_p = sum(count for _, count in n_r), sum(count for _, count in n_p)
        transported = record.phi_cm3 * min(total_r, total_p)
        if transported <= 0:
            continue
        for src, count_r in n_r:
            for dst, count_p in n_p:
                contribution = transported * (count_r / total_r) * (count_p / total_p)
                previous = graph.get_edge_data(src, dst, {}).get("flux_cm3", 0.0)
                graph.add_edge(src, dst, flux_cm3=previous + contribution)
    for src, dst, attrs in graph.edges(data=True):
        attrs["cost"] = -float(__import__("math").log(max(attrs["flux_cm3"], 1e-300)))
    return graph


def npairs(species: Iterable[str], atom: str) -> list[tuple[str, int]]:
    return [(item, atom_count(item, atom)) for item in species if atom_count(item, atom) > 0]


def nh_relevant(record: ReactionFlux) -> bool:
    labels = set(record.reactants) | set(record.products)
    stems = {molecular_formula(item) for item in labels}
    return bool(stems & {"N", "N2", "NH", "NH2", "NH3"})


def display_reactions(records: list[ReactionFlux], max_other: int = 10) -> list[ReactionFlux]:
    """Select a readable NHx subgraph, while always retaining all direct NH sources."""
    forced = [item for item in records if item.route != "other"]
    nhx = [item for item in records if any(molecular_formula(s) in {"NH", "NH2", "NH3"}
                                           for s in item.reactants + item.products)]
    other = [item for item in records if nh_relevant(item) and item not in forced and item not in nhx]
    chosen = {item.reaction: item for item in forced + nhx + other[:max_other]}
    return sorted(chosen.values(), key=lambda item: item.phi_cm3, reverse=True)


def validate_records(records: Iterable[ReactionFlux]) -> dict[str, int]:
    rows = list(records)
    return {
        "reactions": len(rows),
        "n_unbalanced": sum(not row.n_balanced for row in rows),
        "h_unbalanced": sum(not row.h_balanced for row in rows),
        "direct_nh_reactions": sum(row.route != "other" for row in rows),
    }
