#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Audit the corrected-Hong mechanism in a cross-mechanism reaction-family schema.

This script does not claim that a literature mechanism is numerically comparable
until its full kinetic input is locally available and validated.  It inventories
the local corrected-Hong mechanism, constructs stoichiometry-preserving
reaction-family records, and records external candidate mechanisms as
metadata only.  The resulting quotient representation supports
future mechanism-ensemble comparison without conflating reaction labels with
equivalent physics.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import pandas as pd

from _paths import LEGACY_HONG_ROOT
from hong_reaction_hypergraph import parse_reaction


PROJECT = LEGACY_HONG_ROOT
INPUT = PROJECT / "Reproduction" / "2017Hong" / "literature_extract_2017_2018_corrected" / "kinet.inp"
EXTERNAL_FE110_INPUT = (PROJECT / "Reproduction" / "external_models" / "DFT-microkinetic" /
                       "Model_SA_Const_Entropy_base" /
                       "kinet_varyT_metal_auto_sens_DFT_in_entropy_verying_basis.txt")
REPRESENTATIVE_FLUX = PROJECT / "analysis" / "figures_20260829" / "Figure_6_reaction_flux_table.csv"
DEFAULT_OUTPUT = PROJECT / "analysis" / "p11_mechanism_family_audit_20260831"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_reaction(line: str) -> str | None:
    stripped = line.strip()
    if not stripped or stripped.startswith("#") or "=>" not in stripped:
        return None
    return stripped.split("!", 1)[0].strip()


def compact(side: tuple[str, ...]) -> str:
    return "+".join(item.replace(" ", "") for item in side)


def has(side: tuple[str, ...], needle: str) -> bool:
    return any(item.replace(" ", "") == needle for item in side)


def has_prefix(side: tuple[str, ...], prefix: str) -> bool:
    return any(item.replace(" ", "").startswith(prefix) for item in side)


def classify(reactants: tuple[str, ...], products: tuple[str, ...]) -> tuple[str, str]:
    """Return pathway family and layer; direct-NH tests are intentionally first."""
    left, right = compact(reactants), compact(products)
    named = r"H2\((B3SIG|B1SIG|C3PI|A3SIG)\)"
    if re.fullmatch(rf"N\+{named}", left) and right == "H+NH":
        return "Direct NH entry: named electronic H2*", "NH entry"
    if left == "N+H2(RYDBERG_SUM)" and right == "H+NH":
        return "Direct NH entry: Rydberg H2*", "NH entry"
    if re.fullmatch(r"N\+H2\([Vv][0-9]+\)", left) and right == "H+NH":
        return "Direct NH entry: vibrational H2", "NH entry"
    if re.fullmatch(r"N\(2[DP]\)\+H2", left) and right == "H+NH":
        return "Direct NH entry: excited N + H2", "NH entry"
    if left in {"H+N+N2", "H+N+H2", "H+N+@M"} and right in {"NH+N2", "NH+H2", "NH+@M"}:
        return "Direct NH entry: termolecular association", "NH entry"

    if any("Surf" in item or item == "@S" for item in reactants + products):
        return "Surface reaction network", "surface"

    if "NH" in left and "NH2" in right and "NH3" not in right:
        return "NH -> NH2 turnover", "NHx turnover"
    if ("NH" in left or "NH2" in left) and "NH3" in right:
        return "NH/NH2 -> NH3 turnover", "NHx turnover"
    if "NH3" in left and "NH3" not in right:
        return "NH3 loss / ionic conversion", "NHx loss"
    if "NH" in left and "NH" not in right:
        return "NH loss / reverse conversion", "NHx loss"
    if "NH2" in left and "NH2" not in right and "NH3" not in right:
        return "NH2 loss / reverse conversion", "NHx loss"

    if has(reactants, "e") and has(reactants, "N2") and has_prefix(products, "N2("):
        if any("(V" in item for item in products):
            return "N2 vibrational activation", "activation"
        return "N2 electronic activation", "activation"
    if has(reactants, "e") and has(reactants, "H2") and has_prefix(products, "H2("):
        if any("(V" in item for item in products):
            return "H2 vibrational activation", "activation"
        return "H2 electronic activation", "activation"
    if has(reactants, "e") and has(reactants, "N2") and (has(products, "N") or has(products, "N(2D)") or has(products, "N(2P)")):
        return "N2 dissociation / active-N formation", "activation"
    if has(reactants, "e") and has(reactants, "H2") and has(products, "H"):
        return "H2 dissociation / active-H formation", "activation"
    if "WALL" in left or "WALL" in right or (len(reactants) == 1 and len(products) == 1 and ("(" in left or "(" in right)):
        return "State relaxation / wall process", "relaxation"
    if any("^" in item or item.endswith("+") or item.endswith("-") for item in reactants + products):
        return "Ion / electron network", "charged-particle"
    return "Other neutral gas-phase network", "other"


def canonical_species(item: str) -> str:
    value = item.replace(" ", "")
    value = re.sub(r"H2\((B3SIG|B1SIG|C3PI|A3SIG)\)", "H2(electronic)", value)
    value = value.replace("H2(RYDBERG_SUM)", "H2(Rydberg)")
    value = re.sub(r"H2\([Vv][0-9]+\)", "H2(v)", value)
    value = re.sub(r"N2\([Vv][0-9]+\)", "N2(v)", value)
    value = re.sub(r"N2\((A3|B3|a`1|C3)\)", "N2(electronic)", value)
    value = re.sub(r"N\(2[DP]\)", "N(excited)", value)
    return value


def quotient_skeleton(reactants: tuple[str, ...], products: tuple[str, ...]) -> str:
    return f"{' + '.join(canonical_species(x) for x in reactants)} -> {' + '.join(canonical_species(x) for x in products)}"


def read_ledger(path: Path) -> pd.DataFrame:
    rows = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1):
        reaction = clean_reaction(raw)
        if reaction is None:
            continue
        parsed = parse_reaction(reaction)
        if parsed is None:
            continue
        reactants, products = parsed
        family, layer = classify(reactants, products)
        rows.append({
            "mechanism_id": "hong_2017_2018_corrected_local",
            "line_number": line_number,
            "reaction": reaction,
            "reactants": " + ".join(reactants),
            "products": " + ".join(products),
            "family": family,
            "layer": layer,
            "quotient_skeleton": quotient_skeleton(reactants, products),
        })
    if not rows:
        raise RuntimeError(f"No active reactions found in {path}")
    return pd.DataFrame(rows)


def mechanism_registry(input_hash: str) -> pd.DataFrame:
    return pd.DataFrame([
        {
            "mechanism_id": "hong_2017_2018_corrected_local",
            "status": "computable_local",
            "kinetic_input": str(INPUT),
            "sha256": input_hash,
            "scope": "Local corrected Hong reconstruction; input contains a surface subnetwork, while current reported CW/CSTR representative ROP records have no surface rows",
            "comparison_role": "reference model; current results accepted only within stated CW/CSTR boundary",
            "source": "Local literature-traceable 2017/2018 corrected reconstruction",
        },
        {
            "mechanism_id": "elevated_pressure_acs_suschemeng_2025_candidate",
            "status": "not_ported",
            "kinetic_input": "",
            "sha256": "",
            "scope": "Experimental/kinetic study at elevated pressure; reports N2(A), N(2D), H2(v), NHx pathways",
            "comparison_role": "candidate alternative family structure; requires full rate file, provenance audit and harmonized boundary before simulation",
            "source": "doi:10.1021/acssuschemeng.5c06251",
        },
        {
            "mechanism_id": "pulsed_plasma_acs_suschemeng_2022_candidate",
            "status": "not_ported",
            "kinetic_input": "",
            "sha256": "",
            "scope": "Nanosecond-pulsed DBD with OES and kinetic modeling",
            "comparison_role": "candidate waveform/energy-partition comparator; requires full kinetic input and independent P0 validation",
            "source": "doi:10.1021/acssuschemeng.2c06259",
        },
        {
            "mechanism_id": "fe110_dft_microkinetic_jacsau_2024",
            "status": "validated_not_harmonized",
            "kinetic_input": str(EXTERNAL_FE110_INPUT),
            "sha256": sha256(EXTERNAL_FE110_INPUT) if EXTERNAL_FE110_INPUT.is_file() else "",
            "scope": "Public Fe(110) DFT–microkinetic surface branch with a Hong-derived gas-phase backbone; source input contains explicit surface adsorption, dissociation, Eley–Rideal and surface-hydrogenation reactions",
            "comparison_role": "surface-versus-gas competition branch; not an independent validation of the Hong gas-phase NH-entry reactions and not yet boundary-harmonized for numerical comparison",
            "source": "doi:10.1021/jacsau.3c00654; github.com/wwwccttoo/DFT-microkinetic",
        },
    ])


def representative_surface_flux_audit() -> pd.DataFrame:
    """Audit the existing three-state terminal ROP table without extrapolating to all cases."""
    if not REPRESENTATIVE_FLUX.is_file():
        return pd.DataFrame([{"source": str(REPRESENTATIVE_FLUX), "status": "missing", "total_rows": "",
                              "surface_rows": "", "conditions_with_surface_rows": ""}])
    frame = pd.read_csv(REPRESENTATIVE_FLUX)
    surface_mask = frame["reaction"].fillna("").str.contains(r"Surf|@S", regex=True)
    hits = frame[surface_mask]
    return pd.DataFrame([{
        "source": str(REPRESENTATIVE_FLUX),
        "status": "audited_representative_terminal_ROP",
        "total_rows": int(len(frame)),
        "surface_rows": int(len(hits)),
        "conditions_with_surface_rows": int(hits["condition"].nunique()),
    }])


def write_markdown(output: Path, ledger: pd.DataFrame, families: pd.DataFrame, registry: pd.DataFrame, surface_audit: pd.DataFrame, input_hash: str) -> None:
    direct = ledger[ledger["layer"] == "NH entry"]
    skeletons = direct.groupby(["family", "quotient_skeleton"], as_index=False).agg(
        member_reactions=("reaction", lambda values: " | ".join(values)),
        n_members=("reaction", "size"),
    )
    lines = [
        "# Mechanism-family and ensemble-readiness audit",
        "",
        "**Status:** computational inventory; no external candidate mechanism has been simulated.",
        "",
        "## Provenance",
        "",
        f"- Local input: `{INPUT}`",
        f"- SHA-256: `{input_hash}`",
        f"- Active reaction records parsed: {len(ledger)}",
        "- Disabled `#` lines are excluded. A reaction-family label is a comparison ontology, not a kinetic ranking.",
        "",
        "## Core guardrail",
        "",
        "The quotient skeleton maps state-resolved labels to reaction-family classes for cross-mechanism comparison. "
        "It identifies structural correspondence only. It does not assert equality of state energies, cross sections, "
        "rate coefficients, collision selection rules, or pathway fluxes.",
        "",
        "## Active reaction-family census",
        "",
        families.to_markdown(index=False),
        "",
        "## Surface-network boundary audit",
        "",
        f"The input contains {int((ledger['layer'] == 'surface').sum())} active surface-network records. "
        "This does not establish a surface contribution to the reported simulations. The audit below searches only the "
        "existing three representative terminal-ROP tables, not the complete 108-case scan.",
        "",
        surface_audit.to_markdown(index=False),
        "",
        "## Direct-NH entry quotient classes",
        "",
        skeletons.to_markdown(index=False),
        "",
        "## Ensemble registry",
        "",
        registry[["mechanism_id", "status", "scope", "comparison_role", "source"]].to_markdown(index=False),
        "",
        "## Next gate",
        "",
        "A candidate remains unavailable for numerical comparison until its complete kinetic input, rate provenance, species mapping and operating boundary are locally audited and harmonized. "
        "The Fe(110) model is source-acquired and structurally audited, but its runtime dependencies and operating boundary are not yet harmonized; it must not yet enter a flux comparison or a pathway-identifiability result figure.",
        "",
    ]
    (output / "P11_mechanism_family_audit.md").write_text("\n".join(lines), encoding="utf-8")


def draw_figure(output: Path, ledger: pd.DataFrame, families: pd.DataFrame, registry: pd.DataFrame) -> None:
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.2, "axes.labelsize": 8.3,
                         "axes.titlesize": 8.8, "xtick.labelsize": 6.5, "ytick.labelsize": 6.7,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(12.0, 8.4))
    grid = fig.add_gridspec(2, 2, left=0.13, right=0.985, top=0.94, bottom=0.07, hspace=0.48, wspace=0.34,
                            height_ratios=(1.0, 0.97))

    ax = fig.add_subplot(grid[0, 0])
    draw = families.sort_values("n_reactions", ascending=True)
    display_names = {
        "Other neutral gas-phase network": "Other neutral network",
        "Ion / electron network": "Ion/electron network",
        "Surface reaction network": "Surface network",
        "State relaxation / wall process": "State relaxation / wall",
        "N2 electronic activation": "N$_2$ electronic activation",
        "H2 electronic activation": "H$_2$ electronic activation",
        "N2 dissociation / active-N formation": "N$_2$ dissociation / N formation",
        "H2 dissociation / active-H formation": "H$_2$ dissociation / H formation",
        "NH/NH2 -> NH3 turnover": "NH/NH$_2$ → NH$_3$ turnover",
        "NH -> NH2 turnover": "NH → NH$_2$ turnover",
        "NH3 loss / ionic conversion": "NH$_3$ loss / ionic conversion",
        "NH2 loss / reverse conversion": "NH$_2$ loss / reverse conversion",
        "NH loss / reverse conversion": "NH loss / reverse conversion",
        "Direct NH entry: named electronic H2*": "NH entry: named H$_2$*",
        "Direct NH entry: Rydberg H2*": "NH entry: Rydberg H$_2$*",
        "Direct NH entry: vibrational H2": "NH entry: H$_2$(v)",
        "Direct NH entry: excited N + H2": "NH entry: excited N + H$_2$",
        "Direct NH entry: termolecular association": "NH entry: termolecular",
    }
    draw = draw.assign(display=draw["family"].map(display_names).fillna(draw["family"]))
    colors = ["#0072B2" if "Direct NH" in value else "#56B4E9" if value.endswith("activation") else "#999999" for value in draw["family"]]
    ax.barh(draw["display"], draw["n_reactions"], color=colors, edgecolor="white", linewidth=0.35)
    ax.set_xlabel("Active reaction records")
    ax.set_title("Corrected-Hong mechanism census by family")
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(-0.14, 1.08, "A", transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")

    ax = fig.add_subplot(grid[0, 1]); ax.axis("off")
    ax.set_title("Quotient-hypergraph rule for cross-model comparison", loc="left", pad=8)
    source_nodes = [(0.10, 0.75, "H$_2$(B$^3\\Sigma$)"), (0.10, 0.58, "H$_2$(B$^1\\Sigma$)"),
                    (0.10, 0.41, "H$_2$(C$^3\\Pi$)"), (0.10, 0.24, "H$_2$(A$^3\\Sigma$)")]
    for x, y, label in source_nodes:
        ax.add_patch(FancyBboxPatch((x, y - 0.05), 0.25, 0.10, boxstyle="round,pad=0.008", facecolor="#DDEBF7", edgecolor="#555555", lw=0.5))
        ax.text(x + 0.125, y, label, ha="center", va="center", fontsize=6.4)
        ax.add_patch(FancyArrowPatch((x + 0.25, y), (0.48, 0.50), arrowstyle="-|>", mutation_scale=8, lw=0.8, color="#777777"))
    ax.add_patch(FancyBboxPatch((0.48, 0.40), 0.25, 0.20, boxstyle="round,pad=0.012", facecolor="#0072B2", edgecolor="#555555", lw=0.6))
    ax.text(0.605, 0.53, "N + H$_2$(electronic)", ha="center", va="center", fontsize=7.1, color="white", fontweight="bold")
    ax.text(0.605, 0.46, "→ H + NH", ha="center", va="center", fontsize=7.1, color="white", fontweight="bold")
    ax.add_patch(FancyArrowPatch((0.73, 0.50), (0.84, 0.50), arrowstyle="-|>", mutation_scale=9, lw=1.2, color="#333333"))
    ax.add_patch(FancyBboxPatch((0.84, 0.42), 0.12, 0.16, boxstyle="round,pad=0.008", facecolor="#E8F3E5", edgecolor="#555555", lw=0.6))
    ax.text(0.90, 0.50, "NH-entry\nfamily", ha="center", va="center", fontsize=6.5, fontweight="bold")
    ax.text(0.03, 0.04, "State labels are retained in the ledger. The quotient class establishes structural correspondence only;\nfluxes and state-specific rates may break this topological equivalence.", fontsize=6.15, style="italic")
    ax.text(-0.12, 1.08, "B", transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")

    ax = fig.add_subplot(grid[1, 0]); ax.axis("off")
    ax.set_title("NH-entry taxonomy retained for mechanism discrimination", loc="left", pad=8)
    direct = ledger[ledger["layer"] == "NH entry"].groupby("family", as_index=False).agg(n=("reaction", "size"))
    rows = [[row.family.replace("Direct NH entry: ", ""), str(row.n), "candidate observable"] for row in direct.itertuples()]
    observable = {
        "named electronic H2*": "state-sensitive H2*/NH constraint",
        "Rydberg H2*": "lifetime / wall-loss constraint",
        "vibrational H2": "N2/H2 vibrational diagnostic",
        "excited N + H2": "N and excited-N diagnostic",
        "termolecular association": "pressure / density scaling",
    }
    for row in rows:
        row[2] = observable.get(row[0], "multi-observable model discrimination")
    table = ax.table(cellText=rows, colLabels=["Direct-NH family", "Active\nreactions", "Primary discriminating observable"],
                     cellLoc="left", colLoc="left", colWidths=[0.36, 0.17, 0.47], bbox=[0.0, 0.19, 1.0, 0.62])
    table.auto_set_font_size(False); table.set_fontsize(6.35)
    for (row, _col), cell in table.get_celld().items():
        cell.set_linewidth(0.35)
        if row == 0:
            cell.set_facecolor("#DDEBF7"); cell.set_text_props(weight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#F7F7F7")
    ax.text(0.0, 0.05, "The table is a test-design ontology; it contains neither fluxes nor experimental confirmation.", fontsize=6.15, style="italic")
    ax.text(-0.12, 1.08, "C", transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")

    ax = fig.add_subplot(grid[1, 1]); ax.axis("off")
    ax.set_title("Mechanism-ensemble readiness ledger", loc="left", pad=8)
    rows = []
    short_names = {
        "hong_2017_2018_corrected_local": "Hong corrected\n(local)",
        "elevated_pressure_acs_suschemeng_2025_candidate": "Elevated-pressure\nACS 2025 candidate",
        "pulsed_plasma_acs_suschemeng_2022_candidate": "Pulsed-plasma\nACS 2022 candidate",
        "fe110_dft_microkinetic_jacsau_2024": "Fe(110) DFT–microkinetic\nJACS Au 2024",
    }
    provenance_short = {
        "hong_2017_2018_corrected_local": "local archive",
        "elevated_pressure_acs_suschemeng_2025_candidate": "ACS DOI (2025)",
        "pulsed_plasma_acs_suschemeng_2022_candidate": "ACS DOI (2022)",
        "fe110_dft_microkinetic_jacsau_2024": "public code + JACS Au",
    }
    for row in registry.itertuples():
        rows.append([short_names[row.mechanism_id], row.status.replace("_", "\n"),
                     "yes" if row.status == "computable_local" else "no", provenance_short[row.mechanism_id]])
    table = ax.table(cellText=rows, colLabels=["Mechanism", "Status", "Simulation\nallowed?", "Provenance"],
                     cellLoc="left", colLoc="left", colWidths=[0.37, 0.19, 0.16, 0.28], bbox=[0.0, 0.27, 1.0, 0.56])
    table.auto_set_font_size(False); table.set_fontsize(5.55)
    for (row, _col), cell in table.get_celld().items():
        cell.set_linewidth(0.35)
        if row == 0:
            cell.set_facecolor("#DDEBF7"); cell.set_text_props(weight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#F7F7F7")
    ax.text(0.0, 0.12, "The Fe(110) model has reproduced author-boundary data locally, but has not been harmonized\nto the present CW/CSTR operating boundary. No cross-boundary external flux result is claimed.", fontsize=5.45, weight="bold")
    ax.text(0.0, 0.025, "Figure role: readiness candidate only; no cross-mechanism result is claimed before porting.", fontsize=5.65, style="italic")
    ax.text(-0.12, 1.08, "D", transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")

    fig.text(0.5, 0.016, "Mechanism-family audit: a structural comparison framework, not a claim of kinetic equivalence or external-model validation.", ha="center", fontsize=7)
    fig.savefig(output / "Figure_20_mechanism_ensemble_readiness.png", dpi=600)
    fig.savefig(output / "Figure_20_mechanism_ensemble_readiness.pdf")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Output already exists; refusing to overwrite: {args.output}")
    args.output.mkdir(parents=True)
    ledger = read_ledger(args.input)
    family_summary = ledger.groupby(["layer", "family"], as_index=False).agg(n_reactions=("reaction", "size")).sort_values(["layer", "family"])
    direct_skeletons = ledger[ledger["layer"] == "NH entry"].groupby(["family", "quotient_skeleton"], as_index=False).agg(
        n_members=("reaction", "size"), member_reactions=("reaction", lambda values: " | ".join(values))
    )
    input_hash = sha256(args.input)
    registry = mechanism_registry(input_hash)
    ledger.to_csv(args.output / "P11_reaction_family_ledger.csv", index=False)
    family_summary.to_csv(args.output / "P11_family_summary.csv", index=False)
    direct_skeletons.to_csv(args.output / "P11_direct_NH_quotient_classes.csv", index=False)
    registry.to_csv(args.output / "P11_mechanism_ensemble_registry.csv", index=False)
    surface_audit = representative_surface_flux_audit()
    surface_audit.to_csv(args.output / "P11_surface_network_boundary_audit.csv", index=False)
    write_markdown(args.output, ledger, family_summary, registry, surface_audit, input_hash)
    draw_figure(args.output, ledger, family_summary, registry)
    print(f"active_reactions={len(ledger)}")
    print(f"direct_NH_reactions={(ledger['layer'] == 'NH entry').sum()}")
    print(f"output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
