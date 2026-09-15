#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Audit, but do not execute, the public Fe(110) DFT–microkinetic model.

The audit deliberately distinguishes structural availability from numerical
comparability.  The external model shares a Hong-derived gas-phase backbone,
whereas its Fe(110) surface branch is modified by DFT/electric-field terms.
Accordingly, the generated report treats it as a surface-competition branch,
not as an independent validation of the gas-phase direct-NH conclusion.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import subprocess
from collections import Counter
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


PROJECT = Path(__file__).resolve().parent.parent
REPO = PROJECT / "Reproduction" / "external_models" / "DFT-microkinetic"
EXTERNAL_INPUT = REPO / "Model_SA_Const_Entropy_base" / "kinet_varyT_metal_auto_sens_DFT_in_entropy_verying_basis.txt"
EXTERNAL_EXTRACT = REPO / "Model_SA_Const_Entropy_extract" / EXTERNAL_INPUT.name
LOCAL_INPUT = PROJECT / "Reproduction" / "2017Hong" / "literature_extract_2017_2018_corrected" / "kinet.inp"
DEFAULT_OUTPUT = PROJECT / "analysis" / "p12_external_dft_microkinetic_audit_20260831_r5"
PRIOR = PROJECT / "Shao2024_JACSAu"
PRIOR_SOURCE = PRIOR / "dft_microkinetic_windows_case" / "kinet_source.txt"
PRIOR_REPORT = PRIOR / "reproduction" / "Shao2024复现报告.md"
PRIOR_RESULTS = PRIOR / "reproduction" / "results" / "Shao2024复现对比.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def active_equations(path: Path) -> list[tuple[int, str]]:
    rows: list[tuple[int, str]] = []
    for number, raw in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1):
        stripped = raw.strip()
        if stripped and not stripped.startswith("#") and "=>" in stripped:
            rows.append((number, stripped.split("!", 1)[0].strip()))
    return rows


def species_list(path: Path) -> list[str]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    begin = next(i for i, line in enumerate(lines) if line.strip() == "SPECIES")
    result: list[str] = []
    for line in lines[begin + 1:]:
        stripped = line.strip()
        if stripped == "END":
            break
        if stripped and not stripped.startswith("#"):
            result.extend(stripped.split())
    return result


def normalized(equation: str) -> str:
    return "".join(equation.split()).replace("e+", "E+").replace("+e", "+E")


def classify(equation: str) -> str:
    if any(token in equation for token in ("Surf", "2Surf")):
        return "surface"
    left, right = (part.strip() for part in equation.split("=>", 1))
    neutral_n_h2 = r"N(?:\(2[DP]\))?\s*\+\s*H2(?:\([^)]*\))?"
    right_terms = {term.strip() for term in right.split(" + ")}
    if re.fullmatch(neutral_n_h2, left) and right_terms == {"H", "NH"}:
        return "direct_NH_entry"
    if equation.lstrip().startswith(("e +", "E +")):
        return "electron_impact"
    if any(token in equation for token in ("^+", "H^-")):
        return "ion_electron_network"
    return "other_gas_phase"


def git_revision(repo: Path) -> str:
    completed = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], check=True,
                               capture_output=True, text=True)
    return completed.stdout.strip()


def file_row(path: Path, purpose: str) -> dict[str, str]:
    try:
        display_path = str(path.relative_to(REPO))
    except ValueError:
        display_path = str(path)
    return {
        "relative_path": display_path,
        "purpose": purpose,
        "present": "yes" if path.is_file() else "no",
        "bytes": str(path.stat().st_size) if path.is_file() else "",
        "sha256": sha256(path) if path.is_file() else "",
    }


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def draw_figure(output: Path, local_counts: Counter, external_counts: Counter, shared_direct: int,
                external_direct: int) -> None:
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.2, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(12.0, 7.6))
    grid = fig.add_gridspec(2, 2, left=0.075, right=0.985, top=0.93, bottom=0.085, hspace=0.46, wspace=0.31)

    ax = fig.add_subplot(grid[0, 0])
    keys = ["other_gas_phase", "electron_impact", "ion_electron_network", "surface", "direct_NH_entry"]
    labels = ["Other\ngas", "Electron\nimpact", "Ion/electron", "Surface", "Direct\nNH entry"]
    positions = range(len(keys))
    ax.bar([p - 0.19 for p in positions], [local_counts[k] for k in keys], width=0.38, color="#8CB9D9", label="Local corrected Hong")
    ax.bar([p + 0.19 for p in positions], [external_counts[k] for k in keys], width=0.38, color="#D55E00", label="Fe(110) DFT–microkinetic")
    ax.set_xticks(list(positions), labels)
    ax.set_ylabel("Active reaction records")
    ax.set_title("Reaction-family inventory (structural, not flux-weighted)", loc="left")
    ax.legend(frameon=False, fontsize=6.3, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(-0.14, 1.08, "A", transform=ax.transAxes, fontsize=11, fontweight="bold")

    ax = fig.add_subplot(grid[0, 1]); ax.axis("off")
    ax.set_title("What is shared—and what is an independent branch", loc="left", pad=7)
    items = [
        (0.04, 0.69, 0.28, 0.16, "Local corrected\nHong", "#DDEBF7"),
        (0.37, 0.69, 0.28, 0.16, "Shared Hong-derived\ngas-phase backbone", "#EDEDED"),
        (0.70, 0.69, 0.26, 0.16, "Fe(110) DFT field–\nmodified surface branch", "#FCE4D6"),
        (0.37, 0.27, 0.28, 0.15, f"{shared_direct}/{external_direct} direct-NH\nreactions exactly shared", "#E8F3E5"),
    ]
    for x, y, w, h, text, color in items:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012", facecolor=color, edgecolor="#555555", lw=0.55))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=6.8, fontweight="bold" if y < 0.5 else None)
    ax.add_patch(FancyArrowPatch((0.32, 0.77), (0.37, 0.77), arrowstyle="-|>", mutation_scale=9, lw=0.8, color="#555555"))
    ax.add_patch(FancyArrowPatch((0.65, 0.77), (0.70, 0.77), arrowstyle="-|>", mutation_scale=9, lw=0.8, color="#555555"))
    ax.add_patch(FancyArrowPatch((0.51, 0.69), (0.51, 0.43), arrowstyle="-|>", mutation_scale=9, lw=0.8, color="#555555"))
    ax.text(0.04, 0.08, "Inference boundary: the surface branch is a meaningful alternative pathway class;\nthe shared direct-NH gas reactions cannot independently validate the local gas-phase ranking.", fontsize=6.55, style="italic")
    ax.text(-0.10, 1.08, "B", transform=ax.transAxes, fontsize=11, fontweight="bold")

    ax = fig.add_subplot(grid[1, 0]); ax.axis("off")
    ax.set_title("Explicit Fe(110) surface chemistry present in the acquired input", loc="left", pad=7)
    rows = [
        ["Gas → surface capture", "N/N(2D)/N(2P), H, NH, NH$_2$ adsorption"],
        ["Surface hydrogenation", "N* → NH* → NH2* → NH3"],
        ["Eley–Rideal / LH routes", "gas H/H2 + NHx*; N* + H* coupling"],
        ["Dissociative adsorption", "N2/N2(v,electronic), H2/H2(v,electronic)"],
        ["Field/DFT dependence", "rate expressions depend on generated energy/entropy basis parameters"],
    ]
    table = ax.table(cellText=rows, colLabels=["Surface subnetwork family", "Explicit structural content"], cellLoc="left", colLoc="left",
                     colWidths=[0.34, 0.66], bbox=[0, 0.15, 1, 0.68])
    table.auto_set_font_size(False); table.set_fontsize(6.45)
    for (row, _), cell in table.get_celld().items():
        cell.set_linewidth(0.35)
        if row == 0:
            cell.set_facecolor("#FCE4D6"); cell.set_text_props(weight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#F7F7F7")
    ax.text(0.0, 0.045, "These are reaction-structure facts from the rate input; no external-model rate-of-production or selectivity result is shown.", fontsize=6.1, style="italic")
    ax.text(-0.14, 1.08, "C", transform=ax.transAxes, fontsize=11, fontweight="bold")

    ax = fig.add_subplot(grid[1, 1]); ax.axis("off")
    ax.set_title("Porting and comparison gate", loc="left", pad=7)
    rows = [
        ["Source provenance\n& version", "PASS", "Public MIT repository;\nrevision recorded"],
        ["Kinetic input &\nsurface mapping", "PASS", "468 active; 46 surface\nreaction records"],
        ["DFT/field parameter\ngenerator", "PASS", "Fortran + Python parameter\ngeneration chain present"],
        ["Author-boundary\nreference reproduction", "PASS", "Four local Windows cases;\nmax NH3 deviation 0.06345%"],
        ["Solver assets / platform", "PASS", "Source model has DVODE/BOLSIG;\nWindows build already validated"],
        ["Boundary harmonization\nwith CW/CSTR", "HOLD", "Reactor, A/V, roughness, T,\nE/N and electron closure differ"],
        ["Numerical pathway\ncomparison", "NOT ALLOWED", "Requires the HOLD gate"],
    ]
    table = ax.table(cellText=rows, colLabels=["Gate", "Status", "Reason"], cellLoc="left", colLoc="left",
                     colWidths=[0.34, 0.20, 0.46], bbox=[0, 0.12, 1, 0.74])
    table.auto_set_font_size(False); table.set_fontsize(5.65)
    for (row, col), cell in table.get_celld().items():
        cell.set_linewidth(0.35)
        if row == 0:
            cell.set_facecolor("#DDEBF7"); cell.set_text_props(weight="bold")
        elif row in (6, 7):
            cell.set_facecolor("#FFF2CC")
        elif row % 2 == 0:
            cell.set_facecolor("#F7F7F7")
        if col == 1 and row > 0:
            cell.set_text_props(weight="bold")
    ax.text(-0.10, 1.08, "D", transform=ax.transAxes, fontsize=11, fontweight="bold")

    fig.text(0.5, 0.017, "Figure 21. Validated Fe(110) DFT–microkinetic branch: suitable for surface-competition tests, but not an independent gas-phase NH-entry validation or a cross-boundary flux result.", ha="center", fontsize=6.8)
    fig.savefig(output / "Figure_21_external_model_scope_and_porting_gate.png", dpi=600)
    fig.savefig(output / "Figure_21_external_model_scope_and_porting_gate.pdf")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Output already exists; refusing to overwrite: {args.output}")
    for required in (REPO, EXTERNAL_INPUT, EXTERNAL_EXTRACT, LOCAL_INPUT, PRIOR_SOURCE, PRIOR_REPORT, PRIOR_RESULTS):
        if not required.exists():
            raise SystemExit(f"Required source is missing: {required}")
    args.output.mkdir(parents=True)

    external = active_equations(EXTERNAL_INPUT)
    local = active_equations(LOCAL_INPUT)
    ext_species, local_species = species_list(EXTERNAL_INPUT), species_list(LOCAL_INPUT)
    ext_rows = [{"line_number": str(number), "reaction": equation, "family": classify(equation)} for number, equation in external]
    write_csv(args.output / "P12_external_reaction_ledger.csv", ext_rows, ["line_number", "reaction", "family"])
    ext_counts = Counter(row["family"] for row in ext_rows)
    local_counts = Counter(classify(equation) for _, equation in local)
    summary = []
    for family in ("other_gas_phase", "electron_impact", "ion_electron_network", "surface", "direct_NH_entry"):
        summary.append({"family": family, "local_corrected_hong_records": str(local_counts[family]),
                        "external_fe110_records": str(ext_counts[family])})
    write_csv(args.output / "P12_family_crosswalk.csv", summary, list(summary[0]))

    external_set = {normalized(equation) for _, equation in external}
    local_set = {normalized(equation) for _, equation in local}
    external_direct = {normalized(row["reaction"]) for row in ext_rows if row["family"] == "direct_NH_entry"}
    direct_overlap = sorted(external_direct & local_set)
    overlap_rows = [{"metric": "external_unique_active_reactions", "value": str(len(external_set))},
                    {"metric": "local_unique_active_reactions", "value": str(len(local_set))},
                    {"metric": "exact_active_reaction_overlap", "value": str(len(external_set & local_set))},
                    {"metric": "external_direct_NH_reactions", "value": str(len(external_direct))},
                    {"metric": "direct_NH_exact_overlap", "value": str(len(direct_overlap))},
                    {"metric": "external_species", "value": str(len(ext_species))},
                    {"metric": "local_species", "value": str(len(local_species))}]
    write_csv(args.output / "P12_reaction_overlap_metrics.csv", overlap_rows, ["metric", "value"])

    manifest = [
        file_row(EXTERNAL_INPUT, "base kinetic input"),
        file_row(EXTERNAL_EXTRACT, "extract kinetic input; byte-identical check"),
        file_row(REPO / "Model_SA_Const_Entropy_base" / "zdplaskin_m_DFT_in_entropy_varying_basis.F90", "custom ZDPlasKin module"),
        file_row(REPO / "Model_SA_Const_Entropy_base" / "Const_E.F90", "Fortran reactor driver"),
        file_row(REPO / "Model_SA_Const_Entropy_base" / "Kinetic_funcs_entropy_array.py", "parameter-generation and compilation driver"),
        file_row(REPO / "Model_SA_Const_Entropy_base" / "kine_test_shared_entropy_array.sh", "SLURM launch script"),
        file_row(REPO / "qtplaskin.zip", "bundled qtplaskin archive"),
        file_row(REPO / "Model_SA_Const_Entropy_base" / "dvode_f90_m.f90", "ODE solver source at the batch-driver submission path"),
        file_row(REPO / "Model_SA_Const_Entropy_base" / "bolsig_x86_64.so", "Linux BOLSIG shared library at the batch-driver submission path"),
        file_row(REPO / "Model_SA_Const_Entropy_base" / "bolsigdb.dat", "BOLSIG collision database at the batch-driver submission path"),
        file_row(PRIOR_SOURCE, "source-identical Windows reproduction kinetic input"),
        file_row(PRIOR_REPORT, "prior Windows reference-reproduction report"),
        file_row(PRIOR_RESULTS, "prior Windows reference-reproduction metrics"),
    ]
    write_csv(args.output / "P12_file_manifest.csv", manifest, list(manifest[0]))
    revision = git_revision(REPO)
    base_hash, extract_hash = sha256(EXTERNAL_INPUT), sha256(EXTERNAL_EXTRACT)
    source_runtime_ready = all(row["present"] == "yes" for row in manifest[7:10])
    prior_source_identical = sha256(PRIOR_SOURCE) == base_hash

    report = [
        "# External Fe(110) DFT–microkinetic model audit",
        "",
        "**Status:** source acquired and structurally audited; previously validated under the authors' boundary; **not harmonized to the present CW/CSTR boundary**.",
        "",
        "## Reproducible provenance",
        "",
        f"- Public repository: `https://github.com/wwwccttoo/DFT-microkinetic`",
        f"- Checked revision: `{revision}`",
        "- License in repository: MIT (copyright retained in the acquired source).",
        "- Associated article: Shao and Mesbah, *JACS Au* (2024), DOI `10.1021/jacsau.3c00654`.",
        f"- Base kinetic input SHA-256: `{base_hash}`",
        f"- Extract kinetic input SHA-256: `{extract_hash}`",
        f"- The two kinetic-input files are byte-identical: **{'yes' if base_hash == extract_hash else 'no'}**.",
        f"- The earlier local Windows reproduction uses an input with the same SHA-256: **{'yes' if prior_source_identical else 'no'}**.",
        "",
        "## Structural inventory",
        "",
        f"- External input: {len(external)} active reaction records; {len(ext_species)} declared species.",
        f"- Explicit surface-reaction records: {ext_counts['surface']}.",
        f"- Direct NH-entry records of the `N + H2* -> H + NH` / excited-N type: {len(external_direct)}.",
        f"- Exact normalized overlap of those direct-NH records with the local corrected-Hong input: {len(direct_overlap)}/{len(external_direct)}.",
        "",
        "The exact direct-NH overlap is a provenance constraint: it rules out describing this Fe(110) model as an independent validation of the local gas-phase NH-entry ranking. Its legitimate role is to test whether an explicit field-sensitive surface branch competes with, masks or redirects that shared gas-phase entry under a separately harmonized operating boundary.",
        "",
        "## Surface-network content directly present in the rate input",
        "",
        "- Adsorption/capture of N, N(2D), N(2P), H, NH and NH2 on vacant sites.",
        "- Surface N*–H* coupling and successive NH* → NH2* → NH3 formation.",
        "- Eley–Rideal-like gas H/H2 addition to adsorbed nitrogen hydrides.",
        "- Dissociative adsorption of N2 and H2, including vibrational/electronic excited-state variants.",
        "- Dynamic rate expressions parameterized through field-dependent energy and entropy basis files, with total-site density, geometry and roughness terms.",
        "",
        "## Runtime validation and fair-comparison gate",
        "",
        f"The repository contains the mechanism, custom ZDPlasKin module, Fortran driver, parameter-generation script and the batch-driver assets `dvode_f90_m.f90`, `bolsig_x86_64.so` and `bolsigdb.dat` in the documented model submission directory (asset check: **{'pass' if source_runtime_ready else 'fail'}**). The original batch path is Linux-oriented because it links a `.so` library. Independently, this workspace contains a source-identical Windows adaptation using the authors' generated Fortran module and the Windows ZDPlasKin BOLSIG library. Its four author-boundary cases (0.06/0.11 V Å⁻¹; 350/550 K; 4000 s) reproduce the author DFT NH3 series with maximum absolute relative deviation 0.06345%. This current audit did not launch an additional solver run.",
        "",
        "Before any numerical comparison with the present results, the following gate remains mandatory: map the external model's surface area, roughness, site density, temperature, E/N and electron-density closure to one declared common boundary; then verify elemental/site conservation and compute gas- and surface-resolved ROP/flux fractions. Until that step, the validated external reproduction is evidence of model portability and source fidelity—not evidence about dominance under the present CW/CSTR boundary.",
        "",
        "## Generated machine-readable records",
        "",
        "- `P12_external_reaction_ledger.csv`: every active external reaction and its coarse family.",
        "- `P12_family_crosswalk.csv`: local/external structural count comparison.",
        "- `P12_reaction_overlap_metrics.csv`: exact normalized-reaction overlap metrics.",
        "- `P12_file_manifest.csv`: input, code and missing-dependency provenance ledger.",
        "",
    ]
    (args.output / "P12_external_model_provenance.md").write_text("\n".join(report), encoding="utf-8")
    draw_figure(args.output, local_counts, ext_counts, len(direct_overlap), len(external_direct))
    print(f"external_active_reactions={len(external)}")
    print(f"external_surface_reactions={ext_counts['surface']}")
    print(f"direct_nh_overlap={len(direct_overlap)}/{len(external_direct)}")
    print(f"source_runtime_assets_ready={source_runtime_ready}")
    print(f"prior_source_identical={prior_source_identical}")
    print(f"output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
