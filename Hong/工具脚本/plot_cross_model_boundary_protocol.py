#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create a physically guarded protocol for Hong–Shao model comparison.

This is a boundary/identifiability figure, not a simulation.  It makes the
otherwise hidden field-scale mismatch explicit: gas reduced field E/N and
surface local field F_s are linked only through a device-specific enhancement
factor beta; they must never be numerically equated.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = PROJECT / "analysis" / "p13_cross_model_boundary_protocol_20260901_r5"
K_B = 1.380649e-23  # J K-1
P_ATM = 101325.0  # Pa
HONG_T = 300.0  # K
EN_TD = (20.0, 140.0, 240.0)
SURFACE_FIELDS = (0.06, 0.11)  # V Angstrom-1, Shao verified reference settings


def gas_density(temperature_k: float) -> float:
    return P_ATM / (K_B * temperature_k)


def e_bulk_v_per_a(en_td: float, temperature_k: float) -> float:
    return en_td * 1e-21 * gas_density(temperature_k) / 1e10


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def draw(output: Path, ledger: list[dict[str, float]]) -> None:
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.2, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(12.0, 8.15))
    grid = fig.add_gridspec(2, 2, left=0.075, right=0.985, top=0.935, bottom=0.075, hspace=0.49, wspace=0.31)

    ax = fig.add_subplot(grid[0, 0])
    ens = [row["en_td"] for row in ledger]
    bulk = [row["bulk_v_per_a"] for row in ledger]
    ax.plot(ens, bulk, marker="o", ms=4.2, lw=1.5, color="#0072B2", label="Bulk gas field from E/N at 1 atm, 300 K")
    for surface, color in zip(SURFACE_FIELDS, ("#D55E00", "#CC79A7")):
        ax.axhline(surface, color=color, lw=1.15, ls="--", label=f"Shao reference F$_s$ = {surface:.2f} V Å$^{{-1}}$")
    for index, row in enumerate(ledger):
        if index == 0:
            ax.text(34, 1.25e-4, f"β={row['beta_for_0p06']:.0f}–{row['beta_for_0p11']:.0f}",
                    ha="center", va="center", fontsize=6.2, color="#333333")
            continue
        offset = (9, 10) if index == 0 else (0, -17)
        alignment = "left" if index == 0 else "center"
        ax.annotate(f"β={row['beta_for_0p06']:.0f}–{row['beta_for_0p11']:.0f}", (row["en_td"], row["bulk_v_per_a"]),
                    xytext=offset, textcoords="offset points", ha=alignment, fontsize=6.2, color="#333333")
    ax.set_yscale("log")
    ax.set_xlabel("Gas reduced field, E/N (Td)")
    ax.set_ylabel("Field magnitude (V Å$^{-1}$)")
    ax.set_title("Bulk E/N and surface F$_s$ occupy different physical scales", loc="left")
    ax.legend(frameon=False, fontsize=5.9, loc="upper left")
    ax.grid(axis="y", alpha=0.22, lw=0.45)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.02, 0.04, r"Required relation: $F_s=\beta E_{bulk}$; β is a device- and material-dependent\nquantity, not a fitted universal constant and not assumed here.", transform=ax.transAxes, fontsize=6.25, style="italic")
    ax.text(-0.14, 1.08, "A", transform=ax.transAxes, fontsize=11, fontweight="bold")

    ax = fig.add_subplot(grid[0, 1]); ax.axis("off")
    ax.set_title("Five gates before a gas–surface flux comparison is allowed", loc="left", pad=7)
    steps = [
        ("G0", "Frozen sources", "Hong and Shao inputs, hashes, and author-boundary reproduction"),
        ("G1", "Field transfer", r"Obtain $F_s=\beta(E/N)N_g$ from electrostatics or calibrated measurement"),
        ("G2", "Shared reactor", r"Set p, T, x$_{N2}$, n$_e$, A/V, roughness, site density, τ and inlet"),
        ("G3", "Numerical acceptance", "Implement CSTR + site balance; pass common P0 and conservation checks"),
        ("G4", "Identifiability", "Report gas/surface ROP fractions and out-of-sample observables"),
    ]
    y_values = [0.85, 0.68, 0.51, 0.34, 0.17]
    for (gate, title, text), y in zip(steps, y_values):
        color = "#DDEBF7" if gate in ("G0", "G1") else "#FFF2CC" if gate in ("G2", "G3") else "#E8F3E5"
        ax.add_patch(FancyBboxPatch((0.05, y - 0.062), 0.12, 0.105, boxstyle="round,pad=0.007", facecolor=color, edgecolor="#555555", lw=0.5))
        ax.text(0.11, y - 0.01, gate, ha="center", va="center", fontsize=7.3, fontweight="bold")
        ax.add_patch(FancyBboxPatch((0.20, y - 0.062), 0.74, 0.105, boxstyle="round,pad=0.007", facecolor="#FAFAFA", edgecolor="#777777", lw=0.45))
        ax.text(0.23, y + 0.015, title, ha="left", va="center", fontsize=6.7, fontweight="bold")
        ax.text(0.23, y - 0.027, text, ha="left", va="center", fontsize=5.65)
    for y0, y1 in zip(y_values[:-1], y_values[1:]):
        ax.add_patch(FancyArrowPatch((0.11, y0 - 0.065), (0.11, y1 + 0.048), arrowstyle="-|>", mutation_scale=8, lw=0.65, color="#555555"))
    ax.text(-0.10, 1.08, "B", transform=ax.transAxes, fontsize=11, fontweight="bold")

    ax = fig.add_subplot(grid[1, 0]); ax.axis("off")
    ax.set_title("Boundary ledger: variables that can—and cannot—be matched", loc="left", pad=7)
    rows = [
        ["Gas state", "p, T, xN$_2$, n$_e$", "Directly prescribe / measure"],
        ["Electrical forcing", "E/N vs F$_s$", r"Only through $F_s=\beta E_{bulk}$"],
        ["Transport", "τ, inlet/outlet, volume", "Add identical CSTR source/sink terms"],
        ["Surface boundary", "A/V, roughness, site density", "Measure or specify per catalyst"],
        ["Mechanism", "Hong gas + Fe surface branch", "Keep source-specific rate provenance"],
        ["Observables", "NH, NH$_x$, NH$_3$, electric data", "Use independent signals, not NH$_3$ alone"],
    ]
    table = ax.table(cellText=rows, colLabels=["Layer", "Variables", "Valid comparison rule"], cellLoc="left", colLoc="left",
                     colWidths=[0.22, 0.31, 0.47], bbox=[0.0, 0.13, 1.0, 0.72])
    table.auto_set_font_size(False); table.set_fontsize(6.0)
    for (row, _), cell in table.get_celld().items():
        cell.set_linewidth(0.35)
        if row == 0:
            cell.set_facecolor("#DDEBF7"); cell.set_text_props(weight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#F7F7F7")
    ax.text(0.0, 0.045, "A common gas composition alone is insufficient: the surface rate constants embed geometry, site density and local field.", fontsize=6.1, style="italic")
    ax.text(-0.14, 1.08, "C", transform=ax.transAxes, fontsize=11, fontweight="bold")

    ax = fig.add_subplot(grid[1, 1]); ax.axis("off")
    ax.set_title("Pre-registered outputs for pathway identifiability", loc="left", pad=7)
    equations = [
        (r"$\phi_{NH}^{gas}=J_{NH,gas}/(J_{NH,gas}+J_{NH,surf})$", "Gas contribution to first\nN–H formation"),
        (r"$\phi_{NH}^{surf}=J_{NH,surf}/(J_{NH,gas}+J_{NH,surf})$", "Surface contribution to first\nN–H formation"),
        (r"$\mathcal{M}_{gas}=J_{H_2^*\rightarrow NH}/\max(J_{secondary,gas})$", "Gas direct-NH\ndominance margin"),
        (r"$\mathcal{D}_{surf}=\tau\,J_{NH,surf}/n_{N,ref}$", "Surface N–H entry\nDamköhler-like index"),
        (r"$\mathcal{I}=\{NH, NH_x, NH_3, V\! -\! I\! -\! Q, ^{15}NH_3\}$", "Minimum multi-observable\ndiscrimination set"),
    ]
    y = 0.82
    for eq, description in equations:
        ax.add_patch(FancyBboxPatch((0.04, y - 0.065), 0.57, 0.10, boxstyle="round,pad=0.008", facecolor="#F4F7FB", edgecolor="#777777", lw=0.4))
        ax.text(0.325, y - 0.014, eq, ha="center", va="center", fontsize=7.0)
        ax.text(0.65, y - 0.014, description, ha="left", va="center", fontsize=5.75)
        y -= 0.15
    ax.text(0.04, 0.055, "Only after G0–G3: report these quantities with uncertainty and compare a prediction\nthat was not used to calibrate β or surface parameters.", fontsize=6.15, weight="bold")
    ax.text(-0.10, 1.08, "D", transform=ax.transAxes, fontsize=11, fontweight="bold")

    fig.text(0.5, 0.017, "Figure 22. Boundary-harmonization and pathway-identifiability protocol for gas-phase Hong and Fe(110) DFT–microkinetic models. It is a comparison guardrail, not a surface-flux result.", ha="center", fontsize=6.75)
    fig.savefig(output / "Figure_22_boundary_harmonization_protocol.png", dpi=600)
    fig.savefig(output / "Figure_22_boundary_harmonization_protocol.pdf")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Output already exists; refusing to overwrite: {args.output}")
    args.output.mkdir(parents=True)
    ledger = []
    for en in EN_TD:
        bulk = e_bulk_v_per_a(en, HONG_T)
        ledger.append({"en_td": en, "bulk_v_per_a": bulk, "bulk_v_per_m": bulk * 1e10,
                       "beta_for_0p06": 0.06 / bulk, "beta_for_0p11": 0.11 / bulk})
    write_csv(args.output / "P13_field_scale_ledger.csv",
              ["pressure_pa", "temperature_k", "gas_number_density_m-3", "EN_Td", "E_bulk_V_m", "E_bulk_V_A", "beta_for_Fs_0.06_V_A", "beta_for_Fs_0.11_V_A"],
              [[P_ATM, HONG_T, gas_density(HONG_T), row["en_td"], row["bulk_v_per_m"], row["bulk_v_per_a"], row["beta_for_0p06"], row["beta_for_0p11"]] for row in ledger])
    schema = [
        ["pressure", "Pa", "101325", "common prescribed or measured gas pressure", "required"],
        ["gas_temperature", "K", "Hong: 300; Shao reference: 350/550", "must use one common temperature in a comparison", "required"],
        ["x_N2", "mole fraction", "Hong: 0.1–0.9; Shao reference: 1/3", "common feed state", "required"],
        ["electron_density", "cm-3", "Hong: 1.17e8; Shao reference: 8.27e7", "directly measure, calibrate or declare scenario", "required"],
        ["bulk_field", "Td and V m-1", "Hong: 20–240 Td", "E_bulk=(E/N)p/(kBT)", "required"],
        ["surface_field", "V A-1", "Shao reference: 0.06/0.11", "F_s=beta E_bulk; beta needs electrostatic or experimental constraint", "required"],
        ["A_over_V", "model inverse length", "Shao native A*roughness/V = 145.56", "declare whether this or a measured reactor value is shared", "required"],
        ["roughness", "dimensionless", "Shao native: 2.1", "retain source value or define explicit catalyst scenario", "required"],
        ["site_density", "model site-density inputs", "Shao: Tot_sur=1e15; driver Surf=2.87e17", "preserve site balance and do not conflate the two normalizations", "required"],
        ["residence_time", "s", "Hong: 0.010", "add CSTR inlet/outlet to gas and report tau", "required"],
    ]
    write_csv(args.output / "P13_boundary_variable_schema.csv", ["variable", "unit", "current_state", "harmonization_rule", "gate"], schema)
    report = [
        "# Boundary-harmonization and pathway-identifiability protocol",
        "",
        "**Status:** design and dimensional audit. No hybrid CSTR calculation or surface flux result is reported.",
        "",
        "## Why direct field matching is invalid",
        "",
        "The Hong calculation is parameterized by the gas reduced electric field E/N (Td). At pressure p and gas temperature T, the corresponding bulk field is E_bulk = (E/N)p/(k_B T). The Shao Fe(110) model is parameterized by a local surface field F_s (V Å-1) that modifies surface energetics. A physically meaningful relation is F_s = beta E_bulk, where beta contains sheath, geometry, dielectric, packing, material and local-enhancement physics. beta is neither supplied by the present Hong CSTR calculation nor assumed to be one.",
        "",
        "At 1 atm and 300 K, the 20, 140 and 240 Td Hong states correspond to bulk fields of " + ", ".join(f"{row['bulk_v_per_m'] / 1e6:.3f} MV m-1" for row in ledger) + ". To reach the two Shao reference local fields (0.06 and 0.11 V Å-1), required beta spans " + f"{min(row['beta_for_0p06'] for row in ledger):.0f}–{max(row['beta_for_0p11'] for row in ledger):.0f}" + " across those anchors. These are required enhancement factors, not inferred properties of any real device.",
        "",
        "## Comparison protocol",
        "",
        "1. Freeze the source-validated Hong and Shao inputs and retain their separate rate provenance.",
        "2. Obtain beta(E/N, geometry, material) independently through an electrostatic model constrained by device geometry or through calibrated field diagnostics. Do not fit beta only to NH3 and then claim an independent pathway test.",
        "3. Specify a shared reactor state: p, T, x_N2, electron-density closure, A/V, roughness, site density, inlet composition and residence time.",
        "4. Add the identical CSTR inflow/outflow operator to the gas equations while preserving a surface-site balance; require P0, elemental conservation and site conservation before comparing terminal windows.",
        "5. Report gas/surface first-N–H flux fractions, direct-NH margin, NH3 source/sink fluxes and at least one out-of-sample observable such as NH/NHx diagnostics or isotope-resolved NH3.",
        "",
        "## Native surface boundary recovered from the source",
        "",
        "The Shao driver fixes gas pressure to 1 bar and sets `surface_site_density = 2.87e17` for the `Surf` state. The generated surface module sets `THE_V = 48.6`, `THE_AREA = 3370`, `ROUGHNESS = 2.1`, and `TOT_SUR = 1e15`. Thus its model-native effective area-to-volume factor is `THE_AREA*ROUGHNESS/THE_V = 145.556` in the source model's length convention. These values make the source boundary auditable; they are not evidence that the present Hong CSTR has the same geometry, catalyst loading or site normalization.",
        "",
        "## Current gate state",
        "",
        "- G0 source provenance and author-boundary external reproduction: passed.",
        "- G1 field-transfer relation beta: not supplied; numerical mapping is therefore not authorized.",
        "- G2 shared transport/surface geometry: the Shao-native constants are supplied, but no justified declaration maps them to the present Hong CSTR geometry, inlet/outlet operator or catalyst inventory.",
        "- G3 hybrid CSTR numerical acceptance: not attempted.",
        "- G4 cross-model pathway comparison: not allowed until G1–G3 pass.",
        "",
        "The protocol turns a limitation into a falsifiable research object: a surface model may overturn or preserve the gas NH-entry ranking only after the model supplies the physical bridge between bulk discharge forcing and catalyst-local energetics.",
    ]
    (args.output / "P13_boundary_harmonization_protocol.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    draw(args.output, ledger)
    print(f"gas_number_density_300K={gas_density(HONG_T):.8e}")
    print(f"required_beta_range={min(row['beta_for_0p06'] for row in ledger):.3f}-{max(row['beta_for_0p11'] for row in ledger):.3f}")
    print(f"output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
