"""Create and evaluate a no-claim preflight contract for Hong–Fe(110) coupling.

The checker deliberately refuses to label a gas–surface pathway comparison as
ready until five gates are satisfied. It does not calculate a surface flux,
infer beta, alter either source mechanism, or fill missing physical parameters.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import date
from pathlib import Path
from typing import Any

from _paths import LEGACY_HONG_ROOT

ROOT = LEGACY_HONG_ROOT
DEFAULT_OUT = ROOT / "analysis" / "p18_hybrid_comparison_preflight_20260901_r1"
HONG_INPUT = ROOT / "Reproduction" / "2017Hong" / "literature_extract_2017_2018_corrected" / "kinet.inp"
SHAO_INPUT = ROOT / "Shao2024_JACSAu" / "dft_microkinetic_windows_case" / "kinet_source.txt"
EXPECTED_HONG_SHA = "d9b0da025e1eeeeF42fa58963dadcc923366d89ad0621a8ab4516da966794f38".lower()
EXPECTED_SHAO_SHA = "7ce49610345024ce7505433f27bf2dea3ab07d139b6130470c0464fbd0893cc5"
K_B = 1.380649e-23


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def filled(value: Any) -> bool:
    return value is not None and value != "" and value is not False


def template() -> dict[str, Any]:
    """Return a deliberately incomplete comparison request without hidden defaults."""
    return {
        "contract_version": "1.0",
        "purpose": "Preflight only; no gas-surface ROP or performance claim is permitted until all gates pass.",
        "source_locks": {
            "hong_kinetic_input": {
                "path": str(HONG_INPUT.relative_to(ROOT)),
                "expected_sha256": EXPECTED_HONG_SHA,
            },
            "shao_fe110_kinetic_input": {
                "path": str(SHAO_INPUT.relative_to(ROOT)),
                "expected_sha256": EXPECTED_SHAO_SHA,
            },
        },
        "field_transfer": {
            "pressure_pa": None,
            "gas_temperature_k": None,
            "reduced_field_td": None,
            "surface_field_v_per_a": None,
            "beta": None,
            "beta_evidence": None,
            "beta_evidence_options": [
                "electrostatic model constrained by reactor geometry",
                "independently calibrated local-field diagnostic",
            ],
        },
        "shared_reactor_boundary": {
            "n2_mole_fraction": None,
            "h2_mole_fraction": None,
            "electron_density_cm3": None,
            "electron_density_closure": None,
            "residence_time_s": None,
            "inlet_outlet_operator_identical": False,
            "gas_volume_cm3": None,
            "catalyst_area_to_volume_cm_minus1": None,
            "roughness": None,
            "surface_site_normalization": None,
            "catalyst_inventory_basis": None,
        },
        "hybrid_numerical_acceptance": {
            "hybrid_input_frozen": False,
            "terminal_window_p0_passed": False,
            "nitrogen_conservation_passed": False,
            "hydrogen_conservation_passed": False,
            "surface_site_conservation_passed": False,
            "solver_log_or_artifact": None,
        },
        "pathway_identifiability": {
            "gas_and_surface_first_nh_fluxes_reported": False,
            "nh3_source_sink_fluxes_reported": False,
            "out_of_sample_observable": None,
            "out_of_sample_observable_options": [
                "NH or NHx diagnostic",
                "isotope-resolved NH3",
                "electrical V-I-Q measurement",
            ],
            "calibration_did_not_use_only_nh3": False,
        },
    }


def check_source_locks(contract: dict[str, Any]) -> tuple[str, list[str], dict[str, str]]:
    failures: list[str] = []
    actual: dict[str, str] = {}
    for name, expected_path, expected_hash in (
        ("Hong input", HONG_INPUT, EXPECTED_HONG_SHA),
        ("Shao Fe(110) input", SHAO_INPUT, EXPECTED_SHAO_SHA),
    ):
        if not expected_path.exists():
            failures.append(f"{name} missing: {expected_path}")
            continue
        actual_hash = sha256(expected_path)
        actual[name] = actual_hash
        if actual_hash != expected_hash:
            failures.append(f"{name} hash mismatch")
    for key, expected_hash in (("hong_kinetic_input", EXPECTED_HONG_SHA), ("shao_fe110_kinetic_input", EXPECTED_SHAO_SHA)):
        declared = str(contract["source_locks"].get(key, {}).get("expected_sha256", "")).lower()
        if declared != expected_hash:
            failures.append(f"Contract source-lock declaration inconsistent for {key}")
    return ("PASS" if not failures else "BLOCKED", failures, actual)


def check_field(contract: dict[str, Any]) -> tuple[str, list[str], dict[str, float]]:
    field = contract["field_transfer"]
    required = ["pressure_pa", "gas_temperature_k", "reduced_field_td", "surface_field_v_per_a", "beta", "beta_evidence"]
    failures = [f"Missing field-transfer value: {name}" for name in required if not filled(field.get(name))]
    derived: dict[str, float] = {}
    if not failures:
        p, t, en = float(field["pressure_pa"]), float(field["gas_temperature_k"]), float(field["reduced_field_td"])
        e_bulk_v_per_a = en * 1e-21 * p / (K_B * t) / 1e10
        predicted_surface_field = float(field["beta"]) * e_bulk_v_per_a
        derived = {"gas_number_density_m3": p / (K_B * t), "e_bulk_v_per_a": e_bulk_v_per_a,
                   "beta_times_e_bulk_v_per_a": predicted_surface_field,
                   "declared_surface_field_v_per_a": float(field["surface_field_v_per_a"])}
        rel = abs(predicted_surface_field - derived["declared_surface_field_v_per_a"]) / max(abs(derived["declared_surface_field_v_per_a"]), 1e-30)
        derived["relative_field_closure_error"] = rel
        if rel > 0.01:
            failures.append("F_s != beta*E_bulk within 1%: revise field inputs or transfer relation")
    return ("PASS" if not failures else "BLOCKED", failures, derived)


def check_shared_boundary(contract: dict[str, Any]) -> tuple[str, list[str]]:
    boundary = contract["shared_reactor_boundary"]
    required = [
        "n2_mole_fraction", "h2_mole_fraction", "electron_density_cm3", "electron_density_closure",
        "residence_time_s", "inlet_outlet_operator_identical", "gas_volume_cm3",
        "catalyst_area_to_volume_cm_minus1", "roughness", "surface_site_normalization", "catalyst_inventory_basis",
    ]
    failures = [f"Missing shared-boundary item: {name}" for name in required if not filled(boundary.get(name))]
    if filled(boundary.get("n2_mole_fraction")) and filled(boundary.get("h2_mole_fraction")):
        if abs(float(boundary["n2_mole_fraction"]) + float(boundary["h2_mole_fraction"]) - 1.0) > 1e-9:
            failures.append("N2 and H2 mole fractions must close to one for the stated binary feed")
    return ("PASS" if not failures else "BLOCKED", failures)


def check_numerical(contract: dict[str, Any]) -> tuple[str, list[str]]:
    block = contract["hybrid_numerical_acceptance"]
    required = ["hybrid_input_frozen", "terminal_window_p0_passed", "nitrogen_conservation_passed",
                "hydrogen_conservation_passed", "surface_site_conservation_passed", "solver_log_or_artifact"]
    failures = [f"Numerical-acceptance item not supplied or not passed: {name}" for name in required if not filled(block.get(name))]
    return ("PASS" if not failures else "BLOCKED", failures)


def check_identifiability(contract: dict[str, Any]) -> tuple[str, list[str]]:
    block = contract["pathway_identifiability"]
    required = ["gas_and_surface_first_nh_fluxes_reported", "nh3_source_sink_fluxes_reported",
                "out_of_sample_observable", "calibration_did_not_use_only_nh3"]
    failures = [f"Identifiability item not supplied or not passed: {name}" for name in required if not filled(block.get(name))]
    return ("PASS" if not failures else "BLOCKED", failures)


def markdown_report(gates: list[tuple[str, str, list[str]]], derived: dict[str, float], source_hashes: dict[str, str]) -> str:
    all_pass = all(status == "PASS" for _, status, _ in gates)
    lines = [
        "# Hybrid Hong–Fe(110) comparison preflight",
        "",
        "**Status:** " + ("READY FOR A PRE-REGISTERED HYBRID RUN" if all_pass else "NOT AUTHORIZED FOR A GAS–SURFACE FLUX CLAIM"),
        "",
        "This record is a gate evaluator. It reports no hybrid simulation, surface flux, fitted enhancement factor, or catalyst-performance result.",
        "",
        "## Gate results",
        "",
        "| Gate | Status | Requirement |",
        "|---|---|---|",
    ]
    descriptions = {
        "G0": "frozen and verified source inputs",
        "G1": "independently evidenced field transfer with F_s = beta E_bulk closure",
        "G2": "shared CSTR and catalyst boundary",
        "G3": "hybrid numerical, elemental and site-balance acceptance",
        "G4": "pathway-identifying fluxes plus an out-of-sample observable",
    }
    for gate, status, failures in gates:
        detail = descriptions[gate] if not failures else "; ".join(failures)
        lines.append(f"| {gate} | **{status}** | {detail} |")
    lines.extend(["", "## Source-lock hashes", ""])
    for name, digest in source_hashes.items():
        lines.append(f"- {name}: `{digest}`")
    if derived:
        lines.extend(["", "## Field-transfer closure", ""])
        for key, value in derived.items():
            lines.append(f"- {key}: `{value:.12g}`")
    lines.extend([
        "",
        "## Use rule",
        "",
        "A pathway-comparison figure may be produced only after every G0–G4 status is PASS. Filling a missing number with a nominal, literature, or fitted-to-NH3 value is not a substitute for independent field and boundary evidence. Until then, the Fe(110) branch remains a source-validated structural alternative, not a gas–surface flux comparator.",
    ])
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--contract", type=Path, default=None, help="JSON contract; defaults to <output>/hybrid_comparison_contract.json")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    contract_path = args.contract or args.output / "hybrid_comparison_contract.json"
    if not contract_path.exists():
        contract_path.write_text(json.dumps(template(), indent=2) + "\n", encoding="utf-8")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))

    g0, g0_failures, hashes = check_source_locks(contract)
    g1, g1_failures, derived = check_field(contract)
    g2, g2_failures = check_shared_boundary(contract)
    g3, g3_failures = check_numerical(contract)
    g4, g4_failures = check_identifiability(contract)
    gates = [("G0", g0, g0_failures), ("G1", g1, g1_failures), ("G2", g2, g2_failures),
             ("G3", g3, g3_failures), ("G4", g4, g4_failures)]

    (args.output / "P18_hybrid_comparison_preflight.md").write_text(markdown_report(gates, derived, hashes), encoding="utf-8")
    (args.output / "P18_hybrid_comparison_preflight.json").write_text(
        json.dumps({"date": date.today().isoformat(), "gates": [{"gate": gate, "status": status, "failures": failures}
                   for gate, status, failures in gates], "derived": derived, "source_hashes": hashes}, indent=2) + "\n",
        encoding="utf-8",
    )
    with (args.output / "P18_gate_checklist.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["gate", "status", "blocking_items"])
        for gate, status, failures in gates:
            writer.writerow([gate, status, " | ".join(failures)])
    print("; ".join(f"{gate}={status}" for gate, status, _ in gates))
    print(f"contract={contract_path}")
    print(f"report={args.output / 'P18_hybrid_comparison_preflight.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
