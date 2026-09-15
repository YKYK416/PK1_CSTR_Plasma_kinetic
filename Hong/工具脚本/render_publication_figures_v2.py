#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Render a non-destructive publication-style refresh of all manuscript figures.

This runner preserves each existing data-generating script and its accepted
source tables.  It applies only the shared publication theme at render time,
then writes a versioned figure set.  No mechanism input, numerical result or
manuscript image is overwritten.
"""
from __future__ import annotations

import importlib.util
import shutil
import sys
from pathlib import Path

import matplotlib.figure
import pandas as pd

from publication_figure_theme import apply_theme, refine_figure


TOOL_DIR = Path(__file__).resolve().parent
PROJECT = TOOL_DIR.parent
OUT = PROJECT / "analysis" / "publication_figures_20260903_v2"
OLD = PROJECT / "analysis" / "figures_20260829"


def load(filename: str):
    """Load a local generator and ensure its local setup inherits the theme."""
    path = TOOL_DIR / filename
    name = f"publication_refresh_{path.stem}"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    for attr in ("configure", "configure_style"):
        original = getattr(module, attr, None)
        if callable(original):
            def themed(*args, _original=original, **kwargs):
                result = _original(*args, **kwargs)
                apply_theme()
                return result
            setattr(module, attr, themed)
    return module


def call_main(module, args: list[str]) -> None:
    old_argv = sys.argv[:]
    try:
        sys.argv = [str(getattr(module, "__file__", "generator")), *args]
        result = module.main()
        if isinstance(result, int) and result != 0:
            raise RuntimeError(f"{module.__file__} returned {result}")
    finally:
        sys.argv = old_argv


def seed_plot_only_inputs() -> None:
    for name in ("Figure_3_dense_EN_N2_source_data.csv", "Figure_8_9_NHx_CSTR_summary.csv",
                 "Figure_8_NHx_reaction_turnover.csv"):
        source = OLD / name
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, OUT / name)


def render_3_from_accepted_table() -> None:
    """Render the dense map from its frozen accepted table, never rerun CSTR."""
    module = load("plot_hong_dense_cstr_map.py")
    apply_theme()
    table = pd.read_csv(OLD / "Figure_3_dense_EN_N2_source_data.csv")
    if len(table) != 108 or not table.run_status.isin(["success", "retry_success"]).all():
        raise RuntimeError("The Figure 3 source table is not the accepted 108-point CSTR grid")
    table.to_csv(OUT / "Figure_3_dense_EN_N2_source_data.csv", index=False)
    module.make_figure(table, OUT)


def render_17_18() -> None:
    tau = load("hong_tau_robustness.py")
    apply_theme()
    tau_table = pd.read_csv(PROJECT / "analysis" / "p8_tau_robustness_20260831_r3" / "P8_tau_cases.csv")
    tau.plot_results(tau_table, OUT)

    envelope = load("hong_h2star_envelope.py")
    apply_theme()
    envelope_table = pd.read_csv(PROJECT / "analysis" / "p9_h2star_scenario_envelope_20260831_r3" / "P9_scenario_cases.csv")
    envelope.plot_results(envelope_table, OUT)


def render_26() -> None:
    module = load("plot_direct_nh_pathway_concentration.py")
    module.OUT = OUT
    apply_theme()
    module.main()


def main() -> int:
    if OUT.exists() and any(OUT.iterdir()):
        raise SystemExit(f"Refusing to overwrite existing publication refresh: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    apply_theme()

    # Make every source generator use the same final visual cleanup immediately
    # before its normal PDF/PNG export.  This does not alter plotted objects.
    original_savefig = matplotlib.figure.Figure.savefig
    def styled_savefig(fig, *args, **kwargs):
        refine_figure(fig)
        kwargs.setdefault("pad_inches", 0.035)
        return original_savefig(fig, *args, **kwargs)
    matplotlib.figure.Figure.savefig = styled_savefig
    try:
        render_3_from_accepted_table()
        call_main(load("plot_hong_productivity_window.py"), ["--data", str(OUT / "Figure_3_dense_EN_N2_source_data.csv"), "--output", str(OUT)])
        call_main(load("plot_hong_cstr_figures.py"), ["--output", str(OUT)])
        call_main(load("plot_hong_graph_pathways.py"), ["--output", str(OUT)])
        seed_plot_only_inputs()
        call_main(load("plot_hong_nhx_turnover.py"), ["--output", str(OUT), "--plot-only"])
        call_main(load("plot_hong_reaction_control.py"), ["--output", str(OUT)])
        call_main(load("plot_hong_control_atlas_and_regimes.py"), ["--output", str(OUT)])
        call_main(load("plot_hong_extended_figure_candidates.py"), ["--output", str(OUT)])
        render_17_18()
        call_main(load("plot_hong_model_discrimination_atlas.py"), ["--output", str(OUT)])
        call_main(load("hong_mechanism_ensemble_audit.py"), ["--output", str(OUT)])
        call_main(load("audit_external_dft_microkinetic.py"), ["--output", str(OUT)])
        call_main(load("plot_cross_model_boundary_protocol.py"), ["--output", str(OUT)])
        call_main(load("plot_fe110_surface_graph_audit.py"), ["--output", str(OUT)])
        call_main(load("plot_fe110_surface_symmetry_orbits.py"), ["--output", str(OUT)])
        call_main(load("design_fe110_pathway_discrimination.py"), ["--output", str(OUT)])
        render_26()
        call_main(load("plot_graphical_abstract_scope.py"), ["--output", str(OUT)])
    finally:
        matplotlib.figure.Figure.savefig = original_savefig
    print(f"Publication figure refresh written to: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
