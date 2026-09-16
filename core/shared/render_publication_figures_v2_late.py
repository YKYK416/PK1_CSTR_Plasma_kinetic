#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Resume the Fe(110), pathway-concentration and graphical-abstract v2 render."""
from __future__ import annotations

import matplotlib.figure

import render_publication_figures_v2 as refresh
from publication_figure_theme import apply_theme, refine_figure


def main() -> int:
    if not refresh.OUT.is_dir():
        raise SystemExit(f"Missing early v2 render: {refresh.OUT}")
    apply_theme()
    original_savefig = matplotlib.figure.Figure.savefig
    def styled_savefig(fig, *args, **kwargs):
        refine_figure(fig)
        kwargs.setdefault("pad_inches", 0.035)
        return original_savefig(fig, *args, **kwargs)
    matplotlib.figure.Figure.savefig = styled_savefig
    try:
        late = refresh.OUT / "fe110_and_design"
        refresh.call_main(refresh.load("hong_mechanism_ensemble_audit.py"), ["--output", str(late / "figure20")])
        refresh.call_main(refresh.load("audit_external_dft_microkinetic.py"), ["--output", str(late / "figure21")])
        refresh.call_main(refresh.load("plot_cross_model_boundary_protocol.py"), ["--output", str(late / "figure22")])
        refresh.call_main(refresh.load("plot_fe110_surface_graph_audit.py"), ["--output", str(late / "figure23")])
        refresh.call_main(refresh.load("plot_fe110_surface_symmetry_orbits.py"), ["--output", str(late / "figure24")])
        refresh.call_main(refresh.load("design_fe110_pathway_discrimination.py"), ["--output", str(late / "figure25")])
        refresh.render_26()
        refresh.call_main(refresh.load("plot_graphical_abstract_scope.py"), ["--output", str(refresh.OUT)])
    finally:
        matplotlib.figure.Figure.savefig = original_savefig
    print(f"Late publication figure refresh written to: {refresh.OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
