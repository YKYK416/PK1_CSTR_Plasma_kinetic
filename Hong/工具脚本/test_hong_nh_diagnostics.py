#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Focused tests for isolated Hong diagnostic branch transformations."""
import csv
import tempfile
from pathlib import Path

import hong_nh_diagnostics as d

BASE = """H2(A3SIG) => H2 ! 1./(GAMMA_D+H2_WALL_SECOND_PART_E)\nN + H2(RYDBERG_SUM) => H + NH ! k\nN + H2(B3SIG) => H + NH ! k\nN + H2(B1SIG) => H + NH ! k\nN + H2(C3PI) => H + NH ! k\nN + H2(A3SIG) => H + NH ! k\nN(2D) + H2 => H + NH ! k\nN(2P) + H2 => H + NH ! k\n"""

for name, expected in (("no_rydberg_nh", 1), ("no_named_h2star_nh", 4), ("no_excited_n_nh", 2)):
    text, changes = d.branch_kinet_text(BASE, name)
    assert text.count("#DIAGNOSTIC_OFF#") == expected and len(changes) == expected
for factor in (1, 10, 100, 1000):
    text, changes = d.branch_kinet_text(BASE, f"rydberg_loss_x{factor}")
    assert f"{factor}.0d0/(GAMMA_D+H2_WALL_SECOND_PART_E)" in text and len(changes) == 1
    assert "H2(A3SIG) => H2 ! 1./(GAMMA_D+H2_WALL_SECOND_PART_E)" in text

with tempfile.TemporaryDirectory() as td:
    p = Path(td) / "rates.csv"
    with p.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time_s", "dt_s", "cycle", "phase", "reaction", "rate_cm-3s-1"])
        w.writerow(["1e-3", "1e-3", "1", "on", "N+H2(B3SIG)=>H+NH", "2e10"])
        w.writerow(["1e-3", "1e-3", "1", "on", "N(2P)+H2=>H+NH", "1e10"])
        w.writerow(["2e-3", "1e-3", "2", "on", "N+H2(B3SIG)=>H+NH", "4e10"])
    rows, start, end = d.integrate_nh_channels(p)
    by_name = {r["channel"]: r for r in rows}
    assert by_name["Named-H2*"]["phi_cm-3"] == 6e7
    assert by_name["N(2P)+H2"]["share_pct"] == 100.0 / 7.0
    assert start == 0.0 and end == 2e-3
    rows, start, end = d.integrate_nh_channels(p, first_cycle=2)
    by_name = {r["channel"]: r for r in rows}
    assert by_name["Named-H2*"]["phi_cm-3"] == 4e7
    assert start == 1e-3 and end == 2e-3

    cycles = Path(td) / "cycles.csv"
    with cycles.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["cycle", "stable_streak"])
        w.writeheader()
        for i in range(1, 9):
            w.writerow({"cycle": i, "stable_streak": max(0, i - 5)})
    assert d.steady_window_start(cycles) == (6, 8)

print("Hong NH diagnostic branch tests passed ✓")
