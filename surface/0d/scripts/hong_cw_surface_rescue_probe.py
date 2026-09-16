#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""One-case numerical rescue probe for the closed-0D Hong campaign.

This utility does not alter a campaign ledger or replace any production
attempt.  It puts each diagnostic trajectory in ``numerical_rescue_80Td``
under the selected campaign root, together with the exact solver parameters.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import hong_cw_surface_longtime_scan as model


def number_token(value: float) -> str:
    return f"{value:.0e}".replace("+", "").replace("-", "m")


def main() -> int:
    parser = argparse.ArgumentParser(description="Closed-0D Hong numerical rescue probe")
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--temperature-K", type=float, required=True)
    parser.add_argument("--en-Td", type=float, required=True)
    parser.add_argument("--n2-fraction", type=float, required=True)
    parser.add_argument("--atol", type=float, required=True)
    parser.add_argument("--rtol", type=float, default=1.0e-4)
    parser.add_argument("--outer-dt-s", type=float, default=0.1)
    parser.add_argument("--hmax-s", type=float, default=5.0e-3)
    parser.add_argument("--mxstep", type=int, default=100000)
    parser.add_argument("--t-end-s", type=float, default=20.0)
    parser.add_argument("--report-dt-s", type=float, default=5.0)
    parser.add_argument("--timeout-s", type=float, default=1800.0)
    parser.add_argument("--runtime-subdir", default="runtime_build_site_projection_v4",
                        help="campaign-root 下用于本次探针的独立已编译工作树目录")
    parser.add_argument("--label", default="", help="仅用于区分保留的诊断 attempt，不改变物理或数值输入")
    parser.add_argument("--site-projection", action="store_true",
                        help="仅在存档边界投影 5 个守恒表面位点物种；记录最大相对修正量")
    args = parser.parse_args()
    if not 0.0 < args.n2_fraction < 1.0:
        parser.error("n2-fraction 必须位于 0 与 1 之间")
    root = args.campaign_root.resolve()
    # The post-audit runtime includes explicit Fortran completion evidence.
    # ``prepare_worktree`` only touches this named copy; it never changes the
    # original campaign runtime or prior production attempts.
    work = model.prepare_worktree(root / args.runtime_subdir)
    tag = (f"T{args.temperature_K:g}_EN{args.en_Td:g}_N2{args.n2_fraction:g}_"
           f"atol{number_token(args.atol)}_dt{args.outer_dt_s:g}_to{args.t_end_s:g}")
    if args.label:
        tag += "_" + args.label.replace(" ", "_")
    output_dir = root / "numerical_rescue_80Td" / tag
    result = model.run_case(work, root, args.temperature_K, args.en_Td, args.n2_fraction,
                            args.outer_dt_s, args.t_end_s, args.report_dt_s,
                            "numerical_rescue", atol=args.atol, rtol=args.rtol,
                            internal_hmax_s=args.hmax_s, early_diagnostics=False,
                            internal_mxstep=args.mxstep, timeout_s=args.timeout_s,
                            site_projection=args.site_projection,
                            case_output_dir=output_dir)
    (output_dir / "rescue_probe_config.json").write_text(json.dumps({
        "purpose": "Numerical rescue diagnostic only; not a production campaign attempt.",
        "physical_model": "Hong 2017/2018 corrected surface mechanism; closed 0D; CW; no CSTR; no pulse.",
        "input": vars(args),
        "result": result,
    }, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    print(json.dumps({"output_dir": str(output_dir), "result": result}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
