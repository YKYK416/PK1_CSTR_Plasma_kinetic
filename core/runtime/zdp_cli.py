#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ZDPlasKin Python 原生命令行入口（不依赖 Git Bash）。"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import zdp_runtime as runtime


def _build_parser():
    parser = argparse.ArgumentParser(description="ZDPlasKin Python 原生构建/运行入口")
    sub = parser.add_subparsers(dest="command", required=True)
    build = sub.add_parser("build", help="预处理、编译与链接")
    build.add_argument("--directory", default=".", help="算例目录（默认当前目录）")
    build.add_argument("--mode", choices=("full", "main", "clean"), default="full")

    run = sub.add_parser("run", help="运行已编译算例")
    run.add_argument("--directory", default=".", help="算例目录（默认当前目录）")
    for name, default in (("en", "45.1"), ("en_off", "0.1"), ("tg", "300"),
                          ("tend", "1"), ("tag", None),
                          ("freq", "1000"), ("duty", "0.5"), ("cycles", "3"),
                          ("ne_mode", "fix"), ("tau_res", ""), ("atol", ""), ("rtol", ""), ("out_rel", None)):
        run.add_argument("--" + name.replace("_", "-"), default=default)
    run.add_argument("--n2-frac", "--n2frac", dest="n2frac", default="0.3333")
    run.add_argument("--cw", action="store_true")
    run.add_argument("--no-log", action="store_true")
    run.add_argument("--recompile", action="store_true")
    run.add_argument("--restore", action="store_true")
    return parser


def main(argv=None):
    args = _build_parser().parse_args(argv)
    try:
        if args.command == "build":
            return runtime.build(Path(args.directory), mode=args.mode)
        params = vars(args).copy()
        params.pop("command", None)
        params.pop("directory", None)
        return runtime.run(Path(args.directory), params)
    except RuntimeError as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
