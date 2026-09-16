# PK1

本仓库按 Hong 机理的计算边界划分为三个模块，仅保存核心源码、程序包与交接文档。

## 目录说明

| 模块 | 内容 |
| --- | --- |
| `core/` | 必要程序包、ZDPlasKin 求解器源码、QtPlasKin、运行时、共享分析模块、模板和项目文档。 |
| `gas/` | Hong 气相反应计算，按 `0d/` 与 `0d-cstr/` 分开；每个目录包含 `HANDOFF.md`。 |
| `surface/` | Hong 表面反应及 Fe(110) 比较工具，按 `0d/` 与 `0d-cstr/` 分开；每个目录包含 `HANDOFF.md`。 |

每个计算目录的交接文档明确列出模型边界、入口脚本、输入位置和交接检查。研究方案、论文和项目级说明位于 `core/docs/`。

## 版本控制范围

仓库只跟踪集中维护的源码（Python、Fortran、Shell、批处理）、模板、许可证与 Markdown 文档。为避免提交体积过大或混入机器相关内容，`Gas/`、`GasReaction/`、`Hong/analysis/` 等运行工作区，以及案例构建目录中的数值结果、图像、二进制文件、编译产物、日志和缓存默认保留在本地，不纳入 Git。

`Hong/Reproduction/external_models/DFT-microkinetic/` 是一个已独立维护的上游 Git 仓库，因此不嵌套提交到本仓库。
