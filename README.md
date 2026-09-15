# PK1

本仓库用于保存 PK1 项目的核心源码、运行脚本、模板、许可证和 Markdown 文档。

## 目录说明

- `Hong/工具脚本/`：模型运行、分析、绘图与复现相关的 Python、Shell、批处理脚本及模板。
- `Hong/研究方案/`：研究路线、分析方案与结论性文档。
- `Hong/manuscript/`：论文草稿、证据映射、可复现性说明等 Markdown 文档。
- `Hong/Shao2024_JACSAu/`：Shao 2024 JACS Au 工作的复现脚本和说明文档。
- `Hong/1.ZDPlasKin/`、`Hong/qtplaskin-master/`：ZDPlasKin/QtPlasKin 相关源码与模板。

## 版本控制范围

仓库只跟踪集中维护的源码（Python、Fortran、Shell、批处理）、模板、许可证与 Markdown 文档。为避免提交体积过大或混入机器相关内容，`Gas/`、`GasReaction/`、`Hong/analysis/` 等运行工作区，以及案例构建目录中的数值结果、图像、二进制文件、编译产物、日志和缓存默认保留在本地，不纳入 Git。

`Hong/Reproduction/external_models/DFT-microkinetic/` 是一个已独立维护的上游 Git 仓库，因此不嵌套提交到本仓库。
