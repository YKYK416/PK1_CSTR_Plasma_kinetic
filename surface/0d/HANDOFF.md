# Hong 表面反应：0D 交接文档

## 范围

本目录涵盖 Hong 2017/2018 校正机理中金属表面列的封闭 0D 连续波计算；不包含 CSTR 进/出料或停留时间项。

## 入口与职责

- `scripts/hong_cw_surface_longtime_scan.py`：闭合 0D 表面模型长时间扫描主程序。
- `scripts/hong_cw_surface_campaign.py`：可恢复的批量扫描与台账管理。
- `scripts/hong_cw_surface_rescue_probe.py`：受限资源下的诊断性复算。
- `scripts/analyze_hong_closed0d_campaign.py`：已完成 campaign 的汇总与绘图。
- `reproduction/shao2024/`：Shao 2024 JACS Au 的独立表面分支复现与来源说明。

## 运行与输入

求解器和通用运行时位于 `core/`。Hong 原始输入、编译目录和历史结果仍在本机 `Hong/` 工作区，Git 仅保存可维护脚本和说明。

## 交接检查

1. 确认运行边界为 closed 0D，而不是 CSTR。
2. 确认表面位点与元素守恒检查通过。
3. Shao/Fe(110) 分支只能作为表面竞争或结构比较证据，不能替代 Hong 气相路径的独立验证。
