# Hong 表面反应：0D-CSTR 交接文档

## 范围

本目录包含 Hong 表面机理与 CSTR 传输项结合的工作流，以及 Fe(110) 表面支路的边界协调、拓扑审计和判别设计工具。

## 入口与职责

- `scripts/hong_cw_surface_cstr.py`：表面 0D-CSTR 基础模型。
- `scripts/hong_cw_surface_cstr_campaign.py`：表面 CSTR 的可恢复 campaign 与台账。
- `scripts/hybrid_comparison_preflight.py`：混合气相/表面比较的前置检查。
- `scripts/audit_external_dft_microkinetic.py`、`plot_fe110_surface_*.py`：Fe(110) 外部表面支路审计。
- `scripts/design_fe110_pathway_discrimination.py`、`plot_cross_model_boundary_protocol.py`：判别实验与边界协调设计。

## 运行与输入

运行前必须对齐 CSTR 停留时间、进/出料、面积体积比、表面位点密度、温度和外加场定义。外部 DFT-microkinetic 仓库保持独立版本控制，不嵌套提交。

## 交接检查

1. 同时验证元素守恒、表面位点守恒与 CSTR 稳态/末端窗口。
2. 未完成边界协调前，不报告气相与表面通量的可比排序。
3. 将新案例输出、图像和日志保留在本地未跟踪工作区。
