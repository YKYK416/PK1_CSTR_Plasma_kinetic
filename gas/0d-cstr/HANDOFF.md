# Hong 气相反应：0D-CSTR 交接文档

## 范围

本目录是 Hong 气相机理的连续 0D-CSTR 工作流，包含进/出料与停留时间项。表面反应在该分支中明确关闭。

## 入口与职责

- `scripts/hong_cw_longtime_scan.py`：连续波气相 CSTR 长时间收敛与参数扫描主程序。
- `scripts/hong_cw_adaptive_cstr_scan.py`：可恢复的自适应 CSTR 扫描。
- `scripts/hong_cw_cstr_*.py`：收敛复核、统一容差和数值救援工具。
- `scripts/hong_control_switching_atlas.py`、`scripts/hong_tau_robustness.py`、`scripts/hong_h2star_envelope.py`：控制区、停留时间与情景包络分析。
- `scripts/hong_pulse_validation.py`、`scripts/hong_mechanism_ensemble_audit.py`：脉冲数值验收与机理族审计。
- `scripts/plot_hong_*cstr*.py` 及其他 `plot_*.py`：CSTR 数据与论文图后处理。

## 运行与输入

主运行输入和历史案例仍在本地 `GasReaction/` 目录，未提交到 Git。首次交接时应检查脚本中的默认扫描根目录，并使用参数覆盖为接收方机器上的实际位置。

运行前应确认：气相机理输入可编译、表面反应关闭、停留时间与入口条件已声明，并在结果中保留收敛/容差记录。

## 交接检查

1. 不要把 closed 0D 终态误报为 CSTR 出口量。
2. 对数值救援案例保留原始扫描登记与救援原因。
3. 新生成的案例、日志、二进制与图像继续保留本地，不直接提交。
