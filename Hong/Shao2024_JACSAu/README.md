# Shao2024_JACSAu：独立 DFT–微观动力学基准

本目录保存 Ketong Shao 与 Ali Mesbah 公开的 ZDPlasKin DFT–微观动力学算例的 **Windows 可编译副本**，与现有 `Reproduction/2017Hong` 机制相互独立，不会覆盖或混合现有算例。

## 来源与可追溯性

- 论文：K. Shao, A. Mesbah, *A Study on the Role of Electric Field in Low-Temperature Plasma Catalytic Ammonia Synthesis via Integrated Density Functional Theory and Microkinetic Modeling*, **JACS Au** 4 (2024), 525–544；DOI: `10.1021/jacsau.3c00654`。
- 作者公开仓库：<https://github.com/wwwccttoo/DFT-microkinetic>（MIT License，副本见 `upstream_metadata/LICENSE`）。
- 提取版本：上游 `main` 提交 `54f1bfcf3e211e015f760085e6aa5e1bbc303743`。
- 归档 SHA-256：`5C7D78594C228C9B7C8085FAABFC028A069E6DA74336A0F71C633B22CC754603`。
- 提取路径：`Model_SA_Const_Entropy_extract`（作者标注为 DFT–microkinetic model）。

## 内容

`dft_microkinetic_windows_case/` 是可由本项目 Python 原生构建器编译的算例目录：

- `kinet_source.txt`：上游原始机理文件 `kinet_varyT_metal_auto_sens_DFT_in_entropy_verying_basis.txt`，SHA-256 为 `7CE49610345024CE7505433F27BF2DEA3AB07D139B6130470C0464FBD0893CC5`。
- `zdplaskin_m_DFT_in_entropy_varying_basis.F90`：作者随仓库提供的已生成 DFT 模块；`zdplaskin_m.F90` 为其同字节副本，供 Windows 构建器使用。
- `Const_E.F90` 及 `Ele.dat`、`Tgas.dat`、`other_para.dat`、反应能/熵参数文件：作者模型的运行输入。
- `preprocessor.exe`、`bolsig_x86_64_g.lib/.dll`：仅为 Windows 链接补齐的 ZDPlasKin 发行包文件，来源不同于作者仓库。

## 构建

在工作台根目录运行：

```powershell
py -3 工具脚本\zdp_cli.py build --directory Shao2024_JACSAu\dft_microkinetic_windows_case --mode full
```

已验证该命令在本机 CPython + MinGW gfortran 环境下生成 `Const_E.exe`。

## 重要限制

不要将 `kinet_source.txt` 直接改名为 `kinet.inp` 后使用当前 Windows `preprocessor.exe` 重预处理。该原始机理含 515 条反应及多条超过 130 字符的 `$` Fortran 记录；预处理器会截断这些记录，生成不可编译模块。该限制不影响本目录使用作者已生成的 `zdplaskin_m.F90` 编译。

该主程序从 `Ele.dat`、`Tgas.dat`、`other_para.dat` 等文件读取条件，内部默认积分到 `9000 s`；它不是当前 GUI 五种运行参数签名中的任一种。因此本次提取仅验证编译，不把它接入 GUI 的“运行/扫描”按钮，也不擅自改变论文的参数化方式。

与现有完整 Hong 机制的定量差异和新增反应清单见 [机理差异说明.md](机理差异说明.md)。
