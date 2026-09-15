# Python 原生构建与运行

本工作台的 GUI、`zdp_cli.py` 和 `zdp_runtime.py` 直接以 Python 子进程调用：

- `preprocessor.exe`（从 `kinet.inp` 生成 `zdplaskin_m.F90`）；
- `gfortran.exe`（编译并链接）；
- 生成的 `main_*.exe`（执行算例）。

因此，正常使用不需要安装 Git Bash，也不会调用 `build.sh` 或 `run.sh`。这两个文件如仍存在，仅作为旧工作台的人工参考，不参与 GUI 或 Python 命令行流程。

## 复制工作台

复制整个 `Hong` 文件夹即可保留 GUI、模板、发行包和复现实例。复制到新位置后：

1. 双击 `工具脚本\启动GUI.bat`。它相对自身目录启动 GUI，不依赖工作台绝对路径。
2. 在 GUI 的“发行包目录”中重新选择复制后的 `1.ZDPlasKin\ZDPlasKin_2.0a_Windows`（旧路径记录失效时会提示，不会被强制使用）。
3. 若 `gfortran.exe` 不在系统 `PATH`，在 GUI 环境设置中选择其 `bin` 目录，或设置环境变量 `GFORTRAN_DIR`。配置、`GFORTRAN_DIR`、`PATH`、常见 MinGW 目录依次尝试。

Python 负责调度而不是替代 Fortran 求解器：主计算速度仍主要由编译后的 Fortran 程序决定。建议使用 CPython 3.9 或更高版本；本机已用 CPython 3.13 验证。

## 无 GUI 命令

在 `工具脚本` 目录运行，例如：

```powershell
py -3 zdp_cli.py build --directory "..\Reproduction\2017Hong\build" --mode full
py -3 zdp_cli.py run --directory "..\Reproduction\2017Hong\build" --n2-frac 0.5 --tend 1e-6 --tag demo
```

常用可选项包括 `--en`、`--tg`、`--atol`、`--rtol`、`--out-rel`、`--recompile` 与 `--restore`。使用 `zdp_cli.py run --help` 可查看完整列表。
