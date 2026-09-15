# ZDPlasKin 文献复现操作指南（Hong 2017 基线）

> 适用项目：`G:\Kimi_project\ZDPlaskin模拟`
> 本文档帮助你**独立操作**整条复现流水线：改条件 → 重编译 → 跑模拟 → 出图对比。
> 最后更新：2026-08-09

---

## 一、总体架构：一张图看懂流水线

```
                    ┌─────────────────────────────────────────────┐
                    │           输入（你编辑的部分）                │
                    │                                             │
                    │  kinet.inp        bolsigdb.dat   main_hong.F90│
                    │  (反应机理+物种)   (电子碰撞截面)  (运行条件主程序)│
                    └───────┬──────────────────┬──────────────────┘
                            │                  │
                 preprocessor.exe              │
                  (语法检查+代码生成)            │
                            │                  │
                            ▼                  │
                    zdplaskin_m.F90 ◄──────────┘
                    (生成的Fortran模块, 包含物种/反应表 + BOLSIG接口)
                            │
              gfortran 编译 + 链接 bolsig_x86_64_g.lib
              (连同 dvode_f90_m.F90 刚性ODE求解器)
                            │
                            ▼
                      main_hong.exe ──运行──► hong_output_*.csv (时间演化)
                                              hong_rates_*.csv  (反应速率)
                            │
                            ▼
              make_figures.py (Python 后处理) ──► fig_*.png 对比图
```

**核心思想**：ZDPlasKin 把"反应机理"和"求解器"分离。你只需维护三个输入文件，编译一次得到一个可执行程序，之后改条件（N₂比例、时间）只需命令行传参，**不用重新编译**；只有改机理或改物理条件（E/N、温度、压力）才需要重新走一遍流程。

---

## 二、目录结构与文件角色

```
G:\Kimi_project\ZDPlaskin模拟\
├── Reproduction\2017Hong\
│   ├── kinet.inp                    ★ 反应机理（你最常编辑）
│   ├── kinet.inp.bak                  修改前备份
│   ├── 复现结果报告.md                结果与偏差分析
│   ├── make_figures.py              ★ 后处理出图脚本
│   ├── data\
│   │   ├── bolsigdb.dat             ★ 截面数据库（SIGLO 合并，50 过程）
│   │   └── results\                 输出：fig_*.png、summary.json、CSV
│   └── build\                       ★ 构建与运行目录（所有命令在这里执行）
│       ├── preprocessor.exe           ZDPlasKin 预处理器
│       ├── dvode_f90_m.F90            DVODE 刚性 ODE 求解器（不动）
│       ├── zdplaskin_m.F90            ★ 由 preprocessor 生成（不要手改）
│       ├── main_hong.F90            ★ 主程序：设定物理条件
│       ├── main_hong.exe              编译产物，直接运行
│       ├── bolsigdb.dat               从 data\ 复制过来的截面库
│       ├── kinet.inp                  从上一级复制过来的机理
│       ├── bolsig_x86_64_g.dll/.lib   BOLSIG+ 库
│       └── libquadmath-0.dll          gfortran 运行库（必须和 exe 在一起）
└── 1.ZDPlasKin\                     软件原始发行包（不动）
```

**工作约定**：机理和截面库的"母版"在 `Reproduction\2017Hong\` 和 `data\`，改完后复制进 `build\` 再编译。

---

## 三、输入文件详解

### 3.1 kinet.inp —— 反应机理（四个段落）

```fortran
ELEMENTS
E N H S M          ! 元素表：E=电子, S=表面位点(虚拟), M=第三体(虚拟)
END

SPECIES
E N2 H2 N H NH NH2 NH3 N2H M
N2(v1) N2(v2)            ! 振动态：圆括号+数字
N2(A3) N2(B3) N2(a1)     ! N2 电子激发态
H2(b3) H2(B1)            ! H2 激发态
N2^+ N^+ H^+ ...         ! 离子用 ^+
S NS HS NHS NH2S         ! 表面物种（S 后缀=吸附态）
END

BOLSIG
N2                       ! 只列"基础气体物种"，不要列激发态/振动态
H2
N
set dbfile bolsigdb.dat  ! 指定截面数据库文件
END

REACTIONS
反应式 ! 速率表达式
END
```

**速率表达式语法**（Fortran 风格，`Tgas` 是内置变量=气体温度 K）：

```
N + NH => H + N2                                ! 5.0d-11              ← 常数 [cm³/s]
H + NH2 => H2 + NH                              ! 6.6d-11 * exp(-1840.0d0/Tgas)   ← Arrhenius
N + N + M => N2 + M                             ! 1.38d-33 * exp(502.978d0/Tgas)  ← 三体 [cm⁶/s]
```

**BOLSIG 链接反应**：电子碰撞反应不写数字速率，写 `! Bolsig 数据库过程名`，速率由 BOLSIG+ 在线求解：

```
E + N2 => E + N2(A3) ! Bolsig N2->N2(A3)
E + N2 => E + E + N2^+ ! Bolsig N2->N2^+     ← 电离：产物必须 2 个电子
```

**必须遵守的规则**（都是本项目踩过的坑）：

1. **元素守恒**：每个反应左右两边的 E/N/H/S/M 原子数必须相等。电离反应产物写 2 个 `E`。
2. **BOLSIG 过程名必须精确匹配** `bolsigdb.dat` 里的产物名，否则初始化报 `cannot find processes link for <bolsig:XXX>`。改名字要两边一起改。
3. **激发态是"死端"会非物理累积**：每个激发态物种必须有衰减出路（辐射衰减/猝灭/解离），否则会无限堆积（本项目曾出现 N₂(A3) 累积到 0.36% 的 bug）。
4. **表面反应速率**：论文给的是粘附概率 γ，不是速率系数！换算公式：
   - 单位点：`k = γ·(v_th/4)·(A/V)/n_s,vol`  [cm³/s]
   - 双位点：`k = γ·(v_th/4)·(A/V)/n_s,vol²`  [cm⁶/s]
   - 本项目假设：A/V = 10 cm⁻¹，n_s = 1e15 cm⁻² → n_s,vol = 1e16 cm⁻³
5. **数据库没有的电子过程**：写固定速率并加注释说明来源与量级估计依据（见 kinet.inp 末尾的固定速率段）。

### 3.2 bolsigdb.dat —— 截面数据库

纯文本，每个过程一个块：

```
EXCITATION
N2 -> N2(A3)          ← 这就是 kinet.inp 要匹配的名字
6.17                    ← 阈值能 (eV)
SPECIES: e + N2
PROCESS: e + N2 -> e + N2(A3), Excitation
...
0.000000E+00  0.000000E+00    ← 能量(eV) | 截面(m²)
...
```

改过程名时**只改第二行**的 `A->B` 部分；阈值和数据行保持不动。N₂(A3) 原本在 SIGLO 库中分三个振动态分支，本项目已合并为单一过程（截面逐能点求和）。

### 3.3 main_hong.F90 —— 运行条件主程序

关键设定都在文件开头的 `ZDPlasKin_set_*` 调用里：

```fortran
! 总气体密度 [cm-3]：1 atm, 300 K（换压力/温度要重算这个数）
ntot = 2.446d19

! 物理条件：气体温度 + 约化电场（改 E/N 就改这里）
call ZDPlasKin_set_conditions(GAS_TEMPERATURE=300.0d0, REDUCED_FIELD=45.1d0)

! 初始密度
call ZDPlasKin_set_density('N2', n2_frac*ntot)         ! n2_frac 来自命令行
call ZDPlasKin_set_density('H2', (1.0d0-n2_frac)*ntot)
call ZDPlasKin_set_density('M',  ntot, ldens_const=.true.)   ! 第三体锁定
call ZDPlasKin_set_density('S',  1.0d16)                      ! 表面位点
call ZDPlasKin_set_density('E',  1.17d8, ldens_const=.true.)  ! 电子密度固定(论文输入)

! 求解精度
call ZDPlasKin_set_config(ATOL=1.0d3, RTOL=1.0d-4)
```

**命令行参数**（编译一次，随便跑）：

```
main_hong.exe  <n2_fraction>  <t_end_s>  <output_tag>
例：main_hong.exe 0.25 100 N2H2_1_3     → 输出 hong_output_N2H2_1_3.csv
```

**注意**：本复现中电子密度按论文做法**固定**（`ldens_const=.true.`），E/N 也是固定的，所以 Te 由 BOLSIG+ 在该 E/N 下确定（≈1.41 eV），全程不变。要模拟 E/N 随时间变化（脉冲/占空比），需要在时间循环里反复调用 `ZDPlasKin_set_conditions(REDUCED_FIELD=...)`。

---

## 四、完整操作流程（全部命令）

> **快捷方式**：`build\rebuild.sh` 已封装全流程——`./rebuild.sh`（完整重建+冒烟测试）、`./rebuild.sh main`（只重编译主程序）、`./rebuild.sh test`（只冒烟测试）。以下为手动分步命令，理解原理或排查问题时用。

所有命令在 **Git Bash** 中执行，工作目录为 `build\`。注意路径含中文，务必加引号。

### 4.0 环境准备（每次开新终端都要做）

```bash
# gfortran 在 F:\Softwares\Ming64\mingw64\bin，必须先加 PATH（否则静默失败！）
export PATH="/f/Softwares/Ming64/mingw64/bin:$PATH"
gfortran --version   # 验证：应显示 16.1.0
```

### 4.1 修改机理后：重新预处理 + 编译

```bash
cd "/g/Kimi_project/ZDPlaskin模拟/Reproduction/2017Hong"

# 1. 备份并编辑 kinet.inp（用任意编辑器）
cp kinet.inp kinet.inp.bak

# 2. 复制最新机理和截面库到 build/
cp kinet.inp build/kinet.inp
cp data/bolsigdb.dat build/bolsigdb.dat

cd build

# 3. 预处理：语法检查 + 生成 zdplaskin_m.F90
./preprocessor.exe kinet.inp
#    成功会打印物种数/反应数；有语法错误会报行号

# 4. 编译（三步：模块 → 主程序 → 链接）
gfortran -O2 -c dvode_f90_m.F90 zdplaskin_m.F90 main_hong.F90 -static-libgfortran -static-libgcc
gfortran -O2 -o main_hong.exe main_hong.o zdplaskin_m.o dvode_f90_m.o \
  -L. -lbolsig_x86_64_g -static-libgfortran -static-libgcc
```

### 4.2 只改物理条件（E/N/温度/压力/位点密度）：改 main_hong.F90 后只重编译主程序

```bash
cd "/g/Kimi_project/ZDPlaskin模拟/Reproduction/2017Hong/build"
gfortran -O2 -c main_hong.F90 -static-libgfortran -static-libgcc
gfortran -O2 -o main_hong.exe main_hong.o zdplaskin_m.o dvode_f90_m.o \
  -L. -lbolsig_x86_64_g -static-libgfortran -static-libgcc
```

### 4.3 只改 N₂ 比例 / 模拟时长：直接跑，不用编译

```bash
cd "/g/Kimi_project/ZDPlaskin模拟/Reproduction/2017Hong/build"

# 基线（论文 Fig.1：N2:H2=1:2，100 s）
./main_hong.exe 0.3333 100 N2H2_1_2

# 组成扫描（5 组，每组约 0.3 秒）
./main_hong.exe 0.5   10 N2H2_1_1
./main_hong.exe 0.333 10 N2H2_1_2
./main_hong.exe 0.25  10 N2H2_1_3
./main_hong.exe 0.2   10 N2H2_1_4
./main_hong.exe 0.75  10 N2H2_3_1
```

运行成功末尾会打印 `DONE` 和 NH₃/N/H 终态密度。输出两个 CSV：
- `hong_output_<tag>.csv`：时间、Te(eV)、E/N、功率分量 + 全部 35 物种密度（对数时间步，每 10 倍程 10 个点）
- `hong_rates_<tag>.csv`：终态各反应速率（cm⁻³s⁻¹），用于路径分析

### 4.4 出图

```bash
cd "/g/Kimi_project/ZDPlaskin模拟/Reproduction/2017Hong"
python make_figures.py
# 图和 summary.json 写入 data/results/
```

`make_figures.py` 从 `build/` 读 CSV，要新增曲线就改脚本里的 `show = [...]` 物种列表或 `runs = [...]` 运行标签列表。

---

## 五、验证与故障排查清单

### 每次跑完必查

| 检查项 | 正常表现 | 异常提示 |
|---|---|---|
| 初始化链接 | `species link 3 / process link 10` | `cannot find processes link` → bolsigdb.dat 过程名不匹配 |
| Te | 1–5 eV 且全程稳定 | 异常高/低 → E/N 或气体密度设错 |
| NH₃ 量级 | 1e14–1e16 cm⁻³ | 差数量级 → 检查表面速率换算、A/V、n_s |
| 激发态 | N₂(A3) 稳态 ≪ 基态 0.01% | 持续累积 → 该态缺衰减通道（死端） |
| DVODE | 无报错 | 崩溃 → 检查反应守恒、速率量级是否极端 |

### 本项目踩过的坑（避免重蹈）

1. **gfortran 不在 PATH 会静默失败**，必须先 `export PATH`。
2. **缺 libquadmath-0.dll**：exe 双击/运行报缺库，把 mingw64\bin 里的这个 dll 复制到 build/（已做好）。
3. **链接要用 `bolsig_x86_64_g.lib`**（带 `_g` 后缀），不是 `bolsig.lib`。
4. **编译必须加 `-static-libgfortran -static-libgcc`**，否则 gfortran 16 运行时会 0xC0000005 崩溃。
5. **BOLSIG 段只列基础物种**（N2/H2/N），给 N2(v1) 之类加 DENSITY 行会报错。
6. **M（第三体）必须同时声明为元素和物种**，并在主程序里用 `ldens_const=.true.` 锁定为总密度。
7. **BOLSIG 不支持 'a0' 解析拟合格式的截面**（Cosby N₂ 解离），只能用列表截面或固定速率。
8. **警告 "species not configured for BOLSIG+ exceeds 5e-1"** 来自虚拟物种 M，无害，可忽略。

---

## 六、做参数扫描的两种姿势

**E/N 扫描**（下一步推荐）：在 main_hong.F90 外层加一个循环，每次调用 `ZDPlasKin_set_conditions(REDUCED_FIELD=en)` 后重新积分，或每个 E/N 独立跑一遍（把 E/N 也做成命令行参数），然后 Python 汇总。

**占空比/脉冲**（论文 §3.4）：在时间积分循环里按方波切换 `REDUCED_FIELD`（放电相 45.1 Td / 余辉相 ~0 Td），跑多个周期到周期稳态。

---

## 七、关键参考文档

| 文档 | 内容 |
|---|---|
| `Reproduction/2017Hong/复现结果报告.md` | 机制修改清单、结果对比、偏差分析 |
| `复现状态报告_2026-08-09.md` | 项目当前状态与下一步 |
| `ZDPlasKin_BOLSIG_诊断与解决方案报告.md` | BOLSIG 链接问题的完整诊断 |
| `1.ZDPlasKin/ZDPlasKin_manual.txt` | ZDPlasKin 官方手册（API 全集） |
