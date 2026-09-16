#!/usr/bin/env bash
# ============================================================
# build.sh — ZDPlasKin 通用编译脚本（复制到任一构建目录即可用）
#
# 用法：
#   1) 把本脚本复制到含有输入文件的构建目录（该目录应有：
#      kinet.inp + bolsigdb.dat + preprocessor.exe + dvode_f90_m.F90
#      + bolsig 库(.lib/.dll) + 主程序 main_*.F90）
#   2) 在 Git Bash 中执行（可从任意 cwd 调用）：
#        bash build.sh            # 预处理 + 全量编译 + 链接
#        bash build.sh main       # 只重编译主程序（改了物理条件时用）
#        bash build.sh clean      # 删除 .o/.mod 中间产物
#
# 行为：
#   - 自定位：所有操作相对于脚本所在目录，与调用位置无关
#   - gfortran：优先 $GFORTRAN_DIR 环境变量 → PATH → 项目默认
#     F:\Softwares\Ming64\mingw64\bin，都找不到则报错退出
#   - 主程序探测：*.f90/*.F90 中名为 main* 且含 PROGRAM 语句者优先；
#     找不到或存在歧义时给出清晰错误
#   - preprocessor 两次交互自动应答（printf '.\n\n' |）
#   - 预处理后自动应用项目标准补丁（reac_rates use 列表补 density，
#     515 机理必需，65 机理自动跳过）；如目录存在 build_hook.sh 则
#     在编译前执行，用于目录特有的附加补丁
#   - 若目录已有 zdplaskin_m.F90 而无 kinet.inp/preprocessor.exe，
#     跳过预处理直接编译（适合只分发生成模块的场景）
#   - 编译全程带 -ffree-line-length-none（515 反应版内嵌长行必需）
#   - 结束后列出产物及大小；缺产物则非零退出
#
# 注意：路径含中文时请给路径加双引号。
# ============================================================
set -euo pipefail

# 防御（实测坑）：裸 bash（如从 cmd 直接调用、非登录模式）时 PATH 可能不含
# Git 自带的 /usr/bin，导致 dirname/sed/grep 等基础命令"command not found"。
# bash 启动后 /usr/bin 挂载已存在、只是没在 PATH 里——此处仅用内建命令补回。
for _d in /usr/bin /bin /usr/local/bin; do
  if [ -x "$_d/dirname" ]; then
    case ":$PATH:" in *":$_d:"*) ;; *) PATH="$_d:$PATH" ;; esac
  fi
done
export PATH
unset _d
# 另一防御：极简环境（无 TEMP/TMP）下 MinGW 版 gfortran 会退到 C:\Windows
# 写临时文件而报 Permission denied；给它 Windows 形式的临时目录。
# POSIX 形式的 TMPDIR 反而会误导 MinGW 运行时，故一并清除。
if [ -z "${TEMP:-}" ] && [ -z "${TMP:-}" ]; then
  export TEMP="$(cygpath -w /tmp 2>/dev/null || printf 'C:\\Windows\\Temp')"
  export TMP="$TEMP"
fi
case "${TMPDIR:-}" in /*) unset TMPDIR ;; esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
MODE="${1:-full}"

step() { echo ""; echo "=== $1 ==="; }
die()  { echo "错误: $1" >&2; exit 1; }

# ---------- 0. clean 模式 ----------
if [ "$MODE" = "clean" ]; then
  rm -f ./*.o ./*.mod
  echo "已删除 .o / .mod 中间产物"
  exit 0
fi

# ---------- 1. 定位 gfortran ----------
# 注意（本项目实测坑）：这个 MSYS2 版 gfortran 的驱动按 PATH 找配套工具
# (f951/as/ld)，用全路径调用而不把其 bin 目录加进 PATH 会"静默失败"
# （无报错、无产物、rc=1）。因此找到后一律把目录注入 PATH 再按名调用。
# 另一坑：GFORTRAN_DIR 若是 Windows 形式（F:\...），需先转 MSYS 形式——
# bash 启动后对 PATH 的修改不会做 Windows→POSIX 路径转换。
if [ -n "${GFORTRAN_DIR:-}" ] && command -v cygpath >/dev/null 2>&1; then
  case "$GFORTRAN_DIR" in *:\\*|*:/*) GFORTRAN_DIR="$(cygpath -u "$GFORTRAN_DIR")" ;; esac
fi
if [ -n "${GFORTRAN_DIR:-}" ] && [ -x "${GFORTRAN_DIR}/gfortran.exe" ]; then
  export PATH="${GFORTRAN_DIR}:$PATH"
  echo "gfortran: 来自 GFORTRAN_DIR ($GFORTRAN_DIR)"
elif command -v gfortran >/dev/null 2>&1; then
  echo "gfortran: 来自 PATH ($(command -v gfortran))"
elif [ -x "/f/Softwares/Ming64/mingw64/bin/gfortran.exe" ]; then
  export PATH="/f/Softwares/Ming64/mingw64/bin:$PATH"
  echo "gfortran: 来自默认路径 F:\\Softwares\\Ming64\\mingw64\\bin"
else
  die "找不到 gfortran。请安装 MinGW-w64/MSYS2，或设置 GFORTRAN_DIR 环境变量指向其 bin 目录"
fi
gfortran --version | head -1
# 编译器可用性冒烟：试编译一个最小程序，防"静默失败"类问题
_smoke_dir="$(mktemp -d)"
printf 'program t\nprint *,1\nend program t\n' > "$_smoke_dir/t.f90"
if ! ( cd "$_smoke_dir" && gfortran -c t.f90 ) 2>/dev/null; then
  rm -rf "$_smoke_dir"
  die "gfortran 试编译失败（静默失败类问题：请检查 PATH 中 gfortran 配套工具是否齐全）"
fi
rm -rf "$_smoke_dir"

# ---------- 2. 探测源文件 ----------
# dvode（大小写扩展名都接受）
DVODE=""
for f in dvode_f90_m.F90 dvode_f90_m.f90; do [ -f "$f" ] && DVODE="$f" && break; done
[ -n "$DVODE" ] || die "未找到 dvode_f90_m.F90（ODE 求解器源码，应从 ZDPlasKin 发行包复制）"

# bolsig 库：优先 x86_64_g 变体（本项目实测唯一可用），其次任意 bolsig*.lib，再次 *.a
BOLSIG_LINK=""
if [ -f "bolsig_x86_64_g.lib" ]; then
  BOLSIG_LINK="-lbolsig_x86_64_g"
else
  cand="$(ls bolsig*.lib 2>/dev/null | head -1 || true)"
  if [ -n "$cand" ]; then
    base="$(basename "$cand" .lib)"
    BOLSIG_LINK="-l${base}"
    echo "提示: 未找到 bolsig_x86_64_g.lib，改用 $cand（若链接失败请核对库变体）"
  else
    cand="$(ls bolsig*.a 2>/dev/null | head -1 || true)"
    [ -n "$cand" ] && BOLSIG_LINK="-l:$(basename "$cand")"
  fi
fi
[ -n "$BOLSIG_LINK" ] || die "未找到 bolsig 库文件（bolsig_x86_64_g.lib / bolsig*.lib / bolsig*.a）"

# 主程序：main* 且含 PROGRAM 语句者优先；其余含 PROGRAM 的候选报歧义错误
MAIN=""
mapfile -t mains < <(ls main*.F90 main*.f90 2>/dev/null || true)
prog_mains=()
for f in ${mains:+"${mains[@]}"}; do
  grep -qiE '^[[:space:]]*program[[:space:]]' "$f" && prog_mains+=("$f")
done
if [ "${#prog_mains[@]}" -eq 1 ]; then
  MAIN="${prog_mains[0]}"
elif [ "${#prog_mains[@]}" -gt 1 ]; then
  die "发现多个 main* 主程序：${prog_mains[*]}。请移走多余文件或重命名，保留唯一主程序"
else
  # 没有 main*，退而在全部 .f90 里找含 PROGRAM 的文件
  mapfile -t allsrc < <(ls *.F90 *.f90 2>/dev/null || true)
  cands=()
  for f in ${allsrc:+"${allsrc[@]}"}; do
    case "$f" in dvode_f90_m.*|zdplaskin_m.*) continue;; esac
    grep -qiE '^[[:space:]]*program[[:space:]]' "$f" && cands+=("$f")
  done
  if [ "${#cands[@]}" -eq 1 ]; then
    MAIN="${cands[0]}"
    echo "提示: 主程序未以 main 开头，自动选中 $MAIN"
  elif [ "${#cands[@]}" -eq 0 ]; then
    die "未找到含 PROGRAM 语句的主程序 .f90/.F90 文件"
  else
    die "发现多个 PROGRAM 文件：${cands[*]}。请将主程序重命名为 main_*.F90 或移走多余文件"
  fi
fi
EXE="$(basename "$MAIN" | sed 's/\.[Ff]90$//').exe"
echo "主程序: $MAIN → $EXE"

# ---------- 3. 预处理（kinet.inp → zdplaskin_m.F90） ----------
if [ "$MODE" = "full" ]; then
  if [ -f "kinet.inp" ] && [ -f "preprocessor.exe" ]; then
    step "1/2 预处理（preprocessor.exe）"
    # preprocessor 有两次交互：输出文件名("."=默认) + 回车退出
    printf '.\n\n' | ./preprocessor.exe kinet.inp
  elif [ -f "zdplaskin_m.F90" ] || [ -f "zdplaskin_m.f90" ]; then
    echo "提示: 无 kinet.inp 或 preprocessor.exe，跳过预处理，直接编译已有 zdplaskin_m"
  else
    die "既无 kinet.inp+preprocessor.exe，也无已生成的 zdplaskin_m.F90"
  fi
fi

ZD=""
for f in zdplaskin_m.F90 zdplaskin_m.f90; do [ -f "$f" ] && ZD="$f" && break; done
[ -n "$ZD" ] || die "未找到 zdplaskin_m.F90（预处理产物）"

# ---------- 3b. 标准补丁：reac_rates 的 use 列表加入 density ----------
# 515 反应机理的嵌入 Fortran（kinet.inp 里 $ 行）在 reac_rates 子程序内
# 引用 density(...)，而 preprocessor 生成的 use 列表不含 density，
# 直接编译会报 "Function 'density' has no IMPLICIT type"。
# 本项目全部 515 版 rebuild.sh 都带这一补丁；65 版无此行，自动跳过。
if grep -qE "lreaction_block, rrt[[:space:]]*$" "$ZD"; then
  sed -i -E 's/^([[:space:]]*)lreaction_block, rrt[[:space:]]*$/\1lreaction_block, rrt, density/' "$ZD"
  echo "已应用标准补丁: reac_rates use 列表加入 density"
fi

# ---------- 3c. 目录自定义钩子（可选） ----------
# 若目录内存在 build_hook.sh，在编译前执行，用于目录特有的补丁/检查。
if [ -f "build_hook.sh" ]; then
  echo "执行目录钩子: build_hook.sh"
  bash build_hook.sh
fi

# ---------- 4. 编译与链接 ----------
FFLAGS="-O2 -ffree-line-length-none -static-libgfortran -static-libgcc"
if [ "$MODE" = "main" ]; then
  [ -f "zdplaskin_m.o" ] && [ -f "dvode_f90_m.o" ] || die "main 模式需要已有的 zdplaskin_m.o / dvode_f90_m.o，请先跑一次 bash build.sh"
  step "重编译主程序 $MAIN"
  gfortran $FFLAGS -c "$MAIN"
else
  step "2/2 编译（dvode + zdplaskin_m + 主程序）"
  gfortran $FFLAGS -c "$DVODE" "$ZD" "$MAIN"
fi
gfortran -O2 -o "$EXE" "$(basename "$EXE" .exe).o" zdplaskin_m.o dvode_f90_m.o \
  -L. $BOLSIG_LINK -static-libgfortran -static-libgcc
echo "编译完成 ✓"

# ---------- 5. 产物核对 ----------
step "产物清单"
missing=0
for f in "$ZD" "$EXE"; do
  if [ -f "$f" ]; then
    sz=$(stat -c '%s' "$f")
    printf "  %-28s %10d B\n" "$f" "$sz"
  else
    echo "  缺失: $f"
    missing=1
  fi
done
# 提醒运行期随行 dll
for dll in bolsig*.dll libquadmath-0.dll; do
  [ -f "$dll" ] || echo "  提醒: 缺少随行运行库 $dll（exe 运行时可能报缺库）"
done
[ "$missing" -eq 0 ] || die "产物不完整"
echo ""
echo "构建成功：$EXE"
