#!/usr/bin/env bash
# ============================================================
# run.sh — ZDPlasKin 通用运行脚本（复制到任一构建目录即可用）
#
# 用法（Git Bash，可从任意 cwd 调用）：
#   bash run.sh                          # 全默认：EN=45.1 Tg=300 n2=1/3 t=1 s
#   bash run.sh --en 60 --tg 400 --tend 1e-3 --tag test1
#   bash run.sh --en-list "30 60 90" --tag scan1   # E/N 扫描
#   bash run.sh --out runs/batch1/EN60 --tag EN60  # 自定义输出子目录（GUI 用）
#   bash run.sh --restore                # 恢复被 --recompile 补丁过的源码
#   bash run.sh --help
#
# 命名参数（全部有默认值）：
#   --en X        约化场 E/N (Td)              默认 45.1（脉冲版为 EN_on）
#   --en-off X    脉冲 off 相 E/N (Td)         默认 0.1（仅脉冲版）
#   --tg X        气体温度 (K)                 默认 300
#   --n2-frac X   N2 摩尔分数                  默认 0.3333
#   --tend X      模拟时长 (s)                 默认 1（pulse 型时长由 --cycles 决定，此项被忽略）
#   --tag S       输出标签                     默认 r<月日_时分秒>
#   --freq X      脉冲频率 (Hz)                默认 1000（仅脉冲版）
#   --duty X      占空比                       默认 0.5（仅脉冲版）
#   --cycles N    脉冲周期数                   默认 3（仅脉冲版）
#   --ne-mode S   电子密度模式                 默认 fix（仅脉冲版）
#   --cw          连续（恒定 E/N）运行：pulse 型强制 duty=1.0（off 相时长为 0，
#                 主程序跳过余辉相，等效连续放电）；覆盖 --duty。非 pulse 型
#                 目录收到此开关只提示、不影响运行（其它类型本就连续）
#   --atol X      DVODE 绝对容差（程序支持时）
#   --rtol X      DVODE 相对容差（程序支持时）
#   --en-list ".."  以空格分隔的 E/N 列表，逐点运行，tag 自动加 _EN<值>
#   --out PATH    输出子目录（相对构建目录，可多级，如 runs/batch/EN60）；
#                 缺省 runs/<tag>。GUI 用它实现 批次/参数组合 两级目录
#   --no-log      不把控制台输出重定向到 console.log（直接打印到 stdout）
#   --recompile   hong 型主程序专用：sed 补丁源码中的硬编码
#                 REDUCED_FIELD/GAS_TEMPERATURE 并触发重编译（见下）
#   --restore     撤销 --recompile 的源码补丁（从 .runshbak 恢复）
#
# 输出隔离：
#   每次运行在 ./runs/<tag>/（或 --out 指定目录）内进行，自动复制
#   bolsigdb.dat 与 *.DAT 等运行期输入；exe 与 dll 留在构建目录
#   （Windows 按 exe 所在目录解析 dll）。已有结果文件不会被覆盖。
#
# 重要限制（实测结论）：
#   本项目主程序的 E/N、气体温度、配比均不从 kinet.inp 读取，而是：
#   - tscan 型（main_tscan.exe）：全部走命令行参数，本脚本直接可控；
#   - pulse 型（main_pulse.exe）：全部走命令行参数，本脚本直接可控；
#   - autopulse 型（main_auto_pulse.exe，向导生成）：Tgas EN_on EN_off freq duty
#     cycles tag [atol rtol]——组成内置（无 --n2-frac）、电子密度恒定（无 --ne-mode）、
#     时长由 --cycles/--freq 决定（--tend 被忽略）；
#   - hong 型（main_hong.exe / main_full.exe）：只接受
#     n2_frac / t_end / tag 三个参数，E/N 与 Tgas 是 Fortran 硬编码
#     （REDUCED_FIELD=45.1d0、GAS_TEMPERATURE=300.0d0）。
#     对 hong 型使用 --en/--tg 时必须加 --recompile：脚本会先备份
#     源码为 <主程序>.runshbak，再用 sed 补丁并调用 build.sh 重编译。
#     用 --restore 可恢复原始源码。脚本绝不静默修改你的主程序。
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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ---------- 默认值 ----------
EN=45.1; EN_OFF=0.1; TG=300; N2FRAC=0.3333; TEND=1
TAG="r$(date +%m%d_%H%M%S)"
FREQ=1000; DUTY=0.5; CYCLES=3; NEMODE=fix
ATOL=""; RTOL=""; EN_LIST=""; OUTDIR=""
RECOMPILE=0; RESTORE=0; CW=0; NOLOG=0

usage() { sed -n '2,66p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0; }

# ---------- 参数解析 ----------
while [ $# -gt 0 ]; do
  case "$1" in
    --en)       EN="$2"; shift 2;;
    --en-off)   EN_OFF="$2"; shift 2;;
    --tg)       TG="$2"; shift 2;;
    --n2-frac)  N2FRAC="$2"; shift 2;;
    --tend)     TEND="$2"; shift 2;;
    --tag)      TAG="$2"; shift 2;;
    --freq)     FREQ="$2"; shift 2;;
    --duty)     DUTY="$2"; shift 2;;
    --cycles)   CYCLES="$2"; shift 2;;
    --ne-mode)  NEMODE="$2"; shift 2;;
    --cw)       CW=1; shift;;
    --atol)     ATOL="$2"; shift 2;;
    --rtol)     RTOL="$2"; shift 2;;
    --en-list)  EN_LIST="$2"; shift 2;;
    --out)      OUTDIR="$2"; shift 2;;
    --no-log)   NOLOG=1; shift;;
    --recompile) RECOMPILE=1; shift;;
    --restore)  RESTORE=1; shift;;
    --help|-h)  usage;;
    *) echo "未知参数: $1（--help 查看用法）" >&2; exit 2;;
  esac
done

die() { echo "错误: $1" >&2; exit 1; }

# ---------- --restore：恢复源码补丁 ----------
if [ "$RESTORE" -eq 1 ]; then
  found=0
  for bak in *.runshbak; do
    [ -f "$bak" ] || continue
    orig="${bak%.runshbak}"
    cp "$bak" "$orig"
    echo "已恢复: $bak → $orig（如需生效请重新 bash build.sh main）"
    found=1
  done
  [ "$found" -eq 1 ] || echo "没有找到 .runshbak 备份文件"
  exit 0
fi

# ---------- 探测主程序 exe ----------
EXE=""
mapfile -t exes < <(ls main_*.exe 2>/dev/null || true)
if [ "${#exes[@]}" -eq 1 ]; then
  EXE="${exes[0]}"
elif [ "${#exes[@]}" -gt 1 ]; then
  EXE="${exes[0]}"
  echo "提示: 存在多个 main_*.exe（${exes[*]}），使用 $EXE；如需其它请移走多余 exe"
else
  cand="$(ls *.exe 2>/dev/null | grep -vE '^(preprocessor|bolsigplus|bolsigminus|t)\.exe$' | head -1 || true)"
  [ -n "$cand" ] && EXE="$cand"
fi
[ -n "$EXE" ] || die "未找到主程序 exe（请先运行 bash build.sh）"
BASE="$(basename "$EXE" .exe)"
EXE_ABS="$SCRIPT_DIR/$EXE"   # 绝对路径：--out 多级子目录下也能调用

# ---------- Fortran 运行库 PATH ----------
# 实测坑：exe 在 runs/<tag>/ 子目录里运行，Windows 按 PATH 找
# libgfortran-5.dll / libquadmath-0.dll 等运行库；裸环境 PATH 没有
# gfortran bin 时会报 "error while loading shared libraries"。
if [ -n "${GFORTRAN_DIR:-}" ] && command -v cygpath >/dev/null 2>&1; then
  case "$GFORTRAN_DIR" in *:\\*|*:/*) GFORTRAN_DIR="$(cygpath -u "$GFORTRAN_DIR")" ;; esac
fi
if [ -n "${GFORTRAN_DIR:-}" ] && [ -x "${GFORTRAN_DIR}/gfortran.exe" ]; then
  export PATH="${GFORTRAN_DIR}:$PATH"
elif ! command -v gfortran >/dev/null 2>&1 && [ -x "/f/Softwares/Ming64/mingw64/bin/gfortran.exe" ]; then
  export PATH="/f/Softwares/Ming64/mingw64/bin:$PATH"
fi

# 对应源码（用于识别参数签名）
SRC=""
for f in "${BASE}.F90" "${BASE}.f90"; do [ -f "$f" ] && SRC="$f" && break; done

# ---------- 识别主程序类型（参数签名，来自本项目实测） ----------
# hong  : main_hong/main_full —— n2_frac t_end tag（E/N、Tgas 硬编码）
# tscan : main_tscan          —— n2_frac Tgas EN t_end tag [atol rtol]
# pulse : main_pulse          —— n2_frac Tgas EN_on EN_off freq duty cycles ne_mode tag [atol rtol]
# auto  : main_auto（向导生成）—— Tgas EN t_end tag [atol rtol]（组成内置）
# autopulse : main_auto_pulse（向导生成）—— Tgas EN_on EN_off freq duty cycles tag [atol rtol]
#             （组成内置、电子恒定；无 n2/ne_mode/tend）
FLAVOR="hong"
if [ -n "$SRC" ]; then
  # 注意顺序：MAIN_AUTO_PULSE 包含子串 MAIN_AUTO，且 autopulse 源码含 EN_on，
  # 必须最先检查，否则会被误判为 auto 或 pulse
  if grep -q "MAIN_AUTO_PULSE" "$SRC"; then
    FLAVOR="autopulse"
  elif grep -q "MAIN_AUTO" "$SRC"; then
    FLAVOR="auto"
  elif grep -q "EN_on" "$SRC"; then
    FLAVOR="pulse"
  elif grep -qE "EN_set|Tgas_K EN_Td" "$SRC"; then
    FLAVOR="tscan"
  fi
else
  case "$BASE" in
    *auto_pulse*|*autopulse*) FLAVOR="autopulse";;   # 先于 *pulse* / *auto*，防误匹配
    *pulse*) FLAVOR="pulse";;
    *tscan*) FLAVOR="tscan";;
    *auto*) FLAVOR="auto";;
  esac
  echo "提示: 未找到 $BASE 的源码，按 exe 名称猜测类型为 $FLAVOR"
fi
echo "主程序: $EXE（类型: $FLAVOR）"

# ---------- --cw：连续（恒定 E/N）模式 ----------
if [ "$CW" -eq 1 ]; then
  if [ "$FLAVOR" = "pulse" ] || [ "$FLAVOR" = "autopulse" ]; then
    # main_pulse 在 duty=1.0 时 t_off=0，off 相被 if(t_off>0) 守卫整体跳过，
    # 等效连续放电；无除零风险（相位积分仅在时长>0 时调用）。
    DUTY=1.0
    echo "连续模式: --cw 生效，强制 duty=1.0（覆盖 --duty；时长仍为 cycles/freq）"
  else
    echo "提示: --cw 仅对 pulse/autopulse 型主程序有意义（$FLAVOR 型本就是连续场），已忽略"
  fi
fi

# ---------- hong 型的硬编码参数处理 ----------
EN_CHANGED=0; TG_CHANGED=0
[ "$EN" != "45.1" ] && EN_CHANGED=1
[ "$TG" != "300" ] && TG_CHANGED=1
if [ "$FLAVOR" = "hong" ] && { [ "$EN_CHANGED" -eq 1 ] || [ "$TG_CHANGED" -eq 1 ]; }; then
  if [ "$RECOMPILE" -eq 0 ]; then
    echo ""
    echo "注意: $EXE 是 hong 型主程序，E/N 与气体温度是 Fortran 硬编码"
    echo "      （REDUCED_FIELD=45.1d0, GAS_TEMPERATURE=300.0d0），不接受命令行控制。"
    echo "      你指定的 --en $EN / --tg $TG 不会生效。"
    echo "      如需修改：加 --recompile 让脚本备份源码、sed 补丁并重编译；"
    echo "      或直接手改 $SRC 后运行 bash build.sh main。"
    exit 2
  fi
  [ -n "$SRC" ] || die "--recompile 需要主程序源码 ${BASE}.F90 在当前目录"
  [ -f "${SRC}.runshbak" ] || cp "$SRC" "${SRC}.runshbak"
  echo "已备份源码: ${SRC}.runshbak"
  # 补丁值必须"以数字开头"，避免误伤 REDUCED_FIELD=EN_now 这类变量引用
  nmatch=0
  if grep -qE "REDUCED_FIELD=[0-9][0-9.eEdD+-]*" "$SRC"; then
    sed -i -E "s/REDUCED_FIELD=[0-9][0-9.eEdD+-]*/REDUCED_FIELD=${EN}d0/" "$SRC"
    nmatch=$((nmatch+1))
    echo "已补丁: REDUCED_FIELD=${EN}d0"
  fi
  if grep -qE "GAS_TEMPERATURE=[0-9][0-9.eEdD+-]*" "$SRC"; then
    sed -i -E "s/GAS_TEMPERATURE=[0-9][0-9.eEdD+-]*/GAS_TEMPERATURE=${TG}d0/" "$SRC"
    nmatch=$((nmatch+1))
    echo "已补丁: GAS_TEMPERATURE=${TG}d0"
  fi
  [ "$nmatch" -gt 0 ] || die "sed 未匹配到 REDUCED_FIELD=/GAS_TEMPERATURE= 赋值，请手工修改 $SRC"
  echo "触发重编译..."
  if [ -f "build.sh" ]; then
    bash build.sh main
  else
    # 无 build.sh 时的回退编译（gfortran 配套工具按 PATH 解析，需先注入目录）
    if ! command -v gfortran >/dev/null 2>&1; then
      [ -x "${GFORTRAN_DIR:-/f/Softwares/Ming64/mingw64/bin}/gfortran.exe" ] \
        && export PATH="${GFORTRAN_DIR:-/f/Softwares/Ming64/mingw64/bin}:$PATH" \
        || die "找不到 gfortran，无法重编译"
    fi
    gfortran -O2 -ffree-line-length-none -static-libgfortran -static-libgcc -c "$SRC"
    gfortran -O2 -o "$EXE" "${BASE}.o" zdplaskin_m.o dvode_f90_m.o \
      -L. -lbolsig_x86_64_g -static-libgfortran -static-libgcc
  fi
fi

# ---------- 组命令行 ----------
build_args() {
  local en="$1"
  case "$FLAVOR" in
    tscan)
      args=("$N2FRAC" "$TG" "$en" "$TEND")
      ;;
    pulse)
      args=("$N2FRAC" "$TG" "$en" "$EN_OFF" "$FREQ" "$DUTY" "$CYCLES" "$NEMODE")
      ;;
    autopulse)
      # 无 n2/ne_mode/tend：组成内置、电子恒定、时长=cycles/freq
      args=("$TG" "$en" "$EN_OFF" "$FREQ" "$DUTY" "$CYCLES")
      ;;
    auto)
      args=("$TG" "$en" "$TEND")
      ;;
    hong)
      args=("$N2FRAC" "$TEND")
      ;;
  esac
}

run_one() {
  local en="$1" tag="$2"
  build_args "$en"
  # 容差参数是位置参数（tag 之后）：只给一个会错位，要么都给要么都不给
  local extra=()
  if [ -n "$ATOL" ] || [ -n "$RTOL" ]; then
    [ -n "$ATOL" ] && [ -n "$RTOL" ] || die "--atol 与 --rtol 必须同时给出（位置参数不能错位）"
    [ "$FLAVOR" = "hong" ] && die "hong 型主程序不读取容差参数"
    extra=("$ATOL" "$RTOL")
  fi
  # --out 指定时可多级（runs/<批次>/<参数组合>）；缺省 runs/<tag>（旧行为）
  local rundir="${OUTDIR:-runs/${tag}}"
  mkdir -p "$rundir"
  # 运行期输入复制进子目录（ZDPlasKin 按 cwd 读 bolsigdb.dat；熵机制读 *.DAT）
  [ -f "bolsigdb.dat" ] && cp "bolsigdb.dat" "$rundir/"
  for d in *.DAT; do [ -f "$d" ] && cp "$d" "$rundir/"; done
  echo ""
  echo "=== 运行: $EXE ${args[*]} ${tag} ${extra[*]:-}（目录 ${rundir}）==="
  if [ "$NOLOG" -eq 1 ]; then
    ( cd "$rundir" && "$EXE_ABS" "${args[@]}" "$tag" ${extra:+"${extra[@]}"} ) \
      || die "运行失败（--no-log 模式，输出见上方）"
  else
    # tee 同时保留完整 console.log 并把主程序进度实时回传给 GUI。
    # pipefail 已在脚本开头启用，因此求解器失败仍会使本函数失败。
    ( cd "$rundir" && "$EXE_ABS" "${args[@]}" "$tag" ${extra:+"${extra[@]}"} 2>&1 | tee console.log ) \
      || die "运行失败，查看 ${rundir}/console.log"
  fi
  echo "输出文件（${rundir}/）:"
  ( cd "$rundir" && { ls -l *.csv 2>/dev/null || true; ls -l *.txt 2>/dev/null || true; } )
}

# ---------- 执行（单点或 E/N 扫描） ----------
if [ -n "$EN_LIST" ]; then
  [ "$FLAVOR" = "hong" ] && [ "$RECOMPILE" -eq 0 ] && \
    die "hong 型主程序的 E/N 是硬编码，--en-list 需配合 --recompile（但每次只能补丁一个值；请改用循环多次调用 --recompile，或换 tscan/pulse 型主程序）"
  for e in $EN_LIST; do
    run_one "$e" "${TAG}_EN${e}"
  done
else
  run_one "$EN" "$TAG"
fi

echo ""
echo "完成。结果在 runs/ 下按 tag 分目录存放（或 --out 指定目录）；--restore 可撤销源码补丁。"
