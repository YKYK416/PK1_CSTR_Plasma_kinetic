#!/usr/bin/env bash
# run_swarm_all.sh - 对全部 DB 变体跑 BOLSIG- swarm 扫描
# 注意：求解器在 ASCII 暂存目录运行（中文路径会崩），结果复制回工作区
set -uo pipefail
WS_DB="/g/Kimi_project/ZDPlaskin模拟/Reproduction/2017Hong/db_uq"
SCRATCH="/g/Kimi_project/bolsig_uq"
EXE="/g/Kimi_project/ZDPlaskin模拟/1.ZDPlasKin/bolsigminus.exe"

for db in "$WS_DB"/dbs/*.dat; do
  tag=$(basename "$db" .dat)
  dst="$WS_DB/outputs/$tag"
  if [ -f "$dst/output.dat" ] && grep -q FINISHED "$dst/bolsiglog.txt" 2>/dev/null; then
    echo "SKIP $tag"; continue
  fi
  run="$SCRATCH/$tag"
  mkdir -p "$run" "$dst"
  cp "$db" "$run/bolsigdb.dat"
  cp "$WS_DB/scripts/swarm_script.txt" "$run/"
  cp "$EXE" "$run/"
  (cd "$run" && (echo "swarm_script.txt" | timeout 240 ./bolsigminus.exe > console.log 2>&1))
  if grep -q FINISHED "$run/bolsiglog.txt" 2>/dev/null && [ -f "$run/output.dat" ]; then
    cp "$run/output.dat" "$run/bolsiglog.txt" "$dst/"
    echo "OK   $tag"
  else
    cp "$run/console.log" "$run/bolsiglog.txt" "$dst/" 2>/dev/null
    echo "FAIL $tag"
  fi
done
