#!/bin/bash
# 自动更新依赖库
# 用法: ./update_deps.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "============================================"
echo "  更新依赖库 - $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================"

# --- LQCD_Master ---
echo ""
echo "[1/2] 更新 LQCD_Master ..."
if [ -d "LQCD_Master" ]; then
    cd "$SCRIPT_DIR/LQCD_Master"
    git pull --ff-only
    echo "  ✓ LQCD_Master 更新完成"
else
    echo "  ⚠ LQCD_Master 目录不存在，正在克隆..."
    cd "$SCRIPT_DIR"
    git clone https://github.com/sjtu-sai-agents/LQCD_Master.git
    echo "  ✓ LQCD_Master 克隆完成"
fi

# --- lamet-agent ---
echo ""
echo "[2/2] 更新 lamet-agent ..."
if [ -d "lamet-agent" ]; then
    cd "$SCRIPT_DIR/lamet-agent"
    git pull --ff-only
    echo "  ✓ lamet-agent 更新完成"
else
    echo "  ⚠ lamet-agent 目录不存在，正在克隆..."
    cd "$SCRIPT_DIR"
    git clone https://github.com/Greyyy-HJC/lamet-agent.git
    echo "  ✓ lamet-agent 克隆完成"
fi

echo ""
echo "============================================"
echo "  全部依赖库已更新完毕"
echo "============================================"
