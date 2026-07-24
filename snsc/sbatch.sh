#!/bin/bash
# ============================================================================
# Slurm 提交脚本 — 格点 QCD 质子非极化胶子 PDF 计算 (snsc/main.py)
# ============================================================================
#
# GPU 分区: nv100-ins, nv100-sug, dgx2, na100-ins, na100-sug, gpu-debug,
#            na100-40g, na800-sug, na800-pcie, h20-nettr
# CPU 分区: cpu6248R, cpueicc, i72c512g
#
# 若使用 GPU 分区, 需设置 --gpus 参数，并在下方修改 PARTITION 和 XP
#
# Python 环境: miniconda3, conda env "zhangxin-snsc"
#
# 用法:
#   sbatch sbatch.sh                                    # 使用默认参数
#   sbatch --export=ALL,STEPS="1,2,3,6" sbatch.sh       # 覆盖步骤参数
#   sbatch --export=ALL,XP="cupy",PARTITION="na100-sug" sbatch.sh  # GPU 模式
# ============================================================================

# ============================================================================
# Slurm 作业配置 — 默认单 CPU 环境
# ============================================================================
# Slurm 日志写入脚本所在目录的 logs/ 子目录。
# 请在 snsc/ 目录下提交:  cd snsc && sbatch sbatch.sh
# ============================================================================

#SBATCH --job-name=gluon_pdf
#SBATCH --partition=math,cpu6248R,cpueicc,i72c512g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

#SBATCH --output=logs/gluon_pdf_%j.out
#SBATCH --error=logs/gluon_pdf_%j.err

# ============================================================================
# 环境设置
# ============================================================================

# 限制线程数 (单线程作业, 防止 BLAS/LAPACK 多线程干扰误差估计)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# 激活 conda 环境
source /public/home/zhangxin/miniconda3/etc/profile.d/conda.sh && conda activate zhangxin-snsc

# ============================================================================
# 路径配置
# ============================================================================

# 脚本与输出目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN_PY="${SCRIPT_DIR}/main.py"

# 时间戳输出目录
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
OUTPUT_DIR="${SCRIPT_DIR}/output_${TIMESTAMP}"
LOG_FILE="${OUTPUT_DIR}/run.log"

# 创建输出目录
mkdir -p "${OUTPUT_DIR}"

# ============================================================================
# 参数配置 — 通过环境变量或在此处直接修改
# ============================================================================

# --- 后端选择 ---
# numpy (CPU) 或 cupy (GPU)
XP="${XP:-numpy}"

# --- 分析类型 ---
# pdf = 胶子 PDF 全流程, 2pt = 有效质量分析
ANALYSIS_TYPE="${ANALYSIS_TYPE:-pdf}"

# --- 步骤选择 (仅 pdf 模式) ---
# 逗号分隔的步骤编号, 或 "all" 运行全部 13 步
STEPS="${STEPS:-all}"

# --- 系综预设 ---
# L24x72, L32x64, L32x96, L36x108, L48x96, L48x144, L64x128
ENSEMBLE="${ENSEMBLE:-L32x64}"

# --- 数据类型 ---
# complex64 或 complex128
DTYPE="${DTYPE:-complex128}"

# --- 格点参数 (若不使用系综预设) ---
NT="${NT:-64}"
NX="${NX:-32}"
NEV="${NEV:-100}"
NEV1="${NEV1:-100}"
DELTA_Z="${DELTA_Z:-15}"
Z_DIR="${Z_DIR:-2}"

# --- 动量 ---
PX="${PX:-0}"
PY="${PY:-0}"
PZ="${PZ:-6}"
# 动量扫描列表 (逗号分隔), 留空则不扫描
PZ_LIST="${PZ_LIST:-}"

# --- 组态选择 ---
CONF_START="${CONF_START:-20000}"
CONF_STEP="${CONF_STEP:-50}"
CONF_NUM="${CONF_NUM:-1}"

# --- 动量涂抹 ---
MOM_SMEAR="${MOM_SMEAR:-3}"
MOM_SMEAR_PHASE="${MOM_SMEAR_PHASE:--3}"

# --- 规范组态文件路径 (ILDG 二进制) ---
GAUGE_FILE="${GAUGE_FILE:-}"

# --- 数据路径 (读取已有数据) ---
READ_2PT_DIR="${READ_2PT_DIR:-}"
READ_3PT_DIR="${READ_3PT_DIR:-}"
READ_VVV_DIR="${READ_VVV_DIR:-}"
READ_OPE_DIR="${READ_OPE_DIR:-}"

# --- 数据路径 (生成写入) ---
GEN_2PT_DIR="${GEN_2PT_DIR:-}"
GEN_3PT_DIR="${GEN_3PT_DIR:-}"
GEN_VVV_DIR="${GEN_VVV_DIR:-}"
GEN_OPE_DIR="${GEN_OPE_DIR:-}"

# --- 匹配参数 ---
ALPHA_S="${ALPHA_S:-0.2}"
MU_OVER_PZ="${MU_OVER_PZ:-1.0}"

# ============================================================================
# 构建命令行参数
# ============================================================================

ARGS=(
    --analysis-type "${ANALYSIS_TYPE}"
    --xp "${XP}"
    --steps "${STEPS}"
    --dtype "${DTYPE}"
    --output-dir "${OUTPUT_DIR}"
    --Px "${PX}"
    --Py "${PY}"
    --Pz "${PZ}"
    --conf-start "${CONF_START}"
    --conf-step "${CONF_STEP}"
    --conf-num "${CONF_NUM}"
    --delta-z "${DELTA_Z}"
    --z-dir "${Z_DIR}"
    --alpha-s "${ALPHA_S}"
    --mu-over-pz "${MU_OVER_PZ}"
)

# --- 系综预设 (可选) ---
if [ -n "${ENSEMBLE}" ]; then
    ARGS+=(--ensemble "${ENSEMBLE}")
fi

# --- 手动格点参数 (不使用系综时) ---
if [ -z "${ENSEMBLE}" ]; then
    ARGS+=(
        --Nt "${NT}"
        --Nx "${NX}"
        --Nev "${NEV}"
        --Nev1 "${NEV1}"
        --mom-smear "${MOM_SMEAR}"
        --mom-smear-phase "${MOM_SMEAR_PHASE}"
    )
fi

# --- 动量扫描 ---
if [ -n "${PZ_LIST}" ]; then
    ARGS+=(--Pz-list "${PZ_LIST}")
fi

# --- 规范组态文件 ---
if [ -n "${GAUGE_FILE}" ]; then
    ARGS+=(--gauge-file "${GAUGE_FILE}")
fi

# --- 读取路径 ---
[ -n "${READ_2PT_DIR}" ] && ARGS+=(--read-2pt-dir "${READ_2PT_DIR}")
[ -n "${READ_3PT_DIR}" ] && ARGS+=(--read-3pt-dir "${READ_3PT_DIR}")
[ -n "${READ_VVV_DIR}" ] && ARGS+=(--read-VVV-dir "${READ_VVV_DIR}")
[ -n "${READ_OPE_DIR}" ] && ARGS+=(--read-ope-dir "${READ_OPE_DIR}")

# --- 生成路径 ---
[ -n "${GEN_2PT_DIR}" ] && ARGS+=(--gen-2pt-dir "${GEN_2PT_DIR}")
[ -n "${GEN_3PT_DIR}" ] && ARGS+=(--gen-3pt-dir "${GEN_3PT_DIR}")
[ -n "${GEN_VVV_DIR}" ] && ARGS+=(--gen-VVV-dir "${GEN_VVV_DIR}")
[ -n "${GEN_OPE_DIR}" ] && ARGS+=(--gen-ope-dir "${GEN_OPE_DIR}")

# ============================================================================
# 运行主脚本
# ============================================================================

echo "==============================================" | tee -a "${LOG_FILE}"
echo "  格点 QCD 非极化胶子 PDF 计算" | tee -a "${LOG_FILE}"
echo "  开始时间: $(date)" | tee -a "${LOG_FILE}"
echo "  节点: $(hostname)" | tee -a "${LOG_FILE}"
echo "  conda: zhangxin-snsc" | tee -a "${LOG_FILE}"
echo "  输出目录: ${OUTPUT_DIR}" | tee -a "${LOG_FILE}"
echo "  参数: ${ARGS[*]}" | tee -a "${LOG_FILE}"
echo "==============================================" | tee -a "${LOG_FILE}"

# 执行 (无缓冲输出, 便于日志实时查看)
python -u "${MAIN_PY}" "${ARGS[@]}" 2>&1 | tee -a "${LOG_FILE}"

EXIT_CODE=${PIPESTATUS[0]}

echo "==============================================" | tee -a "${LOG_FILE}"
echo "  结束时间: $(date)" | tee -a "${LOG_FILE}"
echo "  退出码: ${EXIT_CODE}" | tee -a "${LOG_FILE}"
echo "  输出目录: ${OUTPUT_DIR}" | tee -a "${LOG_FILE}"
echo "==============================================" | tee -a "${LOG_FILE}"

exit ${EXIT_CODE}
