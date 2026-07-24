#!/bin/bash
# ============================================================================
# Slurm 提交脚本 — LQCD_Master + lamet-agent 联合质子 2pt 蒸馏流水线
# ============================================================================
#
# 架构:
#   Phase 1-2: LQCD_Master Planner→Executor  (工作流规划 + Slurm 代码生成)
#   Phase 3:   snsc/main.py --analysis-type proton-2pt  (蒸馏计算)
#   Phase 4:   HDF5 Bridge  (.npy → lamet-agent 兼容格式)
#   Phase 5:   lamet-agent  (有效质量拟合 + 出版级作图)
#
# 等效于: snsc/sbatch-2pt.sh 的全部功能 + agent 框架集成
#
# 用法:
#   cd agent/snsc && sbatch sbatch.sh
#   sbatch --export=ALL,CONF_ID="46000" sbatch.sh
#   sbatch --export=ALL,PZ_LIST="-2,-3" sbatch.sh
#   sbatch --export=ALL,DISTILLATION_ONLY="1" sbatch.sh
#   sbatch --export=ALL,LAMET_ONLY="1",H5_PATH="artifacts/proton_2pt.h5" sbatch.sh
# ============================================================================

#SBATCH --job-name=agent_2pt
#SBATCH --partition=cpu6248R,cpueicc,i72c512g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

#SBATCH --chdir=/public/home/zhangxin/lattice-pdf/agent/snsc
#SBATCH --output=logs/agent_2pt_%j.out
#SBATCH --error=logs/agent_2pt_%j.err

# ============================================================================
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

source /public/home/zhangxin/miniconda3/etc/profile.d/conda.sh && conda activate zhangxin-snsc

# ============================================================================
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
OUTPUT_DIR="runs/${TIMESTAMP}"
LOG_FILE="${OUTPUT_DIR}/run.log"

mkdir -p "${OUTPUT_DIR}" "${OUTPUT_DIR}/data" "${OUTPUT_DIR}/artifacts" "${OUTPUT_DIR}/plots"

# ============================================================================
# 参数配置
# ============================================================================

CONF_ID="${CONF_ID:-46000}"
PZ_LIST="${PZ_LIST:--2,-3,-4,-5,-6}"
ELEMENT="${ELEMENT:-_Cg5g4}"
MOM_SMEAR="${MOM_SMEAR:--2}"
MOM_SMEAR_PHASE="${MOM_SMEAR_PHASE:-2}"
PZ="${PZ:--2}"
PX="${PX:-0}"
PY="${PY:-0}"

DISTILLATION_ONLY="${DISTILLATION_ONLY:-0}"
LAMET_ONLY="${LAMET_ONLY:-0}"
SKIP_PLAN="${SKIP_PLAN:-1}"  # 默认跳过 LQCD_Master 规划 (仅用蒸馏+lamet)
H5_PATH="${H5_PATH:-}"

# ============================================================================
# 构建参数
# ============================================================================

ARGS=(
    --conf-id="${CONF_ID}"
    --Pz-list="${PZ_LIST}"
    --element="${ELEMENT}"
    --output-dir="${OUTPUT_DIR}"
)

if [ "${SKIP_PLAN}" = "1" ]; then
    ARGS+=(--skip-plan)
fi
if [ "${DISTILLATION_ONLY}" = "1" ]; then
    ARGS+=(--distillation-only)
fi
if [ "${LAMET_ONLY}" = "1" ]; then
    ARGS+=(--lamet-only)
    [ -n "${H5_PATH}" ] && ARGS+=(--h5-path="${H5_PATH}")
fi

# ============================================================================
echo "==============================================" | tee -a "${LOG_FILE}"
echo "  LQCD_Master + lamet-agent 联合流水线" | tee -a "${LOG_FILE}"
echo "  (等效于 snsc/sbatch-2pt.sh + agent 分析)" | tee -a "${LOG_FILE}"
echo "  开始时间: $(date)" | tee -a "${LOG_FILE}"
echo "  节点: $(hostname)" | tee -a "${LOG_FILE}"
echo "  conda: zhangxin-snsc" | tee -a "${LOG_FILE}"
echo "  ------------------------------------" | tee -a "${LOG_FILE}"
echo "  组态: ${CONF_ID}" | tee -a "${LOG_FILE}"
echo "  动量: Pz_list=${PZ_LIST}" | tee -a "${LOG_FILE}"
echo "  涂抹: mom_smear=${MOM_SMEAR}, phase=${MOM_SMEAR_PHASE}" | tee -a "${LOG_FILE}"
echo "  插值: ${ELEMENT}" | tee -a "${LOG_FILE}"
echo "  模式: distillation_only=${DISTILLATION_ONLY}" | tee -a "${LOG_FILE}"
echo "        lamet_only=${LAMET_ONLY}" | tee -a "${LOG_FILE}"
echo "        skip_plan=${SKIP_PLAN}" | tee -a "${LOG_FILE}"
echo "  输出: ${OUTPUT_DIR}" | tee -a "${LOG_FILE}"
echo "==============================================" | tee -a "${LOG_FILE}"

python -u "./run.py" "${ARGS[@]}" 2>&1 | tee -a "${LOG_FILE}"

EXIT_CODE=${PIPESTATUS[0]}

echo "==============================================" | tee -a "${LOG_FILE}"
echo "  结束时间: $(date)" | tee -a "${LOG_FILE}"
echo "  退出码: ${EXIT_CODE}" | tee -a "${LOG_FILE}"
echo "  输出目录: ${OUTPUT_DIR}" | tee -a "${LOG_FILE}"
echo "==============================================" | tee -a "${LOG_FILE}"

exit ${EXIT_CODE}
