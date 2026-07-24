#!/bin/bash
# ============================================================================
# Slurm 提交脚本 — 2pt 有效质量分析 (snsc/main.py --analysis-type 2pt)
# ============================================================================
#
# 目标: 复现 examples/zhangxin/main-2pt.py 的结果
#
# 功能: 读取 Chroma IOG 格式的 2pt 关联函数,
#       计算 cosh 有效质量 (考虑周期性边界条件),
#       进行 Jackknife 误差估计, 并生成出版级有效质量图。
#
# 有效质量公式 (cosh):
#   a·m_eff(t) = arccosh( [C(t-1) + C(t+1)] / [2 C(t)] )
#   m_eff(t) [GeV] = a·m_eff(t) × 0.1973 / a[fm]
#
# Jackknife 误差:
#   C_i^{JK} = (Σ_j C_j - C_i) / (N-1)
#   σ^{JK}(t) = std(C^{JK}(t)) × √(N-1)
#
# 数据源: Chroma IOG 二进制文件
#   依赖: examples/zhangxin/iog_reader/iog.so (已编译的 C 扩展)
#
# Python 环境: miniconda3, conda env "zhangxin-snsc"
#
# 用法:
#   sbatch sbatch-2pt.sh                                      # 默认 (pion, P=0, L24x72)
#   sbatch --export=ALL,HADRON="proton" sbatch-2pt.sh          # 质子分析
#   sbatch --export=ALL,NCONF="100" sbatch-2pt.sh              # 100 个组态
#   sbatch --export=ALL,MEFF_TYPE="log" sbatch-2pt.sh          # 对数有效质量
# ============================================================================

# ============================================================================
# Slurm 作业配置 — 单 CPU 环境
# ============================================================================

#SBATCH --job-name=meff_2pt
#SBATCH --partition=cpu6248R,cpueicc,i72c512g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# ============================================================================
# 环境设置
# ============================================================================

# 限制线程 (单线程作业, 防止 BLAS 多线程干扰 Jackknife 误差)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# 激活 conda 环境
source /public/home/zhangxin/miniconda3/etc/profile.d/conda.sh && conda activate zhangxin-snsc

# ============================================================================
# 路径配置
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN_PY="${SCRIPT_DIR}/main.py"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
OUTPUT_DIR="${SCRIPT_DIR}/output_${TIMESTAMP}"
LOG_FILE="${OUTPUT_DIR}/run.log"

mkdir -p "${OUTPUT_DIR}"

# ============================================================================
# 分析参数 — 匹配 main-2pt.py
# ============================================================================

# --- 格点几何 ---
# main-2pt.py: Nx=24, Nt=72, alttc=0.1053
# 对应系综: beta=6.20, mu_l=-0.2770, mu_s=-0.2400, L=24³×72
NX="${NX:-24}"
NT="${NT:-72}"
ALTTC="${ALTTC:-0.1053}"

# --- 强子与有效质量类型 ---
HADRON="${HADRON:-pion}"
MEFF_TYPE="${MEFF_TYPE:-cosh}"

# --- 动量 ---
PX="${PX:-0}"
PY="${PY:-0}"
PZ="${PZ:-0}"

# --- 组态参数 ---
# main-2pt.py: N_start=10050, gap=50, Ncnfg_iog=52
CONF_START="${CONF_START:-10050}"
CONF_STEP="${CONF_STEP:-50}"
NCONF="${NCONF:-52}"

# --- tsep 和 link_max (文件路径格式参数) ---
TSEP="${TSEP:-36}"
LINK_MAX="${LINK_MAX:-10}"

# --- 时间折叠 ---
TIME_FOLD=""
if [ "${TIME_FOLD_FLAG:-0}" = "1" ]; then
    TIME_FOLD="--time-fold"
fi

# --- 作图范围 (GeV) ---
MEFF_RANGE="${MEFF_RANGE:-0.0,1.0}"

# --- IOG 2pt 文件路径模板 (绝对路径, 来自集群文件系统) ---
# 默认: beta=6.20, L=24×72, sush 的 IOG 数据
# 格式说明符: %d 被 data_analyse 替换为 (Nx, Nt, Px, Py, Pz, ENV, conf)
IOG_2PT_PATH="${IOG_2PT_PATH:-/public/home/sush/3pt_distillation/analyse/pion/24x72/beta6.20_mu-0.2770_ms-0.2400_L%dx%d/sush_iog/pion_2pt_Px%dPy%dPz%d_ENV%d_conf%d_tsep-1_mass-0.2770.iog}"

# ============================================================================
# 构建参数并运行
# ============================================================================

ARGS=(
    --analysis-type 2pt
    --xp numpy
    --output-dir "${OUTPUT_DIR}"
    --Nt "${NT}"
    --Nx "${NX}"
    --alttc "${ALTTC}"
    --hadron "${HADRON}"
    --meff-type "${MEFF_TYPE}"
    --Px "${PX}"
    --Py "${PY}"
    --Pz "${PZ}"
    --conf-start "${CONF_START}"
    --conf-step "${CONF_STEP}"
    --conf-num "${NCONF}"
    --tsep "${TSEP}"
    --link-max "${LINK_MAX}"
    --meff-range "${MEFF_RANGE}"
    --iog-2pt-path "${IOG_2PT_PATH}"
    ${TIME_FOLD}
)

echo "==============================================" | tee -a "${LOG_FILE}"
echo "  2pt 有效质量分析 (匹配 main-2pt.py)" | tee -a "${LOG_FILE}"
echo "  开始时间: $(date)" | tee -a "${LOG_FILE}"
echo "  节点: $(hostname)" | tee -a "${LOG_FILE}"
echo "  conda: zhangxin-snsc" | tee -a "${LOG_FILE}"
echo "  ------------------------------------" | tee -a "${LOG_FILE}"
echo "  格点: ${NT}×${NX}³,  a=${ALTTC} fm" | tee -a "${LOG_FILE}"
echo "  强子: ${HADRON},  有效质量: ${MEFF_TYPE}" | tee -a "${LOG_FILE}"
echo "  动量: P=(${PX},${PY},${PZ})" | tee -a "${LOG_FILE}"
echo "  组态: start=${CONF_START}, step=${CONF_STEP}, N=${NCONF}" | tee -a "${LOG_FILE}"
echo "  tsep=${TSEP}, link_max=${LINK_MAX}" | tee -a "${LOG_FILE}"
echo "  IOG: ${IOG_2PT_PATH}" | tee -a "${LOG_FILE}"
echo "  输出: ${OUTPUT_DIR}" | tee -a "${LOG_FILE}"
echo "==============================================" | tee -a "${LOG_FILE}"

# 运行 (无缓冲输出)
python -u "${MAIN_PY}" "${ARGS[@]}" 2>&1 | tee -a "${LOG_FILE}"

EXIT_CODE=${PIPESTATUS[0]}

echo "==============================================" | tee -a "${LOG_FILE}"
echo "  结束时间: $(date)" | tee -a "${LOG_FILE}"
echo "  退出码: ${EXIT_CODE}" | tee -a "${LOG_FILE}"
echo "  输出目录: ${OUTPUT_DIR}" | tee -a "${LOG_FILE}"
echo "  输出文件:" | tee -a "${LOG_FILE}"
ls -lh "${OUTPUT_DIR}/" | tee -a "${LOG_FILE}"
echo "==============================================" | tee -a "${LOG_FILE}"

exit ${EXIT_CODE}
