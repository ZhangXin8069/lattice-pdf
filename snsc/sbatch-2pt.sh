#!/bin/bash
# ============================================================================
# Slurm 提交脚本 — 质子 2pt 关联函数蒸馏计算
# ============================================================================
#
# 目标: 完全复现 examples/donghx/2pt_proton_Cg5gmu_L24x72_mom2_zdir_dcu.py
#
# 功能: 从本征矢量 + Perambulator 出发,
#       通过蒸馏框架计算质子 2pt 关联函数 C₂(t_snk, t_src),
#       包含动量涂抹、VVV Baryon Block 构造、Wick 收缩、
#       宇称投影和边界条件符号修正。
#
# 输出文件 (命名与参考代码完全一致):
#   Raw 收缩矩阵:
#     twopt_slice_pp_Px{Px}Py{Py}Pz{Pz}_eginphase{mom_smear}{element}_contract_conf{conf_id}.npy
#     形状: (Nt, Nt, 4, 4)
#
#   Parity 投影后:
#     twopt_slice_pp_Px{Px}Py{Py}Pz{Pz}_eginphase{mom_smear}{element}_nopol_ss_conf{conf_id}.npy
#     形状: (Nt, Nt)
#
# Python 环境: miniconda3, conda env "zhangxin-snsc"
#
# 用法:
#   sbatch sbatch-2pt.sh                                          # 默认参数
#   sbatch --export=ALL,CONF_ID="46000" sbatch-2pt.sh              # 指定组态
#   sbatch --export=ALL,PZ_LIST="-2,-3" sbatch-2pt.sh              # 指定动量
#   sbatch --export=ALL,ELEMENT="_Cg5g3" sbatch-2pt.sh             # 指定插值算符
# ============================================================================

# ============================================================================
# Slurm 作业配置 — 单 CPU 环境
# ============================================================================
# Slurm 日志写入脚本所在目录的 logs/ 子目录。
# 请在 snsc/ 目录下提交:  cd snsc && sbatch sbatch-2pt.sh
# ============================================================================

#SBATCH --job-name=proton_2pt
#SBATCH --partition=math,cpu6248R,cpueicc,i72c512g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

#SBATCH --output=logs/proton_2pt_%j.out
#SBATCH --error=logs/proton_2pt_%j.err

# ============================================================================
# 环境设置
# ============================================================================

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
# 分析参数 — 与 2pt_proton_Cg5gmu_L24x72_mom2_zdir_dcu.py 完全一致
# ============================================================================

# --- 格点几何 (硬编码, 与参考代码一致) ---
NX="${NX:-24}"
NT="${NT:-72}"

# --- 本征矢量 — 与参考代码 eig_dir 一致 ---
NEV="${NEV:-100}"
NEV1="${NEV1:-100}"

# --- 组态 — 与参考代码 conf_id = sys.argv[1] 一致 ---
CONF_ID="${CONF_ID:-46000}"

# --- 动量 — 与参考代码 Px=0, Py=0, Pzlist=[-2,-3,-4,-5,-6] 一致 ---
PX="${PX:-0}"
PY="${PY:-0}"
# PZ 作为单动量, PZ_LIST 作为动量扫描 (默认扫描 z 方向)
PZ="${PZ:--2}"
PZ_LIST="${PZ_LIST:--2,-3,-4,-5,-6}"

# --- 动量涂抹 — 与参考代码 mom_smear=-2, momsmear_phase=2 一致 ---
MOM_SMEAR="${MOM_SMEAR:--2}"
MOM_SMEAR_PHASE="${MOM_SMEAR_PHASE:-2}"

# --- 插值算符 — 与参考代码 element="_Cg5g4" 一致 ---
# _Cg5g4 = Cγ₅γ₄ (production default)
# _Cg5g3 = Cγ₅γ₃
# _Cg5   = Cγ₅
ELEMENT="${ELEMENT:-_Cg5g4}"

# --- 数据路径 — 与参考代码一致 ---
# 本征矢量 (读取)
EIG_DIR="${EIG_DIR:-/public/group/lqcd/eigensystem/beta6.20_mu-0.2770_ms-0.2400_L24x72/${CONF_ID}/}"
# Perambulator (读取)
PERAM_U_DIR="${PERAM_U_DIR:-/public/home/sunp/sunpeng_new_disk/mom_smear_perambulators/beta6.20_mu-0.2770_ms-0.2400_L24x72/output_dir_data/mz2_my0_mx0/${CONF_ID}/}"
# 参考数据 (已有 donghx 结果, 用于对比验证)
CORR_NUCL_DIR="${CORR_NUCL_DIR:-/public/group/lqcd/donghx/2pt_Result/beta6.20_mu-0.2770_ms-0.2400_L24x72/momsmear-2z/${CONF_ID}/}"

# ============================================================================
# 构建参数并运行
# ============================================================================

ARGS=(
    --analysis-type proton-2pt
    --xp numpy
    --output-dir "${OUTPUT_DIR}"
    --Nt "${NT}"
    --Nx "${NX}"
    --Nev "${NEV}"
    --Nev1 "${NEV1}"
    --conf-start "${CONF_ID}"
    --conf-step 1
    --conf-num 1
    --Px "${PX}"
    --Py "${PY}"
    --Pz "${PZ}"
    --Pz-list "${PZ_LIST}"
    --mom-smear "${MOM_SMEAR}"
    --mom-smear-phase "${MOM_SMEAR_PHASE}"
    --element "${ELEMENT}"
    --eig-dir "${EIG_DIR}"
    --peram-u-dir "${PERAM_U_DIR}"
    --corr-nucl-dir "${CORR_NUCL_DIR}"
)

echo "==============================================" | tee -a "${LOG_FILE}"
echo "  质子 2pt 关联函数蒸馏计算" | tee -a "${LOG_FILE}"
echo "  (复现 2pt_proton_Cg5gmu_L24x72_mom2_zdir_dcu.py)" | tee -a "${LOG_FILE}"
echo "  开始时间: $(date)" | tee -a "${LOG_FILE}"
echo "  节点: $(hostname)" | tee -a "${LOG_FILE}"
echo "  conda: zhangxin-snsc" | tee -a "${LOG_FILE}"
echo "  ------------------------------------" | tee -a "${LOG_FILE}"
echo "  格点: ${NT}×${NX}³, Nev=${NEV}, Nev1=${NEV1}" | tee -a "${LOG_FILE}"
echo "  组态: ${CONF_ID}" | tee -a "${LOG_FILE}"
echo "  动量: Pz_list=${PZ_LIST}" | tee -a "${LOG_FILE}"
echo "  涂抹: mom_smear=${MOM_SMEAR}, phase=${MOM_SMEAR_PHASE}" | tee -a "${LOG_FILE}"
echo "  插值: ${ELEMENT}" | tee -a "${LOG_FILE}"
echo "  本征矢量: ${EIG_DIR}" | tee -a "${LOG_FILE}"
echo "  Peramb:    ${PERAM_U_DIR}" | tee -a "${LOG_FILE}"
echo "  参考对比:  ${CORR_NUCL_DIR}" | tee -a "${LOG_FILE}"
echo "  生成输出:  ${OUTPUT_DIR}/data/" | tee -a "${LOG_FILE}"
echo "  日志:      ${OUTPUT_DIR}" | tee -a "${LOG_FILE}"
echo "==============================================" | tee -a "${LOG_FILE}"

# 运行 (无缓冲输出)
python -u "${MAIN_PY}" "${ARGS[@]}" 2>&1 | tee -a "${LOG_FILE}"

EXIT_CODE=${PIPESTATUS[0]}

echo "==============================================" | tee -a "${LOG_FILE}"
echo "  结束时间: $(date)" | tee -a "${LOG_FILE}"
echo "  退出码: ${EXIT_CODE}" | tee -a "${LOG_FILE}"
echo "  输出目录: ${CORR_NUCL_DIR}" | tee -a "${LOG_FILE}"
echo "  日志目录: ${OUTPUT_DIR}" | tee -a "${LOG_FILE}"
echo "==============================================" | tee -a "${LOG_FILE}"

exit ${EXIT_CODE}
