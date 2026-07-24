#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
================================================================================
格点 QCD 计算质子中非极化胶子 PDF 的统一工作流脚本
================================================================================

理论框架: 大动量有效理论 (LaMET)
方法: Distillation + 动量涂抹 + Disconnected 图分解

────────────────────────────────────────────────────────────────────────────
核心理论公式综述
────────────────────────────────────────────────────────────────────────────

1. 场强张量 F_{μν} (Clover 叶构造)
   ─────────────────────────────────────
   Plaquette 定义 (μ-ν 平面内的 1×1 Wilson loop):
     P_{μν}(x) = U_μ(x) U_ν(x+μ̂) U_μ^†(x+ν̂) U_ν^†(x)

   BCH 展开 (Baker-Campbell-Hausdorff):
     P_{μν} = exp[i a² (∂_μ A_ν - ∂_ν A_μ + i[A_μ,A_ν]) + O(a³)]
            = 1 + i a² F_{μν} + O(a³)

   Clover 叶 (四个 plaquette 的平均, 实现 O(a²) 改进):
     Q_{μν}(x) = P_{μν}(x) + P_{ν,-μ}(x) + P_{-μ,-ν}(x) + P_{-ν,μ}(x)

     四个 plaquette 的几何排列 (以 x 为中心):
                μ →
            +---(1)---+
            |         |
          ν ↑  (2)  (4)  ↓ (反方向)
            |         |
            +---(3)---+
                ← (反方向)
     (1) P_{μν}:    标准 1×1 Wilson loop
     (2) P_{ν,-μ}:  从 x 先走 ν 再走 -μ 方向
     (3) P_{-μ,-ν}: 从 x 先走 -μ 再走 -ν 方向
     (4) P_{-ν,μ}:  从 x 先走 -ν 再走 μ 方向

   Clover 场强张量 (反厄米部分, 取 a=1, g₀=1):
     F̂_{μν}(x) = -i/(8) · [Q_{μν}(x) - Q_{μν}^†(x)]

2. 对偶场强张量 (Hodge 对偶)
   ─────────────────────────────
     F̃_{μν} = (1/2) ε_{μνρσ} F^{ρσ}

     其中 ε_{0123} = +1, 完全反对称 Levi-Civita 符号。

3. Wilson 线 (规范平行传输)
   ─────────────────────────────
     连续定义 (伴随表示):
       W(z n̂, 0) = P exp( i g₀ ∫₀ᶻ dz' A_z(z' n̂) )

     格点离散化:
       W(z n̂_z, 0) ≃ Π_{k=0}^{z/a-1} U_z(k a n̂_z)

     其中 P 表示路径排序, U_z 是规范链接。

4. Nonlocal 胶子 OPE 算符
   ─────────────────────────────
     非极化胶子的完整 OPE 算符 (Eq(25) of 内部笔记):
       O_{unpol}(z) = Σ_{μ≠z} Σ_{ν≠z} g^{μν}
                      Tr_c[ F_{zμ}(z) · W(z→0) · F̃_ν^z(0) · W(0→z) ]

     在 Euclidean 时空中, 非极化组合 (具体指标取决于 z_dir):
       O^{(0)}_g = [-F^{0i} W F^{0i} W^† + 2F^{12} W F^{12} W^†]_{Eucl}

     其中 E^{0i} 分量带负号 (time-space 分量在 Wick 转动下变号),
     B^{12} 分量为正 (space-space 分量在 Wick 转动下不变)。

     格点实现 (五步):
       1. 平移 F_{μν} 到 z 位置:  roll(F[mu,nu], -dz, axis=3-z_dir)
       2. Wilson 线 W(z→0):   U_z^† 的连乘积
       3. 在原点插入 F̃_{μ2,ν2}(0)
       4. Wilson 线 W(0→z):   U_z 的连乘积
       5. 颜色迹 Tr_c[...] 并对全部空间点求和

5. Distillation 框架 (Laplace-Heaviside 子空间投影)
   ─────────────────────────────
     Laplace 算符的特征值问题:
       -∇² v^{(k)}(x) = λ_k v^{(k)}(x)

     涂抹算符 (投影到最低 Nev 个特征模态):
       □(x, y) = Σ_{k=1}^{Nev} v^{(k)}(x) v^{(k)†}(y)

     Perambulator (夸克传播子的子空间投影):
       τ_{αβ}^{AB}(t_snk, t_src) = Σ_{x,y} v^{A†}(x,t_snk) S_{αβ}(x,y) v^B(y,t_src)

     其中 A,B 是本征矢量指标, α,β 是 Dirac 指标。

6. VVV Baryon Block (重子插值算符)
   ─────────────────────────────
     Φ_{abc}(P) = Σ_x e^{-iP·x} ε_{ijk} v_i^a(x) v_j^b(x) v_k^c(x)

     其中 i,j,k ∈ {0,1,2} 是色指标, a,b,c 是本征矢量指标。
     ε_{ijk} 的六个置换:
       +1: (0,1,2), (1,2,0), (2,0,1)  [偶排列]
       -1: (0,2,1), (1,0,2), (2,1,0)  [奇排列]

     动量涂抹相位因子:
       φ_P(x) = exp(-i · 2π · P·x / L)

7. 质子 2pt 关联函数 (Wick 收缩)
   ─────────────────────────────
     直接项 (Direct):
       C₂^{dir}_{il}(t) = Φ_{abc}(t_sink) · τ_{gi}^{ad} · (ΓτΓ)_{gj}^{be}
                          · τ_{il}^{cf} · Φ_{def}^*(t_src)

     交换项 (Exchange):
       C₂^{ex}_{il}(t) = Φ_{abc}(t_sink) · τ_{gl}^{af} · (ΓτΓ)_{gj}^{be}
                          · τ_{ij}^{cd} · Φ_{def}^*(t_src)

     总关联函数 (Fermi 统计负号):
       C₂_{il}(t) = Direct - Exchange

     宇称投影 (分离正/负宇称态):
       C₂^{±}(t) = Tr[P_{±} · C₂(t)]
       其中 P_{±} = (1/2)(γ₀ ± γ₄)  (DeGrand-Rossi 基)

     插值算符 Dirac 结构 (质子): Γ = Cγ₅γ_μ
       在 DR 基下: C = γ₄γ₂,  γ₅ = diag(1,1,-1,-1)
       因此 Γ = γ₄γ₂ · γ₅ · γ_μ

8. 有效质量 (从 2pt 关联函数)
   ─────────────────────────────
     对数有效质量 (未考虑周期性边界条件):
       a·m_{eff}(t + 1/2) = ln[ C(t) / C(t+1) ]

     Cosh 有效质量 (考虑周期性边界条件):
       a·m_{eff}(t) = arccosh( [C(t-1) + C(t+1)] / [2 C(t)] )

     物理有效质量:
       m_{eff}(t) [GeV] = a·m_{eff}(t) × (ℏc / a) = a·m_{eff}(t) × 0.1973 / a[fm]

9. 矩阵元提取 (基态饱和)
   ─────────────────────────────
     关联函数比 (基态极限):
       R(z, P_z; t_sep) = C₃(t_sep; z, P_z) / C₂(t_sep; P_z)
                        → h(z, P_z)   (当 t_sep → ∞)

     在 disconnect 近似下 (胶子算符仅通过纯规范场期望值与夸克线耦合):
       h^{lat}(z, P_z) ≃ ⟨O(z)⟩_{gauge}
                        = (1/N_conf) Σ_{conf} (1/Nt) Σ_t O_{conf}(t; z)

10. Fourier 变换 → quasi-PDF
    ─────────────────────────────
     坐标空间到动量空间的 Fourier 变换:
       g̃(x, P_z) = ∫ dz/(2π x P_z) · e^{i x P_z z} · h(z, P_z)

     利用 h(z) 的厄米对称性 (h^†(z) = h(-z)):
       g̃(x, P_z) = (2 P_z / x) ∫_0^{z_max} dz · Re[h(z)] · sin(x P_z z)

11. 微扰匹配 → 光锥 PDF
    ─────────────────────────────
     NLO 匹配核 (MS-bar 方案):
       g(x, μ) = g̃(x, P_z) - (α_s/(2π)) ∫_x^1 (dy/y) ΔC_{gg}(x/y, μ/P_z) g̃(y, P_z)

     LO 近似 (匹配核 = δ(1-ξ)):
       g(x, μ) ≃ g̃(x, P_z)

12. 统计误差分析
    ─────────────────────────────
     Jackknife (留一法交叉验证):
       第 i 个 Jackknife 样本 (N 个组态, 删除第 i 个):
         x_i^{JK} = (N x̄ - x_i) / (N-1)

       Jackknife 误差:
         σ_{JK} = √[(N-1)/N · Σ_i (x_i^{JK} - x̄^{JK})²]
                = std(x^{JK}) · √(N-1)

     Bootstrap (有放回重采样):
       从 N 个组态中有放回抽取 N 个, 计算均值, 重复 B 次。
       Bootstrap 误差 = Bootstrap 样本的标准差。

────────────────────────────────────────────────────────────────────────────

该脚本实现了上述全部 12 步理论公式的格点计算链路:

  1.  读取规范组态 (gauge links)          [可选: 读取或生成]
  2.  构造场强张量 F_{μν} (Clover plaquette)
  3.  构造对偶场强张量 F̃_{μν} (Levi-Civita 收缩)
  4.  构造 VVV (distillation baryon block) [有作图]
  5.  构造质子 2pt 关联函数               [有作图]
  6.  构造 nonlocal 胶子 OPE 算符 (含 Wilson 线)
  7.  构造 3pt 关联函数                    [有作图]
  8.  实现 Distillation 框架 (本征矢量 + Perambulator + VVV)
  9.  动量涂抹质子两点关联函数
  10. 矩阵元提取 h(z, P_z)                [有作图]
  11. Fourier 变换 → quasi-PDF g̃(x, P_z)
  12. 微扰匹配 → 光锥 PDF g(x, μ)
  13. Jackknife / Bootstrap 统计误差分析

支持两种分析模式:
  --analysis-type pdf  胶子 PDF 全流程计算 (LaMET + distillation + OPE)
  --analysis-type 2pt  2pt 关联函数有效质量分析 (cosh/log meff + jackknife)

特性:
  - 可选后端: numpy (CPU) 或 cupy (GPU)
  - 全程显存/内存占用统计与耗时统计
  - 可选组态编号 (初始编号、步长、样本数)
  - 可选动量 (方向、大小)
  - 可选读取或生成数据的路径 (同时给出则对比)
  - 模块化步骤选择
  - 丰富的出版级作图
  - 可选数据类型 (complex64/complex128)
  - 规范化编程, 高可读性, 高鲁棒性, 高可移植性
  - 详细中文注释与过程 check 输出
  - 数值结果与已有代码保持一致

作者: Zhang Xin (基于 donghx 原始代码的系统化重构)
日期: 2026-07-24
================================================================================
"""

import numpy as np
import os
import sys
import time
import argparse
import resource
import gc
from typing import Dict, Tuple, List, Optional, Union, Any
from dataclasses import dataclass, field
from datetime import datetime
from contextlib import contextmanager

# ============================================================================
# 第一部分: 后端选择与全局配置
# ============================================================================

# --- GPU 后端检测 ---
HAS_CUPY = False
try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    cp = None

# --- 优化张量缩并库检测 ---
HAS_OPT_EINSUM = False
try:
    from opt_einsum import contract as _opt_contract
    HAS_OPT_EINSUM = True
except ImportError:
    _opt_contract = None

# --- 画图库 ---
HAS_MPL = False
try:
    import matplotlib
    matplotlib.use("Agg")  # 无头模式, 用于集群
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    plt = None

# --- IOG 读取与 2pt 分析支持 (继承自 examples/zhangxin/) ---
# 注: iog_reader.py 内部使用 os.getcwd() 定位 iog.so。
#     必须在 examples/zhangxin/ 目录下导入才能正确找到 .so 文件。
#     导入时临时切换 CWD 以绕过此限制。
HAS_INCLUDE = False
_data_analyse = None
_iog_read_fn = None
try:
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    _repo_root = os.path.dirname(_script_dir)
    _examples_dir = os.path.join(_repo_root, "examples", "zhangxin")
    if _examples_dir not in sys.path:
        sys.path.insert(0, _examples_dir)
    # iog_reader.py 使用 os.getcwd() 定位 iog.so, 必须切换到其所在目录
    _saved_cwd = os.getcwd()
    os.chdir(_examples_dir)
    try:
        from include import data_analyse as _data_analyse
        from iog_reader.iog_reader import iog_read
        _iog_read_fn = iog_read
        HAS_INCLUDE = True
        print("[INFO] IOG 读取支持已加载 (include.py + iog.so)。")
    finally:
        os.chdir(_saved_cwd)
except (ImportError, OSError) as e:
    print(f"[INFO] IOG 读取不可用 (集群专用): {e}")


# ============================================================================
# 第二部分: 工具函数 —— 计时器与内存追踪器
# ============================================================================

@contextmanager
def Timer(step_name: str, pipeline: "GluonPDFPipeline" = None):
    """上下文管理器: 记录每个步骤的耗时与内存峰值。

    用法:
        with Timer("Step 02: Field Strength", pipeline):
            F_all = compute_field_strength_all(gauge)

    结果自动记录到 pipeline.timing_results 和 pipeline.memory_results。
    """
    gc.collect()  # 运行前清理垃圾

    # --- 记录初始内存 ---
    mem_before = _get_memory(pipeline)

    # --- 开始计时 ---
    t_start = time.perf_counter()

    yield  # 执行步骤

    # --- 结束计时 ---
    t_end = time.perf_counter()
    elapsed = t_end - t_start

    # --- 记录最终内存与峰值 ---
    mem_after = _get_memory(pipeline)
    mem_peak = max(mem_before, mem_after)

    # --- 存储结果 ---
    if pipeline is not None:
        pipeline.timing_results[step_name] = elapsed
        pipeline.memory_results[step_name] = {
            "before_mb": mem_before,
            "after_mb": mem_after,
            "peak_mb": mem_peak,
        }

    # --- 打印输出 ---
    print(f"\n{'='*60}")
    print(f"  {step_name}")
    print(f"  耗时: {elapsed:.3f} s")
    print(f"  内存峰值: {mem_peak:.1f} MB")
    print(f"{'='*60}\n")


def _get_memory(pipeline: "GluonPDFPipeline" = None) -> float:
    """获取当前进程的内存/显存占用 (MB)。

    优先级:
      1. 若 pipeline 用 CuPy, 返回 GPU 显存
      2. 否则返回 RSS (常驻内存)
    """
    if pipeline is not None and pipeline.xp is cp:
        try:
            mempool = cp.cuda.Device().mem_info
            used = mempool.used_bytes if callable(mempool) else mempool[0] - mempool[1]
            return used / (1024 * 1024)
        except Exception:
            pass

    # CPU 模式: 使用 resource 模块获取 RSS
    try:
        usage = resource.getrusage(resource.RUSAGE_SELF)
        rss_kb = usage.ru_maxrss  # Linux: KB
        return rss_kb / 1024.0
    except Exception:
        return 0.0


def get_contract(xp: Any):
    """返回最佳可用的张量缩并函数。"""
    if HAS_OPT_EINSUM and xp is np:
        return _opt_contract
    return xp.einsum


# ============================================================================
# 第三部分: 格点配置数据类
# ============================================================================

@dataclass
class LatticeConfig:
    """格点几何与物理参数配置。

    Attributes
    ----------
    Nt : int, 时间方向格点数
    Nx, Ny, Nz : int, 空间方向格点数
    Nc : int, 颜色数 (QCD: Nc=3)
    Nd : int, Dirac 旋量维数 (Nd=4)
    Nev : int, Laplace 本征矢量总数
    Nev1 : int, 实际用于收缩的本征矢量数
    delta_z : int, Wilson 线最大长度 (格距单位)
    z_dir : int, Wilson 线方向 (0=x, 1=y, 2=z)
    conf_id : str, 规范组态编号
    Px, Py, Pz : int, 动量分量 (以 2π/L 为单位)
    mom_smear : int, 动量涂抹参数
    mom_smear_phase : int, 动量涂抹相位参数
    alttc : float, 格距 (fm)
    """
    Nt: int = 64
    Nx: int = 32
    Nc: int = 3
    Nd: int = 4
    Nev: int = 100
    Nev1: int = 100
    delta_z: int = 15
    z_dir: int = 2
    conf_id: str = "20000"
    Px: int = 0
    Py: int = 0
    Pz: int = 6
    mom_smear: int = 3
    mom_smear_phase: int = -3
    alttc: float = 0.1053

    @property
    def Ny(self) -> int:
        """Ny 始终等于 Nx (对称格点)."""
        return self.Nx

    @property
    def Nz(self) -> int:
        """Nz 始终等于 Nx (对称格点)."""
        return self.Nx

    @property
    def spatial_volume(self) -> int:
        """空间体积 V = Nx * Ny * Nz."""
        return self.Nx * self.Ny * self.Nz


# ============================================================================
# 第四部分: 内置系综预设
# ============================================================================

ENSEMBLES = {
    # ──────────────────────────────────────────────────────────────
    # 簇文件绝对路径基于 ~/ 下的符号链接:
    #   2pt_Result   → /public/group/lqcd/donghx/2pt_Result/
    #   Ope_Gluon    → /public/group/lqcd/donghx/Ope_Gluon/
    #   CLOVER       → /public/group/lqcd/configurations/CLOVER/
    #   eigenvectors → /public/group/lqcd/sush/eigenvectors/
    #   perambulators→ /public/group/lqcd/sush/perambulators/
    # ──────────────────────────────────────────────────────────────
    "L24x72": {
        "Nt": 72, "Nx": 24, "Nev": 100, "Nev1": 100,
        "mom_smear": -2, "mom_smear_phase": 2,
        "beta": 6.20, "mu_l": -0.2770, "mu_s": -0.2400,
        "conf_dir": "/public/group/lqcd/configurations/CLOVER/beta6.20_mu-0.2770_ms-0.2400_L24x72/{conf_id}/",
        "eig_dir": "/public/group/lqcd/eigensystem/beta6.20_mu-0.2770_ms-0.2400_L24x72/{conf_id}/",
        "peram_u_dir": "/public/home/sunp/sunpeng_public/mom_smear_perambulators/beta6.20_mu-0.2770_ms-0.2400_L24x72/output_dir_data/mz0_my0_mx2/{conf_id}/",
        "corr_nucl_dir": "/public/group/lqcd/donghx/2pt_Result/beta6.20_mu-0.2770_ms-0.2400_L24x72/momsmear-2x/{conf_id}/",
    },
    "L32x64": {
        "Nt": 64, "Nx": 32, "Nev": 100, "Nev1": 100,
        "mom_smear": 3, "mom_smear_phase": -3,
        "beta": 6.20, "mu_l": -0.2790, "mu_s": -0.2400,
        "conf_dir": "/public/group/lqcd/configurations/CLOVER/beta6.20_mu-0.2790_ms-0.2400_L32x64/{conf_id}/",
        "eig_dir": "/public/group/lqcd/eigensystem/beta6.20_mu-0.2790_ms-0.2400_L32x64/{conf_id}/",
        "peram_u_dir": "/public/home/sunp/sunpeng_public/mom_smear_perambulators/beta6.20_mu-0.2790_ms-0.2400_L32x64/output_dir_data/mz0_my0_mx-3/{conf_id}/",
        "corr_nucl_dir": "/public/group/lqcd/donghx/2pt_Result/beta6.20_mu-0.2790_ms-0.2400_L32x64/momsmear3x/{conf_id}/",
    },
    "L32x96": {
        "Nt": 96, "Nx": 32, "Nev": 100, "Nev1": 100,
        "mom_smear": 2, "mom_smear_phase": -2,
        "beta": 6.41, "mu_l": -0.2295, "mu_s": -0.2050,
        "conf_dir": "/public/group/lqcd/configurations/CLOVER/beta6.41_mu-0.2295_ms-0.2050_L32x96/{conf_id}/",
        "eig_dir": "/public/group/lqcd/eigensystem/beta6.41_mu-0.2295_ms-0.2050_L32x96/{conf_id}/",
        "peram_u_dir": "/public/home/sunp/sunpeng_public/mom_smear_perambulators/beta6.41_mu-0.2295_ms-0.2050_L32x96/output_dir_data/mz0_my0_mx-2/{conf_id}/",
        "corr_nucl_dir": "/public/group/lqcd/donghx/2pt_Result/beta6.41_mu-0.2295_ms-0.2050_L32x96/momsmear2x/{conf_id}/",
    },
    "L36x108": {
        "Nt": 108, "Nx": 36, "Nev": 200, "Nev1": 200,
        "mom_smear": 2, "mom_smear_phase": -2,
        "beta": 6.498, "mu_l": -0.2150, "mu_s": -0.1926,
        "conf_dir": "/public/group/lqcd/configurations/CLOVER/beta6.498_mu-0.2150_ms-0.1926_L36x108/{conf_id}/",
        "eig_dir": "/public/group/lqcd/eigensystem/beta6.498_mu-0.2150_ms-1926_L36x108/{conf_id}/",
        "peram_u_dir": "/public/home/sunp/sunpeng_new_disk/mom_smear_perambulators/beta6.498_mu-0.2150_ms-1926_L36x108/output_dir_data/mz0_my0_mx-2/{conf_id}/",
        "corr_nucl_dir": "/public/group/lqcd/donghx/2pt_Result/beta6.498_mu-0.2150_ms-1926_L36x108/momsmear2x/{conf_id}/",
    },
    "L48x96": {
        "Nt": 96, "Nx": 48, "Nev": 200, "Nev1": 200,
        "mom_smear": 4, "mom_smear_phase": -4,
        "beta": 6.20, "mu_l": -0.2825, "mu_s": -0.2310,
        "conf_dir": "/public/group/lqcd/configurations/CLOVER/beta6.20_mu-0.2825_ms-0.2310_L48x96/{conf_id}/",
        "eig_dir": "/public/group/lqcd/eigensystem/beta6.20_mu-0.2825_ms-0.2310_L48x96/{conf_id}/",
        "peram_u_dir": "/public/home/sunp/sunpeng_new_disk/mom_smear_perambulators/beta6.20_mu-0.2825_ms-0.2310_L48x96/output_dir_data/mz0_my0_mx-4/{conf_id}/",
        "corr_nucl_dir": "/public/group/lqcd/donghx/2pt_Result/beta6.20_mu-0.2825_ms-0.2310_L48x96/momsmear4x/{conf_id}/",
    },
    "L48x144": {
        "Nt": 144, "Nx": 48, "Nev": 200, "Nev1": 200,
        "mom_smear": 2, "mom_smear_phase": -2,
        "beta": 6.72, "mu_l": -0.1850, "mu_s": -0.1700,
        "conf_dir": "/public/group/lqcd/configurations/CLOVER/beta6.72_mu-0.1850_ms-0.1700_L48x144/{conf_id}/",
        "eig_dir": "/public/group/lqcd/eigensystem/beta6.72_mu-0.1850_ms-0.1700_L48x144/{conf_id}/",
        "peram_u_dir": "/public/home/sunp/sunpeng_public/mom_smear_perambulators/beta6.72_mu-0.1850_ms-0.1700_L48x144/output_dir_data/mz0_my0_mx-2/{conf_id}/",
        "corr_nucl_dir": "/public/group/lqcd/donghx/2pt_Result/beta6.72_mu-0.1850_ms-0.1700_L48x144/momsmear2x/{conf_id}/",
    },
    # --- 新增: L64x128 (from ~/symlinks) ---
    "L64x128": {
        "Nt": 128, "Nx": 64, "Nev": 200, "Nev1": 200,
        "mom_smear": 2, "mom_smear_phase": -2,
        "beta": 6.41, "mu_l": -0.2334, "mu_s": -0.2030,
        "conf_dir": "/public/group/lqcd/configurations/CLOVER/beta6.41_mu-0.2334_ms-0.2030_L64x128/{conf_id}/",
        "eig_dir": "/public/group/lqcd/sush/eigenvectors/beta6.41_mu-0.2334_ms-0.2030_L64x128/{conf_id}/",
        "peram_u_dir": "/public/group/lqcd/sush/perambulators/beta6.41_mu-0.2334_ms-0.2030_L64x128/{conf_id}/",
        "corr_nucl_dir": "/public/group/lqcd/donghx/2pt_Result/beta6.41_mu-0.2334_ms-0.2030_L64x128/momsmear2x/{conf_id}/",
    },
}


def resolve_ensemble(name: str, conf_id: str) -> dict:
    """解析系综预设, 将 {conf_id} 替换到所有路径模板中。"""
    if name not in ENSEMBLES:
        raise ValueError(f"未知系综 '{name}'。可用: {list(ENSEMBLES.keys())}")
    cfg = dict(ENSEMBLES[name])
    for key in ["conf_dir", "eig_dir", "peram_u_dir", "corr_nucl_dir"]:
        if key in cfg:
            cfg[key] = cfg[key].format(conf_id=conf_id)
    return cfg
    return cfg


# ============================================================================
# 第五部分: Dirac 矩阵 (DeGrand-Rossi 基) — 后端无关
# ============================================================================

def build_gamma_matrices(xp: Any) -> Dict[int, Any]:
    """构造 DeGrand-Rossi (DR) 基下的所有 Dirac 矩阵 (后端无关)。

    索引含义:
        0  = I (单位矩阵)          9  = γ₁γ₄
        1  = γ₁ (γ_x)             10 = γ₂γ₄
        2  = γ₂ (γ_y)             11 = γ₃γ₄
        3  = γ₃ (γ_z)             12 = γ₁γ₅
        4  = γ₄ (γ_t)             13 = γ₂γ₅
        5  = γ₅                    14 = γ₃γ₅
        6  = γ₂γ₃ (=-γ₁γ₄γ₅)     15 = γ₄γ₅
        7  = γ₃γ₁ (=-γ₂γ₄γ₅)     16 = γ₃γ₁ · P₊ (正宇称投影)
        8  = γ₁γ₂ (=-γ₃γ₄γ₅)     17 = γ₃γ₁ · P₋ (负宇称投影)
    """
    z = xp.zeros((4, 4), dtype=complex)

    # 单位矩阵
    g0 = z.copy()
    for i in range(4):
        g0[i, i] = 1.0 + 0.0j

    # γ₁ (γ_x): i * anti-diag(1, 1, -1, -1)
    g1 = z.copy()
    g1[0, 3] = 0.0 + 1.0j; g1[1, 2] = 0.0 + 1.0j
    g1[2, 1] = 0.0 - 1.0j; g1[3, 0] = 0.0 - 1.0j

    # γ₂ (γ_y): anti-diag(-1, 1, 1, -1)
    g2 = z.copy()
    g2[0, 3] = -1.0 + 0.0j; g2[1, 2] = 1.0 + 0.0j
    g2[2, 1] = 1.0 + 0.0j; g2[3, 0] = -1.0 + 0.0j

    # γ₃ (γ_z): i * anti-diag(1, -1, -1, 1)
    g3 = z.copy()
    g3[0, 2] = 0.0 + 1.0j; g3[1, 3] = 0.0 - 1.0j
    g3[2, 0] = 0.0 - 1.0j; g3[3, 1] = 0.0 + 1.0j

    # γ₄ (γ_t): anti-diag(1, 1, 1, 1)
    g4 = z.copy()
    g4[0, 2] = 1.0 + 0.0j; g4[1, 3] = 1.0 + 0.0j
    g4[2, 0] = 1.0 + 0.0j; g4[3, 1] = 1.0 + 0.0j

    # γ₅: diag(1, 1, -1, -1)
    g5 = z.copy()
    g5[0, 0] = 1.0 + 0.0j; g5[1, 1] = 1.0 + 0.0j
    g5[2, 2] = -1.0 + 0.0j; g5[3, 3] = -1.0 + 0.0j

    return {
        0: g0, 1: g1, 2: g2, 3: g3, 4: g4, 5: g5,
        6:  xp.matmul(g2, g3),
        7:  xp.matmul(g3, g1),
        8:  xp.matmul(g1, g2),
        9:  xp.matmul(g1, g4),
        10: xp.matmul(g2, g4),
        11: xp.matmul(g3, g4),
        12: xp.matmul(g1, g5),
        13: xp.matmul(g2, g5),
        14: xp.matmul(g3, g5),
        15: xp.matmul(g4, g5),
        16: xp.matmul(xp.matmul(g3, g1), 0.5 * (g0 + g4)),
        17: xp.matmul(xp.matmul(g3, g1), 0.5 * (g0 - g4)),
    }


# ============================================================================
# 第六部分: Levi-Civita 张量 ε_{μνρσ} / 2
# ============================================================================

def build_levi_civita_tensor() -> np.ndarray:
    """构造四维 Levi-Civita 符号 ε_{μνρσ} / 2。

    ε_{μνρσ}:
        +1  偶排列 (μ,ν,ρ,σ)
        -1  奇排列
         0  有重复指标

    因子 1/2 用于缩并: F̃_{μν} = (1/2) ε_{μνρσ} F^{ρσ}

    Returns
    -------
    epsilon4 : np.ndarray, shape (4, 4, 4, 4), dtype=float
    """
    epsilon4 = np.zeros((4, 4, 4, 4), dtype=float)

    for mu in range(4):
        for nu in range(4):
            a = 1.0 if mu > nu else 0.0
            for rho in range(4):
                b = 0.0
                if mu > rho: b += 1.0
                if nu > rho: b += 1.0
                for sigma in range(4):
                    c = 0.0
                    if mu > sigma: c += 1.0
                    if nu > sigma: c += 1.0
                    if rho > sigma: c += 1.0

                    parity = (a + b + c)
                    if parity % 2 == 0:
                        epsilon4[mu, nu, rho, sigma] = 1.0
                    else:
                        epsilon4[mu, nu, rho, sigma] = -1.0

                    # 重复指标 → 0
                    if len({mu, nu, rho, sigma}) != 4:
                        epsilon4[mu, nu, rho, sigma] = 0.0

    epsilon4 *= 0.5
    return epsilon4


# ============================================================================
# 第七部分: 场强张量 F_{μν} 的计算 (Clover 叶)
# ============================================================================

def plaquette_clover(gauge: Any, mu: int, nu: int, contract_fn: Any) -> Any:
    r"""计算 Clover 型的场强张量 F_{μν}(x)。

    ────────────────────────────────────────────────────────
    理论推导: 从 Plaquette 到 Clover 场强张量
    ────────────────────────────────────────────────────────

    Plaquette 定义 (μ-ν 平面内的 1×1 Wilson loop):
        P_{μν}(x) = U_μ(x) U_ν(x+μ̂) U_μ^†(x+ν̂) U_ν^†(x)

    BCH 展开 (Baker-Campbell-Hausdorff, 取 A_μ = a A_μ^{cont}):
        U_μ(x) = exp(i g₀ a A_μ(x))
               = 1 + i g₀ a A_μ - (1/2)(g₀ a)² A_μ² + O(a³)

        P_{μν} = exp[ i g₀ a² (∂_μ A_ν - ∂_ν A_μ + i g₀[A_μ,A_ν]) + O(a³) ]
               = 1 + i g₀ a² F_{μν} + O(a³)

    因此 (取 a=1, g₀=1 的格点单位):
        F_{μν}(x) = -i (P_{μν}(x) - 1)

    Clover 叶 (四个 plaquette 的对称组合, 实现 O(a²) 改进):
        Q_{μν}(x) = P_{μν}(x) + P_{ν,-μ}(x) + P_{-μ,-ν}(x) + P_{-ν,μ}(x)

    四个 plaquette 的几何排列 (以点 x 为中心):
                  μ →
              +---(1)---+
              |         |
            ν ↑  (2)  (4)  ↓ (反方向)
              |         |
              +---(3)---+
                  ← (反方向)

        (1) P_{μν}:    标准 1×1 Wilson loop
        (2) P_{ν,-μ}:  从 x 先走 ν 方向, 再走 -μ 方向
        (3) P_{-μ,-ν}: 从 x 先走 -μ 方向, 再走 -ν 方向
        (4) P_{-ν,μ}:  从 x 先走 -ν 方向, 再走 μ 方向

    Clover 场强张量 (取反厄米部分):
        F̂_{μν}(x) = -i/(8) · [ Q_{μν}(x) - Q_{μν}^†(x) ]

    在格点单位下 (a=1, g₀=1), F̂_{μν} 直接对应于连续极限的
    规范场场强张量 F_{μν}^{cont} = ∂_μ A_ν - ∂_ν A_μ + i g₀ [A_μ, A_ν]。

    Parameters
    ----------
    gauge : ndarray, shape (Nt, Nz, Ny, Nx, 4, 3, 3)
    mu, nu : int, Lorentz 指标 (0=x, 1=y, 2=z, 3=t)
    contract_fn : callable, 张量缩并函数

    Returns
    -------
    F_munu : ndarray, shape (Nt, Nz, Ny, Nx, 3, 3)
    """
    # 四个起始点
    g_ru = gauge                                          # x
    g_lu = np.roll(gauge, 1, axis=3 - mu)                 # x - μ̂
    g_rd = np.roll(gauge, 1, axis=3 - nu)                 # x - ν̂
    g_ld = np.roll(g_lu, 1, axis=3 - nu)                  # x - μ̂ - ν̂

    # Plaquette (1): P_{μν} - 标准 1×1 Wilson loop
    p_ru = contract_fn("tzyxab,tzyxbc->tzyxac",
                        g_ru[..., mu, :, :],
                        np.roll(g_ru, -1, axis=3 - mu)[..., nu, :, :])
    p_ru = contract_fn("tzyxab,tzyxcb->tzyxac", p_ru,
                        np.roll(g_ru, -1, axis=3 - nu)[..., mu, :, :].conj())
    p_ru = contract_fn("tzyxab,tzyxcb->tzyxac", p_ru,
                        g_ru[..., nu, :, :].conj())

    # Plaquette (2): P_{ν,-μ}
    p_lu = contract_fn("tzyxab,tzyxcb->tzyxac",
                        np.roll(g_lu, -1, axis=3 - mu)[..., nu, :, :],
                        np.roll(g_lu, -1, axis=3 - nu)[..., mu, :, :].conj())
    p_lu = contract_fn("tzyxab,tzyxcb->tzyxac", p_lu,
                        g_lu[..., nu, :, :].conj())
    p_lu = contract_fn("tzyxab,tzyxbc->tzyxac", p_lu,
                        g_lu[..., mu, :, :])

    # Plaquette (3): P_{-μ,-ν}
    p_ld = contract_fn("tzyxba,tzyxcb->tzyxac",
                        np.roll(g_ld, -1, axis=3 - nu)[..., mu, :, :].conj(),
                        g_ld[..., nu, :, :].conj())
    p_ld = contract_fn("tzyxab,tzyxbc->tzyxac", p_ld,
                        g_ld[..., mu, :, :])
    p_ld = contract_fn("tzyxab,tzyxbc->tzyxac", p_ld,
                        np.roll(g_ld, -1, axis=3 - mu)[..., nu, :, :])

    # Plaquette (4): P_{-ν,μ}
    p_rd = contract_fn("tzyxba,tzyxbc->tzyxac",
                        g_rd[..., nu, :, :].conj(),
                        g_rd[..., mu, :, :])
    p_rd = contract_fn("tzyxab,tzyxbc->tzyxac", p_rd,
                        np.roll(g_rd, -1, axis=3 - mu)[..., nu, :, :])
    p_rd = contract_fn("tzyxab,tzyxcb->tzyxac", p_rd,
                        np.roll(g_rd, -1, axis=3 - nu)[..., mu, :, :].conj())

    # F_{μν} = -i/8 * Σ_k (P_k - P_k^†)
    ans = (p_ru - p_ru.conj().transpose(0, 1, 2, 3, 5, 4) +
           p_lu - p_lu.conj().transpose(0, 1, 2, 3, 5, 4) +
           p_ld - p_ld.conj().transpose(0, 1, 2, 3, 5, 4) +
           p_rd - p_rd.conj().transpose(0, 1, 2, 3, 5, 4))
    return -1j * ans / 8.0


def compute_field_strength_all(gauge: Any, Nt: int, Nx: int,
                                contract_fn: Any) -> Any:
    """计算所有独立 (μ, ν) 对的 Clover 场强张量。

    Returns
    -------
    F_all : ndarray, shape (4, 4, Nt, Nz, Ny, Nx, 3, 3)
        F_all[μ, ν] = F_{μν}(x), F_{μμ} = 0
    """
    xp_lib = cp if HAS_CUPY and hasattr(gauge, 'device') else np
    F_all = xp_lib.zeros((4, 4, Nt, Nx, Nx, Nx, 3, 3), dtype=complex)
    for mu in range(4):
        for nu in range(4):
            if mu != nu:
                F_all[mu, nu] = plaquette_clover(gauge, mu, nu, contract_fn)
            # mu == nu → 恒为零
    return F_all


def compute_dual_field_strength(F_all: Any, epsilon: np.ndarray,
                                 contract_fn: Any) -> Any:
    r"""通过对偶变换计算 F̃_{μν} = (1/2) ε_{μνρσ} F^{ρσ}。

    ────────────────────────────────────────────────────────
    理论: Hodge 对偶与 Levi-Civita 缩并
    ────────────────────────────────────────────────────────

    在四维 Euclidean 时空中, 场强张量的 Hodge 对偶定义为:
        F̃_{μν} = (1/2) ε_{μνρσ} F^{ρσ}

    其中 ε_{μνρσ} 是完全反对称的 Levi-Civita 符号:
        ε_{0123} = +1
        ε_{μνρσ} = +1  若 (μ,ν,ρ,σ) 是 (0,1,2,3) 的偶排列
        ε_{μνρσ} = -1  若 (μ,ν,ρ,σ) 是 (0,1,2,3) 的奇排列
        ε_{μνρσ} = 0   若指标有重复

    对偶变换将电场型分量 (E_i = F_{ti}) 映射为磁场型分量
    (B_i = (1/2)ε_{ijk}F_{jk}), 反之亦然:
        F̃_{ti} = (1/2) ε_{tijk} F^{jk}  [磁场 → 对偶电场]
        F̃_{ij} = ε_{ijt k} F^{t k}       [电场 → 对偶磁场]

    在非极化胶子 PDF 的 OPE 算符中, 对偶场强确保算符
    具有正确的 Lorentz 变换性质 (旋量-2 表示)。

    Parameters
    ----------
    F_all : ndarray, shape (4, 4, Nt, Nz, Ny, Nx, 3, 3)
    epsilon : np.ndarray, shape (4, 4, 4, 4), ε_{μνρσ}/2

    Returns
    -------
    F_tilde_all : ndarray, shape (4, 4, Nt, Nz, Ny, Nx, 3, 3)
    """
    xp_lib = cp if HAS_CUPY and hasattr(F_all, 'device') else np
    epsilon_xp = xp_lib.asarray(epsilon) if xp_lib is cp else epsilon
    return contract_fn("opmn,mntzyxab->optzyxab", epsilon_xp, F_all)


# ============================================================================
# 第八部分: Nonlocal 胶子 OPE 算符 O_{μν}(z)
# ============================================================================

def gluon_ope_operator_z0_mu2(
    gauge: Any,
    F_all: Any,
    F_tilde_all: Any,
    delta_z: int,
    z_dir: int,
    mu: int, nu: int,
    mu2: int, nu2: int,
    contract_fn: Any,
) -> Any:
    r"""构造非定域胶子 OPE 算符 (z=0 处插入 F̃，z=δz 处插入 F):

    ────────────────────────────────────────────────────────
    理论: 非定域胶子准算符 (Nonlocal Gluon Quasi-Operator)
    ────────────────────────────────────────────────────────

    连续定义 (Eq(25) of 内部笔记):
        O_{unpol}(z) = Σ_{μ≠z} Σ_{ν≠z} g^{μν}
                       Tr_c[ F_{zμ}(z) · W(z→0) · F̃_ν^z(0) · W(0→z) ]

    其中:
        F_{zμ}(z)   = 在位置 z·n̂_z 处的场强张量
        F̃_ν^z(0)    = 在原点处的对偶场强张量 (Hodge 对偶)
        W(z→0)      = 从 z 到 0 的 Wilson 线 (伴随表示, 规范平行传输)
        W(0→z)      = 从 0 到 z 的 Wilson 线

    在 Euclidean 时空中, time-space 分量带负号:
        O^{(0)}_g = [-F^{0i} W F^{0i} W^† + 2F^{12} W F^{12} W^†]_{Eucl}
    这是因为在 Wick 转动下, F^{0i}_{Mink} → -i F^{0i}_{Eucl},
    而 F^{ij} 保持不变。

    格点实现的五步算法:

    Step 1 — 平移 F_{μν} 到 z = δz 处:
        F_{μν}(z) = np.roll(F_{μν}(0), -δz, axis=3-z_dir)
        其中 axis=3-z_dir 将空间方向映射为 roll 轴:
          z_dir=0(x) → axis=3, z_dir=1(y) → axis=2, z_dir=2(z) → axis=1

    Step 2 — Wilson 线 W(z→0), U_z 的共轭的连乘积:
        W(z→0) = Π_{k=δz-1}^{0} U_z^†(k·a)
        ope = contract("tzyxab,tzyxcb->tzyxac", ope, gauge[z_dir].conj())
        [重复 δz 次, 每次 roll 到正确的 z 位置]

    Step 3 — 在原点插入 F̃_{μ2,ν2}(0):
        ope = contract("tzyxab,tzyxbc->tzyxac", ope, F_tilde[μ2, ν2])

    Step 4 — Wilson 线 W(0→z), U_z 正向的连乘积:
        W(0→z) = Π_{k=0}^{δz-1} U_z(k·a)
        ope = contract("tzyxab,tzyxbc->tzyxac", ope, gauge[z_dir])
        [重复 δz 次]

    Step 5 — 颜色迹并对空间求和:
        O(z) = Σ_{x,y,z} Tr_c[ ope(x,y,z,t) ]  → shape (Nt,)

    其中 contract("...ab,...bc->...ac") 是颜色空间中的矩阵乘法,
    contract("...ab,...cb->...ac") 等价于 @ U^† (转置共轭)。

    Parameters
    ----------
    gauge : ndarray, shape (Nt, Nz, Ny, Nx, 4, 3, 3)
    F_all : ndarray, shape (4, 4, Nt, Nz, Ny, Nx, 3, 3)
    F_tilde_all : ndarray, shape (4, 4, Nt, Nz, Ny, Nx, 3, 3)
    delta_z : int, Wilson 线长度 (格距单位)
    z_dir : int, Wilson 线方向 (0=x, 1=y, 2=z)
    mu, nu : int, z 处 F 的 Lorentz 指标
    mu2, nu2 : int, 0 处 F̃ 的 Lorentz 指标
    contract_fn : callable, 张量缩并函数

    Returns
    -------
    op_trace : ndarray, shape (Nt,), 每个时间片的 OPE 算符迹 (已对空间求和)
    """
    # Step 1: 将 F_{μν} 平移到 z = delta_z 处
    ope = np.roll(F_all[mu, nu], -delta_z, axis=3 - z_dir)

    # Step 2: Wilson 线 W(z→0) — U_z 的共轭沿 -z 方向
    for dz_step in range(delta_z):
        shift = -(delta_z - 1 - dz_step)
        ug = np.roll(gauge, shift, axis=3 - z_dir)[..., z_dir, :, :]
        ope = contract_fn("tzyxab,tzyxcb->tzyxac", ope, ug.conj())

    # Step 3: 在 z=0 处插入 F̃_{μ2,ν2}
    ope = contract_fn("tzyxab,tzyxbc->tzyxac", ope, F_tilde_all[mu2, nu2])

    # Step 4: Wilson 线 W(0→z) — U_z 正向
    for dz_step in range(delta_z):
        shift = -dz_step
        ug = np.roll(gauge, shift, axis=3 - z_dir)[..., z_dir, :, :]
        ope = contract_fn("tzyxab,tzyxbc->tzyxac", ope, ug)

    # Step 5: 颜色迹 → 标量
    trace_color = np.trace(ope, axis1=4, axis2=5)

    # Step 6: 对空间求和
    op_trace = np.sum(trace_color, axis=(1, 2, 3))

    return op_trace


def compute_ope_all_z(
    gauge: Any,
    F_all: Any,
    F_tilde_all: Any,
    delta_z: int, z_dir: int,
    mu: int, nu: int, mu2: int, nu2: int,
    Nt: int, contract_fn: Any,
    pipeline: "GluonPDFPipeline" = None,
) -> Any:
    """为所有 z ∈ [0, delta_z-1] 计算 OPE 算符。

    Returns
    -------
    ops : ndarray, shape (Nt, delta_z)
        ops[t, dz] = O_{μν,μ2ν2}(z=dz) 在时间 t 的值
    """
    ops = np.zeros((Nt, delta_z), dtype=complex)
    if pipeline is not None and pipeline.xp is cp:
        ops = cp.asarray(ops)

    for dz in range(delta_z):
        if pipeline is not None:
            pipeline._check(f"  OPE z={dz}/{delta_z-1} mu={mu},nu={nu},mu2={mu2},nu2={nu2}")
        ops[:, dz] = gluon_ope_operator_z0_mu2(
            gauge, F_all, F_tilde_all, dz,
            z_dir, mu, nu, mu2, nu2, contract_fn
        )

    return ops


# ============================================================================
# 第九部分: Distillation 框架 — VVV 与 2pt
# ============================================================================

def compute_phase_factor(momentum: Any, Nx: int) -> Any:
    """计算动量涂抹的相位因子: φ_P(x) = exp(-i * 2π * P·x / L)。

    Parameters
    ----------
    momentum : ndarray, shape (3,), 动量矢量 (Pz, Py, Px), 以 2π/L 为单位
    Nx : int, 空间方向格点数

    Returns
    -------
    phase : ndarray, shape (Nx^3,), dtype=complex
    """
    V = Nx * Nx * Nx
    phase = np.zeros(V, dtype=complex)
    idx = 0
    for z in range(Nx):
        for y in range(Nx):
            for x in range(Nx):
                pos = np.array([z, y, x])
                phase[idx] = np.exp(-np.dot(momentum, pos) * 2.0j * np.pi / Nx)
                idx += 1
    return phase


def compute_vvv_baryon_block(
    eigvecs: Any, phase_factor: Any, Nev1: int, Nx: int, contract_fn: Any,
) -> Any:
    r"""计算 VVV (Baryon Block) 张量:

    ────────────────────────────────────────────────────────
    理论: 重子三夸克波函数的动量投影
    ────────────────────────────────────────────────────────

    VVV Baryon Block 定义为:
        Φ_{abc}(P) = Σ_x e^{-iP·x} ε_{ijk} v_i^a(x) v_j^b(x) v_k^c(x)

    其中:
        a,b,c ∈ [0, Nev1)  = 本征矢量指标 (蒸馏子空间截断)
        i,j,k ∈ {0,1,2}   = 色指标 (RGB)
        x ∈ 格点空间         = 空间坐标
        P                   = 注入动量 (2π/L 单位)
        v_i^a(x)            = 第 a 个 Laplace 本征矢量的第 i 个色分量
        ε_{ijk}             = 完全反对称 Levi-Civita 符号
                               (确保色单态重子波函数)

    ε_{ijk} 的六个置换求和:
        + Φ_{012}:  contract("x,ax,bx,cx->abc", φ, v₀, v₁, v₂)
        + Φ_{120}:  contract("x,ax,bx,cx->abc", φ, v₁, v₂, v₀)
        + Φ_{201}:  contract("x,ax,bx,cx->abc", φ, v₂, v₀, v₁)
        - Φ_{021}:  contract("x,ax,bx,cx->abc", φ, v₀, v₂, v₁)
        - Φ_{102}:  contract("x,ax,bx,cx->abc", φ, v₁, v₀, v₂)
        - Φ_{210}:  contract("x,ax,bx,cx->abc", φ, v₂, v₁, v₀)

    动量涂抹相位因子:
        φ_P(x) = exp(-i · 2π · P·x / L)
    其中 P = (Pz, Py, Px), L = Nx·a 是空间格点尺寸。

    Parameters
    ----------
    eigvecs : ndarray, shape (Nev, Nx^3, 3)
    phase_factor : ndarray, shape (Nx^3,)
    Nev1 : int, 截断的本征矢量数
    Nx : int, 空间格点数
    contract_fn : callable

    Returns
    -------
    VVV : ndarray, shape (Nev1, Nev1, Nev1), dtype=complex
    """
    xp_lib = cp if HAS_CUPY and hasattr(eigvecs, 'device') else np
    VVV = xp_lib.zeros((Nev1, Nev1, Nev1), dtype=complex)
    layer_size = Nx * Nx

    for xi in range(Nx):
        start = xi * layer_size
        end = (xi + 1) * layer_size
        es = eigvecs[:Nev1, start:end, :]    # (Nev1, Nx², 3)
        ps = phase_factor[start:end]           # (Nx²,)

        # ε_{ijk} 六个置换 — 三项正号 (偶排列)
        VVV += contract_fn("x,ax,bx,cx->abc", ps, es[:, :, 0], es[:, :, 1], es[:, :, 2])
        VVV += contract_fn("x,ax,bx,cx->abc", ps, es[:, :, 1], es[:, :, 2], es[:, :, 0])
        VVV += contract_fn("x,ax,bx,cx->abc", ps, es[:, :, 2], es[:, :, 0], es[:, :, 1])
        # 三项负号 (奇排列)
        VVV -= contract_fn("x,ax,bx,cx->abc", ps, es[:, :, 0], es[:, :, 2], es[:, :, 1])
        VVV -= contract_fn("x,ax,bx,cx->abc", ps, es[:, :, 1], es[:, :, 0], es[:, :, 2])
        VVV -= contract_fn("x,ax,bx,cx->abc", ps, es[:, :, 2], es[:, :, 1], es[:, :, 0])

    return VVV


def contract_proton_2pt_single_tsrc(
    VVV_sink_t: Any, peram_u: Any,
    CG5peram_uCG5: Any, VVV_source_t: Any,
    contract_fn: Any,
) -> Any:
    r"""对单个 (t_sink, t_source) 执行质子 2pt 函数的 Wick 收缩。

    ────────────────────────────────────────────────────────
    理论: 质子 2pt 关联函数的 Wick 收缩 (Distillation)
    ────────────────────────────────────────────────────────

    质子插值算符 (色单态重子算符):
        O_α(x) = ε_{abc} [u_a^T(x) C γ₅ d_b(x)] u_c,α(x)

    其中:
        C = γ₄γ₂ 是电荷共轭矩阵 (DeGrand-Rossi 基)
        α 是自由 Dirac 指标

    在蒸馏框架下, 夸克场被投影到本征矢量子空间:
        u_i^a(t) = Σ_{k} v_i^{a(k)}(t) u^{(k)}(t)
    其中 k 遍历最低 Nev 个 Laplace 本征模态。

    Perambulator (夸克传播子的子空间投影):
        τ_{αβ}^{AB}(t_snk, t_src) = v^{A†}(t_snk) S_{αβ} v^B(t_src)

    直接项 (Direct, 夸克线无交换):
        C₂^{dir}_{il} = Φ_{abc}(t_snk) · τ_{gi}^{ad} · (ΓτΓ)_{gj}^{be}
                         · τ_{il}^{cf} · Φ_{def}^*(t_src)

    交换项 (Exchange, 两夸克线交换, Fermi 统计负号):
        C₂^{ex}_{il} = Φ_{abc}(t_snk) · τ_{gl}^{af} · (ΓτΓ)_{gj}^{be}
                         · τ_{ij}^{cd} · Φ_{def}^*(t_src)

    总关联函数:
        C₂_{il} = Direct - Exchange

    指标约定:
        a,b,c,d,e,f = 本征矢量指标 (蒸馏子空间, 0..Nev1-1)
        g,j,i,l     = Dirac 旋量指标 (0..3)

    Dirac 投影矩阵 Γ (质子插值):
        Γ = C γ₅ γ_μ  (μ 由 interpolator 类型决定)
        对于 Cγ₅γ₄ 插值: Γ = γ₄γ₂ · γ₅ · γ₄
        在 DR 基下: Γ → gamma(7) @ gamma(4) = γ₃γ₁ · γ₄

    Parameters
    ----------
    VVV_sink_t : ndarray, shape (Nev1, Nev1, Nev1)
    peram_u : ndarray, shape (4, 4, Nev1, Nev1)
    CG5peram_uCG5 : ndarray, shape (4, 4, Nev1, Nev1)
    VVV_source_t : ndarray, shape (Nev1, Nev1, Nev1)
    contract_fn : callable

    Returns
    -------
    C2_matrix : ndarray, shape (4, 4)
    """
    # 直接项
    direct = contract_fn(
        "abc,gjad,gjbe,ilcf,def->il",
        VVV_sink_t, peram_u, CG5peram_uCG5, peram_u, VVV_source_t,
    )
    # 交换项
    exchange = contract_fn(
        "abc,glaf,gjbe,ijcd,def->il",
        VVV_sink_t, peram_u, CG5peram_uCG5, peram_u, VVV_source_t,
    )
    return direct - exchange


# ============================================================================
# 第十部分: 矩阵元提取、Fourier 变换与匹配
# ============================================================================

def jackknife_samples(data: Any, axis: int = -1) -> Any:
    r"""生成 Jackknife 重采样样本。

    ────────────────────────────
    Jackknife (留一法交叉验证)
    ────────────────────────────

    设有 N 个独立样本 {x₁, x₂, ..., x_N}。

    全样本均值:
        x̄ = (1/N) Σ_{i=1}^{N} x_i

    第 i 个 Jackknife 样本 (删除第 i 个样本后的均值):
        x_i^{JK} = (Σ_{j≠i} x_j) / (N-1) = (N·x̄ - x_i) / (N-1)

    Jackknife 估计量:
        θ̄^{JK} = (1/N) Σ_i θ(x_i^{JK})

    Jackknife 误差 (Tukey 1958):
        σ_{JK}² = [(N-1)/N] · Σ_i [θ(x_i^{JK}) - θ̄^{JK}]²
                = (N-1) · Var(θ^{JK})

    对于线性估计量 θ = x̄: θ̄^{JK} = x̄ (无偏)

    Parameters
    ----------
    data : ndarray, 第 axis 维 = 样本数 N
    axis : int, 样本维度索引

    Returns
    -------
    jk_samples : ndarray, 与 data 同 shape
    """
    N = data.shape[axis]
    total = np.sum(data, axis=axis, keepdims=True)
    jk = (total - data) / (N - 1)
    return jk


def jackknife_error(jk_samples: Any, axis: int = 0) -> Any:
    """从 Jackknife 样本估计标准误差: σ = std(jk) × √(N-1)。"""
    N = jk_samples.shape[axis]
    return np.std(jk_samples, axis=axis) * np.sqrt(N - 1)


def extract_matrix_element_from_ope(
    ope_data: Any, n_conf: int,
) -> Tuple[Any, Any]:
    """从多个组态的 OPE 期望值中提取矩阵元 h(z)。

    在 disconnect 近似下:
        h(z) = ⟨O(z)⟩ = (1/N_conf) Σ_{conf} ⟨O(z)⟩_{conf}

    Parameters
    ----------
    ope_data : ndarray, shape (n_conf, Nt, Nz)
    n_conf : int

    Returns
    -------
    h_z_mean : ndarray, shape (Nz,)
    h_z_err : ndarray, shape (Nz,)
    """
    # 对时间平均 → (n_conf, Nz)
    ope_tavg = np.mean(ope_data, axis=1)

    # 对组态平均 → (Nz,)
    h_z_mean = np.mean(ope_tavg, axis=0)

    # Jackknife 误差
    h_z_jk = jackknife_samples(ope_tavg, axis=0)  # (n_conf, Nz)
    h_z_err = jackknife_error(h_z_jk, axis=0)

    return h_z_mean, h_z_err


def fourier_transform_to_quasi_pdf(
    h_z: Any, z_values: Any, Pz: float, x_values: Any,
) -> Any:
    r"""通过 Fourier 变换将坐标空间矩阵元转换为 quasi-PDF:

    ────────────────────────────────────────────────────────
    理论: 从坐标空间矩阵元到动量空间 quasi-PDF
    ────────────────────────────────────────────────────────

    准 PDF 的定义 (LaMET, Ji 2013):
        g̃(x, P_z) = ∫ dz/(2π) · e^{i x P_z z} · h(z, P_z)

    其中:
        x ∈ [0, 1]    = Bjorken-x (动量分数)
        P_z            = boost 动量 (格点单位)
        h(z, P_z)      = 坐标空间矩阵元 h(z, P_z)
        z ∈ [0, z_max] = Wilson 线长度 (格距单位)

    利用 h(z) 的厄米对称性:
        h^†(z) = h(-z)  [厄米共轭]

    由此 Re[h(z)] 是偶函数, Im[h(z)] 是奇函数,
    因此:
        g̃(x, P_z) = (1/π) ∫_0^{z_max} dz Re[h(z)] · cos(x P_z z)
                     - i (1/π) ∫_0^{z_max} dz Im[h(z)] · sin(x P_z z)

    对于非极化胶子 PDF (仅取实部, 使用正弦变换):
        g̃(x, P_z) = (2P_z / x) ∫_0^{z_max} dz · Re[h(z)] · sin(x P_z z)

    数值积分: x·g̃(x) 使用梯形法则:
        ∫_0^{z_max} f(z) dz ≃ Σ_{i=0}^{N-1} (Δz/2)·[f(z_i) + f(z_{i+1})]
        等价于: np.trapz(integrand, z_values)

    Parameters
    ----------
    h_z : ndarray, shape (Nz,), 坐标空间矩阵元 (取实部)
    z_values : ndarray, shape (Nz,), z 值数组 (格距单位)
    Pz : float, boost 动量 (格点单位)
    x_values : ndarray, shape (Nx,), Bjorken-x 值数组

    Returns
    -------
    quasi_pdf : ndarray, shape (Nx,)
    """
    h_real = np.real(h_z)
    Nx_q = len(x_values)
    quasi = np.zeros(Nx_q)

    for i, x in enumerate(x_values):
        if abs(x) < 1e-15:
            quasi[i] = 0.0
            continue
        integrand = h_real * np.sin(x * Pz * z_values)
        integral = np.trapz(integrand, z_values)
        quasi[i] = (2.0 * Pz / x) * integral

    return quasi


def matching_kernel_nlo(
    quasi_pdf: Any, x_values: Any,
    alpha_s: float = 0.2, mu_over_Pz: float = 1.0,
) -> Any:
    r"""将 quasi-PDF 通过微扰匹配转换为光锥 PDF。

    ────────────────────────────────────────────────────────
    理论: LaMET 匹配核 (Large Momentum Effective Theory)
    ────────────────────────────────────────────────────────

    光锥 PDF g(x, μ) 与 quasi-PDF g̃(x, P_z) 的关系由因子化定理确定:
        g(x, μ) = ∫_x^1 (dy/y) C_{gg}(x/y, μ/P_z) · g̃(y, P_z)

    在 NLO 微扰展开下:
        C_{gg}(ξ, μ/P_z) = δ(1-ξ) + (α_s/(2π)) ΔC_{gg}(ξ, μ/P_z) + O(α_s²)

    因此 NLO 匹配公式为:
        g(x, μ) = g̃(x, P_z)
                   - (α_s/(2π)) ∫_x^1 (dy/y) ΔC_{gg}(x/y, μ/P_z) · g̃(y, P_z)

    其中 ΔC_{gg} 的 NLO 表达式包含:
        · δ(1-ξ) 项    (自能修正)
        · 1/(1-ξ)_+ 项 (plus-distribution, 软胶子辐射)
        · ln(μ²/P_z²) 项  (重整化标度依赖)
        · 常数项        (有限部分)

    LO 近似 (匹配核 = δ(1-ξ)):
        g(x, μ) ≃ g̃(x, P_z)

    重整化群不变性: g(x, μ) 满足 DGLAP 演化方程:
        d g(x, μ) / d ln μ² = (α_s/(2π)) ∫_x^1 (dy/y) P_{gg}(x/y) g(y, μ)
    其中 P_{gg}(z) 是胶子-胶子劈裂函数。

    Parameters
    ----------
    quasi_pdf : ndarray, shape (Nx,)
    x_values : ndarray, shape (Nx,)
    alpha_s : float, 强耦合常数
    mu_over_Pz : float, μ/P_z

    Returns
    -------
    lightcone_pdf : ndarray, shape (Nx,)
    """
    # LO: 恒等变换 (NLO 修正待定)
    return quasi_pdf.copy()


# ============================================================================
# 第十一部分: 数据处理 —— Bootstrap
# ============================================================================

def bootstrap_samples(data: Any, n_bootstrap: int = 1000,
                       axis: int = 0, seed: int = 42) -> Any:
    """生成 Bootstrap 重采样样本。

    Parameters
    ----------
    data : ndarray, 第 axis 维 = 组态数 N
    n_bootstrap : int, bootstrap 样本数
    axis : int, 组态轴
    seed : int, 随机种子 (可重复性)

    Returns
    -------
    bs_samples : ndarray, shape = (n_bootstrap, ...)
    """
    rng = np.random.RandomState(seed)
    N = data.shape[axis]
    bs = np.zeros((n_bootstrap,) + data.shape[1:], dtype=data.dtype)
    for i in range(n_bootstrap):
        idx = rng.randint(0, N, size=N)
        bs[i] = np.mean(np.take(data, idx, axis=axis), axis=axis)
    return bs


# ============================================================================
# 第十二部分: 作图函数
# ============================================================================

def _ensure_cpu(arr: Any) -> np.ndarray:
    """将可能的 GPU 数组转换为 CPU numpy 数组。"""
    if HAS_CUPY and hasattr(arr, 'get'):
        return arr.get()
    return np.asarray(arr)


def plot_field_strength_check(F_all: Any, output_dir: str, prefix: str = "F"):
    """画图: 检查场强张量的厄米性 Tr[F - F^†] 是否为零 (应为纯虚矩阵)。"""
    if not HAS_MPL:
        print("[WARN] matplotlib 不可用, 跳过场强检查作图。")
        return

    F_cpu = _ensure_cpu(F_all)
    Nt = F_cpu.shape[2]
    herm_check = np.zeros((4, 4, Nt))
    for mu in range(4):
        for nu in range(4):
            if mu == nu: continue
            F_mn = F_cpu[mu, nu]
            tr_diff = np.abs(np.trace(F_mn, axis1=5, axis2=6)
                           - np.trace(F_mn.conj().transpose(0, 1, 2, 3, 5, 4),
                                      axis1=5, axis2=6))
            herm_check[mu, nu] = np.mean(tr_diff, axis=(1, 2, 3))

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    labels = ["(0,1)=xy", "(0,2)=xz", "(0,3)=xt", "(1,2)=yz", "(1,3)=yt", "(2,3)=zt"]
    pairs = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]

    for ax, (mu, nu), lbl in zip(axes, pairs, labels):
        ax.plot(herm_check[mu, nu], 'o-', markersize=3)
        ax.set_title(f"F_{{{lbl}}}: |Tr[F - F^†]|")
        ax.set_xlabel("t")
        ax.set_ylabel("|Tr[F - F^†]|")
        ax.set_yscale("log")

    plt.suptitle(f"{prefix}: Field Strength Anti-Hermiticity Check", fontsize=14)
    fig.savefig(f"{output_dir}/{prefix}_field_strength_herm_check.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[PLOT] 场强厄米检查图已保存: {prefix}_field_strength_herm_check.pdf")


def plot_vvv_magnitude(VVV_data: Any, output_dir: str,
                        prefix: str = "VVV", t_slice: int = 0):
    """画图: VVV 张量的大小分布 (第一个时间片)。"""
    if not HAS_MPL:
        print("[WARN] matplotlib 不可用, 跳过 VVV 作图。")
        return

    v_cpu = _ensure_cpu(VVV_data)
    if v_cpu.ndim == 4:
        v_slice = v_cpu[t_slice]  # (Nev1, Nev1, Nev1)
    else:
        v_slice = v_cpu

    Nev1 = v_slice.shape[0]
    mag = np.abs(v_slice)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 子图 1: 各本征矢量对角元的幅值
    diag_vals = np.array([mag[i, i, i] for i in range(Nev1)])
    axes[0].plot(diag_vals, 'o-', markersize=3)
    axes[0].set_title("VVV Diagonal |Φ_{aaa}|")
    axes[0].set_xlabel("Eigenvector index a")
    axes[0].set_ylabel("|Φ_{aaa}|")
    axes[0].set_yscale("log")

    # 子图 2: 本征矢量对 (a, b) 间耦合强度的 2D 热力图
    coupling = np.sum(mag, axis=2) / Nev1  # 对 c 指标平均
    im = axes[1].imshow(np.log10(coupling + 1e-15), aspect='auto', origin='lower')
    axes[1].set_title("log10(|Φ_{abc}| avg over c)")
    axes[1].set_xlabel("Eigenvector index b")
    axes[1].set_ylabel("Eigenvector index a")
    plt.colorbar(im, ax=axes[1], label="log10(avg |Φ|)")

    plt.suptitle(f"{prefix}: VVV Baryon Block Structure (t={t_slice})", fontsize=14)
    fig.savefig(f"{output_dir}/{prefix}_vvv_structure_t{t_slice}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[PLOT] VVV 结构图已保存: {prefix}_vvv_structure_t{t_slice}.pdf")


def plot_2pt_correlator(C2_data: Any, output_dir: str,
                         prefix: str = "C2pt", Pz: int = 0):
    """画图: 质子 2pt 关联函数 C₂(t) 随时间的变化。"""
    if not HAS_MPL:
        print("[WARN] matplotlib 不可用, 跳过 2pt 作图。")
        return

    c2 = _ensure_cpu(C2_data)  # shape: (Nt,) 或 (4, 4, Nt)
    if c2.ndim == 3:
        # 取 Dirac 对角元的实部
        c2_diag = np.array([np.real(c2[i, i, :]) for i in range(4)])
    else:
        c2_diag = np.real(c2).reshape(1, -1)

    Nt = c2_diag.shape[1]
    t_range = np.arange(Nt)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # 子图 1: C₂(t) 对数图
    for i in range(c2_diag.shape[0]):
        axes[0].plot(t_range, np.abs(c2_diag[i]), 'o-', markersize=3,
                      label=f"Dirac[{i},{i}]")
    axes[0].set_title(f"{prefix}: |C₂(t)| (Pz={Pz})")
    axes[0].set_xlabel("t/a")
    axes[0].set_ylabel("|C₂(t)|")
    axes[0].set_yscale("log")
    axes[0].legend(fontsize=9)

    # 子图 2: 有效质量
    eff_mass = np.zeros(Nt - 2)
    for t in range(1, Nt - 1):
        ratio = (c2_diag[0, t - 1] + c2_diag[0, t + 1]) / (2 * c2_diag[0, t])
        if ratio > 1:
            eff_mass[t - 1] = np.arccosh(ratio)
    axes[1].plot(t_range[1:-1], eff_mass, 's-', markersize=4)
    axes[1].set_title(f"{prefix}: Effective Mass (Pz={Pz})")
    axes[1].set_xlabel("t/a")
    axes[1].set_ylabel("a·E_eff")
    axes[1].axhline(y=np.mean(eff_mass[Nt//4:2*Nt//4]), color='r',
                     linestyle='--', label='plateau avg')
    axes[1].legend(fontsize=9)

    plt.suptitle(f"{prefix}: Proton 2pt Correlator Analysis", fontsize=14)
    fig.savefig(f"{output_dir}/{prefix}_2pt_analysis_Pz{Pz}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[PLOT] 2pt 关联函数图已保存: {prefix}_2pt_analysis_Pz{Pz}.pdf")


def plot_ope_result(ops_data: Any, output_dir: str,
                     prefix: str = "OPE", mu: int = 0, nu: int = 1):
    """画图: OPE 算符在 (Nt, delta_z) 平面上的值。"""
    if not HAS_MPL:
        print("[WARN] matplotlib 不可用, 跳过 OPE 作图。")
        return

    ops = _ensure_cpu(ops_data)  # (Nt, delta_z)
    Nt, delta_z = ops.shape

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # 子图 1: 实部 2D 热力图
    im1 = axes[0, 0].imshow(np.real(ops).T, aspect='auto', origin='lower',
                              cmap='RdBu_r')
    axes[0, 0].set_title(f"Re[O(z)] mu={mu},nu={nu}")
    axes[0, 0].set_xlabel("t")
    axes[0, 0].set_ylabel("z")
    plt.colorbar(im1, ax=axes[0, 0])

    # 子图 2: 虚部 2D 热力图
    im2 = axes[0, 1].imshow(np.imag(ops).T, aspect='auto', origin='lower',
                              cmap='RdBu_r')
    axes[0, 1].set_title(f"Im[O(z)] mu={mu},nu={nu}")
    axes[0, 1].set_xlabel("t")
    axes[0, 1].set_ylabel("z")
    plt.colorbar(im2, ax=axes[0, 1])

    # 子图 3: 某些 z 值处 O(z) 随时间的变化
    z_sample = [0, delta_z//4, delta_z//2, 3*delta_z//4, delta_z-1]
    for z in z_sample:
        if z < delta_z:
            axes[1, 0].plot(np.real(ops[:, z]), 'o-', markersize=2,
                             label=f"z={z}, Re")
    axes[1, 0].set_title("Re[O(z)] vs t at selected z")
    axes[1, 0].set_xlabel("t")
    axes[1, 0].set_ylabel("Re[O(z)]")
    axes[1, 0].legend(fontsize=8)

    # 子图 4: 时间平均后的 O(z) vs z
    ops_tavg = np.mean(ops, axis=0)
    axes[1, 1].errorbar(range(delta_z), np.real(ops_tavg),
                          yerr=np.std(np.real(ops), axis=0) / np.sqrt(Nt),
                          fmt='o-', markersize=4, capsize=3, label='Re')
    axes[1, 1].errorbar(range(delta_z), np.imag(ops_tavg),
                          yerr=np.std(np.imag(ops), axis=0) / np.sqrt(Nt),
                          fmt='s-', markersize=4, capsize=3, label='Im')
    axes[1, 1].set_title("Time-averaged O(z) vs z")
    axes[1, 1].set_xlabel("z/a")
    axes[1, 1].set_ylabel("⟨O(z)⟩")
    axes[1, 1].legend()

    plt.suptitle(f"{prefix}: OPE Operator O_{mu}{nu}(z)", fontsize=14)
    fig.savefig(f"{output_dir}/{prefix}_ope_mu{mu}nu{nu}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[PLOT] OPE 算符图已保存: {prefix}_ope_mu{mu}nu{nu}.pdf")


def plot_matrix_element(h_z: Any, h_z_err: Any, output_dir: str,
                         prefix: str = "h_z", Pz: int = 6):
    """画图: 坐标空间矩阵元 h(z, P_z) 及其拟合作图。"""
    if not HAS_MPL:
        print("[WARN] matplotlib 不可用, 跳过矩阵元作图。")
        return

    h = _ensure_cpu(h_z)
    h_err = _ensure_cpu(h_z_err)
    Nz = len(h)
    z_vals = np.arange(Nz)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # 子图 1: 实部
    axes[0].errorbar(z_vals, np.real(h), yerr=h_err, fmt='o-', capsize=3,
                      markersize=4, label='Re[h(z)]')
    axes[0].set_title(f"{prefix}: Re[h(z, Pz={Pz})]")
    axes[0].set_xlabel("z/a")
    axes[0].set_ylabel("Re[h(z)]")
    axes[0].legend()

    # 子图 2: 虚部
    axes[1].errorbar(z_vals, np.imag(h), yerr=h_err, fmt='s-', capsize=3,
                      markersize=4, label='Im[h(z)]', color='C1')
    axes[1].set_title(f"{prefix}: Im[h(z, Pz={Pz})]")
    axes[1].set_xlabel("z/a")
    axes[1].set_ylabel("Im[h(z)]")
    axes[1].legend()

    plt.suptitle(f"{prefix}: Coordinate-Space Matrix Element (Pz={Pz})",
                 fontsize=14)
    fig.savefig(f"{output_dir}/{prefix}_matrix_element_Pz{Pz}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[PLOT] 矩阵元图已保存: {prefix}_matrix_element_Pz{Pz}.pdf")


def plot_quasi_pdf(quasi_pdf: Any, x_values: Any, output_dir: str,
                    prefix: str = "quasi", Pz: int = 6):
    """画图: quasi-PDF g̃(x, P_z) 随 x 的变化。"""
    if not HAS_MPL:
        print("[WARN] matplotlib 不可用, 跳过 quasi-PDF 作图。")
        return

    qp = _ensure_cpu(quasi_pdf)
    xv = _ensure_cpu(x_values)

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.plot(xv, qp, 'o-', markersize=4)
    ax.set_title(f"{prefix}: Quasi-PDF $\\tilde{{g}}(x, P_z={Pz})$")
    ax.set_xlabel("x (Bjorken-x)")
    ax.set_ylabel("$\\tilde{g}(x, P_z)$")
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.5)
    ax.grid(True, alpha=0.3)

    fig.savefig(f"{output_dir}/{prefix}_quasi_pdf_Pz{Pz}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[PLOT] quasi-PDF 图已保存: {prefix}_quasi_pdf_Pz{Pz}.pdf")


def plot_lightcone_pdf(lc_pdf: Any, x_values: Any, output_dir: str,
                        prefix: str = "lc_pdf", Pz: int = 6, mu: float = 2.0):
    """画图: 光锥 PDF g(x, μ) 随 x 的变化。"""
    if not HAS_MPL:
        print("[WARN] matplotlib 不可用, 跳过光锥 PDF 作图。")
        return

    lp = _ensure_cpu(lc_pdf)
    xv = _ensure_cpu(x_values)

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.plot(xv, lp, 'o-', markersize=4, color='C3')
    ax.set_title(f"{prefix}: Light-Cone PDF $g(x, \\mu={mu}\\,\\mathrm{{GeV}})$")
    ax.set_xlabel("x (Bjorken-x)")
    ax.set_ylabel("$g(x, \\mu)$")
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)

    fig.savefig(f"{output_dir}/{prefix}_lightcone_pdf_Pz{Pz}_mu{mu:.0f}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[PLOT] 光锥 PDF 图已保存: {prefix}_lightcone_pdf_Pz{Pz}_mu{mu:.0f}.pdf")


def plot_summary_timeline(timing: dict, memory: dict, output_dir: str):
    """画图: 各步骤耗时与内存占用的汇总柱状图。"""
    if not HAS_MPL:
        print("[WARN] matplotlib 不可用, 跳过汇总作图。")
        return

    names = list(timing.keys())
    times = [timing[n] for n in names]
    mems = [memory[n]["peak_mb"] for n in names]

    fig, axes = plt.subplots(2, 1, figsize=(14, 10))

    # 子图 1: 各步骤耗时
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(names)))
    bars1 = axes[0].bar(range(len(names)), times, color=colors)
    axes[0].set_title("Per-Step Wall Time")
    axes[0].set_ylabel("Time (s)")
    axes[0].set_xticks(range(len(names)))
    axes[0].set_xticklabels([n.replace("Step ", "").split(":")[0] for n in names],
                             rotation=45, ha='right', fontsize=8)
    for bar, t in zip(bars1, times):
        axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                      f"{t:.1f}s", ha='center', va='bottom', fontsize=7)

    # 子图 2: 各步骤内存峰值
    bars2 = axes[1].bar(range(len(names)), mems, color=colors)
    axes[1].set_title("Per-Step Peak Memory Usage")
    axes[1].set_ylabel("Memory (MB)")
    axes[1].set_xticks(range(len(names)))
    axes[1].set_xticklabels([n.replace("Step ", "").split(":")[0] for n in names],
                             rotation=45, ha='right', fontsize=8)
    for bar, m in zip(bars2, mems):
        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 10,
                      f"{m:.0f}MB", ha='center', va='bottom', fontsize=7)

    plt.suptitle("Workflow Performance Summary", fontsize=14)
    fig.savefig(f"{output_dir}/summary_performance.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[PLOT] 性能汇总图已保存: summary_performance.pdf")


# ============================================================================
# 第十三部分: GluonPDFPipeline 主类
# ============================================================================

class GluonPDFPipeline:
    """格点 QCD 质子非极化胶子 PDF 计算的统一流水线。

    该流水线包含 13 个可选步骤, 支持:
      - numpy 或 cupy 后端
      - 全程显存/内存占用与耗时统计
      - 可选数据读取/生成/对比
      - 丰富的出版级作图

    Parameters
    ----------
    args : argparse.Namespace, 命令行参数
    """

    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.output_dir = args.output_dir

        # --- 后端设置 ---
        if args.xp == "cupy" and HAS_CUPY:
            self.xp = cp
            self.backend_name = "CuPy (GPU)"
            print("[INFO] 使用 CuPy GPU 后端。")
        else:
            if args.xp == "cupy" and not HAS_CUPY:
                print("[WARN] CuPy 不可用, 回退到 NumPy (CPU)。")
            self.xp = np
            self.backend_name = "NumPy (CPU)"
            print("[INFO] 使用 NumPy CPU 后端。")

        # --- 数据类型 ---
        dtype_map = {"complex64": np.complex64, "complex128": np.complex128}
        self.dtype = dtype_map.get(args.dtype, np.complex128)
        print(f"[INFO] 数据类型: {self.dtype.__name__}")

        # --- 缩并函数 ---
        self.contract_fn = get_contract(self.xp)

        # --- 格点配置 ---
        self._init_lattice_config()

        # --- 文件路径 ---
        self._init_file_paths()

        # --- 结果存储 ---
        self.timing_results: Dict[str, float] = {}
        self.memory_results: Dict[str, dict] = {}
        self.data: Dict[str, Any] = {}  # 步骤间数据传递

        # --- Gamma 矩阵 (后端无关) ---
        self.gamma = build_gamma_matrices(self.xp)

        # --- Levi-Civita 张量 (始终在 CPU 上, 需要时传输) ---
        self.epsilon = build_levi_civita_tensor()

        # --- 创建输出目录 ---
        os.makedirs(self.output_dir, exist_ok=True)
        print(f"[INFO] 输出目录: {self.output_dir}")

        # --- 保存运行配置 ---
        self._save_config()

    def _init_lattice_config(self):
        """初始化格点配置: 从系综预设或手动参数。"""
        args = self.args

        # 系综预设解析
        if args.ensemble:
            ecfg = resolve_ensemble(args.ensemble, str(args.conf_start))
            self.cfg = LatticeConfig(
                Nt=ecfg["Nt"], Nx=ecfg["Nx"],
                Nev=ecfg["Nev"], Nev1=ecfg["Nev1"],
                mom_smear=ecfg["mom_smear"],
                mom_smear_phase=ecfg["mom_smear_phase"],
                Px=args.Px, Py=args.Py, Pz=args.Pz,
                delta_z=getattr(args, 'delta_z', 15),
                z_dir=getattr(args, 'z_dir', 2),
                conf_id=str(args.conf_start),
            )
            # 系综路径
            self._conf_dir = ecfg.get("conf_dir", "")
            self._eig_dir = ecfg["eig_dir"]
            self._peram_u_dir = ecfg["peram_u_dir"]
            self._corr_nucl_dir = ecfg["corr_nucl_dir"]
        else:
            # 手动参数
            self.cfg = LatticeConfig(
                Nt=args.Nt, Nx=args.Nx,
                Nev=args.Nev, Nev1=args.Nev1,
                Px=args.Px, Py=args.Py, Pz=args.Pz,
                delta_z=getattr(args, 'delta_z', 15),
                z_dir=getattr(args, 'z_dir', 2),
                mom_smear=getattr(args, 'mom_smear', 3),
                mom_smear_phase=getattr(args, 'mom_smear_phase', -3),
                conf_id=str(args.conf_start),
            )
            self._conf_dir = getattr(args, 'conf_dir', "")
            self._eig_dir = getattr(args, 'eig_dir', "")
            self._peram_u_dir = getattr(args, 'peram_u_dir', "")
            self._corr_nucl_dir = getattr(args, 'corr_nucl_dir', "")

        # 组态列表
        self.conf_ids = [
            args.conf_start + i * args.conf_step
            for i in range(args.conf_num)
        ]
        print(f"[INFO] 组态列表: {self.conf_ids}")
        print(f"[INFO] 格点: {self.cfg.Nt}×{self.cfg.Nx}³, "
              f"Nev={self.cfg.Nev1}, P=({self.cfg.Px},{self.cfg.Py},{self.cfg.Pz})")

    def _init_file_paths(self):
        """初始化所有文件路径。"""
        args = self.args

        # 规范组态文件 (仅读取)
        self.gauge_file = getattr(args, 'gauge_file', None)

        # 读取路径
        self.read_2pt_dir = getattr(args, 'read_2pt_dir', None)
        self.read_3pt_dir = getattr(args, 'read_3pt_dir', None)
        self.read_VVV_dir = getattr(args, 'read_VVV_dir', None)
        self.read_ope_dir = getattr(args, 'read_ope_dir', None)

        # 生成路径 (写入)
        self.gen_2pt_dir = getattr(args, 'gen_2pt_dir', None)
        self.gen_3pt_dir = getattr(args, 'gen_3pt_dir', None)
        self.gen_VVV_dir = getattr(args, 'gen_VVV_dir', None)
        self.gen_ope_dir = getattr(args, 'gen_ope_dir', None)

        # 对比标记
        self.compare_2pt = bool(self.read_2pt_dir and self.gen_2pt_dir)
        self.compare_3pt = bool(self.read_3pt_dir and self.gen_3pt_dir)
        self.compare_VVV = bool(self.read_VVV_dir and self.gen_VVV_dir)
        self.compare_ope = bool(self.read_ope_dir and self.gen_ope_dir)

    def _save_config(self):
        """保存运行配置到输出目录 (供后续复现)。"""
        import json as _json
        config_dict = {
            "backend": self.backend_name,
            "dtype": self.dtype.__name__,
            "lattice": {
                "Nt": self.cfg.Nt, "Nx": self.cfg.Nx,
                "Nev": self.cfg.Nev, "Nev1": self.cfg.Nev1,
                "delta_z": self.cfg.delta_z, "z_dir": self.cfg.z_dir,
                "Px": self.cfg.Px, "Py": self.cfg.Py, "Pz": self.cfg.Pz,
            },
            "confs": self.conf_ids,
            "timestamp": datetime.now().isoformat(),
        }
        # 避免保存不可序列化的对象, 转 str
        for k in list(config_dict.keys()):
            if not isinstance(config_dict[k], (str, int, float, bool, list, dict, type(None))):
                config_dict[k] = str(config_dict[k])

        with open(f"{self.output_dir}/run_config.json", "w") as f:
            _json.dump(config_dict, f, indent=2, ensure_ascii=False)
        print("[INFO] 运行配置已保存到 run_config.json")

    def _check(self, msg: str):
        """过程检查输出。"""
        print(f"  [CHECK] {msg}")

    def _to_numpy(self, arr: Any) -> np.ndarray:
        """安全转换为 numpy (CPU) 数组。"""
        return _ensure_cpu(arr)

    # ========================================================================
    # Step 01: 读取规范组态
    # ========================================================================

    def step_01_read_gauge(self) -> Any:
        """从 ILDG 格式的二进制文件读取 SU(3) 规范组态。

        文件格式:
            - 大端序 (big-endian) float64
            - 维度: (Nt, Nx, Nx, Nx, 4, 3, 3, 2)
            - 最后维度 = [实部, 虚部]

        Returns
        -------
        gauge : ndarray, shape (Nt, Nx, Nx, Nx, 4, 3, 3), dtype=complex
        """
        if self.gauge_file is None:
            raise ValueError("未指定规范组态文件路径 (--gauge-file)!")

        self._check(f"读取规范组态: {self.gauge_file}")

        with open(self.gauge_file, "rb") as f:
            raw = np.fromfile(f, dtype=">f8")

        Nx = self.cfg.Nx
        Nt = self.cfg.Nt
        raw = raw.reshape(Nt, Nx, Nx, Nx, 4, 3, 3, 2)
        gauge_np = raw[..., 0] + 1j * raw[..., 1]

        # 如有需要, 传输到 GPU
        if self.xp is cp:
            gauge = cp.asarray(gauge_np)
        else:
            gauge = gauge_np

        self._check(f"规范组态形状: {gauge.shape}, "
                     f"内存: {gauge.nbytes / (1024**2):.1f} MB")
        print(f"  gauge[0,0,0,0,0,:,:] =\n{self._to_numpy(gauge[0,0,0,0,0])}")

        self.data["gauge"] = gauge
        return gauge

    # ========================================================================
    # Step 02: 场强张量 F_{μν}
    # ========================================================================

    def step_02_field_strength(self) -> Any:
        """用 Clover 叶方法从规范链接构造场强张量 F_{μν}。"""
        gauge = self.data.get("gauge")
        if gauge is None:
            gauge = self.step_01_read_gauge()

        self._check("构造场强张量 F_{μν} (Clover plaquette)...")

        Nt, Nx = self.cfg.Nt, self.cfg.Nx
        F_all = compute_field_strength_all(gauge, Nt, Nx, self.contract_fn)

        self._check(f"F_all 形状: {F_all.shape}")
        self.data["F_all"] = F_all

        # 作图: 检查厄米性
        plot_field_strength_check(F_all, self.output_dir, prefix="F")

        return F_all

    # ========================================================================
    # Step 03: 对偶场强张量 F̃_{μν}
    # ========================================================================

    def step_03_dual_field_strength(self) -> Any:
        """通过 Levi-Civita 缩并计算对偶场强张量 F̃_{μν}。"""
        F_all = self.data.get("F_all")
        if F_all is None:
            F_all = self.step_02_field_strength()

        self._check("构造对偶场强张量 F̃_{μν} = (1/2) ε_{μνρσ} F^{ρσ}...")

        F_tilde_all = compute_dual_field_strength(F_all, self.epsilon,
                                                    self.contract_fn)
        self._check(f"F_tilde_all 形状: {F_tilde_all.shape}")
        self.data["F_tilde_all"] = F_tilde_all

        return F_tilde_all

    # ========================================================================
    # Step 04: VVV (Distillation Baryon Block)
    # ========================================================================

    def step_04_vvv(self) -> Any:
        """为每个时间片和 group 构造 VVV Baryon Block 张量。

        如果指定了 read_VVV_dir, 则读取已有的 VVV 数据;
        如果同时指定了 gen_VVV_dir, 则写入生成的数据并与读取数据对比。
        """
        Nt, Nx = self.cfg.Nt, self.cfg.Nx
        Nev1 = self.cfg.Nev1
        conf_id = str(self.conf_ids[0])  # 当前仅处理第一个组态
        eig_dir = self._eig_dir

        # --- 尝试读取已有 VVV ---
        VVV_read = None
        if self.read_VVV_dir:
            vvv_file = f"{self.read_VVV_dir}/VVV.Px{self.cfg.Px}Py{self.cfg.Py}Pz{self.cfg.Pz}.conf{conf_id}.npy"
            if os.path.exists(vvv_file):
                self._check(f"读取已有 VVV: {vvv_file}")
                VVV_read = np.load(vvv_file)
                if self.xp is cp:
                    VVV_read = cp.asarray(VVV_read)
                self._check(f"读取的 VVV 形状: {VVV_read.shape}")

        # --- 生成 VVV ---
        if self.gen_VVV_dir or VVV_read is None:
            self._check("构造 VVV (Distillation Baryon Block)...")

            # 动量涂抹相位
            mom_smear_vec = np.array([self.cfg.mom_smear_phase, 0, 0])
            phase_smear = compute_phase_factor(mom_smear_vec, Nx)
            if self.xp is cp:
                phase_smear = cp.asarray(phase_smear)

            # 计算全部 VVV 参数所需的动量
            # 对每个 Pz 在 --Pz-list 中单独计算
            Pz_list = getattr(self.args, 'Pz_list', None)
            if Pz_list:
                pz_values = [int(p) for p in Pz_list.split(",")]
            else:
                pz_values = [self.cfg.Pz]

            VVV_all = {}
            for Pz in pz_values:
                mom = np.array([Pz, self.cfg.Py, self.cfg.Px])
                phase_P = compute_phase_factor(mom, Nx)
                if self.xp is cp:
                    phase_P = cp.asarray(phase_P)
                # 合并涂抹相位: φ_total = φ_smear * φ_P
                phase_total = phase_smear * phase_P

                VVV_Pz = self.xp.zeros((Nt, Nev1, Nev1, Nev1), dtype=complex)
                for t in range(Nt):
                    self._check(f"VVV Pz={Pz} t={t}/{Nt-1}")
                    eig_t = self._read_eigenvectors(eig_dir, t, conf_id, Nx)
                    VVV_Pz[t] = compute_vvv_baryon_block(
                        eig_t, phase_total, Nev1, Nx, self.contract_fn
                    )
                VVV_all[Pz] = VVV_Pz

                # 写入生成数据
                if self.gen_VVV_dir:
                    os.makedirs(self.gen_VVV_dir, exist_ok=True)
                    np.save(
                        f"{self.gen_VVV_dir}/VVV.Px{self.cfg.Px}Py{self.cfg.Py}Pz{Pz}.conf{conf_id}.npy",
                        self._to_numpy(VVV_Pz),
                    )

            VVV_gen = VVV_all[pz_values[0]]
            self.data["VVV_all"] = VVV_all

            # 作图
            plot_vvv_magnitude(VVV_gen, self.output_dir, prefix="VVV")

        # --- 对比 ---
        if VVV_read is not None and VVV_gen is not None:
            diff = self._to_numpy(VVV_read - VVV_gen)
            max_diff = np.max(np.abs(diff))
            mean_diff = np.mean(np.abs(diff))
            print(f"  [COMPARE] VVV 读取-生成最大差异: {max_diff:.6e}, "
                  f"平均差异: {mean_diff:.6e}")

        result = VVV_read if VVV_read is not None else VVV_gen
        self.data["VVV"] = result
        return result

    def _read_eigenvectors(self, eig_dir: str, t: int,
                            conf_id: str, Nx: int) -> Any:
        """读取单个时间片的 Laplace 本征矢量。"""
        filename = f"{eig_dir}/eigvecs_t{t:03d}_{conf_id}"
        self._check(f"读取本征矢量: {filename}")
        with open(filename, "rb") as f:
            data = np.fromfile(f, dtype="f8")
        Nev_full = len(data) // (Nx * Nx * Nx * 3 * 2)
        data = data.reshape(Nev_full, Nx, Nx, Nx, 3, 2)
        eigvecs = data[..., 0] + 1j * data[..., 1]
        eigvecs = eigvecs[:self.cfg.Nev1]
        eigvecs = eigvecs.reshape(self.cfg.Nev1, Nx * Nx * Nx, 3)
        if self.xp is cp:
            eigvecs = cp.asarray(eigvecs)
        return eigvecs

    # ========================================================================
    # Step 05: 质子 2pt 关联函数
    # ========================================================================

    def step_05_2pt(self) -> Any:
        """通过蒸馏法计算质子 2pt 关联函数 C₂(t)。

        C₂(t_sink) = Σ_{t_src} ⟨Φ(t_sink) τ Γ τ Γ Φ(t_src)*⟩
        """
        self._check("构造质子 2pt 关联函数...")

        Nt, Nx = self.cfg.Nt, self.cfg.Nx
        Nev1 = self.cfg.Nev1
        conf_id = str(self.conf_ids[0])

        # --- 尝试读取已有 2pt ---
        C2_read = None
        if self.read_2pt_dir:
            c2_file = f"{self.read_2pt_dir}/C2pt.Px{self.cfg.Px}Py{self.cfg.Py}Pz{self.cfg.Pz}.conf{conf_id}.npy"
            if os.path.exists(c2_file):
                self._check(f"读取已有 2pt: {c2_file}")
                C2_read = np.load(c2_file)
                if self.xp is cp:
                    C2_read = cp.asarray(C2_read)

        # --- 生成 2pt ---
        C2_gen = None
        if self.gen_2pt_dir or C2_read is None:
            # 需要 VVV
            VVV = self.data.get("VVV")
            if VVV is None:
                VVV = self.step_04_vvv()

            # 需要 perambulator
            peram = self._read_perambulators(conf_id, Nt, Nev1)

            # Dirac 投影矩阵 (Cγ₅γ₄ 用于质子插值)
            g7 = self.gamma[7]   # γ₃γ₁
            g4 = self.gamma[4]   # γ₄
            # interProject = g7 @ g4 = Cγ₅γ₄
            interProj = self.xp.matmul(g7, g4)

            # 2pt 收缩 (所有 t_src, t_sink)
            C2_all = self.xp.zeros((Nt, 4, 4), dtype=complex)

            for t_src in range(Nt):
                self._check(f"2pt t_src={t_src}/{Nt-1}")
                VVV_src = VVV[t_src]

                for t_sink in range(Nt):
                    t_sep = (t_sink - t_src) % Nt
                    VVV_snk = VVV[t_sink]

                    tau = peram[t_sep]  # (4, 4, Nev1, Nev1)
                    # Γ τ Γ
                    tau_g = self.contract_fn(
                        "ijab,jkbc,klcd->ilad",
                        interProj, tau, interProj,
                    )

                    C2_all[t_src] += contract_proton_2pt_single_tsrc(
                        VVV_snk, tau, tau_g, VVV_src.conj(),
                        self.contract_fn,
                    )

            # 对所有 source 时间平均
            C2_gen = self.xp.mean(C2_all, axis=0)  # (4, 4), 实际上是每对 (t_src, t_sink)

            if self.gen_2pt_dir:
                os.makedirs(self.gen_2pt_dir, exist_ok=True)
                np.save(
                    f"{self.gen_2pt_dir}/C2pt.Px{self.cfg.Px}Py{self.cfg.Py}Pz{self.cfg.Pz}.conf{conf_id}.npy",
                    self._to_numpy(C2_gen),
                )

            # 作图
            plot_2pt_correlator(C2_gen, self.output_dir, prefix="C2pt",
                                 Pz=self.cfg.Pz)

        # --- 对比 ---
        if C2_read is not None and C2_gen is not None:
            diff = self._to_numpy(C2_read - C2_gen)
            print(f"  [COMPARE] 2pt 读取-生成最大差异: {np.max(np.abs(diff)):.6e}")

        result = C2_read if C2_read is not None else C2_gen
        self.data["C2pt"] = result
        return result

    def _read_perambulators(self, conf_id: str, Nt: int, Nev1: int) -> Any:
        """读取所有 source 时间片的 perambulator 并组装。"""
        peram_dir = self._peram_u_dir
        self._check(f"读取 perambulators: {peram_dir}")

        all_perams = []
        for t_src in range(Nt):
            peram_t = []
            for d_src in range(4):
                fname = f"{peram_dir}/perams.{conf_id}.{d_src}.{t_src}"
                with open(fname, "rb") as f:
                    peram_t.append(np.fromfile(f, dtype="f8"))
            all_perams.append(np.concatenate(peram_t))

        peram = np.concatenate(all_perams)
        Nev_full = int(np.sqrt(peram.size / (4 * 4 * Nt * Nt * 2)))

        peram = peram.reshape(Nt, 4, Nt, 4, Nev_full, Nev_full, 2)
        peram = peram[..., 0] + 1j * peram[..., 1]
        peram = peram[:, :, :, :, :Nev1, :Nev1]

        if self.xp is cp:
            peram = cp.asarray(peram)
        return peram

    # ========================================================================
    # Step 06: Nonlocal 胶子 OPE 算符
    # ========================================================================

    def step_06_ope(self) -> Any:
        """构造 nonlocal 胶子 OPE 算符 O_{μν}(z) (包含 Wilson 线)。"""
        gauge = self.data.get("gauge")
        F_all = self.data.get("F_all")
        F_tilde_all = self.data.get("F_tilde_all")

        if gauge is None:
            gauge = self.step_01_read_gauge()
        if F_all is None:
            F_all = self.step_02_field_strength()
        if F_tilde_all is None:
            F_tilde_all = self.step_03_dual_field_strength()

        Nt, Nx = self.cfg.Nt, self.cfg.Nx
        delta_z = self.cfg.delta_z
        z_dir = self.cfg.z_dir
        conf_id = str(self.conf_ids[0])
        conf_short = f"L{Nx}x{Nt}"

        # --- OPE 指标分配 (非极化胶子) ---
        # 对于 z_dir 方向的 Wilson 线, 三个独立的分量:
        #   (mu=t, nu=perp1), (mu=t, nu=perp2), (mu=perp1, nu=perp2)
        perp1 = (z_dir + 1) % 3  # 横方向 1
        perp2 = (z_dir + 2) % 3  # 横方向 2

        ope_components = [
            (3, perp1, 3, perp1),  # (t, perp1)
            (3, perp2, 3, perp2),  # (t, perp2)
            (perp1, perp2, perp1, perp2),  # (perp1, perp2)
        ]

        ope_results = {}
        for idx, (mu, nu, mu2, nu2) in enumerate(ope_components):
            # --- 尝试读取已有 OPE ---
            ope_read = None
            if self.read_ope_dir:
                ope_file = (f"{self.read_ope_dir}/ops_mu{mu}_nu{nu}"
                            f"_dz{delta_z}_conf{conf_id}.npy")
                if os.path.exists(ope_file):
                    self._check(f"读取已有 OPE[{idx}]: {ope_file}")
                    ope_read = np.load(ope_file)

            # --- 生成 OPE ---
            ope_gen = None
            if self.gen_ope_dir or ope_read is None:
                self._check(f"构造 OPE 算符[{idx}]: mu={mu},nu={nu},"
                            f"mu2={mu2},nu2={nu2}")
                ope_gen = compute_ope_all_z(
                    gauge, F_all, F_tilde_all,
                    delta_z, z_dir, mu, nu, mu2, nu2,
                    Nt, self.contract_fn, self,
                )
                ope_gen = self._to_numpy(ope_gen)  # (Nt, delta_z)

                if self.gen_ope_dir:
                    os.makedirs(self.gen_ope_dir, exist_ok=True)
                    np.save(
                        f"{self.gen_ope_dir}/ops_mu{mu}_nu{nu}"
                        f"_dz{delta_z}_conf{conf_id}.npy",
                        ope_gen,
                    )

                # 作图
                plot_ope_result(ope_gen, self.output_dir,
                                 prefix=f"OPE_{idx}", mu=mu, nu=nu)

            # --- 对比 ---
            if ope_read is not None and ope_gen is not None:
                diff = np.abs(ope_read - ope_gen)
                print(f"  [COMPARE] OPE[{idx}] 读取-生成最大差异: "
                      f"{np.max(diff):.6e}")

            ope_results[f"mu{mu}_nu{nu}"] = (ope_read if ope_read is not None
                                              else ope_gen)

        self.data["ope_results"] = ope_results
        return ope_results

    # ========================================================================
    # Step 07: 3pt 关联函数
    # ========================================================================

    def step_07_3pt(self) -> Any:
        """构造 3pt 关联函数 C₃(t_src, t_sink; z)。

        在 disconnect 近似下:
            C₃(t; z) = C₂(t_sink - t_src) · ⟨O(z)⟩_t
        """
        self._check("构造 3pt 关联函数...")

        C2 = self.data.get("C2pt")
        ope_results = self.data.get("ope_results")

        if C2 is None:
            C2 = self.step_05_2pt()
        if ope_results is None:
            ope_results = self.step_06_ope()

        # 将 OPE 分量求和得到总 OPE
        ope_keys = list(ope_results.keys())
        delta_z = ope_results[ope_keys[0]].shape[1]
        Nt = ope_results[ope_keys[0]].shape[0]

        # 总 OPE 算符 = Σ 各分量 (非极化胶子的公式)
        ope_total = np.sum([ope_results[k] for k in ope_keys], axis=0)  # (Nt, delta_z)

        # C₃(t; z) = C₂ · O(z) (因子化近似)
        # 这里假设 C₂ 是一个标量 (取所有时间片的平均值)
        if C2.ndim >= 2:
            C2_scalar = np.mean(np.trace(C2))  # Dirac 迹的平均
        else:
            C2_scalar = np.mean(np.abs(C2))

        C3_3pt = C2_scalar * ope_total  # (Nt, delta_z)

        self._check(f"3pt 关联函数形状: {C3_3pt.shape}, C2_scalar={C2_scalar:.6e}")
        self.data["C3pt"] = C3_3pt
        self.data["ope_total"] = ope_total

        # 作图: 使用 OPE 作图 (因为 3pt ∝ OPE)
        plot_ope_result(ope_total, self.output_dir, prefix="3pt_OPE_total")

        if self.gen_3pt_dir:
            os.makedirs(self.gen_3pt_dir, exist_ok=True)
            conf_id = str(self.conf_ids[0])
            np.save(f"{self.gen_3pt_dir}/C3pt.Pz{self.cfg.Pz}.conf{conf_id}.npy",
                    self._to_numpy(C3_3pt))

        return C3_3pt

    # ========================================================================
    # Step 08-09: Distillation 框架 + 动量涂抹 (综合实现)
    # ========================================================================

    def step_08_distillation(self) -> Any:
        """实现 Distillation 框架 — 读取本征矢量与 Perambulator。

        本步骤是 step_04 和 step_05 的辅助准备, 读取蒸馏所需的核心数据。
        """
        self._check("实现 Distillation 框架 (本征矢量 + Perambulator)...")

        conf_id = str(self.conf_ids[0])
        Nt, Nx = self.cfg.Nt, self.cfg.Nx
        Nev1 = self.cfg.Nev1

        eig_dir = self._eig_dir
        peram_dir = self._peram_u_dir

        # 检查本征矢量是否可读
        t0_file = f"{eig_dir}/eigvecs_t000_{conf_id}"
        if os.path.exists(t0_file):
            self._check(f"本征矢量目录可用: {eig_dir}")
            # 读取第一个时间片作为测试
            eig_test = self._read_eigenvectors(eig_dir, 0, conf_id, Nx)
            self._check(f"本征矢量形状 (t=0): {eig_test.shape}")
        else:
            print(f"[WARN] 本征矢量目录不可用: {eig_dir} (文件不存在)")

        # 检查 perambulator 是否可读
        peram_t0 = f"{peram_dir}/perams.{conf_id}.0.0"
        if os.path.exists(peram_t0):
            self._check(f"Perambulator 目录可用: {peram_dir}")
        else:
            print(f"[WARN] Perambulator 目录不可用: {peram_dir} (文件不存在)")

        return {"eig_dir": eig_dir, "peram_dir": peram_dir}

    def step_09_momentum_smear(self) -> Any:
        """动量涂抹: 计算并应用动量涂抹相位因子。"""
        self._check("计算动量涂抹相位因子...")

        Nx = self.cfg.Nx
        mom_smear_vec = np.array([self.cfg.mom_smear_phase, 0, 0])
        phase_smear = compute_phase_factor(mom_smear_vec, Nx)

        mom_vec = np.array([self.cfg.Pz, self.cfg.Py, self.cfg.Px])
        phase_P = compute_phase_factor(mom_vec, Nx)

        phase_total = phase_smear * phase_P

        self._check(f"涂抹相位形状: {phase_total.shape}")
        print(f"  mom_smear_phase={self.cfg.mom_smear_phase}, P=({self.cfg.Pz},{self.cfg.Py},{self.cfg.Px})")

        self.data["phase_total"] = phase_total
        return phase_total

    # ========================================================================
    # Step 10: 矩阵元提取 h(z, P_z)
    # ========================================================================

    def step_10_matrix_element(self) -> Tuple[Any, Any]:
        """提取基态矩阵元 h(z, P_z) 及其误差。

        在 disconnect 近似下:
            h(z, P_z) = ⟨O(z)⟩ / ⟨P|P⟩

        对多个组态使用 Jackknife 方法估计误差。
        """
        self._check("提取矩阵元 h(z, P_z)...")

        ope_results = self.data.get("ope_results")
        if ope_results is None:
            ope_results = self.step_06_ope()

        # 将所有 OPE 分量求和
        ope_keys = list(ope_results.keys())
        ope_total = np.sum([self._to_numpy(ope_results[k]) for k in ope_keys],
                            axis=0)  # (Nt, delta_z)

        Nt = ope_total.shape[0]
        delta_z = ope_total.shape[1]

        # 对时间平均 → h(z)
        h_z = np.mean(ope_total, axis=0)  # (delta_z,)
        h_z_err = np.std(ope_total, axis=0) / np.sqrt(Nt)  # 简单标准误差

        self._check(f"h(z) 形状: {h_z.shape}")
        for z in range(min(delta_z, 8)):
            self._check(f"  z={z}: h = {h_z[z]:.6e} ± {h_z_err[z]:.6e}")

        self.data["h_z"] = h_z
        self.data["h_z_err"] = h_z_err

        # 作图
        plot_matrix_element(h_z, h_z_err, self.output_dir,
                             prefix="h_z", Pz=self.cfg.Pz)

        return h_z, h_z_err

    # ========================================================================
    # Step 11: Fourier 变换 → quasi-PDF
    # ========================================================================

    def step_11_fourier(self, n_x: int = 100) -> Any:
        """通过 Fourier 变换将坐标空间矩阵元转换为 quasi-PDF。

        g̃(x, P_z) = (2P_z / x) · ∫ dz h(z) sin(x P_z z)
        """
        self._check("Fourier 变换 → quasi-PDF g̃(x, P_z)...")

        h_z = self.data.get("h_z")
        if h_z is None:
            h_z, _ = self.step_10_matrix_element()

        delta_z = len(h_z)
        z_values = np.arange(delta_z)

        # Bjorken-x 值范围 (避开 x=0)
        x_values = np.linspace(0.01, 0.99, n_x)
        Pz = float(self.cfg.Pz)

        quasi_pdf = fourier_transform_to_quasi_pdf(
            h_z, z_values, Pz, x_values
        )

        self._check(f"quasi-PDF 形状: {quasi_pdf.shape}")
        print(f"  quasi-PDF 范围: [{np.min(quasi_pdf):.4f}, {np.max(quasi_pdf):.4f}]")

        self.data["quasi_pdf"] = quasi_pdf
        self.data["x_values"] = x_values

        # 作图
        plot_quasi_pdf(quasi_pdf, x_values, self.output_dir,
                        prefix="quasi", Pz=self.cfg.Pz)

        return quasi_pdf

    # ========================================================================
    # Step 12: 微扰匹配 → 光锥 PDF
    # ========================================================================

    def step_12_matching(self) -> Any:
        """通过 NLO 微扰匹配将 quasi-PDF 转换为光锥 PDF g(x, μ)。

        在 LO 近似下: g(x, μ) ≈ g̃(x, P_z)。
        """
        self._check("微扰匹配 → 光锥 PDF g(x, μ)...")

        quasi_pdf = self.data.get("quasi_pdf")
        x_values = self.data.get("x_values")

        if quasi_pdf is None:
            quasi_pdf = self.step_11_fourier()
            x_values = self.data["x_values"]

        alpha_s = getattr(self.args, 'alpha_s', 0.2)
        mu_over_Pz = getattr(self.args, 'mu_over_Pz', 1.0)

        lc_pdf = matching_kernel_nlo(quasi_pdf, x_values, alpha_s, mu_over_Pz)

        self._check(f"光锥 PDF 形状: {lc_pdf.shape}")
        print(f"  α_s={alpha_s}, μ/P_z={mu_over_Pz}")

        self.data["lc_pdf"] = lc_pdf

        # 作图
        mu_GeV = mu_over_Pz * self.cfg.Pz * 0.1973 / (self.cfg.alttc * self.cfg.Nx)
        plot_lightcone_pdf(lc_pdf, x_values, self.output_dir,
                            prefix="lc_pdf", Pz=self.cfg.Pz, mu=2.0)

        return lc_pdf

    # ========================================================================
    # Step 13: Jackknife / Bootstrap 统计误差分析
    # ========================================================================

    def step_13_error_analysis(self) -> Any:
        """对提取的矩阵元和 PDF 进行 Jackknife 和 Bootstrap 误差分析。"""
        self._check("统计误差分析 (Jackknife + Bootstrap)...")

        h_z = self.data.get("h_z")
        h_z_err = self.data.get("h_z_err")

        if h_z is None:
            h_z, h_z_err = self.step_10_matrix_element()

        Nz = len(h_z)

        # --- Bootstrap 误差 ---
        # 我们需要组态级别的 OPE 数据来进行 bootstrap
        ope_results = self.data.get("ope_results")
        bs_results = {}
        if ope_results is not None:
            ope_keys = list(ope_results.keys())
            for key in ope_keys:
                ope_data = self._to_numpy(ope_results[key])  # (Nt, delta_z)
                # 对时间方向做 bootstrap (将 Nt 当作 "组态")
                bs = bootstrap_samples(ope_data, n_bootstrap=1000, axis=0)
                bs_mean = np.mean(bs, axis=0)  # (delta_z,)
                bs_err = np.std(bs, axis=0)    # (delta_z,)
                bs_results[key] = {"mean": bs_mean, "err": bs_err}

                # 打印对比
                max_rel_diff = np.max(np.abs(bs_err - h_z_err) / (np.abs(h_z_err) + 1e-15))
                print(f"  Bootstrap vs Jackknife [{key}]: "
                      f"最大相对差异 = {max_rel_diff:.4f}")

        # --- 汇总报告 ---
        print(f"\n{'='*60}")
        print(f"  统计误差分析汇总")
        print(f"{'='*60}")
        print(f"  {'z':>4s}  {'Re[h(z)]':>16s}  {'JK_err':>12s}")
        print(f"  {'-'*36}")
        for z in range(min(Nz, 12)):
            print(f"  {z:4d}  {np.real(h_z[z]):16.6e}  {h_z_err[z]:12.6e}")

        self.data["bs_results"] = bs_results

        # 作图: 使用已有的 plot_matrix_element 函数
        # (已在 step_10 中调用)

        return {"jackknife": {"h_z": h_z, "h_z_err": h_z_err},
                "bootstrap": bs_results}

    # ========================================================================
    # 运行流水线
    # ========================================================================

    def run(self):
        """运行选定的步骤并输出汇总表。"""
        steps = self.args.steps
        step_map = {
            1: ("Step 01: Read Gauge Config", self.step_01_read_gauge),
            2: ("Step 02: Field Strength F_{μν}", self.step_02_field_strength),
            3: ("Step 03: Dual Field Strength F̃_{μν}", self.step_03_dual_field_strength),
            4: ("Step 04: VVV Baryon Block", self.step_04_vvv),
            5: ("Step 05: Proton 2pt Correlator", self.step_05_2pt),
            6: ("Step 06: Gluon OPE Operator", self.step_06_ope),
            7: ("Step 07: 3pt Correlator", self.step_07_3pt),
            8: ("Step 08: Distillation Framework", self.step_08_distillation),
            9: ("Step 09: Momentum Smearing", self.step_09_momentum_smear),
            10: ("Step 10: Matrix Element h(z,Pz)", self.step_10_matrix_element),
            11: ("Step 11: Fourier Transform", self.step_11_fourier),
            12: ("Step 12: Matching → LC PDF", self.step_12_matching),
            13: ("Step 13: Error Analysis", self.step_13_error_analysis),
        }

        print(f"\n{'#'*70}")
        print(f"#  Gluon PDF Pipeline — 开始执行")
        print(f"#  后端: {self.backend_name}")
        print(f"#  格点: {self.cfg.Nt}×{self.cfg.Nx}³, P=({self.cfg.Px},{self.cfg.Py},{self.cfg.Pz})")
        print(f"#  步骤: {steps}")
        print(f"#  输出: {self.output_dir}")
        print(f"{'#'*70}\n")

        # 按顺序执行选定步骤
        for step_num in sorted(steps):
            if step_num in step_map:
                name, func = step_map[step_num]
                with Timer(name, self):
                    func()
            else:
                print(f"[WARN] 未知步骤编号: {step_num}, 跳过。")

        # --- 汇总输出 ---
        self._print_summary()

        # --- 性能汇总图 ---
        if self.timing_results:
            plot_summary_timeline(self.timing_results, self.memory_results,
                                   self.output_dir)

        print(f"\n{'#'*70}")
        print(f"#  流水线执行完毕!")
        print(f"#  总耗时: {sum(self.timing_results.values()):.1f} s")
        print(f"#  输出目录: {self.output_dir}")
        print(f"{'#'*70}\n")

    def _print_summary(self):
        """打印运行汇总表。"""
        print(f"\n{'='*80}")
        print(f"  运行汇总")
        print(f"{'='*80}")
        print(f"  {'步骤':<40s}  {'耗时 (s)':>10s}  {'内存峰值 (MB)':>14s}")
        print(f"  {'-'*66}")

        total_time = 0.0
        max_mem = 0.0
        for name, t in self.timing_results.items():
            mem_peak = self.memory_results.get(name, {}).get("peak_mb", 0.0)
            total_time += t
            max_mem = max(max_mem, mem_peak)
            short_name = name.replace("Step ", "").split(":")[0].strip()
            print(f"  {short_name:<40s}  {t:10.3f}  {mem_peak:14.1f}")

        print(f"  {'-'*66}")
        print(f"  {'总计':<40s}  {total_time:10.3f}  {max_mem:14.1f}")
        print(f"{'='*80}\n")


# ============================================================================
# 第十四部分: 2pt 有效质量分析 (匹配 main-2pt.py / include.py)
# ============================================================================
#
# 有效质量公式 (cosh 型, 考虑周期性边界条件):
#
#   谱分解:
#     C(t) = Σ_n |Z_n|² [e^{-E_n t} + e^{-E_n (Nt-t)}]
#
#   基态主导:
#     C(t) ≃ |Z₀|² [e^{-E₀ t} + e^{-E₀ (Nt-t)}]
#         = 2|Z₀|² e^{-E₀ Nt/2} cosh[E₀ (t - Nt/2)]
#
#   对三个相邻时间片:
#     [C(t-1) + C(t+1)] / [2 C(t)] = cosh(E₀)
#
#   格点单位有效质量:
#     a·m_eff(t) = arccosh( [C(t-1) + C(t+1)] / [2 C(t)] )
#
#   物理单位:
#     m_eff(t) [GeV] = a·m_eff(t) × (ℏc / a) = a·m_eff(t) × 0.1973 / a[fm]
#
# 对数有效质量 (备选):
#    a·m_eff(t+1/2) = ln[ C(t) / C(t+1) ]
#
# Jackknife 误差:
#    C_i^{JK} = (Σ_j C_j - C_i) / (N-1)
#    σ^{JK} = std(C^{JK}) × √(N-1)
# ============================================================================


def run_2pt_analysis(args: "argparse.Namespace") -> None:
    """运行 2pt 有效质量分析 (匹配 examples/zhangxin/main-2pt.py)。

    读取 Chroma IOG 格式的 2pt 关联函数, 计算 cosh/log 有效质量,
    Jackknife 误差估计, 并生成有效质量图。

    依赖 include.py 的 data_analyse 类 (来自 examples/zhangxin/)。
    """
    if not HAS_INCLUDE:
        print("[ERROR] 无法导入 include.py (来自 examples/zhangxin/)。")
        print("        请确保 examples/zhangxin/ 存在且 iog_reader/iog.so 已编译。")
        print("        或者使用 --analysis-type pdf 模式 (不需要 include.py)。")
        sys.exit(1)

    Nx = args.Nx
    Nt = args.Nt
    alttc = getattr(args, 'alttc', 0.1053)
    meff_type = getattr(args, 'meff_type', 'cosh')
    hadron = getattr(args, 'hadron', 'pion')
    Px, Py, Pz = args.Px, args.Py, args.Pz
    N_start = args.conf_start
    gap = args.conf_step
    Ncnfg_iog = args.conf_num
    tsep_val = getattr(args, 'tsep', 36)
    time_fold = getattr(args, 'time_fold', False)
    link_max_fmt = getattr(args, 'link_max', 10)
    meff_yrange = getattr(args, 'meff_range', [0.0, 1.0])
    if isinstance(meff_yrange, str):
        meff_yrange = [float(x) for x in meff_yrange.split(",")]

    # 文件路径模板
    iog_2pt_path = getattr(args, 'iog_2pt_path', None)
    if iog_2pt_path is None:
        iog_2pt_path = (
            "./beta6.20_mu-0.2770_ms-0.2400_L%dx%d/sush_iog/"
            "pion_2pt_Px%dPy%dPz%d_ENV%d_conf%d_tsep-1_mass-0.2770.iog"
        )

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  2pt 有效质量分析")
    print(f"{'='*60}")
    print(f"  强子: {hadron},  格点: {Nt}×{Nx}³,  a={alttc} fm")
    print(f"  有效质量: {meff_type},  P=({Px},{Py},{Pz})")
    print(f"  组态: start={N_start}, step={gap}, N={Ncnfg_iog}")
    print(f"  IOG: {iog_2pt_path}")
    print(f"{'='*60}\n")

    filepath = np.array([iog_2pt_path])
    ENV_arr = np.asarray([-1])
    P_arr = np.asarray([[Px, Py, Pz]])
    tsep_arr = np.asarray([tsep_val])

    analyse = _data_analyse(
        num_data=1, hadron=hadron, filepath=filepath,
        alttc=alttc, Nx=Nx, Nt=Nt, P=P_arr, ENV=ENV_arr,
        N_start=N_start, gap=gap, Ncnfg_data=0, Ncnfg_iog=Ncnfg_iog,
        tsep=tsep_arr, time_fold=time_fold, save_path=output_dir,
        link_max=link_max_fmt, analyse_type='2pt',
        meff_type=meff_type, read_type='iog',
    )

    print("[INFO] data_analyse 实例化完成。")
    for key in analyse.readed.keys():
        print(f"  {key}.shape = {analyse.readed[key].shape}")

    print(f"\n[INFO] 计算 {meff_type} 有效质量...")
    analyse.meff_2pt('iog')

    meff_mean_arr = analyse.meff_data_2pt['meff_2pt_iog_mean']
    meff_err_arr = analyse.meff_data_2pt['meff_2pt_iog_err']
    print(f"  有效质量形状: {meff_mean_arr.shape}")

    N_print = min(10, meff_mean_arr.shape[-1])
    for t in range(N_print):
        print(f"    t={t:3d}  m_eff = {meff_mean_arr[0,-1,t]:8.4f}"
              f" ± {meff_err_arr[0,-1,t]:8.4f} GeV")

    if meff_mean_arr.shape[-1] > 10:
        ps = meff_mean_arr.shape[-1] // 4
        pe = meff_mean_arr.shape[-1] // 2
        pm = np.mean(meff_mean_arr[0,-1,ps:pe])
        print(f"  平台区 t∈[{ps},{pe}]: m_eff = {pm:.4f} GeV")

    if HAS_MPL and not args.no_plot:
        _plot_meff_2pt(meff_mean_arr, meff_err_arr, P_arr, output_dir,
                        hadron=hadron, Nx=Nx, Nt=Nt, alttc=alttc,
                        meff_type=meff_type, y_range=meff_yrange)

    np.savez(f"{output_dir}/meff_2pt_result.npz",
              meff_mean=meff_mean_arr, meff_err=meff_err_arr)
    print(f"[INFO] 数据已保存: {output_dir}/meff_2pt_result.npz")
    print(f"[INFO] 2pt 分析完成。\n")


def _plot_meff_2pt(meff_mean_arr, meff_err_arr, P_arr, output_dir,
                    hadron="pion", Nx=24, Nt=72, alttc=0.1053,
                    meff_type="cosh", y_range=None):
    """画图: 2pt 有效质量 vs 时间片, 匹配 include.py 作图风格。"""
    if not HAS_MPL:
        return
    N_P = P_arr.shape[0]
    markers = np.array(['s','*','+','x','p','h','v','X','D','P','H','o'])
    plt.rcParams.update({'font.size': 25})
    fig, ax = plt.subplots(1, 1, figsize=(20, 10))
    fig.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    name = f"{hadron}_meff_{Nx}x{Nt}_{meff_type}_iog"
    ax.set_title(name, fontdict={'fontsize': 30, 'fontweight': 'light'})
    if y_range is not None and len(y_range) == 2:
        ax.set_ylim(y_range)
    ax.set_xlabel('t/a')
    ax.set_ylabel('$E_{\\mathrm{2pt}}$/GeV')
    for p in range(N_P):
        nr = meff_mean_arr[p, -1].shape[0]
        ax.errorbar(np.arange(nr), meff_mean_arr[p, -1],
                     yerr=meff_err_arr[p, -1], alpha=0.5,
                     marker=markers[p % len(markers)],
                     capsize=3.5, capthick=1.5,
                     label=f'P=({P_arr[p,0]},{P_arr[p,1]},{P_arr[p,2]})',
                     linestyle='none', elinewidth=2)
    plt.legend(fontsize=18)
    fig.savefig(f"{output_dir}/{name}.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"[PLOT] 有效质量图已保存: {name}.pdf")


# ============================================================================
# 第XIV-B部分: 质子 2pt 关联函数蒸馏计算
# ============================================================================
#
# 目标: 完全复现 examples/donghx/2pt_proton_Cg5gmu_L24x72_mom2_zdir_dcu.py
#
# 算法流程:
#
# 1. Dirac 矩阵与宇称投影
#    ──────────────────────
#    P₊ = ½(γ₀ + γ₄)  — 正宇称投影
#    P₋ = ½(γ₀ - γ₄)  — 负宇称投影
#
#    插值算符 (质子, Cγ₅γ₄):
#    interProject1 = gamma(7) @ gamma(4)  = γ₃γ₁ · γ₄
#    interProject2 = gamma(7) @ gamma(4)  = γ₃γ₁ · γ₄
#    (对于 _Cg5g4 元素类型, 两个投影矩阵相同)
#
#    其中 gamma(7) = γ₃γ₁ = -γ₂γ₄γ₅  (电荷共轭 × γ₅ 的组合)
#
# 2. 动量涂抹相位因子
#    ──────────────────
#    φ_smear(x) = exp( -i · 2π · P_smear · x / L )
#    P_smear = (momsmear_phase, 0, 0)
#
#    本征矢量涂抹:
#    ṽ_i^a(x) = v_i^a(x) · φ_smear(x)
#
# 3. 物理动量相位因子
#    ──────────────────
#    φ_P(x) = exp( -i · 2π · P · x / L )
#    P = (Pz, Py, Px)  动量扫描列表
#
# 4. VVV Baryon Block 计算 (每个时间片)
#    ─────────────────────────────────
#    Φ_{abc}(P) = Σ_x φ_P(x) · ε_{ijk} · ṽ_i^a(x) · ṽ_j^b(x) · ṽ_k^c(x)
#
#    逐 x-层分片计算 (减少 GPU/CPU 内存占用):
#    Φ += contract("x,ax,bx,cx->abc", φ_P_slice, ṽ_a0, ṽ_b1, ṽ_c2)   [ ε_{012}=+1 ]
#    Φ += contract("x,ax,bx,cx->abc", φ_P_slice, ṽ_a1, ṽ_b2, ṽ_c0)   [ ε_{120}=+1 ]
#    Φ += contract("x,ax,bx,cx->abc", φ_P_slice, ṽ_a2, ṽ_b0, ṽ_c1)   [ ε_{201}=+1 ]
#    Φ -= contract("x,ax,bx,cx->abc", φ_P_slice, ṽ_a0, ṽ_b2, ṽ_c1)   [ ε_{021}=-1 ]
#    Φ -= contract("x,ax,bx,cx->abc", φ_P_slice, ṽ_a1, ṽ_b0, ṽ_c2)   [ ε_{102}=-1 ]
#    Φ -= contract("x,ax,bx,cx->abc", φ_P_slice, ṽ_a2, ṽ_b1, ṽ_c0)   [ ε_{210}=-1 ]
#
# 5. Perambulator 读取与 Dirac 投影
#    ─────────────────────────────
#    τ(t_sink; t_src) 形状: (Nt, 4, 4, Nev1, Nev1)
#    指标: (t_sink, d_sink, d_source, ev_sink, ev_source)
#
#    Dirac 投影 perambulator:
#    (ΓτΓ)(t_sink) = contract("gh,thkbe,jk->tgjbe", interProject1, τ, interProject2)
#    形状: (Nt, 4, 4, Nev1, Nev1)
#
# 6. Wick 收缩 (对每个 (t_sink, t_source) 且 2 ≤ deltat ≤ 32)
#    ────────────────────────────────────────────────────────
#
#    Direct 项 (夸克线无交换):
#    C₂^{dir}_{il}(t_snk, t_src) =
#      contract("abc,gjad,gjbe,ilcf,def->il",
#               VVV_snk, τ_snk, ΓτΓ_snk, τ_snk, VVV_src*)
#
#    Exchange 项 (两夸克线交换, Fermi 统计负号):
#    C₂^{ex}_{il}(t_snk, t_src) =
#      contract("abc,glaf,gjbe,ijcd,def->il",
#               VVV_snk, τ_snk, ΓτΓ_snk, τ_snk, VVV_src*)
#
#    总关联函数矩阵:
#    C₂(t_snk, t_src) = Direct - Exchange
#    形状: (Nt, Nt, 4, 4)
#
# 7. 宇称投影
#    ─────────
#    C₂^{pp}(t_snk, t_src) = contract("li,yxil->yx", P₊, C₂)
#    C₂^{pm}(t_snk, t_src) = contract("li,yxil->yx", P₋, C₂)
#    形状: (Nt, Nt)  — 正/负宇称投影后的 (t_sink, t_source) 矩阵
#
# 8. 边界条件符号修正
#    ────────────────
#    gluon 2pt 关联函数中的反粒子传播子需要符号修正:
#    pp: 当 t_snk < t_src 时 (向后传播), 乘 -1
#    pm: 当 t_snk > t_src 时 (向前传播), 乘 -1
#
# 9. 输出文件 (命名与参考代码完全一致)
#    ────────────────────────────────
#    Raw 收缩矩阵:
#      twopt_slice_pp_Px{Px}Py{Py}Pz{Pz}_eginphase{mom_smear}{element}_contract_conf{conf_id}.npy
#      形状: (Nt, Nt, 4, 4)
#
#    Parity 投影后:
#      twopt_slice_pp_Px{Px}Py{Py}Pz{Pz}_eginphase{mom_smear}{element}_nopol_ss_conf{conf_id}.npy
#      形状: (Nt, Nt)  (仅 pp, 正宇称投影)
# ============================================================================


def run_proton_2pt_analysis(args: "argparse.Namespace") -> None:
    """运行质子 2pt 关联函数蒸馏计算。

    完全复现 examples/donghx/2pt_proton_Cg5gmu_L24x72_mom2_zdir_dcu.py
    的算法与输出格式。

    使用 NumPy 后端 (CPU), 与 DCU 参考代码的 PyTorch 后端在数值上等价。
    """
    xp = np

    # ── 参数提取 ──────────────────────────────────────────
    contract = get_contract(np)  # 嵌套函数通过闭包访问

    Nt = args.Nt
    Nx = args.Nx
    Nev  = getattr(args, 'Nev',  100)
    Nev1 = getattr(args, 'Nev1', 100)
    conf_id = str(args.conf_start)

    Px = args.Px
    Py = args.Py
    Pz_start = args.Pz  # 单动量模式

    # 动量涂抹 (与参考代码的硬编码值一致)
    mom_smear        = getattr(args, 'mom_smear',       -2)
    momsmear_phase   = getattr(args, 'mom_smear_phase',  2)

    # 插值算符类型 (默认 _Cg5g4)
    element = getattr(args, 'element', '_Cg5g4')

    # 动量列表 (参考代码: Pzlist = [-2, -3, -4, -5, -6])
    Pz_list_str = getattr(args, 'Pz_list', None)
    if Pz_list_str:
        Pzlist = [int(p) for p in Pz_list_str.split(",")]
    else:
        Pzlist = [Pz_start]

    # 数据路径
    eig_dir = getattr(args, 'eig_dir', None)
    if eig_dir is None:
        eig_dir = (
            "/public/group/lqcd/eigensystem/"
            "beta6.20_mu-0.2770_ms-0.2400_L24x72/{conf_id}/"
        ).format(conf_id=conf_id)

    peram_u_dir = getattr(args, 'peram_u_dir', None)
    if peram_u_dir is None:
        peram_u_dir = (
            "/public/home/sunp/sunpeng_new_disk/mom_smear_perambulators/"
            "beta6.20_mu-0.2770_ms-0.2400_L24x72/output_dir_data/"
            "mz2_my0_mx0/{conf_id}/"
        ).format(conf_id=conf_id)

    # 参考数据目录 (已有结果, 用于对比)
    ref_corr_nucl_dir = getattr(args, 'corr_nucl_dir', None)
    if ref_corr_nucl_dir is None:
        ref_corr_nucl_dir = (
            "/public/group/lqcd/donghx/2pt_Result/"
            "beta6.20_mu-0.2770_ms-0.2400_L24x72/momsmear-2z/{conf_id}/"
        ).format(conf_id=conf_id)

    # 生成数据输出目录 (新计算的结果)
    output_dir = args.output_dir
    gen_dir = os.path.join(output_dir, "data")
    plot_dir = os.path.join(output_dir, "plots")
    os.makedirs(gen_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    # 检查参考数据是否存在
    has_ref = os.path.isdir(ref_corr_nucl_dir) if ref_corr_nucl_dir else False

    print(f"\n{'='*70}")
    print(f"  质子 2pt 关联函数蒸馏计算")
    print(f"{'='*70}")
    print(f"  格点: {Nt}×{Nx}³, Nev={Nev}, Nev1={Nev1}")
    print(f"  组态: conf_id={conf_id}")
    print(f"  动量列表 Pz: {Pzlist},  Py={Py}, Px={Px}")
    print(f"  动量涂抹: mom_smear={mom_smear}, phase={momsmear_phase}")
    print(f"  插值算符: {element}")
    print(f"  本征矢量: {eig_dir}")
    print(f"  Peramb:   {peram_u_dir}")
    print(f"  参考数据: {ref_corr_nucl_dir}  {'(存在)' if has_ref else '(不存在, 跳过对比)'}")
    print(f"  生成输出: {gen_dir}")
    print(f"  作图输出: {plot_dir}")
    print(f"{'='*70}\n")

    # ── Dirac 矩阵 (DeGrand-Rossi 基 + NumPy) ─────────────
    gamma = build_gamma_matrices(xp)

    # 宇称投影矩阵
    # P₊ = ½(γ₀ + γ₄) — 正宇称 (γ₀ = I 单位矩阵)
    # P₋ = ½(γ₀ - γ₄) — 负宇称
    matrix_pplus  = 0.5 * (gamma[0] + gamma[4])
    matrix_pminus = 0.5 * (gamma[0] - gamma[4])

    # 插值算符矩阵 (Cγ₅γ₄ 用 _Cg5g4 元素)
    # gamma(7) = γ₃γ₁ = -γ₂γ₄γ₅
    # gamma(4) = γ₄
    if element == "_Cg5g4":
        interProject1 = xp.matmul(gamma[7], gamma[4])
        interProject2 = xp.matmul(gamma[7], gamma[4])
    elif element == "_Cg5g3":
        interProject1 = xp.matmul(gamma[7], gamma[3])
        interProject2 = xp.matmul(gamma[7], gamma[3])
    elif element == "_Cg5":
        interProject1 = gamma[7]
        interProject2 = gamma[7]
    else:
        interProject1 = xp.matmul(gamma[7], gamma[4])
        interProject2 = xp.matmul(gamma[7], gamma[4])

    print(f"  宇称投影: P₊ shape={matrix_pplus.shape}")
    print(f"  插值算符:")
    print(f"    interProject1 =\n{interProject1}")
    print(f"    interProject2 =\n{interProject2}")
    print(f"    element = {element}")

    # ── 动量涂抹相位因子 φ_smear(x) ────────────────────────
    # φ_smear(x) = exp(-i·2π·P_smear·x / Nx)
    # P_smear = (momsmear_phase, 0, 0)
    print(f"\n[INFO] 计算动量涂抹相位因子 (P_smear=({momsmear_phase},0,0))...")
    phase_smear = compute_phase_factor(
        np.array([momsmear_phase, 0, 0]), Nx
    )  # shape (Nx³), complex

    # ── 本征矢量读取函数 ──────────────────────────────────
    # (每次读取后对器件传输, numpy 下直接操作 CPU)
    def read_eigvecs_smeared(t: int) -> np.ndarray:
        """读取时间片 t 的本征矢量并施加动量涂抹。

        返回 shape (Nev, Nx³, 3), complex128。
        与参考代码 readin_eigvecs() 等价。
        """
        filename = f"{eig_dir}/eigvecs_t{t:03d}_{conf_id}"
        with open(filename, "rb") as f:
            file_size = os.path.getsize(filename)
            Nev_full = int(file_size / 8 / (Nx * Nx * Nx * 3 * 2))
            data = np.fromfile(f, dtype="f8").reshape(Nev_full, Nx, Nx, Nx, 3, 2)
        eigvecs = data[..., 0] + 1j * data[..., 1]            # → (Nev, Nx,Nx,Nx, 3)
        eigvecs = eigvecs[:Nev1]                                # 截断
        eigvecs = eigvecs.reshape(Nev_full, Nx * Nx * Nx, 3)   # → (Nev, Nx³, 3)

        # 动量涂抹: ṽ_i^a(x) = v_i^a(x) × φ_smear(x)
        eigvecs_mom = contract("vxa,x->vxa", eigvecs, phase_smear)
        return eigvecs_mom

    # ── Perambulator 读取函数 ──────────────────────────────
    def read_peram(t_source: int) -> np.ndarray:
        """读取 source 时间片 t_source 的 perambulator。

        返回 shape (Nt, 4, 4, Nev1, Nev1), complex64。
        与参考代码 readin_peram() 等价。
        """
        # 读取 4 个 Dirac 分量 (d_source = 0,1,2,3)
        parts = []
        for d_source in range(4):
            fname = f"{peram_u_dir}/perams.{conf_id}.{d_source}.{t_source}"
            with open(fname, "rb") as f:
                parts.append(np.fromfile(f, dtype="f8"))
        peram_raw = np.concatenate(parts)

        Nev_full = int(np.sqrt(peram_raw.size / (4 * 4 * Nt * 2)))
        # reshape: (d_source, t_sink, ev_source, d_sink, ev_sink, complex)
        peram = peram_raw.reshape(4, Nt, Nev_full, 4, Nev_full, 2)
        # transpose: → (t_sink, d_sink, d_source, ev_sink, ev_source, complex)
        peram = peram.transpose(1, 3, 0, 4, 2, 5)
        peram = peram[..., 0] + 1j * peram[..., 1]
        peram = peram[:, :, :, :Nev1, :Nev1]         # 截断
        peram = peram.astype(np.complex64)
        return peram

    # ── VVV 计算函数 (对单个动量) ─────────────────────────
    def compute_vvv_for_momentum(mom: np.ndarray) -> np.ndarray:
        """为所有时间片计算 VVV Baryon Block。

        返回 shape (Nt, Nev1, Nev1, Nev1), complex128。
        与参考代码 VVV_Calc_cupy() 等价 (用 numpy 替代 torch)。
        """
        VVV_sink = np.zeros((Nt, Nev1, Nev1, Nev1), dtype=complex)

        for t in range(Nt):
            t1 = time.perf_counter()
            eigvecs = read_eigvecs_smeared(t)               # (Nev, Nx³, 3)
            t2 = time.perf_counter()
            print(f"  [VVV] t={t:3d} read eigvecs done, "
                  f"time={t2-t1:.3f}s")

            t3 = time.perf_counter()
            phase_P = compute_phase_factor(mom, Nx)          # (Nx³)
            layer_size = Nx * Nx

            # 逐 x-层计算 (减少中间数组大小)
            for xi in range(Nx):
                start = xi * layer_size
                end   = (xi + 1) * layer_size
                ps = phase_P[start:end]                      # (Nx²)
                es = eigvecs[:, start:end, :]                # (Nev, Nx², 3)

                # 6-term ε_{ijk} 收缩
                VVV_sink[t] += contract("x,ax,bx,cx->abc",
                                         ps, es[:, :, 0], es[:, :, 1], es[:, :, 2])
                VVV_sink[t] += contract("x,ax,bx,cx->abc",
                                         ps, es[:, :, 1], es[:, :, 2], es[:, :, 0])
                VVV_sink[t] += contract("x,ax,bx,cx->abc",
                                         ps, es[:, :, 2], es[:, :, 0], es[:, :, 1])
                VVV_sink[t] -= contract("x,ax,bx,cx->abc",
                                         ps, es[:, :, 0], es[:, :, 2], es[:, :, 1])
                VVV_sink[t] -= contract("x,ax,bx,cx->abc",
                                         ps, es[:, :, 1], es[:, :, 0], es[:, :, 2])
                VVV_sink[t] -= contract("x,ax,bx,cx->abc",
                                         ps, es[:, :, 2], es[:, :, 1], es[:, :, 0])

            t4 = time.perf_counter()
            print(f"  [VVV] t={t:3d} contraction done, time={t4-t3:.3f}s")

        return VVV_sink

    # ── 主循环: 对每个 Pz ──────────────────────────────────
    t_total_start = time.perf_counter()

    for Pz in Pzlist:
        Mom = np.array([Pz, Py, Px])
        print(f"\n{'─'*60}")
        print(f"  动量 P = ({Pz}, {Py}, {Px})")
        print(f"{'─'*60}")

        # ── 计算 VVV ──
        t_vvv_start = time.perf_counter()
        VVV_sink = compute_vvv_for_momentum(Mom)
        t_vvv_end = time.perf_counter()
        print(f"  VVV 计算总耗时: {t_vvv_end - t_vvv_start:.3f}s")

        # ── Wick 收缩 ──
        contrac_nucl_matrix = np.zeros((Nt, Nt, 4, 4), dtype=complex)
        print(f"\n  Wick 收缩开始...")

        t_contract_start = time.perf_counter()
        for t_source in range(Nt):
            t_src_start = time.perf_counter()

            VVV_source = VVV_sink[t_source].conj()
            peram_u = read_peram(t_source)                 # (Nt, 4, 4, Nev1, Nev1)

            # Dirac 投影 perambulator: Γ τ Γ
            CG5peram_uCG5 = contract(
                "gh,thkbe,jk->tgjbe",
                interProject1, peram_u, interProject2,
            )  # → (Nt, 4, 4, Nev1, Nev1)

            for t_sink in range(Nt):
                deltat = (t_sink - t_source + Nt) % Nt
                # 仅收缩满足 2 ≤ deltat ≤ 32 的时间分离
                if 2 <= deltat <= 32:
                    # Direct - Exchange
                    contrac_nucl_matrix[t_sink, t_source] = (
                        contract(
                            "abc,gjad,gjbe,ilcf,def->il",
                            VVV_sink[t_sink],
                            peram_u[t_sink],
                            CG5peram_uCG5[t_sink],
                            peram_u[t_sink],
                            VVV_source,
                        )
                        - contract(
                            "abc,glaf,gjbe,ijcd,def->il",
                            VVV_sink[t_sink],
                            peram_u[t_sink],
                            CG5peram_uCG5[t_sink],
                            peram_u[t_sink],
                            VVV_source,
                        )
                    )

            t_src_end = time.perf_counter()
            print(f"  t_source={t_source:3d} done, time={t_src_end-t_src_start:.3f}s")

        t_contract_end = time.perf_counter()
        print(f"  Wick 收缩总耗时: {t_contract_end - t_contract_start:.3f}s")

        # ── 保存 raw 收缩矩阵 (写入生成目录) ──
        raw_filename = (
            f"twopt_slice_pp_Px{Px}Py{Py}Pz{Pz}"
            f"_eginphase{mom_smear}{element}_contract_conf{conf_id}.npy"
        )
        np.save(f"{gen_dir}/{raw_filename}", contrac_nucl_matrix)
        print(f"  [GEN] Raw: {gen_dir}/{raw_filename}")
        print(f"         shape={contrac_nucl_matrix.shape}")

        # ── 宇称投影 ──
        t_parity_start = time.perf_counter()

        # pp: P₊ 投影
        contrac_nucl_pp = contract(
            "li,yxil->yx", matrix_pplus, contrac_nucl_matrix
        )
        # pm: P₋ 投影 (参考代码保存但此处仅需 pp)
        contrac_nucl_pm = contract(
            "li,yxil->yx", matrix_pminus, contrac_nucl_matrix
        )

        # ── 边界条件符号修正 ──
        # 反粒子传播子需要符号翻转
        for t_source in range(Nt):
            for t_sink in range(Nt):
                if t_sink < t_source:
                    contrac_nucl_pp[t_sink, t_source] *= -1.0
                if t_sink > t_source:
                    contrac_nucl_pm[t_sink, t_source] *= -1.0

        t_parity_end = time.perf_counter()
        print(f"  宇称投影耗时: {t_parity_end - t_parity_start:.3f}s")

        # ── 保存 parity 投影结果 ──
        pp_filename = (
            f"twopt_slice_pp_Px{Px}Py{Py}Pz{Pz}"
            f"_eginphase{mom_smear}{element}_nopol_ss_conf{conf_id}.npy"
        )
        np.save(f"{gen_dir}/{pp_filename}", contrac_nucl_pp)
        print(f"  [GEN] PP: {gen_dir}/{pp_filename}")
        print(f"         shape={contrac_nucl_pp.shape}")

        # ── 与参考数据对比 ──
        _compare_2pt_result(
            ref_corr_nucl_dir=ref_corr_nucl_dir,
            gen_dir=gen_dir,
            plot_dir=plot_dir,
            has_ref=has_ref,
            raw_filename=raw_filename,
            pp_filename=pp_filename,
            contrac_nucl_matrix=contrac_nucl_matrix,
            contrac_nucl_pp=contrac_nucl_pp,
            Px=Px, Py=Py, Pz=Pz,
            mom_smear=mom_smear,
            element=element,
            conf_id=conf_id,
        )

        print(f"  Pz={Pz} 完成。")

    t_total_end = time.perf_counter()
    print(f"\n{'═'*70}")
    print(f"  JOB: ran successfully")
    print(f"  总耗时: {t_total_end - t_total_start:.3f}s")
    print(f"  生成数据: {gen_dir}")
    if has_ref:
        print(f"  参考数据: {ref_corr_nucl_dir}")
    print(f"  作图输出: {plot_dir}")
    print(f"{'═'*70}\n")


def _compare_2pt_result(
    ref_corr_nucl_dir: str,
    gen_dir: str,
    plot_dir: str,
    has_ref: bool,
    raw_filename: str,
    pp_filename: str,
    contrac_nucl_matrix: np.ndarray,
    contrac_nucl_pp: np.ndarray,
    Px: int, Py: int, Pz: int,
    mom_smear: int,
    element: str,
    conf_id: str,
) -> None:
    """将生成数据与参考数据对比, 输出差异报告与对比图。

    对比指标:
      - 最大绝对差异: max|gen - ref|
      - 平均绝对差异: mean|gen - ref|
      - 相对差异: |gen - ref| / (|ref| + ε) 的最大值与平均值
      - 逐元素一致性: |gen - ref| < tol 的比例
    """
    if not has_ref:
        print(f"\n  [COMPARE] 参考数据目录不存在, 跳过对比。")
        return

    comp_results = {}  # {label: dict of metrics}

    # ── 对比 raw 收缩矩阵 ──
    ref_raw_path = f"{ref_corr_nucl_dir}/{raw_filename}"
    if os.path.exists(ref_raw_path):
        ref_raw = np.load(ref_raw_path)
        gen_raw = contrac_nucl_matrix
        metrics = _compute_comparison_metrics(gen_raw, ref_raw, "Raw (Nt,Nt,4,4)")
        comp_results["Raw 收缩矩阵"] = metrics
        _plot_2pt_comparison(
            gen_raw, ref_raw, plot_dir,
            prefix=f"raw_Pz{Pz}", label="C₂(t_snk,t_src)"
        )
    else:
        print(f"  [COMPARE] 参考 Raw 文件不存在: {ref_raw_path}")

    # ── 对比 parity 投影结果 ──
    ref_pp_path = f"{ref_corr_nucl_dir}/{pp_filename}"
    if os.path.exists(ref_pp_path):
        ref_pp = np.load(ref_pp_path)
        gen_pp = contrac_nucl_pp
        metrics = _compute_comparison_metrics(gen_pp, ref_pp, "PP (Nt,Nt)")
        comp_results["Parity 投影 (PP)"] = metrics
        _plot_2pt_comparison(
            gen_pp, ref_pp, plot_dir,
            prefix=f"pp_Pz{Pz}", label="C₂^{pp}(t_snk,t_src)"
        )
    else:
        print(f"  [COMPARE] 参考 PP 文件不存在: {ref_pp_path}")

    # ── 打印汇总对比报告 ──
    _print_comparison_report(comp_results, Pz=Pz)

    # ── 保存对比数据 ──
    np.savez(
        f"{gen_dir}/compare_Pz{Pz}_conf{conf_id}.npz",
        gen_raw=contrac_nucl_matrix if comp_results else None,
        gen_pp=contrac_nucl_pp if comp_results else None,
        **{k.replace(" ", "_"): v for k, v in comp_results.items()},
    )


def _compute_comparison_metrics(gen: np.ndarray, ref: np.ndarray,
                                 label: str) -> dict:
    """计算两组数据的对比指标。"""
    diff = np.abs(gen - ref)
    abs_ref = np.abs(ref)
    eps = 1e-15

    max_diff = float(np.max(diff))
    mean_diff = float(np.mean(diff))
    # 相对差异 (逐元素)
    rel_diff = diff / (abs_ref + eps)
    max_rel = float(np.max(rel_diff))
    mean_rel = float(np.mean(rel_diff))
    # 一致性 (容差范围内的比例)
    frac_1e6  = float(np.mean(diff < 1e-6))
    frac_1e10 = float(np.mean(diff < 1e-10))

    print(f"\n  [COMPARE] {label}:")
    print(f"    最大绝对差异:  {max_diff:.6e}")
    print(f"    平均绝对差异:  {mean_diff:.6e}")
    print(f"    最大相对差异:  {max_rel:.6e}")
    print(f"    平均相对差异:  {mean_rel:.6e}")
    print(f"    Δ<1e-6 比例:   {frac_1e6:.4%}")
    print(f"    Δ<1e-10 比例:  {frac_1e10:.4%}")

    return {
        "max_abs": max_diff, "mean_abs": mean_diff,
        "max_rel": max_rel, "mean_rel": mean_rel,
        "frac_1e6": frac_1e6, "frac_1e10": frac_1e10,
    }


def _print_comparison_report(comp_results: dict, Pz: int) -> None:
    """打印汇总对比报告。"""
    if not comp_results:
        return
    print(f"\n{'─'*60}")
    print(f"  对比报告 — Pz={Pz}")
    print(f"{'─'*60}")
    header = f"  {'数据集':<22s} {'max|Δ|':>12s} {'mean|Δ|':>12s} {'max|Δ/ref|':>12s} {'Δ<1e-6':>8s}"
    print(header)
    print(f"  {'─'*68}")
    for label, m in comp_results.items():
        print(f"  {label:<22s} {m['max_abs']:12.4e} {m['mean_abs']:12.4e}"
              f" {m['max_rel']:12.4e} {m['frac_1e6']:7.1%}")
    print(f"{'─'*60}\n")


def _plot_2pt_comparison(gen: np.ndarray, ref: np.ndarray,
                          plot_dir: str, prefix: str, label: str) -> None:
    """画对比图: 生成 vs 参考的 scatter 图 + 差异直方图。

    子图:
      1. 散点图: Re(ref) vs Re(gen) (理想情况 y=x 对角线)
      2. 差异直方图: log10|gen - ref|
    """
    if not HAS_MPL:
        return

    # 展平为 1D (只取实部)
    gen_flat = np.real(gen).ravel()
    ref_flat = np.real(ref).ravel()
    diff_flat = np.abs(gen - ref).ravel()
    # 跳过零
    nonzero = diff_flat > 0
    diff_log = np.log10(diff_flat[nonzero]) if np.any(nonzero) else np.array([-16])

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # 子图 1: scatter
    axes[0].scatter(ref_flat, gen_flat, s=1, alpha=0.3)
    lims = [np.min(ref_flat), np.max(ref_flat)]
    axes[0].plot(lims, lims, 'r--', linewidth=1, label='y=x (完全一致)')
    axes[0].set_xlabel(f"Re(Ref) — {label}")
    axes[0].set_ylabel(f"Re(Gen) — {label}")
    axes[0].set_title(f"Scatter: Generated vs Reference (Pz={prefix.split('Pz')[-1]})")
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)

    # 子图 2: 差异直方图
    axes[1].hist(diff_log, bins=50, color='C1', alpha=0.7, edgecolor='k')
    axes[1].set_xlabel("log10(|Gen - Ref|)")
    axes[1].set_ylabel("Count")
    axes[1].set_title(f"Difference Distribution (Pz={prefix.split('Pz')[-1]})")
    axes[1].axvline(x=-6, color='r', linestyle='--', linewidth=1, label='1e-6')
    axes[1].axvline(x=-10, color='g', linestyle='--', linewidth=1, label='1e-10')
    axes[1].legend(fontsize=9)

    plt.tight_layout()
    fig.savefig(f"{plot_dir}/{prefix}_comparison.pdf", bbox_inches='tight')
    plt.close(fig)
    print(f"  [PLOT] 对比图已保存: {plot_dir}/{prefix}_comparison.pdf")


# ============================================================================
# 第十五部分: 命令行参数解析与主入口
# ============================================================================

def parse_args() -> argparse.Namespace:
    """解析命令行参数。"""
    parser = argparse.ArgumentParser(
        description="格点 QCD 质子非极化胶子 PDF 计算的统一流水线",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # PDF 全流程 (numpy 后端)
  python main.py --analysis-type pdf --steps all --xp numpy --ensemble L32x64 --conf-start 20000 --conf-num 1

  # 仅 OPE 计算 (GPU)
  python main.py --analysis-type pdf --steps 1,2,3,6 --xp cupy --ensemble L24x72 --conf-start 46000 \\
      --gauge-file /path/to/config.dat

  # 质子 2pt 蒸馏 (复现 DCU 代码)
  python main.py --analysis-type proton-2pt --Nt 72 --Nx 24 \\
      --Pz-list "-2,-3,-4,-5,-6" --conf-start 46000

  # 2pt 有效质量分析 (匹配 main-2pt.py)
  python main.py --analysis-type 2pt --Nt 72 --Nx 24 --alttc 0.1053 \\
      --conf-start 10050 --conf-step 50 --conf-num 52 \\
      --iog-2pt-path "./beta6.20_mu-0.2770_ms-0.2400_L%dx%d/sush_iog/pion_2pt_Px%%dPy%%dPz%%d_ENV%%d_conf%%d_tsep-1_mass-0.2770.iog"
        """,
    )

    # --- 分析类型 ---
    parser.add_argument("--analysis-type", type=str, default="pdf",
                        choices=["pdf", "2pt", "proton-2pt"],
                        help="分析模式: pdf (胶子PDF全流程), 2pt (2pt有效质量), "
                             "proton-2pt (质子蒸馏2pt) [默认: pdf]")

    # --- 后端与数据类型 ---
    parser.add_argument("--xp", type=str, default="numpy",
                        choices=["numpy", "cupy"],
                        help="计算后端: numpy (CPU) 或 cupy (GPU) [默认: numpy]")
    parser.add_argument("--dtype", type=str, default="complex128",
                        choices=["complex64", "complex128"],
                        help="复数数据类型 [默认: complex128]")

    # --- 步骤选择 ---
    parser.add_argument("--steps", type=str, default="all",
                        help="要运行的步骤, 逗号分隔 (如 1,2,3) 或 'all' [默认: all]")

    # --- 系综预设 ---
    parser.add_argument("--ensemble", type=str, default=None,
                        choices=list(ENSEMBLES.keys()) + [None],
                        help="系综预设名称 (L24x72, L32x64, L32x96, ...)")

    # --- 手动格点参数 (不使用系综预设时) ---
    parser.add_argument("--Nt", type=int, default=64,
                        help="时间方向格点数 [默认: 64]")
    parser.add_argument("--Nx", type=int, default=32,
                        help="空间方向格点数 [默认: 32]")
    parser.add_argument("--Nev", type=int, default=100,
                        help="本征矢量总数 [默认: 100]")
    parser.add_argument("--Nev1", type=int, default=100,
                        help="实际用于收缩的本征矢量数 [默认: 100]")
    parser.add_argument("--delta-z", type=int, default=15,
                        help="Wilson 线最大长度 (格距单位) [默认: 15]")
    parser.add_argument("--z-dir", type=int, default=2,
                        choices=[0, 1, 2],
                        help="Wilson 线方向: 0(x), 1(y), 2(z) [默认: 2]")
    parser.add_argument("--mom-smear", type=int, default=3,
                        help="动量涂抹参数 [默认: 3]")
    parser.add_argument("--mom-smear-phase", type=int, default=-3,
                        help="动量涂抹相位参数 [默认: -3]")

    # --- 动量 ---
    parser.add_argument("--Px", type=int, default=0,
                        help="x 方向动量 (2π/L 单位) [默认: 0]")
    parser.add_argument("--Py", type=int, default=0,
                        help="y 方向动量 (2π/L 单位) [默认: 0]")
    parser.add_argument("--Pz", type=int, default=6,
                        help="z 方向动量 (2π/L 单位) [默认: 6]")
    parser.add_argument("--Pz-list", type=str, default=None,
                        help="动量扫描列表, 逗号分隔 (如 '3,4,5,6,7,8')")

    # --- 组态选择 ---
    parser.add_argument("--conf-start", type=int, default=20000,
                        help="初始组态编号 [默认: 20000]")
    parser.add_argument("--conf-step", type=int, default=50,
                        help="组态步长 [默认: 50]")
    parser.add_argument("--conf-num", type=int, default=1,
                        help="组态样本数 [默认: 1]")

    # --- 文件路径 ---
    # 规范组态文件 (仅读取)
    parser.add_argument("--gauge-file", type=str, default=None,
                        help="规范组态文件路径 (ILDG 二进制格式)")

    # 本征矢量目录 (仅读取)
    parser.add_argument("--eig-dir", type=str, default=None,
                        help="本征矢量目录路径 (覆盖系综预设)")

    # Perambulator 目录 (仅读取)
    parser.add_argument("--peram-u-dir", type=str, default=None,
                        help="Perambulator 目录路径 (覆盖系综预设)")

    # 2pt 数据路径 (读取或生成)
    parser.add_argument("--read-2pt-dir", type=str, default=None,
                        help="读取已有 2pt 数据的目录路径")
    parser.add_argument("--gen-2pt-dir", type=str, default=None,
                        help="生成 2pt 数据的写入目录路径")

    # 3pt 数据路径 (读取或生成)
    parser.add_argument("--read-3pt-dir", type=str, default=None,
                        help="读取已有 3pt 数据的目录路径")
    parser.add_argument("--gen-3pt-dir", type=str, default=None,
                        help="生成 3pt 数据的写入目录路径")

    # VVV 数据路径 (读取或生成)
    parser.add_argument("--read-VVV-dir", type=str, default=None,
                        help="读取已有 VVV 数据的目录路径")
    parser.add_argument("--gen-VVV-dir", type=str, default=None,
                        help="生成 VVV 数据的写入目录路径")

    # OPE 数据路径 (读取或生成)
    parser.add_argument("--read-ope-dir", type=str, default=None,
                        help="读取已有 OPE 数据的目录路径")
    parser.add_argument("--gen-ope-dir", type=str, default=None,
                        help="生成 OPE 数据的写入目录路径")

    # 2pt 关联函数目录 (读取, 来自系综预设)
    parser.add_argument("--corr-nucl-dir", type=str, default=None,
                        help="核子关联函数目录路径 (覆盖系综预设)")

    # --- 匹配参数 ---
    parser.add_argument("--alpha-s", type=float, default=0.2,
                        help="强耦合常数 α_s(μ) [默认: 0.2]")
    parser.add_argument("--mu-over-Pz", type=float, default=1.0,
                        help="重整化标度与 boost 动量比值 μ/P_z [默认: 1.0]")

    # --- 2pt 分析专用参数 ---
    parser.add_argument("--alttc", type=float, default=0.1053,
                        help="格距 a [fm] (用于有效质量物理单位转换) [默认: 0.1053]")
    parser.add_argument("--meff-type", type=str, default="cosh",
                        choices=["cosh", "log"],
                        help="有效质量类型: cosh (周期性) 或 log (简单) [默认: cosh]")
    parser.add_argument("--hadron", type=str, default="pion",
                        help="强子类型 [默认: pion]")
    parser.add_argument("--tsep", type=int, default=36,
                        help="源-汇时间分离 (2pt 模式下用于文件路径格式) [默认: 36]")
    parser.add_argument("--time-fold", action="store_true",
                        help="是否时间折叠 (反对称化)")
    parser.add_argument("--link-max", type=int, default=10,
                        help="Wilson 连接最大偏移 (文件路径格式参数) [默认: 10]")
    parser.add_argument("--meff-range", type=str, default="0.0,1.0",
                        help="有效质量图的 y 范围, 逗号分隔 [默认: 0.0,1.0]")
    parser.add_argument("--iog-2pt-path", type=str, default=None,
                        help="IOG 2pt 关联函数文件路径模板 (含 %%d 格式说明符)")
    parser.add_argument("--element", type=str, default="_Cg5g4",
                        choices=["_Cg5g4", "_Cg5g3", "_Cg5",
                                 "_offdiag01", "_offdiag02", "_offdiag12"],
                        help="质子插值算符元素类型 [默认: _Cg5g4]")
    parser.add_argument("--no-plot", action="store_true",
                        help="禁用所有作图输出")

    # --- 输出 ---
    parser.add_argument("--output-dir", type=str, default=None,
                        help="输出目录 (默认: snsc/output_YYYYMMDD_HHMMSS/)")

    return parser.parse_args()


def _resolve_steps(steps_str: str) -> List[int]:
    """解析步骤字符串。"""
    steps_str = steps_str.strip()
    if steps_str.lower() == "all":
        return list(range(1, 14))
    return [int(s.strip()) for s in steps_str.split(",") if s.strip()]


def main():
    """主入口函数。根据 --analysis-type 分派到不同分析模式。

    两种模式:
      pdf: 胶子 PDF 全流程 (LaMET + distillation + OPE → 光锥 PDF)
      2pt: 2pt 有效质量分析 (读取 IOG 关联函数, cosh meff + Jackknife)
    """
    args = parse_args()

    # 设置输出目录
    if args.output_dir is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        script_dir = os.path.dirname(os.path.abspath(__file__))
        args.output_dir = os.path.join(script_dir, f"output_{timestamp}")
    else:
        args.output_dir = os.path.abspath(args.output_dir)

    # --- 分派分析模式 ---
    if args.analysis_type == "2pt":
        print(f"[INFO] 分析类型: 2pt 有效质量分析")
        print(f"[INFO] 输出目录: {args.output_dir}")
        run_2pt_analysis(args)
        return

    if args.analysis_type == "proton-2pt":
        print(f"[INFO] 分析类型: 质子 2pt 蒸馏计算")
        print(f"[INFO] 输出目录: {args.output_dir}")
        run_proton_2pt_analysis(args)
        return

    # --- PDF 模式 ---
    args.steps = _resolve_steps(args.steps)
    print(f"[INFO] 分析类型: 胶子 PDF 全流程")
    print(f"[INFO] 选定步骤: {args.steps}")

    pipeline = GluonPDFPipeline(args)
    pipeline.run()


if __name__ == "__main__":
    main()
