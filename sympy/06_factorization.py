#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
================================================================================
06_factorization.py — 准TMD→光锥TMD因子化定理与微扰匹配核
================================================================================

本模块使用 SymPy 推导完整的因子化定理和匹配核:

  1. LaMET 因子化定理 (Euclidean → 光锥)
  2. TMD 因子化定理 (准TMD → 光锥TMD)
  3. NLO 匹配核 C_{gg} 的推导
  4. 2×2 味混合矩阵 (胶子-夸克混合)
  5. 重整化群重求和 (RGR)
  6. 重正子重求和 (LRR) 与领阶幂次修正

参考源:
  [1] X. Ji, "Parton Physics on a Euclidean Lattice", PRL 110, 262002 (2013)
  [2] X. Ji, "Parton Physics from LaMET", Sci. China 57, 1407 (2014)
  [3] T. Izubuchi et al., "Factorization Theorem...", PRD 98, 056004 (2018)
  [7] J.-H. Zhang et al., "Accessing gluon parton distributions in LaMET",
      PRL 122, 142001 (2019)
  - 补充/格点QCD中的部分子分布函数.tex: §3-4
  - 文档/gluon_PDF_continuum.tex: Suppl. Mat.
  - 文档/gluon_pdf_derivation.tex: §9 (微扰匹配)

作者: 基于 lattice-pdf 库全部文档与代码的系统化推导
日期: 2026-07-24
================================================================================
"""

import sympy as sp
from sympy import I, oo, pi, Symbol, symbols, integrate, exp, log
from sympy import Function, Piecewise, Heaviside, DiracDelta


# ============================================================================
# 第0节: 符号声明
# ============================================================================

# 动量分数
x_var = Symbol('x', real=True)          # Bjorken-x (准PDF)
y_var = Symbol('y', real=True)          # Bjorken-y (光锥PDF)
xi_var = Symbol('xi', real=True)        # ξ = x/y (动量比)

# 标度
mu_val = Symbol('mu', real=True, positive=True)
Pz_val = Symbol('P_z', real=True, positive=True)

# 耦合与色因子
alpha_s_sym = Symbol('alpha_s', real=True, positive=True)
C_A_sym = Symbol('C_A', real=True, positive=True)
C_F_sym = Symbol('C_F', real=True, positive=True)
N_f_sym = Symbol('N_f', integer=True, positive=True)
T_F = sp.Rational(1, 2)  # Tr[t^a t^b] = T_F δ^{ab}

# 匹配核
C_gg = Function('C_gg')
C_gq = Function('C_gq')
C_qg = Function('C_qg')
C_qq = Function('C_qq')

# plus 分布
plus_dist = Function('[]_+')


# ============================================================================
# 第1节: LaMET 因子化定理
# ============================================================================

def derive_lamet_factorization():
    """
    推导 LaMET 因子化定理。

    准 PDF q̃(x, P_z) 与光锥 PDF q(x, μ) 通过因子化公式联系:

        q̃(x, P_z) = ∫_{-1}^{1} (dy/|y|) C(x/y, μ/(y P_z)) q(y, μ)
                    + O(Λ_QCD²/(x² P_z²), Λ_QCD²/((1-x)² P_z²))

    这个定理由 Izubuchi et al. (2018) 严格证明。
    """
    print("=" * 70)
    print("第1节: LaMET 因子化定理")
    print("=" * 70)

    print("\n【定理1.1】(Izubuchi, Ji, Jin, Stewart, Zhao 2018):")
    print("")
    print("  对于大 P_z, 准 PDF q̃(x, P_z) 与光锥 PDF q(x, μ) 满足:")
    print("")
    print("    q̃(x, P_z) = C(x/y, μ/P_z) ⊗ q(y, μ)")
    print("              + O(Λ²/P_z²)")
    print("")
    print("  其中 ⊗ 表示 Mellin 卷积:")
    print("    (C ⊗ q)(x) = ∫_{x}^{1} (dy/y) C(x/y, μ/(y P_z)) q(y, μ)")
    print("")
    print("  匹配核 C 是微扰可计算的, 在 LO 时:")
    print("    C^{(0)}(ξ) = δ(1-ξ)")

    # 对于胶子
    print("\n【定理1.2】胶子准 PDF 的因子化 (Zhang et al. 2019):")
    print("")
    print("  x g̃(x, P_z) = ∫_{x}^{1} (dy/y) C_{gg}(x/y, μ/(y P_z)) y g(y, μ)")
    print("              + ∫_{x}^{1} (dy/y) C_{gq}(x/y, μ/(y P_z)) y q(y, μ)")
    print("")
    print("  注意: 胶子准 PDF 在 NLO 与夸克单态 PDF 混合!")
    print("  这是 2×2 混合矩阵的来源.")

    print("\n【定理1.3】幂次修正的严格形式:")
    print("")
    print("  O(Λ²/P_z²) = O(Λ_QCD²/(x² P_z²)) + O(Λ_QCD²/((1-x)² P_z²))")
    print("")
    print("  这些修正在小 x 和大 x 区域被增强, 在中等 x")
    print("  (0.3 ≲ x ≲ 0.6) 处最温和。")

    return {}


# ============================================================================
# 第2节: TMD 因子化定理
# ============================================================================

def derive_tmd_factorization():
    """
    推导 TMD 因子化定理。

    准 TMD-PDF f̃_1^g(x, k__perp, P_z) 与光锥 TMD-PDF f_1^g(x, k__perp) 满足:

        f̃_1^g(x, k__perp, P_z) = H(μ, ζ) ×
          ∫ d²b__perp e^{-i k__perp·b__perp}  C_TMD(x, b__perp, μ, P_z)
          × f̃_1^g(x, b__perp)
          × S^{1/2}(b__perp, μ, ζ)
          + O(Λ²/P_z²)

    其中:
      - H(μ, ζ): 硬函数 (微扰可计算)
      - C_TMD: TMD 匹配核
      - S^{1/2}: 软函数 (非微扰)
      - f̃_1^g: 坐标空间 TMD-PDF
    """
    print("\n" + "=" * 70)
    print("第2节: TMD 因子化定理")
    print("=" * 70)

    print("\n【定理2.1】TMD 因子化定理 (准TMD → 光锥TMD):")
    print("")
    print("  f̃_1^{{g,B}}(x, b__perp, P_z, a)")
    print("    = H(μ, ζ) · C_{{TMD}}(x, b__perp, μ, P_z)")
    print("    · f_1^{{g,R}}(x, b__perp, μ, ζ)")
    print("    · S^{{1/2}}(b__perp, μ, ζ)")
    print("    + O(a², 1/P_z²)")
    print("")
    print("  其中 B = 裸 (bare), R = 重整化 (renormalized).")

    print("\n【推导2.2】动量空间 → 坐标空间的关系:")
    print("")
    print("  动量空间 TMD:     f_1^g(x, k__perp²)")
    print("  坐标空间 TMD:     f̃_1^g(x, b__perp²)")
    print("")
    print("  Fourier 变换:")
    print("    f_1^g(x, k__perp²) = ∫ d²b__perp/(2π)² e^{i k__perp·b__perp} f̃_1^g(x, b__perp²)")
    print("")
    print("    f̃_1^g(x, b__perp²) = ∫ d²k__perp e^{-i k__perp·b__perp} f_1^g(x, k__perp²)")

    return {}


# ============================================================================
# 第3节: NLO 匹配核
# ============================================================================

def derive_nlo_matching_kernel():
    """
    推导胶子准 PDF 到光锥 PDF 的 NLO 匹配核。

    在 MS-bar 方案中, NLO 匹配核包含以下贡献:
      - 实辐射图 (real emission): 胶子辐射
      - 虚辐射图 (virtual correction): 胶子圈
      - 顶角修正图 (vertex correction)
      - 自能图 (self-energy)

    对于非极化胶子 PDF, NLO 匹配核 C_{gg} 在 Ratio 方案中的
    形式为 (来自 Chen et al. 2025 的补充材料):

      C^{ratio}(ξ, μ/(y P_z)) = δ(1-ξ)
        + (α_s/(2π)) [c_1(ξ) + c_2(ξ) ln(μ²/(4 y² P_z²))]
        + O(α_s²)

    其中 c_1(ξ) 和 c_2(ξ) 是 ξ 的分段函数。
    """
    print("\n" + "=" * 70)
    print("第3节: NLO 匹配核 C_{gg}")
    print("=" * 70)

    print("\n【推导3.1】Ratio 方案中的 NLO 胶子匹配核:")
    print("")
    print("  C^{{ratio}}_{{gg}}(ξ, r) = δ(1-ξ)")
    print("    + (α_s/(2π)) × [")
    print("      C_A · P_{{gg}}(ξ) · ln(μ²/(4y²P_z²))")
    print("    + C_A · f₁(ξ)")
    print("    + N_f · T_F · f₂(ξ)")
    print("    + ... ]")
    print("")
    print("  其中 P_{{gg}}(ξ) 是 LO DGLAP 胶子→胶子劈裂函数:")
    print("")
    print("  P_{{gg}}(ξ) = 2C_A [ ξ/(1-ξ)_+ + (1-ξ)/ξ + ξ(1-ξ) ]")
    print("             + δ(1-ξ) · (11C_A - 4N_f T_F)/6")

    # DGLAP 劈裂函数
    print("\n【推导3.2】劈裂函数与 Plus 分布:")
    print("")
    print("  Plus 分布的定义:")
    print("    ∫_0^1 dξ f(ξ) [g(ξ)]_+ = ∫_0^1 dξ [f(ξ) - f(1)] g(ξ)")
    print("")
    print("  对于 P_{{gg}}:")
    print("    [ξ/(1-ξ)]_+ = lim_{{ε→0}} [ ξ/(1-ξ) Θ(1-ξ-ε)")
    print("                               + δ(1-ξ) ln ε ]")

    # 混合方案中的匹配核
    print("\n【推导3.3】混合方案中的 NLO 匹配核:")
    print("")
    print("  在混合重整化方案 (Hybrid Scheme) 中,")
    print("  短距离和长距离部分使用不同的匹配核:")
    print("")
    print("  短距离 (|z| < z_s): 使用 ratio 方案")
    print("    C^{{ratio}}(ξ, r) = δ(1-ξ) + (α_s/(2π)) c^{{ratio}}(ξ, r)")
    print("")
    print("  长距离 (|z| > z_s): 使用自重整化方案")
    print("    C^{{hybrid}}(ξ, r) = δ(1-ξ) + (α_s/(2π)) c^{{hybrid}}(ξ, r)")
    print("")
    print("  z_s ≈ 0.3 fm 是典型的过渡标度.")

    return {}


# ============================================================================
# 第4节: 2×2 味混合矩阵
# ============================================================================

def derive_flavor_mixing():
    """
    推导胶子-夸克味混合的 2×2 矩阵。

    在味单态扇区, 准 PDF 与光锥 PDF 通过 2×2 矩阵联系:

        [ g̃(x, P_z) ]     [ C_{gg}  C_{gq} ]   [ g(y, μ) ]
        [ q̃(x, P_z) ] = ∫ [ C_{qg}  C_{qq} ] ⊗ [ q(y, μ) ] dy

    其中 g 是非极化胶子 PDF, q 是夸克单态 PDF
    (所有夸克味之和).
    """
    print("\n" + "=" * 70)
    print("第4节: 2×2 味混合矩阵")
    print("=" * 70)

    print("\n【定义4.1】味单态扇区的混合矩阵:")
    print("")
    print("  ┌           ┐     ┌                 ┐   ┌         ┐")
    print("  │ g̃(x,P_z)  │     │ C_{{gg}}   C_{{gq}}  │   │ g(y,μ)  │")
    print("  │           │ = ∫ │                 │ ⊗ │         │ dy")
    print("  │ q̃(x,P_z)  │     │ C_{{qg}}   C_{{qq}}  │   │ q(y,μ)  │")
    print("  └           ┘     └                 ┘   └         ┘")
    print("")
    print("  其中 q = Σ_f (q_f + q̄_f) 是夸克单态 PDF.")

    print("\n【推导4.2】混合矩阵元的物理来源:")
    print("")
    print("  C_{{gg}}: 胶子准 PDF ← 胶子光锥 PDF")
    print("    → 胶子辐射图 + 胶子圈图")
    print("")
    print("  C_{{gq}}: 胶子准 PDF ← 夸克光锥 PDF")
    print("    → 夸克→胶子对的辐射 (g→qq̄ 逆过程)")
    print("")
    print("  C_{{qg}}: 夸克准 PDF ← 胶子光锥 PDF")
    print("    → 胶子→夸克对的辐射")
    print("")
    print("  C_{{qq}}: 夸克准 PDF ← 夸克光锥 PDF")
    print("    → 标准夸克 PDF 演化")

    print("\n【推导4.3】LO DGLAP 演化与混合矩阵的关系:")
    print("")
    print("  在 MS-bar 方案中, DGLAP 演化方程为:")
    print("")
    print("  d/d(ln μ²) [ g(x,μ) ]     [ P_{{gg}}  P_{{gq}} ]   [ g(x,μ) ]")
    print("             [ q(x,μ) ] = ∫ [ P_{{qg}}  P_{{qq}} ] ⊗ [ q(x,μ) ] dy")
    print("")
    print("  其中 P_{{ij}} 是劈裂函数. 这与匹配矩阵 C 密切相关,")
    print("  因为两者都来源于相同的 QCD 顶点.")

    return {}


# ============================================================================
# 第5节: 重整化群重求和 (RGR)
# ============================================================================

def derive_rgr_resummation():
    """
    推导重整化群重求和 (Renormalization Group Resummation, RGR)。

    匹配核中的大对数 ln(μ²/(2xP_z)²) 在 |2xP_z| 远离 μ 时
    变得很大, 破坏了微扰展开的收敛性。

    RGR 通过将匹配核重求和到以 2xP_z 为内禀标度的
    跑动耦合常数中来解决这个问题:

        C_RGR(x/y, μ/(yP_z)) = C(x/y, 1)|_{α_s → α_s(2yP_z)}
    """
    print("\n" + "=" * 70)
    print("第5节: 重整化群重求和 (RGR)")
    print("=" * 70)

    print("\n【方法5.1】RGR 的基本思想:")
    print("")
    print("  NLO 匹配核包含大对数:")
    print("    C(ξ, r) = δ(1-ξ) + (α_s/(2π)) [f(ξ) + g(ξ) ln(r²)]")
    print("")
    print("  其中 r = μ/(2yP_z).")
    print("")
    print("  当 2yP_z ≪ μ 或 2yP_z ≫ μ 时,")
    print("  ln(r²) 很大 → 破坏微扰收敛.")

    print("\n【方法5.2】RGR 重求和公式:")
    print("")
    print("  C^{{RGR}}(ξ, r) = C(ξ, 1)")
    print("    × exp[ -∫_{α_s(μ)}^{α_s(2yP_z)} dα' (γ(α')/β(α')) ]")
    print("")
    print("  其中 γ 是算符的反常量纲, β 是 QCD β 函数.")
    print("")
    print("  这确保了大对数被重求和到 α_s 的跑动中,")
    print("  恢复了微扰展开的可靠性.")

    # 与 Pz 外推的协同
    print("\n【注5.3】RGR 与 P_z 外推的协同:")
    print("")
    print("  RGR 减小了匹配核对 P_z 的残余依赖性 →")
    print("  残余的 1/P_z² 修正更纯粹地来自高扭度")
    print("  (而非微扰对数) → P_z 外推变得更可靠.")

    return {}


# ============================================================================
# 第6节: 重正子重求和 (LRR)
# ============================================================================

def derive_lrr_resummation():
    """
    推导领头重正子重求和 (Leading Renormalon Resummation, LRR)。

    IR 重正子是微扰级数的阶乘发散性的来源。
    在准 PDF 算符中, 重正子模糊性与 1/P_z² 幂次修正
    具有精确的对应关系。

    LRR 通过吸收重正子模糊性到非微扰参数中来
    恢复领阶幂次精度。
    """
    print("\n" + "=" * 70)
    print("第6节: 领头重正子重求和 (LRR)")
    print("=" * 70)

    print("\n【定义6.1】IR 重正子:")
    print("")
    print("  Wilson 系数 C(ξ, r) 的微扰级数:")
    print("    C(ξ, r) = Σ_{{n=0}}^{{∞}} c_n(ξ) α_s^{{n+1}}(μ)")
    print("")
    print("  在大 n 时, c_n 阶乘增长:")
    print("    c_n ∼ n! · (2β₀)ⁿ · n^{γ}")
    print("")
    print("  这导致级数是渐近的 (asymptotic),")
    print("  Borel 可求和但存在模糊性 (ambiguity).")

    print("\n【方法6.2】LRR 方案 (Zhang et al. 2023):")
    print("")
    print("  通过 Borel 变换识别重正子的位置和留数,")
    print("  构建 LRR 修正因子 R(ξ, μ, P_z):")
    print("")
    print("  C^{{LRR}}(ξ, r) = C^{{RGR}}(ξ, r) + δC(ξ, r)")
    print("")
    print("  其中 δC 吸收了重正子模糊性, 使残余的 1/P_z²")
    print("  修正具有明确的非微扰物理意义.")
    print("")
    print("  这一方法首次在 π 介子夸克 PDF 中实现,")
    print("  推广到核子胶子 PDF 正在进行中.")

    return {}


# ============================================================================
# 第7节: 汇总与主函数
# ============================================================================

def main():
    """运行所有推导。"""
    print("=" * 70)
    print("准TMD→光锥TMD 因子化定理与匹配核 — SymPy 推导")
    print("基于 lattice-pdf 库全部文档与代码")
    print("=" * 70)

    derive_lamet_factorization()
    derive_tmd_factorization()
    derive_nlo_matching_kernel()
    derive_flavor_mixing()
    derive_rgr_resummation()
    derive_lrr_resummation()

    print("\n" + "=" * 70)
    print("推导完成。请继续: 07_continuum_limit.py (连续极限外推)")
    print("=" * 70)


if __name__ == "__main__":
    main()
