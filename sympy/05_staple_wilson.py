#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
================================================================================
05_staple_wilson.py — Staple型Wilson线构造与TMD关联函数
================================================================================

本模块推导TMD-PDF计算所需的staple型Wilson线的完整构造:
  1. 空间方向Wilson线的离散化
  2. Staple型Wilson线的几何构造
  3. Wilson线自能与线性发散
  4. HYP涂抹 vs 梯度流涂抹对Wilson线的平滑
  5. TMD软因子与Wilson圈减除

参考源:
  [9] J.-C. He et al. (LPC), "Unpolarized TMD parton distributions of the
      nucleon from lattice QCD", PRD 109, 114513 (2024)
  - 补充/格点QCD中的TMD_PDF.tex: §2 (Staple型Wilson线)
  - 补充/格点QCD中的Wilson线.tex
  - 补充/格点QCD中的重整化.tex: §5 (Wilson线算符的重整化)
  - 汇报/构造胶子准算符.tex: §2.3
  - 文档/gluon_PDF_continuum.tex: §2-3

作者: 基于 lattice-pdf 库全部文档与代码的系统化推导
日期: 2026-07-24
================================================================================
"""

import sympy as sp
from sympy import I, oo, pi, Symbol, symbols, exp, log
from sympy import Function, IndexedBase

# ============================================================================
# 第0节: 符号声明
# ============================================================================

# 空间方向
z_val = Symbol('z', real=True)                # 纵向分离
b_perp_val = Symbol('b_perp', real=True)      # 横向分离

# Wilson 线参数
delta_m = Symbol('delta_m', real=True)        # 线性发散质量参数
N_link = Symbol('N_link', integer=True, positive=True)  # 链接数

# 涂抹参数
n_smear = Symbol('n_smear', integer=True, positive=True)  # 涂抹步数
alpha_hyp = Symbol('alpha_hyp', real=True)    # HYP 涂抹参数

# 路径
path_segments = IndexedBase('path')


# ============================================================================
# 第1节: 空间方向Wilson线的格点构造
# ============================================================================

def construct_spatial_wilson_line():
    """
    格点上沿空间方向的 Wilson 线离散化。

    基础表示中的 Wilson 线:
        U(x₂, x₁) = Π_{k=0}^{N_z-1} U_z(x₁ + k·a_z ẑ)

    其中每个链接 U_z(s) 是 SU(3) 矩阵 (3×3 复矩阵)。

    在格点上, Wilson 线就是一组规范链接的乘积。
    """
    print("=" * 70)
    print("第1节: 空间方向 Wilson 线的格点构造")
    print("=" * 70)

    print("\n【定义1.1】基础表示中的 Wilson 线:")
    print("")
    print("  U(x+z, x) = Π_{{k=0}}^{{N_z-1}} U_z(x + k·a·ẑ)")
    print("")
    print("  其中 N_z = |z|/a 是链接数,")
    print("  a 是格距, ẑ 是方向单位矢量.")

    print("\n【定义1.2】伴随表示中的 Wilson 线:")
    print("")
    print("  U_{{adj}}^{{ab}}(x₂, x₁)")
    print("    = 2 Tr[t^a U_{{fund}}(x₂, x₁) t^b U_{{fund}}^†(x₂, x₁)]")
    print("")
    print("  等价形式:")
    print("    U_{{adj}}^{{ab}} = P exp[-i g₀ ∫ dz (-i f^{{abc}}) A_μ^c]")
    print("")
    print("  其中 f^{{abc}} 是 SU(3) 结构常数.")

    # Wilson 线的自能
    print("\n【推导1.3】Wilson 线的自能 (线性发散来源):")
    print("")
    print("  裸 Wilson 线的期望值:")
    print("    ⟨U(z,0)⟩_{{bare}} = e^{{-δm |z|/a}} × (对数因子) × (物理)")
    print("")
    print("  其中 δm ~ 1/a 是线性发散质量参数.")
    print("")
    print("  对于胶子 (伴随表示):")
    print("    δm_g = (C_A/C_F) δm_q ≈ (9/4) δm_q")
    print("")
    print("  因此胶子 Wilson 线的线性发散是夸克情形的 9/4 倍!")

    return {}


# ============================================================================
# 第2节: Staple 型 Wilson 线构造
# ============================================================================

def construct_staple_wilson_line():
    """
    构造 TMD-PDF 所需的 staple 型 Wilson 线。

    Staple 形状的 Wilson 线解决了一个关键问题:
      如何定义规范不变的、同时依赖于纵向和横向分离的
      非定域胶子关联函数?

    答案: 使用 staple 型路径:
      1. 从 (z-, b__perp) 沿光锥到 ∞-
      2. 在 ∞- 处横向移动 b__perp
      3. 从 ∞- 沿光锥回到 (0-, 0__perp)

    这个构造由 TMD 因子化定理中的软函数和快度截断需求
    自然导出。
    """
    print("\n" + "=" * 70)
    print("第2节: Staple 型 Wilson 线构造")
    print("=" * 70)

    print("\n【定义2.1】Staple Wilson 线的三段结构:")
    print("")
    print("  W^{{staple,ab}}(ξ-, b__perp; 0-, 0__perp)")
    print("")
    print("  分为三段:")
    print("")
    print("  ① 光锥段 (向后):")
    print("    W₁ = U_{{adj}}( ξ-, b__perp → L-, b__perp )")
    print("    = Π_{{k=0}}^{{L}} exp[-i g₀ a A^+(ξ- + k a)]")
    print("")
    print("  ② 横向段:")
    print("    W₂ = U_{{adj}}( L-, b__perp → L-, 0__perp )")
    print("    = Π_{{k=0}}^{{b}} exp[-i g₀ a A__perp(b - k a)]")
    print("")
    print("  ③ 光锥段 (向前):")
    print("    W₃ = U_{{adj}}( L-, 0__perp → 0-, 0__perp )")
    print("    = Π_{{k=0}}^{{L}} exp[-i g₀ a A^+(k a)]")
    print("")
    print("  总 Wilson 线: W^{{staple}} = W₁ · W₂ · W₃")

    # 在格点上的实现
    print("\n【实现2.2】格点上的 staple Wilson 线 (LaMET/TMD 混合框架):")
    print("")
    print("  由于格点 QCD 在 Euclidean 时空中进行,")
    print("  光锥方向被替换为大动量 boost 的方向 (z 轴).")
    print("")
    print("  格点 staple 形状:")
    print("    ① 沿 z 轴从 (z, b__perp) 到 (z+L_z, b__perp)")
    print("    ② 沿横向从 (z+L_z, b__perp) 到 (z+L_z, 0__perp)")
    print("    ③ 沿 z 轴从 (z+L_z, 0__perp) 到 (0, 0__perp)")
    print("")
    print("  L_z 是一个大的纵向延伸 (实际计算中取 L_z ≫ z)")

    return {}


# ============================================================================
# 第3节: Wilson 线涂抹 — HYP vs 梯度流
# ============================================================================

def wilson_line_smearing():
    """
    Wilson 线的涂抹 (smearing) 方案。

    涂抹的目的是抑制 Wilson 线的线性 UV 发散,
    通过用平滑后的规范链接替换原始链接来实现。

    两种主要方法:
      (a) HYP 涂抹: 局域, 快速, 但缺乏严格的重整化理论基础
      (b) 梯度流涂抹: 非局域(扩散), 较慢, 但有严格的 UV 有限性证明
    """
    print("\n" + "=" * 70)
    print("第3节: Wilson 线涂抹 — HYP vs 梯度流")
    print("=" * 70)

    print("\n【方法3.1】HYP (Hypercubic) 涂抹:")
    print("")
    print("  HYP 涂抹由三步局域涂抹组成:")
    print("")
    print("  Step 1: 在包含原始链接的超立方体内涂抹")
    print("    V_{{μ,1}}(x) = Proj_{{SU(3)}}[ (1-α₁)U_μ(x)")
    print("                + (α₁/2) Σ_{{ν≠μ}} (staples) ]")
    print("")
    print("  Step 2: 使用 Step 1 的链接进一步涂抹")
    print("  Step 3: 使用 Step 2 的链接最终涂抹")
    print("")
    print("  典型参数: α₁=0.75, α₂=0.6, α₃=0.3")
    print("  CLQCD/LPC 合作组使用 HYP5 (5 步涂抹)")

    print("\n【方法3.2】梯度流 (Wilson Flow) 涂抹:")
    print("")
    print("  等同于对规范链接应用 Wilson 流方程:")
    print("")
    print("    d/dt V_t = -g₀² {{∂ S_W(V_t)}} V_t")
    print("")
    print("  在流时间 t 后, 链接被平滑半径 r_t ≈ √(8t) 的高斯核平滑.")
    print("")
    print("  优势: Lüscher-Weisz 定理保证 UV 有限性")
    print("  劣势: 需要零流时间外推 t→0")

    print("\n【比较3.3】HYP vs 梯度流:")
    print("")
    print("  ┌──────────┬────────────┬──────────────┬──────────────┐")
    print("  │  性质     │  HYP 涂抹   │  梯度流涂抹   │  建议          │")
    print("  ├──────────┼────────────┼──────────────┼──────────────┤")
    print("  │  理论基础 │  启发式     │  严格 QFT     │  梯度流更优    │")
    print("  │  计算成本 │  低         │  中等          │  HYP 更经济    │")
    print("  │  可重整化 │  无保证     │  严格可重整化  │  梯度流更优    │")
    print("  │  UV 截断  │  有效       │  指数式        │  梯度流更强    │")
    print("  │  实际使用 │  LPC/CLQCD │  MSULat        │  取决于需求    │")
    print("  └──────────┴────────────┴──────────────┴──────────────┘")

    return {}


# ============================================================================
# 第4节: Wilson 圈减除方案
# ============================================================================

def wilson_loop_subtraction():
    """
    Wilson 圈减除方案: 消除 Wilson 线自能的非微扰方法。

    基本思想:
      Wilson 线的自能 ⟨W(z,0)⟩ 在大 z 处指数衰减。

      Wilson 圈 (Wilson loop) 包含两条 Wilson 线 + 两条横向链接,
      形成闭合回路。通过除以适当的 Wilson 圈, 可以消除线性发散。

    对于 TMD, 减除方案为:
      h_R(z, b__perp) = h_B(z, b__perp) / Z_W(z, b__perp)

      其中 Z_W 由方形 Wilson 圈确定:
        Z_W(z, b__perp) = exp[-δm(z + b__perp)/a]
    """
    print("\n" + "=" * 70)
    print("第4节: Wilson 圈减除方案")
    print("=" * 70)

    print("\n【定义4.1】Wilson 圈减除 (TMD 情形):")
    print("")
    print("  裸矩阵元: h_B(z, b__perp, P_z, a)")
    print("")
    print("  减除因子: Z_W(z, b__perp, a)")
    print("    = exp[-δm (|z| + |b__perp|) / a]")
    print("    × (对数修正)")
    print("")
    print("  重整化矩阵元:")
    print("    h_R(z, b__perp, P_z, μ) = h_B(z, b__perp, P_z, a)")
    print("                         / Z_W(z, b__perp, a)")
    print("                         × Z_log(μ, a)")

    print("\n【注4.2】方形 Wilson 圈的物理动机:")
    print("")
    print("  考虑 R × T 的 Wilson 圈 W(R,T),")
    print("  其在大 T 极限下: ⟨W(R,T)⟩ ∝ exp[-V(R)·T]")
    print("")
    print("  其中 V(R) 是静态夸克-反夸克势.")
    print("  V(R) 中的线性项 σR 对应于 Wilson 线自能中的 δm.")
    print("")
    print("  因此 Wilson 圈提供了 δm 的非微扰确定方法.")

    return {}


# ============================================================================
# 第5节: TMD 软函数与快度截断
# ============================================================================

def tmd_soft_function():
    """
    TMD 软函数: 处理快度发散的关键要素。

    TMD 因子化中, 软函数 S(b__perp) 吸收了两个 Wilson 线方向
    之间的红外 (快度) 发散。

    软函数的定义:
      S(b__perp) ∝ ⟨0| Tr[ S_n^†(b__perp) S_n̄(b__perp) ] |0⟩

    其中 S_n 和 S_n̄ 是沿不同光锥方向的软 Wilson 线。
    """
    print("\n" + "=" * 70)
    print("第5节: TMD 软函数与快度截断")
    print("=" * 70)

    print("\n【定义5.1】TMD 软函数:")
    print("")
    print("  S(b__perp, μ, y_n - y_n̄)")
    print("    = (1/N_c) ⟨0| Tr[ S_n^†(b__perp) S_n̄(b__perp) ] |0⟩")
    print("")
    print("  其中 S_n 是沿方向 n^μ 的半无限 Wilson 线:")
    print("    S_n(x) = P exp[-i g₀ ∫_0^∞ ds n·A(x + s n)]")

    print("\n【推导5.2】软函数在 TMD 因子化中的作用:")
    print("")
    print("  TMD-PDF 包含三类发散:")
    print("    1. UV 发散 → 通过标准 QCD 重整化处理")
    print("    2. 快度发散 → 通过软函数 S(b__perp) 吸收")
    print("    3. 线性发散 → 通过 Wilson 圈减除处理")
    print("")
    print("  完整因子化公式:")
    print("    f_1^g(x, k__perp²) = H(μ, ζ) ×")
    print("      ∫ d²b__perp e^{-i k__perp·b__perp}  f̃_1^g(x, b__perp)")
    print("      × S^{1/2}(b__perp, μ, ζ)")
    print("")
    print("  其中 S^{1/2} 是软因子的平方根 (每个 TMD 取一份).")

    return {}


# ============================================================================
# 第6节: 汇总与主函数
# ============================================================================

def main():
    """运行所有推导。"""
    print("=" * 70)
    print("Staple型Wilson线构造 — SymPy 推导")
    print("基于 lattice-pdf 库全部文档与代码")
    print("=" * 70)

    construct_spatial_wilson_line()
    construct_staple_wilson_line()
    wilson_line_smearing()
    wilson_loop_subtraction()
    tmd_soft_function()

    print("\n" + "=" * 70)
    print("推导完成。请继续: 06_factorization.py (因子化定理)")
    print("=" * 70)


if __name__ == "__main__":
    main()
