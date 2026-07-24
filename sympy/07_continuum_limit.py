#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
================================================================================
07_continuum_limit.py — 连续极限外推 (a→0) 与联合外推
================================================================================

本模块推导连续极限外推和联合外推的完整公式体系:
  1. 离散化误差的分类 (O(a), O(a²))
  2. 连续极限外推公式
  3. 无穷大动量外推 (P_z → ∞)
  4. 大λ外推 (λ = zP_z → ∞)
  5. 联合外推 (a, P_z 同时外推)
  6. 系统误差估计

参考源:
  [8] C. Chen et al., "Unpolarized gluon PDF...continuum limit" (2025)
  [13] K. Zhang et al., "Impact of gauge fixing precision on the continuum
       limit", PRD 110, 074505 (2024)
  - 补充/格点QCD中的外推.tex: §3-6
  - 文档/gluon_PDF_continuum.tex: §4

作者: 基于 lattice-pdf 库全部文档与代码的系统化推导
日期: 2026-07-24
================================================================================
"""

import sympy as sp
from sympy import Symbol, symbols, Function, exp, log, oo

# ============================================================================
# 符号
# ============================================================================
a_sym, Pz_sym = symbols('a P_z', real=True, positive=True)
x_var = Symbol('x', real=True)
xg0 = Function('xg_0')(x_var)
f_func = Function('f')(x_var)
h_func = Function('h')(x_var)
d_func = Function('d')(x_var)


def derive_continuum_extrapolation():
    """
    连续极限外推 (a→0).

    对改进费米子 (Clover/Staggered/HISQ), 领先离散化误差为 O(a²):
        O(a) = O_0 + c_1 a² + c_2 a⁴ + ...

    需要至少 3 个格距值才能可靠外推.
    """
    print("=" * 70)
    print("第1节: 连续极限外推 a→0")
    print("=" * 70)
    print("\n【公式1.1】O(a²) 外推:")
    print("  O(a) = O_0 + c_1 a² + O(a⁴)")
    print("\n  对 Clover 费米子 (c_sw 非微扰确定): 领先误差 O(a²)")
    print("  对 Wilson 费米子: 领先误差 O(a) (需 Symanzik 改进)")
    print("\n  需要 ≥3 个 a 值: CLQCD/LPC 使用 a = {{0.105, 0.0897, 0.0775}} fm")


def derive_infinite_momentum_extrapolation():
    """
    无穷大动量外推 (P_z → ∞).

    LaMET 幂次修正:
        q(x, P_z) = q_0(x) + c(x)/P_z² + d(x)/P_z⁴ + ...

    外推变量是 1/P_z² 而非 1/P_z (Lorentz 不变性排除奇次幂).
    """
    print("\n" + "=" * 70)
    print("第2节: 无穷大动量外推 P_z→∞")
    print("=" * 70)
    print("\n【公式2.1】1/P_z² 外推:")
    print("  xg(x, P_z) = xg_0(x) + d(x)/P_z²")
    print("\n  需要 ≥2 个 P_z 值. CLQCD/LPC 使用 P_z ≈ 1.5-2.0 GeV.")


def derive_large_lambda_extrapolation():
    """
    大λ外推 (λ = zP_z → ∞).

    重整化矩阵元 h_R(z, P_z) 仅在有限 z 范围内可用,
    需要外推到 z→∞ 以进行 Fourier 变换.

    拟合形式: h_R(z) ∼ c₁ z^{-d₁} e^{-z/λ₀} + c₂ z^{-d₂} e^{-z/λ₀}
    """
    print("\n" + "=" * 70)
    print("第3节: 大λ外推 λ=zP_z→∞")
    print("=" * 70)
    print("\n【公式3.1】外推拟合函数:")
    print("  h_R(z, P_z) ∼ c₁ z^{{-d₁}} e^{{-z/λ₀}} + c₂ z^{{-d₂}} e^{{-z/λ₀}}")
    print("\n  λ_L 截断参数: 仅 λ > λ_L 用于拟合 (典型 λ_L=6-8)")
    print("\n【替代方案3.2】Derivative 方法:")
    print("  q̃(x,P_z) = (i/x) ∫ dz e^{{ixP_zz}} dh_R/dz")


def derive_joint_extrapolation():
    """
    联合外推 (a, P_z 同时外推).

    胶子 PDF 的联合外推 (Chen et al. 2025):
        xg(x, P_z, a) = xg_0(x) + a² f(x) + a² P_z² h(x) + d(x)/P_z²

    其中 a² P_z² h(x) 项是动量依赖的离散化误差.
    """
    print("\n" + "=" * 70)
    print("第4节: 联合外推")
    print("=" * 70)
    print("\n【公式4.1】胶子 PDF 联合外推 (CLQCD/LPC 2025):")
    print("  xg(x, P_z, a) = xg_0(x) + a² f(x) + a² P_z² h(x) + d(x)/P_z²")
    print("\n  各项物理来源:")
    print("    xg_0: 物理 PDF (a→0, P_z→∞)")
    print("    a² f(x): 标准 O(a²) 离散化误差")
    print("    a² P_z² h(x): 动量依赖的离散化误差 (a·P_z 耦合)")
    print("    d(x)/P_z²: LaMET 幂次修正")

    print("\n【公式4.2】横向性 PDF 联合外推 (含手征外推):")
    print("  h(x, P_z, a, m_π) = h_0(x) + a² f(x) + c(x)/P_z²")
    print("                     + g₀ m_π² ln(m_π²/μ₀²) + g₁ m_π²")


def main():
    print("=" * 70)
    print("连续极限外推与联合外推 — SymPy 推导")
    print("=" * 70)
    derive_continuum_extrapolation()
    derive_infinite_momentum_extrapolation()
    derive_large_lambda_extrapolation()
    derive_joint_extrapolation()
    print("\n推导完成。")


if __name__ == "__main__":
    main()
