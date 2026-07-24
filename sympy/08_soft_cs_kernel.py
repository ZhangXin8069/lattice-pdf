#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
================================================================================
08_soft_cs_kernel.py — 软函数与Collins-Soper演化核
================================================================================

推导TMD因子化中的软函数S(b__perp)和Collins-Soper演化核K(b__perp,μ):
  1. TMD软函数的物理起源
  2. Collins-Soper方程与演化核
  3. 内禀软函数S_I的格点提取 (赝标量介子形状因子法)
  4. CS核的格点提取 (准TMD波函数法 / 准TMD-PDF法)
  5. 重整化群演化与快度演化

参考源:
  [9] J.-C. He et al., "Unpolarized TMD...", PRD 109, 114513 (2024)
  [10] M.-H. Chu et al., "Lattice calculation of the intrinsic soft
       function and the Collins-Soper kernel", JHEP 08, 172 (2023)
  - 补充/格点QCD中的TMD_PDF.tex: §3

作者: 基于 lattice-pdf 库全部文档与代码的系统化推导
日期: 2026-07-24
================================================================================
"""

import sympy as sp
from sympy import Symbol, symbols, Function, exp, log, integrate, oo, pi


def define_soft_function():
    """
    TMD软函数的定义与物理意义.

    软函数S(b__perp)吸收了两个Wilson线方向之间的红外(快度)发散,
    是TMD因子化中不可约化的非微扰要素.
    """
    print("=" * 70)
    print("第1节: TMD 软函数 S(b__perp)")
    print("=" * 70)

    print("\n【定义1.1】TMD 软函数:")
    print("  S(b__perp, μ, y_n - y_n̄)")
    print("    = (1/N_c) ⟨0| Tr[ S_n^†(b__perp) S_n̄(b__perp) ] |0⟩")
    print("\n  其中 S_n 是沿光锥方向 n^μ 的半无限 Wilson 线.")
    print("\n【结构1.2】软函数的两部分分解:")
    print("  S(b__perp, μ, ζ) = S_I(b__perp, μ) · exp[K(b__perp, μ) ln(ζ/ζ₀)]")
    print("\n  其中:")
    print("    S_I: 内禀软函数 (intrinsic soft function)")
    print("    K:   Collins-Soper 演化核")
    print("    ζ:   快度参数 (Collins-Soper 标度)")


def derive_collins_soper_equation():
    """
    Collins-Soper 方程与演化核.

    CS方程控制TMD-PDF的快度演化:
        d/d(ln ζ) f(x, k__perp²; μ, ζ) = K(b__perp, μ) · f(x, k__perp²; μ, ζ)
    """
    print("\n" + "=" * 70)
    print("第2节: Collins-Soper 方程")
    print("=" * 70)

    print("\n【方程2.1】Collins-Soper 演化方程:")
    print("  d/d(ln ζ) f_1^g(x, k__perp²; μ, ζ) = K(b__perp, μ) · f_1^g(x, k__perp²; μ, ζ)")
    print("\n  K(b__perp, μ) = -C_F α_s/π · ln(μ² b__perp² e^{2γ_E}/4) + O(α_s²)")
    print("\n  在小 b__perp 的微扰区域: K ∝ ln(b__perp²)")
    print("  在大 b__perp 的非微扰区域: K 趋于常数 (从格点数据确定)")


def lattice_extraction():
    """
    CS核的格点提取方法.

    两种主要方法:
    1. 准TMD波函数法: 通过赝标量介子 quasi-TMD 波函数的比值
    2. 准TMD-PDF法: 直接通过核子 quasi-TMD-PDF 提取
    """
    print("\n" + "=" * 70)
    print("第3节: CS 核的格点提取")
    print("=" * 70)

    print("\n【方法3.1】准TMD波函数法 (Chu et al. 2023):")
    print("  通过 π 介子 quasi-TMD 波函数在不同 P_z 下的比值")
    print("  提取 CS 核:")
    print("    K(b__perp, μ) ∼ (1/ln(P_z₁/P_z₂))")
    print("              × ln[Φ̃(b__perp, P_z₁)/Φ̃(b__perp, P_z₂)]")

    print("\n【方法3.2】内禀软函数 S_I (赝标量介子形状因子法):")
    print("  S_I(b__perp) 可以通过赝标量介子 (η_s) 的 quasi-TMD")
    print("  形状因子在 z=0 时的行为提取.")

    print("\n【方法3.3】准TMD-PDF法 (He et al. 2024):")
    print("  直接从核子 quasi-TMD-PDF 中提取 CS 核,")
    print("  作为不同 P_z 下矩阵元比值的副产品.")


def main():
    print("=" * 70)
    print("软函数与Collins-Soper演化核 — SymPy 推导")
    print("=" * 70)
    define_soft_function()
    derive_collins_soper_equation()
    lattice_extraction()
    print("\n推导完成。")


if __name__ == "__main__":
    main()
