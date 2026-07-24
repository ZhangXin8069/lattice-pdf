#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
================================================================================
09_correlation_fn.py — 两点/三点关联函数与矩阵元提取
================================================================================

推导格点上关联函数的构造与矩阵元提取:
  1. 质子两点关联函数 (Distillation + 动量涂抹)
  2. 非连通三点关联函数 (Disconnected 3pt)
  3. 矩阵元提取 (比值法, 双态拟合)
  4. VVV 重子块构造
  5. Wick 收缩 (直接项 + 交换项)

参考源:
  - code/gluon_pdf_full_workflow.py: §5, §7
  - examples/2pt_proton_Cg5gmu_*.py: 质子2pt实现
  - 补充/格点QCD中的关联函数.tex
  - 补充/格点QCD中的维克收缩与傅里叶变换.tex
  - 补充/格点QCD中的连通图与非连通图.tex
  - Note of gluon PDFs (内部笔记): §1

作者: 基于 lattice-pdf 库全部文档与代码的系统化推导
日期: 2026-07-24
================================================================================
"""

import sympy as sp
from sympy import I, Symbol, symbols, Function, exp, pi, oo


def derive_two_point_function():
    """
    质子两点关联函数 (Distillation 框架).

    C₂(Δt, P) = ⟨N(P, t_sink) N̄(P, t_src)⟩

    其中 N 是质子插值算符:
        N = P⁺ ε^{abc} (u^a^T C γ₅ γₜ d^b) u^c
    """
    print("=" * 70)
    print("第1节: 质子两点关联函数")
    print("=" * 70)

    print("\n【定义1.1】质子插值算符:")
    print("  N_α(P,t) = Σ_x e^{-iP·x} ε^{{abc}}")
    print("             [u^{{aT}}(x,t) C γ₅ γ_t d^b(x,t)] u_α^c(x,t)")
    print("\n【定义1.2】两点函数:")
    print("  C₂(Δt, P) = ⟨ Tr[ P₊ N(P, t_sink) N̄(P, t_src) ] ⟩")
    print("\n  其中 P₊ = (1+γ₄)/2 是正宇称投影算符.")

    print("\n【蒸馏因子化1.3】:")
    print("  C₂ = Φ_{abc}(t_sink)")
    print("     × τ^{ad}(t_sink, t_src)")
    print("     × (Cγ₅γ_t) τ^{be}(t_sink, t_src) (Cγ₅γ_t)^T")
    print("     × τ^{cf}(t_sink, t_src)")
    print("     × Φ_{def}*(t_src)")
    print("\n  其中 Φ 是 VVV 重子块, τ 是 perambulator.")


def derive_disconnected_three_point():
    """
    非连通三点关联函数.

    胶子流插入与核子无夸克线连接 → Disconnected 图.

    C₃(z, P_z, t_sep) = ⟨O_g(z) C₂(P_z)⟩ - ⟨O_g(z)⟩⟨C₂(P_z)⟩

    需要真空期望值减除.
    """
    print("\n" + "=" * 70)
    print("第2节: 非连通三点关联函数")
    print("=" * 70)

    print("\n【定义2.1】非连通三点函数:")
    print("  C₃(z, P_z, t_sep)")
    print("    = ⟨(O_g(z,t_i) - ⟨O_g⟩)(C₂^N(P_z,t_snk,t_src) - ⟨C₂^N⟩)⟩")
    print("\n  关键: 胶子流和质子2pt可独立计算 → 灵活性")
    print("  但信噪比远差于连通图 → 需要大量组态")


def derive_matrix_element_extraction():
    """
    矩阵元提取方法.

    比值法:
        h(z, P_z) = lim_{Δt→∞} C₃(z, P_z, Δt) / C₂(P_z, Δt)

    在 disconnect 近似下:
        h(z) = ⟨O(z)⟩ (对时间和组态平均)

    双态拟合:
        R(Δt) = c₀ + c₁ e^{-ΔE·Δt}
    """
    print("\n" + "=" * 70)
    print("第3节: 矩阵元提取")
    print("=" * 70)

    print("\n【方法3.1】比值法 (Plateau 拟合):")
    print("  h(z, P_z) = average_{{Δt∈[t_min, t_max]}} C₃(Δt,z)/C₂(Δt)")
    print("\n【方法3.2】双态拟合:")
    print("  R(Δt, z) = c₀(z) + c₁(z) e^{-ΔE(z) Δt}")
    print("\n【方法3.3】Jackknife 误差:")
    print("  σ_JK² = (N-1)/N · Σ_i (θ_i - θ̄)²")


def derive_vvv_and_wick():
    """
    VVV 重子块与 Wick 收缩.

    VVV: Φ_{abc}(P) = Σ_x e^{-iP·x} ε_{ijk} v_i^a(x) v_j^b(x) v_k^c(x)

    Wick 收缩: 直接项 - 交换项
    """
    print("\n" + "=" * 70)
    print("第4节: VVV 重子块与 Wick 收缩")
    print("=" * 70)

    print("\n【定义4.1】VVV 重子块:")
    print("  Φ_{{abc}}(P) = Σ_x e^{{-iP·x}} ε_{{ijk}} v_i^a(x) v_j^b(x) v_k^c(x)")
    print("\n  其中 v_i^a(x) 是 Laplace 本征矢量 (蒸馏).")
    print("\n【定义4.2】Wick 收缩 (质子):")
    print("  C₂ = Tr[Φ(t_snk)·τ·(Cγ₅γ_t)·τ·(Cγ₅γ_t)^T·τ·Φ(t_src)*]  [直接项]")
    print("      - Tr[Φ(t_snk)·τ·(Cγ₅γ_t)·τ·(Cγ₅γ_t)^T·τ·Φ(t_src)*]  [交换项]")


def main():
    print("=" * 70)
    print("关联函数与矩阵元提取 — SymPy 推导")
    print("=" * 70)
    derive_two_point_function()
    derive_disconnected_three_point()
    derive_matrix_element_extraction()
    derive_vvv_and_wick()
    print("\n推导完成。")


if __name__ == "__main__":
    main()
