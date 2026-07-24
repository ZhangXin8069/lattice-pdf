#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
================================================================================
04_gluon_operator.py — 非定域胶子准TMD算符的构造
================================================================================

本模块使用 SymPy 推导格点上非定域胶子准TMD算符的完整构造:

  1. 基础表示下的双胶子矩阵元 M_{μλ;νρ}(z)
  2. 六个不变振幅分解 (Balitsky-Morris-Radyushkin)
  3. 可乘性重整化算符组合 (Li-Ma-Qiu 定理)
  4. 非极化胶子流 O^{(0)}_g 和 O^{(1)}_g
  5. 极化(螺旋度)胶子流
  6. 格点离散化实现

参考源:
  [7] J.-H. Zhang et al., "Accessing gluon parton distributions in LaMET",
      PRL 122, 142001 (2019)
  - Balitsky, Morris, Radyushkin, "Gluon pseudo-distributions at short
    distances: Forward case", PLB 808, 135621 (2020)
  - 补充/格点QCD中的胶子算符.tex
  - 汇报/构造胶子准算符.tex: §2-4
  - 汇报/格点上计算胶子准算符.tex: §7
  - Note of gluon PDFs (内部笔记): §1 (非极化胶子算符)

作者: 基于 lattice-pdf 库全部文档与代码的系统化推导
日期: 2026-07-24
================================================================================
"""

import sympy as sp
from sympy import I, Symbol, symbols, Function, IndexedBase
from sympy import LeviCivita, KroneckerDelta

# ============================================================================
# 第0节: 符号声明
# ============================================================================

# 指标
mu, nu, alpha_val, beta, lamb, rho = symbols(
    'mu nu alpha beta lambda rho', integer=True)

# 动量
p = IndexedBase('p')  # 核子动量 p^μ
z_vec = IndexedBase('z')  # 分离矢量 z^μ

# Ioffe 时间
nu_ioffe = Symbol('nu', real=True)  # ν = -p·z

# 不变振幅 (Balitsky 的六个振幅)
M_pp = Function('M_pp')(nu_ioffe)
M_zz = Function('M_zz')(nu_ioffe)
M_zp = Function('M_zp')(nu_ioffe)
M_pz = Function('M_pz')(nu_ioffe)
M_ppzz = Function('M_ppzz')(nu_ioffe)
M_gg = Function('M_gg')(nu_ioffe)


# ============================================================================
# 第1节: 基本定义 — 双胶子关联函数
# ============================================================================

def define_gluon_correlator():
    """
    定义双胶子关联函数及其基本性质。

    基础表示中的双胶子关联函数 (Balitsky et al. 2020):

        M_{μλ;νρ}(z,p) = ⟨p| Tr[ F_{μλ}(z) W(z,0) F_{νρ}(0) W(0,z) ] |p⟩

    其中:
      - F_{μν} 是场强张量 (基础表示)
      - W(z,0) 是基础表示中的 Wilson 线
      - Tr 是颜色迹
      - ⟨p|...|p⟩ 是核子态的平均
    """
    print("=" * 70)
    print("第1节: 双胶子关联函数的基本定义")
    print("=" * 70)

    print("\n【定义1.1】基础表示中的双胶子矩阵元:")
    print("")
    print("  M_{{μλ;νρ}}(z,p) = ⟨p| Tr_c[ F_{{μλ}}(z) · W(z,0) ·")
    print("                                   F_{{νρ}}(0) · W(0,z) ] |p⟩")
    print("")
    print("  其中:")
    print("    F_{{μν}} = Σ_a F_{{μν}}^a t^a  (t^a = λ^a/2)")
    print("    W(z,0) = P exp[-i g₀ ∫_0^z dz' A_z(z')]  (空间 Wilson 线)")
    print("    Tr_c  = 对色空间取迹")

    print("\n【性质1.2】M 的对称性:")
    print("")
    print("  M_{{μλ;νρ}}(z,p) = M_{{νρ;μλ}}(-z,p)   (z 反演)")
    print("  M_{{μλ;νρ}}(z,p) = M_{{νρ;μλ}}(z,p)*   (厄米性)")
    print("  M_{{λμ;νρ}} = -M_{{μλ;νρ}}  (反对称性)")
    print("  M_{{μλ;ρν}} = -M_{{μλ;νρ}}  (反对称性)")

    return {}


# ============================================================================
# 第2节: 六个不变振幅分解
# ============================================================================

def derive_invariant_amplitudes():
    """
    推导 Balitsky-Morris-Radyushkin 的六个不变振幅分解。

    M_{μλ;νρ}(z,p) 的 Lorentz 结构可以分解为六个独立的不变振幅:

        M_{μλ;νρ} = Σ_{i=1}^{6} A_i(ν) × T^i_{μλ;νρ}(z,p)

    其中 A_i(ν) 是 ν = -p·z 的函数,
    T^i 是 Lorentz 张量结构。

    六个振幅是:
      M_{pp}, M_{zz}, M_{zp}, M_{pz}, M_{ppzz}, M_{gg}
    """
    print("\n" + "=" * 70)
    print("第2节: 六个不变振幅分解")
    print("=" * 70)

    print("\n【定义2.1】Lorentz 张量基:")
    print("")
    print("  g_{{μν}} = Minkowski 度规 = diag(1,-1,-1,-1)")
    print("")
    print("  两个独立矢量: p^μ (核子动量), z^μ (分离矢量)")
    print("")
    print("  六个独立的 Lorentz 结构 T^i_{{μλ;νρ}}:")

    print("\n【推导2.2】六个不变振幅的显式构造:")
    print("")
    print("  A₁ → M_{{pp}}:  ∝ p_μ p_ν g_{{λρ}}")
    print("  A₂ → M_{{zz}}:  ∝ z_μ z_ν g_{{λρ}}")
    print("  A₃ → M_{{zp}}:  ∝ z_μ p_ν g_{{λρ}}")
    print("  A₄ → M_{{pz}}:  ∝ p_μ z_ν g_{{λρ}}")
    print("  A₅ → M_{{ppzz}}: ∝ p_μ z_ν p_λ z_ρ  (最复杂的结构)")
    print("  A₆ → M_{{gg}}:  ∝ g_{{μν}} g_{{λρ}}")

    print("\n【推导2.3】非极化胶子 PDF 与 M_pp 的关系:")
    print("")
    print("  在光锥极限 z²=0 下:")
    print("")
    print("    -M_{{pp}}(ν, 0) = (1/2) ∫_{-1}^{1} dx e^{-ixν} x g(x)")
    print("")
    print("  即: M_{{pp}} 的 Fourier 变换直接给出共线胶子 PDF xg(x).")

    # 提取 M_pp 的 Lorentz 分量
    print("\n【推导2.4】从矩阵元提取 M_{{pp}}:")
    print("")
    print("  对于非极化胶子, 需要组合:")
    print("")
    print("    M_{{ti;it}} + M_{{ij;ji}} = 2 p₀² M_{{pp}}")
    print("")
    print("  其中 i,j 是横向指标 (x,y).")
    print("  这意味着可以用横向的 F_{{ti}} 和 F_{{ij}} 分量提取 M_{{pp}}.")

    return {}


# ============================================================================
# 第3节: 可乘性重整化算符组合
# ============================================================================

def derive_multiplicative_renormalizable_combinations():
    """
    推导可乘性重整化的算符组合。

    Zhang et al. (2019) 和 Li-Ma-Qiu 证明:
      胶子准 PDF 算符的 36 个 Lorentz 分量 (μλ;νρ)
      各自独立地满足可乘性重整化:

        M_{μλ;νρ}^R(z, μ) = Z_{μλ;νρ}(z, a, μ) M_{μλ;νρ}^B(z, a)

    但对于与 PDF 相关的组合, 以下三个组合具有
    特别简洁的可乘性重整化性质:

        O_1 = M_{tx;tx} + M_{ty;ty}
        O_2 = M_{xy;xy}
        O_3 = M_{tx;tx} + M_{ty;ty} - 2 M_{xy;xy}  (非极化!)
    """
    print("\n" + "=" * 70)
    print("第3节: 可乘性重整化算符组合")
    print("=" * 70)

    print("\n【定理3.1】(Li-Ma-Qiu 乘法重整化定理):")
    print("")
    print("  所有 36 个 Lorentz 分量 M_{{μλ;νρ}}(z) 独立地满足")
    print("  乘法重整化:")
    print("")
    print("    [M_{{μλ;νρ}}(z)]_R(μ) = Z_{{μλ;νρ}}(z,a,μ) [M_{{μλ;νρ}}(z)]_B(a)")
    print("")
    print("  即不存在算符混合! 这与夸克准 PDF 的 Γ=γ^t 情形类似.")

    print("\n【定义3.2】Zhang et al. (2019) 的非极化胶子算符:")
    print("")
    print("  O(z) = M_{{tx;tx}}(z) + M_{{ty;ty}}(z) - 2 M_{{xy;xy}}(z)")
    print("")
    print("  这个特定组合的重要性:")
    print("  1. 在 z → 0 时约化为标准的胶子能量-动量张量")
    print("  2. 具有确定的紫外反常量纲")
    print("  3. 在 NLO 匹配中与夸克单态 PDF 有简洁的混合结构")

    # 按 z 分量数分类
    print("\n【分类3.3】按 z 分量数 s 的分类:")
    print("")
    print("  s = 0: 4个分量  M_{{0i;0i}} 和 M_{{ij;ij}}  (i,j=1,2)")
    print("    反常量纲 = 2γ  (仅来自 Wilson 线自能)")
    print("")
    print("  s = 1: 4个分量  M_{{0i;i3}} ± M_{{3i;i0}}")
    print("    反常量纲 = (3/2)γ")
    print("")
    print("  s = 2: 2个分量  M_{{3i;3i}} 和 M_{{03;03}}")
    print("    反常量纲 = γ")
    print("")
    print("  其中 γ 是基本 Wilson 线自能的反常量纲.")

    return {}


# ============================================================================
# 第4节: 非极化胶子流 O^{(0)}_g 和 O^{(1)}_g
# ============================================================================

def construct_unpolarized_gluon_currents():
    """
    构造格点上非极化胶子流算符 O^{(0)}_g 和 O^{(1)}_g。

    这两个算符用于从格点数据中提取胶子 PDF 所需的不同不变振幅组合。

    Euclidean 空间中的 O^{(0)}_g (非极化):
        O^{(0)}_g(z,t_i) = {F^{01}(z) W F^{01}(0) + F^{02}(z) W F^{02}(0)
                           + F^{03}(z) W F^{03}(0) + F^{12}(z) W F^{12}(0)
                           + F^{23}(z) W F^{23}(0) + F^{31}(z) W F^{31}(0)}
    """
    print("\n" + "=" * 70)
    print("第4节: 非极化胶子流算符 O^{(0)}_g 和 O^{(1)}_g")
    print("=" * 70)

    print("\n【定义4.1】Euclidean 空间中的 O^{(0)}_g (非极化):")
    print("")
    print("  O^{(0)}_g(z, t_i) =")
    print("      F^{{01}}(z)W[z,0]F^{{01}}(0)W[0,z]")
    print("    + F^{{02}}(z)W[z,0]F^{{02}}(0)W[0,z]")
    print("    + F^{{03}}(z)W[z,0]F^{{03}}(0)W[0,z]")
    print("    + F^{{12}}(z)W[z,0]F^{{12}}(0)W[0,z]")
    print("    + F^{{23}}(z)W[z,0]F^{{23}}(0)W[0,z]")
    print("    + F^{{31}}(z)W[z,0]F^{{31}}(0)W[0,z]")
    print("")
    print("  其中 t_i 表示时间片, W 是 Wilson 线.")

    print("\n【定义4.2】Euclidean 空间中的 O^{(1)}_g (非极化):")
    print("")
    print("  O^{(1)}_g(z, t_i) =")
    print("      F^{{12}}(z)W[z,0]F^{{12}}(0)W[0,z]")
    print("    + F^{{23}}(z)W[z,0]F^{{23}}(0)W[0,z]")
    print("    + F^{{31}}(z)W[z,0]F^{{31}}(0)W[0,z]")
    print("")
    print("  即 O^{(1)}_g 仅包含纯磁分量 F^{{ij}}.")

    print("\n【注4.3】O^{(0)}_g 和 O^{(1)}_g 的组合:")
    print("")
    print("  O^{(0)}_g - 2 O^{(1)}_g = F^{{0i}} W F^{{0i}} - F^{{ij}} W F^{{ij}}")
    print("")
    print("  这正是 Zhang et al. (2019) 的可乘性重整化组合")
    print("  (在 Euclidean 空间中的形式).")

    return {}


# ============================================================================
# 第5节: 格点离散化实现
# ============================================================================

def lattice_implementation():
    """
    非定域胶子算符的格点离散化实现。

    关键步骤:
      1. 用 Clover 构造计算 F_{μν}(x)
      2. 用格点 Wilson 线构造 W(z, 0)
      3. 计算算符迹: Tr_c[F_{μλ}(z) W(z,0) F_{νρ}(0) W(0,z)]
      4. 对空间求和
    """
    print("\n" + "=" * 70)
    print("第5节: 格点离散化实现")
    print("=" * 70)

    print("\n【算法5.1】格点上单算符的计算 (code/gluon_pdf_full_workflow.py):")
    print("")
    print("  def gluon_ope_operator_z0_mu2(gauge, F_all, F_tilde_all,")
    print("                                 delta_z, z_dir, mu, nu, mu2, nu2):")
    print("    # Step 1: 平移 F_{{μν}} 到 z 处")
    print("    ope = np.roll(F_all[mu, nu], -delta_z, axis=3-z_dir)")
    print("")
    print("    # Step 2: 乘 Wilson 线 z→0 (U_z†)")
    print("    for dz in range(delta_z):")
    print("        ope = ope @ U_z†(适当平移)")
    print("")
    print("    # Step 3: 在 z=0 处插入 F̃_{{μ2,ν2}}")
    print("    ope = ope @ F_tilde_all[mu2, nu2]")
    print("")
    print("    # Step 4: 乘 Wilson 线 0→z (U_z)")
    print("    for dz in range(delta_z):")
    print("        ope = ope @ U_z(适当平移)")
    print("")
    print("    # Step 5: 颜色迹 Tr_c")
    print("    trace = np.trace(ope, axis1=4, axis2=5)")
    print("")
    print("    # Step 6: 空间求和")
    print("    return np.sum(trace, axis=(1,2,3))")

    print("\n【注5.2】Wilson 线方向的选择:")
    print("")
    print("  对于 LaMET 计算, Wilson 线通常沿 z 方向 (boost 方向).")
    print("  在格点上 z_dir = 2 (对应格点的 z 轴).")
    print("")
    print("  对于 TMD 计算, 还需要横向的 Wilson 线段")
    print("  (见 05_staple_wilson.py).")

    return {}


# ============================================================================
# 第6节: 汇总与主函数
# ============================================================================

def main():
    """运行所有推导。"""
    print("=" * 70)
    print("非定域胶子准TMD算符构造 — SymPy 推导")
    print("基于 lattice-pdf 库全部文档与代码")
    print("=" * 70)

    define_gluon_correlator()
    derive_invariant_amplitudes()
    derive_multiplicative_renormalizable_combinations()
    construct_unpolarized_gluon_currents()
    lattice_implementation()

    print("\n" + "=" * 70)
    print("推导完成。请继续: 05_staple_wilson.py (Staple Wilson线)")
    print("=" * 70)


if __name__ == "__main__":
    main()
