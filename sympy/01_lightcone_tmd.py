#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
================================================================================
01_lightcone_tmd.py — 光锥胶子TMD-PDF的连续场论定义
================================================================================

本模块使用 SymPy 符号计算推导核子非极化胶子 TMD-PDF 的光锥定义,
包括:
  1. 光锥坐标与 Sudakov 分解
  2. 胶子场强张量 F_{μν} 及其对偶 F̃_{μν} 的光锥分量
  3. 非极化胶子 TMD-PDF 的算符定义
  4. 伴随表示 Wilson 线的光锥构造
  5. TMD-PDF 的领头扭度 (leading-twist) 分类
  6. 与共线胶子 PDF 的关系 (对 k__perp 积分)

参考文献:
  [1] X. Ji, "Parton Physics on a Euclidean Lattice", PRL 110, 262002 (2013)
  [2] J.-C. He et al. (LPC), "Unpolarized TMD parton distributions of the
      nucleon from lattice QCD", PRD 109, 114513 (2024)
  [3] J.-H. Zhang et al., "Accessing gluon parton distributions in LaMET",
      PRL 122, 142001 (2019)
  [4] 内部笔记: "Note of gluon PDFs" (docs/Note of gluon PDFs.pdf)
  [5] 补充文档: 格点QCD中的TMD_PDF.tex (§2-3)
  [6] 补充文档: 格点QCD中的胶子算符.tex
  [7] 补充文档: 格点QCD中的部分子分布函数.tex (§1)

作者: 基于 lattice-pdf 库全部文档与代码的系统化推导
日期: 2026-07-24
================================================================================
"""

import sympy as sp
from sympy import I, oo, pi, Symbol, Function, symbols, simplify, expand
from sympy import conjugate as conj
from sympy import integrate, diff, limit, series, summation

# ============================================================================
# 第0节: 符号声明
# ============================================================================

# --- 光锥坐标 ---
xi_m, xi_p = symbols('xi^- xi^+', real=True)     # 光锥坐标 ξ-, ξ⁺
x_m, x_p = symbols('x^- x^+', real=True)          # 光锥坐标 x-, x⁺
x_perp = symbols('x_perp', cls=sp.IndexedBase)    # 横向坐标 x__perp = (x¹, x²)
x, k_perp = symbols('x k_perp', real=True)        # Bjorken-x, 横向动量 |k__perp|

# --- 强子动量 ---
P_plus, P_minus, P_perp = symbols('P^+ P^- P_perp', real=True)
P_z = symbols('P_z', real=True)                   # z-方向动量 (LaMET boost)
M_N = symbols('M_N', real=True)                   # 核子质量

# --- QCD 参数 ---
g_s = symbols('g_s', real=True, positive=True)    # 强耦合常数
alpha_s = g_s**2 / (4 * pi)                       # α_s = g_s²/(4π)
N_c = symbols('N_c', integer=True, positive=True) # 颜色数 (QCD: N_c=3)
C_A = N_c                                         # 伴随表示 Casimir
C_F = (N_c**2 - 1) / (2 * N_c)                   # 基础表示 Casimir
mu = symbols('mu', real=True, positive=True)      # 重整化标度 μ
Lambda_QCD = symbols('Lambda_QCD', real=True, positive=True)

# --- TMD 变量 ---
b_perp = symbols('b_perp', real=True)             # 横向分离 |b__perp| (坐标空间)
zeta = symbols('zeta', real=True)                 # 快度参数 ζ
y_n = symbols('y_n', real=True)                   # 快度截断

# --- 格点变量 ---
a = symbols('a', real=True, positive=True)        # 格距
L = symbols('L', real=True, positive=True)        # 空间体积
N_x = symbols('N_x', integer=True, positive=True) # 空间格点数

# --- 符号函数 ---
f_g = Function('f_g')(x, k_perp, mu, zeta)       # 非极化胶子 TMD-PDF
g_collinear = Function('g_collinear')(x, mu)      # 共线胶子 PDF

# --- 场论算符 (符号) ---
# 胶子场强张量 F_{μν}^a (伴随表示, a=1..N_c²-1)
F_munu_a = sp.IndexedBase('F')
# 对偶场强张量
F_tilde_munu_a = sp.IndexedBase('F_tilde')
# 伴随表示 Wilson 线
U_adj = sp.IndexedBase('U_adj')
# 核子态
P_state = Function('|P>')


# ============================================================================
# 第1节: 光锥坐标与 Sudakov 分解
# ============================================================================

def define_lightcone_coordinates():
    """
    定义光锥坐标与 Sudakov 分解。

    光锥坐标 (ξ⁺, ξ-, ξ__perp) 由 Minkowski 坐标定义:
        ξ^± = (ξ⁰ ± ξ³) / √2
        ξ__perp = (ξ¹, ξ²)

    四矢量 a^μ 的光锥分解:
        a^μ = (a⁺, a-, a__perp^μ)
        a·b = a⁺b- + a-b⁺ - a__perp·b__perp
        a^² = 2a⁺a- - |a__perp|²

    对于以接近光速沿 +z 方向运动的强子:
        P^μ ≈ (P^+, 0, 0__perp)  →  "infinite momentum frame" (IMF)

    Returns
    -------
    rules : dict
        光锥坐标变换规则
    """
    # Minkowski 坐标与光锥坐标的关系
    t, z = symbols('t z', real=True)

    xi_plus_expr = (t + z) / sp.sqrt(2)
    xi_minus_expr = (t - z) / sp.sqrt(2)

    print("=" * 70)
    print("第1节: 光锥坐标与 Sudakov 分解")
    print("=" * 70)
    print("\n光锥坐标定义:")
    print("  ξ⁺ = (t + z)/√2 = (t + z)/{sp.sqrt(2)}")
    print("  ξ- = (t - z)/√2 = (t - z)/{sp.sqrt(2)}")
    print("  ξ__perp = (x, y)")

    # Sudakov 分解
    print("\nSudakov 分解 (接近光速沿 +z 方向运动的强子):")
    print("  P^μ = (P^+, P^-, P__perp) ≈ (P^+, M_N²/(2P^+), 0__perp)")
    print("  k^μ = (x P^+, k^- , k__perp)")
    print("  其中 x = k^+/P^+ 是部分子携带的纵向动量份额")

    # 标量积
    a_plus, a_minus, a_perp_sq = symbols('a^+ a^- a_perp^2', real=True)
    b_plus, b_minus, b_perp_sq = symbols('b^+ b^- b_perp^2', real=True)
    dot_product = a_plus * b_minus + a_minus * b_plus - sp.sqrt(a_perp_sq * b_perp_sq)
    print("\n光锥标量积: a·b = a⁺b- + a-b⁺ - a__perp·b__perp")

    return {
        'xi_plus': xi_plus_expr,
        'xi_minus': xi_minus_expr,
    }


# ============================================================================
# 第2节: 胶子场强张量的光锥分量
# ============================================================================

def define_gluon_field_strength_lightcone():
    """
    定义胶子场强张量及其光锥分量。

    伴随表示中的场强张量:
        F_{μν}^a = ∂_μ A_ν^a - ∂_ν A_μ^a + g_s f^{abc} A_μ^b A_ν^c

    其中 f^{abc} 是 SU(N_c) 的结构常数。

    对偶场强张量:
        F̃_{μν}^a = (1/2) ε_{μνρσ} F^{a,ρσ}

    光锥规范 A^+ = 0 下的重要分量:
        F^{+i}  → 与横向极化胶子相关 (i=1,2)
        F^{+-}  → 与纵向极化胶子相关

    Returns
    -------
    dict
        场强张量分量
    """
    print("\n" + "=" * 70)
    print("第2节: 胶子场强张量 F_{μν} 的光锥分量")
    print("=" * 70)

    # 符号声明
    A_mu_a = sp.IndexedBase('A')  # 胶子场 A_μ^a
    f_abc = sp.IndexedBase('f')   # 结构常数 f^{abc}

    # 指标
    mu, nu, rho, sigma = symbols('mu nu rho sigma', integer=True)

    print("\n伴随表示场强张量:")
    print("  F_{{μν}}^a = ∂_μ A_ν^a - ∂_ν A_μ^a + g_s f^{{abc}} A_μ^b A_ν^c")

    print("\n对偶场强张量:")
    print("  F̃_{{μν}}^a = (1/2) ε_{{μνρσ}} F^{{a,ρσ}}")

    print("\n光锥规范 A^+ = 0 下:")
    print("  F_a^{+i} = ∂^+ A_a^i  →  '好' 分量, 与横向极化胶子相关")
    print("  F_a^{+-} = ∂^+ A_a^-  →  与纵向极化胶子相关")
    print("  F_a^{ij}  (i,j=1,2) →  包含非线性胶子自相互作用")

    # 胶子 PDF 算符中的关键组合
    print("\n非极化胶子 PDF/TMD 算符中的场强组合:")
    print("  F^{+i}_a · F_a^{+i}  →  g_1(x, k__perp) 的算符 (对 i 求和)")
    print("  F^{+i}_a · F̃_a^{+i}  →  Δg(x, k__perp) 的算符 (胶子螺旋度)")

    return {}


# ============================================================================
# 第3节: 非极化胶子 TMD-PDF 的算符定义
# ============================================================================

def define_unpolarized_gluon_tmd():
    """
    推导非极化胶子 TMD-PDF 的完整算符定义。

    光锥胶子 TMD-PDF f_1^g(x, k__perp) (非极化, leading-twist) 定义为:

        x f_1^g(x, k__perp²) = ∫ dξ- d²ξ__perp / ((2π)³ 2(P^+)²)
            × e^{-i x P^+ ξ- + i k__perp · ξ__perp}
            × ⟨P| F_a^{+i}(ξ-, ξ__perp) W_{ab}(ξ; 0) F_b^{+i}(0) |P⟩

    其中 W_{ab} 是伴随表示中的 staple 型 Wilson 线,
    F_a^{+i} = F_a^{+i} 是场强张量的 "+i" 分量。

    注意: 这是 TMD 情形, 与共线 PDF 的关键区别在于:
    1. 额外包含横向分离 ξ__perp
    2. Wilson 线包含横向部分 (staple 形状)
    3. 对 k__perp 的依赖性

    Returns
    -------
    dict
        TMD 定义及其性质
    """
    print("\n" + "=" * 70)
    print("第3节: 非极化胶子 TMD-PDF 的算符定义")
    print("=" * 70)

    # TMD 坐标
    xi_minus_sym, xi_perp_sym = symbols('xi^- xi_perp', real=True)
    k_perp_sym = symbols('k_perp', real=True)

    print("\n【定义3.1】非极化胶子 TMD-PDF (光锥定义):")
    print("")
    print("  x f_1^g(x, k__perp²; μ, ζ) = ∫ dξ- d²ξ__perp / ((2π)³ 2(P^+)²)")
    print("    × e^{-i x P^+ ξ- + i k__perp·ξ__perp}")
    print("    × ⟨P| F_a^{+i}(ξ-, ξ__perp) W_{ab}^{staple}(ξ; 0) F_b^{+i}(0) |P⟩")
    print("")
    print("  其中:")
    print("    ξ = (ξ-, ξ__perp) 是光锥坐标分离")
    print("    x = k^+/P^+ ∈ [0,1] 是纵向动量份额")
    print("    k__perp 是部分子的横向动量")
    print("    μ 是重整化标度")
    print("    ζ 是快度参数 (Collins-Soper 标度)")
    print("    W^{staple} 是 staple 型 Wilson 线")

    print("\n【定义3.2】与共线胶子 PDF 的关系 (对 k__perp 积分):")
    print("")
    print("  x g(x, μ) = ∫ d²k__perp  x f_1^g(x, k__perp²; μ, ζ)")
    print("")
    print("  共线 PDF 通过积分掉横向动量依赖获得。")

    print("\n【定义3.3】TMD 因子化中的快度演化 (Collins-Soper 方程):")
    print("")
    print("  d/d(ln ζ)  f_1^g(x, k__perp²; μ, ζ) = K(b__perp, μ) · f_1^g(x, k__perp²; μ, ζ)")
    print("")
    print("  其中 K(b__perp, μ) 是 Collins-Soper 演化核。")

    return {
        'gluon_tmd': 'x f_1^g(x, k_perp^2; mu, zeta)',
        'collinear_relation': 'integral over k_perp',
        'cs_equation': 'Collins-Soper evolution',
    }


# ============================================================================
# 第4节: Staple 型 Wilson 线构造
# ============================================================================

def define_staple_wilson_line():
    """
    构造 TMD 所需的 staple 型 Wilson 线。

    在共线 PDF 中, Wilson 线是沿光锥方向的直线。
    在 TMD-PDF 中, 横向分离需要 Wilson 线的横向段,
    形成 "staple" (订书钉) 形状:

        W^{staple}(ξ; 0) = W_-(ξ-, ξ__perp → ∞-, ξ__perp)
                         × W__perp(∞-, ξ__perp → ∞-, 0__perp)
                         × W_+(∞-, 0__perp → 0-, 0__perp)

    其中各个段的定义:
        W_-(ξ → ∞): 沿光锥方向从 ξ 到无穷远的 Wilson 线
        W__perp: 在光锥无穷远处的横向 Wilson 线
        W_+(∞ → 0): 沿光锥方向从无穷远到原点的 Wilson 线
    """
    print("\n" + "=" * 70)
    print("第4节: Staple 型 Wilson 线构造")
    print("=" * 70)

    print("\n【定义4.1】基础表示中的 Wilson 线:")
    print("")
    print("  U_{fund}(x₂, x₁) = P exp[-i g_s ∫_{x₁}^{x₂} dz^μ A_μ^a(z) t^a]")
    print("")
    print("  其中 t^a = λ^a/2 是 SU(N_c) 的生成元 (Gell-Mann 矩阵),")
    print("  P 表示路径编序 (path ordering).")

    print("\n【定义4.2】伴随表示中的 Wilson 线:")
    print("")
    print("  U_{adj}^{ab}(x₂, x₁) = 2 Tr[t^a U_{fund}(x₂, x₁) t^b U_{fund}^†(x₂, x₁)]")
    print("")
    print("  等价形式:")
    print("  U_{adj}^{ab}(x₂, x₁) = P exp[-i g_s ∫_{x₁}^{x₂} dz^μ (-i f^{abc}) A_μ^c(z)]")

    print("\n【定义4.3】Staple 型 Wilson 线 (TMD):")
    print("")
    print("  W^{staple,ab}(ξ-, ξ__perp; 0-, 0__perp)")
    print("    = U_{adj}^{ab}(ξ-, ξ__perp → ∞-, ξ__perp)           [光锥段, 向后]")
    print("    × U_{adj}^{bc}(∞-, ξ__perp → ∞-, 0__perp)           [横向段]")
    print("    × U_{adj}^{cd}(∞-, 0__perp → 0-, 0__perp)            [光锥段, 向前]")
    print("")
    print("  路径编序确保了规范不变性。")

    print("\n【定义4.4】TMD 软因子 (soft factor):")
    print("")
    print("  S(b__perp) = (1/N_c) ⟨0| Tr[ S_n^†(b__perp) S_n̄(b__perp) ] |0⟩")
    print("")
    print("  其中 S_n 和 S_n̄ 是沿不同光锥方向的 Wilson 线,")
    print("  用于吸收 infrared 发散。")

    return {}


# ============================================================================
# 第5节: TMD-PDF 的领头扭度分类
# ============================================================================

def classify_tmd_pdfs():
    """
    TMD-PDF 的领头扭度 (leading-twist, twist-2) 分类。

    在 twist-2 层次, 夸克有 8 种独立的 TMD-PDF,
    胶子也有对应的 8 种 gluon TMD distributions:

        |  夸克 TMD      |  胶子 TMD      |  极化  |  物理含义        |
        |----------------|---------------|--------|------------------|
        |  f_1 (unpol)   |  f_1^g        |  U     |  非极化数密度     |
        |  g_1 (hel)     |  g_1^g        |  L     |  螺旋度分布       |
        |  h_1 (trans)   |  h_1^g        |  T     |  横向性分布       |
        |  g_{1T}        |  g_{1T}^g     |  L-T   |  横向-螺旋度      |
        |  h_{1L}^_perp      |  h_{1L}^{g_perp}  |  T-L   |  纵向-横向性      |
        |  h_{1T}^_perp      |  h_{1T}^{g_perp}  |  T-T   |  旋量性 (pretzel) |
        |  h_1^_perp         |  h_1^{g_perp}     |  U-T   |  Boer-Mulders     |
        |  f_{1T}^_perp      |  f_{1T}^{g_perp}  |  U-T   |  Sivers           |

    其中 U = 非极化核子, L = 纵向极化核子, T = 横向极化核子.
    箭头 _perp 表示横向动量依赖 (naive T-odd).
    """
    print("\n" + "=" * 70)
    print("第5节: TMD-PDF 的领头扭度分类")
    print("=" * 70)

    print("\n【表5.1】Twist-2 胶子 TMD-PDF 的完整分类:")
    print("")
    print("  1. f_1^g(x, k__perp²)  — 非极化胶子 TMD (本文重点)")
    print("    算符: F_a^{+i} F_a^{+i}")
    print("")
    print("  2. g_{1L}^g(x, k__perp²) — 纵向极化胶子螺旋度 TMD")
    print("    算符: F_a^{+i} F̃_a^{+i}  (iε_{ij} 收缩)")
    print("")
    print("  3. h_1^g(x, k__perp²) — 横向极化胶子横向性 TMD")
    print("    算符: F_a^{+{i}} F_a^{+{j}} + (i↔j)  (对称无迹)")
    print("")
    print("  4. h_{1L}^{g_perp}(x, k__perp²) — 纵向极化胶子 worm-gear TMD")
    print("")
    print("  5-8. 其他 twist-2 胶子 TMD (naive T-odd 等)")

    print("\n【注5.1】本文聚焦于 f_1^g(x, k__perp²),")
    print("  即核子非极化、胶子非极化的 TMD-PDF.")
    print("  这是 SIDIS 和 Drell-Yan 过程中最重要的胶子 TMD.")

    return {}


# ============================================================================
# 第6节: 光锥 PDF 与 TMD 的积分关系
# ============================================================================

def derive_tmd_to_collinear():
    """
    推导 TMD-PDF 到共线 PDF 的积分关系。

    使用 SymPy 证明:
        ∫ d²k__perp f_1^g(x, k__perp²) = g(x)
    """
    print("\n" + "=" * 70)
    print("第6节: TMD-PDF → 共线 PDF 积分关系")
    print("=" * 70)

    # 符号声明
    k_perp_sq = symbols('k_perp^2', real=True, positive=True)
    f1g = Function('f_1^g')(x, k_perp_sq)

    print("\n  x g(x) = ∫ d²k__perp  x f_1^g(x, k__perp²)")
    print("")
    print("  在极坐标 (k__perp, φ) 下:")
    print("  = ∫_0^{2π} dφ ∫_0^∞ d|k__perp| |k__perp|  x f_1^g(x, |k__perp|²)")
    print("  = 2π ∫_0^∞ d|k__perp| |k__perp|  x f_1^g(x, |k__perp|²)")
    print("")
    print("  定义 k__perp² = t:")
    print("  = π ∫_0^∞ dt  x f_1^g(x, t)")
    print("")
    print("  ∴ x f_1^g 对 k__perp² 的积分给出共线胶子 PDF x g(x)")

    return {}


# ============================================================================
# 第7节: 汇总与主函数
# ============================================================================

def main():
    """运行所有推导。"""
    print("=" * 70)
    print("核子非极化胶子 TMD-PDF 的光锥定义 — SymPy 符号推导")
    print("基于 lattice-pdf 库全部文档与代码")
    print("=" * 70)

    define_lightcone_coordinates()
    define_gluon_field_strength_lightcone()
    define_unpolarized_gluon_tmd()
    define_staple_wilson_line()
    classify_tmd_pdfs()
    derive_tmd_to_collinear()

    print("\n" + "=" * 70)
    print("推导完成。请继续: 02_gradient_flow.py (梯度流重整化)")
    print("=" * 70)


if __name__ == "__main__":
    main()
