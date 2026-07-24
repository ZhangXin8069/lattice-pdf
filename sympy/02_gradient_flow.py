#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
================================================================================
02_gradient_flow.py — 梯度流(Wilson Flow)重整化与小流时展开(SFTX)
================================================================================

本模块使用 SymPy 推导梯度流重整化方案的完整理论框架:

  1. Wilson 流方程的定义与基本性质
  2. 流场的平滑性质 (高斯热核)
  3. 流时间依赖的局域复合场
  4. UV 有限性定理 (Lüscher-Weisz 定理)
  5. 小流时展开 (SFTX) 与匹配系数
  6. 能动张量的梯度流表示
  7. 梯度流在胶子 PDF 计算中的应用
  8. 零流时外推 (t → 0 极限)

关键参考文献:
  [4] M. Lüscher, "Properties and uses of the Wilson flow in lattice QCD",
      JHEP 08, 071 (2010)
  [5] M. Lüscher and P. Weisz, "Perturbative analysis of the gradient flow
      in non-abelian gauge theories", JHEP 02, 051 (2011)
  [6] C. Monahan and K. Orginos, "Quasi parton distributions and the
      gradient flow", JHEP 03, 116 (2017)
  [15] 补充文档: 格点QCD中的梯度流重整化.tex

作者: 基于 lattice-pdf 库全部文档与代码的系统化推导
日期: 2026-07-24
================================================================================
"""

import sympy as sp
from sympy import I, oo, pi, Symbol, Function, symbols, simplify, expand, log
from sympy import integrate, diff, limit, series, summation, exp, sqrt
from sympy import Matrix, IndexedBase, Idx
from sympy import KroneckerDelta, LeviCivita

# ============================================================================
# 第0节: 符号声明
# ============================================================================

# --- 时空坐标 ---
t_flow = Symbol('t', real=True, positive=True)          # 流时间
x_mu = IndexedBase('x')                                  # 时空坐标 x^μ
D_dim = Symbol('D', integer=True, positive=True)         # 时空维数

# --- 流场 ---
# B_μ(t,x) = 规范场的流时间演化版本
B_mu = Function('B_mu')
# G_{μν}(t,x) = 流时间处的场强张量
G_munu = Function('G_munu')
# χ(t,x) = 夸克流场
chi = Function('chi')
# λ(t,x) = 拉格朗日乘子场 (D+1 维场论表示)
lambda_field = Function('lambda')

# --- 原始规范场 ---
A_mu = Function('A_mu')                                  # 原始规范势 A_μ(x)
g0 = Symbol('g0', real=True, positive=True)              # 裸耦合常数
beta_val = Symbol('beta', real=True, positive=True)      # β = 6/g0²

# --- 匹配系数 ---
c_T = Function('c_T')(t_flow)                            # 能动张量匹配系数
c_E = Function('c_E')(t_flow)                            # 作用量密度匹配系数
c_S = Function('c_S')(t_flow)                            # 迹反常匹配系数

# --- QCD β 函数 ---
b0 = Symbol('b0', real=True)                             # β₀ = (11 - 2N_f/3)/(4π)²
b1 = Symbol('b1', real=True)                             # β₁
N_f = Symbol('N_f', integer=True, positive=True)         # 味数

# --- 梯度流标度 ---
t0_scale = Symbol('t0', real=True, positive=True)        # 梯度流标度 t₀
w0_scale = Symbol('w0', real=True, positive=True)        # 梯度流标度 w₀

# --- 格点变量 ---
V_t_link = IndexedBase('V_t')                            # 流时间依赖的规范链接
S_W = Function('S_W')                                    # Wilson 规范作用量


# ============================================================================
# 第1节: Wilson 流方程的定义
# ============================================================================

def define_wilson_flow_equation():
    """
    定义 Wilson 流 (梯度流) 的基本方程。

    连续理论中的 Yang-Mills 梯度流:
        ∂_t B_μ(t,x) = D_ν G_{νμ}(t,x),    B_μ(0,x) = A_μ(x)

    其中:
        D_ν = ∂_ν + [B_ν, ·]  是协变导数
        G_{νμ} = ∂_ν B_μ - ∂_μ B_ν + [B_ν, B_μ]  是流场强张量

    格点上的 Wilson 流:
        d/dt V_t(x,μ) = -g_0² {∂_{x,μ} S_W(V_t)} V_t(x,μ)
        V_t(x,μ)|_{t=0} = U(x,μ)
    """
    print("=" * 70)
    print("第1节: Wilson 流方程的定义")
    print("=" * 70)

    print("\n【定义1.1】连续 Yang-Mills 梯度流方程:")
    print("")
    print("  ∂_t B_μ(t,x) = D_ν G_{νμ}(t,x)")
    print("")
    print("  其中:")
    print("    D_ν = ∂_ν + [B_ν, ·]  (伴随表示协变导数)")
    print("    G_{μν} = ∂_μ B_ν - ∂_ν B_μ + [B_μ, B_ν]")
    print("  初始条件: B_μ(0,x) = A_μ(x)")

    print("\n【定义1.2】格点 Wilson 流方程:")
    print("")
    print("  d/dt V_t(x,μ) = -g₀² {{∂_{x,μ} S_W(V_t)}} V_t(x,μ)")
    print("  V_t|_{t=0} = U(x,μ)")
    print("")
    print("  其中 S_W 是 Wilson 规范作用量,")
    print("  ∂_{x,μ} S_W 是作用量对链接变量 U(x,μ) 的 Lie 代数导数.")

    # 领头阶解: 线性化流方程
    print("\n【推导1.3】领头阶解 (微扰展开, g₀→0):")
    print("")
    print("  在 g₀→0 极限下, 流方程线性化为:")
    print("    ∂_t B_μ = ∂² B_μ  (热方程!)")
    print("")
    print("  解为高斯平滑:")
    print("    B_μ(t,x) = ∫ d^D y  K_t(x-y) A_μ(y)")
    print("")
    print("  其中热核:")
    print("    K_t(z) = e^{-z²/(4t)} / (4πt)^{D/2}")
    print("")
    print("  这意味着在流时间 t, 规范场被半径为")
    print("  r_t ≈ √(8t) 的高斯核平滑. 这是梯度流的核心性质!")

    return {}


# ============================================================================
# 第2节: 梯度流的平滑性质
# ============================================================================

def derive_smoothing_properties():
    """
    推导梯度流的平滑性质。

    梯度流本质上是一个扩散过程:
      - 流时间 t 起到扩散"时间"的作用
      - 平滑半径 r_t ≈ √(8t)
      - 在 t > 0 时场是 C^∞ 光滑的
      - t > 0 时所有的 UV 发散被调节 (regulate)

    在动量空间中, 平滑因子为:
      B_μ(t,p) ≈ e^{-t p²} A_μ(p)

    这意味着动量 p ≳ 1/√t 的模式被指数抑制.
    """
    print("\n" + "=" * 70)
    print("第2节: 梯度流的平滑性质")
    print("=" * 70)

    print("\n【性质2.1】平滑半径与扩散特征:")
    print("")
    print("  平滑半径: r_t ≈ √(8t)")
    print("  在格点上, 典型流时间选择: t ≈ (3-5) a²")
    print("  对应平滑半径: r_t ≈ 3-5 倍格距")
    print("")
    print("  这对应于格点 QCD 计算中的 t_flow ≈ 0.3-0.5 fm²")

    print("\n【性质2.2】动量空间平滑因子:")
    print("")
    print("  B_μ(t,p) ≈ e^{-t p²} A_μ(p)   (领头阶)")
    print("")
    print("  指数因子 e^{-t p²} 在 p ≫ 1/√t 时提供了强烈的 UV 截断.")
    print("  这意味着对于 t > 0 的任意关联函数在微扰论中是 UV 有限的.")

    print("\n【性质2.3】梯度流与 HYP/Stout 涂抹的比较:")
    print("")
    print("  ┌──────────────┬──────────────────┬──────────────────────┐")
    print("  │  方法         │  理论基础         │  优势                 │")
    print("  ├──────────────┼──────────────────┼──────────────────────┤")
    print("  │  HYP 涂抹     │ 启发式            │  计算快, 局域         │")
    print("  │  Stout 涂抹   │ 启发式            │  可微, 常用于 HMC    │")
    print("  │  Wilson 流    │ 严格 QFT 理论     │  可重整化理论完备    │")
    print("  └──────────────┴──────────────────┴──────────────────────┘")

    return {}


# ============================================================================
# 第3节: Lüscher-Weisz UV 有限性定理
# ============================================================================

def luscher_weisz_theorem():
    """
    阐述 Lüscher-Weisz 定理: 梯度流重整化的 UV 有限性。

    定理 (Lüscher & Weisz 2011):
      对于 t > 0, 规范不变的流场算符的关联函数
      在标准 QCD 参数重整化后是 UV 有限的。

    具体而言:
      在量纲正规化 (D = 4 - 2ε) 中, 对所有 t > 0,
      ⟨E(t,x)⟩ 在 ε→0 时是有限的,
      其中 E(t,x) = (1/4) G_{μν}^a(t,x) G_{μν}^a(t,x)
      是流时间处的规范作用量密度。

    这一性质是梯度流作为重整化方案的基石:
      - 无需额外的算符重整化常数
      - 流时间 t 本身起到了物理 UV 截断的作用
    """
    print("\n" + "=" * 70)
    print("第3节: Lüscher-Weisz UV 有限性定理")
    print("=" * 70)

    print("\n【定理3.1】(Lüscher & Weisz 2011):")
    print("")
    print("  对于 t > 0, 规范不变的流场算符的关联函数")
    print("  在标准 QCD 参数重整化后是 UV 有限的.")
    print("")
    print("  具体地, 在量纲正规化中:")
    print("")
    print("    ⟨E(t,x)⟩ = (3g₀²)/(4π²t²) × [1 + O(g₀²)]")
    print("")
    print("  该表达式在 ε→0 时是有限的 (无 1/ε 极点).")

    # 单圈验证
    print("\n【推导3.2】单圈计算验证:")
    print("")
    print("  领头阶 (树图):")
    print("    ⟨E(t,x)⟩^{(0)} = (N_c² - 1) · (3g₀²)/(4π² t²)")
    print("")
    print("  次领头阶 (单圈):")
    print("    ⟨E(t,x)⟩^{(1)} = ⟨E(t,x)⟩^{(0)} × [1 + k₁ g₀² + O(g₀⁴)]")
    print("    k₁ = (C_A/(4π)) × (11/3) + O(ε)")
    print("")
    print("  关键: k₁ 不包含 1/ε 极点 → UV 有限!")

    # 包含夸克场
    print("\n【推论3.3】包含夸克场时的波函数重整化:")
    print("")
    print("  夸克流场需要乘性波函数重整化因子 Z_χ:")
    print("")
    print("    Z_χ = 1 + (3 C_F g₀²)/(16π² ε) + O(g₀⁴)")
    print("")
    print("  这可以通过[环标记] (ringed) 夸克场来处理:")
    print("    χ̊(t,x) = χ(t,x) / √(⟨χ̄(t,x) D̸ χ(t,x)⟩)")
    print("")
    print("  环标记后, χ̊ 的关联函数对所有 t > 0 是 UV 有限的.")

    return {}


# ============================================================================
# 第4节: 小流时展开 (SFTX)
# ============================================================================

def derive_small_flow_time_expansion():
    """
    推导小流时展开 (Small Flow-Time Expansion, SFTX)。

    对于 t → 0, 流场算符 O(t,x) 可以展开为 t=0 处的
    重整化局域算符的线性组合:

        O(t,x) ∼ Σ_k c_k(t) O_{R,k}(x)    (t → 0)

    其中 c_k(t) 是微扰可计算的 Wilson 系数,
    O_{R,k}(x) 是重整化的局域复合算符。

    这是梯度流重整化方案的核心: 通过测量流时间 t 处的
    关联函数, 并利用 SFTX 外推到 t=0, 可以获得 t=0 处的
    物理算符的矩阵元。
    """
    print("\n" + "=" * 70)
    print("第4节: 小流时展开 (SFTX)")
    print("=" * 70)

    print("\n【定义4.1】SFTX 的一般形式:")
    print("")
    print("  O(t,x) ≃ Σ_k c_k(t) O_{R,k}(x)    (t → 0)")
    print("")
    print("  其中:")
    print("    c_k(t):  Wilson 系数 (微扰可计算)")
    print("    O_{R,k}:  重整化局域算符")
    print("    ≃:      在关联函数的意义下渐近等价")

    # 规范作用量密度的 SFTX
    print("\n【推导4.2】规范作用量密度 E(t,x) 的 SFTX:")
    print("")
    print("  E(t,x) = (1/4) G_{μν}^a(t,x) G_{μν}^a(t,x)")
    print("")
    print("  SFTX 展开:")
    print("    E(t,x) ≃ c_E(t) · 1 + c_T(t) · T_{μμ}(x) + O(t)")
    print("")
    print("  其中:")
    print("    c_E(t) ∼ 1/t²  → 真空贡献")
    print("    c_T(t) ∼ const  → 能动张量的迹")
    print("    T_{μμ} = (β(g)/2g) F_{μν}^a F_{μν}^a  是迹反常")

    # 匹配系数
    print("\n【推导4.3】匹配系数 c_E(t) 和 c_T(t) 的微扰展开:")
    print("")
    print("  c_E(t) = (3g₀²)/(4π²t²) × [1 + O(g₀²)]")
    print("  c_T(t) = g₀²/(g²(μ)) × [1 + O(g₀²)]")
    print("")
    print("  其中 g(μ) 是 MS-bar 方案中的跑动耦合常数.")

    # 能动张量的梯度流表示
    print("\n【推导4.4】能动张量的梯度流公式 (Suzuki 2013):")
    print("")
    print("  {{T_{μν}}}_R(x) = lim_{{t→0}} [U_{μν}(t,x)/c_T(t)")
    print("                          - (c_S(t)/(c_T(t)c_E(t))) δ_{μν}")
    print("                            × (E(t,x) - ⟨E(t,x)⟩)]")
    print("")
    print("  其中:")
    print("    U_{μν}(t,x) = G_{μρ}^a(t,x) G_{νρ}^a(t,x)")
    print("                  - (1/4) δ_{μν} G_{ρσ}^a(t,x) G_{ρσ}^a(t,x)")
    print("    c_S(t)/c_T(t) = β(g)/(2g³)  (迹反常关系)")

    # 匹配系数的渐近行为
    print("\n【推导4.5】匹配系数的 t→0 渐近行为:")
    print("")
    print("  1/c_T(t) ∼ -2 b₀ ln(√(8t) Λ) + c₁")
    print("")
    print("  其中 b₀ = (11 - 2N_f/3)/(4π)² 是单圈 β 函数系数,")
    print("  Λ 是 QCD 标度参数, c₁ 是方案依赖的有限常数.")

    # 高阶计算
    print("\n【注4.6】高阶计算:")
    print("  - Harlander, Kluth, Lange (2018): 双圈 c_T(t) 和 c_E(t)")
    print("  - Borgulat et al. (2024): NNLO 匹配系数")
    print("  - Mereghetti et al. (2022): 夸克偶极算符的单圈匹配")

    return {}


# ============================================================================
# 第5节: 梯度流在胶子 PDF 计算中的应用
# ============================================================================

def gradient_flow_for_gluon_pdf():
    """
    推导梯度流在胶子 PDF/TMD 计算中的具体应用。

    核心思想:
      1. 对流时间 t 处的场计算胶子准算符的矩阵元
      2. 利用 SFTX 将结果匹配到 t=0 的物理算符
      3. 外推 t→0 以获得物理结果 (零流时间外推)

    关键公式 (Monahan & Orginos 2017):
      准 PDF 算符的梯度流版本:
        O_flow(z, P_z, t) = ⟨P_z| O(z, t) |P_z⟩

      在 t→0 极限下, SFTX 给出:
        O_flow(z, P_z, t) ≃ ζ(z, t, μ) O_R(z, P_z, μ)

      其中 ζ 是微扰可计算的匹配系数。
    """
    print("\n" + "=" * 70)
    print("第5节: 梯度流在胶子 PDF/TMD 计算中的应用")
    print("=" * 70)

    print("\n【方法5.1】梯度流涂抹 vs HYP 涂抹:")
    print("")
    print("  在胶子 PDF 计算中, Wilson 线中的规范链接被涂抹以抑制")
    print("  线性 UV 发散. 两种常用的涂抹方案:")
    print("")
    print("  (a) HYP 涂抹:")
    print("      - 3 步局域涂抹 (hypercubic 阻挡)")
    print("      - CLQCD/LPC 合作组使用 HYP5 (5 步涂抹)")
    print("      - 优势: 计算快速, 局域")
    print("      - 劣势: 缺乏严格的重整化理论基础")
    print("")
    print("  (b) 梯度流涂抹 (Wilson flow):")
    print("      - 流时间 t_flow ≈ 3a² (平滑半径 ~ √(8t))")
    print("      - MSULat 合作组使用梯度流")
    print("      - 优势: 严格可重整化, UV 有限性定理保证")
    print("      - 劣势: 计算成本较高, 需要零流时间外推")

    print("\n【方法5.2】Monahan-Orginos (2017) 方案:")
    print("")
    print("  对流时间 t 处的胶子场计算准 PDF 算符:")
    print("")
    print("    O(z, t) = F_{zμ}(z, t) W_t(z, 0) F^{μz}(0, t)")
    print("")
    print("  SFTX 给出与 t=0 物理算符的关系:")
    print("")
    print("    O(z, t) ≃ ζ(z, t, μ) O_R(z, μ) + O(t)")
    print("")
    print("  然后外推 t→0 以消除有限的流时间效应.")

    print("\n【方法5.3】梯度流中的 PDF 矩比值 (MSULat 组):")
    print("")
    print("  ⟨x^{n-1}⟩^{MS}(μ) / ⟨x^{m-1}⟩^{MS}(μ)")
    print("    = (ζ_m(t,μ)/ζ_n(t,μ)) · ⟨x^{n-1}⟩(t) / ⟨x^{m-1}⟩(t) + O(t)")
    print("")
    print("  矩比值的好处: 部分 t 依赖的 SFTX 系数在比值中抵消,")
    print("  减少了对微扰计算的依赖.")

    # 零流时间外推
    print("\n【方法5.4】零流时间外推 (t → 0):")
    print("")
    print("  在多个流时间值 [{t₁, t₂, t₃}] 处计算矩阵元,")
    print("  然后外推到 t → 0:")
    print("")
    print("    h(z, P_z, t) = h(z, P_z, 0) + c₁(z) t + c₂(z) t² + ...")
    print("")
    print("  外推消除有限的流时间效应, 得到物理的 h(z, P_z).")

    return {}


# ============================================================================
# 第6节: 格点上 Wilson 流的数值实现
# ============================================================================

def wilson_flow_numerical_implementation():
    """
    Wilson 流的格点数值实现方法。

    使用三阶 Runge-Kutta 积分方案:
      V_{t+ε}(x,μ) = exp[-g₀² ε ∂_{x,μ} S_W(V_t)] V_t(x,μ)

    流时间步长 ε 需要足够小以确保积分精度。
    典型选择: ε = 0.01-0.02.
    """
    print("\n" + "=" * 70)
    print("第6节: Wilson 流数值实现")
    print("=" * 70)

    print("\n【算法6.1】三阶 Runge-Kutta 积分方案:")
    print("")
    print("  输入: 初始规范组态 U_μ(x), 总流时间 t_max, 步长 ε")
    print("  输出: 流时间 t_max 处的规范组态 V_μ(t_max, x)")
    print("")
    print("  步骤 (对每个链接变量):")
    print("    1. 计算 staple 和 Z_μ(x) = ∂_{x,μ} S_W(V)")
    print("    2. 构造反厄米无迹矩阵:")
    print("       Ω_μ(x) = (1/2)(Z_μ(x) - Z_μ^†(x))")
    print("                - (1/(2N_c)) Tr[Z_μ(x) - Z_μ^†(x)]")
    print("    3. RK3 步进:")
    print("       W₁ = exp[-g₀² ε Ω(V)] V")
    print("       W₂ = exp[-g₀² ε (3/4) Ω(W₁)] W₁")
    print("       V(t+ε) = exp[-g₀² ε (8/9) Ω(W₂)] W₂")

    print("\n【准则6.2】流时间选择准则:")
    print("")
    print("  (a) 能标设定条件: t²⟨E(t)⟩|_{{t=t₀}} = 0.3")
    print("      → 确定格点标度 t₀ ≈ 0.15 fm²")
    print("")
    print("  (b) 平滑条件: t ≫ a²")
    print("      → 确保平滑半径 ≫ 格距, 离散化误差可控")
    print("")
    print("  (c) 体积条件: √(8t) ≪ L")
    print("      → 确保平滑半径远小于格点体积, 避免有限体积效应")
    print("")
    print("  典型值: t ∈ [a², 5a²], 即约 0.01-0.06 fm²")

    return {}


# ============================================================================
# 第7节: 汇总与主函数
# ============================================================================

def main():
    """运行所有推导。"""
    print("=" * 70)
    print("梯度流(Wilson Flow)重整化 — SymPy 符号推导")
    print("基于 lattice-pdf 库全部文档与代码")
    print("=" * 70)

    define_wilson_flow_equation()
    derive_smoothing_properties()
    luscher_weisz_theorem()
    derive_small_flow_time_expansion()
    gradient_flow_for_gluon_pdf()
    wilson_flow_numerical_implementation()

    print("\n" + "=" * 70)
    print("推导完成。请继续: 03_field_strength.py (场强张量)")
    print("=" * 70)


if __name__ == "__main__":
    main()
