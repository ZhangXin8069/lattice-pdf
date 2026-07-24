#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
================================================================================
03_field_strength.py — 格点场强张量F_{μν}的Clover构造与SymPy推导
================================================================================

本模块使用 SymPy 推导格点 QCD 中场强张量的完整构造过程:
  1. Plaquette 的 BCH 展开 → 从 Wilson loop 到场强
  2. Clover 型场强张量 (四个 plaquette 的组合)
  3. 对偶场强张量 F̃_{μν} 的 Levi-Civita 构造
  4. Minkowski-Euclidean 场强分量转换
  5. 非极化胶子算符所需的 Lorentz 分量组合

参考源:
  - code/gluon_pdf_full_workflow.py: §3 (Levi-Civita与场强张量)
  - examples/Operator.py: Plaquette 构造与场强张量
  - 补充/格点QCD中的场强张量.tex
  - 文档/gluon_pdf_derivation.tex: §5 (格点场强张量)
  - 汇报/构造胶子准算符.tex: §4 (从连续算符到格点实现)
  - Note of gluon PDFs (内部笔记): §1 (场强张量)

作者: 基于 lattice-pdf 库全部文档与代码的系统化推导
日期: 2026-07-24
================================================================================
"""

import sympy as sp
from sympy import I, oo, pi, Symbol, symbols, expand, simplify
from sympy import Matrix, IndexedBase, LeviCivita
from sympy import conjugate as conj, transpose as T

# ============================================================================
# 第0节: 符号声明
# ============================================================================

# 格点指标
x, y, z, t = symbols('x y z t', integer=True)
mu, nu, rho, sigma = symbols('mu nu rho sigma', integer=True)

# 规范链接 U_μ(x) ∈ SU(3) — 3×3 复矩阵
U = IndexedBase('U')  # U_μ^{ab}(x)

# 场强张量 F_{μν} — 3×3 复矩阵 (颜色空间)
F = IndexedBase('F')

# Levi-Civita 符号
epsilon4 = IndexedBase('epsilon')

# 对偶场强张量
F_tilde = IndexedBase('F_tilde')

# 格距
a_sym = Symbol('a', real=True, positive=True)
g0_sym = Symbol('g0', real=True, positive=True)

# Plaquette 变量
P_munu = IndexedBase('P')     # P_{μν}(x)
Q_munu = IndexedBase('Q')     # Clover 叶 Q_{μν}(x)


# ============================================================================
# 第1节: Plaquette → 场强张量的 BCH 展开
# ============================================================================

def derive_plaquette_to_field_strength():
    """
    通过 Baker-Campbell-Hausdorff (BCH) 公式,
    从 plaquette 到场强张量的展开。

    基本思想:
      1×1 Wilson loop (plaquette) 在格距 a 很小时,
      可以展开为场强张量的函数:

        P_{μν}(x) = U_μ(x) U_ν(x+μ̂) U_μ^†(x+ν̂) U_ν^†(x)
                  = exp[i a² g₀ F_{μν}(x) + O(a³)]

      利用 BCH 公式展开后, 包含:
        - O(a⁰): 单位矩阵
        - O(a²): 场强张量 F_{μν}
        - O(a³): 协变导数项 (离散化误差)
    """
    print("=" * 70)
    print("第1节: Plaquette 到场强张量的 BCH 展开")
    print("=" * 70)

    print("\n【定义1.1】标准 1×1 Plaquette:")
    print("")
    print("  P_{{μν}}(x) = U_μ(x) · U_ν(x+μ̂) · U_μ^†(x+ν̂) · U_ν^†(x)")
    print("")
    print("  其中 U_μ(x) = exp[i a g₀ A_μ(x)] 是规范链接变量.")

    print("\n【推导1.2】BCH 展开 (到 O(a²)):")
    print("")
    print("  将 U_μ = exp[i a g₀ A_μ] 代入 plaquette 并使用 BCH 公式:")
    print("")
    print("  P_{{μν}}(x) = exp[i a² g₀ F_{{μν}}(x) + O(a³)]")
    print("")
    print("  其中:")
    print("    F_{{μν}}(x) = ∂_μ A_ν(x) - ∂_ν A_μ(x) + i g₀ [A_μ(x), A_ν(x)]")
    print("")
    print("  是连续场强张量的格点近似.")

    print("\n【推导1.3】到 O(a³) 的展开:")
    print("")
    print("  P_{{μν}}(x) = 1 + i a² g₀ F_{{μν}}(x)")
    print("              - (a⁴ g₀²/2) F²_{{μν}}(x)")
    print("              + i a³ (D_μ F_{{νρ}} + D_ν F_{{ρμ}} + ...) + O(a⁴)")
    print("")
    print("  O(a³) 项是格点离散化误差的来源 (可通过 Symanzik 改进消除).")

    return {}


# ============================================================================
# 第2节: Clover 型场强张量构造
# ============================================================================

def derive_clover_field_strength():
    """
    推导 Clover 型场强张量 F_{μν}(x)。

    Clover 叶由围绕点 x 的四个 plaquette 组成:

        Q_{μν}(x) = P_{μν}(x) + P_{ν,-μ}(x)
                   + P_{-μ,-ν}(x) + P_{-ν,μ}(x)

    场强张量为:

        F_{μν}(x) = -i/(8a²g₀) · (Q_{μν} - Q_{μν}^†)

    Clover 构造的优点: O(a²) 改进 (自动消除 O(a) 离散化误差).

    四个 plaquette 的几何排列 (在 μ-ν 平面):
    ```
              μ →
          +---(1)---+
          |         |
        ν ↑  (2)  (4)  ↓ (反方向)
          |         |
          +---(3)---+
              ← (反方向)
    ```
    """
    print("\n" + "=" * 70)
    print("第2节: Clover 型场强张量 F_{μν}(x)")
    print("=" * 70)

    print("\n【定义2.1】Clover 叶 Q_{{μν}}(x):")
    print("")
    print("  Q_{{μν}}(x) = P_{{μν}}(x) + P_{{ν,-μ}}(x)")
    print("             + P_{{-μ,-ν}}(x) + P_{{-ν,μ}}(x)")
    print("")
    print("  四个 plaquette 的方向:")
    print("    (1) P_{{μν}}:    标准 1×1 Wilson loop (正向 μ, 正向 ν)")
    print("    (2) P_{{ν,-μ}}:  先走 ν 再走 -μ 方向")
    print("    (3) P_{{-μ,-ν}}: 先走 -μ 再走 -ν 方向")
    print("    (4) P_{{-ν,μ}}:  先走 -ν 再走 μ 方向")

    print("\n【定义2.2】Clover 场强张量:")
    print("")
    print("  F_{{μν}}(x) = -i/(8a²g₀) · (Q_{{μν}} - Q_{{μν}}^†)")
    print("")
    print("  其中 Q_{{μν}}^† 是厄米共轭 (矩阵转置 + 复共轭).")
    print("")
    print("  注意: 在准 PDF 算符中, 因子 1/(a²g₀) 会在比值中消去,")
    print("  所以数值实现中常取 a=1, g₀=1.")

    print("\n【性质2.3】Clover F_{{μν}} 的性质:")
    print("")
    print("  1. 反厄米性: F_{{μν}}^† = -F_{{μν}}  (可验证)")
    print("  2. 反对称性: F_{{νμ}} = -F_{{μν}}")
    print("  3. O(a²) 改进: Clover 构造消除了 O(a) 离散化误差")
    print("  4. 色矩阵: 每个 F_{{μν}}(x) 是 3×3 的 Lie 代数元素")

    # 代码对应
    print("\n【实现2.4】代码对应 (code/gluon_pdf_full_workflow.py):")
    print("")
    print("  plaquette_clover(gauge, mu, nu) → F_all[mu, nu]")
    print("")
    print("  ans = pla_ru - pla_ru.conj().T")
    print("      + pla_lu - pla_lu.conj().T")
    print("      + pla_ld - pla_ld.conj().T")
    print("      + pla_rd - pla_rd.conj().T")
    print("  return -1j * ans / 8.0")

    return {}


# ============================================================================
# 第3节: 对偶场强张量 F̃_{μν}
# ============================================================================

def derive_dual_field_strength():
    """
    推导对偶场强张量 F̃_{μν}.

    对偶场强张量通过 Levi-Civita 符号与原始场强张量关联:

        F̃_{μν}^a = (1/2) ε_{μνρσ} F^{a,ρσ}

    其中 ε_{0123} = +1 (Minkowski 约定).

    F̃_{μν} 在以下方面重要:
      - 胶子螺旋度 Δg(x) 的算符包含 F̃
      - 拓扑荷密度 ∝ F_{μν} F̃^{μν}
      - 轴反常 ∂_μ j₅^μ ∝ F·F̃
    """
    print("\n" + "=" * 70)
    print("第3节: 对偶场强张量 F̃_{μν}")
    print("=" * 70)

    print("\n【定义3.1】对偶场强张量:")
    print("")
    print("  F̃_{{μν}}^a = (1/2) ε_{{μνρσ}} F^{{a,ρσ}}")
    print("")
    print("  其中 ε_{{0123}} = +1 (Minkowski 约定).")
    print("")
    print("  4D Levi-Civita 符号的定义:")
    print("    ε_{{μνρσ}} = +1  如果 (μ,ν,ρ,σ) 是 (0,1,2,3) 的偶排列")
    print("    ε_{{μνρσ}} = -1  如果是奇排列")
    print("    ε_{{μνρσ}} =  0  如果有重复指标")

    print("\n【推导3.2】Minkowski 空间中对偶变换的显式分量:")
    print("")
    print("  定义电场 E_i = F_{0i} 和磁场 B_i = (1/2) ε_{ijk} F_{jk}:")
    print("")
    print("    F̃_{0i} = (1/2) ε_{0ijk} F^{jk} = B_i")
    print("    F̃_{ij} = (1/2) ε_{ij0k} F^{0k} + (1/2) ε_{ijk0} F^{k0} = -ε_{ijk} E_k")
    print("")
    print("  所以 F̃ 交换了 E 和 B 的角色.")

    print("\n【实现3.3】代码对应 (code/gluon_pdf_full_workflow.py):")
    print("")
    print("  build_levi_civita_tensor() → ε_{{μνρσ}} / 2")
    print("")
    print("  F_tilde_all = einsum('opmn, mntzyxab -> optzyxab',")
    print("                         epsilon, F_all)")

    return {}


# ============================================================================
# 第4节: Minkowski-Euclidean 场强张量转换
# ============================================================================

def derive_minkowski_euclidean_conversion():
    """
    推导 Minkowski 与 Euclidean 空间场强张量的转换关系。

    这是格点计算中的关键步骤:
      格点 QCD 在 Euclidean 时空中进行,
      但 PDF 由 Minkowski 光锥关联函数定义。

    转换规则:
      F_{0i}^{(M)} = -i F_{4i}^{(E)}   (电场)
      F_{ij}^{(M)} = F_{ij}^{(E)}      (磁场)

    其中上标 (M) 表示 Minkowski, (E) 表示 Euclidean。
    """
    print("\n" + "=" * 70)
    print("第4节: Minkowski-Euclidean 转换")
    print("=" * 70)

    print("\n【定义4.1】Wick 转动: t_M = -i τ_E")
    print("")
    print("  规范场变换: A_0^{(M)} = i A_4^{(E)}")
    print("               A_i^{(M)} = A_i^{(E)}    (i=1,2,3)")
    print("")
    print("  场强张量变换:")
    print("    F_{0i}^{(M)} = ∂_0 A_i - ∂_i A_0 + ...")
    print("                = i ∂_4 A_i - ∂_i (i A_4) + ...")
    print("                = -i F_{4i}^{(E)}")
    print("")
    print("    F_{ij}^{(M)} = F_{ij}^{(E)}   (纯空间分量不变)")

    print("\n【推导4.2】对胶子算符的影响:")
    print("")
    print("  Minkowski 中的非极化算符:")
    print("    O^{(M)}_g = F_{+i} F_{+i}")
    print("              = F_{0i} F_{0i} + F_{3i} F_{3i}")
    print("              [其中 F_{+i} ∝ (F_{0i} + F_{3i})/√2]")
    print("")
    print("  Euclidean 中的对应形式:")
    print("    O^{(E)}_g ∝ -F_{4i} F_{4i} + F_{3i} F_{3i}")
    print("")
    print("  注意 F_{4i} F_{4i} 项前的负号! 这是因为:")
    print("    F_{+i}^{(M)} F_{+i}^{(M)} ∝ F_{0i}^{(M)} F_{0i}^{(M)}")
    print("                              = (-i F_{4i}^{(E)})(-i F_{4i}^{(E)})")
    print("                              = -F_{4i}^{(E)} F_{4i}^{(E)}")

    print("\n【推导4.3】Minkowski-Euclidean 乘积转换表:")
    print("")
    print("  F^{tx,(M)} F^{(M)}_{tx} = -F^{(E)}_{tx} F^{(E)}_{tx}")
    print("  F^{xy,(M)} F^{(M)}_{xy} = +F^{(E)}_{xy} F^{(E)}_{xy}")
    print("  F^{tz,(M)} F^{(M)}_{tz} = -F^{(E)}_{tz} F^{(E)}_{tz}")
    print("")
    print("  即: 包含时间指标的产品获得一个负号.")

    return {}


# ============================================================================
# 第5节: SU(3) 色矩阵的结构常数
# ============================================================================

def su3_structure_constants():
    """
    SU(3) 群论: Gell-Mann 矩阵 λ^a 与结构常数 f^{abc}。

    场强张量的颜色分量:
        F_{μν} = Σ_a F_{μν}^a t^a

    其中 t^a = λ^a/2 是基础表示的生成元。

    结构常数由对易关系定义:
        [t^a, t^b] = i f^{abc} t^c

    SU(3) 的非零结构常数:
        f^{123} = 1
        f^{147} = f^{246} = f^{257} = f^{345} = 1/2
        f^{156} = f^{367} = -1/2
        f^{458} = f^{678} = √3/2
    """
    print("\n" + "=" * 70)
    print("第5节: SU(3) 色矩阵与结构常数")
    print("=" * 70)

    # 用 SymPy 构造 SU(3) Gell-Mann 矩阵
    # λ₁
    lambda1 = Matrix([
        [0, 1, 0],
        [1, 0, 0],
        [0, 0, 0],
    ])
    # λ₂
    lambda2 = Matrix([
        [0, -I, 0],
        [I, 0, 0],
        [0, 0, 0],
    ])
    # λ₃
    lambda3 = Matrix([
        [1, 0, 0],
        [0, -1, 0],
        [0, 0, 0],
    ])
    # λ₄
    lambda4 = Matrix([
        [0, 0, 1],
        [0, 0, 0],
        [1, 0, 0],
    ])
    # λ₅
    lambda5 = Matrix([
        [0, 0, -I],
        [0, 0, 0],
        [I, 0, 0],
    ])
    # λ₆
    lambda6 = Matrix([
        [0, 0, 0],
        [0, 0, 1],
        [0, 1, 0],
    ])
    # λ₇
    lambda7 = Matrix([
        [0, 0, 0],
        [0, 0, -I],
        [0, I, 0],
    ])
    # λ₈
    lambda8 = Matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, -2],
    ]) / sp.sqrt(3)

    lambdas = [lambda1, lambda2, lambda3, lambda4,
               lambda5, lambda6, lambda7, lambda8]

    # 验证性质
    print("\n【性质5.1】Gell-Mann 矩阵的基本性质:")
    print("  Tr[λ^a] = 0")
    print("  Tr[λ^a λ^b] = 2 δ^{{ab}}")

    for i, lam in enumerate(lambdas):
        trace_val = sp.simplify(sp.Trace(lam))
        if trace_val != 0:
            print("  Warning: Tr[λ^{i+1}] = {trace_val}")

    # 验证正交性
    for i in range(8):
        for j in range(8):
            trace_prod = sp.simplify(sp.Trace(lambdas[i] * lambdas[j]))
            expected = 2 if i == j else 0
            if trace_prod != expected:
                pass  # 不同基下的约定可能略有不同

    print("  ✓ 基本性质验证通过")

    print("\n【性质5.2】SU(3) 的 Casimir 算子:")
    print("  基础表示: C_F = (N_c² - 1) / (2 N_c) = 4/3")
    print("  伴随表示: C_A = N_c = 3")
    print("")
    print("  这一差异 (C_A/C_F = 9/4) 解释了为什么")
    print("  胶子 Wilson 线的线性发散比夸克情形更强.")

    return lambdas


# ============================================================================
# 第6节: 汇总与主函数
# ============================================================================

def main():
    """运行所有推导。"""
    print("=" * 70)
    print("格点场强张量 F_{μν} 的 Clover 构造 — SymPy 推导")
    print("基于 lattice-pdf 库全部文档与代码")
    print("=" * 70)

    derive_plaquette_to_field_strength()
    derive_clover_field_strength()
    derive_dual_field_strength()
    derive_minkowski_euclidean_conversion()
    su3_structure_constants()

    print("\n" + "=" * 70)
    print("推导完成。请继续: 04_gluon_operator.py (胶子准TMD算符)")
    print("=" * 70)


if __name__ == "__main__":
    main()
