#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
10_full_workflow.py — 端到端完整计算工作流
================================================================================

整合前面所有模块, 构建从格点规范组态到物理胶子TMD-PDF的完整计算工作流.

10步计算管线 (对应 code/gluon_pdf_full_workflow.py):

  Step 1:  读取规范组态 (gauge links)
  Step 2:  构造场强张量 F_{μν} (Clover plaquette)
  Step 3:  构造对偶场强张量 F̃_{μν} (Levi-Civita)
  Step 4:  构造非定域胶子准TMD算符 O(z, b__perp) (Wilson线 + staple)
  Step 5:  梯度流涂抹 (Wilson flow renormalization)
  Step 6:  Distillation框架 (本征矢量 + Perambulator + VVV)
  Step 7:  质子两点关联函数 + 非连通三点关联函数
  Step 8:  矩阵元提取 h(z, b__perp, P_z)
  Step 9:  Fourier变换 → quasi-TMD → 重整化 → 匹配 → 光锥TMD
  Step 10: 连续极限外推 (a→0) + 无穷大动量外推 (P_z→∞) + 软函数

完整的公式链 (Eq. 链接):

  规范组态 U_μ(x)
    → [Step 2-3] F_{μν}(x), F̃_{μν}(x)
    → [Step 4] O(z, b__perp) = Tr[F(z) W^{staple} F̃(0) W^{staple†}]
    → [Step 5] O_flow(t) → O_R (SFTX, t→0)
    → [Step 6-7] C₂(P_z, Δt), C₃(z, b__perp, P_z, Δt)
    → [Step 8] h(z, b__perp, P_z) = C₃ / C₂
    → [Step 9] Fourier: g̃(x, k__perp, P_z)
              Matching: g(x, k__perp, μ) = C ⊗ g̃ + O(1/P_z²)
              TMD factor: f₁^g(x, k__perp²) = H · g̃ · S^{1/2}
    → [Step 10] Joint extrap: (a→0, P_z→∞) → 物理结果

参考源:
  [1-20] 全部前文所列参考文献
  - code/gluon_pdf_full_workflow.py: 完整实现
  - code/gluon_pdf_workflow.py: 生产级实现 (MPI + 系综预设)

作者: 基于 lattice-pdf 库全部文档与代码的系统化推导
日期: 2026-07-24
================================================================================
"""


def print_full_workflow():
    """打印完整的10步计算工作流。"""
    print("=" * 80)
    print("  核子非极化胶子 TMD-PDF 的格点 QCD 计算")
    print("  梯度流重整化 + 连续极限 + LaMET/TMD 因子化")
    print("  完整 10 步计算工作流")
    print("=" * 80)

    steps = [
        ("Step 1: 读取规范组态",
         "U_μ(x) ∈ SU(3), 形状 [Nt, Nz, Ny, Nx, 4, 3, 3]\n"
         "CLQCD 系综: TITLS+TITLC, a ∈ {0.105, 0.0897, 0.0775} fm\n"
         "m_π ≈ 300 MeV, N_conf > 700"),

        ("Step 2: 构造场强张量 F_{μν}",
         "Clover 构造: Q_{μν} = P_{μν} + P_{ν,-μ} + P_{-μ,-ν} + P_{-ν,μ}\n"
         "F_{μν} = -i/(8a²g₀) · (Q_{μν} - Q_{μν}†)\n"
         "O(a²) 改进, 6 个独立分量"),

        ("Step 3: 对偶场强张量 F̃_{μν}",
         "F̃_{μν} = (1/2) ε_{μνρσ} F^{ρσ}\n"
         "用于非极化算符中的 F·F̃ 组合"),

        ("Step 4: 非定域胶子准TMD算符",
         "O(z, b__perp) = Tr_c[ F_{zμ}(z, b__perp) W^{staple}(z,b__perp;0,0)\n"
         "                    F̃^{μz}(0,0) W^{staple†} ]\n"
         "可乘性重整化组合: M_{tx;tx} + M_{ty;ty} - 2M_{xy;xy}"),

        ("Step 5: 梯度流涂抹与重整化",
         "Wilson flow: ∂_t B_μ = D_ν G_{νμ}, t_flow ≈ 3a²\n"
         "SFTX: O(t) ≃ Σ_k c_k(t) O_{R,k}  (t→0)\n"
         "零流时间外推: h(t) = h(0) + c₁t + c₂t² + ...\n"
         "HYP5 涂抹 (替代方案, CLQCD/LPC 使用)"),

        ("Step 6: Distillation 框架",
         "Laplace 本征矢量: ∇² v_k = λ_k v_k, N_ev = 100\n"
         "Perambulator: τ(t_snk, t_src) = V† D-¹ V\n"
         "VVV: Φ_{abc}(P) = Σ_x e^{-iP·x} ε_{ijk} v_i^a v_j^b v_k^c"),

        ("Step 7: 关联函数",
         "质子2pt: C₂(Δt, P_z) (宇称投影 P₊)\n"
         "非连通3pt: C₃(z, b__perp, P_z, Δt)\n"
         "  = ⟨(O_g - ⟨O_g⟩)(C₂^N - ⟨C₂^N⟩)⟩\n"
         "真空期望值减除 + 多源平均"),

        ("Step 8: 矩阵元提取",
         "h(z, b__perp, P_z) = lim_{Δt→∞} C₃(Δt) / C₂(Δt)\n"
         "在 disconnect 近似: h = ⟨O⟩ (组态+时间平均)\n"
         "双态拟合: R(Δt) = c₀ + c₁ e^{-ΔE Δt}\n"
         "Jackknife 误差分析"),

        ("Step 9: Fourier变换 + 匹配 + TMD因子化",
         "Fourier: g̃(x, k__perp, P_z) = ∫ dz d²b__perp e^{ixP_zz - ik__perp·b__perp} h(z,b__perp)\n"
         "重整化: h_R = h_B / Z_R (混合方案: ratio + 自重整化)\n"
         "匹配: g(x,μ) = C_{gg} ⊗ g̃ + C_{gq} ⊗ q̃  (2×2 NLO)\n"
         "TMD: f₁^g(x,k__perp²;μ,ζ) = H(μ,ζ) × ∫db__perp J₀(k__perpb__perp) g̃(x,b__perp) × S^{1/2}"),

        ("Step 10: 连续极限外推",
         "联合外推:\n"
         "  xg(x, P_z, a) = xg_0(x) + a² f(x)\n"
         "                 + a² P_z² h(x) + d(x)/P_z²\n"
         "需要 ≥3 个 a 值 + ≥2 个 P_z 值\n"
         "最终结果: 物理胶子 TMD-PDF f₁^g(x, k__perp²; μ, ζ)"),
    ]

    for i, (title, desc) in enumerate(steps):
        print("\n{'─' * 70}")
        print("  {title}")
        print("{'─' * 70}")
        for line in desc.split('\n'):
            print("  {line}")

    print("\n{'=' * 80}")
    print("  关键数值参数 (CLQCD/LPC 2025):")
    print("    格距: a = {{0.105, 0.0897, 0.0775}} fm")
    print("    π质量: m_π ≈ 300 MeV")
    print("    核子动量: P_z ≤ 1.97 GeV  (n_z ≤ 6 × 2π/L)")
    print("    蒸馏矢量: N_ev = 100")
    print("    涂抹: HYP5 (10步)")
    print("    重整化: 混合方案 (z_s ≈ 0.3 fm)")
    print("    匹配: NLO, 2×2 味混合矩阵")
    print("    误差: Jackknife, 系统误差 4 源")
    print("{'=' * 80}")


def print_reference_chain():
    """打印完整的参考文献链。"""
    print("\n{'=' * 80}")
    print("  关键公式的参考源映射")
    print("{'=' * 80}")

    refs = [
        ("光锥胶子PDF定义",
         "Eq.(20)-type from Note of gluon PDFs; Balitsky 2020"),
        ("准胶子PDF算符",
         "Zhang et al. PRL 122, 142001 (2019); Li-Ma-Qiu PRL 122, 062002"),
        ("LaMET因子化定理",
         "Izubuchi et al. PRD 98, 056004 (2018)"),
        ("梯度流 (Wilson Flow)",
         "Lüscher JHEP 08, 071 (2010); Lüscher-Weisz JHEP 02, 051 (2011)"),
        ("Monahan-Orginos梯度流PDF",
         "Monahan & Orginos JHEP 03, 116 (2017)"),
        ("混合重整化方案",
         "Ji et al. Nucl. Phys. B 964, 115311 (2021)"),
        ("自重整化方案",
         "Huo et al. Nucl. Phys. B 969, 115443 (2021)"),
        ("连续极限胶子PDF",
         "Chen et al. (CLQCD/LPC) arXiv:2510.26425 (2025)"),
        ("核子非极化TMD",
         "He et al. (LPC) PRD 109, 114513 (2024)"),
        ("CS核 + 软函数",
         "Chu et al. JHEP 08, 172 (2023)"),
    ]

    for formula, ref in refs:
        print("  • {formula}")
        print("    → {ref}")


def main():
    print_full_workflow()
    print_reference_chain()

    print("\n{'=' * 80}")
    print("  完整推导完成。")
    print("  所有 sympy 模块位于 ./sympy/ 目录:")
    print("    01_lightcone_tmd.py    - 光锥胶子TMD-PDF定义")
    print("    02_gradient_flow.py    - 梯度流(Wilson flow)重整化")
    print("    03_field_strength.py   - 场强张量Clover构造")
    print("    04_gluon_operator.py   - 非定域胶子准TMD算符")
    print("    05_staple_wilson.py    - Staple型Wilson线")
    print("    06_factorization.py    - 因子化定理与匹配核")
    print("    07_continuum_limit.py  - 连续极限外推")
    print("    08_soft_cs_kernel.py   - 软函数与CS演化核")
    print("    09_correlation_fn.py   - 关联函数与矩阵元提取")
    print("    10_full_workflow.py    - 完整10步工作流 (本文件)")
    print("    references.py          - 完整参考文献数据库")
    print("{'=' * 80}")


if __name__ == "__main__":
    main()
