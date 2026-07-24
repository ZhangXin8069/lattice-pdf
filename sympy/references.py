#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
references.py — 完整参考文献数据库
================================================================================

本模块包含本项目中所有sympy推导引用的完整参考文献列表,
按主题分类组织。

参考文献分类:
  [A] LaMET 框架与因子化定理
  [B] 胶子算符与可乘性重整化
  [C] 梯度流重整化
  [D] TMD 因子化与软函数
  [E] 连续极限与格点计算
  [F] 格点方法与 Distillation
  [G] 本库补充文档
  [H] 教材

作者: 基于 lattice-pdf 库全部文档与代码
日期: 2026-07-24
================================================================================
"""

REFERENCES = {
    # =========================================================================
    # [A] LaMET 框架与因子化定理
    # =========================================================================
    "Ji2013": {
        "cite": "X. Ji, PRL 110, 262002 (2013)",
        "title": "Parton Physics on a Euclidean Lattice",
        "arxiv": "1305.1539",
        "file": "docs/Parton Physics on Euclidean Lattice.pdf",
        "key_contribution": "LaMET 框架的奠基性论文: 准PDF 概念, Lorentz boost",
    },
    "Ji2014": {
        "cite": "X. Ji, Sci. China Phys. Mech. Astron. 57, 1407 (2014)",
        "title": "Parton Physics from Large-Momentum Effective Field Theory",
        "arxiv": "1404.6680",
        "file": "docs/Parton Physics from Large-Momentum Effective Field Theory.pdf",
        "key_contribution": "LaMET 系统表述, HQET 类比, 普适性类",
    },
    "Izubuchi2018": {
        "cite": "T. Izubuchi et al., PRD 98, 056004 (2018)",
        "title": "Factorization Theorem Relating Euclidean and Light-Cone Parton Distributions",
        "arxiv": "1801.03917",
        "file": "docs/Factorization Theorem Relating Euclidean and Light-Cone Parton Distributions.pdf",
        "key_contribution": "LaMET因子化定理的严格 OPE 证明",
    },
    "Ji2017More": {
        "cite": "X. Ji, J.-H. Zhang, Y. Zhao, Nucl. Phys. B 924, 366 (2017)",
        "title": "More on Large-Momentum Effective Theory Approach to Parton Physics",
        "arxiv": "1706.07416",
        "file": "docs/More On Large-Momentum Effective Theory Approach to Parton Physics.pdf",
        "key_contribution": "解析延拓, 胶子极化, 普适性类扩展",
    },

    # =========================================================================
    # [B] 胶子算符与可乘性重整化
    # =========================================================================
    "Zhang2019": {
        "cite": "J.-H. Zhang et al., PRL 122, 142001 (2019)",
        "title": "Accessing Gluon Parton Distributions in Large Momentum Effective Theory",
        "arxiv": "1808.10824",
        "file": "docs/Accessing gluon parton distributions in large momentum effective theory.pdf",
        "key_contribution": "辅助重夸克场方法, 4个可乘性重整化构建块, 2×2混合矩阵",
    },
    "Li2018": {
        "cite": "Z.-Y. Li, Y.-Q. Ma, J.-W. Qiu, PRL 122, 062002 (2019)",
        "title": "Multiplicative Renormalizability of Quasi-Parton Operators",
        "arxiv": "1809.01836",
        "file": "docs/Multiplicative renormalizability of quasi-parton operators.pdf",
        "key_contribution": "36个胶子准算符可乘性重整化的严格证明",
    },
    "Balitsky2020": {
        "cite": "I. Balitsky, W. Morris, A. Radyushkin, PLB 808, 135621 (2020)",
        "title": "Gluon Pseudo-Distributions at Short Distances: Forward Case",
        "arxiv": "1910.13963",
        "file": "docs/Short-Distance Structure of Unpolarized Gluon Pseudodistributions.pdf",
        "key_contribution": "六个不变振幅分解, 单圈紫外/红外分离",
    },
    "Yao2023": {
        "cite": "F. Yao, Y. Ji, J.-H. Zhang, JHEP 11, 021 (2023)",
        "title": "Connecting Euclidean to Light-Cone Correlations",
        "arxiv": "2212.14415",
        "file": "docs/Connecting Euclidean to light-cone correlations...pdf",
        "key_contribution": "味单态完备 NLO 匹配核, 非极化/极化参数对比",
    },

    # =========================================================================
    # [C] 梯度流重整化
    # =========================================================================
    "Luscher2010": {
        "cite": "M. Lüscher, JHEP 08, 071 (2010)",
        "title": "Properties and Uses of the Wilson Flow in Lattice QCD",
        "arxiv": "1006.4518",
        "file": "docs/Properties_and_uses_of_the_Wilson_flow_in_lattice_QCD_Luscher_2010.pdf",
        "key_contribution": "Wilson 流的基本性质, 平滑半径",
    },
    "Luscher2011": {
        "cite": "M. Lüscher, P. Weisz, JHEP 02, 051 (2011)",
        "title": "Perturbative Analysis of the Gradient Flow in Non-Abelian Gauge Theories",
        "arxiv": "1101.0963",
        "file": "docs/Perturbative_analysis_of_the_gradient_flow_in_non-abelian_gauge_theories_Luscher_Weisz_2011.pdf",
        "key_contribution": "UV 有限性定理 (Lüscher-Weisz 定理)",
    },
    "Luscher2013": {
        "cite": "M. Lüscher, JHEP 04, 123 (2013)",
        "title": "Chiral Symmetry and the Yang-Mills Gradient Flow",
        "arxiv": "1302.5246",
        "file": "docs/Chiral_symmetry_and_the_Yang-Mills_gradient_flow_Luscher_2013.pdf",
    },
    "Luscher2014": {
        "cite": "M. Lüscher, PoS LATTICE2013, 016 (2014)",
        "title": "Future Applications of the Yang-Mills Gradient Flow in Lattice QCD",
        "arxiv": "1308.5598",
        "file": "docs/Future_applications_of_the_Yang-Mills_gradient_flow_in_lattice_QCD_Luscher_2014.pdf",
    },
    "Suzuki2013": {
        "cite": "H. Suzuki, PTEP 2013, 083B03 (2013)",
        "title": "Energy-Momentum Tensor from the Yang-Mills Gradient Flow",
        "arxiv": "1304.0533",
        "file": "docs/Energy-momentum_tensor_from_the_Yang-Mills_gradient_flow_Suzuki_2013.pdf",
        "key_contribution": "能动张量的梯度流表示, c_T(t), c_S(t) 系数",
    },
    "Monahan2017": {
        "cite": "C. Monahan, K. Orginos, JHEP 03, 116 (2017)",
        "title": "Quasi Parton Distributions and the Gradient Flow",
        "arxiv": "1612.01584",
        "file": "docs/Quasi_parton_distributions_and_the_gradient_flow_Monahan_Orginos_2017.pdf",
        "key_contribution": "梯度流用于准 PDF 算符, SFTX 匹配",
    },
    "Monahan2018": {
        "cite": "C. Monahan, PRD 97, 054507 (2018)",
        "title": "Smeared Quasidistributions in Perturbation Theory",
        "arxiv": "1710.04607",
        "file": "docs/Smeared_quasidistributions_in_perturbation_theory_Monahan_2018.pdf",
    },
    "Shindler2018": {
        "cite": "A. Shindler, PRD 99, 054505 (2019)",
        "title": "Off-Lightcone Wilson-Line Operators in Gradient Flow",
        "arxiv": "1809.04474",
        "file": "docs/Off-lightcone_Wilson-line_operators_in_gradient_flow.pdf",
    },

    # =========================================================================
    # [D] TMD 因子化与软函数
    # =========================================================================
    "He2024": {
        "cite": "J.-C. He et al. (LPC), PRD 109, 114513 (2024)",
        "title": "Unpolarized TMD Parton Distributions of the Nucleon from Lattice QCD",
        "arxiv": "2211.02340",
        "file": "docs/Unpolarized Transverse-Momentum-Dependent Parton Distributions of the Nucleon from Lattice QCD.pdf",
        "key_contribution": "核子非极化TMD-PDF首次格点计算",
    },
    "Zhang2022renormTMD": {
        "cite": "K. Zhang et al. (LPC), PRD 106, 094510 (2022)",
        "title": "Renormalization of TMD Parton Distribution on the Lattice",
        "arxiv": "2205.13402",
        "file": "docs/Renormalization of transverse-momentum-dependent parton distribution on the lattice.pdf",
        "key_contribution": "Wilson圈减除方案, 5格距系统验证, RI/MOM失败分析",
    },
    "Chu2023soft": {
        "cite": "M.-H. Chu et al. (LPC), JHEP 08, 172 (2023)",
        "title": "Lattice Calculation of the Intrinsic Soft Function and the Collins-Soper Kernel",
        "arxiv": "2306.06488",
        "file": "docs/Lattice Calculation of the Intrinsic Soft Function and the Collins-Soper Kernel.pdf",
        "key_contribution": "内禀软函数 S_I 和 CS 核的格点提取",
    },
    "Chu2023TMDWF": {
        "cite": "M.-H. Chu et al. (LPC), PRD 108, 094513 (2023)",
        "title": "Transverse-Momentum-Dependent Wave Functions of Pion from Lattice QCD",
        "arxiv": "2302.09961",
        "file": "docs/Transverse-Momentum-Dependent Wave Functions of Pion from Lattice QCD.pdf",
    },
    "LPC2020soft": {
        "cite": "Q.-A. Zhang et al. (LPC), PRL 125, 192001 (2020)",
        "title": "Lattice-QCD Calculations of TMD Soft Function Through LaMET",
        "arxiv": "2005.14572",
        "file": "docs/Lattice-QCD Calculations of TMD Soft Function Through Large-Momentum Effective Theory.pdf",
    },
    "Chu2022CS": {
        "cite": "M.-H. Chu et al. (LPC), PRD 106, 034509 (2022)",
        "title": "Nonperturbative Determination of Collins-Soper Kernel from Quasi TMD Wave Functions",
        "arxiv": "2204.00200",
        "file": "docs/Nonperturbative Determination of Collins-Soper Kernel from Quasi Transverse-Momentum Dependent Wave Functions.pdf",
    },
    "Ma2025BoerMulders": {
        "cite": "L. Ma et al. (LPC), (2025)",
        "title": "Quark Transverse Spin-Momentum Correlation...The Boer-Mulders Function",
        "arxiv": "2502.11807",
        "file": "docs/Quark Transverse Spin-Momentum Correlation of the Nucleon from Lattice QCD The Boer-Mulders Function.pdf",
    },

    # =========================================================================
    # [E] 连续极限与格点计算 (胶子PDF)
    # =========================================================================
    "Chen2025gluon": {
        "cite": "C. Chen et al. (CLQCD, LPC), arXiv:2510.26425 (2025)",
        "title": "Unpolarized Gluon PDF of the Nucleon from Lattice QCD in the Continuum Limit",
        "file": "docs/Unpolarized gluon PDF of the nucleon from lattice QCD in the continuum limit.pdf",
        "key_contribution": "首次三格距连续极限胶子PDF, 联合外推",
    },
    "Chen2025first": {
        "cite": "C. Chen et al., PRD 111, 074506 (2025)",
        "title": "First Self-Renormalized Gluon PDF of Nucleon from LaMET in the Continuum Limit",
        "arxiv": "2408.12819",
        "file": "docs/First Self-Renormalized Gluon PDF of Nucleon from Large-Momentum Effective Theory in the Continuum Limit.pdf",
    },
    "NieMiera2025": {
        "cite": "A. NieMiera et al., arXiv:2510.17758 (2025)",
        "title": "First Self-Renormalized Gluon PDF...Continuum Limit",
        "file": "docs/First Self-Renormalized Gluon PDF of Nucleon from Large-Momentum Effective Theory in the Continuum Limit.pdf",
        "key_contribution": "MSULat 梯度流自重整化, 三格距连续极限",
    },
    "Fan2018gluon": {
        "cite": "Z.-Y. Fan et al., PRL 121, 242001 (2018)",
        "title": "Gluon Quasi-PDF from Lattice QCD",
        "arxiv": "1808.02077",
        "file": "docs/Gluon Quasi-PDF From Lattice QCD.pdf",
        "key_contribution": "首次胶子准PDF格点计算 (探路性)",
    },
    "Fan2023": {
        "cite": "Z. Fan, W. Good, H.-W. Lin, PRD 108, 014508 (2023)",
        "title": "Physical-Continuum Limit of the Nucleon Gluon Parton Distribution from Lattice QCD",
        "arxiv": "2210.09985",
        "file": "docs/Physical-Continuum Limit of the Nucleon Gluon Parton Distribution from Lattice QCD.pdf",
    },
    "Good2025a": {
        "cite": "W. Good, F. Yao, H.-W. Lin, arXiv:2505.13321 (2025)",
        "title": "First Nucleon Gluon PDF from Large Momentum Effective Theory",
        "file": "docs/First Nucleon Gluon PDF from Large Momentum Effective Theory.pdf",
    },
    "Lin2025gluon": {
        "cite": "W. Good, K. Hasan, H.-W. Lin, J. Phys. G 52, 035105 (2025)",
        "title": "Gluon Parton Distribution of the Nucleon from 2+1+1-Flavor Lattice QCD in the Physical-Continuum Limit",
        "arxiv": "2409.02750",
        "file": "docs/Gluon Parton Distribution of the Nucleon from 2+1+1-Flavor Lattice QCD in the Physical-Continuum Limit.pdf",
    },
    "Egerer2022": {
        "cite": "C. Egerer et al. (HadStruc), PRD 106, 094511 (2022)",
        "title": "Towards the Determination of the Gluon Helicity Distribution in the Nucleon from Lattice QCD",
        "arxiv": "2207.08733",
        "file": "docs/Towards the determination of the gluon helicity distribution in the nucleon from lattice quantum chromodynamics.pdf",
    },

    # =========================================================================
    # [F] 重整化方法
    # =========================================================================
    "Ji2021hybrid": {
        "cite": "X. Ji et al., Nucl. Phys. B 964, 115311 (2021)",
        "title": "A Hybrid Renormalization Scheme for Quasi Light-Front Correlations in LaMET",
        "arxiv": "2008.03886",
        "file": "docs/A Hybrid Renormalization Scheme for Quasi Light-Front Correlations in Large-Momentum Effective Theory.pdf",
        "key_contribution": "混合重整化方案 (ratio + 自重整化)",
    },
    "Huo2021self": {
        "cite": "Y.-K. Huo et al. (LPC), Nucl. Phys. B 969, 115443 (2021)",
        "title": "Self-Renormalization of Quasi-Light-Front Correlators on the Lattice",
        "arxiv": "2103.02965",
        "file": "docs/Self-Renormalization of Quasi-Light-Front Correlators on the Lattice.pdf",
        "key_contribution": "自重整化方案: z-比值消除线性发散和重整子",
    },
    "Zhang2023RGR": {
        "cite": "R. Zhang et al., Phys. Lett. B 844, 138081 (2023)",
        "title": "Leading Power Accuracy in Lattice Calculations of Parton Distributions",
        "arxiv": "2305.05212",
        "file": "docs/Leading Power Accuracy in Lattice Calculations of Parton Distributions.pdf",
    },
    "Su2023LRR": {
        "cite": "Y. Su et al., Nucl. Phys. B 991, 116201 (2023)",
        "title": "Power Corrections and Renormalons in Parton Quasi-Distributions",
        "arxiv": "2209.01236",
        "file": "docs/Power corrections and renormalons in parton quasi-distributions.pdf",
    },
    "Liu2019isovector": {
        "cite": "Y.-S. Liu et al. (LPC), PRD 101, 034020 (2020)",
        "title": "Unpolarized Isovector Quark Distribution Function from Lattice QCD",
        "arxiv": "1807.06566",
        "file": "docs/Unpolarized isovector quark distribution function from Lattice QCD A systematic analysis of renormalization and matching.pdf",
    },
    "Zhang2024gauge_fixing": {
        "cite": "K. Zhang et al. (LPC), PRD 110, 074505 (2024)",
        "title": "Impact of Gauge Fixing Precision on the Continuum Limit of Non-Local Quark-Bilinear Lattice Operators",
        "arxiv": "2405.14097",
        "file": "docs/Impact of gauge fixing precision on the continuum limit of non-local quark-bilinear lattice operators.pdf",
    },

    # =========================================================================
    # [G] 格点方法与 Distillation
    # =========================================================================
    "Peardon2009": {
        "cite": "M. Peardon et al. (Hadron Spectrum), PRD 80, 054506 (2009)",
        "title": "A Novel Quark-Field Creation Operator Construction for Hadronic Physics in Lattice QCD",
        "arxiv": "0905.2160",
        "file": "docs/A novel quark-field creation operator construction for hadronic physics in lattice QCD.pdf",
        "key_contribution": "Distillation 方法的原始提出",
    },
    "Egerer2021a": {
        "cite": "C. Egerer et al., PRD 103, 034502 (2021)",
        "title": "Distillation at High-Momentum",
        "arxiv": "2009.10691",
        "file": "docs/Distillation at High-Momentum.pdf",
        "key_contribution": "动量涂抹蒸馏, 高动量核子态",
    },
    "Zhang2025kinematic": {
        "cite": "R. Zhang et al., arXiv:2501.00729 (2025)",
        "title": "Kinematically-Enhanced Interpolating Operators for Boosted Hadrons",
        "file": "docs/Kinematically-enhanced interpolating operators for boosted hadrons.pdf",
    },

    # =========================================================================
    # [H] 组态生成 (CLQCD)
    # =========================================================================
    "CLQCD2024": {
        "cite": "Z.-C. Hu et al. (CLQCD), PRD 109, 054507 (2024)",
        "title": "Charmed Meson Masses and Decay Constants in the Continuum",
        "arxiv": "2310.00814",
        "file": "docs/Charmed meson masses and decay constants in the continuum from the tadpole improved clover ensembles.pdf",
    },
    "CLQCD2025": {
        "cite": "H.-Y. Du et al. (CLQCD), PRD 111, 054504 (2025)",
        "title": "Quark Masses and Low Energy Constants in the Continuum",
        "arxiv": "2408.03548",
        "file": "docs/Quark masses and low energy constants in the continuum from the tadpole improved clover ensembles.pdf",
    },

    # =========================================================================
    # [I] 本库内部文档
    # =========================================================================
    "NoteGluonPDFs": {
        "cite": "内部笔记",
        "title": "Note of Gluon PDFs",
        "file": "docs/Note of gluon PDFs.pdf, 文档/note_of_gluon_PDFs.tex",
        "key_contribution": "Eq.(20)矩阵元, Eq.(25)OPE算符, 非极化/极化胶子流构造",
    },
    "GluonPDFDerivation": {
        "cite": "内部文档",
        "title": "胶子PDF推导 (gluon_pdf_derivation.tex)",
        "file": "文档/gluon_pdf_derivation.tex",
        "key_contribution": "Jaffe-Manohar分解, 完整格点推导",
    },
    "GluonPDFContinuum": {
        "cite": "内部文档",
        "title": "胶子PDF连续极限 (gluon_PDF_continuum.tex)",
        "file": "文档/gluon_PDF_continuum.tex",
        "key_contribution": "三格距联合外推细节, NLO匹配核完整形式",
    },
    "SupplementGradientFlow": {
        "cite": "补充文档",
        "title": "格点QCD中的梯度流重整化",
        "file": "补充/格点QCD中的梯度流重整化.tex",
    },
    "SupplementTMD": {
        "cite": "补充文档",
        "title": "格点QCD中的TMD_PDF",
        "file": "补充/格点QCD中的TMD_PDF.tex",
    },
    "SupplementExtrapolation": {
        "cite": "补充文档",
        "title": "格点QCD中的外推",
        "file": "补充/格点QCD中的外推.tex",
    },
    "SupplementRenormalization": {
        "cite": "补充文档",
        "title": "格点QCD中的重整化",
        "file": "补充/格点QCD中的重整化.tex",
    },
    "SupplementGluonOperator": {
        "cite": "补充文档",
        "title": "格点QCD中的胶子算符",
        "file": "补充/格点QCD中的胶子算符.tex",
    },
    "SupplementPDF": {
        "cite": "补充文档",
        "title": "格点QCD中的部分子分布函数",
        "file": "补充/格点QCD中的部分子分布函数.tex",
    },
    "SupplementWilsonLine": {
        "cite": "补充文档",
        "title": "格点QCD中的Wilson线",
        "file": "补充/格点QCD中的Wilson线.tex",
    },
    "SupplementFieldStrength": {
        "cite": "补充文档",
        "title": "格点QCD中的场强张量",
        "file": "补充/格点QCD中的场强张量.tex",
    },
    "SupplementLaMET": {
        "cite": "补充文档",
        "title": "格点QCD中的大动量有效理论",
        "file": "补充/格点QCD中的大动量有效理论.tex",
    },
    "SupplementGluonPolarization": {
        "cite": "补充文档",
        "title": "格点QCD中的胶子极化",
        "file": "补充/格点QCD中的胶子极化.tex",
    },
    "ReportGluonQuasiOp": {
        "cite": "汇报文档",
        "title": "格点上计算胶子准算符",
        "file": "汇报/格点上计算胶子准算符.tex",
    },
    "ReportConstructGluonOp": {
        "cite": "汇报文档",
        "title": "构造胶子准算符",
        "file": "汇报/构造胶子准算符.tex",
    },
    "ReferTheoryWorkflow": {
        "cite": "参考文档",
        "title": "理论解析与工作流",
        "file": "refer/理论解析与工作流.tex",
    },

    # =========================================================================
    # [J] 教材
    # =========================================================================
    "GattringerLang": {
        "cite": "C. Gattringer, C. B. Lang, Lect. Notes Phys. 788, Springer (2010)",
        "title": "Quantum Chromodynamics on the Lattice",
        "file": "books/Quantum Chromodynamics on the Lattice.pdf",
    },
    "Gupta": {
        "cite": "R. Gupta, arXiv:hep-lat/9807028 (1997)",
        "title": "Introduction to Lattice QCD",
        "file": "books/INTRODUCTION TO LATTICE QCD.pdf",
    },
    "PeskinSchroeder": {
        "cite": "M. E. Peskin, D. V. Schroeder, Westview Press (1995)",
        "title": "An Introduction to Quantum Field Theory",
        "file": "books/An Introduction to Quantum Field Theory.pdf",
    },

    # =========================================================================
    # [K] 代码
    # =========================================================================
    "CodeFullWorkflow": {
        "cite": "code/gluon_pdf_full_workflow.py",
        "title": "胶子PDF完整工作流 (~1900行)",
        "key_contribution": "10步管线单脚本实现",
    },
    "CodeWorkflow": {
        "cite": "code/gluon_pdf_workflow.py",
        "title": "胶子PDF工作流 (生产级, ~1580行)",
        "key_contribution": "子命令架构, MPI支持, 系综预设",
    },
    "ExamplesOperators": {
        "cite": "examples/Operator.py, examples/gamma_matrix_cupy_DR.py",
        "title": "GPU算符实现 (donghx)",
        "key_contribution": "Plaquette, Clover F, DR基Dirac矩阵 (CuPy)",
    },
}


def print_references_by_category():
    """按类别打印参考文献。"""
    categories = {
        "[A] LaMET框架": ["Ji2013", "Ji2014", "Izubuchi2018", "Ji2017More"],
        "[B] 胶子算符": ["Zhang2019", "Li2018", "Balitsky2020", "Yao2023"],
        "[C] 梯度流": ["Luscher2010", "Luscher2011", "Luscher2013",
                       "Luscher2014", "Suzuki2013", "Monahan2017",
                       "Monahan2018", "Shindler2018"],
        "[D] TMD/软函数": ["He2024", "Zhang2022renormTMD", "Chu2023soft",
                          "Chu2023TMDWF", "LPC2020soft", "Chu2022CS",
                          "Ma2025BoerMulders"],
        "[E] 胶子PDF格点": ["Chen2025gluon", "Chen2025first", "NieMiera2025",
                           "Fan2018gluon", "Fan2023", "Good2025a",
                           "Lin2025gluon", "Egerer2022"],
        "[F] 重整化方法": ["Ji2021hybrid", "Huo2021self", "Zhang2023RGR",
                          "Su2023LRR", "Liu2019isovector",
                          "Zhang2024gauge_fixing"],
        "[G] 格点方法": ["Peardon2009", "Egerer2021a", "Zhang2025kinematic"],
        "[H] 组态": ["CLQCD2024", "CLQCD2025"],
        "[I] 内部文档": ["NoteGluonPDFs", "GluonPDFDerivation",
                        "GluonPDFContinuum", "SupplementGradientFlow",
                        "SupplementTMD", "SupplementExtrapolation",
                        "SupplementRenormalization", "SupplementGluonOperator",
                        "SupplementPDF", "SupplementWilsonLine",
                        "SupplementFieldStrength", "SupplementLaMET",
                        "SupplementGluonPolarization", "ReportGluonQuasiOp",
                        "ReportConstructGluonOp", "ReferTheoryWorkflow"],
        "[J] 教材": ["GattringerLang", "Gupta", "PeskinSchroeder"],
        "[K] 代码": ["CodeFullWorkflow", "CodeWorkflow", "ExamplesOperators"],
    }

    for cat_name, keys in categories.items():
        print("\n" + "─" * 70)
        print("  " + cat_name)
        print("─" * 70)
        for key in keys:
            ref = REFERENCES[key]
            print("  [%s] %s" % (key, ref['cite']))
            print("         %s" % ref['title'])
            if 'arxiv' in ref:
                print("         arXiv:%s" % ref['arxiv'])
            if 'file' in ref:
                print("         文件: %s" % ref['file'])


if __name__ == "__main__":
    print("=" * 70)
    print("完整参考文献数据库")
    print("共 {len(REFERENCES)} 条参考文献")
    print("=" * 70)
    print_references_by_category()
