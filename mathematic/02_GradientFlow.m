(* ::Package:: *)

(* ============================================================================
   第2部分: 梯度流 (Gradient Flow / Wilson Flow) 重整化理论
   ============================================================================

   参考文献:
   - Luescher, "Properties and uses of the Wilson flow in lattice QCD",
     JHEP 08, 071 (2010)
   - Luescher & Weisz, "Perturbative analysis of the gradient flow in
     non-abelian gauge theories", JHEP 02, 051 (2011)
   - Luescher, "Chiral symmetry and the Yang-Mills gradient flow",
     JHEP 04, 123 (2013)
   - Suzuki, "Energy-momentum tensor from the Yang-Mills gradient flow",
     PTEP 2013, 083B03 (2013)
   - Monahan & Orginos, "Quasi parton distributions and the gradient flow",
     JHEP 03, 116 (2017)
   - 补充/格点QCD中的梯度流重整化.tex
   ============================================================================ *)

BeginPackage["GradientFlow`"];

GradientFlowDerivation;
WilsonFlowEquation;
FlowTimeDependentFieldStrength;
GradientFlowUVFiniteness;
ShortFlowTimeExpansion;
GradientFlowRenormalizedGluonOperator;

Begin["`Private`"];

(* --------------------------------------------------------------------------
   2.1 Wilson流方程 (纯规范场)
   --------------------------------------------------------------------------
   规范场B_mu(t,x)沿流时间t的演化方程:
   \partial_t B_mu(t,x) = D_nu G_{nu,mu}(t,x),   B_mu(0,x) = A_mu(x)

   其中流时间依赖的场强张量和协变导数:
   G_{mu,nu}(t,x) = \partial_mu B_nu - \partial_nu B_mu + [B_mu, B_nu]
   D_mu = \partial_mu + [B_mu, ·]

   格点上的Wilson流:
   \dot{V}_t(x,mu) = -g_0^2 {\partial_{x,mu} S_W(V_t)} V_t(x,mu),  V_0 = U
   -------------------------------------------------------------------------- *)

WilsonFlowEquation[verbose_:False] := Module[{flow},
  flow = {
    "ContinuumEquation" ->
      "\\partial_t B_\\mu(t,x) = D_\\nu G_{\\nu\\mu}(t,x),  B_\\mu(0,x) = A_\\mu(x)",
    "FlowFieldStrength" ->
      "G_{\\mu\\nu}(t,x) = \\partial_\\mu B_\\nu - \\partial_\\nu B_\\mu + [B_\\mu, B_\\nu]",
    "LatticeEquation" ->
      "\\dot{V}_t(x,\\mu) = -g_0^2 \\{\\partial_{x,\\mu} S_W(V_t)\\} V_t(x,\\mu)",
    "GaugeFixing" ->
      "\\partial_t B_\\mu = D_\\nu G_{\\nu\\mu} + \\alpha_0 D_\\mu \\partial_\\nu B_\\nu  (\\alpha_0=1: Feynman)"
  };

  If[verbose,
    Print["=== Wilson流方程 ==="];
    Print["连续形式:"];
    Print["  \\partial_t B_\\mu(t,x) = D_\\nu G_{\\nu\\mu}(t,x)"];
    Print["  B_\\mu(0,x) = A_\\mu(x)  (初始条件)"];
    Print["格点形式 (Wilson flow):"];
    Print["  \\dot{V}_t(x,\\mu) = -g_0^2 {\\partial_{x,\\mu} S_W(V_t)} V_t(x,\\mu)"];
  ];

  flow
];

(* --------------------------------------------------------------------------
   2.2 领头阶解: 热核与平滑性质
   --------------------------------------------------------------------------
   在Feynman规范 (alpha_0=1) 下, 线性化的流方程可通过热核求解:

   B_{mu,1}(t,x) = \int d^D y  K_t(x-y) A_mu(y)

   热核: K_t(z) = exp(-z^2/(4t)) / (4\pi t)^{D/2}

   流场的平滑范围: RMS半径 = sqrt(8t)  (在D=4维中)

   这明确地展示了梯度流的核心性质:
   - 对于t>0, 规范场被高斯核在范围sqrt(8t)内平滑
   - 高动量模式被指数级抑制: ~exp(-t p^2)
   -------------------------------------------------------------------------- *)

HeatKernelSolution[verbose_:False] := Module[{},
  If[verbose,
    Print["=== 热核解 (领头阶) ==="];
    Print["B_{\\mu,1}(t,x) = \\int d^D y K_t(x-y) A_\\mu(y)"];
    Print["K_t(z) = e^{-z^2/(4t)} / (4\\pi t)^{D/2}"];
    Print["平滑范围: RMS半径 = \\sqrt{8t}  (D=4)"];
    Print["动量空间: \\tilde{B}_\\mu(t,p) = e^{-t p^2} \\tilde{A}_\\mu(p)"];
    Print["高动量抑制因子: e^{-t p^2}"];
  ];

  {
   "LeadingOrderSolution" ->
     "B_{\\mu,1}(t,x) = \\int d^D y K_t(x-y) A_\\mu(y)",
   "HeatKernel" -> "K_t(z) = e^{-z^2/(4t)} / (4\\pi t)^{D/2}",
   "SmoothingRange" -> "\\sqrt{8t}",
   "MomentumSpace" -> "e^{-t p^2}"
  }
];

(* --------------------------------------------------------------------------
   2.3 夸克场的梯度流
   --------------------------------------------------------------------------
   夸克场沿规范协变热方程演化:
   \partial_t \chi(t,x) = \Delta \chi(t,x)
   \partial_t \bar{\chi}(t,x) = \bar{\chi}(t,x) \overleftarrow{\Delta}

   其中 \Delta = D_mu D_mu 是规范协变Laplace算符

   注意: 夸克流方程保持手征对称性 - 流夸克场在手征转动下的变换
   方式与基本夸克场完全相同。

   夸克传播子 (领头阶):
   <\chi(t,x) \bar{\chi}(s,y)> = \int d^Dp/(2\pi)^D e^{ip(x-y)} e^{-(t+s)p^2}/(m_0+i\slashed{p})
   指数因子 e^{-(t+s)p^2} 对t+s>0保证了UV有限性
   -------------------------------------------------------------------------- *)

QuarkGradientFlow[verbose_:False] := Module[{},
  If[verbose,
    Print["=== 夸克场梯度流 ==="];
    Print["\\partial_t \\chi = D_\\mu D_\\mu \\chi,  \\chi(0) = \\psi"];
    Print["手征对称性保持: \\chi 与 \\psi 手征变换方式相同"];
    Print["传播子: e^{-(t+s)p^2}/(m_0+i\\slashed{p})"];
    Print["t+s>0时无奇异性"];
  ];

  {
   "QuarkFlowEquation" ->
     "\\partial_t \\chi = D_\\mu D_\\mu \\chi, \\partial_t \\bar{\\chi} = \\bar{\\chi} \\overleftarrow{D_\\mu D_\\mu}",
   "ChiralSymmetry" -> "保持手征对称性",
   "Propagator" -> "e^{-(t+s)p^2}/[m_0 + i\\slashed{p}]  + O(g_0^2)"
  }
];

(* --------------------------------------------------------------------------
   2.4 梯度流的核心定理: UV有限性
   --------------------------------------------------------------------------
   Luescher-Weisz定理 (2011):

   对于由梯度流在流时间 t>0 产生的规范场 B_mu(t,x), 其任意规范不变
   局域乘积的关联函数在标准QCD参数 (耦合常数和质量) 重整化后是UV有限的。
   流场 B_mu(t,x) 本身不需要任何波函数重整化。

   证明思路:
   1. 梯度流关联函数等价于一个 D+1 维局域场论的关联函数
      (额外维度 = 流时间 t)
   2. 在 D+1 维理论中, 体传播子包含因子 e^{-t p^2}
   3. 对于 t>0, 所有涉及体传播子的Feynman图都UV收敛
   4. 唯一的发散来自 t=0 边界 - 即标准4维QCD发散
   5. BRS对称性保证体部分不接受任何无穷大的抵消项

   包含夸克场时: 流夸克场需要乘性重整化因子 Z_chi,
   但 Z_chi 不依赖于流时间t。

   对胶子算符的应用: 纯胶子算符 E_t(x) = 1/4 G_{mu,nu}^a G_{mu,nu}^a
   根本不需要重整化! (n=nbar=0)
   -------------------------------------------------------------------------- *)

GradientFlowUVFiniteness[verbose_:False] := Module[{},
  If[verbose,
    Print["=== 梯度流的UV有限性 ==="];
    Print["Luescher-Weisz定理 (2011):"];
    Print["  对于t>0, 规范不变算符的关联函数在"];
    Print["  标准QCD参数重整化后是UV有限的"];
    Print[""];
    Print["关键点:"];
    Print["  1. 流规范场 B_mu 不需要波函数重整化"];
    Print["  2. 流夸克场需要 Z_chi, 但Z_chi不依赖于t"];
    Print["  3. 纯胶子算符 (如 E_t) 完全不需要重整化"];
    Print[""];
    Print["D+1维场论证明:"];
    Print["  额外维度 = 流时间t (质量量纲 -2)"];
    Print["  体传播子 ∝ e^{-t p^2} / p^2"];
    Print["  t>0时所有体图UV收敛"];
    Print["  仅t=0边界需要标准QCD抵消项"];
  ];

  {
   "Theorem" -> "Luescher-Weisz (2011): 梯度流关联函数在t>0处UV有限",
   "GaugeField" -> "B_\\mu 不需要波函数重整化",
   "QuarkField" -> "\\chi 需要 Z_\\chi, Z_\\chi不依赖于t",
   "GluonOperator" -> "纯胶子算符 (n=nbar=0) 完全不需要重整化",
   "ProofMethod" -> "D+1维局域场论 + 幂次计数 + BRS对称性"
  }
];

(* --------------------------------------------------------------------------
   2.5 小流时展开 (Short Flow-Time Expansion, SFTX)
   --------------------------------------------------------------------------
   对于任意规范不变的局域算符 O(t,x), 在 t->0 极限下:

   O(t,x) ~ \sum_k c_k(t) O_{R,k}(x)

   其中 O_{R,k}(x) 是 t=0 处重整化的局域算符,
   c_k(t) 是有限的展开系数 (Wilson系数)。

   对规范作用量密度:
   E(t,x) ~ <E(t,x)> + c_E(t) [1/4 F_{rho,sigma}^a F_{rho,sigma}^a]_R(x) + O(t)

   对能动张量 (Suzuki 2013):
   T_{mu,nu}_R(x) = lim_{t->0} [1/c_T(t) U_{mu,nu}(t,x)
                    - c_S(t)/(c_T(t)c_E(t)) delta_{mu,nu} (E(t,x) - <E(t,x)>)]
   -------------------------------------------------------------------------- *)

ShortFlowTimeExpansion[verbose_:False] := Module[{},
  If[verbose,
    Print["=== 小流时展开 (SFTX) ==="];
    Print["O(t,x) ~ \\sum_k c_k(t) O_{R,k}(x)   (t -> 0)"];
    Print[""];
    Print["能动张量梯度流公式 (Suzuki 2013):"];
    Print["  T_{mu,nu}_R(x) = lim_{t->0} [1/c_T(t) U_{mu,nu}(t,x)"];
    Print["    - c_S/(c_T c_E) delta_{mu,nu} (E - <E>)]"];
    Print[""];
    Print["系数渐近行为 (t->0+):"];
    Print["  1/c_T(t) ~ -2 b_0 ln(\\sqrt{8t}\\Lambda) + c_1"];
  ];

  {
   "SFTX" -> "O(t,x) ~ \\sum_k c_k(t) O_{R,k}(x)",
   "EMT_Formula" -> "Suzuki (2013): 能动张量的梯度流表示",
   "AsymptoticCoeffs" -> "1/c_T(t) ~ -2b_0 ln(\\sqrt{8t}\\Lambda)"
  }
];

(* --------------------------------------------------------------------------
   2.6 梯度流在部分子分布函数中的应用
   --------------------------------------------------------------------------
   Monahan & Orginos (2017) 的关键思想:

   对夸克场和规范场同时施加梯度流, 并将流时间 tau 固定在物理单位中,
   则所有涉及涂抹后Wilson线的矩阵元在连续极限下都是有限的。

   涂抹后的准PDF矩阵元:
   h^{(s)}(z/sqrt{tau}, sqrt{tau} P_z, sqrt{tau} Lambda_QCD, sqrt{tau} M_N)

   关键优势:
   1. 梯度流作为规范不变的UV正规化, 替代格距a的角色
   2. 涂抹后矩阵元在a->0 (保持sqrt{tau}固定) 下有限
   3. 转动对称性保持: twist-2算符之间仅对数混合, 无幂次混合
   4. 在连续理论中独立于格点细节进行匹配

   对胶子PDF的实际应用:
   - CLQCD/LPC: HYP5涂抹 (等价于离散梯度流步骤)
   - MSULat: 梯度流涂抹 (tau_W = 3 a^2)
   - 两种方法给出定性一致的结果

   梯度流中胶子算符的自重整化 (MSULat 2025):
   使用梯度流涂抹后自重整化保持稳定:
   k = 0.61(27), Lambda_QCD = 0.286(85) GeV, m_0 = 0.13(30) GeV
   -------------------------------------------------------------------------- *)

GradientFlowRenormalizedGluonOperator[verbose_:False] := Module[{},
  If[verbose,
    Print["=== 梯度流在胶子PDF计算中的应用 ==="];
    Print["Monahan-Orginos (2017):"];
    Print["  梯度流正规化准PDF算符中的Wilson线"];
    Print["  涂抹后矩阵元在a->0下有限"];
    Print[""];
    Print["HYP涂抹 vs 梯度流 (内部笔记对比):"];
    Print["  HYP5涂抹 与 Wilson flow (tau=1.4) 定性一致"];
    Print[""];
    Print["自重整化 + 梯度流 (MSULat 2025):"];
    Print["  k = 0.61(27), Lambda_QCD = 0.286(85) GeV"];
    Print["  m_0 = 0.13(30) GeV"];
  ];

  {
   "MonahanOrginos2017" ->
     "梯度流作为UV正规化, 替代格距a的角色",
   "KeyProperty" ->
     "涂抹后矩阵元在连续极限a->0 (保持sqrt{tau}固定) 下有限",
   "PracticalImplementations" -> {
     "CLQCD/LPC: 10步HYP涂抹 (HYP5)",
     "MSULat: 梯度流涂抹 (tau_W = 3 a^2)"
   },
   "SelfRenormalizationParams" -> {
     "k" -> 0.61, "k_uncertainty" -> 0.27,
     "LambdaQCD" -> 0.286, "LambdaQCD_uncertainty" -> 0.085,
     "m0" -> 0.13, "m0_uncertainty" -> 0.30
   }
  }
];

(* --------------------------------------------------------------------------
   主推导函数
   -------------------------------------------------------------------------- *)
GradientFlowDerivation[] := Module[{},
  Print["=== Part 2: 梯度流 (Gradient Flow) 重整化 ==="];

  Print["\n2.1 Wilson流方程 (纯规范场):"];
  Print["  \\partial_t B_\\mu = D_\\nu G_{\\nu\\mu},  B_\\mu(0) = A_\\mu"];
  Print["  格点: \\dot{V}_t = -g_0^2 {\\partial S_W(V_t)} V_t"];

  Print["\n2.2 热核解 (领头阶):"];
  Print["  B_{\\mu,1}(t,x) = \\int d^D y K_t(x-y) A_\\mu(y)"];
  Print["  K_t(z) = e^{-z^2/(4t)} / (4\\pi t)^{D/2}"];
  Print["  平滑范围: RMS半径 = \\sqrt{8t}"];

  Print["\n2.3 Luescher-Weisz UV有限性定理:"];
  Print["  对于t>0, 梯度流规范不变算符 = UV有限"];
  Print["  B_\\mu 不需要波函数重整化"];
  Print["  纯胶子算符完全不需要重整化!"];

  Print["\n2.4 小流时展开 (SFTX):"];
  Print["  O(t,x) ~ \\sum_k c_k(t) O_{R,k}(x)  (t->0)"];

  Print["\n2.5 Monahan-Orginos (2017) PDF应用:"];
  Print["  梯度流正规化 = 替代格距a"];
  Print["  涂抹后矩阵元在连续极限下有限"];
  Print["  仅对数混合, 无幂次混合"];

  Print["\n2.6 参考源:"];
  Print["  - Luescher, JHEP 08, 071 (2010)"];
  Print["  - Luescher & Weisz, JHEP 02, 051 (2011)"];
  Print["  - Monahan & Orginos, JHEP 03, 116 (2017)"];
  Print["  - 补充/格点QCD中的梯度流重整化.tex"];
];

End[];
EndPackage[];
GradientFlowDerivation[];
