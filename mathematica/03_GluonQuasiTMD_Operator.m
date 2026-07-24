(* ::Package:: *)

(* ============================================================================
   第3部分: 胶子准TMD-PDF算符
   Gluon Quasi-TMD-PDF Operator Construction
   ============================================================================

   参考文献:
   - Balitsky, Morris, Radyushkin, PLB 808, 135621 (2020)
   - Zhang, Ji, Schaefer, Wang, Zhao, PRL 122, 142001 (2019)
   - Li, Ma, Qiu, PRL 122, 062002 (2019)
   - He et al. (LPC), PRD 109, 114513 (2024)
   - 文档/gluon_pdf_derivation.tex
   - 补充/格点QCD中的胶子算符.tex
   ============================================================================ *)

BeginPackage["GluonQuasiTMD`"];

GluonQuasiTMDDerivation;
InvariantAmplitudeDecomposition;
GluonQuasiTMDOperator;
GluonCurrentUnpolarized;
GluonCurrentHelicity;
MultiplicativeRenormalizability;

Begin["`Private`"];

(* --------------------------------------------------------------------------
   3.1 伴随表示下的胶子TMD算符定义
   --------------------------------------------------------------------------
   光锥胶子PDF (伴随表示):
   x g(x, mu) = \int dz^-/(2\pi P^+) e^{-ixP^+z^-}
                <P| F^{+mu}_a(z^-) U_{ab}(z^-,0) F^{b mu}_+(0) |P>

   其中:
   - z^{pm} = (z^0 pm z^3)/sqrt{2} 是光锥坐标
   - F^{+mu} = F^{0mu} + F^{3mu} (光锥规范下的场强张量分量)
   - U_{ab} 是伴随表示中的光锥Wilson线:
     U_{ab}(z^-,0) = P exp[i g \int_0^{z^-} d\eta^- A^+(eta^-)]_{ab}
   - A^{bc}_mu = i f^{abc} A^a_mu 是伴随表示规范场
   -------------------------------------------------------------------------- *)

(* --------------------------------------------------------------------------
   3.2 基础表示下的矩阵元一般形式
   --------------------------------------------------------------------------
   在格点计算中, 胶子算符在基础表示下构造 (计算更方便):

   M^{mu,lambda; nu,rho}(z) = \sum_{\vec{x}} Tr[
     F^{mu,lambda}(\vec{x} + z\hat{z})
     U(\vec{x} + z\hat{z}, \vec{x})
     F^{nu,rho}(\vec{x})
     U(\vec{x}, \vec{x} + z\hat{z})
   ]

   其中:
   - F^{mu,nu} 由Clover项构造
   - U 是基础表示中的规范链接
   - 迹 (trace) 取遍颜色指标

   利用恒等式 2 Tr[T^a U T^b U^dag] = U^{ab} 将伴随表示转换为
   基础表示。

   矩阵元分解为不变振幅 (Balitsky-Morris-Radyushkin):
   M_{mu,alpha; lambda,beta}(z,p) = ...
   -------------------------------------------------------------------------- *)

(* --------------------------------------------------------------------------
   3.3 不变振幅分解
   --------------------------------------------------------------------------
   利用Lorentz不变性, 矩阵元 M_{mu,alpha; lambda,beta}(z,p) 可分解为
   六个不变振幅:

   M_{pp}:   张量结构 (g_{mu,lambda} p_alpha p_beta - ...)
   M_{zz}:   张量结构 (g_{mu,lambda} z_alpha z_beta - ...)
   M_{zp}:   张量结构 (g_{mu,lambda} z_alpha p_beta - ...)
   M_{pz}:   张量结构 (g_{mu,lambda} p_alpha z_beta - ...)
   M_{ppzz}: 张量结构 (p_mu z_alpha - p_alpha z_mu)(p_lambda z_beta - p_beta z_lambda)
   M_{gg}:   张量结构 (g_{mu,lambda} g_{alpha,beta} - g_{mu,beta} g_{alpha,lambda})

   每个振幅 M_i 都是不变间隔 z^2 和 Ioffe时间 nu = -(p·z) 的函数。

   对格点运动学:
   p^mu = (E, 0, 0, P_z),  z^mu = (0, 0, 0, z_3)
   有: nu = z_3 P_z,  z^2 = -z_3^2

   非极化胶子PDF由 M_{pp} 振幅决定:
   在光锥极限下: g^{alpha,beta} M_{+alpha; beta+}(z^-, p) = -2 p^2_+ M_{pp}(nu, 0)
   -M_{pp}(nu, 0) = 1/2 \int_{-1}^{1} dx e^{-ix nu} x g(x)
   -------------------------------------------------------------------------- *)

InvariantAmplitudeDecomposition[verbose_:False] := Module[{amps, relations},
  amps = {
    "Mpp" -> "(g_{\\mu\\lambda}p_\\alpha p_\\beta - g_{\\mu\\beta}p_\\alpha p_\\lambda - g_{\\alpha\\lambda}p_\\mu p_\\beta + g_{\\alpha\\beta}p_\\mu p_\\lambda) M_{pp}",
    "Mzz" -> "(g_{\\mu\\lambda}z_\\alpha z_\\beta - g_{\\mu\\beta}z_\\alpha z_\\lambda - g_{\\alpha\\lambda}z_\\mu z_\\beta + g_{\\alpha\\beta}z_\\mu z_\\lambda) M_{zz}",
    "Mzp" -> "(g_{\\mu\\lambda}z_\\alpha p_\\beta - g_{\\mu\\beta}z_\\alpha p_\\lambda - g_{\\alpha\\lambda}z_\\mu p_\\beta + g_{\\alpha\\beta}z_\\mu p_\\lambda) M_{zp}",
    "Mpz" -> "(g_{\\mu\\lambda}p_\\alpha z_\\beta - g_{\\mu\\beta}p_\\alpha z_\\lambda - g_{\\alpha\\lambda}p_\\mu z_\\beta + g_{\\alpha\\beta}p_\\mu z_\\lambda) M_{pz}",
    "Mppzz" -> "(p_\\mu z_\\alpha - p_\\alpha z_\\mu)(p_\\lambda z_\\beta - p_\\beta z_\\lambda) M_{ppzz}",
    "Mgg" -> "(g_{\\mu\\lambda}g_{\\alpha\\beta} - g_{\\mu\\beta}g_{\\alpha\\lambda}) M_{gg}"
  };

  relations = {
    "Mpp_to_PDF" ->
      "-M_{pp}(\\nu, 0) = \\frac{1}{2} \\int_{-1}^{1} dx e^{-i x \\nu} x g(x)",
    "IoffeTime" -> "\\nu = z_3 P_z",
    "MppExtraction" ->
      "M_{ti;it} + M_{ij;ji} = 2 p_0^2 M_{pp}  (关键组合!)",
    "Symmetry" ->
      "M_{pp}, M_{zz}, M_{gg}, M_{ppzz}, M_{pz}-M_{zp} 是 \\nu 的偶函数"
  };

  If[verbose,
    Print["=== 不变振幅分解 (Balitsky-Morris-Radyushkin) ==="];
    Print["六个不变振幅: M_{pp}, M_{zz}, M_{zp}, M_{pz}, M_{ppzz}, M_{gg}"];
    Print["非极化胶子PDF <-> M_{pp}:"];
    Print["  -M_{pp}(\\nu) = 1/2 \\int dx e^{-ix\\nu} x g(x)"];
    Print["关键组合: M_{ti;it} + M_{ij;ji} = 2 p_0^2 M_{pp}"];
  ];

  {amps, relations}
];

(* --------------------------------------------------------------------------
   3.4 非极化胶子流算符
   --------------------------------------------------------------------------
   在基础表示下, 使用 M_{pp} 提取关系定义的胶子流:

   乘法可重整化的胶子算符组合 (Zhang et al. 2019):
   O(z) = M^{tx;tx}(z) + M^{ty;ty}(z) - 2 M^{xy;xy}(z)

   该组合的关键性质: M^{ti;it} 和 M^{ij;ji} 具有相同的单圈UV反常量纲,
   使得整个组合在一圈水平上 (至少) 是乘法可重整化的。

   Euclidean空间中的非极化胶子流:
   O^{(0)}_g(z, t_i) = [-F^{0i}(z) W F^{0i}(0) W^dag
                       + 2 F^{12}(z) W F^{12}(0) W^dag]_{Eucl}

   时间分量前的 "-" 号来自 Minkowski<->Euclidean 转换:
   F^{0i}_{(M)} F^{0i}_{(M)} = -F^{0i}_{(E)} F^{0i}_{(E)}

   辅助算符 (用于交叉检验):
   O^{(1)}_g(z, t_i) = F^{0i}(z) W F^{0i}(0) W^dag
   -------------------------------------------------------------------------- *)

GluonCurrentUnpolarized[verbose_:False] := Module[{},
  If[verbose,
    Print["=== 非极化胶子流算符 ==="];
    Print["乘法可重整化组合:"];
    Print["  O(z) = M_{tx;tx}(z) + M_{ty;ty}(z) - 2 M_{xy;xy}(z)"];
    Print[""];
    Print["Euclidean空间形式:"];
    Print["  O^{(0)}_g = -F^{0i} W F^{0i} W^dag + 2 F^{12} W F^{12} W^dag"];
    Print["  O^{(1)}_g =  F^{0i} W F^{0i} W^dag  (辅助)"];
    Print[""];
    Print["Minkowski/Euclidean关键转换:"];
    Print["  F^{0i}_{(M)}F^{0i}_{(M)} = -F^{0i}_{(E)}F^{0i}_{(E)}"];
    Print["  F^{ij}_{(M)}F^{ij}_{(M)} =  F^{ij}_{(E)}F^{ij}_{(E)}"];
  ];

  {
   "MultiplicativelyRenormalizable" ->
     "O(z) = M_{tx;tx} + M_{ty;ty} - 2 M_{xy;xy}",
   "UnpolarizedCurrent_Eucl" ->
     "O^{(0)}_g = -F^{0i}WF^{0i}W^\\dagger + 2F^{12}WF^{12}W^\\dagger",
   "AuxiliaryCurrent_Eucl" ->
     "O^{(1)}_g = F^{0i}WF^{0i}W^\\dagger",
   "KeyRelation" ->
     "M_{ti;it} + M_{ij;ji} = 2 p_0^2 M_{pp}"
  }
];

(* --------------------------------------------------------------------------
   3.5 极化 (螺旋度) 胶子流算符
   --------------------------------------------------------------------------
   极化胶子PDF由对偶场强张量的关联定义:

   Euclidean空间中的螺旋度胶子流:
   O^{(0)}_g(z) = -{ [F^{01} W F^{23} - F^{02} W F^{13} + 2 F^{12} W F^{03}]
                   - (z <-> -z) }_{Eucl}

   其中对偶关系已被展开:
   \tilde{F}^{mu,nu} = 1/2 \epsilon^{mu,nu,rho,sigma} F_{rho,sigma}

   辅助螺旋度算符:
   O^{(1)}_g(z) = {F^{01} W F^{23} + F^{02} W F^{13}} - (z <-> -z)

   螺旋度PDF的Ioffe时间分布:
   -i \tilde{I}_p(\nu) = \tilde{M}^{(+)}_{(ps)}(\nu) - \nu \tilde{M}_{(pp)}(\nu)
   \Delta G = \int_0^\infty d\nu \tilde{I}_p(\nu)
   -------------------------------------------------------------------------- *)

GluonCurrentHelicity[verbose_:False] := Module[{},
  If[verbose,
    Print["=== 极化(螺旋度)胶子流算符 ==="];
    Print["Euclidean螺旋度胶子流:"];
    Print["  O^{(0)}_g = -{[F^{01}WF^{23}-F^{02}WF^{13}+2F^{12}WF^{03}] - (z<->-z)}"];
    Print["螺旋度PDF Ioffe时间分布:"];
    Print["  -i\\tilde{I}_p(\\nu) = \\tilde{M}^{(+)}_{(ps)}(\\nu) - \\nu\\tilde{M}_{(pp)}(\\nu)"];
  ];

  {
   "HelicityCurrent" ->
     "O^{(0)}_g = -{[F^{01}WF^{23}-F^{02}WF^{13}+2F^{12}WF^{03}] - (z<->-z)}",
   "IoffeTimeDistribution" ->
     "-i\\tilde{I}_p(\\nu) = \\tilde{M}^{(+)}_{(ps)} - \\nu\\tilde{M}_{(pp)}",
   "HelicityPDFMoment" ->
     "\\Delta G = \\int_0^\\infty d\\nu \\tilde{I}_p(\\nu)"
  }
];

(* --------------------------------------------------------------------------
   3.6 准胶子算符的乘法可重整性 (Li-Ma-Qiu定理)
   --------------------------------------------------------------------------
   Li, Ma, Qiu (2019) 严格证明:

   在规范不变的UV正规化方案 (如维数正规化 DR) 下,
   所有36个纯准胶子算符的UV发散都局域在时空中, 且可以不互相混合地
   被乘性重整化。

   O^{mu,nu,rho,sigma}_{g,R}(xi) = e^{-C_g|xi_z|} Z^{-1}_{wg}
     Z^{-s/2}_{vg1} Z^{-(2-s)/2}_{vg2} O^{mu,nu,rho,sigma}_{bg}(xi)

   其中 s 是从 {mu,nu,rho,sigma} 中选择z分量的数目。

   分类:
   - s=0 (无z分量): F^{ti}F^{ti}, F^{ij}F^{ij} - 仅需 C_g, Z_{wg}
   - s=1 (一个z分量): F^{ti}F^{zi} - 额外需 Z^{-1/2}_{vg1}
   - s=2 (两个z分量): F^{zi}F^{zi} - 需 Z^{-1}_{vg2}

   关键结论: 准胶子算符在UV重整化阶段不与夸克准PDF混合。
   但在微扰匹配到光锥PDF的阶段, 通过硬匹配核发生 2x2 混合。
   -------------------------------------------------------------------------- *)

MultiplicativeRenormalizability[verbose_:False] := Module[{},
  If[verbose,
    Print["=== 乘法可重整性 (Li-Ma-Qiu 2019) ==="];
    Print["36个准胶子算符均乘法可重整化"];
    Print["重整化常数仅依赖于 z分量数 s"];
    Print["s=0: C_g + Z_{wg}"];
    Print["s=1: C_g + Z_{wg} + Z^{-1/2}_{vg1}"];
    Print["s=2: C_g + Z_{wg} + Z^{-1}_{vg2}"];
  ];

  {
   "Theorem" -> "所有36个纯准胶子算符乘法可重整化 (Li-Ma-Qiu 2019)",
   "RenormFactors" -> {
     "C_g" -> "线性发散 (Wilson线自能)",
     "Z_{wg}" -> "规范链端点对数发散",
     "Z_{vg1}" -> "顶点重整化 (一个z分量)",
     "Z_{vg2}" -> "顶点重整化 (两个z分量)"
   },
   "Classification" -> {
     "s=0" -> "F^{ti}F^{ti}, F^{ij}F^{ij} - 无 Z_{vg}",
     "s=1" -> "F^{ti}F^{zi} - Z^{-1/2}_{vg1}",
     "s=2" -> "F^{zi}F^{zi} - Z^{-1}_{vg2}"
   },
   "FlavorSingletMixing" ->
     "UV重整化阶段无混合; 微扰匹配阶段发生 2x2 混合"
  }
];

(* --------------------------------------------------------------------------
   3.7 准TMD-PDF算符 (区别于准PDF算符)
   --------------------------------------------------------------------------
   准TMD-PDF算符需要额外的横向分离 b_perp, 使用 staple型 Wilson线:

   \tilde{h}_{chi,Gamma}(b_perp, z, L, P_z; 1/a) =
     <chi(P_z)| \bar{psi}(0_perp,0) Gamma W_{sqsubset}(b_perp, L, z)
               psi(b_perp, z) |chi(P_z)>

   对胶子TMD-PDF, 算符涉及场强张量而非夸克场:
   M^{mu,lambda; nu,rho}_{TMD}(z, b_perp) = <P|
     F^{mu,lambda}(b_perp, z) U_{staple}(...) F^{nu,rho}(0, 0) |P>

   其中 U_{staple} 是连接两个空间分离胶子场的staple型伴随表示Wilson线。
   -------------------------------------------------------------------------- *)

(* --------------------------------------------------------------------------
   主推导函数
   -------------------------------------------------------------------------- *)
GluonQuasiTMDDerivation[] := Module[{},
  Print["=== Part 3: 胶子准TMD-PDF算符 ==="];

  Print["\n3.1 光锥胶子PDF定义 (伴随表示):"];
  Print["  x g(x,\\mu) = \\int dz^-/(2\\pi P^+) e^{-ixP^+z^-}"];
  Print["    <P| F^{+\\mu}_a(z^-) U_{ab}(z^-,0) F^{b\\mu}_+(0) |P>"];

  Print["\n3.2 基础表示格点算符:"];
  Print["  M^{\\mu\\lambda;\\nu\\rho}(z) = \\sum_{\\vec{x}} Tr["];
  Print["    F^{\\mu\\lambda}(\\vec{x}+z\\hat{z}) U F^{\\nu\\rho}(\\vec{x}) U^\\dagger]"];

  Print["\n3.3 不变振幅分解 (Balitsky-Morris-Radyushkin):"];
  Print["  6个振幅: M_{pp}, M_{zz}, M_{zp}, M_{pz}, M_{ppzz}, M_{gg}"];
  Print["  非极化胶子PDF 由 M_{pp} 决定"];
  Print["  M_{ti;it} + M_{ij;ji} = 2 p_0^2 M_{pp}"];

  Print["\n3.4 非极化胶子流 (乘法可重整化):"];
  Print["  O(z) = M_{tx;tx} + M_{ty;ty} - 2 M_{xy;xy}"];

  Print["\n3.5 Euclidean非极化胶子流:"];
  Print["  O^{(0)}_g = -F^{0i}WF^{0i}W^dag + 2F^{12}WF^{12}W^dag"];

  Print["\n3.6 极化胶子流:"];
  Print["  O^{(0)}_g = -{[F^{01}WF^{23}-F^{02}WF^{13}+2F^{12}WF^{03}] - (z<->-z)}"];

  Print["\n3.7 乘法可重整性 (Li-Ma-Qiu 2019):"];
  Print["  36个准胶子算符均乘法可重整化"];
  Print["  分类: s=0,1,2 (z分量数) - 不同重整化常数"];

  Print["\n3.8 参考源:"];
  Print["  - Balitsky, Morris, Radyushkin, PLB 808, 135621 (2020)"];
  Print["  - Zhang, Ji et al., PRL 122, 142001 (2019)"];
  Print["  - Li, Ma, Qiu, PRL 122, 062002 (2019)"];
  Print["  - 文档/gluon_pdf_derivation.tex, Sec.2-3, 6"];
  Print["  - 补充/格点QCD中的胶子算符.tex"];
];

End[];
EndPackage[];
GluonQuasiTMDDerivation[];
