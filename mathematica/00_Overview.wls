#!/usr/bin/env wolframscript
(* ::Title:: *)
(*\[ | \[ | *)

(* ::Section:: *)
(*\[ | \[ | *)

(* ::Text:: *)
(*\[ | *)

(* ::Text:: *)
(*\[ | *)

(* ::Text:: *)
(*\[ | *)

(* ::Text:: *)
(*\[ | *)

(* ::Text:: *)
(*\[ | *)

(* ::Text:: *)
(*\[ | *)

(* ::Text:: *)
(*\[ | *)

(* ::Text:: *)
(*\[ | *)

(* ::Text:: *)
(*\[ | *)

Print["========================================================================"];
Print["  \[ | "];
Print["  \[ | "];
Print["  \[ | "];
Print["========================================================================"];
Print["  Author: Based on CLQCD/LPC collaboration documented derivations"];
Print["  Date: 2026-07-21"];
Print["  Repository: /root/lattice-pdf/"];
Print["========================================================================"];
Print[""];

(* ::Section:: *)
(*\[ | \[ | *)

Print["--- Physical Motivation ---"];
Print["The proton spin crisis motivates the Jaffe-Manohar sum rule:"];
Print["  J = (1/2)\[CapitalDelta]\[CapitalSigma] + L_q + L_G + \[CapitalDelta]G"];
Print["where \[CapitalDelta]G is the gluon spin contribution determined by"];
Print["the gluon helicity distribution \[CapitalDelta]g(x)."];
Print[""];
Print["The unpolarized gluon TMD PDF g(x, k_perp) generalizes the collinear PDF"];
Print["to include transverse momentum dependence, giving a 3D tomographic"];
Print["picture of the nucleon."];
Print[""];

(* ::Section:: *)
(*\[ | \[ | *)

Print["--- LaMET (Large Momentum Effective Theory) Framework ---"];
Print["Key insight (Ji, 2013): In a nucleon boosted to large momentum P_z,"];
Print["equal-time spatial correlations approach light-cone correlations."];
Print[""];
Print["Quasi-TMD-PDF can be computed on the Euclidean lattice:"];
Print["  f_tilde(x, b_perp, P_z, \[Mu]) = \[Integral] dz/(2\[Pi]) e^{-i x P_z z}"];
Print["    \[Times] h_R(b_perp, z, P_z, \[Mu])"];
Print[""];
Print["Factorization theorem connects quasi-TMD to light-cone TMD:"];
Print["  f_tilde / Sqrt[S_I] = H(\[Zeta]_z/\[Mu]^2) exp[(1/2) ln(\[Zeta]_z/\[Zeta]) K] f"];
Print[""];

(* ::Section:: *)
(*\[ | \[ | *)

Print["--- Gradient Flow Renormalization ---"];
Print["Gradient flow (Wilson flow) provides a rigorous renormalization framework:"];
Print["  \[PartialD]_t B_\[Mu] = D_\[Nu] G_\[Nu]\[Mu],   B_\[Mu](0,x) = A_\[Mu](x)"];
Print["  \[PartialD]_t \[Chi] = \[CapitalDelta] \[Chi],    \[Chi](0,x) = \[Psi](x)"];
Print[""];
Print["Key theorem (L\[UDoubleDot]scher-Weisz, 2011): For t > 0, all gauge-invariant"];
Print["correlation functions are UV finite after standard QCD renormalization."];
Print["No additional renormalization of the flowed gauge field is needed."];
Print[""];
Print["For the gluon TMD PDF, gradient flow:"];
Print["  (1) Regulates Wilson line power divergences"];
Print["  (2) Preserves rotational symmetry -> only logarithmic mixing"];
Print["  (3) Provides a theoretically rigorous continuum limit"];
Print[""];

(* ::Section:: *)
(*\[ | \[ | *)

Print["--- Continuum Limit Extrapolation ---"];
Print["Joint extrapolation in lattice spacing a and momentum P_z:"];
Print["  x g(x, P_z, a) = x g_0(x) + a^2 F(x) + a^2 P_z^2 H(x) + D(x)/P_z^2"];
Print[""];
Print["Three ensembles used (CLQCD):"];
Print["  C24P29: a = 0.105 fm, m_\[Pi] = 292 MeV"];
Print["  E32P29: a = 0.0897 fm, m_\[Pi] = 287 MeV"];
Print["  F32P30: a = 0.0775 fm, m_\[Pi] = 300 MeV"];
Print[""];

(* ::Section:: *)
(*\[ | \[ | *)

Print["--- Key References ---"];
Print["[1] Note of gluon PDFs (internal note, docs/Note of gluon PDFs.pdf)"];
Print["[2] Chen et al., Unpolarized gluon PDF...continuum limit (2025), arXiv:2510.26425"];
Print["    (docs/Unpolarized gluon PDF of the nucleon from lattice QCD in the continuum limit.pdf)"];
Print["[3] He et al., Unpolarized TMD parton distributions... (2024), PRD 109, 114513"];
Print["    (docs/Unpolarized Transverse-Momentum-Dependent Parton Distributions...pdf)"];
Print["[4] Zhang et al., Renormalization of TMD parton distribution... (2022), PRD 106, 094510"];
Print["    (docs/Renormalization of transverse-momentum-dependent parton distribution...pdf)"];
Print["[5] L\[UDoubleDot]scher, Properties and uses of the Wilson flow (2010), JHEP 08, 071"];
Print["    (docs/Properties and uses of the Wilson flow in lattice QCD_Luscher_2010.pdf)"];
Print["[6] L\[UDoubleDot]scher & Weisz, Perturbative analysis of the gradient flow (2011) JHEP 02, 051"];
Print["    (docs/Perturbative analysis of the gradient flow...pdf)"];
Print["[7] Monahan & Orginos, Quasi parton distributions and the gradient flow (2017) JHEP 03, 116"];
Print["    (docs/Quasi parton distributions and the gradient flow_Monahan_Orginos_2017.pdf)"];
Print["[8] Balitsky, Morris, Radyushkin, PLB 808, 135621 (2020), arXiv:1910.13963"];
Print["[9] Ji, PRL 110, 262002 (2013), LaMET founding paper"];
Print["[10] Ji et al., Hybrid renormalization scheme... (2021), NPB 964, 115311"];
Print["[11] Fan et al., Gluon Quasi-PDF From Lattice QCD (2018), PRL 121, 242001"];
Print["[12] 格点QCD中的梯度流重整化.tex (补充/)"];
Print["[13] 格点QCD中的TMD_PDF.tex (补充/)"];
Print["[14] 格点QCD中的光锥PDF、准PDF、赝PDF.tex (补充/)"];
Print["[15] gluon_pdf_derivation.tex (文档/)"];
Print["[16] gluon_PDF_continuum.tex (文档/)"];
Print[""];

Print["========================================================================"];
Print["  The following files contain the complete Mathematica derivation:"];
Print["    01_FieldStrength_GluonOperator.wls"];
Print["    02_GradientFlow_Renormalization.wls"];
Print["    03_TMD_WilsonLine_SoftFunction.wls"];
Print["    04_QuasiPDF_Matching_Factorization.wls"];
Print["    05_ContinuumLimit_FullWorkflow.wls"];
Print["========================================================================"];

(* ::Section:: *)
(*\[ | \[ | *)

(* ::Text:: *)
(*\[ | *)
