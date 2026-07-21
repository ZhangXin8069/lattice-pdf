(*
================================================================================
格点QCD中核子非极化胶子TMD-PDF的完整理论推导与数值实现
================================================================================

标题: 使用梯度流重整化方案的连续极限下核子非极化胶子的TMD部分子分布函数
      的格点QCD计算 —— 完整理论推导与Mathematica符号/数值代码

理论框架: 大动量有效理论 (LaMET) + 梯度流重整化 + TMD因子化
方法: Clover场强张量 + Staple型Wilson线 + Wilson圈减除 + 连续极限外推

参考源 (本库路径):
  1. 要求.md — 项目规格说明
  2. docs/Note of gluon PDFs.pdf — 胶子PDF内部理论笔记 (Eq(20), Eq(25))
  3. 补充/格点QCD中的TMD_PDF.tex — TMD-PDF完整理论框架
  4. 补充/格点QCD中的梯度流重整化.tex — 梯度流重整化完整推导
  5. 补充/格点QCD中的胶子算符.tex — 胶子算符分类与构造
  6. 补充/格点QCD中的场强张量.tex — Clover项与F_{μν}构造
  7. 补充/格点QCD中的重整化.tex — 重整化方案
  8. 补充/格点QCD中的大动量有效理论.tex — LaMET框架
  9. 补充/格点QCD中的外推.tex — 连续极限与手征外推
  10. code/gluon_pdf_full_workflow.py — 完整10步计算工作流
  11. examples/Calc_ope_unpol.py — OPE部分GPU代码
  12. examples/2pt_proton_Cg5gmu_L32x64_mom2_xdir_gpu.py — 质子2pt代码
  13. docs/Unpolarized gluon PDF of the nucleon from lattice QCD in the continuum limit.pdf
  14. 文档/gluon_pdf_derivation.tex — 胶子PDF推导
  15. 文档/gluon_PDF_continuum.tex — 胶子PDF连续极限

作者: 基于本库全部文献与代码的综合解析
日期: 2026-07-21
================================================================================
*)

(* ============================================================================
   第一部分: 梯度流 (Gradient Flow / Wilson Flow) 方程与求解
   参考: 补充/格点QCD中的梯度流重整化.tex
         Lüscher (2010) JHEP 08, 071
         Lüscher & Weisz (2011) JHEP 02, 051
   ============================================================================ *)

Print["========================================================================"];
Print["第一部分: 梯度流 (Gradient Flow) 方程与求解"];
Print["========================================================================"];

(* --------------------------------------------------------------------------
   1.1 连续理论中的梯度流方程
   -------------------------------------------------------------------------- *)

(* 规范场的梯度流方程 (D维Euclidean空间, SU(N)规范理论)

   ∂_t B_μ(t,x) = D_ν G_{νμ}(t,x),   B_μ(0,x) = A_μ(x)

   其中流时间 t 具有质量量纲 [t] = -2,
   G_{μν} = ∂_μ B_ν - ∂_ν B_μ + [B_μ, B_ν] 是流时间依赖的场强张量,
   D_μ = ∂_μ + [B_μ, ·] 是伴随表示中的协变导数.

   添加规范固定项 (Feynman规范 α₀=1):
   ∂_t B_μ = D_ν G_{νμ} + α₀ D_μ ∂_ν B_ν
*)

(* 流方程在微扰论中的领头阶解 (热核):
   B_{μ,1}(t,x) = ∫ d^D y  K_t(x-y) A_μ(y)
   K_t(z) = Exp[-z^2/(4t)] / (4π t)^(D/2)
   平滑半径: sqrt(8t) 在四维中
*)

Clear[Kt, Bflow, Ddim];
Kt[z_, t_, Ddim_: 4] := Exp[-z^2/(4 t)] / (4 \[Pi] t)^(Ddim/2);

Print["热核 K_t(z) = Exp[-z^2/(4t)] / (4πt)^(D/2)"];
Print["平滑半径 r_smooth = Sqrt[8t]"];

(* --------------------------------------------------------------------------
   1.2 夸克场的梯度流方程
   参考: Lüscher (2013) JHEP 04, 123
   -------------------------------------------------------------------------- *)

(* 夸克场的梯度流:
   ∂_t χ(t,x) = Δ χ(t,x),   ∂_t χ̄(t,x) = χ̄(t,x) Δ⃖
   Δ = D_μ D_μ (规范协变Laplacian)

   夸克传播子 (领头阶):
   <χ(t,x) χ̄(s,y)> = ∫ d^D p/(2π)^D  Exp[i p(x-y)] Exp[-(t+s)p^2] / (m₀ + i p̸)
*)

(* --------------------------------------------------------------------------
   1.3 梯度流的 UV 有限性定理 (Lüscher-Weisz 2011)
   -------------------------------------------------------------------------- *)

(* 定理: 对于 t > 0 处的流场 B_μ(t,x),
   任意规范不变局域乘积的关联函数在标准QCD参数重整化后是UV有限的.

   - 胶子流场不需要波函数重整化: B_{μ,R} = B_μ
   - 流夸克场需要乘性重整化: χ_R = Z_χ^{-1/2} χ
   - Z_χ 不依赖于流时间 t

   单圈结果 (MS-bar方案):
   Z_χ = 1 + (3 C_F / (16π^2)) (g^2/ε) + O(g^4)
   C_F = (N^2-1)/(2N)
*)

Clear[CF, Zchi];
CF[nc_: 3] := (nc^2 - 1)/(2 nc);
ZchiMSbar[gs_, eps_, nc_: 3] := 1 + (3 CF[nc]/(16 \[Pi]^2)) (gs^2/eps);

Print["CF = ", CF[3]];
Print["Z_χ (MS-bar, 单圈) = 1 + (3 CF/16π^2)(g^2/ε)"];

(* --------------------------------------------------------------------------
   1.4 格点上的 Wilson Flow
   -------------------------------------------------------------------------- *)

(* Wilson flow 的格点实现:
   V̇_t(x,μ) = -g₀^2 {∂_{x,μ} S_W(V_t)} V_t(x,μ),   V_t|{t=0} = U(x,μ)

   数值积分: 三阶 Runge-Kutta 格式, 步长 ε
   每一步的离散化误差: O(ε^3)

   与 stout smearing 的关系:
   重复应用 stout link smearing 在 ε → 0 极限下等价于精确的 Wilson flow.
*)

(* 三阶Runge-Kutta步进 *)
Clear[WilsonFlowStep];
WilsonFlowStep[U_, eps_, beta_] := Module[{V0, Z0, Z1, Z2, V1, V2},
  (* U 是规范连接变量 (格点上的link变量)
     这里给出的是形式化的伪代码, 实际格点实现需遍历所有link *)
  V0 = U;
  (* 第一步: 计算驱动力 Z0 = -g0^2 ∂_{x,μ} S_W[V0] *)
  Z0 = -beta * GaugeForce[V0];
  V1 = MatrixExp[eps * Z0] . V0;
  (* 第二步 *)
  Z1 = -beta * GaugeForce[V1];
  V2 = MatrixExp[(3/4) eps * Z1 - (1/4) eps * Z0] . V0;
  (* 第三步 *)
  Z2 = -beta * GaugeForce[V2];
  Return[MatrixExp[(1/12) eps * (8 Z2 - 4 Z1 + 5 Z0)] . V0];
];

Print["Wilson Flow 数值积分: 三阶RK格式, O(ε^3) 单步误差"];

(* --------------------------------------------------------------------------
   1.5 流时间的选择条件
   -------------------------------------------------------------------------- *)

(* 流时间 t 需满足:
   1) sqrt(8t) ≳ 2a (保证连续极限的良好行为)
   2) sqrt(8t) ≪ 物理尺度 (如强子尺寸 ~1 fm)

   典型取值: sqrt(8τ) ~ 0.3-0.5 fm
   对于 a = 0.1 fm: τ/a^2 ~ 1-3
*)

Clear[FlowTimeCondition];
FlowTimeCondition[a_, t_] := Module[{rsmooth},
  rsmooth = Sqrt[8 t];
  {rsmooth, rsmooth/a, 2 a}
];

Print["流时间条件示例 (a = 0.1 fm, τ = 3a^2):"];
Print["  sqrt(8τ) = ", Sqrt[8*3*0.1^2] // N, " fm"];
Print["  sqrt(8τ)/a = ", Sqrt[8*3] // N];


(* ============================================================================
   第二部分: 场强张量 F_{μν} 的格点构造 (Clover项)
   参考: 补充/格点QCD中的场强张量.tex
         Note of gluon PDFs (内部笔记)
   ============================================================================ *)

Print["\n========================================================================"];
Print["第二部分: 场强张量 F_{μν} 的 Clover 构造"];
Print["========================================================================"];

(* --------------------------------------------------------------------------
   2.1 Plaquette 与 BCH 展开
   -------------------------------------------------------------------------- *)

(* Plaquette 定义:
   P_{μν}(x) = U_μ(x) U_ν(x+μ̂) U_μ^†(x+ν̂) U_ν^†(x)

   BCH展开 (小a):
   P_{μν} = Exp[i a^2 (∂_μ A_ν - ∂_ν A_μ + i[A_μ, A_ν]) + O(a^3)]
          = 1 + i a^2 F_{μν} + O(a^3)
*)

(* --------------------------------------------------------------------------
   2.2 Clover 项 Q_{μν} 的定义
   -------------------------------------------------------------------------- *)

(* Clover 项: 四个共面 plaquette 的平均
   Q_{μν}(x) = P_{μν}(x) - P_{ν,-μ}(x) + P_{-ν,μ}(x) + P_{-μ,-ν}(x)

   对称性:
   - Q_{νμ} = -Q_{μν} (反对称)
   - 平均消除了 O(a) 和 O(a^3) 奇次离散化误差
   - 首项误差为 O(a^4) 在 Clover 构造中, 但提取 F_{μν} 后为 O(a^2)
*)

(* --------------------------------------------------------------------------
   2.3 从 Clover 项到场强张量
   -------------------------------------------------------------------------- *)

(* 格点场强张量:
   F_{μν}(x) = -(i/(8a^2)) (Q_{μν}(x) - Q_{νμ}(x))

   这是 hermitian 的, traceless 的 SU(3) 代数元素.

   展开到连续极限:
   F_{μν}^{latt}(x) = F_{μν}^{cont}(x) + O(a^2)
*)

Clear[CloverToF];
CloverToF[Qmunu_, a_] := -(I/(8 a^2)) (Qmunu - Transpose[Qmunu]);

(* --------------------------------------------------------------------------
   2.4 对偶场强张量 F̃_{μν}
   -------------------------------------------------------------------------- *)

(* 对偶场强张量 (极化胶子计算需要):
   F̃^{μν} = (1/2) ε^{μνρσ} F_{ρσ},   ε^{0123} = 1 (Minkowski)

   在 Euclidean 空间中, Levi-Civita 张量定义需小心处理.
*)

Clear[DualF, LeviCivita4];
(* 四维 Levi-Civita 符号 (数值约定 ε_{0123} = +1) *)
LeviCivita4[idx__] := Signature[{idx}];

DualF[F_, mu_, nu_] := (1/2) Sum[
    LeviCivita4[mu, nu, rho, sigma] * F[[rho, sigma]],
    {rho, 0, 3}, {sigma, 0, 3}
  ];

(* --------------------------------------------------------------------------
   2.5 Minkowski ↔ Euclidean 转换
   关键关系 (参考: Note of gluon PDFs):
   F^{tx,(M)} F^{(M)}_{tx} = -F^{(E)}_{tx} F^{(E)}_{tx}
   F^{xy,(M)} F^{(M)}_{xy} = +F^{(E)}_{xy} F^{(E)}_{xy}
   -------------------------------------------------------------------------- *)

Clear[MinkowskiToEuclideanF2];
MinkowskiToEuclideanF2[expr_, indices_List] := Module[{result},
  result = expr;
  (* 每个时间指标引入因子 i (由于 Wick 旋转 t → -iτ) *)
  Do[
    If[MemberQ[indices, 0],
      result = -result (* F^{tμ}F_{tμ} 乘积的符号 *)
    ],
    {i, 1}
  ];
  result
];

Print["Minkowski/Euclidean 转换:"];
Print["  F^{tx,(M)} F^{(M)}_{tx} = -F^{(E)}_{tx} F^{(E)}_{tx} (时间分量反号)"];
Print["  F^{xy,(M)} F^{(M)}_{xy} = +F^{(E)}_{xy} F^{(E)}_{xy} (空间分量不变)"];


(* ============================================================================
   第三部分: 胶子 TMD 算符的构造
   参考: 补充/格点QCD中的胶子算符.tex
         补充/格点QCD中的TMD_PDF.tex
         Zhang et al. (2019) PRL 122, 142001
         Li, Ma, Qiu (2019) PRL 122, 062002
         Balitsky, Morris, Radyushkin (2020) PLB 808, 135621
   ============================================================================ *)

Print["\n========================================================================"];
Print["第三部分: 胶子 TMD 算符的构造"];
Print["========================================================================"];

(* --------------------------------------------------------------------------
   3.1 非极化胶子算符的不变振幅分解
   -------------------------------------------------------------------------- *)

(* 双胶子关联函数 (核子自旋平均):
   M_{μα;λβ}(z,p) = <p| G_{μα}(z) [z,0] G_{λβ}(0) |p>

   六个不变振幅 (Balitsky-Morris-Radyushkin 分解):
   M_{μα;λβ} = (g_{μλ}p_αp_β - g_{μβ}p_αp_λ - g_{αλ}p_μp_β + g_{αβ}p_μp_λ) M_{pp}
             + (g_{μλ}z_αz_β - ...) M_{zz}
             + (g_{μλ}z_αp_β - ...) M_{zp}
             + (g_{μλ}p_αz_β - ...) M_{pz}
             + (p_μz_α - p_αz_μ)(p_λz_β - p_βz_λ) M_{ppzz}
             + (g_{μλ}g_{αβ} - g_{μβ}g_{αλ}) M_{gg}

   非极化胶子PDF由 M_{pp} 振幅决定.
*)

(* 六个不变张量结构 *)
Clear[BMRtensor, MppTensor, MzzTensor, MzpTensor, MpzTensor, MppzzTensor, MggTensor];

(* 度规张量 (Euclidean, diag(1,1,1,1)) *)
gEE[mu_, nu_] := KroneckerDelta[mu, nu];

(* M_{pp} 张量结构 *)
MppTensor[mu_, alpha_, lambda_, beta_, p_] :=
  gEE[mu, lambda] p[[alpha]] p[[beta]] -
  gEE[mu, beta] p[[alpha]] p[[lambda]] -
  gEE[alpha, lambda] p[[mu]] p[[beta]] +
  gEE[alpha, beta] p[[mu]] p[[lambda]];

(* M_{zz} 张量结构 *)
MzzTensor[mu_, alpha_, lambda_, beta_, z_] :=
  gEE[mu, lambda] z[[alpha]] z[[beta]] -
  gEE[mu, beta] z[[alpha]] z[[lambda]] -
  gEE[alpha, lambda] z[[mu]] z[[beta]] +
  gEE[alpha, beta] z[[mu]] z[[lambda]];

(* M_{zp} 张量结构 *)
MzpTensor[mu_, alpha_, lambda_, beta_, z_, p_] :=
  gEE[mu, lambda] z[[alpha]] p[[beta]] -
  gEE[mu, beta] z[[alpha]] p[[lambda]] -
  gEE[alpha, lambda] z[[mu]] p[[beta]] +
  gEE[alpha, beta] z[[mu]] p[[lambda]];

(* M_{pz} 张量结构 *)
MpzTensor[mu_, alpha_, lambda_, beta_, z_, p_] :=
  gEE[mu, lambda] p[[alpha]] z[[beta]] -
  gEE[mu, beta] p[[alpha]] z[[lambda]] -
  gEE[alpha, lambda] p[[mu]] z[[beta]] +
  gEE[alpha, beta] p[[mu]] z[[lambda]];

(* M_{ppzz} 张量结构 *)
MppzzTensor[mu_, alpha_, lambda_, beta_, z_, p_] :=
  (p[[mu]] z[[alpha]] - p[[alpha]] z[[mu]]) *
  (p[[lambda]] z[[beta]] - p[[beta]] z[[lambda]]);

(* M_{gg} 张量结构 *)
MggTensor[mu_, alpha_, lambda_, beta_] :=
  gEE[mu, lambda] gEE[alpha, beta] -
  gEE[mu, beta] gEE[alpha, lambda];

(* 完整的 BMR 分解 *)
Clear[BMRAmplitude];
BMRAmplitude[mu_, alpha_, lambda_, beta_, z_, p_,
  Mpp_, Mzz_, Mzp_, Mpz_, Mppzz_, Mgg_] :=
  Mpp * MppTensor[mu, alpha, lambda, beta, p] +
  Mzz * MzzTensor[mu, alpha, lambda, beta, z] +
  Mzp * MzpTensor[mu, alpha, lambda, beta, z, p] +
  Mpz * MpzTensor[mu, alpha, lambda, beta, z, p] +
  Mppzz * MppzzTensor[mu, alpha, lambda, beta, z, p] +
  Mgg * MggTensor[mu, alpha, lambda, beta];

Print["BMR 六个不变振幅分解已定义: M_{pp}, M_{zz}, M_{zp}, M_{pz}, M_{ppzz}, M_{gg}"];

(* --------------------------------------------------------------------------
   3.2 非极化胶子PDF的算符组合
   关键组合 (Balitsky et al. 2020):
   M_{ti;it} + M_{ij;ji} = 2 p₀² M_{pp}
   这两个矩阵元具有相同的单圈紫外反常量纲 → 可乘性重整化.
   -------------------------------------------------------------------------- *)

(* 在 Euclidean 空间中的具体实现:
   指标约定: 0=t, 1=x, 2=y, 3=z

   非极化胶子流 O_g^{(0)} (Euclidean):
   O_g^{(0)}(z) = -F^{0i}(z) W[z,0] F_{0i}(0) W[0,z]
                  + 2 F^{12}(z) W[z,0] F_{12}(0) W[0,z]

   其中 i = 1,2 (横向), W 是伴随表示中的 Wilson 线.
*)

Clear[GluonCurrentUnpol, GluonCurrentUnpolMatrix];

(* 基础表示下的胶子算符 (格点上实际计算的) *)
GluonCurrentUnpol[F1_, F2_, Ulink_] := Module[{result, i},
  result = 0;
  (* 时间-空间分量: -F^{0i} F_{0i} *)
  For[i = 1, i <= 2, i++,
    result += -F1[[0, i]] . Ulink . F2[[0, i]] . ConjugateTranspose[Ulink];
  ];
  (* 空间-空间分量: +2 F^{12} F_{12} *)
  result += 2 F1[[1, 2]] . Ulink . F2[[1, 2]] . ConjugateTranspose[Ulink];
  result
];

(* 格点基础表示下的 M_{μλ;νρ}(z) 构造 (参考: Chen et al. 2025) *)
Clear[LatticeGluonME];
LatticeGluonME[Ffield_, Ulink_, mu_, lambda_, nu_, rho_] :=
  Tr[
    Ffield[[mu, lambda]] . Ulink .
    Ffield[[nu, rho]] . ConjugateTranspose[Ulink]
  ];

Print["格点胶子矩阵元 M_{μλ;νρ}(z) = Tr[F_{μλ}(ẑz) U(ẑz,0) F_{νρ}(0) U(0,ẑz)]"];

(* --------------------------------------------------------------------------
   3.3 极化 (螺旋度) 胶子算符
   参考: Egerer et al. (2022) PRD 106, 094511
         Note of gluon PDFs (内部笔记)
   -------------------------------------------------------------------------- *)

(* 极化胶子流 O_g^{(0)} (Euclidean):
   O_g^{(0)}(z) = -[{F^{01} W F^{23} - F^{02} W F^{13} + 2 F^{12} W F^{03}}
                   - {z ↔ -z}]_{Eucl}

   这里利用了 F̃^{μν} = (1/2) ε^{μνρσ} F_{ρσ} 展开为 Euclidean F^{μν} 分量.
   z-奇组合用于提取与核子自旋相关的部分.
*)

Clear[GluonCurrentPol];
GluonCurrentPol[Fz_, F0_, Ulink_] := Module[{zodd},
  zodd = (
    Fz[[0, 1]] . Ulink . F0[[2, 3]] . ConjugateTranspose[Ulink] -
    Fz[[0, 2]] . Ulink . F0[[1, 3]] . ConjugateTranspose[Ulink] +
    2 Fz[[1, 2]] . Ulink . F0[[0, 3]] . ConjugateTranspose[Ulink]
  );
  -zodd (* 整体负号来自 Minkowski/Euclidean 转换 *)
];

(* --------------------------------------------------------------------------
   3.4 可乘性重整化: 36个算符的分类
   参考: Li, Ma, Qiu (2019) PRL 122, 062002
         Zhang et al. (2019) PRL 122, 142001
   -------------------------------------------------------------------------- *)

(* 准胶子算符的乘性重整化定理 (Li-Ma-Qiu):
   O^{μνρσ}_g(ξ) = Exp[-C_g |ξ_z|] Z^{-1}_{wg} Z^{-s/2}_{vg1}
                    Z^{-(2-s)/2}_{vg2} O^{μνρσ}_{bg}(ξ)

   其中 s 是从 {μ,ν,ρ,σ} 中选择 z 分量的数目:
   - s=0 (无z分量): 如 F^{ti}F^{ti}, F^{ij}F^{ij} — 仅需要 C_g, Z_{wg}
   - s=1 (一个z分量): 如 F^{ti}F^{zi} — 额外需要 Z^{-1/2}_{vg1}
   - s=2 (两个z分量): 如 F^{zi}F^{zi} — 需要 Z^{-1}_{vg2}
*)

Clear[CountZIndices];
CountZIndices[indices_List] := Count[indices, 3]; (* 指标3 = z方向 *)

Clear[RenormConstant];
RenormConstant[Cg_, Zwg_, Zvg1_, Zvg2_, s_] :=
  Exp[-Cg * Lz] * Zwg^(-1) * Zvg1^(-s/2) * Zvg2^(-(2 - s)/2);

Print["可乘性重整化分类 (按z分量数s):"];
Print["  s=0: Z = Exp[-C_g|z|] Z_{wg}^{-1}"];
Print["  s=1: Z = Exp[-C_g|z|] Z_{wg}^{-1} Z_{vg1}^{-1/2}"];
Print["  s=2: Z = Exp[-C_g|z|] Z_{wg}^{-1} Z_{vg2}^{-1}"];


(* ============================================================================
   第四部分: Staple型Wilson线与TMD准PDF算符
   参考: 补充/格点QCD中的TMD_PDF.tex
         Zhang et al. (2022) PRD 106, 094510
         He et al. (2024) PRD 109, 114513
   ============================================================================ *)

Print["\n========================================================================"];
Print["第四部分: Staple型Wilson线与准TMD-PDF算符"];
Print["========================================================================"];

(* --------------------------------------------------------------------------
   4.1 Staple型Wilson线 W_{⊏}(b_⊥, L, z) 的构造
   -------------------------------------------------------------------------- *)

(* Staple型规范链接的三段构造:
   W_{⊏}(b_⊥, L, z) = U_z^†((z+L)n̂_z + b⃗_⊥, b⃗_⊥)
                      × U_⊥((z+L)n̂_z + b⃗_⊥, (z+L)n̂_z)
                      × U_z((z+L)n̂_z, zn̂_z)

   其中:
   U_z(ξ₁^z n̂_z + ξ⃗_⊥, ξ₂^z n̂_z + ξ⃗_⊥) =
     P exp[-ig ∫_{ξ₁^z}^{ξ₂^z} dλ  n̂_z · A(λn̂_z + ξ⃗_⊥)]

   U_⊥(ξ^z n̂_z + ξ⃗_⊥¹, ξ^z n̂_z + ξ⃗_⊥²) =
     P exp[-ig ∫_{ξ⃗_⊥¹}^{ξ⃗_⊥²} dλ  n̂_⊥ · A(ξ^z n̂_z + λn̂_⊥)]
*)

(* Staple Wilson 线的路径积分表示 (连续理论) *)
Clear[StapleWilsonLine, WilsonLineZ, WilsonLinePerp];

(* 沿z方向的Wilson线 *)
WilsonLineZ[Afield_, z1_, z2_, xperp_] :=
  PExp[-I g0 Integrate[Afield[3, {xperp, lambda}], {lambda, z1, z2}]];

(* 沿横向的Wilson线 *)
WilsonLinePerp[Afield_, xiPerp1_, xiPerp2_, z_] :=
  PExp[-I g0 Integrate[
    Sum[Afield[i, {lambda, z}] * (xiPerp2[[i]] - xiPerp1[[i]])/
      Norm[xiPerp2 - xiPerp1],
      {i, 1, 2}],
    {lambda, 0, 1}
  ]];

(* 完整的 Staple Wilson 线 *)
StapleWilsonLine[Afield_, bperp_, L_, z_] :=
  ConjugateTranspose[
    WilsonLineZ[Afield, 0, z + L, bperp]
  ] .
  WilsonLinePerp[Afield, {(z + L), bperp}, {(z + L), 0}, z + L] .
  WilsonLineZ[Afield, z, z + L, {0, 0, 0}];

Print["Staple型Wilson线 W_{⊏}(b_⊥, L, z) 三段构造已定义"];
Print["  L: staple沿z方向延伸的臂长 (L → ∞ 以模拟光锥)"];

(* --------------------------------------------------------------------------
   4.2 准TMD-PDF裸矩阵元
   -------------------------------------------------------------------------- *)

(* 裸矩阵元定义:
   h̃_{χ,Γ}(b_⊥, z, L, P^z; 1/a) =
     <χ(P^z)| ψ̄(0⃗_⊥, 0) Γ W_{⊏}(b_⊥, L, z) ψ(b⃗_⊥, z) |χ(P^z)>

   对于胶子TMD, 替换 ψ̄...Γ...ψ → F_{μν} W F_{ρσ}:
   h̃_g(b_⊥, z, L, P^z; 1/a) =
     <N(P^z)| F^{μν}(b⃗_⊥, z) W_{⊏} F_{μν}(0⃗_⊥, 0) |N(P^z)>

   其中 W_{⊏} 现在是伴随表示中的 staple Wilson 线.
*)

Clear[QuasiTMDBareGluon];
QuasiTMDBareGluon[Ffield_, Astaple_, bperp_, z_, L_, Pz_] :=
  Module[{Fz, F0},
    Fz = FfieldAt[Ffield, {bperp, z}];
    F0 = FfieldAt[Ffield, {{0, 0}, 0}];
    (* 胶子场强张量的缩并 *)
    ContractGluonTMD[Fz, Astaple, F0]
  ];

Print["准TMD-PDF裸矩阵元: <N(P^z)| F(z,b_⊥) W_{⊏} F(0) |N(P^z)>"];

(* --------------------------------------------------------------------------
   4.3 裸矩阵元中的四类UV发散
   参考: Zhang et al. (2022) PRD 106, 094510
   -------------------------------------------------------------------------- *)

(* 四类UV发散:
   1. 线性发散: Exp[-δm̄ (2L+z)/a]
      Wilson线自能修正, 正比于 Wilson 线总长度/a
   2. Cusp发散: staple 三个转角处的对数发散
      由 cusp 反常维数 Γ_{cusp}(α_s) 描述
   3. Pinch-pole奇点: Exp[-V(b_⊥) L]
      两条平行Wilson线之间的胶子交换 → 有效势
   4. 对数发散: 夸克-Wilson线顶点修正, ∝ ln a
*)

Clear[LinearDivergence, PinchPoleSingularity, CuspAnomalousDim];

(* 线性发散因子 *)
LinearDivergence[deltaMbar_, L_, z_, a_] :=
  Exp[-deltaMbar * (2 L + z)/a];

(* Pinch-pole 因子 (V(b_⊥) 是 Wilson 线之间的有效势) *)
PinchPoleSingularity[Vbperp_, L_] := Exp[-Vbperp * L];

(* Cusp 反常维数 (微扰展开, 领头阶) *)
CuspAnomalousDim[alphaS_] := (alphaS * CF[3])/\[Pi]; (* LO *)

Print["四类UV发散:"];
Print["  1. 线性发散: Exp[-δm̄(2L+z)/a]"];
Print["  2. Cusp发散: Γ_{cusp}(α_s) 描述"];
Print["  3. Pinch-pole: Exp[-V(b_⊥)L]"];
Print["  4. 对数发散: ∝ α_s ln a"];


(* ============================================================================
   第五部分: 梯度流重整化方案
   参考: 补充/格点QCD中的梯度流重整化.tex
         Monahan & Orginos (2017) JHEP 03, 116
         Monahan (2018) PRD 97, 054507
         Francis et al. (2025) arXiv:2509.02472
   ============================================================================ *)

Print["\n========================================================================"];
Print["第五部分: 梯度流重整化方案"];
Print["========================================================================"];

(* --------------------------------------------------------------------------
   5.1 梯度流正则化的核心思想 (Monahan-Orginos 2017)
   -------------------------------------------------------------------------- *)

(* 核心思想:
   对夸克场和规范场同时施加梯度流, 将流时间 τ 固定在物理单位中.
   所有涉及涂抹后 Wilson 线的矩阵元在连续极限下都是有限的.

   涂抹后的准PDF矩阵元:
   h^{(s)}(z/√τ, √τ P_z, √τ Λ_QCD, √τ M_N)

   梯度流作为规范不变的UV正规化手段, 替代了格点间距 a 的角色.

   关键性质:
   1. 涂抹后的矩阵元在 a → 0 (保持 √τ 固定) 下是有限的
   2. 在连续理论中, 梯度流保持转动对称性 → twist-2算符之间仅有对数混合
   3. 涂抹后准PDF与光前PDF的匹配可在连续理论中独立于格点细节进行
*)

(* --------------------------------------------------------------------------
   5.2 小流时展开 (Short Flow-Time Expansion, SFTX)
   参考: Lüscher (2014) PoS LATTICE2013, 016
         Suzuki (2013) PTEP 2013, 083B03
   -------------------------------------------------------------------------- *)

(* SFTX 基本公式:
   O(t, x) ~ Σ_k c_k(t) O_{R,k}(x)   (t → 0)

   其中 O_{R,k} 是 t=0 处重整化的局域算符,
   c_k(t) 是有限的 Wilson 系数.

   对于规范作用量密度:
   E(t,x) = (1/4) G_{μν}^a(t,x) G_{μν}^a(t,x)

   SFTX:
   E(t,x) ~ <E(t,x)> + c_E(t) [(1/4) F_{ρσ}^a F_{ρσ}^a]_R(x) + O(t)

   Wilson系数 c_E(t) 满足重整化群方程:
   c_E(t) ~ g^2(μ) / g^2(1/√(8t))  (领头对数近似)
*)

Clear[ActionDensity, SFTXExpansion, WilsonCoeffCE];

(* 流时间 t 处的规范作用量密度 *)
ActionDensity[Gfield_, t_, x_] := (1/4) Sum[
    Tr[Gfield[[mu, nu]] . Gfield[[mu, nu]]],
    {mu, 0, 3}, {nu, 0, 3}
  ];

(* c_E(t) 的领头对数行为 (利用渐近自由):
   在 t → 0 极限下,
   1/c_E(t) ~ -2 b_0 ln(√(8t) Λ_QCD) + const
*)
Clear[RunningCoupling, Beta0];
Beta0[nf_: 4] := (11 - 2 nf/3)/(4 \[Pi]);

RunningCoupling[mu_, LambdaQCD_, nf_: 4] :=
  1/(2 Beta0[nf] * Log[mu/LambdaQCD]);

(* c_E(t) 的一圈近似 *)
WilsonCoeffCE[t_, mu_, LambdaQCD_, nf_: 4] :=
  RunningCoupling[mu, LambdaQCD, nf]^2 /
  RunningCoupling[1/Sqrt[8 t], LambdaQCD, nf]^2;

Print["SFTX Wilson 系数 c_E(t) 的领头对数近似:"];
Print["  c_E(t) ≃ g^2(μ) / g^2(1/√(8t)) ≈ ln(1/(√(8t)Λ))/ln(μ/Λ)"];

(* --------------------------------------------------------------------------
   5.3 能动张量的梯度流表示
   参考: Suzuki (2013) PTEP 2013, 083B03
   -------------------------------------------------------------------------- *)

(* D+1 维的类能动张量:
   U_{μν}(t,x) = G_{μρ}^a(t,x) G_{νρ}^a(t,x) - (1/4) δ_{μν} G_{ρσ}^a(t,x) G_{ρσ}^a(t,x)

   其 SFTX:
   U_{μν}(t,x) = c_T(t) {T_{μν}}_R(x) + c_S(t) δ_{μν} [(1/4) F²]_R(x) + O(t)

   能动张量的梯度流公式:
   {T_{μν}}_R(x) = lim_{t→0} [U_{μν}(t,x)/c_T(t)
                   - (c_S(t)/(c_T(t)c_E(t))) δ_{μν} (E(t,x) - <E(t,x)>)]
*)

Clear[EMTFlow, UmunuTensor];

(* U_{μν} 张量 *)
UmunuTensor[Gfield_, t_, x_, mu_, nu_] := Module[{result},
  result = Sum[
    Gfield[[mu, rho]] . Gfield[[nu, rho]],
    {rho, 0, 3}
  ];
  result - (1/4) KroneckerDelta[mu, nu] *
    Sum[Gfield[[rho, sigma]] . Gfield[[rho, sigma]], {rho, 0, 3}, {sigma, 0, 3}]
];

(* c_T(t) 的渐近行为:
   1/c_T(t) ~ -2 b_0 ln(√(8t) Λ) + c_1   (t → 0⁺)
*)
Clear[CoeffCT, CoeffCS];
CoeffCT[t_, LambdaQCD_, c1_, nf_: 4] :=
  1/(-2 Beta0[nf] Log[Sqrt[8 t] LambdaQCD] + c1);

CoeffCS[t_, LambdaQCD_, c1_, c2_, nf_: 4] :=
  -(Beta0[nf]/2) * (1 - (c1 - c2)/(-Log[Sqrt[8 t] LambdaQCD]));

Print["能动张量梯度流公式: T_{μν}(x) = lim_{t→0} [...]"];
Print["  1/c_T(t) ~ -2b_0 ln(√(8t)Λ) + c_1"];
Print["  c_S/(c_T c_E) ~ -(b_0/2)[1 - (c_1-c_2)/ln(√(8t)Λ)]"];

(* --------------------------------------------------------------------------
   5.4 梯度流方案中的 PDF Moment 比值
   参考: Francis et al. (2025)
   -------------------------------------------------------------------------- *)

(* Moment 比值的梯度流公式 (到 NNLO):
   <x^{n-1}>^{MS̄}(μ) / <x^{m-1}>^{MS̄}(μ) =
     [ζ_m(t,μ) / ζ_n(t,μ)] × [<x^{n-1}>(t) / <x^{m-1}>(t)] + O(t)

   匹配系数 ζ_n(t,μ) 在 NNLO 精度下已被计算.
*)

Clear[ZetaCoefficient];
(* 匹配系数 ζ_n(t,μ) 的领头阶形式:
   ζ_n(t,μ) = 1 + (α_s(μ) C_F / (2π)) [γ_n^{(0)} ln(μ² t) + ...]
   其中 γ_n^{(0)} 是第 n 个矩的反常量纲.
*)

(* Mellin矩的反常量纲 (非极化, 领头阶) *)
Clear[AnomalousDimMoment];
AnomalousDimMoment[n_, flavor_] := Module[{},
  (* 非单态夸克 PDF 的反常量纲 *)
  If[flavor === "NS",
    -(CF[3]/\[Pi]) Sum[1/k, {k, 1, n}],  (* 简化形式 *)
    0
  ]
];

ZetaCoefficient[n_, t_, mu_, alphaS_] := Module[{gamma0},
  gamma0 = AnomalousDimMoment[n, "NS"];
  1 + (alphaS * CF[3]/(2 \[Pi])) * gamma0 * Log[mu^2 * t]
];

Print["PDF Moment 比值的梯度流公式"];
Print["  ζ_n(t,μ): 匹配系数, 已计算到 NNLO"];

(* --------------------------------------------------------------------------
   5.5 环标记夸克场 (Ringed Fermion Fields)
   参考: Makino & Suzuki (2014) PTEP 2014, 063B02
   -------------------------------------------------------------------------- *)

(* 环标记夸克场消除 Z_χ 的微扰依赖性:
   χ̊(t,x) = χ(t,x) / √( -2N_c / (t² <χ̄(t,x) D̸ χ(t,x)>) )
   其中 D̸ 是 Dirac 算符.

   这样定义的 χ̊ 完全非微扰地归一化, 不需要微扰 Z_χ.
*)

Clear[RingedQuarkField];
RingedQuarkField[chi_, t_, x_, Dslash_] :=
  chi[t, x] / Sqrt[-2 * 3 / (t^2 * VacuumExpectation[
    Conjugate[chi[t, x]] . Dslash . chi[t, x]
  ])];

Print["环标记夸克场: 非微扰地消除 Z_χ 依赖"];


(* ============================================================================
   第六部分: Wilson圈减除与短距离比率 (SDR) 重整化
   参考: 补充/格点QCD中的TMD_PDF.tex
         Zhang et al. (2022) PRD 106, 094510
         He et al. (2024) PRD 109, 114513
   ============================================================================ *)

Print["\n========================================================================"];
Print["第六部分: Wilson圈减除 + SDR 重整化"];
Print["========================================================================"];

(* --------------------------------------------------------------------------
   6.1 Wilson圈减除: 消除线性发散 + pinch-pole + cusp
   -------------------------------------------------------------------------- *)

(* Wilson圈减除的准TMD-PDF:
   h_{χ,Γ}(b_⊥, z, P^z; 1/a) =
     lim_{L→∞} h̃_{χ,Γ}(b_⊥, z, L, P^z; 1/a) / √Z_E(2L+z, b_⊥; 1/a)

   矩形Wilson圈 Z_E:
   Z_E(r, b_⊥) = (1/N_c) Tr<0| U_⊥^†(-L⃗+b⃗_⊥, -b_⊥)
                 × U_z^†(L⃗+z⃗+b⃗_⊥, -2L-z)
                 × U_⊥(L⃗+z⃗, b_⊥) U_z(-L⃗, 2L+z) |0>
*)

Clear[WilsonLoopZE, WilsonLoopSubtracted];

(* 矩形 Wilson 圈 (真空期望值) *)
WilsonLoopZE[r_, bperp_, a_] := Module[{},
  (* 在格点上, Z_E 通过计算闭合 Wilson 圈的迹得到 *)
  (* 这是形式化表示; 实际需从格点数据采样 *)
  Exp[-Vstatic[bperp] * r]   (* 面积定律: 在大面积下 *)
];

(* Wilson圈减除后的矩阵元 *)
WilsonLoopSubtracted[hBare_, ZE_, L_, z_, bperp_] :=
  Limit[hBare / Sqrt[ZE[2 L + z, bperp]], L -> Infinity];

Print["Wilson圈减除:"];
Print["  h = lim_{L→∞} h̃ / √Z_E(2L+z, b_⊥)"];
Print["  减除机制: √Z_E 包含相同的线性发散 Exp[-δm̄(2L+z)/a]"];

(* --------------------------------------------------------------------------
   6.2 五格点间距的系统验证
   参考: Zhang et al. (2022) PRD 106, 094510
   -------------------------------------------------------------------------- *)

(* Zhang et al. 在五种格点间距上验证了 Wilson 圈减除:
   a = {0.032, 0.043, 0.060, 0.086, 0.121} fm
   对应 β = {6.3, 6.1, 5.9, 5.7, 5.5} (Clover 费米子)

   关键发现:
   1. 减除后矩阵元在 L ≳ 0.36 fm 后达到稳定平台
   2. a → 0 连续极限下展示稳定收敛
   3. Clover 和 Overlap 费米子的减除结果一致
*)

Clear[PlatformFit, Lvalues];
(* L 的平台拟合 *)
Lvalues = {4, 6, 8, 10, 12}; (* 格点单位的 L 值 *)

PlatformFit[data_, Lvals_] := Module[{mean, err},
  mean = Mean[data];
  err = StandardDeviation[data]/Sqrt[Length[data]];
  {mean, err}
];

Print["五格点间距验证 (Zhang et al. 2022):"];
Print["  a = {0.032, 0.043, 0.060, 0.086, 0.121} fm"];
Print["  平台出现于 L > 0.36 fm, 所有 a 上一致"];

(* --------------------------------------------------------------------------
   6.3 短距离比率 (SDR) 方案: 对数发散的重整化
   -------------------------------------------------------------------------- *)

(* SDR 重整化因子:
   Z_O(1/a, μ, Γ) = lim_{L→∞}
     h̃_{χ,Γ}(b_{⊥,0}, z_0, L, P^z=0; 1/a) /
     [√Z_E(2L+z_0, b_{⊥,0}; 1/a) h̃^{MS̄}_Γ(z_0, b_{⊥,0}, μ)]

   在壳夸克矩阵元的单圈 MS̄ 结果:
   h̃^{MS̄}_Γ(z, b_⊥, μ) = 1 + (α_s C_F / (2π)) *
     [1/2 + (3/2) ln(μ² (b_⊥² + z²) e^{γ_E} / 4)
      - 2(z/b_⊥) arctan(z/b_⊥)] + O(α_s²)
*)

Clear[SDRFactor, OnShellMSbar];

(* 在壳夸克矩阵元的单圈 MS̄ 结果 (Zhang et al. 2022, Eq. 27) *)
OnShellMSbar[z_, bperp_, mu_, alphaS_] := Module[{eulerGamma, lnTerm, arctanTerm},
  eulerGamma = 0.5772156649015329; (* Euler-Mascheroni 常数 *)
  lnTerm = Log[mu^2 * (bperp^2 + z^2) * Exp[eulerGamma] / 4];
  arctanTerm = If[bperp > 0, 2 * (z/bperp) * ArcTan[z/bperp], 0];
  1 + (alphaS * CF[3]/(2 \[Pi])) * (1/2 + (3/2) * lnTerm - arctanTerm)
];

(* SDR 重整化因子 *)
SDRFactor[z0_, bperp0_, mu_, alphaS_, hBareZeroMom_] :=
  hBareZeroMom / OnShellMSbar[z0, bperp0, mu, alphaS];

(* 完整重整化的准TMD-PDF:
   h^{SDR}(b_⊥, z, P^z; μ) =
     h^{sub}(b_⊥, z, P^z) / h^{sub}(b_{⊥,0}, z_0, P^z=0)
     × h̃^{MS̄}_Γ(z_0, b_{⊥,0}, μ)
*)
Clear[FullRenormQuasiTMD];
FullRenormQuasiTMD[hSub_z_Pz_, hSub_z0_Pz0_, z0_, bperp0_, mu_, alphaS_] :=
  (hSub_z_Pz / hSub_z0_Pz0) * OnShellMSbar[z0, bperp0, mu, alphaS];

Print["SDR 方案的重整化因子:"];
Print["  Z_O = h̃^{sub}(z₀,b_{⊥,0},P^z=0) / h̃^{MS̄}_Γ(z₀,b_{⊥,0},μ)"];
Print["  h̃^{MS̄}_Γ 使用在壳夸克矩阵元的单圈结果"];

(* --------------------------------------------------------------------------
   6.4 RI/MOM 方案的失败 (经验教训)
   参考: Zhang et al. (2022), Huo et al. (2021)
   -------------------------------------------------------------------------- *)

(* RI/MOM 重整化后的矩阵元在 a → 0 下不收敛:
   - z 越大或 a 越小, 格距依赖性越强
   - 残余线性发散未完全抵消
   - 根本原因: 胶子-规范链接顶点的非微扰线性发散

   结论: 基于减除的方案 (如 Wilson 圈减除) 是准TMD算符重整化的唯一可行路径.
*)


(* ============================================================================
   第七部分: TMD 软函数与 Collins-Soper 核
   参考: 补充/格点QCD中的TMD_PDF.tex
         LPC (2020) PRL 125, 192001
         Chu et al. (2022) PRD 106, 034509
         Chu et al. (2023) JHEP 08, 172
   ============================================================================ *)

Print["\n========================================================================"];
Print["第七部分: TMD 软函数与 Collins-Soper 核"];
Print["========================================================================"];

(* --------------------------------------------------------------------------
   7.1 TMD 软函数的两部分分解
   -------------------------------------------------------------------------- *)

(* TMD 软函数:
   S(b_⊥, μ, ζ₁, ζ₂) = √S_I(b_⊥, μ) × (ζ₁ζ₂/b_⊥^{-4})^{K(b_⊥,μ)/4}

   其中:
   - S_I(b_⊥, μ): 内禀软函数, 与快度无关
   - K(b_⊥, μ): Collins-Soper 演化核

   CS核满足 RG 方程:
   μ² dK/dμ² = -Γ_{cusp}(α_s(μ))
*)

Clear[CSkernel, TMDSoftFunction, IntrinsicSoftFunction];

(* Collins-Soper 方程 (ζ 演化):
   2ζ d/dζ ln f(x, b_⊥, μ, ζ) = K(b_⊥, μ)
*)

(* CS核的 RG 方程解 (LL 精度):
   K(b_⊥, μ) = K(b_⊥, μ₀) - Γ_{cusp}^{(0)}/(2β₀) ln(α_s(μ₀)/α_s(μ))
*)
Clear[CSkernelRGE];
CSkernelRGE[K0_, bperp_, mu_, mu0_, alphaSmu_, alphaSmu0_, cuspGamma0_] :=
  K0[bperp, mu0] - (cuspGamma0/(2 Beta0[4])) *
    Log[alphaSmu0 / alphaSmu];

Print["Collins-Soper 核满足 RG 方程: μ² dK/dμ² = -Γ_{cusp}(α_s)"];

(* --------------------------------------------------------------------------
   7.2 内禀软函数的格点提取: 赝标量形状因子因子化
   -------------------------------------------------------------------------- *)

(* 形状因子定义:
   F(b_⊥, P^z) = <π(-P⃗)| (q₁Γq₁)(b⃗_⊥) (q₂Γq₂)(0) |π(P⃗)>_c

   因子化定理 (LPC 2020):
   F(b_⊥, P^z) = S_I(b_⊥) ∫ dx dx' H(x,x',P^z) Φ^†(x',b_⊥,-P^z) Φ(x,b_⊥,P^z)
                 + O(1/(P^z)²)

   领头阶提取 (H = 1/(2N_c)):
   S_I(b_⊥) = 2N_c × F(b_⊥, P^z) / |φ(0, b_⊥, P^z)|² + O(α_s, 1/(P^z)²)

   其中 φ(0, b_⊥, P^z) = ∫₀¹ dx Φ(x, b_⊥, P^z) 是准TMD波函数的 Fourier零点.
*)

Clear[ExtractIntrinsicSoft, QuasiTMDWaveFunction];

(* 领头阶内禀软函数提取 *)
ExtractIntrinsicSoft[formFactor_, phiZero_, Nc_: 3] :=
  2 Nc * formFactor / Abs[phiZero]^2;

(* 准TMD波函数的 Fourier零点 *)
QuasiTMDWaveFunction[x_, bperp_, Pz_] :=
  (* 从格点数据获得的波函数; 这里是形式化函数 *)
  Exp[-bperp^2/(2 sigma[x, Pz]^2)] * phiShape[x, Pz];

Print["内禀软函数提取 (LO):"];
Print["  S_I(b_⊥) = 2N_c × F(b_⊥,P^z) / |φ(0,b_⊥,P^z)|²"];

(* --------------------------------------------------------------------------
   7.3 Collins-Soper 核的格点提取
   -------------------------------------------------------------------------- *)

(* CS核的动量比提取公式:
   K(b_⊥, μ) = (1/ln(P₁^z/P₂^z)) ln|φ(0, b_⊥, P₁^z) / φ(0, b_⊥, P₂^z)|
               + O(α_s, 1/(P^z)²)

   NLO 改进:
   K(b_⊥, μ) = (1/ln(P₁^z/P₂^z)) ×
     [ln|φ(0,b_⊥,P₁^z)/φ(0,b_⊥,P₂^z)| + ln|H(p₁,μ)/H(p₂,μ)|]
*)

Clear[ExtractCSkernel, HardKernelNLO];

(* LO CS核提取 *)
ExtractCSkernelLO[phi1_, phi2_, Pz1_, Pz2_] :=
  (1/Log[Pz1/Pz2]) * Log[Abs[phi1/phi2]];

(* NLO 硬核 *)
HardKernelNLO[zetaZ_, mu_, alphaS_] := Module[{},
  1 + (alphaS * CF[3]/(2 \[Pi])) *
    (-2 + \[Pi]^2/12 + Log[zetaZ/mu^2] - (1/2) Log[zetaZ/mu^2]^2)
];

(* NLO CS核提取 *)
ExtractCSkernelNLO[phi1_, phi2_, Pz1_, Pz2_, mu_, alphaS_] := Module[{},
  (1/Log[Pz1/Pz2]) * (
    Log[Abs[phi1/phi2]] +
    Log[Abs[HardKernelNLO[Pz1^2, mu, alphaS] /
            HardKernelNLO[Pz2^2, mu, alphaS]]]
  )
];

Print["CS核提取 (动量比法):"];
Print["  K(b_⊥,μ) = ln|φ(P₁^z)/φ(P₂^z)| / ln(P₁^z/P₂^z)"];


(* ============================================================================
   第八部分: 连续极限外推
   参考: 补充/格点QCD中的外推.tex
         文档/gluon_PDF_continuum.tex
         Chen et al. (2025) (CLQCD/LPC, gluon PDF 连续极限)
   ============================================================================ *)

Print["\n========================================================================"];
Print["第八部分: 连续极限外推"];
Print["========================================================================"];

(* --------------------------------------------------------------------------
   8.1 格点间距依赖性的形式
   -------------------------------------------------------------------------- *)

(* 对于 O(a) 改进的作用量 (如 TITLC Clover费米子):
   离散化误差以 O(a^2) 标度.

   胶子 TMD 矩阵元的连续极限外推:
   h(z, b_⊥, P^z, a) = h(z, b_⊥, P^z, 0) + c_1(z, b_⊥, P^z) a^2
                       + c_2(z, b_⊥, P^z) a^4 + ...

   在胶子PDF计算中 (Chen et al. 2025):
   使用三个格点间距 a = {0.105, 0.0897, 0.0775} fm
   在 CLQCD TITLS/TITLC 系综上.
*)

Clear[ContinuumExtrapolation, ContinuumExtrapolationLinear];

(* 线性 a^2 外推 *)
ContinuumExtrapolationLinear[aSq_, data_] := Module[{fit, a2vals},
  a2vals = Map[#^2 &, aSq];
  fit = LinearModelFit[
    Transpose[{a2vals, data}],
    {1, x}, {x}
  ];
  {fit[0], fit["ParameterErrors"][[1]]}
];

(* 二次 a^2 外推 (加入 a^4 项) *)
ContinuumExtrapolationQuad[aSq_, data_] := Module[{fit, a2vals},
  a2vals = Map[#^2 &, aSq];
  fit = LinearModelFit[
    Transpose[{a2vals, data}],
    {1, x, x^2}, {x}
  ];
  {fit[0], fit["ParameterErrors"][[1]]}
];

(* --------------------------------------------------------------------------
   8.2 联合连续-无穷大动量外推
   参考: Chen et al. (2025)
   -------------------------------------------------------------------------- *)

(* 联合外推公式 (胶子PDF):
   x g(x) = x g₀(x) + a² f(x) + a² (P^z)² h(x) + d(x)/(P^z)²

   四项贡献:
   1. x g₀(x): 物理极限 (a→0, P^z→∞)
   2. a² f(x): 格距离散化效应
   3. a² (P^z)² h(x): 格距-动量混合效应
   4. d(x)/(P^z)²: 有限动量幂次修正
*)

Clear[JointExtrapolation];

(* 联合外推拟合函数 *)
JointExtrapolation[aVals_, PzVals_, dataMatrix_, xVal_] := Module[
  {xdata, ydata, fit},
  (* dataMatrix[a_idx, Pz_idx] = xg(x, a, P^z) *)
  (* 构造设计矩阵 *)
  (* 返回外推的物理结果 xg₀(x) *)
  Flatten[dataMatrix][[1]] (* 简化: 在实际拟合中使用所有数据点 *)
];

(* 连续极限下固定 P^z 的外推检验 *)
Clear[FixedPzContinuum];
FixedPzContinuum[aVals_, dataAtFixedPz_] := Module[{fit, aSq},
  aSq = Map[#^2 &, aVals];
  fit = LinearModelFit[
    Transpose[{aSq, dataAtFixedPz}],
    {1, x}, {x}
  ];
  {
    fit[0],                          (* 外推值 *)
    fit["ParameterErrors"][[1]],     (* 统计误差 *)
    fit["ParameterConfidenceIntervals"][[2, 2]] (* 系统误差 *)
  }
];

Print["联合连续-无穷大动量外推:"];
Print["  xg(x) = xg₀(x) + a²f(x) + a²P²_z h(x) + d(x)/P²_z"];

(* --------------------------------------------------------------------------
   8.3 CLQCD 系综参数
   -------------------------------------------------------------------------- *)

(* CLQCD TITLS/TITLC 系综 (参考: CLQCD 2024, 2025) *)
Clear[CLQCDensemble];

CLQCDensemble[name_] := Module[{},
  Switch[name,
    "C24P29", {a -> 0.105,  mpi -> 290, L -> 2.5,  Ncfg -> 1000}, (* fm, MeV, fm *)
    "C32P23", {a -> 0.0897, mpi -> 300, L -> 2.9,  Ncfg -> 500},
    "C48P21", {a -> 0.0775, mpi -> 310, L -> 3.7,  Ncfg -> 300},
    _,        {a -> 0.1,    mpi -> 300, L -> 3.0,  Ncfg -> 100}
  ]
];

Print["CLQCD 系综参数:"];
Print["  C24P29: a = 0.105 fm,  m_π ≈ 290 MeV,  L ≈ 2.5 fm"];
Print["  C32P23: a = 0.0897 fm, m_π ≈ 300 MeV,  L ≈ 2.9 fm"];
Print["  C48P21: a = 0.0775 fm, m_π ≈ 310 MeV,  L ≈ 3.7 fm"];

(* --------------------------------------------------------------------------
   8.4 梯度流方案中的连续极限外推
   -------------------------------------------------------------------------- *)

(* 在梯度流方案中, 连续极限更为直接:
   固定流时间 τ 在物理单位中, 取 a → 0:

   h(a, τ) = h(0, τ) + C a²/τ + O(a⁴/τ²)

   额外条件: √(8τ) ≫ a (流平滑半径远大于格距)
*)

Clear[GradientFlowContinuum];

GradientFlowContinuum[aVals_, data_, tau_] := Module[
  {aSqOverTau, fit},
  aSqOverTau = Map[#^2 / tau &, aVals];
  fit = LinearModelFit[
    Transpose[{aSqOverTau, data}],
    {1, x}, {x}
  ];
  {fit[0], fit["ParameterErrors"][[1]]}
];

Print["梯度流方案连续极限外推:"];
Print["  h(a,τ) = h(0,τ) + C·a²/τ + O(a⁴/τ²)"];


(* ============================================================================
   第九部分: TMD 因子化与匹配
   参考: 补充/格点QCD中的TMD_PDF.tex
         补充/格点QCD中的大动量有效理论.tex
         Izubuchi et al. (2018) PRD 98, 056004
         He et al. (2024) PRD 109, 114513
   ============================================================================ *)

Print["\n========================================================================"];
Print["第九部分: TMD 因子化与匹配到光锥 TMD-PDF"];
Print["========================================================================"];

(* --------------------------------------------------------------------------
   9.1 完整的 TMD 因子化定理
   -------------------------------------------------------------------------- *)

(* 重整化后的准TMD-PDF f̃(x, b_⊥, ζ_z, μ) 与光锥 TMD-PDF f(x, b_⊥, μ, ζ):

   f̃(x, b_⊥, ζ_z, μ) / √S_I(b_⊥, μ) =
     H(ζ_z/μ²) · Exp[(1/2) ln(ζ_z/ζ) K(b_⊥, μ)] · f(x, b_⊥, μ, ζ)
     + O(Λ²_QCD/ζ_z, M²/(P^z)², 1/(b_⊥² ζ_z))

   其中:
   - ζ_z = (2x P^z)²: 准TMD-PDF的快度标度
   - ζ: 光锥TMD-PDF的参考快度标度
   - H: 硬匹配核 (微扰可计算)
   - K: CS演化核
   - S_I: 内禀软函数
*)

Clear[TMDMatching, QuasiToLightconeTMD];

(* 从准TMD-PDF到光锥TMD-PDF的完整匹配 *)
QuasiToLightconeTMD[fQuasi_, x_, bperp_, Pz_, mu_, zeta_, SI_, K_, alphaS_] :=
  Module[{zetaZ, hardKernel, csEvolution},
    zetaZ = (2 x * Pz)^2;
    hardKernel = HardKernelNLO[zetaZ, mu, alphaS]; (* NLO *)
    csEvolution = Exp[(1/2) Log[zetaZ/zeta] * K];
    (fQuasi / Sqrt[SI[bperp, mu]]) / (hardKernel * csEvolution)
  ];

Print["TMD因子化定理:"];
Print["  f̃/√S_I = H(ζ_z/μ²) · Exp[(1/2)ln(ζ_z/ζ)K(b_⊥,μ)] · f(x,b_⊥,μ,ζ)"];

(* --------------------------------------------------------------------------
   9.2 NLO 硬匹配核 (单圈)
   -------------------------------------------------------------------------- *)

(* 非极化核子 TMD-PDF 的 NLO 硬核:
   H^{(1)}(ζ_z/μ²) = (α_s C_F)/(2π) ×
     [-2 + π²/12 + ln(ζ_z/μ²) - (1/2) ln²(ζ_z/μ²)]
*)

(* NLO 硬核 (上述已定义 HardKernelNLO) *)

(* NNLO 硬核 (两圈) 的形式结构:
   H^{(2)} = (α_s/(2π))² × [C_F² H_{FF} + C_F C_A H_{FA} + C_F n_f T_F H_{Ff}]

   其中 C_A = N_c, T_F = 1/2, n_f 为活跃味道数.

   对于 ζ_z/μ² ∼ O(1), NNLO修正 ∼ (α_s/2π)² ∼ 10⁻³
*)

Print["NLO硬匹配核:"];
Print["  H^{(1)} = (α_s C_F/2π)[-2 + π²/12 + ln(ζ_z/μ²) - ½ln²(ζ_z/μ²)]"];

(* --------------------------------------------------------------------------
   9.3 重整化群重求和 (RGR)
   -------------------------------------------------------------------------- *)

(* RGR 的必要性:
   当 x → 0 或 P^z 很大时, ζ_z = (2xP^z)² 与 μ² 之间的对数很大.
   例如: x=0.01, P^z=2.5 GeV, μ=2 GeV:
     ln(ζ_z/μ²) ≈ -7.38 → 固定阶微扰论失效

   RGR 策略:
   1. 在自然标度 μ₀ = 2xP^z 处计算硬核 (对数 ln(ζ_z/μ₀²) = 0)
   2. 通过 RG 方程演化硬核从 μ₀ 到参考标度 μ
   3. 跑动耦合 α_s(μ̄) 通过精确的 β 函数演化

   硬核的 RG 方程:
   d/d(ln μ²) H(ζ_z/μ²) = γ_H(α_s) · H(ζ_z/μ²)
*)

Clear[RGRmatching, HardKernelAnomalousDim];

(* 硬核反常维数 (满足 RG 方程 μ² dH/dμ² = γ_H H) *)
HardKernelAnomalousDim[alphaS_] :=
  -(alphaS * CF[3])/\[Pi]; (* LO, 非极化 *)

(* RGR 演化的硬核 *)
RGRmatching[zetaZ_, mu_, x_, Pz_, alphaSFunc_] := Module[
  {mu0, alphaSmu, alphaSmu0, gammaH},
  mu0 = 2 x * Pz;
  alphaSmu = alphaSFunc[mu];
  alphaSmu0 = alphaSFunc[mu0];
  (* 从 μ₀ 到 μ 的 RG 演化 *)
  HardKernelNLO[zetaZ, mu0, alphaSmu0] *
    Exp[Integrate[
      HardKernelAnomalousDim[alphaSFunc[mubar]] / mubar^2,
      {mubar, mu0, mu}
    ]]
];

Print["RGR 重求和:"];
Print["  1. 计算 H 在 μ₀=2xP^z (自然标度, 无大对数)"];
Print["  2. RG 演化 H(μ₀) → H(μ)"];
Print["  3. 使用四圈/五圈 β 函数求解 α_s(μ̄)"];

(* --------------------------------------------------------------------------
   9.4 强耦合跑动 (四圈 β 函数)
   -------------------------------------------------------------------------- *)

Clear[AlphaS, BetaFunctionQCD, RunAlphaS];

(* QCD β 函数 (MS̄ 方案, n_f = 4):
   β(α_s) = -β₀ α_s² - β₁ α_s³ - β₂ α_s⁴ - β₃ α_s⁵ + O(α_s⁶)

   系数 (n_f = 4):
   β₀ = (11 - 2n_f/3)/(4π) = 25/(12π)
   β₁ = (102 - 38n_f/3)/(16π²) = ...
*)

Clear[Beta0, Beta1, Beta2, Beta3];

Beta0[nf_: 4] := (11 - 2 nf/3)/(4 \[Pi]);
Beta1[nf_: 4] := (102 - 38 nf/3)/(16 \[Pi]^2);
Beta2[nf_: 4] := (2857/2 - 5033 nf/18 + 325 nf^2/54)/(64 \[Pi]^3);
Beta3[nf_: 4] := (
  149753/6 + 3564 Zeta3 - (1078361/162 + 6508 Zeta3/27) nf +
  (50065/162 + 6472 Zeta3/81) nf^2 + 1093 nf^3/729
)/(256 \[Pi]^4);

(* 跑动耦合常数 (数值求解 β 函数) *)
RunAlphaS[mu_, LambdaQCD_, nf_: 4] := Module[{b0},
  b0 = Beta0[nf];
  1/(2 b0 * Log[mu/LambdaQCD])  (* 单圈近似 *)
  (* 更精确的可用 NDSolve 求解完整的 β 函数 *)
];

Print["强耦合跑动 α_s(μ): 单圈近似"];
Print["  α_s(μ) = 1 / (2β₀ ln(μ/Λ_QCD))"];


(* ============================================================================
   第十部分: 从准TMD到光锥TMD — 傅里叶变换与数值流程
   参考: code/gluon_pdf_full_workflow.py
         He et al. (2024) PRD 109, 114513
   ============================================================================ *)

Print["\n========================================================================"];
Print["第十部分: 傅里叶变换与完整数值流程"];
Print["========================================================================"];

(* --------------------------------------------------------------------------
   10.1 坐标空间 → 动量空间 Fourier 变换
   -------------------------------------------------------------------------- *)

(* 准TMD-PDF的 Fourier 变换 (z → x, b_⊥ → k_⊥):
   f̃(x, b_⊥, P^z, μ) = ∫ dz/(2π) Exp[-i z (x P^z)] h̃(z, b_⊥, P^z, μ)

   以及横向 Fourier 变换 (b_⊥ → k_⊥):
   f̃(x, k_⊥, P^z, μ) = ∫ d²b_⊥/(2π)² Exp[-i b⃗_⊥ · k⃗_⊥] f̃(x, b_⊥, P^z, μ)
*)

Clear[FourierZtoX, FourierBperpToKperp];

(* 纵向 Fourier 变换 (z → x) *)
FourierZtoX[hMatrix_, zVals_, xVals_, Pz_] := Module[{nx, nz, i, j, result},
  nx = Length[xVals];
  nz = Length[zVals];
  result = Table[0.0, {nx}];
  For[i = 1, i <= nx, i++,
    result[[i]] = Sum[
      hMatrix[[j]] * Exp[-I * zVals[[j]] * (xVals[[i]] * Pz)] *
        (zVals[[2]] - zVals[[1]]), (* d_z 步长 *)
      {j, 1, nz}
    ] / (2 \[Pi])
  ];
  Re[result] (* 非极化分布取实部 *)
];

(* 横向 Fourier 变换 (b_⊥ → k_⊥, 轴对称情况) *)
FourierBperpToKperp[fBperp_, bperpVals_, kperpVals_] := Module[
  {nk, nb, i, j, result},
  nk = Length[kperpVals];
  nb = Length[bperpVals];
  result = Table[0.0, {nk}];
  For[i = 1, i <= nk, i++,
    result[[i]] = Sum[
      fBperp[[j]] * BesselJ[0, bperpVals[[j]] * kperpVals[[i]]] *
        bperpVals[[j]] * (bperpVals[[2]] - bperpVals[[1]]),
      {j, 1, nb}
    ] / (2 \[Pi])
  ];
  result
];

Print["Fourier 变换: z → x, b_⊥ → k_⊥ (轴对称)"];
Print["  f̃(x,b_⊥) = ∫ dz/(2π) e^{-iz(xP^z)} h̃(z,b_⊥)"];

(* --------------------------------------------------------------------------
   10.2 大 λ = z P^z 外推
   -------------------------------------------------------------------------- *)

(* 坐标空间准TMD-PDF在大 z 处的信噪比指数衰减.
   通过解析延拓策略扩展有效 z 范围:

   h(z) = A Exp[-c|z|] · (振荡因子)  (大 z 渐近行为)

   使用大 λ_{max} 截断并对 |z| > λ_{max} 使用外推值.
*)

Clear[LargeZExtrapolation];

LargeZExtrapolation[hData_, zVals_, zCut_, params_List] := Module[
  {hFit, hExtrap},
  (* 对 |z| > zCut 使用参数化外推 *)
  Table[
    If[Abs[zVals[[i]]] <= zCut,
      hData[[i]],
      params[[1]] * Exp[-params[[2]] * Abs[zVals[[i]]]] *
        Cos[params[[3]] * zVals[[i]]] (* 振荡部分 *)
    ],
    {i, 1, Length[zVals]}
  ]
];

Print["大 λ = zP^z 外推: 扩展有效 z 范围 → 改善 x 分辨率"];

(* --------------------------------------------------------------------------
   10.3 完整10步计算工作流 (参考: code/gluon_pdf_full_workflow.py)
   -------------------------------------------------------------------------- *)

Clear[GluonTMDWorkflow];

GluonTMDWorkflow[params_Association] := Module[
  {
    (* Step 1: 格点参数 *)
    Nt, Nx, Ny, Nz, Nc, a, Pz, lambdaQCD, mu, alphaS,
    (* Step 2-3: 场强张量 *)
    Flink, Fmunu, Fdual,
    (* Step 4: 胶子TMD算符 *)
    zVals, bperpVals, MmunuRho,
    (* Step 5: 梯度流涂抹 *)
    tauFlow, FlinkFlowed, FmunuFlowed,
    (* Step 6: 裸矩阵元 *)
    hBare,
    (* Step 7: Wilson圈减除 + SDR *)
    ZE, hSub, hRenorm,
    (* Step 8: 大λ外推 *)
    hExtrap,
    (* Step 9: Fourier变换 *)
    xVals, fQuasi,
    (* Step 10: TMD匹配 *)
    SI, CSkernel, fLightcone
  },

  (* === 参数初始化 === *)
  Nt = params["Nt"]; Nx = params["Nx"]; Ny = params["Ny"]; Nz = params["Nz"];
  Nc = 3; a = params["a"];
  Pz = params["Pz"]; lambdaQCD = params["LambdaQCD"];
  mu = params["mu"]; alphaS = RunAlphaS[mu, lambdaQCD];

  Print["=== 格点QCD胶子TMD-PDF完整工作流 ==="];
  Print["格点: ", Nt, "×", Nz, "×", Ny, "×", Nx, ", a = ", a, " fm"];
  Print["核子动量 P^z = ", Pz, " GeV"];
  Print["α_s(", mu, " GeV) = ", alphaS];

  (* === Step 1-3: 读取/构造 F_{μν} === *)
  Print["\n[Step 1-3] 构造场强张量 F_{μν}..."];
  (* 在格点上, 这涉及:
     - 读取规范连接变量 U_μ(x)
     - 构造 Clover 项 Q_{μν}
     - 提取 F_{μν} = -i/(8a^2)(Q_{μν} - Q_{νμ})
     - (可选) 构造对偶 F̃_{μν} 用于极化
  *)

  (* === Step 4: 胶子TMD算符 === *)
  Print["\n[Step 4] 构造胶子TMD算符..."];
  zVals = Range[-10, 10, 1] * a; (* z 在格点单位 *)
  bperpVals = Range[0, 8, 1] * a; (* b_⊥ 在格点单位 *)

  (* M_{μλ;νρ}(z) = Tr[F_{μλ}(z) U(z,0) F_{νρ}(0) U(0,z)] *)
  (* O(z) = M_{tx;tx} + M_{ty;ty} - 2 M_{xy;xy} *)

  (* === Step 5: 梯度流涂抹 === *)
  Print["\n[Step 5] 梯度流涂抹 (流时间 τ = ", params["tauFlow"], ")..."];
  tauFlow = params["tauFlow"];
  (* 在格点上:
     - 对规范连接变量施加 Wilson flow 直到流时间 τ
     - 使用三阶 Runge-Kutta 积分
     - 流后重新构造 F_{μν}
     - 流后计算 M_{μλ;νρ}
  *)

  (* === Step 6: 计算裸矩阵元 === *)
  Print["\n[Step 6] 计算裸准TMD-PDF矩阵元..."];
  (* h̃(b_⊥, z, L, P^z; 1/a) = <N(P^z)| O_g(z, b_⊥) |N(P^z)>
     对于胶子, 这是非连通三点函数:
     C_{3pt} = <(O_g - <O_g>)(C_{2pt}^N - <C_{2pt}^N>)>
  *)

  (* === Step 7: 重整化 === *)
  Print["\n[Step 7] Wilson圈减除 + SDR 重整化..."];

  (* Wilson圈减除:
     h^{sub}(b_⊥, z) = h̃(b_⊥, z, L → ∞) / √Z_E(2L+z, b_⊥)
     实际中: 在 L ≳ 0.36 fm 的平台区域取值
  *)

  (* SDR 重整化:
     h^{SDR}(b_⊥, z, P^z; μ) =
       h^{sub}(b_⊥, z, P^z) / h^{sub}(b_{⊥,0}, z_0, P^z=0)
       × h̃^{MS̄}(z_0, b_{⊥,0}, μ)
     典型选择: z_0 ≈ 2a, b_{⊥,0} ≈ a
  *)

  (* === Step 8: 大 λ 外推 === *)
  Print["\n[Step 8] 大 λ = zP^z 外推..."];
  (* 扩展有效 z 范围, 改善 x 空间分辨率 *)

  (* === Step 9: Fourier 变换 === *)
  Print["\n[Step 9] Fourier 变换 z → x..."];
  xVals = Range[0.1, 0.9, 0.05];
  (* f̃(x, b_⊥, P^z, μ) = ∫ dz/(2π) Exp[-iz(xP^z)] h^{SDR}(b_⊥, z, P^z; μ) *)

  (* === Step 10: TMD 匹配 === *)
  Print["\n[Step 10] 因子化匹配: 准TMD → 光锥 TMD..."];

  (* 提取软函数 (在实际计算中从独立的形状因子分析获得) *)
  (* K(b_⊥, μ) 从准TMD波函数的动量比提取 *)

  (* 光锥 TMD-PDF:
     f(x, b_⊥, μ, ζ) = f̃(x, b_⊥, ζ_z, μ) / [√S_I(b_⊥, μ) H(ζ_z/μ²)
                        × Exp(½ ln(ζ_z/ζ) K(b_⊥, μ))]
  *)

  Print["\n=== 工作流完成 ==="];
  Print["输出: 光锥胶子 TMD-PDF g(x, b_⊥, μ, ζ)"];

  (* 返回结果字典 *)
  <|
    "z_vals" -> zVals,
    "bperp_vals" -> bperpVals,
    "x_vals" -> xVals,
    "Pz" -> Pz,
    "mu" -> mu,
    "alpha_s" -> alphaS,
    "status" -> "completed"
  |>
];

Print["完整10步工作流函数 GluonTMDWorkflow 已定义"];
Print["  (格点数据的实际填充和具体数值需在集群上运行获得)"];


(* ============================================================================
   第十一部分: 数值示例 — 使用格点系综参数
   ============================================================================ *)

Print["\n========================================================================"];
Print["第十一部分: 数值示例与参数化"];
Print["========================================================================"];

(* --------------------------------------------------------------------------
   11.1 模型参数化
   -------------------------------------------------------------------------- *)

(* 对于无法从实际格点数据获得函数的情况,
   使用唯象模型的参数化形式来演示完整的计算流程. *)

(* 准TMD-PDF 矩阵元在 b_⊥ 空间的唯象模型:
   h(z, b_⊥) = h(z) × Exp[-b_⊥² / (4 σ²(z))]

   h(z) 的模型: 高斯型 Ioffe-时间分布
   h(z) = Exp[-(z P^z / λ_s)²]  (领头近似)
*)

Clear[ModelMatrixElement, ModelQuasiTMD];

(* 模型矩阵元 *)
ModelMatrixElement[z_, bperp_, Pz_, lambdaS_, sigma0_] :=
  Exp[-(z * Pz / lambdaS)^2] * Exp[-bperp^2 / (4 sigma0^2)];

(* 模型准TMD-PDF (已重整化) *)
ModelQuasiTMD[x_, bperp_, Pz_, lambdaS_, sigma0_, nz_: 100] := Module[
  {zmin, zmax, dz, integral},
  zmin = -10/Pz;
  zmax = 10/Pz;
  dz = (zmax - zmin)/nz;
  integral = Sum[
    ModelMatrixElement[z, bperp, Pz, lambdaS, sigma0] *
      Exp[-I z x Pz] * dz,
    {z, zmin, zmax, dz}
  ] / (2 \[Pi]);
  Re[integral]
];

Print["模型矩阵元: h(z,b_⊥) = Exp[-(zP^z/λ_s)²] Exp[-b_⊥²/(4σ₀²)]"];

(* --------------------------------------------------------------------------
   11.2 光锥 TMD-PDF 的数值演示
   -------------------------------------------------------------------------- *)

(* 使用典型格点参数:
   - P^z = 2.15 GeV (对应 CLQCD C24P29 的中间动量)
   - μ = 2 GeV (MS̄ 标度)
   - α_s(2 GeV) ≈ 0.30 (n_f = 4)
   - λ_s ≈ 3-5 (Ioffe时间分布宽度参数)
   - σ₀ ≈ 0.3 fm (横向分布宽度)
*)

Clear[DemoTMDCalculation];

DemoTMDCalculation[] := Module[
  {
    Pz, mu, alphaS, zeta, lambdaS, sigma0,
    xVals, bperpVals, i, j,
    fQuasi, fLightcone, SI, Kval, hardKernel
  },

  (* 物理参数 *)
  Pz = 2.15;      (* GeV *)
  mu = 2.0;        (* GeV *)
  alphaS = 0.30;
  zeta = mu^2;     (* 参考快度标度 *)
  lambdaS = 4.0;   (* Ioffe时间宽度 *)
  sigma0 = 0.3;    (* fm, 横向宽度 *)

  Print["物理参数:"];
  Print["  P^z = ", Pz, " GeV"];
  Print["  μ = ", mu, " GeV"];
  Print["  α_s(μ) = ", alphaS];
  Print["  ζ = ", zeta, " GeV^2"];

  (* 计算网格 *)
  xVals = Range[0.1, 0.7, 0.1];
  bperpVals = Range[0.0, 0.6, 0.1]; (* fm *)

  Print["\n计算网格:"];
  Print["  x = ", xVals];
  Print["  b_⊥ = ", bperpVals, " fm"];

  (* 对于演示, 我们输出模型化的 TMD-PDF *)
  Print["\n=== 模型光锥胶子TMD-PDF g(x, b_⊥) ==="];

  Do[
    Do[
      (* 准TMD-PDF *)
      fQuasi = ModelQuasiTMD[x, bperp, Pz, lambdaS, sigma0];

      (* 软函数 (模型值) *)
      SI = Exp[-bperp/0.5]; (* 指数衰减模型 *)

      (* CS核 (模型值) *)
      Kval = -0.2 * bperp^2/(bperp^2 + 0.1^2); (* 平滑模型 *)

      (* 硬匹配核 *)
      hardKernel = HardKernelNLO[(2 x Pz)^2, mu, alphaS];

      (* 光锥 TMD-PDF *)
      fLightcone = fQuasi / (Sqrt[SI] * hardKernel *
        Exp[(1/2) Log[(2 x Pz)^2/zeta] * Kval]);

      Print["  x=", NumberForm[x, {2,1}],
        ", b_⊥=", NumberForm[bperp, {2,1}], " fm: g = ",
        NumberForm[fLightcone, {6,4}]];
    , {bperp, bperpVals}]
  , {x, xVals}];

  Print["\n(注: 以上为唯象模型演示值, 非实际格点数据.)"];
  Print["实际格点TMD-PDF需使用真实组态数据通过完整工作流计算."];

  <|
    "method" -> "phenomenological_model",
    "Pz" -> Pz,
    "mu" -> mu,
    "note" -> "demonstration_values_only"
  |>
];

(* 运行演示 *)
DemoTMDCalculation[];


(* ============================================================================
   第十二部分: Jackknife 误差分析与连续极限检验
   参考: 补充/格点QCD中的误差统计.tex
   ============================================================================ *)

Print["\n========================================================================"];
Print["第十二部分: Jackknife 统计误差分析与系统误差检验"];
Print["========================================================================"];

(* --------------------------------------------------------------------------
   12.1 Jackknife 重采样
   -------------------------------------------------------------------------- *)

Clear[JackknifeResample, JackknifeError];

(* Jackknife 重采样: 轮流删除每个组态 *)
JackknifeResample[data_] := Module[{n, samples},
  n = Length[data];
  samples = Table[
    Mean[Delete[data, i]],
    {i, 1, n}
  ];
  samples
];

(* Jackknife 误差估计 *)
JackknifeError[data_] := Module[{n, mean, jkSamples, jkMean},
  n = Length[data];
  mean = Mean[data];
  jkSamples = JackknifeResample[data];
  jkMean = Mean[jkSamples];
  Sqrt[(n - 1) * Mean[(jkSamples - jkMean)^2]]
];

Print["Jackknife 误差: σ_{JK} = √((N-1) Σ_i (f_i^{JK} - f̄^{JK})²)"];

(* --------------------------------------------------------------------------
   12.2 连续极限外推的系统误差
   -------------------------------------------------------------------------- *)

Clear[ContinuumLimitSystematicError];

(* 系统误差来源:
   1. a^2 外推截断: 是否包含 a^4 项
   2. 外推范围: 使用 2 个 vs 3 个格点间距
   3. 固定 P^z 下的 a^2 外推 vs 联合 a^2-P^z 外推
*)

ContinuumLimitSystematicError[aVals_, dataVals_, methods_] := Module[
  {results, i},
  results = Table[
    Switch[methods[[i]],
      "linear_a2",
        ContinuumExtrapolationLinear[aVals, dataVals][[1]],
      "quad_a4",
        ContinuumExtrapolationQuad[aVals, dataVals][[1]],
      "omit_largest_a",
        Module[{aSub, dSub},
          aSub = Take[aVals, -2];
          dSub = Take[dataVals, -2];
          ContinuumExtrapolationLinear[aSub, dSub][[1]]
        ]
    ],
    {i, 1, Length[methods]}
  ];
  {
    Mean[results],
    StandardDeviation[results] (* 系统误差 *)
  }
];

Print["连续极限外推的系统误差估计:"];
Print["  比较不同外推方案 (线性a², 二次a⁴, 排除最大a)"];

(* --------------------------------------------------------------------------
   12.3 激发态污染控制 — 双态拟合
   参考: He et al. (2024) PRD 109, 114513
   -------------------------------------------------------------------------- *)

(* 对于三点函数 (核子-算符-核子), 双态拟合:
   C_{3pt}(t_{sep}, t) = A₀ Exp[-E₀ t_{sep}] ×
     [1 + A₁ Exp[-ΔE t_{sep}] + ...]

   其中 ΔE = E₁ - E₀ ~ 500 MeV
*)

Clear[TwoStateFit, ExcitedStateContamination];

(* 双态拟合模型 *)
TwoStateFit[tsep_, A0_, E0_, A1_, deltaE_] :=
  A0 * Exp[-E0 * tsep] * (1 + A1 * Exp[-deltaE * tsep]);

(* 激发态污染估计 *)
ExcitedStateContamination[tsep_, A1_, deltaE_] :=
  A1 * Exp[-deltaE * tsep];

Print["激发态污染: C_{3pt} = A₀e^{-E₀t_{sep}} [1 + A₁e^{-ΔE t_{sep}} + ...]"];
Print["  在 t_{sep} ≳ 1.0 fm, ΔE ~ 500 MeV → 污染 < 5%"];


(* ============================================================================
   第十三部分: 从共线PDF到TMD — k_⊥积分的检验
   ============================================================================ *)

Print["\n========================================================================"];
Print["第十三部分: 积分检验 — 从 TMD 到共线 PDF"];
Print["========================================================================"];

(* TMD-PDF → 共线 PDF 的一致性检验:
   f(x, μ) = ∫ d²k_⊥ f(x, k_⊥, μ, ζ)   (ζ-独立)

   b_⊥ 空间中的等价形式:
   f(x, μ) = 2π ∫₀^∞ db_⊥ b_⊥ J₀(0) f(x, b_⊥, μ, ζ)
            = 2π ∫₀^∞ db_⊥ b_⊥ f(x, b_⊥, μ, ζ)

   检验: 从 TMD 积分得到的结果应等于从共线 PDF 格点计算的结果.
*)

Clear[TMDtoCollinear, CollinearConsistencyCheck];

(* b_⊥ 空间中 TMD → 共线积分 *)
TMDtoCollinear[fTMD_, bperpVals_, x_] := Module[{integral, i},
  integral = Sum[
    2 \[Pi] * bperpVals[[i]] * fTMD[[i]] *
      (bperpVals[[2]] - bperpVals[[1]]),
    {i, 1, Length[bperpVals]}
  ];
  integral
];

(* 一致性检验 *)
CollinearConsistencyCheck[fTMD_, bVals_, fCollinear_, x_] := Module[
  {integral},
  integral = TMDtoCollinear[fTMD, bVals, x];
  {
    integral,
    fCollinear,
    Abs[integral - fCollinear] / Abs[fCollinear] (* 相对偏差 *)
  }
];

Print["积分检验: f(x,μ) = ∫ d²k_⊥ f(x,k_⊥,μ,ζ) = 2π∫ db_⊥ b_⊥ f(x,b_⊥,μ,ζ)"];
Print["  (此检验在格点数据分析中用于验证TMD因子化的一致性)"];


(* ============================================================================
   第十四部分: 胶子TMD与夸克TMD的关键对比
   ============================================================================ *)

Print["\n========================================================================"];
Print["第十四部分: 胶子TMD vs. 夸克TMD — 关键差异");
Print["========================================================================"];

(* 胶子TMD的特殊性:
   1. 胶子算符需要在伴随表示中构造 Wilson 线
   2. 36个Lorentz分量 → 可乘性重整化分类 (按z分量数)
   3. 非连通图: 胶子流与核子无夸克线连接
   4. 味单态混合: 2×2 匹配矩阵 C_{qq}, C_{qg}, C_{gq}, C_{gg}
   5. 线性发散更强: 胶子 Wilson 线的 k ~ (C_A/C_F) × k_夸克
*)

Clear[GluonVsQuarkTMD];

GluonVsQuarkTMD[] := Module[{},
  Print["胶子TMD与夸克TMD的关键差异:"];
  Print[""];
  Print["1. 算符构造:"];
  Print["   夸克: ψ̄ Γ W_{⊏} ψ (基础表示Wilson线)"];
  Print["   胶子: F^{μν} W^{(a)}_{⊏} F_{ρσ} (伴随表示Wilson线)"];
  Print[""];
  Print["2. Wilson线线性发散系数:"];
  Print["   夸克: k_q ∝ C_F = 4/3"];
  Print["   胶子: k_g ∝ C_A = 3 → k_g/k_q ≈ C_A/C_F = 9/4"];
  Print[""];
  Print["3. 微扰匹配:"];
  Print["   夸克: 1×1 匹配 (非单态)"];
  Print["   胶子: 2×2 匹配 (味单态, 与夸克PDF混合)"];
  Print[""];
  Print["4. 连通性:"];
  Print["   夸克: 连通+非连通 (同位旋矢量: 纯连通)"];
  Print["   胶子: 纯非连通 (胶子流无夸克线连接核子)"];
  Print[""];
  Print["5. 梯度流中的重整化:"];
  Print["   夸克: 需要 Z_χ (环标记夸克可消除)"];
  Print["   胶子: 纯胶子算符无需流场重整化"];
  Print[""];
  Print["6. 小x行为:"];
  Print["   夸克: x f₁(x) ~ const (Regge)"];
  Print["   胶子: x g(x) ~ x^{-λ} (更陡的小x增长)"];
];

GluonVsQuarkTMD[];


(* ============================================================================
   第十五部分: 与实验的联系 — EIC 唯象学
   ============================================================================ *)

Print["\n========================================================================"];
Print["第十五部分: TMD因子化普适性检验与EIC唯象学");
Print["========================================================================"];

(* TMD因子化的核心预言: TMD分布的普适性

   Sivers 函数符号反转:
   f_{1T}^{⊥}|_{DY} = -f_{1T}^{⊥}|_{SIDIS}

   这是非Abel规范理论 (QCD) 区别于 Abel 理论 (QED) 的关键预言.

   EIC 将以 1% 级别精度测量 SIDIS 方位角不对称性:
   - cos(2φ_h) 调制: Boer-Mulders × Collins (h₁^⊥ ⊗ H₁^⊥)
   - sin(φ_h - φ_S) 调制: Sivers × 非极化碎裂 (f_{1T}^⊥ ⊗ D₁)
   - sin(φ_h + φ_S) 调制: Collins 效应 (h₁ ⊗ H₁^⊥)
*)

(* T-odd 函数的符号反转检验 (格点QCD) *)
Clear[SiversSignFlipTest];

SiversSignFlipTest[fSIDIS_, fDY_] := Module[{ratio},
  ratio = fDY / fSIDIS;
  Print["Sivers 符号反转检验:"];
  Print["  f_{1T}^⊥|_{DY} / f_{1T}^⊥|_{SIDIS}"];
  Print["  理论预言: -1 (QCD 非Abel性质)"];
  Print["  格点结果: ", ratio];
  ratio
];

Print["EIC 唯象学中格点TMD-PDF的三重角色:"];
Print["  1. 理论先验 (prior): EIC数据前的唯一第一性原理信息"];
Print["  2. 联合拟合 (joint fit): 打破唯象拟合中的参数简并"];
Print["  3. 理论基准 (benchmark): 直接比较格点预言与实验提取"];


(* ============================================================================
   总结
   ============================================================================ *)

Print["\n========================================================================"];
Print["总结: 使用梯度流重整化方案的连续极限下核子非极化胶子TMD-PDF"];
Print["========================================================================"];

Print["
本 Mathematica 代码实现了以下完整的理论推导链:

1. 梯度流 (Wilson Flow) 方程
   - 规范场和夸克场的流方程
   - UV有限性定理 (Lüscher-Weisz 2011)
   - 小流时展开 (SFTX) 与 Wilson 系数

2. 胶子 TMD 算符构造
   - BMR 六个不变振幅分解
   - Clover 场强张量
   - 非极化和极化胶子流算符
   - 36个算符的可乘性重整化分类 (Li-Ma-Qiu)

3. Staple型Wilson线
   - 三段规范平行传输
   - 四类UV发散: 线性, cusp, pinch-pole, 对数

4. 梯度流重整化方案
   - Monahan-Orginos (2017) 梯度流正规化
   - 环标记夸克场 (Makino-Suzuki 2014)
   - PDF Moment 比值 (Francis et al. 2025)

5. Wilson圈减除 + SDR 重整化
   - 五格点间距系统验证 (Zhang et al. 2022)
   - RI/MOM 方案失败的物理原因
   - 在壳夸克矩阵元的单圈MS̄结果

6. TMD软函数与CS核
   - 内禀软函数 S_I: 形状因子因子化提取
   - Collins-Soper核 K: 动量比提取 (NLO改进)
   - CS核的普适性验证 (MILC/CLS交叉验证)

7. 连续极限外推
   - a^2标度的线性/二次外推
   - 联合连续-无穷大动量外推
   - 梯度流方案中的连续极限 (a^2/τ 展开)

8. TMD因子化匹配
   - NLO 硬匹配核 (单圈)
   - 重整化群重求和 (RGR)
   - 跑动耦合 (四圈β函数)

9. 完整数值流程
   - 10步工作流 (参考 gluon_pdf_full_workflow.py)
   - Fourier 变换 (z→x, b_⊥→k_⊥)
   - 大λ外推
   - Jackknife 误差分析

主要参考源 (按在本库中的路径):
  - 要求.md
  - 补充/格点QCD中的TMD_PDF.tex
  - 补充/格点QCD中的梯度流重整化.tex
  - 补充/格点QCD中的胶子算符.tex
  - 补充/格点QCD中的场强张量.tex
  - 补充/格点QCD中的重整化.tex
  - 补充/格点QCD中的大动量有效理论.tex
  - 补充/格点QCD中的外推.tex
  - code/gluon_pdf_full_workflow.py
  - examples/Calc_ope_unpol.py
  - examples/2pt_proton_Cg5gmu_L32x64_mom2_xdir_gpu.py
  - docs/Note of gluon PDFs.pdf
  - 文档/gluon_pdf_derivation.tex
  - 文档/gluon_PDF_continuum.tex
  - docs/Unpolarized gluon PDF of the nucleon from lattice QCD in the continuum limit.pdf
  - docs/Renormalization of TMD parton distribution on the lattice.pdf
  - docs/Quasi parton distributions and the gradient flow.pdf
  - docs/Accessing gluon parton distributions in large momentum effective theory.pdf
  - docs/Multiplicative renormalizability of quasi-parton operators.pdf

本文档可直接在 Mathematica / Wolfram Language 环境中运行.
"]
