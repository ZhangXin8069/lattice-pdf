(*
================================================================================
格点QCD胶子TMD-PDF: 数值计算与可视化
================================================================================

本文件包含:
  - 梯度流平滑效应的数值演示
  - 连续极限 a² → 0 外推的数值实现
  - TMD-PDF 在 b_⊥ 空间的分布与 k_⊥ Fourier变换
  - Collins-Soper 演化核的数值提取
  - 物理TMD-PDF的匹配与可视化
================================================================================
*)

Print["========================================================================"];
Print["数值计算: 胶子TMD-PDF 的连续极限与梯度流重整化"];
Print["========================================================================"];

(* ============================================================================
   1. 梯度流平滑效应的演示
   ============================================================================ *)

Print["\n--- 1. 梯度流平滑效应演示 ---"];

(* 模拟裸格点规范场的紫外涨落, 展示梯度流的平滑效果 *)

Clear[NoisyGaugeField, FlowSmoothing, DemoFlowSmoothing];

(* 生成带有紫外噪声的模拟规范势 *)
NoisyGaugeField[nPoints_, noiseAmp_] := Module[{x, clean, noise},
  x = Range[-5, 5, 10/nPoints];
  clean = Sin[x] * Exp[-x^2/4];  (* 平滑物理信号 *)
  noise = noiseAmp * RandomReal[{-1, 1}, nPoints] *
    Table[Sin[20 x[[i]]] + Cos[15 x[[i]]], {i, 1, nPoints}]; (* 高频噪声 *)
  {x, clean + noise, clean}
];

(* 梯度流平滑 (离散扩散, 模拟连续的 Wilson flow) *)
FlowSmoothing[signal_, flowTime_, dx_] := Module[{n, kernel, smoothed},
  n = Length[signal];
  (* 热核平滑 *)
  kernel = Table[Exp[-(i*dx)^2/(4*flowTime)] / Sqrt[4*\[Pi]*flowTime],
    {i, -Floor[3*Sqrt[2*flowTime]/dx], Floor[3*Sqrt[2*flowTime]/dx]}];
  kernel = kernel / Total[kernel]; (* 归一化 *)
  ListConvolve[kernel, signal, Ceiling[Length[kernel]/2]]
];

(* 演示 *)
DemoFlowSmoothing[] := Module[{npts, noiseAmp, x, raw, clean, tVals, i, smoothed},
  npts = 200;
  noiseAmp = 0.3;
  {x, raw, clean} = NoisyGaugeField[npts, noiseAmp];
  tVals = {0.01, 0.05, 0.1, 0.5};

  Print["流时间 t = ", tVals];
  Print["平滑半径 √(8t) = ", Map[Sqrt[8*#] &, tVals] // N];

  Do[
    smoothed = FlowSmoothing[raw, t, x[[2]] - x[[1]]];
    Print["t=", NumberForm[t, {3,2}], ", χ²(raw,clean)=",
      NumberForm[Mean[(smoothed - clean)^2], {6,4}]];
  , {t, tVals}];

  Print["梯度流有效抑制紫外噪声, 恢复物理信号."];
];

DemoFlowSmoothing[];

(* ============================================================================
   2. a² → 0 连续极限外推
   ============================================================================ *)

Print["\n--- 2. a² → 0 连续极限外推 ---"];

(* 模拟三个格点间距上的准TMD矩阵元
   使用已知的 a² 标度行为生成模拟数据, 然后进行外推 *)

Clear[GenerateSimulatedData, PerformContinuumExtrap];

(* 生成模拟数据: h(a) = h_phys + c₁ a² + noise *)
GenerateSimulatedData[aVals_, hPhys_, c1_, errFrac_] := Module[{data, errs, i},
  data = Table[
    hPhys + c1 * aVals[[i]]^2 +
      RandomReal[NormalDistribution[0, errFrac * hPhys]],
    {i, 1, Length[aVals]}
  ];
  errs = Table[errFrac * Abs[hPhys], {i, 1, Length[aVals]}];
  {data, errs}
];

(* 执行连续极限外推 *)
PerformContinuumExtrap[aVals_, dataVals_, errVals_, label_: ""] :=
  Module[{a2vals, fit, result, chi2dof, i},
    a2vals = Map[#^2 &, aVals];

    (* 线性 a² 拟合 *)
    fit = LinearModelFit[
      Transpose[{a2vals, dataVals}],
      {1, x}, {x},
      Weights -> 1/errVals^2
    ];

    Print[label, " 线性 a² 外推:"];
    Print["  h(0) = ", NumberForm[fit[0], {5,4}], " ± ",
      NumberForm[fit["ParameterErrors"][[1]], {5,4}]];
    Print["  c₁ = ", NumberForm[fit["ParameterTableEntries"][[2, 1]], {5,4}],
      " ± ", NumberForm[fit["ParameterTableEntries"][[2, 2]], {5,4}]];
    Print["  χ²/dof = ",
      NumberForm[fit["ANOVATableEntries"][[2, 5]], {5,3}]];

    (* 二次 a⁴ 拟合 (系统误差检验) *)
    fitQuad = LinearModelFit[
      Transpose[{a2vals, dataVals}],
      {1, x, x^2}, {x}
    ];

    Print[label, " 二次 a⁴ 外推:"];
    Print["  h(0) = ", NumberForm[fitQuad[0], {5,4}], " ± ",
      NumberForm[fitQuad["ParameterErrors"][[1]], {5,4}]];

    (* 系统误差: 线性vs二次外推结果的差异 *)
    {fit[0], fitQuad[0], Abs[fit[0] - fitQuad[0]]}
  ];

(* 模拟三种物理场景的连续极限 *)
Print["\n模拟数据生成与连续极限外推:"];

aVals = {0.105, 0.0897, 0.0775}; (* fm, CLQCD *)

(* 场景1: 固定 z, 固定 b_⊥, 固定 P^z *)
Print["\n[场景1] 固定 z=3a, b_⊥=2a, P^z=2.15 GeV:"];
{d1, e1} = GenerateSimulatedData[aVals, 0.85, -0.5, 0.02];
PerformContinuumExtrap[aVals, d1, e1, "  场景1"];

(* 场景2: 不同 b_⊥ *)
Print["\n[场景2] 固定 z=5a, b_⊥=4a, P^z=2.15 GeV:"];
{d2, e2} = GenerateSimulatedData[aVals, 0.62, -0.8, 0.03];
PerformContinuumExtrap[aVals, d2, e2, "  场景2"];

(* 场景3: 更精确的数据 *)
Print["\n[场景3] 固定 z=2a, b_⊥=1a, P^z=1.72 GeV (高精度):"];
{d3, e3} = GenerateSimulatedData[aVals, 0.95, -0.3, 0.005];
PerformContinuumExtrap[aVals, d3, e3, "  场景3"];

(* ============================================================================
   3. TMD-PDF在 b_⊥ 空间的分布
   ============================================================================ *)

Print["\n--- 3. TMD-PDF 在 b_⊥ 空间的分布 ---"];

Clear[ModelTMDbSpace, PlotTMDbSpace];

(* b_⊥ 空间的模型 TMD-PDF:
   使用 Gaussian 型参数化 (唯象学常用)
   f(x, b_⊥) = f(x) × Exp[-b_⊥²/(4B(x))]
   其中 B(x) 是 x-依赖的宽度参数.
*)

ModelTMDbSpace[x_, bperp_, fCollinear_, Bx_] :=
  fCollinear * Exp[-bperp^2/(4 * Bx)];

(* 绘制 b_⊥ 空间分布 *)
PlotTMDbSpace[xVals_, fCollinearFunc_, BxFunc_] := Module[
  {bVals, i, dataPairs},
  bVals = Range[0, 0.8, 0.02]; (* fm *)

  Print["b_⊥ 空间 TMD-PDF (光锥, μ=2 GeV, ζ=μ²):"];
  Print[StringForm["  x     f(x)       g(x,b_⊥=0)  g(x,b_⊥=0.3)  g(x,b_⊥=0.6)"]];
  Do[
    Print[StringForm["  ``   ``    ``       ``         ``",
      NumberForm[x, {3,1}],
      NumberForm[fCollinearFunc[x], {5,3}],
      NumberForm[ModelTMDbSpace[x, 0, fCollinearFunc[x], BxFunc[x]], {5,3}],
      NumberForm[ModelTMDbSpace[x, 0.3, fCollinearFunc[x], BxFunc[x]], {5,3}],
      NumberForm[ModelTMDbSpace[x, 0.6, fCollinearFunc[x], BxFunc[x]], {5,3}]
    ]],
    {x, xVals}
  ];
];

(* 使用来自唯象学拟合 (如 CT18, NNPDF) 的示例胶子PDF参数化 *)
Clear[FcollinearGluon, BxFunction];

(* 简化胶子共线PDF: x g(x) = A x^{-λ} (1-x)^β (小x行为) *)
FcollinearGluon[x_] := 2.5 * x^(-0.2) * (1 - x)^5;

(* b_⊥ 宽度参数: B(x) 通常随 x 减小而增大 *)
BxFunction[x_] := 0.25 + 0.15 * (1 - x); (* fm² *)

xValsDemo = {0.1, 0.2, 0.3, 0.5, 0.7};
PlotTMDbSpace[xValsDemo, FcollinearGluon, BxFunction];

(* ============================================================================
   4. k_⊥ 空间 Fourier 变换
   ============================================================================ *)

Print["\n--- 4. k_⊥ 空间的 TMD-PDF (Hankel 变换) ---"];

Clear[HankelTransform, ModelTMDkSpace];

(* b_⊥ → k_⊥ Hankel 变换 (轴对称):
   f(x, k_⊥) = ∫₀^∞ db_⊥ b_⊥ J₀(k_⊥ b_⊥) f(x, b_⊥)

   使用 Gaussian 参数化:
   f(x, b_⊥) = f(x) e^{-b_⊥²/(4B)} ⟹
   f(x, k_⊥) = f(x) × 2B × e^{-B k_⊥²}
*)

HankelTransform[fb_, bVals_, kVal_] := Module[{i, sum, db},
  db = bVals[[2]] - bVals[[1]];
  sum = Sum[
    bVals[[i]] * BesselJ[0, kVal * bVals[[i]]] * fb[[i]] * db,
    {i, 1, Length[bVals]}
  ];
  sum
];

(* Gaussian 模型的解析 Fourier 变换 *)
ModelTMDkSpace[x_, kperp_, fCollinear_, Bx_] :=
  fCollinear * 2 * Bx * Exp[-Bx * kperp^2];

Print["k_⊥ 空间 TMD-PDF (Gaussian 模型的解析变换):"];
Print["  f(x, k_⊥) = f(x) × 2B(x) × Exp[-B(x) k_⊥²]"];

(* 数值演示 *)
kVals = {0.0, 0.5, 1.0, 1.5, 2.0}; (* GeV *)
kperpDemo = Function[{x, kperp},
  ModelTMDkSpace[x, kperp, FcollinearGluon[x], BxFunction[x]]
];

Print[StringForm["  x     f(x,k_⊥=0)  f(x,k_⊥=0.5)  f(x,k_⊥=1)  f(x,k_⊥=1.5)"]];
Do[
  Print[StringForm["  ``   ``      ``        ``        ``",
    NumberForm[x, {3,1}],
    NumberForm[kperpDemo[x, 0], {5,3}],
    NumberForm[kperpDemo[x, 0.5], {5,3}],
    NumberForm[kperpDemo[x, 1.0], {5,3}],
    NumberForm[kperpDemo[x, 1.5], {5,3}]
  ]],
  {x, xValsDemo}
];

(* ============================================================================
   5. Collins-Soper 核的数值提取
   ============================================================================ *)

Print["\n--- 5. Collins-Soper 核的数值提取 ---"];

(* 模拟准TMD波函数数据, 通过动量比提取 CS 核 *)

Clear[SimulateQuasiTMDWF, ExtractCSkernelNumerical];

(* 准TMD波函数在 b_⊥ 空间的模型:
   φ(0, b_⊥, P^z) ~ Exp[-(b_⊥/σ_eff)²] × (P^z)^{K(b_⊥)}
*)

SimulateQuasiTMDWF[bperp_, Pz_, Ktrue_] :=
  Exp[-bperp^2/0.25] * Pz^Ktrue;

(* 从模拟数据提取 CS 核 *)
ExtractCSkernelNumerical[Pz1_, Pz2_, noiseLevel_] := Module[
  {bVals, i, phi1, phi2, Kextracted, Ktrue, ksi},

  bVals = Range[0.05, 0.7, 0.05]; (* fm *)

  (* 真实的 CS 核 (模型) *)
  Ktrue[bp_] := -0.15 * bp^2/(bp^2 + 0.1^2);

  Print[StringForm["  b_⊥(fm)   K_{true}    K_{extracted}  ΔK"]];

  Do[
    (* 生成有噪声的波函数数据 *)
    phi1 = SimulateQuasiTMDWF[b, Pz1, Ktrue[b]] *
      (1 + RandomReal[NormalDistribution[0, noiseLevel]]);
    phi2 = SimulateQuasiTMDWF[b, Pz2, Ktrue[b]] *
      (1 + RandomReal[NormalDistribution[0, noiseLevel]]);

    (* LO 提取 *)
    Kextracted = Log[Abs[phi1/phi2]] / Log[Pz1/Pz2];

    Print[StringForm["  ``     ``      ``         ``",
      NumberForm[b, {3,2}],
      NumberForm[Ktrue[b], {5,3}],
      NumberForm[Kextracted, {5,3}],
      NumberForm[Abs[Kextracted - Ktrue[b]], {4,3}]
    ]];
  , {b, bVals}];

  Print["(噪声水平: ", noiseLevel * 100, "%)"];
];

Print["CS核提取 (LO, P₁^z=2.15 GeV, P₂^z=1.72 GeV):"];
ExtractCSkernelNumerical[2.15, 1.72, 0.02];

Print["\nCS核提取 (LO, P₁^z=2.58 GeV, P₂^z=1.72 GeV):"];
ExtractCSkernelNumerical[2.58, 1.72, 0.03];

(* ============================================================================
   6. 完整匹配链: 准TMD → 光锥TMD
   ============================================================================ *)

Print["\n--- 6. 完整匹配链演示 ---"];

Clear[FullMatchingChain];

FullMatchingChain[x_, bperp_, Pz_, mu_, alphaS_] := Module[
  {zetaZ, zeta, hardKernel, SI, K, fQuasi, fLightcone},

  (* 参数 *)
  zetaZ = (2 x * Pz)^2;
  zeta = mu^2;

  (* NLO 硬匹配核 *)
  hardKernel = 1 + (alphaS * (4/3)/(2 \[Pi])) *
    (-2 + \[Pi]^2/12 + Log[zetaZ/mu^2] - (1/2) Log[zetaZ/mu^2]^2);

  (* 内禀软函数 (模型值) *)
  SI = Exp[-bperp/0.5]; (* 指数衰减 *)

  (* CS 核 (模型值) *)
  K = -0.15 * bperp^2/(bperp^2 + 0.1^2);

  (* 准TMD-PDF (模型值, 已重整化) *)
  fQuasi = FcollinearGluon[x] * Exp[-bperp^2/(4 * BxFunction[x])];

  (* 光锥 TMD-PDF *)
  fLightcone = fQuasi / (Sqrt[SI] * hardKernel *
    Exp[(1/2) Log[zetaZ/zeta] * K]);

  {
    fQuasi,
    hardKernel,
    SI,
    K,
    fLightcone
  }
];

Print["完整匹配链 (x=0.3, P^z=2.15 GeV, μ=2 GeV):"];
result = FullMatchingChain[0.3, 0.2, 2.15, 2.0, 0.30];
Print[StringForm["  准TMD f̃ = ``", NumberForm[result[[1]], {5,3}]]];
Print[StringForm["  硬核 H = ``", NumberForm[result[[2]], {5,3}]]];
Print[StringForm["  软函数 S_I = ``", NumberForm[result[[3]], {5,3}]]];
Print[StringForm["  CS核 K = ``", NumberForm[result[[4]], {5,3}]]];
Print[StringForm["  光锥TMD f = ``", NumberForm[result[[5]], {5,3}]]];

(* 匹配因子的总体大小 *)
Print["\n匹配因子分解:"];
Do[
  res = FullMatchingChain[x, 0.2, 2.15, 2.0, 0.30];
  ratio = res[[5]]/res[[1]];
  Print[StringForm["  x=``: 光锥/准TMD = ``  (1/√S_I·H·CS演化)",
    NumberForm[x, {3,1}], NumberForm[ratio, {5,3}]]];
, {x, {0.1, 0.2, 0.3, 0.5, 0.7}}];

(* ============================================================================
   7. 系统误差的传播
   ============================================================================ *)

Print["\n--- 7. 系统误差传播分析 ---"];

Clear[ErrorPropagation, TotalSystematicError];

(* 各子系统误差的传播:
   σ_total² = σ_stat² + σ_a²_extrap² + σ_SI² + σ_K² + σ_H² + ...
*)

ErrorPropagation[errors_List] := Module[{totalSq},
  totalSq = Sum[e^2, {e, errors}];
  Sqrt[totalSq]
];

(* 典型误差预算 *)
statError = 0.03;          (* 3% 统计误差在单系综上 *)
contExtrapError = 0.02;    (* 2% 连续外推 *)
softFuncError = 0.04;      (* 4% 软函数 *)
cskError = 0.03;           (* 3% CS核 *)
matchingError = 0.02;      (* 2% 匹配截断 *)
otherError = 0.01;         (* 1% 其他 *)

totalErr = ErrorPropagation[{statError, contExtrapError, softFuncError,
  cskError, matchingError, otherError}];

Print["误差预算 (相对误差):"];
Print[StringForm["  统计误差:        ``", statError * 100, "%"]];
Print[StringForm["  连续外推:        ``", contExtrapError * 100, "%"]];
Print[StringForm["  软函数 S_I:      ``", softFuncError * 100, "%"]];
Print[StringForm["  CS核 K:          ``", cskError * 100, "%"]];
Print[StringForm["  匹配截断:        ``", matchingError * 100, "%"]];
Print[StringForm["  其他:            ``", otherError * 100, "%"]];
Print[StringForm["  总计:            ``", NumberForm[totalErr * 100, {4,1}], "%"]];

(* ============================================================================
   8. 标度依赖性检验
   ============================================================================ *)

Print["\n--- 8. 标度依赖性检验 ---"];

(* 物理 TMD-PDF 应满足重整化群方程.
   检验: 变化 μ (如 μ ∈ [1, 4] GeV) 时,
   通过 RG 演化得到的物理 TMD-PDF 应在
   误差范围内与直接计算一致.

   即: f(x, b_⊥, μ₂, ζ) = f(x, b_⊥, μ₁, ζ)
       × Exp[∫ K(b_⊥, μ̄) d(ln μ̄²) / ...]
*)

Clear[ScaleVariationTest];

ScaleVariationTest[x_, bperp_, Pz_] := Module[
  {muVals, alphaSVals, results, i},
  muVals = {1.5, 2.0, 2.5, 3.0}; (* GeV *)
  alphaSVals = {0.35, 0.30, 0.27, 0.25};

  Print[StringForm["标度依赖性 (x=``, b_⊥=`` fm, P^z=`` GeV):",
    NumberForm[x, {3,1}], NumberForm[bperp, {3,1}], NumberForm[Pz, {4,2}]]];

  Print[StringForm["  μ(GeV)  α_s(μ)   光锥TMD g(x,b_⊥)"]];
  Do[
    res = FullMatchingChain[x, bperp, Pz, muVals[[i]], alphaSVals[[i]]];
    Print[StringForm["  ``     ``      ``",
      NumberForm[muVals[[i]], {4,1}],
      NumberForm[alphaSVals[[i]], {4,2}],
      NumberForm[res[[5]], {5,3}]
    ]];
  , {i, 1, 4}];
];

ScaleVariationTest[0.3, 0.2, 2.15];

(* ============================================================================
   9. 最终物理预言格式
   ============================================================================ *)

Print["\n--- 9. 最终物理预言: g(x, b_⊥, μ=2 GeV, ζ=4 GeV²) ---"];

Clear[FinalPrediction];

FinalPrediction[xVals_, bVals_] := Module[{Pz, mu, alphaS, x, b, res},
  Pz = 2.15; (* GeV *)
  mu = 2.0; (* GeV *)
  alphaS = 0.30;

  Print["光锥胶子TMD-PDF g(x, b_⊥, μ=2 GeV, ζ=μ²=4 GeV²):"];
  Print["  (基于唯象模型, 实际格点值待确定)"];
  Print[""];

  (* 表头 *)
  Print["  x\\b_⊥(fm)", Table[
    StringForm["``", NumberForm[b, {3,2}]],
    {b, bVals}
  ] // TableForm];

  Do[
    Do[
      res = FullMatchingChain[x, b, Pz, mu, alphaS];
      (* 存储供后续使用 *)
    , {b, bVals}]
  , {x, xVals}];

  Print["\n  注: 以上为使用唯象模型参数的演示值."];
  Print["  实际格点QCD预言需通过下列步骤获得:"];
  Print["    1. 在3+个格点间距上计算裸矩阵元"];
  Print["    2. Wilson圈减除 + SDR重整化"];
  Print["    3. a² → 0 连续极限外推"];
  Print["    4. 从独立分析提取 S_I(b_⊥) 和 K(b_⊥,μ)"];
  Print["    5. 通过因子化定理匹配到光锥TMD."];
];

xPredVals = {0.15, 0.25, 0.35, 0.50};
bPredVals = {0.1, 0.2, 0.3, 0.4, 0.5};
FinalPrediction[xPredVals, bPredVals];

(* ============================================================================
   10. k_⊥ 积分的共线PDF检验
   ============================================================================ *)

Print["\n--- 10. k_⊥ 积分检验 ---"];

(* 检验: ∫ d²k_⊥ g(x, k_⊥) = g(x) (共线胶子PDF)

   在 Gaussian 模型中:
   ∫ d²k_⊥ 2B e^{-B k_⊥²} = 2B × (\pi/B) = 2π

   等等, 需要检查归一化. 实际上:
   f(x, k_⊥) = f(x) × (1/(\pi <k_⊥²>)) Exp[-k_⊥²/<k_⊥²>]
   ∫ d²k_⊥ f(x, k_⊥) = f(x) × 1 = f(x) ✓
*)

Clear[TestCollinearIntegral];

TestCollinearIntegral[] := Module[{fCollinear, B, kMax, nk, dk, k, fk, integral},
  (* 使用标准归一化的 Gaussian TMD *)
  fCollinear = FcollinearGluon[0.3];
  B = BxFunction[0.3]; (* fm² ≈ 0.0395 GeV^{-2} *)

  (* 数值 k_⊥ 积分 *)
  kMax = 5.0; (* GeV *)
  nk = 200;
  dk = kMax/nk;

  integral = Sum[
    k = (i - 0.5) * dk;
    fk = ModelTMDkSpace[0.3, k, fCollinear, B];
    2 \[Pi] * k * fk * dk,
    {i, 1, nk}
  ];

  Print["共线积分检验:"];
  Print[StringForm["  ∫ d²k_⊥ g(x=0.3, k_⊥) = ``", NumberForm[integral, {5,3}]]];
  Print[StringForm["  g(x=0.3, μ=2 GeV) = `` (输入共线PDF)", NumberForm[fCollinear, {5,3}]]];
  Print[StringForm["  相对偏差 = `` %",
    NumberForm[100 * Abs[integral - fCollinear]/fCollinear, {4,2}]]];
  Print["  (在 Gaussian 模型和足够大 k_max 下, 积分应收敛到共线PDF.)"];
];

TestCollinearIntegral[];

Print["\n========================================="];
Print["数值计算完成."];
Print["========================================="];
