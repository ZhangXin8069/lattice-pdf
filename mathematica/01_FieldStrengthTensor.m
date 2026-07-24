(* ::Package:: *)

(* ============================================================================
   第1部分: 格点上的场强张量 (Field Strength Tensor on the Lattice)
   ============================================================================

   参考文献:
   - Gattringer & Lang, "Quantum Chromodynamics on the Lattice", Ch.2-3
   - 补充/格点QCD中的场强张量.tex
   - 文档/gluon_pdf_derivation.tex, Sec.5
   - code/gluon_pdf_full_workflow.py, Steps 2-3
   ============================================================================ *)

BeginPackage["FieldStrengthTensor`"];

FieldStrengthTensorDerivation;
GaugeLinkToFieldStrength;
CloverTerm;
PlaquetteExpansion;
MinkowskiEuclideanRelation;
DualFieldStrength;
FmuNuMatrix;

Begin["`Private`"];

(* --------------------------------------------------------------------------
   1.1 格点规范链接与规范势的关系
   --------------------------------------------------------------------------
   格点规范链接 U_mu(n) 与连续规范势 A_mu(x) 的关系:
   U_mu(n) = P exp[ig \int_0^a ds A_mu(n a + s \hat{mu})]
           = exp[ig a A_mu(n a) + O(a^2)]

   在微扰展开中: U_mu(n) = 1 + i g a A_mu(n) - 1/2 g^2 a^2 A_mu(n)^2 + ...
   -------------------------------------------------------------------------- *)

(* Plaquette: 最基本的Wilson圈, 1x1方格 *)
Plaquette[U_, x_, mu_, nu_] := Module[{},
  (* P_{mu,nu}(x) = U_mu(x) U_nu(x+mu) U_mu^dag(x+nu) U_nu^dag(x) *)
  U[mu, x] . U[nu, x + UnitVector[mu]] .
  ConjugateTranspose[U[mu, x + UnitVector[nu]]] .
  ConjugateTranspose[U[nu, x]]
];

(* --------------------------------------------------------------------------
   1.2 Baker-Campbell-Hausdorff展开: 从Plaquette到场强张量
   --------------------------------------------------------------------------
   将plaquette展开到O(a^2):
   P_{mu,nu} = 1 + i a^2 F_{mu,nu}(x) + O(a^3)

   其中场强张量:
   F_{mu,nu}^a = \partial_mu A_nu^a - \partial_nu A_mu^a - g f^{abc} A_mu^b A_nu^c

   等价地: F_{mu,nu} = -i [D_mu, D_nu],  D_mu = \partial_mu + i g A_mu
   -------------------------------------------------------------------------- *)

(* BCH展开: 验证Plaquette展开到O(a^2) *)
PlaquetteExpansion[verbose_:False] := Module[{result},
  (* 在连续极限下:
     P_{mu,nu} = exp[i a^2 (\partial_mu A_nu - \partial_nu A_mu + i[A_mu, A_nu]) + O(a^3)]
               = 1 + i a^2 F_{mu,nu} + O(a^3)
  *)
  result = {
    "Plaquette" -> "P_{mu,nu}(x) = U_mu(x) U_nu(x+mu) U_mu^dag(x+nu) U_nu^dag(x)",
    "BCH_Expansion" -> "P_{mu,nu} = exp[i a^2 F_{mu,nu}(x) + O(a^3)]",
    "FieldStrength" -> "F_{mu,nu}(x) = \partial_mu A_nu - \partial_nu A_mu - g f^{abc} A_mu^b A_nu^c",
    "CovariantForm" -> "F_{mu,nu} = -i [D_mu, D_nu], D_mu = \partial_mu + i g A_mu"
  };

  If[verbose,
    Print["=== Plaquette展开 ==="];
    Print["P_{mu,nu} = U_mu(x) U_nu(x+mu) U_mu^dag(x+nu) U_nu^dag(x)"];
    Print["使用BCH公式展开:"];
    Print["  U_mu(x) = exp(i a A_mu(x))"];
    Print["  U_nu(x+mu) = exp(i a A_nu(x) + i a^2 \partial_mu A_nu(x) + ...)"];
    Print["  => P_{mu,nu} = 1 + i a^2 F_{mu,nu}(x) + O(a^3)"];
  ];

  result
];

(* --------------------------------------------------------------------------
   1.3 Clover项: 数值稳定的场强张量构造
   --------------------------------------------------------------------------
   为减小统计涨落, 通过平均四个定向plaquette来构造clover项:
   Q_{mu,nu}(x) = P_{mu,nu}(x) - P_{nu,-mu}(x) + P_{-nu,mu}(x) + P_{-mu,-nu}(x)

   场强张量由clover项的反称部分给出:
   F_{mu,nu} = -i/(8 a^2) (Q_{mu,nu} - Q_{nu,mu})

   这一构造自动保证了F_{mu,nu}的反对称性: F_{mu,nu} = -F_{nu,mu}
   -------------------------------------------------------------------------- *)

CloverTerm[verbose_:False] := Module[{clover, fieldStrength},
  clover = {
    "Qmunu" -> "Q_{mu,nu}(x) = P_{mu,nu} - P_{nu,-mu} + P_{-nu,mu} + P_{-mu,-nu}",
    "FourPlaquettes" -> {
      "P_{mu,nu}"   -> "(mu, nu) 正向plaquette",
      "P_{nu,-mu}"  -> "(nu, -mu) plaquette",
      "P_{-nu,mu}"  -> "(-nu, mu) plaquette",
      "P_{-mu,-nu}" -> "(-mu, -nu) plaquette"
    }
  };

  fieldStrength = {
    "Fmunu" -> "F_{mu,nu} = -i/(8 a^2) (Q_{mu,nu} - Q_{nu,mu})",
    "Antisymmetry" -> "F_{mu,nu} = -F_{nu,mu}  (自动满足)",
    "NumberComponents" -> "6个独立分量: F_{01},F_{02},F_{03},F_{12},F_{13},F_{23}"
  };

  If[verbose,
    Print["=== Clover项构造 ==="];
    Print["Q_{mu,nu} = P_{mu,nu} - P_{nu,-mu} + P_{-nu,mu} + P_{-mu,-nu}"];
    Print["F_{mu,nu} = -i/(8 a^2) (Q_{mu,nu} - Q_{nu,mu})"];
  ];

  {clover, fieldStrength}
];

(* --------------------------------------------------------------------------
   1.4 Minkowski与Euclidean空间的场强张量关系
   --------------------------------------------------------------------------
   Minkowski空间: 度规 g_{mu,nu} = diag(1, -1, -1, -1)
   Euclidean空间: 度规 g_{mu,nu} = diag(1, 1, 1, 1)   (已旋转后)

   Wick旋转: t -> -i tau 导致:
   A_0^(M) = i A_0^(E)    (时间分量获得因子 i)
   A_i^(M) = A_i^(E)      (空间分量不变)

   因此场强张量的转换关系:
   F_{0i}^{(M)} = i F_{0i}^{(E)}
   F_{ij}^{(M)} = F_{ij}^{(E)}

   乘积关系 (关键于胶子流构造):
   F_{0i}^{(M)} F_{0i}^{(M)} = -F_{0i}^{(E)} F_{0i}^{(E)}
   F_{ij}^{(M)} F_{ij}^{(M)} =  F_{ij}^{(E)} F_{ij}^{(E)}
   -------------------------------------------------------------------------- *)

MinkowskiEuclideanRelation[verbose_:False] := Module[{relation},
  relation = {
    "WickRotation" -> "t -> -i\\tau,  A_0^(M) = i A_0^(E)",
    "FieldStrengthTransform" -> {
      "F_{0i}^{(M)} = i F_{0i}^{(E)}",
      "F_{ij}^{(M)} = F_{ij}^{(E)}"
    },
    "ProductRelations" -> {
      "F_{0i}^{(M)} F_{0i}^{(M)} = -F_{0i}^{(E)} F_{0i}^{(E)}",
      "F_{ij}^{(M)} F_{ij}^{(M)} =  F_{ij}^{(E)} F_{ij}^{(E)}"
    }
  };

  If[verbose,
    Print["=== Minkowski <-> Euclidean 关系 ==="];
    Print["Wick旋转: t -> -i\\tau"];
    Print["F_{0i}^{(M)} = i F_{0i}^{(E)},  F_{ij}^{(M)} = F_{ij}^{(E)}"];
    Print["F_{0i}^{(M)}F_{0i}^{(M)} = -F_{0i}^{(E)}F_{0i}^{(E)}"];
    Print["F_{ij}^{(M)}F_{ij}^{(M)} =  F_{ij}^{(E)}F_{ij}^{(E)}"];
  ];

  relation
];

(* --------------------------------------------------------------------------
   1.5 场强张量矩阵的显式形式
   --------------------------------------------------------------------------
   在Minkowski空间中 (F_{mu,nu} = \partial_mu A_nu - \partial_nu A_mu - g[A_mu, A_nu]):

   F_{mu,nu}^{(M)} = [[   0,   E_x,   E_y,   E_z ],
                      [ -E_x,     0,  -B_z,   B_y ],
                      [ -E_y,   B_z,     0,  -B_x ],
                      [ -E_z,  -B_y,   B_x,     0 ]]

   在Euclidean空间中 (所有分量实):
   F_{mu,nu}^{(E)} = [[   0,  E_x,  E_y,  E_z ],
                      [ -E_x,    0,  B_z, -B_y ],
                      [ -E_y, -B_z,    0,  B_x ],
                      [ -E_z,  B_y, -B_x,    0 ]]
   -------------------------------------------------------------------------- *)

FmuNuMatrix[] := Module[{},
  (* Euclidean空间中的场强张量矩阵 *)
  {
   "Minkowski" -> {
     {"0", "E_x", "E_y", "E_z"},
     {"-E_x", "0", "-B_z", "B_y"},
     {"-E_y", "B_z", "0", "-B_x"},
     {"-E_z", "-B_y", "B_x", "0"}
   },
   "Euclidean" -> {
     {"0", "E_x", "E_y", "E_z"},
     {"-E_x", "0", "B_z", "-B_y"},
     {"-E_y", "-B_z", "0", "B_x"},
     {"-E_z", "B_y", "-B_x", "0"}
   },
   "Note" -> "在Euclidean空间中, 电场分量获得因子i: E_i^(E) = i E_i^(M)"
  }
];

(* --------------------------------------------------------------------------
   1.6 对偶场强张量
   --------------------------------------------------------------------------
   对偶场强张量: \tilde{F}_{mu,nu} = 1/2 \epsilon_{mu,nu,rho,sigma} F^{rho,sigma}

   约定: \epsilon_{0123} = +1
   在Euclidean空间中, Levi-Civita张量的定义保持一致

   显式关系:
   \tilde{F}_{01} = F_{23},  \tilde{F}_{02} = -F_{13},  \tilde{F}_{03} = F_{12}
   \tilde{F}_{23} = F_{01},  \tilde{F}_{13} = -F_{02},  \tilde{F}_{12} = F_{03}
   -------------------------------------------------------------------------- *)

DualFieldStrength[verbose_:False] := Module[{dual},
  dual = {
    "Definition" -> "\\tilde{F}_{mu,nu} = 1/2 \\epsilon_{mu,nu,rho,sigma} F^{rho,sigma}",
    "Convention" -> "\\epsilon_{0123} = +1",
    "ExplicitRelations" -> {
      "\\tilde{F}_{01} = F_{23},  \\tilde{F}_{23} = F_{01}",
      "\\tilde{F}_{02} = -F_{13}, \\tilde{F}_{13} = -F_{02}",
      "\\tilde{F}_{03} = F_{12},  \\tilde{F}_{12} = F_{03}"
    }
  };

  If[verbose,
    Print["=== 对偶场强张量 ==="];
    Print["\\tilde{F}_{mu,nu} = 1/2 \\epsilon_{mu,nu,rho,sigma} F^{rho,sigma}"];
    Print["\\tilde{F}_{01}=F_{23}, \\tilde{F}_{02}=-F_{13}, \\tilde{F}_{03}=F_{12}"];
  ];

  dual
];

(* --------------------------------------------------------------------------
   主推导函数
   -------------------------------------------------------------------------- *)
FieldStrengthTensorDerivation[] := Module[{},
  Print["=== Part 1: 格点场强张量 ==="];

  Print["\n1.1 格点规范链接:"];
  Print["  U_mu(n) = P exp(i g \\int_0^a ds A_mu(n a + s \\hat{mu}))"];

  Print["\n1.2 Plaquette构造:"];
  Print["  P_{mu,nu}(x) = U_mu(x) U_nu(x+\\hat{mu}) U_mu^dag(x+\\hat{nu}) U_nu^dag(x)"];
  Print["  BCH展开: P_{mu,nu} = 1 + i a^2 F_{mu,nu}(x) + O(a^3)"];

  Print["\n1.3 Clover项 (场强张量格点定义):"];
  Print["  Q_{mu,nu} = P_{mu,nu} - P_{nu,-mu} + P_{-nu,mu} + P_{-mu,-nu}"];
  Print["  F_{mu,nu} = -\\frac{i}{8a^2}(Q_{mu,nu} - Q_{nu,mu})"];
  Print["  该构造包含4个plaquette的平均, 减小统计涨落"];

  Print["\n1.4 Minkowski <-> Euclidean 关键转换:"];
  Print["  F_{0i}^{(M)} F_{0i}^{(M)} = -F_{0i}^{(E)} F_{0i}^{(E)}"];
  Print["  F_{ij}^{(M)} F_{ij}^{(M)} =  F_{ij}^{(E)} F_{ij}^{(E)}"];

  Print["\n1.5 参考源:"];
  Print["  - Gattringer & Lang, 'QCD on the Lattice', Sec.2-3"];
  Print["  - 补充/格点QCD中的场强张量.tex"];
  Print["  - 文档/gluon_pdf_derivation.tex, Eq.(16-19)"];
  Print["  - code/gluon_pdf_full_workflow.py, Steps 2-3"];
];

End[];
EndPackage[];
FieldStrengthTensorDerivation[];
