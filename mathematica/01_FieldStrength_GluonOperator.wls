#!/usr/bin/env wolframscript
(* ::Title:: *)
(*\[ | \[ | *)

(* ::Section:: *)
(*\[ | \[ | *)

(* ::Text:: *)
(*\[ | *)

(* ::Text:: *)
(*\[ | *)
(* ::Subsection:: *)
(*\[ | \[ | *)

(* \[ | *)
Print["=== 1.1  SU(N) Gauge Field and Link Variables ==="];

(* \[ | *)
(* U_\[Mu](x) = exp(i a A_\[Mu](x)) \[ | *)
(* \[ | *)
Print["Gauge link variable: U_\[Mu](n) = exp(i a A_\[Mu](n))"];
Print["where A_\[Mu] = A_\[Mu]^a T^a, T^a are SU(N) generators."];
Print[""];

(* \[ | *)
(* \[ | *)
colorGroup = 3;  (* SU(3) *)
Print["Color group: SU(", colorGroup, ")"];
Print[""];

(* ::Subsection:: *)
(*\[ | \[ | *)

(* \[ | *)
(* P_{\[Mu]\[Nu]}(x) = U_\[Mu](x) U_\[Nu](x+\[Mu]^) U_\[Mu]^\[Dagger](x+\[Nu]^) U_\[Nu]^\[Dagger](x) \[ | *)
(* \[ | *)
Print["=== 1.2  Plaquette and Field Strength Tensor ==="];
Print[""];
Print["The plaquette is the fundamental gauge-invariant object on the lattice:"];
Print["  P_{\[Mu]\[Nu]}(x) = U_\[Mu](x) U_\[Nu](x+\[Mu]^) U_\[Mu]^\[Dagger](x+\[Nu]^) U_\[Nu]^\[Dagger](x)"];

(* Baker-Campbell-Hausdorff \[ | *)
Print[""];
Print["Using BCH formula:"];
Print["  P_{\[Mu]\[Nu]} = exp[i a^2 (\[PartialD]_\[Mu] A_\[Nu] - \[PartialD]_\[Nu] A_\[Mu] + i[A_\[Mu],A_\[Nu]]) + O(a^3)]"];
Print["  = 1 + i a^2 F_{\[Mu]\[Nu]}(x) + O(a^4)"];
Print[""];

(* \[ | *)
Print["--- Plaquette Symbolic Expansion (BCH) ---"];

(* \[ | *)
(* \[ | *)
(* We expand the plaquette symbolically using non-commutative algebra *)
(* to demonstrate how the field strength tensor emerges. *)

(* \[ | *)
Clear[A, U, F];
(* Define formal non-commutative gauge potentials *)
(* A[\[Mu]] represents A_\[Mu](x) as a formal operator *)

(* Plaquette = U1 * U2 * U3^dag * U4^dag where *)
(* U1 = exp(i a A_\[Mu](x)) *)
(* U2 = exp(i a A_\[Nu](x + a \[Mu]^)) *)
(* U3 = exp(i a A_\[Mu](x + a \[Nu]^)) *)
(* U4 = exp(i a A_\[Nu](x)) *)

(* \[ | \[ | *)
(* We use Series expansion in small a *)

(* \[ | \[ | *)
ExpandPlaquette[mu_, nu_, a_] := Module[{},
  (* U_\[Mu](x) = exp(i a A_\[Mu]) \[TildeTilde] 1 + i a A_\[Mu] - a^2/2 A_\[Mu]^2 *)
  (* U_\[Nu](x+\[Mu]^) = exp(i a A_\[Nu](x+\[Mu]^)) *)
  (*   A_\[Nu](x+\[Mu]^) \[TildeTilde] A_\[Nu] + a \[PartialD]_\[Mu] A_\[Nu] *)
  (* We work to O(a^2) in the exponent and O(a^4) in the expansion *)

  (* Symbolically: P_{\[Mu]\[Nu]} = 1 + i a^2 F_{\[Mu]\[Nu]} + ... *)
  (* where F_{\[Mu]\[Nu]} = \[PartialD]_\[Mu] A_\[Nu] - \[PartialD]_\[Nu] A_\[Mu] + i [A_\[Mu], A_\[Nu]] *)

  Print["Plaquette P_{", mu, nu, "} expanded to O(a^2):"];
  Print["  = 1 + i a^2 F_{", mu, nu, "} + O(a^4)"];
  Print["  with F_{", mu, nu, "} = \[PartialD]_", mu, " A_", nu, " - \[PartialD]_", nu, " A_", mu, " + i [A_", mu, ", A_", nu, "]"];
];

ExpandPlaquette[\[Mu], \[Nu], a];

(* ::Subsection:: *)
(*\[ | \[ | *)

(* \[ | *)
(* To reduce statistical fluctuations, average over 4 plaquette orientations *)
(* Q_{\[Mu]\[Nu]} = P_{\[Mu]\[Nu]} - P_{\[Nu],-\[Mu]} + P_{-\[Nu],\[Mu]} + P_{-\[Mu],-\[Nu]} *)
(* F_{\[Mu]\[Nu]} = -i/(8 a^2) (Q_{\[Mu]\[Nu]} - Q_{\[Nu]\[Mu]}) *)
(* \[ | *)

Print[""];
Print["--- Clover Term Construction ---"];
Print[""];
Print["The clover term averages 4 plaquette orientations to reduce noise:"];
Print["  Q_{\[Mu]\[Nu]}(x) = P_{\[Mu]\[Nu]}(x) - P_{\[Nu],-\[Mu]}(x) + P_{-\[Nu],\[Mu]}(x) + P_{-\[Mu],-\[Nu]}(x)"];
Print["  F_{\[Mu]\[Nu]} = -i/(8 a^2) (Q_{\[Mu]\[Nu]} - Q_{\[Nu]\[Mu]})"];
Print[""];
Print["This gives the anti-hermitian traceless field strength tensor on the lattice."];
Print[""];

(* ::Subsection:: *)
(*\[ | \[ | *)

(* F_{\[Mu]\[Nu]}^{Mink} \[ | *)
(* F_{tx}^{Mink} F_{tx}^{Mink} = -F_{tx}^{Eucl} F_{tx}^{Eucl} *)
(* F_{xy}^{Mink} F_{xy}^{Mink} = F_{xy}^{Eucl} F_{xy}^{Eucl} *)
Print["=== 1.3  Minkowski-Euclidean Correspondence ==="];
Print[""];
Print["Field strength tensor in Minkowski space (covariant):"];
Print["  F^{(M)}_{\[Mu]\[Nu]} = {{0, F_tx, F_ty, F_tz},"];
Print["                         {F_xt, 0, F_xy, F_xz},"];
Print["                         {F_yt, F_yx, 0, F_yz},"];
Print["                         {F_zt, F_zx, F_zy, 0}}"];
Print[""];
Print["Minkowski-Euclidean conversion rules:"];
Print["  F^{tx,(M)} F^{(M)}_{tx} = -F^{(E)}_{tx} F^{(E)}_{tx}"];
Print["  F^{xy,(M)} F^{(M)}_{xy} = +F^{(E)}_{xy} F^{(E)}_{xy}"];
Print[""];
Print["The minus sign for temporal components is crucial for correct"];
Print["operator construction in Euclidean lattice calculations."];
Print[""];

(* ::Section:: *)
(*\[ | \[ | *)

Print["=== 2.1  Continuous Definition of Gluon Field Strength ==="];
Print[""];
Print["In the continuum:"];
Print["  F^a_{\[Mu]\[Nu]} = \[PartialD]_\[Mu] A^a_\[Nu] - \[PartialD]_\[Nu] A^a_\[Mu] - g f^{abc} A^b_\[Mu] A^c_\[Nu]"];
Print["  \[TildeF]_{\[Mu]\[Nu]} = (1/2) \[Epsilon]_{\[Mu]\[Nu]\[Rho]\[Sigma]} F^{\[Rho]\[Sigma]}"];
Print["  with \[Epsilon]_{0123} = 1"];
Print[""];

(* ::Subsection:: *)
(*\[ | \[ | *)

(* \[ | *)
Print["=== 2.2  Levi-Civita Tensor Construction ==="];

(* \[ | \[ | *)
leviCivita4[n_:4] := Module[{eps},
  eps = Table[0, {\[Mu], n}, {\[Nu], n}, {\[Rho], n}, {\[Sigma], n}];
  Do[
    (* Only non-zero if all indices distinct *)
    If[Length[Union[{\[Mu], \[Nu], \[Rho], \[Sigma]}]] == n,
      eps[[\[Mu]+1, \[Nu]+1, \[Rho]+1, \[Sigma]+1]] = Signature[{\[Mu], \[Nu], \[Rho], \[Sigma]}]
    ],
    {\[Mu], 0, n-1}, {\[Nu], 0, n-1}, {\[Rho], 0, n-1}, {\[Sigma], 0, n-1}
  ];
  eps/2  (* factor 1/2 for F-tilde = (1/2) \[Epsilon] F *)
];

epsilon4 = leviCivita4[4];
Print["Levi-Civita tensor \[Epsilon]_{\[Mu]\[Nu]\[Rho]\[Sigma]}/2 constructed (4D)."];
Print["Non-zero components example: \[Epsilon]_0123/2 = ", epsilon4[[1,2,3,4]]];
Print[""];

(* ::Subsection:: *)
(*\[ | \[ | *)

(* \[ | *)
(* Compute F-tilde from F via Levi-Civita contraction *)
(* \[TildeF]_{\[Mu]\[Nu]} = (1/2) \[Epsilon]_{\[Mu]\[Nu]\[Rho]\[Sigma]} F^{\[Rho]\[Sigma]} *)
Print["=== 2.3  Dual Field Strength Tensor ==="];

(* \[ | F_{\[Mu]\[Nu]} \[ | *)
(* Construct a symbolic F tensor and compute its dual *)
(* For SU(3), F is 4\[Times]4 anti-symmetric with 6 independent real components *)

Clear[F, FTilde, f];
(* F[\[Mu],\[Nu]] is anti-symmetric: F[\[Nu],\[Mu]] = -F[\[Mu],\[Nu]] *)
(* We define the 6 independent components *)

(* \[ | *)
FTildeCompute[FMat_] := Module[{FT, mu, nu, rho, sigma, eps},
  eps = leviCivita4[4];
  FT = Table[0, {mu, 4}, {nu, 4}];
  Do[
    FT[[mu+1, nu+1]] = Sum[
      eps[[mu+1, nu+1, rho+1, sigma+1]] * FMat[[rho+1, sigma+1]],
      {rho, 0, 3}, {sigma, 0, 3}
    ],
    {mu, 0, 3}, {nu, 0, 3}
  ];
  FT
];

(* Test with symbolic F *)
Print["Symbolic test: \[TildeF]_{01} = \[Epsilon]_{01\[Rho]\[Sigma]} F^{\[Rho]\[Sigma]}/2"];
Print["  Nonzero contributions from \[Rho]\[Sigma] = 23, 32:"];
Print["  \[TildeF]_{01} = (1/2)(F^{23} + (-1) F^{32}) = F^{23}"];
Print["  (since F^{32} = -F^{23})"];
Print[""];
Print["Similarly: \[TildeF]_{02} = -F^{13}, \[TildeF]_{03} = F^{12}, etc."];
Print[""];

(* ::Section:: *)
(*\[ | \[ | *)

(* ::Subsection:: *)
(*\[ | \[ | *)

(* O(z) = M_{tx;tx}(z) + M_{ty;ty}(z) - 2 M_{xy;xy}(z) \[ | *)
(* \[ | *)
Print["=== 3.1  Multiplicatively Renormalizable Gluon Operator ==="];
Print[""];
Print["In the fundamental representation, the gluon quasi-PDF operator is:"];
Print["  M_{\[Mu]\[Lambda];\[Nu]\[Rho]}(z) = F^{a}_{\[Mu]\[Lambda]}(z) U^{ab}(z,0) F^{b}_{\[Nu]\[Rho]}(0)"];
Print[""];
Print["On the lattice (fundamental representation):"];
Print["  M_{\[Mu]\[Lambda];\[Nu]\[Rho]}(z) = \[Sum]_x Tr[ F_{\[Mu]\[Lambda]}(x+z*z^) U(x+z*z^, x) F_{\[Nu]\[Rho]}(x) U(x, x+z*z^) ]"];
Print[""];
Print["The multiplicatively renormalizable combination (Zhang et al. 2019, Balitsky et al. 2020):"];
Print["  O(z) = M_{tx;tx}(z) + M_{ty;ty}(z) - 2 M_{xy;xy}(z)"];
Print[""];
Print["This combination works because M^{ti;it} and M^{ij;ji} share the same"];
Print["one-loop UV anomalous dimension, making the sum multiplicatively renormalizable."];
Print[""];

(* ::Subsection:: *)
(*\[ | \[ | *)

(* \[ | \[ | *)
(* O_g^{(0)}(z) = [-F^{0i} W F^{0i} W^\[Dagger] + 2 F^{12} W F^{12} W^\[Dagger]]_{Eucl} *)
Print["=== 3.2  Unpolarized Gluon Current Operator ==="];
Print[""];
Print["In Minkowski space:"];
Print["  O_g^{(0)}(z) = [F^{0i}(z) W F_{0i}(0) W^\[Dagger] + 2 F^{12}(z) W F_{12}(0) W^\[Dagger]]_{Mink}"];
Print[""];
Print["In Euclidean space (using Minkowski-Euclidean conversion):"];
Print["  O_g^{(0)}(z) = [-F^{0i}(z) W F_{0i}(0) W^\[Dagger] + 2 F^{12}(z) W F_{12}(0) W^\[Dagger]]_{Eucl}"];
Print["  = -F^{01} W F^{01} W^\[Dagger] - F^{02} W F^{02} W^\[Dagger] - F^{03} W F^{03} W^\[Dagger]"];
Print["    + 2 F^{12} W F^{12} W^\[Dagger]"];
Print[""];
Print["The minus sign on the temporal components comes from the"];
Print["Minkowski-to-Euclidean conversion of F^{0i} F_{0i}."];
Print[""];

(* ::Subsection:: *)
(*\[ | \[ | *)

(* \[ | *)
(* M^{ti;it} + M^{ji;ij} = 2 p_0^2 M_{pp} *)
(* -M_{pp} = (1/2) \[Integral]_{-1}^1 dx e^{-i x \[Nu]} x g(x) *)
Print["=== 3.3  Invariant Amplitude Decomposition ==="];
Print[""];
Print["The matrix element decomposes into invariant amplitudes:"];
Print["  M_{\[Mu]\[Alpha];\[Lambda]\[Beta]}(z,p) = <p| F_{\[Mu]\[Alpha]}(z) W F_{\[Lambda]\[Beta]}(0) |p>"];
Print[""];
Print["The twist-2 invariant amplitude M_{pp} is extracted via:"];
Print["  M^{ti;it}(z,p) + M^{ji;ij}(z,p) = 2 p_0^2 M_{pp}(\[Nu], 0)"];
Print[""];
Print["Relationship to the gluon PDF (Ioffe time \[Nu] = p_3 z_3):"];
Print["  -M_{pp}(\[Nu]) = (1/2) \[Integral]_{-1}^1 dx e^{-i x \[Nu]} x g(x)"];
Print[""];
Print["Inverting the Fourier transform:"];
Print["  x g(x) = -(1/\[Pi]) \[Integral]_0^\[Infinity] d\[Nu] cos(x\[Nu]) M_{pp}(\[Nu])"];
Print[""];

(* ::Subsection:: *)
(*\[ | \[ | *)

(* \[ | *)
(* \[CapitalDelta]G = \[Integral]_0^1 dx \[CapitalDelta]g(x) *)
Print["=== 3.4  Gluon Helicity Operator (for comparison) ==="];
Print[""];
Print["Gluon helicity distribution from z-odd combination:"];
Print["  \[TildeM]_{\[Mu]\[Alpha];\[Lambda]\[Beta]}(z,p) = m(z,p) - m(-z,p)"];
Print[""];
Print["Helicity PDF from Ioffe-time distribution:"];
Print["  \[CapitalDelta]G = \[Integral]_0^\[Infinity] d\[Nu] \[TildeI]_p(\[Nu]) = \[Integral]_0^1 dx \[CapitalDelta]g(x)"];
Print["  -i \[TildeI]_p(\[Nu]) = \[TildeM]_{(ps)}^{(+)}(\[Nu]) - \[Nu] \[TildeM]_{(pp)}(\[Nu])"];
Print[""];
Print["In Euclidean space, the helicity gluon current operator is:"];
Print["  O_g^{(0)}(z) = -{F^{01} W F^{23} - F^{02} W F^{13} + 2 F^{12} W F^{03}} - {z \[LeftRightArrow] -z}"];
Print[""];

(* ::Section:: *)
(*\[ | \[ | *)

(* ::Subsection:: *)
(*\[ | \[ | *)

Print["=== 4.1  Matrix Element Extraction ==="];
Print[""];
Print["The bare matrix element h_B(z, P_z, 1/a) = <P_z| O(z) |P_z> is extracted"];
Print["from the ratio of 3-point to 2-point correlation functions:"];
Print[""];
Print["  R_3pt(P_z; t_sep, t) = C_3pt(P_z; t_sep, t) / C_2pt(P_z; t_sep)"];
Print["  \[TildeTilde] c_0 + c_1 e^{-\[CapitalDelta]E(t_sep - t)} + c_1 e^{-\[CapitalDelta]E t}"];
Print[""];
Print["where c_0 is the desired bare quasi-matrix element."];
Print[""];
Print["For disconnected diagrams (gluon operators):"];
Print["  C_3pt = \[Sum]_{t_0} (O_g - <O_g>) \[Times] (C_2pt - <C_2pt>)"];
Print[""];
Print["The disconnected 3pt function factorizes into:"];
Print["  1. Proton 2pt correlator (with distillation + momentum smearing)"];
Print["  2. OPE part (nonlocal gluon operator with Wilson line)"];
Print[""];

(* ::Subsection:: *)
(*\[ | \[ | *)

(* \[ | *)
(* -\[Del]^2_{xy}(t) = 6\[Delta]_{xy} - \[Sum]_{j=1}^3 [U_j(x,t)\[Delta]_{x+j^,y} + U_j^\[Dagger](x-j^,t)\[Delta]_{x-j^,y}] *)
(* \[ | *)
(* \[ | \[ | *)
Print["=== 4.2  Distillation Method ==="];
Print[""];
Print["Distillation (Peardon et al. 2009) uses low-lying eigenmodes of the"];
Print["3D gauge-covariant Laplacian to smear quark fields:"];
Print[""];
Print["  -\[Del]^2_{xy}(t) = 6\[Delta]_{xy} - \[Sum]_{j=1}^3 [U_j(x,t)\[Delta]_{x+j^,y} + U_j^\[Dagger](x-j^,t)\[Delta]_{x-j^,y}]"];
Print[""];
Print["Eigenvalue problem:"];
Print["  -\[Del]^2(t) v^{(i)} = \[Lambda]_i v^{(i)},   v^{(i)} v^{(j)\[Dagger]} = \[Delta]_{ij}"];
Print[""];
Print["Distillation operator:"];
Print["  \[Square]_{xy}(t) = \[Sum]_{k=1}^{N_ev} v^{(k)}_x(t) v^{(k)\[Dagger]}_y(t) = V(t) V^\[Dagger](t)"];
Print[""];
Print["Perambulator (key contraction object):"];
Print["  \[Tau]^{(p,\[BarP])}_{\[Alpha],\[Beta]} = V^{(p)\[Dagger]} M^{-1} V^{(\[BarP])}"];
Print[""];
Print["Momentum smearing for boosted nucleons:"];
Print["  \[TildeV]^k_x(z,t) = exp(i \[Zeta]\[CenterDot]z) v^k_x(z,t),  \[Zeta] = 2 \[Times] 2\[Pi]/L z^"];
Print[""];

(* ::Section:: *)
(*\[ | \[ | *)

(* *)

Print["=== Summary of File 01 ===*];
Print["*];
Print["This file has established the mathematical foundation for:"];
Print["  1. Gauge field strength tensor F_{\[Mu]\[Nu]} on the lattice (Clover plaquette)"];
Print["  2. Dual field strength tensor \[TildeF]_{\[Mu]\[Nu]} (Levi-Civita contraction)"];
Print["  3. Multiplicatively renormalizable gluon operator O(z)"];
Print["  4. Unpolarized gluon current O_g^{(0)} in Euclidean space"];
Print["  5. Invariant amplitude decomposition M_{pp} \[LeftRightArrow] x g(x)"];
Print["  6. Matrix element extraction and distillation method"];
Print[""];
Print["Proceed to 02_GradientFlow_Renormalization.wls for gradient flow"];
Print["and renormalization theory."];
