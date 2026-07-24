## Core Physics Objective

Compute the proton two-point correlation function at five non-zero boost momenta (Pz = −2,−3,−4,−5,−6 in units of 2π/Ls) using the **distillation (Laplacian-Heaviside) method** on configuration 46000 of the C24P29 ensemble (24³×72, β=6.20, a≈0.1052 fm). The interpolator uses the **Cg5g4 Dirac structure** (Γ = γ₃γ₁γ₄) for the third quark at the sink, with its conjugate Γ̄ = −Cg5g4 at the source.

## Key Revisions from Critique

### 1. Wick Contraction with External Spin Indices (Critical Fix)
The original formula contracted sink and source perambulator spin indices directly through `(Cg5g4)_{γ'γ}`, producing a scalar. The corrected formula introduces **external spin indices σ' (sink) and σ (source)**:

```
C_{σ'σ} = Σ Φ Φ* [τττ − exchange] (Cγ₅)(Cγ₅) (Cg5g4)_{σ'γ} (−Cg5g4)_{γ'σ}
```

This yields a genuine (Nt, Nt, 4, 4) complex matrix where σ' is the external sink spin and σ the external source spin. The projector acts as an outer product — not a trace — preserving the 4×4 structure.

### 2. Source Projector (Critical Fix)
The original plan omitted the source projector entirely. Since (Cg5g4)^† = −Cg5g4 in the DeGrand-Rossi basis, the conjugate operator at the source carries **Γ̄ = γ₄Γ^†γ₄ = −Cg5g4**. Both sink and source projectors are now explicitly included: `(Cg5g4)_{σ'γ}` at the sink and `(−Cg5g4)_{γ'σ}` at the source.

### 3. Exchange Term Index Consistency (Critical Fix)
After the i↔k index relabeling in the exchange term (which flips ε_{ijk} → −ε_{kji} and cancels the explicit minus sign), the "free" quark at both sink and source has the same perambulator spin indices (γ at sink, γ' at source) as the direct term. Therefore the same projector outer product applies consistently to both terms.

### 4. VVV Block Without Extra Smearing Phase (Substantial Fix)
The momentum-smearing phase e^{iζ·x} (ζ=−2×2π/24) is already baked into the perambulators during their construction. The VVV baryon block now uses **only the momentum-projection phase**:

```
Φ_{ijk}(Pz, t) = Σ_x e^{−i Pz·x} ε_{abc} v_i^a v_j^b v_k^c
```

Applying the smearing phase again at the sink would effectively shift the momentum to Pz−ζ, producing wrong correlators for all Pz except the coincidental Pz=−2.

### 5. Two-Sided Parity Projection (Substantial Fix)
Replaced single-sided `Tr[P₊ C]` with the correct two-sided form:

```
C_pp = Re Tr[P₊ @ C_raw @ P₊^†]   where P₊ = (I + γ₄)/2
```

A single-sided trace mixes parity channels when the source already carries a non-trivial spin structure (the −Cg5g4 projector). Two-sided projection correctly isolates positive parity at both source and sink.

### 6. Concrete Contraction Algorithm
Provided a feasible staged contraction using **numpy einsum with `optimize='optimal'`**, which automatically finds the memory-minimizing contraction path. Two formulations are given:
- **Single-einsum** (recommended for ≥128 GB RAM): direct + exchange terms as separate einsum calls.
- **Two-stage fallback** (for <64 GB RAM): pre-sum over i,i' before completing the diquark and third-quark contractions.

### 7. Data Validation at Load Time
Added explicit **shape and dtype checks** for eigenvectors and perambulators, with clear abort messages on mismatch. Added a **gamma basis verification step**: if perambulators were computed in a non-DGR convention, basis-rotation matrices must be inserted before all Dirac algebra.

## Data Flow

| Input | Source | Format |
|-------|--------|--------|
| Laplace eigenvectors (Nev=100) | /public/group/lqcd/eigensystem/.../46000/ | .npy, (13824, 100) per time slice |
| Momentum-smeared perambulators | /public/home/sunp/.../mz2_my0_mx0/46000/ | .npy/HDF5, (72,72,4,4,100,100) complex |
| Reference 2pt data | /public/group/lqcd/donghx/.../momsmear-2z/46000/ | .npy |

| Output | Description |
|--------|-------------|
| proton_C2pt_raw_Pz{nP}.npy | (72,72,4,4) complex — raw spin matrix with sink/source projectors |
| proton_C2pt_pp_Pz{nP}.npy | (72,72) real — two-sided parity-projected correlator |
| proton_meff_Pz{nP}.npy | (72,) float64 — effective mass from arccosh on folded correlator |
| plots/*.pdf | Log-scale C2pt, effective mass plateau, scatter comparison |

## DeGrand-Rossi Gamma Convention (Reference)

| Matrix | Definition | Properties |
|--------|-----------|------------|
| γ_i (i=1,2,3) | [[0, iσ_i], [−iσ_i, 0]] | Hermitian |
| γ₄ | [[0, I], [I, 0]] | Hermitian |
| γ₅ | [[I, 0], [0, −I]] | Hermitian, {γ₅,γ_μ}=0 |
| γ₀ | I_{4×4} | Identity |
| C = γ₂γ₄ | diag(−iσ₂, iσ₂) | C^† = −C, C^T = −C |
| Cγ₅ | diag(−iσ₂, −iσ₂) | Symmetric: (Cγ₅)^T = Cγ₅ |
| Cg5g4 = γ₃γ₁γ₄ | antidiag(−iσ₃σ₁, iσ₃σ₁) | (Cg5g4)^† = −Cg5g4; γ₄(Cg5g4)γ₄ = Cg5g4 |