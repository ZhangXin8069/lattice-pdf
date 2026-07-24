# Proton 2pt Distillation — Planner Summary

## Computational Strategy

The proton 2pt correlator is computed via the distillation method, which
projects quark fields onto the Nev=100 lowest eigenmodes of the 3D gauge-covariant
Laplace operator -∇².  This reduces the computational cost from O(V³) to O(Nev³)
while preserving the low-energy physics relevant for ground-state extraction.

### Three-Stage Computation

**Stage 1 — VVV Baryon Block Construction (per timeslice, per momentum)**

The VVV tensor Φ_{abc}(P) encodes the momentum-projected three-quark wavefunction
in the distillation subspace.  For each timeslice t ∈ [0, 71], we:
1. Read eigenvectors v_i^a(x) (shape Nev×Nx³×3) from binary files
2. Apply momentum smearing: ṽ = v × exp(-i·2π·P_smear·x/Nx)
3. For the target momentum P, compute the 6-term ε_{ijk} Levi-Civita contraction:
   Φ_{abc}(P) = Σ_x e^{-iP·x} ε_{ijk} ṽ_i^a(x) ṽ_j^b(x) ṽ_k^c(x)

The contraction uses x-layer slicing (24 layers) to manage memory for the
intermediate (Nev×Nx²×3) arrays.  Results are cached as .npy for reuse across
momentum values.

**Stage 2 — Wick Contraction (per source timeslice)**

The 2pt correlator is constructed from two Feynman diagrams:
- **Direct term**: quarks propagate from source to sink without crossing
- **Exchange term**: two quark lines are swapped (Fermi minus sign)

For each t_source ∈ [0, 71]:
1. Read perambulator τ(t_snk; t_src) — shape (Nt, 4, 4, Nev, Nev)
2. Compute Dirac-projected perambulator: ΓτΓ = contract("gh,thkbe,jk→tgjbe", Γ, τ, Γ)
3. For t_sink where 2 ≤ deltat ≤ 32:
   C₂ = Direct(Φ_snk, τ_snk, ΓτΓ_snk, τ_snk, Φ_src*) - Exchange(...)
   Output shape: (4, 4) — Dirac spin indices

Time separation cut [2, 32] ensures ground-state dominance while maintaining
adequate statistics for effective mass extraction.

**Stage 3 — Parity Projection & Effective Mass**

1. Parity projection: C₂^{pp} = Tr[P₊ · C₂], boundary sign correction
2. 1D correlator: C(deltat) = mean_{t_src} C₂^{pp}(t_src+deltat, t_src)
3. Cosh effective mass: a·m_eff(t) = arccosh((C(t-1)+C(t+1))/(2·C(t)))
4. Physical units: m_eff [GeV] = a·m_eff × 0.1973 / 0.1053

### Reference Validation

Generated results are compared against the reference data from:
`/public/group/lqcd/donghx/2pt_Result/beta6.20_mu-0.2770_ms-0.2400_L24x72/momsmear-2z/`

Comparison metrics: max absolute diff, mean absolute diff, relative diff,
fraction of elements within 1e-6 and 1e-10 tolerances.

### Resource Estimate

- VVV: ~25s per timeslice (CPU, 24×Nev³×Nx contractions) → ~30 min per momentum
- Wick contraction: ~5 min per t_source (72 sources) → ~6 h per momentum
- Total for 5 momenta: ~32 h (with VVV caching, subsequent momenta ~6 h each)
- Memory: ~4 GB (VVV block: 72×100³×16 bytes ~ 1.1 GB; perambulator ~ 0.9 GB)
