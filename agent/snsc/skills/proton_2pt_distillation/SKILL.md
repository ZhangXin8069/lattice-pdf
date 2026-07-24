# Proton 2pt Correlator via Distillation

Compute the proton two-point correlation function using the distillation
(Laplacian-Heaviside) method with momentum smearing on momentum-smeared
perambulators.

## Physics

The proton 2pt correlator in the distillation framework is:

  C₂(t_snk, t_src) = Σ_{x_snk,x_src} ⟨0| χ(t_snk, x_snk) χ̄(t_src, x_src) |0⟩

where the interpolator χ = ε_{abc} [u_a^T Cγ₅ d_b] u_c (proton operator).

In the distillation subspace (Nev lowest eigenmodes of the 3D Laplace
operator -∇²):

  C₂^{il}(t_snk, t_src) = Φ_{abc}(t_snk) · τ^{ad}_{gi} · (ΓτΓ)^{be}_{gj}
                         · τ^{cf}_{il} · Φ^*_{def}(t_src)    [Direct]
                       - Φ_{abc}(t_snk) · τ^{af}_{gl} · (ΓτΓ)^{be}_{gj}
                         · τ^{cd}_{ij} · Φ^*_{def}(t_src)    [Exchange]

where:
  Φ_{abc}(P) = Σ_x e^{-iP·x} ε_{ijk} v_i^a(x) v_j^b(x) v_k^c(x)  [VVV block]
  τ(t_snk, t_src) = perambulator (quark propagator in distillation subspace)
  Γ = Cγ₅γ₄ (interpolator Dirac structure, element="_Cg5g4")

Parity projection separates positive/negative parity states:
  P₊ = ½(γ₀ + γ₄),  P₋ = ½(γ₀ - γ₄)
  C₂^{pp}(t_snk, t_src) = Tr[P₊ · C₂(t_snk, t_src)]

Boundary condition sign correction:
  pp: C₂^{pp}(t_snk, t_src) *= -1  when t_snk < t_src
  pm: C₂^{pm}(t_snk, t_src) *= -1  when t_snk > t_src

Effective mass (cosh method, periodic BC):
  a·m_eff(t) = arccosh( [C(t-1) + C(t+1)] / [2 C(t)] )
  m_eff(t) [GeV] = a·m_eff(t) × 0.1973 / a[fm]

## Data Dependencies

Required input files (all per-conf_id, cluster paths):

| Data | Path Pattern | Format |
|------|-------------|--------|
| eigenvectors | `eig_dir/eigvecs_t{t:03d}_{conf_id}` | binary float64, (Nev, Nx³, 3) complex |
| perambulators | `peram_dir/perams.{conf_id}.{d_src}.{t_src}` | binary float64, (Nt, 4, 4, Nev1, Nev1) |
| gauge configs | `conf_dir/cfg_{conf_id}.lime` (optional, for OPE) | ILDG binary |

Optional reference data (for validation):
| Data | Path Pattern |
|------|-------------|
| reference raw | `ref_dir/twopt_slice_pp_..._contract_conf{id}.npy` |
| reference pp | `ref_dir/twopt_slice_pp_..._nopol_ss_conf{id}.npy` |

## Output Files

| File | Shape | Description |
|------|-------|-------------|
| `VVV_Nev1{Nev1}_Px{Px}Py{Py}Pz{Pz}_conf{id}.npy` | (Nt, Nev1, Nev1, Nev1) | VVV baryon block (cached) |
| `twopt_slice_pp_Px{Px}Py{Py}Pz{Pz}_eginphase{ms}{elem}_contract_conf{id}.npy` | (Nt, Nt, 4, 4) | raw contraction matrix |
| `twopt_slice_pp_Px{Px}Py{Py}Pz{Pz}_eginphase{ms}{elem}_nopol_ss_conf{id}.npy` | (Nt, Nt) | parity-projected pp correlator |
| `meff_Pz{Pz}_conf{id}.npz` | — | effective mass data + 1D C2pt |
| `meff_Pz{Pz}_{elem}_conf{id}.pdf` | — | effective mass plot |

## Ensemble Presets

All paths use absolute cluster paths (matching ~/ symlinks):

| Ensemble | Nt | Nx | Nev1 | mom_smear | phase |
|----------|----|----|------|-----------|-------|
| L24x72 | 72 | 24 | 100 | -2 | 2 |
| L32x64 | 64 | 32 | 100 | 3 | -3 |
| L32x96 | 96 | 32 | 100 | 2 | -2 |
| L36x108 | 108 | 36 | 200 | 2 | -2 |
| L48x96 | 96 | 48 | 200 | 4 | -4 |
| L48x144 | 144 | 48 | 200 | 2 | -2 |
| L64x128 | 128 | 64 | 200 | 2 | -2 |

## References

- Peardon et al., Phys.Rev.D 80 (2009) 054506 — distillation method
- 2pt_proton_Cg5gmu_L24x72_mom2_zdir_dcu.py — reference DCU implementation
- snsc/main.py —arp-unified pipeline (see --analysis-type proton-2pt)
- agent/lamet-agent — correlator_analysis stage for effective mass fitting
