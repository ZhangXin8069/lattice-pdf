# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Lattice QCD research project by Zhang Xin. The goal is to compute the **unpolarized gluon Parton Distribution Function (PDF) of the proton** from first-principles lattice QCD calculations.

**Theoretical framework**: Large Momentum Effective Theory (LaMET), which connects Euclidean lattice correlators to light-cone PDFs via quasi-PDFs and perturbative matching.

**Method**: The gluon PDF involves **disconnected diagrams**. The three-point correlation function is decoupled into two independent parts:
1. **Proton 2pt correlator** — computed with momentum-smeared distillation
2. **OPE (Operator Product Expansion) part** — the nonlocal gluon operator

The matrix element to compute is Eq(20) from the internal note; the OPE decomposition follows Eq(25).

## Key References

**Primary** (in `docs/`):
- `Note of gluon PDFs.pdf` — internal theory note with the core equations
- `Unpolarized gluon PDF of the nucleon from lattice QCD in the continuum limit.pdf`
- `First Nucleon Gluon PDF from Large Momentum Effective Theory.pdf`

**Textbooks** (in `books/`):
- `Quantum Chromodynamics on the Lattice.pdf` — Gattringer & Lang
- `An Introduction to Quantum Field Theory.pdf` — Peskin & Schroeder
- `INTRODUCTION TO LATTICE QCD.pdf` — lecture notes

**All other papers in `docs/`** cover: quasi-PDFs, pseudo-PDFs, renormalization, distillation, TMDs, gluon operators, Wilson lines, and related topics.

## Repository Structure

| Directory | Content |
|---|---|
| `code/` | **Main analysis code**: `gluon_pdf_full_workflow.py` (~1900 lines) — a systematic reconstruction of the full 10-step computation pipeline from gauge configurations to light-cone PDF |
| `examples/` | **Original GPU code** by donghx: `2pt_proton_Cg5gmu_L32x64_mom2_xdir_gpu.py` (proton 2pt), `Calc_ope_unpol.py` (OPE), with CuPy GPU utilities. Contains symlinks to remote data on the cluster (`/public/group/lqcd/donghx/`) |
| `refer/` | **Theory analysis & workflow document**: `理论解析与工作流.tex` — comprehensive LaTeX document with the theoretical derivation and computational workflow |
| `reports/` | **Beamer slides**: `gluon_pdf_slides.tex` — presentation slides for the gluon PDF results |
| `补充/` | **Supplementary LaTeX notes** (~35 topics in Chinese) covering all aspects of lattice QCD: fermion actions, gauge fields, Wilson lines, correlation functions, renormalization, distillation, Monte Carlo methods, error analysis, Feynman diagrams, etc. |
| `文档/` | **Documentation**: gluon PDF derivation (`gluon_pdf_derivation.tex`), continuum extrapolation (`gluon_PDF_continuum.tex`), internal note (`note_of_gluon_PDFs.tex`) |
| `汇报/` | **Reports**: gluon quasi-operator construction on the lattice, lattice calculation methods |
| `代码/` | **Code analysis**: `code_analysis.tex` — analysis of the GPU computation code |
| `要求.md` | Project requirements/specification document (the starting point) |

## Computational Workflow (`code/gluon_pdf_full_workflow.py`)

The 10-step pipeline from gauge configurations to light-cone PDF:

1. Read gauge configurations (gauge links U_μ(x))
2. Construct field strength tensor F_{μν} (Clover plaquette)
3. Construct dual field strength tensor F̃_{μν} (Levi-Civita contraction)
4. Construct nonlocal gluon OPE operator O_{μν}(z) with Wilson lines
5. Implement Distillation framework (eigenvectors + perambulator + VVV)
6. Momentum-smeared proton two-point correlation function
7. Matrix element extraction h(z, P_z)
8. Fourier transform → quasi-PDF g̃(x, P_z)
9. Perturbative matching → light-cone PDF g(x, μ)
10. Jackknife statistical error analysis

### Python Dependencies

```python
numpy, cupy        # GPU acceleration (falls back to numpy if cupy unavailable)
opt_einsum          # Optimized tensor contractions (falls back to numpy.einsum)
```

The code is designed to run on HPC clusters with GPUs. Conditional imports handle environments where CuPy/opt_einsum are unavailable.

## LaTeX Compilation

All LaTeX documents use **XeLaTeX** with Chinese support (ctex package). The project uses MiKTeX/TeX Live with these key packages:
`ctex`, `xeCJK`, `fontspec`, `physics`, `braket`, `beamer` (with metropolis/ Madrid/ CambridgeUS themes)

**Compile a document**:
```bash
# Single pass
xelatex -interaction=nonstopmode -halt-on-error <file>.tex

# With latexmk (auto-reruns for cross-references)
latexmk -xelatex -interaction=nonstopmode -halt-on-error <file>.tex

# Clean auxiliary files
latexmk -c <file>.tex
```

**Check if required packages are available**:
```bash
kpsewhich <package>.sty   # e.g., kpsewhich ctex.sty
```

**Compile the slides** and check page count:
```bash
cd reports && latexmk -xelatex -interaction=nonstopmode gluon_pdf_slides.tex
cd reports && mdls -name kMDItemNumberOfPages gluon_pdf_slides.pdf
```

## Important Project Conventions

- **Chinese is the primary language** for documentation, notes, and comments. Code docstrings and variable names are in English. The `要求.md` file and all supplementary notes (`补充/`) are in Chinese.
- **Python 3** is the standard. The `.vscode/settings.json` configures `conda` as the package/environment manager.
- **Data files** (HDF5, npy, npz, h5) live on the HPC cluster, not in this repo — they are gitignored. The `examples/` directory has symlinks to cluster paths.
- **GPU code** in `examples/` is the original (donghx's); `code/gluon_pdf_full_workflow.py` is a systematic refactored version.
- **The `.gitignore`** excludes: build artifacts, LaTeX auxiliary files (`.aux`, `.log`, `.toc`, `.nav`, `.snm`, `.xdv`, `.fdb_latexmk`, `.fls`), compiled binaries, HDF5/npy data, Python cache, and editor config.
- **Placeholder files named `PASS`** in `code/` and `log/` mark directories that need to be tracked by git but were otherwise empty.
