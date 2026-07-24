# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is Zhang Xin's Lattice QCD research project. The goal is to compute the **unpolarized gluon Parton Distribution Function (PDF) of the proton** from first-principles lattice QCD calculations.

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

**Paper library** (`docs/`): 80+ PDFs covering quasi-PDFs, pseudo-PDFs, renormalization, distillation, TMDs, gluon operators, Wilson lines, gradient flow, and related topics. Use `pdftotext` to extract text for reading.

## Repository Structure

| Directory | Content |
|---|---|
| `code/` | **Two workflow implementations**: `gluon_pdf_full_workflow.py` (~1900 lines, single-script pipeline) and `gluon_pdf_workflow.py` (~1580 lines, subcommand-based with MPI and ensemble presets) |
| `examples/` | **64 GPU/DCU scripts** by donghx: proton 2pt correlators across 6 lattice ensembles, 3 momentum directions, and 2 GPU backends. Plus utility modules: `Operator.py`, `gamma_matrix_cupy_DR.py`, `input_output_4_cupy.py`. Contains symlinks to remote data on the cluster (`/public/group/lqcd/donghx/`). |
| `docs/` | **80+ reference papers** (PDF) on LaMET, quasi-PDFs, distillation, renormalization, gradient flow, etc. |
| `agent/` | **AI agent submodules** and documentation: three git submodules for agentic lattice QCD workflows |
| `refer/` | `理论解析与工作流.tex` — comprehensive LaTeX document with theoretical derivation and computational workflow |
| `reports/` | **Beamer slides**: `gluon_pdf_slides.tex` (theoretical basis & workflow) and `gluon_pdf_continuum_beamer.tex` (continuum-limit results) |
| `补充/` | **34 supplementary LaTeX notes** (Chinese) covering all aspects of lattice QCD: fermion actions, gauge fields, Wilson lines, correlation functions, renormalization, distillation, Monte Carlo methods, error analysis, etc. |
| `文档/` | **Core documents**: gluon PDF derivation (`gluon_pdf_derivation.tex`), continuum extrapolation (`gluon_PDF_continuum.tex`), internal note (`note_of_gluon_PDFs.tex`) |
| `汇报/` | **Reports**: gluon quasi-operator construction on the lattice, lattice calculation methods |
| `代码/` | **Code analysis**: `code_analysis.tex` — analysis of the GPU computation code |

## Two Workflow Implementations in `code/`

The `code/` directory contains two independent implementations of the gluon PDF workflow. They compute the same physics but differ in architecture and intended use.

### `gluon_pdf_full_workflow.py` — Self-contained pipeline

A single-script implementation of all 10 steps with dataclass-based configuration. Best for understanding the full algorithm or running on a single machine.

```bash
python code/gluon_pdf_full_workflow.py                                    # default params
python code/gluon_pdf_full_workflow.py --Pz 6 --conf 20000 --delta_z 15  # custom
```

Uses `argparse`. All logic is in one file: `LatticeConfig` + `CalcParams` dataclasses, then standalone functions for each step (plaquette, field strength, OPE operator, distillation, Fourier transform, matching, jackknife).

### `gluon_pdf_workflow.py` — Subcommand-based with MPI

A modular implementation designed for production HPC runs. Provides subcommands, built-in ensemble presets, MPI support via `mpi4py`, and backward-compatible stdin parameter parsing.

```bash
# Subcommands: 2pt, ope, pla, full
python code/gluon_pdf_workflow.py 2pt --ensemble L32x64 --conf-id 20000
python code/gluon_pdf_workflow.py ope --ensemble L32x64 --conf-id 20000 --gauge-file config.dat
python code/gluon_pdf_workflow.py pla --ensemble L32x64 --conf-id 20000 --gauge-file config.dat
python code/gluon_pdf_workflow.py full --ensemble L32x64 --conf-id 20000 --gauge-file config.dat

# Backward-compatible: read params from file (same format as examples/)
python code/gluon_pdf_workflow.py 2pt --params-file params.txt

# Multi-GPU via MPI
mpirun -np 4 python code/gluon_pdf_workflow.py ope --ensemble L32x64 --conf-id 20000 --gauge-file config.dat
```

**Built-in ensemble presets** (in `ENSEMBLES` dict, lines 107–151):

| Preset | β | L³×T | Nev | mom_smear | Notes |
|--------|---|------|-----|-----------|-------|
| `L24x72` | 6.20 | 24³×72 | 100 | -2 | Coarse, a≈0.105 fm |
| `L32x64` | 6.20 | 32³×64 | 100 | 3 | Coarse, a≈0.105 fm |
| `L32x96` | 6.41 | 32³×96 | 100 | 2 | Fine, a≈0.083 fm |
| `L36x108` | 6.498 | 36³×108 | 200 | 2 | Finer, a≈0.071 fm |
| `L48x96` | 6.20 | 48³×96 | 200 | 4 | Coarse large-volume |
| `L48x144` | 6.72 | 48³×144 | 200 | 2 | Physical pion mass |

Each preset hardcodes cluster paths for eigenvectors, perambulators, and correlator output directories (with `{conf_id}` substitution). Override individual paths with CLI arguments.

**When to use which**: Use `gluon_pdf_workflow.py` for production runs on the cluster (it handles MPI, file I/O, and ensemble presets). Use `gluon_pdf_full_workflow.py` for understanding the algorithm, testing on small lattices, or single-machine runs without MPI.

## GPU Code in `examples/`

The `examples/` directory contains 64 scripts for proton 2pt correlator computation. They follow a strict naming convention:

```
2pt_proton_Cg5gmu_{L}x{T}_mom{m}_{dir}_{backend}.py
```

- **{L}x{T}**: lattice size, e.g. `L32x64`, `L24x72`, `L32x96`, `L36x108`, `L48x144`
- **mom{m}**: momentum in 2π/L units, e.g. `mom0` (rest), `mom2`
- **{dir}**: momentum direction: `xdir`, `ydir`, `zdir`
- **{backend}**: `gpu` = NVIDIA CUDA (CuPy), `dcu` = AMD/Hygon DCU (ROCm/HIP via a PyTorch compatibility layer)

**Parameter passing**: Examples use **stdin redirection** with key-value pairs (parsed by `fileinput.input()`), NOT argparse:

```bash
python examples/2pt_proton_Cg5gmu_L32x64_mom2_xdir_gpu.py < params.txt
```

Example params.txt format:
```
Nt 64
Nx 32
Nev 100
Px 0
Py 0
Pz 6
```

**Key utility modules**:
- `Operator.py` — plaquette construction, field strength tensor F_{μν} (CuPy)
- `gamma_matrix_cupy_DR.py` — gamma matrices in DeGrand-Rossi (chiral-variant) basis
- `input_output_4_cupy.py` — I/O utilities for gauge configs and correlators

**Data paths**: The examples contain symlinks to cluster storage (`/public/group/lqcd/donghx/`). These will not resolve on a local machine — the actual data lives on the HPC cluster filesystem.

**Tensor conventions in examples**: Gauge field shape `[Nt, Nz, Ny, Nx, dir, color, color]` for OPE code (differs from the `code/` convention of `[color, color, dir, x, y, z, t]`).

## LaTeX Compilation

All LaTeX documents use **XeLaTeX** with Chinese support (ctex package). Key packages: `ctex`, `xeCJK`, `fontspec`, `physics`, `braket`, `beamer` (with metropolis/Madrid/CambridgeUS themes).

**Compile a document**:
```bash
# Single document (two-pass for cross-references)
cd 文档 && xelatex -interaction=nonstopmode -halt-on-error gluon_pdf_derivation.tex
xelatex -interaction=nonstopmode -halt-on-error gluon_pdf_derivation.tex

# With latexmk (auto-reruns)
latexmk -xelatex -interaction=nonstopmode -halt-on-error <file>.tex

# Compile from any directory with output in source dir
xelatex -interaction=nonstopmode -output-directory=/root/lattice-pdf/补充 \
  /root/lattice-pdf/补充/格点QCD蒸馏方法解析.tex

# Clean auxiliary files
latexmk -c <file>.tex
```

**Use the `/compile-latex` skill** for batch operations:
```
/compile-latex 文档/gluon_pdf_derivation.tex   # compile one file
/compile-latex 补充                              # compile all 34 supplement notes
/compile-latex all                                # compile all 46 .tex files
/compile-latex status                             # check PDF status
/compile-latex clean                              # clean aux files
```

(Full details in `.claude/skills/compile-latex.md`)

## Project Skills

Skills are slash commands defined in `.claude/skills/`. Root-level skills (work anywhere in the repo):

| Skill | Purpose |
|-------|---------|
| `compile-latex` | Compile any LaTeX document — single file, whole directory, or all 46 .tex files |
| `run-gluon-pdf` | Run `gluon_pdf_full_workflow.py` with various options (all, custom, check) |
| `update-agent-docs` | Recompile three agent submodule PDFs (LQCD_Master, lamet-agent, PyQUDA) |

Agent-level skills in `agent/.claude/skills/`:
| Skill | Purpose |
|-------|---------|
| `generate-and-update` | Unified compile of all 46 LaTeX docs + run Python workflow |
| `update-deps-docs` | Update three dependency-library PDFs |

## AI Agent Submodules (`agent/`)

The `agent/` directory contains three git submodules for AI-driven lattice QCD research:

| Submodule | Path | Description |
|-----------|------|-------------|
| **LQCD_Master** | `agent/LQCD_Master` | Translates natural-language physics requests into executable PyQUDA programs and Slurm submission scripts. Planner→Executor pipeline with human-in-the-loop checkpoints. |
| **lamet-agent** | `agent/lamet-agent` | Python-first LaMET workflow framework. Runs the full 5-stage analysis pipeline (correlator analysis → renormalization → Fourier transform → perturbative matching → extrapolation) via LLM-driven agentic stages. Has its own detailed `AGENTS.md`. |
| **PyQUDA** | `agent/PyQUDA` | Python wrapper for QUDA (HMC, gauge/fermion smearing, gradient flow, multigrid). |

```bash
# Clone with submodules
git clone --recurse-submodules <repo-url>

# After cloning without --recurse-submodules
git submodule update --init --recursive

# Update all submodules to latest remote
cd agent && ./update_deps.sh

# Update submodule documentation (compile PDFs)
cd agent && ./update_docs.sh             # compile all three docs
cd agent && ./update_docs.sh status       # check PDF status
cd agent && ./update_docs.sh deps         # update submodules + compile docs
```

## Important Project Conventions

- **Chinese is the primary language** for documentation, notes, and comments. Code docstrings and variable names are in English.
- **Python 3** is the standard. This project has **no `env.sh`** (unlike sibling projects PyQCU/PyQUDA) — dependencies are managed via the system Python or conda.
- **Data files** (HDF5, npy, npz, h5) live on the HPC cluster at `/public/group/lqcd/` and `/public/home/sunp/`, not in this repo. They are gitignored. The `examples/` directory has symlinks to cluster paths.
- **Two parameter-passing modes coexist**: `code/` workflow files use `argparse` (modern); `examples/` scripts use stdin redirection (original, backward-compatible via `fileinput.input()`).
- **GPU backends**: `_gpu.py` = NVIDIA CUDA/CuPy; `_dcu.py` = AMD/Hygon DCU via ROCm/HIP. The DCU scripts use a PyTorch compatibility layer rather than raw CuPy.
- **CuPy/numpy fallback pattern**: `try: import cupy as cp; HAS_CUPY = True` / `except ImportError: cp = np; HAS_CUPY = False`
- **opt_einsum/numpy fallback pattern**: `try: from opt_einsum import contract` / `except ImportError: contract = np.einsum`
- **Tensor conventions**: gauge fields `[color, color, direction, x, y, z, t]` (code/) or `[t, z, y, x, dir, color, color]` (examples/OPE); fermion fields `[spin, color, x, y, z, t]`; parity-split tensors prepend a `[2]` dimension (even/odd sites).
- **Gamma matrices**: DeGrand-Rossi (DR, chiral-variant) basis, defined in `examples/gamma_matrix_cupy_DR.py`.
- **Placeholder files named `PASS`** in `code/` and `log/` mark directories that need to be tracked by git but were otherwise empty.
- **The `.gitignore`** excludes: build artifacts, LaTeX auxiliary files, compiled binaries, HDF5/npy data, Python cache, and editor config.
