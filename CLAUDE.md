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
| `examples/` | **All Python code**, organized by contributor: `donghx/` (proton 2pt + OPE), `zhangxin/` (gluon PDF workflow + data analysis), `huangcl/` (Chroma-based analysis) |
| `docs/` | **80+ reference papers** (PDF) on LaMET, quasi-PDFs, distillation, renormalization, gradient flow, etc. |
| `agent/` | **AI agent submodules**: LQCD_Master, lamet-agent, PyQUDA (git submodules; see `agent/CLAUDE.md` in each for details) |
| `refer/` | `理论解析与工作流.tex` — comprehensive LaTeX document with theoretical derivation and computational workflow |
| `reports/` | **Beamer slides**: `gluon_pdf_slides.tex` (theoretical basis & workflow) and `gluon_pdf_continuum_beamer.tex` (continuum-limit results) |
| `补充/` | **35 supplementary LaTeX notes** (Chinese) covering all aspects of lattice QCD: fermion actions, gauge fields, Wilson lines, correlation functions, renormalization, distillation, Monte Carlo methods, error analysis, etc. |
| `文档/` | **Core documents**: gluon PDF derivation (`gluon_pdf_derivation.tex`), continuum extrapolation (`gluon_PDF_continuum.tex`), internal note (`note_of_gluon_PDFs.tex`) |
| `汇报/` | **Reports**: gluon quasi-operator construction on the lattice, lattice calculation methods |
| `代码/` | **Code analysis**: `code_analysis.tex` — analysis of the GPU computation code |

## GPU Code in `examples/donghx/`

donghx's 64 scripts for proton 2pt correlator computation follow a strict naming convention:

```
2pt_proton_Cg5gmu_{L}x{T}_mom{m}_{dir}_{backend}.py
```

- **{L}x{T}**: lattice size, e.g. `L32x64`, `L24x72`, `L32x96`, `L36x108`, `L48x144`, `L48x96`
- **mom{m}**: momentum in 2π/L units, e.g. `mom0` (rest), `mom2`
- **{dir}**: momentum direction: `xdir`, `ydir`, `zdir`
- **{backend}**: `gpu` = NVIDIA CUDA (CuPy), `dcu` = AMD/Hygon DCU (ROCm/HIP via a PyTorch compatibility layer)
- Some L32x96 variants have no backend suffix (`.py` only) — these are CPU versions

**Parameter passing**: These scripts use **stdin redirection** with key-value pairs (parsed by `fileinput.input()`), NOT argparse:

```bash
python examples/donghx/2pt_proton_Cg5gmu_L32x64_mom2_xdir_gpu.py < params.txt
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

**Key utility modules** (in `examples/donghx/`):
- `Operator.py` — plaquette construction, field strength tensor F_{μν} (CuPy)
- `gamma_matrix_cupy_DR.py` — gamma matrices in DeGrand-Rossi (chiral-variant) basis
- `input_output_4_cupy.py` — I/O utilities for gauge configs and correlators

**OPE calculation scripts** (in `examples/donghx/`):
- `Calc_ope_unpol.py` — unpolarized gluon OPE operator calculation (CuPy, non-MPI)
- `Calc_ope_unpol_new.py` — updated version with MPI support
- `Calc_ope_gauge_fix_unpol.py` — gauge-fixed variant
- `Calc_ope_helicity.py` / `Calc_ope_helicity_new.py` — helicity (polarized) gluon OPE
- `Calc_ope_gauge_fix_helicity.py` — gauge-fixed helicity OPE
- `Calc_pla.py` — plaquette calculation
- `Calc_VVV.py` — VVV (distillation perambulator contraction) calculation
- `check_exist.py` / `check_exist_TMD.py` — data existence checkers

**Data paths**: The directory contains symlinks to cluster storage (`/public/group/lqcd/donghx/`). These will not resolve on a local machine — the actual data lives on the HPC cluster filesystem.

**Tensor conventions in donghx code**: Gauge field shape `[Nt, Nz, Ny, Nx, dir, color, color]` (differs from the zhangxin convention of `[color, color, dir, x, y, z, t]`).

## Zhangxin's Workflow & Analysis in `examples/zhangxin/`

### Two Gluon PDF Workflow Implementations

Two independent implementations compute the same physics but differ in architecture:

#### `gluon_pdf_full_workflow.py` — Self-contained pipeline

A single-script implementation (~1900 lines) of all 10 steps with dataclass-based configuration. Best for understanding the full algorithm or running on a single machine.

```bash
python examples/zhangxin/gluon_pdf_full_workflow.py                                    # default params
python examples/zhangxin/gluon_pdf_full_workflow.py --Pz 6 --conf 20000 --delta_z 15  # custom
```

Uses `argparse`. All logic is in one file: `LatticeConfig` + `CalcParams` dataclasses, then standalone functions for each step (plaquette, field strength, OPE operator, distillation, Fourier transform, matching, jackknife).

#### `gluon_pdf_workflow.py` — Subcommand-based with MPI

A modular implementation (~1600 lines) designed for production HPC runs. Provides subcommands, built-in ensemble presets, MPI support via `mpi4py`, and backward-compatible stdin parameter parsing.

```bash
# Subcommands: 2pt, ope, pla, full
python examples/zhangxin/gluon_pdf_workflow.py 2pt --ensemble L32x64 --conf-id 20000
python examples/zhangxin/gluon_pdf_workflow.py ope --ensemble L32x64 --conf-id 20000 --gauge-file config.dat
python examples/zhangxin/gluon_pdf_workflow.py pla --ensemble L32x64 --conf-id 20000 --gauge-file config.dat
python examples/zhangxin/gluon_pdf_workflow.py full --ensemble L32x64 --conf-id 20000 --gauge-file config.dat

# Backward-compatible: read params from file (same format as examples/donghx/)
python examples/zhangxin/gluon_pdf_workflow.py 2pt --params-file params.txt

# Multi-GPU via MPI
mpirun -np 4 python examples/zhangxin/gluon_pdf_workflow.py ope --ensemble L32x64 --conf-id 20000 --gauge-file config.dat
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

### Data Analysis Framework

Zhangxin's analysis framework for post-processing distillation and IOG correlator data:

- **`include.py`** — Core module (~900 lines). Provides the `data_analyse` class for reading, analyzing, and plotting distillation 2pt/3pt correlator data. Supports both raw `.dat` files and Chroma IOG binary files. Handles jackknife resampling, effective mass extraction, ratio analysis, and link (Wilson line) analysis.
- **`main.py`** — Pion 3pt ratio analysis: loads `_3pt.npy` and `_2pt.npy`, computes the 3pt/2pt ratio R(z) with jackknife errors, and plots `pion_3pt.png`.
- **`main-2pt.py`** — Pion 2pt effective mass analysis from IOG files. Computes and plots `cosh` effective masses.
- **`main-3pt.py`** — Pion 3pt data extraction from IOG files. Reads 3pt and 2pt correlators, saves them as `.npy` for further analysis.
- **`main_iog.py`** — Full IOG-based analysis pipeline: reads Chroma IOG files, extracts 2pt/3pt correlators, computes ratios, fits link data, and plots results.
- **`_main.py`** — Extended IOG analysis with multi-source data (sush distillation + Chroma IOG), link ratio fitting, and effective mass comparison.
- **`lsq_tools.py`** — Utility functions for reading `.dat` correlator files and least-squares fitting with `lsqfit` + `gvar`.

### IOG Reader (`iog_reader/`)

A C extension module for reading Chroma's IOG binary format:
- `iog_reader.py` — Python wrapper using `ctypes` to call the compiled `.so`
- `iog.so` — compiled shared library (gitignored, must be built on the cluster)
- Reads IOG files into pandas DataFrames with integer labels and real/imaginary components

### Ensemble Data (`beta6.20_mu-0.2770_ms-0.2400_L24x72/`)

Contains per-ensemble analysis subdirectories with IOG files, Slurm error logs, and Chroma XML parameter files.

## Huangcl's Code in `examples/huangcl/`

- **`code.py`** — Chroma-based analysis script for L24x72 ensemble. Reads Chroma 2pt/3pt correlators, computes effective masses and ratios, supports jackknife resampling. Uses matplotlib for plotting.
- **`submit.sh`** — Slurm submission script for running `code.py` on the cluster.
- **`result/`** — Output directory for analysis results.

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

## Project Skills (`.claude/skills/`)

Skills are slash commands defined in `.claude/skills/`. Root-level skills (work anywhere in the repo):

| Skill | Purpose |
|-------|---------|
| `compile-latex` | Compile any LaTeX document — single file, whole directory, or all 35+ .tex files |
| `run-gluon-pdf` | Run `gluon_pdf_full_workflow.py` with various options (all, custom, check) |
| `update-agent-docs` | Recompile three agent submodule PDFs (LQCD_Master, lamet-agent, PyQUDA) |

Agent-level skills in `agent/.claude/skills/`:
| Skill | Purpose |
|-------|---------|
| `generate-and-update` | Unified compile of all LaTeX docs + run Python workflow |
| `update-deps-docs` | Update three dependency-library PDFs |

## Unified Workflow Scripts (`snsc/`)

The `snsc/` directory contains unified, production-quality entry points for the gluon PDF pipeline, consolidating code from `examples/donghx/`, `examples/zhangxin/`, and `examples/huangcl/` into a single configurable CLI with three analysis modes.

```bash
# All three modes via main.py
python snsc/main.py --help

# Mode 1: PDF (full LaMET pipeline, 13 optional steps)
python snsc/main.py --analysis-type pdf --steps 1,2,3,6 --xp numpy --ensemble L32x64 \
    --conf-start 20000 --conf-num 1 --gauge-file /path/to/config.dat

# Mode 2: proton-2pt (distillation, VVV + Wick contraction, matches donghx DCU code)
python snsc/main.py --analysis-type proton-2pt --Nt 72 --Nx 24 \
    --Pz-list "-2,-3,-4,-5,-6" --conf-start 46000 --element _Cg5g4

# Mode 3: 2pt (effective mass from IOG data, matches main-2pt.py)
python snsc/main.py --analysis-type 2pt --Nt 72 --Nx 24 --alttc 0.1053 \
    --conf-start 10050 --conf-step 50 --conf-num 52

# Slurm submission (always submit from snsc/ directory)
cd snsc
sbatch sbatch.sh          # PDF full pipeline (CPU)
sbatch sbatch-2pt.sh      # Proton 2pt distillation (CPU)
```

Key features:
- `--xp numpy|cupy` for backend selection, `--dtype complex64|complex128`
- Per-step wall time + peak memory tracking via `Timer` context manager
- Automatic VVV caching (`.npy` in output dir — skip recomputation)
- Read/generate/compare paths for 2pt, 3pt, VVV, OPE data
- Reference data comparison with max/mean absolute+relative diff metrics + scatter/histogram plots
- Cosh effective mass extraction from 2pt with plateau estimation
- Ensemble presets for all 6 lattice geometries + L64x128
- Timestamped output: `snsc/output_YYYYMMDD_HHMMSS/{data/, plots/, run_config.json}`
- Default conda env: `source /public/home/zhangxin/miniconda3/etc/profile.d/conda.sh && conda activate zhangxin-snsc`

**Slurm**: Both `sbatch.sh` and `sbatch-2pt.sh` use `#SBATCH --chdir=/public/home/zhangxin/lattice-pdf/snsc` so the working directory is always correct regardless of where `sbatch` is invoked. Slurm stdout/err goes to `snsc/logs/*_%j.{out,err}`; application logs go to `output_*/run.log`.

**Dependencies** (proton-2pt mode): requires `opt_einsum`, `numpy`, `matplotlib`. The 2pt mode additionally requires the compiled `examples/zhangxin/iog_reader/iog.so` and `include.py` from `examples/zhangxin/`.

## AI Agent Submodules (`agent/`)

Three git submodules for AI-driven lattice QCD research:

| Submodule | Path | Description |
|-----------|------|-------------|
| **LQCD_Master** | `agent/LQCD_Master` | Natural language → PyQUDA programs + Slurm scripts. Planner→Executor pipeline with checkpoints. |
| **lamet-agent** | `agent/lamet-agent` | LaMET analysis framework. 5-stage agentic pipeline (correlator→renorm→Fourier→matching→extrapolation). Uses manifest JSON to define workflows, pydantic-validated. |
| **PyQUDA** | `agent/PyQUDA` | Python wrapper for QUDA (HMC, gauge/fermion smearing, gradient flow, multigrid). |

### LQCD_Master Usage

```bash
cd agent/LQCD_Master

# Interactive mode: describe the physics task in natural language
python run.py

# Non-interactive with task file
python run.py --task task.txt --non-interactive
python run.py --task "Compute proton 2pt with momentum Pz=3 on L32x64 ensemble"

# Test mode: generate artifacts without submitting jobs
python run.py --task "..." --test

# List available skills
python run.py --list-skills
```

**Architecture**: `run.py` → `WorkflowOrchestrator` → `PlannerAgent` (generates plan from NL task) → `ExecutorAgent` (generates PyQUDA code + Slurm scripts). Output goes to `runs/<timestamp>/`. Uses `core_architecture/` and `utils/` modules; relies on `prompt_toolkit` for interactive input and a `.env` file for LLM API keys.

### lamet-agent Usage

```bash
cd agent/lamet-agent

# Validate a manifest
lamet-agent validate path/to/manifest.json

# Interactive planning (review/repair a manifest before running)
lamet-agent plan manifest.json --backend api --model deepseek/deepseek-chat

# Run a workflow (staged agent loop)
lamet-agent run manifest.json --backend api --model deepseek/deepseek-chat \
    --max-tool-steps 40 --report-language en
```

**Manifest structure**: JSON with `metadata` (run_id, root_directory, stages), `inputs` (correlators/artifacts/kernels with data_paths, ensemble, momenta, volume), and `stages` (per-stage defaults + jobs). Schema enforced by pydantic in `manifest.py`. Correlator data is read as HDF5/NetCDF via `core/data.py`. Effective mass uses `core/plotting.py` (cosh meff, publication-style figures with Times New Roman + STIX math).

**Stages**: `correlator_analysis` (2pt ground-state fitting, meff plots), `renormalization` (hybrid self-renormalization), `fourier_transform`, `perturbative_matching`, `extrapolation`, `review`. Each stage has `functions.py` (computation), `prompts.py` (LLM prompts), `skills.py` (tool definitions).

### agent/snsc/ — LQCD_Master + lamet-agent 联合流水线

Bridge directory running the proton-2pt distillation pipeline via both agent frameworks. Equivalent functionality to `snsc/sbatch-2pt.sh` but with agent-driven workflow orchestration and publication-quality analysis.

```
agent/snsc/
├── run.py                              # 主入口 (5-phase orchestration)
├── manifest.json                       # lamet-agent correlator manifest
├── sbatch.sh                           # Slurm 提交脚本
├── skills/proton_2pt_distillation/     # LQCD_Master 技能定义
│   └── SKILL.md                        # 蒸馏工作流描述
├── artifacts/                          # HDF5 数据 (lamet-agent 兼容)
├── runs/                               # 时间戳输出目录
└── logs/                               # Slurm 输出日志
```

**5-Phase Pipeline:**

| Phase | Framework | Description |
|-------|-----------|-------------|
| 1-2 | LQCD_Master Planner→Executor | NL task → structured plan + code generation |
| 3 | snsc/main.py | Distillation computation (eigvecs→VVV→Wick→2pt) |
| 4 | HDF5 Bridge | .npy output → lamet-agent-compatible HDF5 |
| 5 | lamet-agent | Effective mass fitting + publication plots |

```bash
cd agent/snsc

# 完整流水线 (蒸馏 + lamet meff 分析)
python run.py --conf-id 46000 --Pz-list "-2,-3,-4,-5,-6"

# 仅蒸馏
python run.py --distillation-only --conf-id 46000

# 仅 lamet 分析 (需已有 HDF5)
python run.py --lamet-only --h5-path artifacts/proton_2pt.h5

# 含 LQCD_Master 规划 (需要 .env 配置 LLM API)
python run.py --task "Compute proton 2pt on L24x72 ensemble, Pz=-2..-6"

# Slurm 提交
sbatch sbatch.sh
sbatch --export=ALL,CONF_ID="46000" sbatch.sh
```

**Key architectural decisions:**
- Reuses `snsc/main.py` for distillation computation (not duplicating code)
- lamet-agent manifest uses `resample_mode: "jk"` (jackknife) for error analysis
- Default `--skip-plan` on Slurm to avoid LLM costs in production
- HDF5 Bridge converts (Nt,Nt) parity-projected 2pt → (Lt, n_cfg) complex format lamet-agent expects
- lamet-agent does correlator_analysis stage with `fitting_form: "NonBreit"` for multi-momentum independent fits

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

## Python Dependencies

This project has **no `env.sh`** (unlike sibling projects PyQCU/PyQUDA). Standard dependencies:

- **GPU**: `cupy` (NVIDIA CUDA) or PyTorch (AMD DCU via ROCm/HIP compatibility layer)
- **Core**: `numpy`, `scipy`, `opt_einsum` (optional, falls back to `numpy.einsum`)
- **Distillation/analysis**: `lsqfit`, `gvar`, `sympy`, `matplotlib`, `pandas`, `proplot`
- **I/O**: `h5py` (HDF5), `fileinput` (stdlib, for stdin parameter parsing)
- **MPI**: `mpi4py` (optional, for multi-GPU distributed runs)
- **Chroma IOG**: `iog_reader/iog.so` (compiled C extension, gitignored — build on cluster)

## Important Project Conventions

- **Chinese is the primary language** for documentation, notes, and comments. Code docstrings and variable names are in English.
- **Python 3** is the standard. Dependencies are managed via the system Python or conda.
- **Data files** (HDF5, npy, npz, h5) live on the HPC cluster at `/public/group/lqcd/`, `/public/home/sunp/`, and `/public/home/sush/`, not in this repo. They are gitignored. The `examples/donghx/` directory has symlinks to cluster paths.
- **Two parameter-passing modes coexist**: `examples/zhangxin/` workflow files use `argparse` (modern); `examples/donghx/` scripts use stdin redirection (original, via `fileinput.input()`).
- **GPU backends**: `_gpu.py` = NVIDIA CUDA/CuPy; `_dcu.py` = AMD/Hygon DCU via ROCm/HIP with PyTorch compatibility layer.
- **CuPy/numpy fallback pattern**: `try: import cupy as cp; HAS_CUPY = True` / `except ImportError: cp = np; HAS_CUPY = False`
- **opt_einsum/numpy fallback pattern**: `try: from opt_einsum import contract` / `except ImportError: contract = np.einsum`
- **Tensor conventions**:
  - zhangxin code: gauge fields `[color, color, direction, x, y, z, t]`; fermion fields `[spin, color, x, y, z, t]`
  - donghx code: gauge fields `[t, z, y, x, dir, color, color]`
  - Parity-split tensors prepend a `[2]` dimension (even/odd sites).
- **Gamma matrices**: DeGrand-Rossi (DR, chiral-variant) basis, defined in `examples/donghx/gamma_matrix_cupy_DR.py`.
- **Placeholder files named `PASS`** mark directories that need to be tracked by git but were otherwise empty.
- **The `.gitignore`** excludes: build artifacts, LaTeX auxiliary files, compiled binaries, HDF5/npy data, Python cache, IOG `.so` files, and editor config.
