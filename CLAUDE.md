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

**All other papers in `docs/`** cover: quasi-PDFs, pseudo-PDFs, renormalization, distillation, TMDs, gluon operators, Wilson lines, and related topics.

## Repository Structure

| Directory | Content |
|---|---|
| `code/` | **Main analysis code**: `gluon_pdf_full_workflow.py` (~1900 lines) — a systematic reconstruction of the full 10-step computation pipeline from gauge configurations to light-cone PDF |
| `examples/` | **Original GPU code** by donghx: `2pt_proton_Cg5gmu_L32x64_mom2_xdir_gpu.py` (proton 2pt), `Calc_ope_unpol.py` (OPE, MPI-parallel), with CuPy GPU utilities (`Operator.py`, `gamma_matrix_cupy_DR.py`, `input_output_4_cupy.py`). Contains symlinks to remote data on the cluster (`/public/group/lqcd/donghx/`) |
| `agent/` | **AI agent submodules** and documentation: three git submodules for agentic lattice QCD workflows |
| `refer/` | **Theory analysis & workflow document**: `理论解析与工作流.tex` — comprehensive LaTeX document with the theoretical derivation and computational workflow |
| `reports/` | **Beamer slides**: `gluon_pdf_slides.tex` (theoretical basis & workflow) and `gluon_pdf_continuum_beamer.tex` (continuum-limit results) |
| `补充/` | **Supplementary LaTeX notes** (~35 topics in Chinese) covering all aspects of lattice QCD: fermion actions, gauge fields, Wilson lines, correlation functions, renormalization, distillation, Monte Carlo methods, error analysis, Feynman diagrams, etc. |
| `文档/` | **Documentation**: gluon PDF derivation (`gluon_pdf_derivation.tex`), continuum extrapolation (`gluon_PDF_continuum.tex`), internal note (`note_of_gluon_PDFs.tex`) |
| `汇报/` | **Reports**: gluon quasi-operator construction on the lattice, lattice calculation methods |
| `代码/` | **Code analysis**: `code_analysis.tex` — analysis of the GPU computation code |
| `要求.md` | Project requirements/specification document (the starting point) |

## AI Agent Submodules (`agent/`)

The `agent/` directory contains three git submodules for AI-driven lattice QCD research:

| Submodule | Path | Description |
|-----------|------|-------------|
| **LQCD_Master** | `agent/LQCD_Master` | Translates natural-language physics requests into executable PyQUDA programs and Slurm submission scripts. Uses a Planner→Executor pipeline with human-in-the-loop checkpoints. Benchmarked at 90% accuracy on 70 validation tasks. |
| **lamet-agent** | `agent/lamet-agent` | Python-first LaMET workflow framework. Runs the full 5-stage analysis pipeline (correlator analysis → renormalization → Fourier transform → perturbative matching → extrapolation) via LLM-driven agentic stages. Reads JSON manifests, writes NetCDF intermediate artifacts. |
| **PyQUDA** | `agent/PyQUDA` | Python wrapper for QUDA (HMC, gauge/fermion smearing, gradient flow, multigrid). Required by LQCD_Master and lamet-agent. |

### Working with submodules

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
cd agent && ./update_docs.sh lamet-agent  # compile single doc
```

### Agent scripts

- `agent/update_deps.sh` — `git pull --ff-only` on all three submodules (clones if missing)
- `agent/update_docs.sh` — compiles PDF documentation from LaTeX sources in `agent/文档/` using xelatex (two-pass). Subcommands: `all`, `LQCD_Master`, `lamet-agent`, `PyQUDA`, `status`, `deps`
- `agent/文档/` — LaTeX source and compiled PDFs for each submodule's documentation

### lamet-agent quick reference

```bash
cd agent/lamet-agent
pip install -e ".[dev,analysis]"
lamet-agent validate examples/pion_pdf_cg_manifest.json
lamet-agent run examples/pion_pdf_cg_manifest.json --backend api --model provider/model_id --verbose
```

The lamet-agent has its own detailed `AGENTS.md` at `agent/lamet-agent/AGENTS.md` — consult it for stage architecture, manifest conventions, and development rules.

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

### Running the workflow

```bash
# Basic run with defaults
python code/gluon_pdf_full_workflow.py

# With parameters
python code/gluon_pdf_full_workflow.py \
  --conf 20000 --Pz 6 --delta_z 15 \
  --gauge_file gauge_config.dat \
  --eig_dir ./eigvecs/ --peram_dir ./perams/ \
  --output_dir ./results/
```

### Python Dependencies

```python
numpy, cupy        # GPU acceleration (falls back to numpy if cupy unavailable)
opt_einsum          # Optimized tensor contractions (falls back to numpy.einsum)
```

The code is designed to run on HPC clusters with GPUs. Conditional imports handle environments where CuPy/opt_einsum are unavailable.

## GPU Code in `examples/`

The original GPU code uses different conventions from the refactored `code/` version:

- **Parameter passing**: via stdin redirection with key-value pairs (parsed by `fileinput.input()`), not argparse
- **Multi-GPU (OPE)**: uses `mpi4py` for MPI distribution across processes
- **Data paths**: contain symlinks to cluster storage (`/public/group/lqcd/donghx/`)
- **DeGrand-Rossi basis**: gamma matrices use the DR (chiral-variant) basis, defined in `gamma_matrix_cupy_DR.py`
- **Tensor conventions**: gauge field shape `[Nt, Nz, Ny, Nx, dir, color, color]` for OPE code; field strength tensor uses plaquette construction in `Operator.py`

## LaTeX Compilation

All LaTeX documents use **XeLaTeX** with Chinese support (ctex package). Key packages: `ctex`, `xeCJK`, `fontspec`, `physics`, `braket`, `beamer` (with metropolis/Madrid/CambridgeUS themes).

**Compile a document**:
```bash
# Single document (from its directory)
cd 文档 && xelatex -interaction=nonstopmode -halt-on-error gluon_pdf_derivation.tex

# Two-pass for cross-references
xelatex -interaction=nonstopmode gluon_pdf_derivation.tex
xelatex -interaction=nonstopmode gluon_pdf_derivation.tex

# With latexmk (auto-reruns)
latexmk -xelatex -interaction=nonstopmode -halt-on-error <file>.tex

# Compile from any directory with output in source dir
xelatex -interaction=nonstopmode -output-directory=/root/lattice-pdf/补充 \
  /root/lattice-pdf/补充/格点QCD蒸馏方法解析.tex

# Clean auxiliary files
latexmk -c <file>.tex
```

**Compile beamer slides**:
```bash
cd reports
latexmk -xelatex -interaction=nonstopmode gluon_pdf_slides.tex
latexmk -xelatex -interaction=nonstopmode gluon_pdf_continuum_beamer.tex
```

**Check if required packages are available**:
```bash
kpsewhich <package>.sty   # e.g., kpsewhich ctex.sty
```

**Extract text from PDFs** for reading/analysis:
```bash
pdftotext "path/to/document.pdf" /tmp/output.txt
```

## Project Skills

Skills are slash commands for common project tasks, defined as `.claude/skills/<name>.md` files.
Skills are scoped to the directory they live in — root-level skills work anywhere in the repo;
agent-level skills only trigger when working under `agent/`.

### Root-level skills (`.claude/skills/`)

| Skill | Purpose |
|-------|---------|
| `compile-latex` | Compile any LaTeX document — single file, whole directory, or all 46 .tex files |
| `run-gluon-pdf` | Run the gluon PDF full workflow (`code/gluon_pdf_full_workflow.py`) |
| `update-agent-docs` | Recompile three agent submodule PDFs (LQCD_Master, lamet-agent, PyQUDA) |

```
/compile-latex 文档/gluon_pdf_derivation.tex   # compile one file
/compile-latex 补充                              # compile all 34 supplement notes
/compile-latex all                                # compile all 46 .tex files
/compile-latex status                             # check PDF status
/compile-latex clean                              # clean aux files

/run-gluon-pdf all                                # run full 10-step pipeline
/run-gluon-pdf check                              # check Python dependencies
/run-gluon-pdf custom --Pz 6 --conf 20000        # custom parameters

/update-agent-docs all                            # recompile all three agent PDFs
/update-agent-docs deps                           # pull source + compile
/update-agent-docs status                         # check agent PDF status
```

### Agent-level skills (`agent/.claude/skills/`)

These duplicate the three root skills above but are scoped to `agent/`:

| Skill | Root equivalent | Purpose |
|-------|----------------|---------|
| `generate-and-update` | `compile-latex` + `run-gluon-pdf` | 统一编译全部 46 个 LaTeX 文档和运行 Python 工作流 |
| `update-deps-docs` | `update-agent-docs` | 更新三个依赖库 PDF 文档 |

### Website skills (`../lattice-qcd-at-imp.top/.claude/skills/`)

The group website is a sibling project with its own skills:

| Skill | Purpose |
|-------|---------|
| `build-website` | 从零构建整个课题组网站（SPA, GitHub Pages） |
| `update-website` | 更新网站内容（论文/会议/讲习班/学生/导师/翻译/主题） |

## Important Project Conventions

- **Chinese is the primary language** for documentation, notes, and comments. Code docstrings and variable names are in English.
- **Python 3** is the standard. The `.vscode/settings.json` configures `conda` as the package/environment manager.
- **Data files** (HDF5, npy, npz, h5) live on the HPC cluster, not in this repo — they are gitignored. The `examples/` directory has symlinks to cluster paths.
- **GPU code** in `examples/` is the original (donghx's); `code/gluon_pdf_full_workflow.py` is a systematic refactored version.
- **The `.gitignore`** excludes: build artifacts, LaTeX auxiliary files (`.aux`, `.log`, `.toc`, `.nav`, `.snm`, `.xdv`, `.fdb_latexmk`, `.fls`), compiled binaries, HDF5/npy data, Python cache, and editor config.
- **Placeholder files named `PASS`** in `code/` and `log/` mark directories that need to be tracked by git but were otherwise empty.
- **Tensor conventions**: gauge fields `[color, color, direction, x, y, z, t]` or `[t, z, y, x, dir, color, color]` depending on the code module; fermion fields `[spin, color, x, y, z, t]`; parity-split tensors prepend a `[2]` dimension (even/odd sites).
- **CuPy/numpy fallback pattern**: `try: import cupy as cp; HAS_CUPY = True` / `except ImportError: cp = np; HAS_CUPY = False`
- **opt_einsum/numpy fallback pattern**: `try: from opt_einsum import contract` / `except ImportError: contract = np.einsum`
