# Mathematica 代码: 格点QCD胶子TMD-PDF推导

本目录包含使用梯度流重整化方案的连续极限下核子非极化胶子TMD部分子分布函数的格点QCD计算的完整 Mathematica/Wolfram Language 推导与数值实现代码。

## 文件说明

| 文件 | 内容 | 行数 |
|------|------|------|
| `gluon_tmd_gradient_flow_continuum.wl` | **主文件**: 完整理论推导 (15个部分), 涵盖梯度流方程 → TMD因子化 → 连续极限 → 数值工作流 | ~1500 |
| `symbolic_derivations.wl` | **符号推导**: 热核展开, BCH展开, BMR张量验证, SFTX RG解, 硬核推导, CS核RG方程 | ~600 |
| `numerical_computations.wl` | **数值计算**: 梯度流平滑演示, a²连续极限外推, b_⊥/k_⊥分布, CS核提取, 匹配链, 误差分析 | ~600 |
| `README.md` | 本文件 | — |

## 使用方法

在 Mathematica 或 Wolfram Engine 中运行:

```mathematica
(* 加载并运行主推导文件 *)
<< "/root/lattice-pdf/mathematic/gluon_tmd_gradient_flow_continuum.wl"

(* 加载符号推导 *)
<< "/root/lattice-pdf/mathematic/symbolic_derivations.wl"

(* 加载数值计算 *)
<< "/root/lattice-pdf/mathematic/numerical_computations.wl"
```

或在命令行:
```bash
wolframscript -f gluon_tmd_gradient_flow_continuum.wl
```

## 理论结构 (主文件 15 个部分)

### 第一部分: 梯度流 (Gradient Flow) 方程与求解
- 规范场和夸克场的流方程 (Lüscher 2010, 2013)
- 热核与平滑半径 √(8t)
- UV有限性定理 (Lüscher-Weisz 2011)
- 格点 Wilson Flow (三阶Runge-Kutta)

### 第二部分: 场强张量 F_{μν} 的 Clover 构造
- Plaquette BCH展开
- Clover项 Q_{μν} 定义与对称性
- Minkowski ↔ Euclidean 转换

### 第三部分: 胶子 TMD 算符的构造
- BMR 六个不变振幅分解 (Balitsky-Morris-Radyushkin 2020)
- 非极化和极化胶子流算符
- 可乘性重整化分类 (Li-Ma-Qiu 2019, Zhang et al. 2019)

### 第四部分: Staple型Wilson线与准TMD-PDF算符
- 三段规范平行传输构造
- 四类UV发散 (线性, cusp, pinch-pole, 对数)

### 第五部分: 梯度流重整化方案
- Monahan-Orginos (2017) 梯度流正规化
- 小流时展开 (SFTX)
- 环标记夸克场 (Makino-Suzuki 2014)
- PDF Moment 比值 (Francis et al. 2025)

### 第六部分: Wilson圈减除 + SDR 重整化
- 五格点间距系统验证 (Zhang et al. 2022)
- RI/MOM方案失败分析
- 在壳夸克矩阵元单圈MS̄结果

### 第七部分: TMD 软函数与 Collins-Soper 核
- 内禀软函数 S_I: 形状因子因子化提取
- CS核 K: 动量比提取 (NLO改进)
- MILC/CLS 双系综交叉验证

### 第八部分: 连续极限外推
- a² 标度线性/二次外推
- 联合连续-无穷大动量外推
- CLQCD系综参数 (C24P29, C32P23, C48P21)

### 第九部分: TMD 因子化与匹配
- 完整因子化定理 (Izubuchi et al. 2018)
- NLO 硬匹配核 (单圈)
- 重整化群重求和 (RGR)
- 四圈 β 函数与强耦合跑动

### 第十部分: 傅里叶变换与完整数值流程
- z → x Fourier变换 (纵向)
- b_⊥ → k_⊥ Hankel变换 (横向)
- 大 λ = zP^z 外推
- 完整10步工作流 (参考 code/gluon_pdf_full_workflow.py)

### 第十一部分: 数值示例与参数化
- 唯象模型 TMD-PDF
- 典型格点参数 (P^z=2.15 GeV, μ=2 GeV)

### 第十二部分: Jackknife 误差分析
- Jackknife重采样
- 系统误差传播
- 激发态污染双态拟合

### 第十三部分: 积分检验
- TMD → 共线PDF k_⊥积分
- 一致性检验

### 第十四部分: 胶子TMD vs 夸克TMD 对比
- 六项关键差异 (算符, 发散, 匹配, 连通性, 梯度流, 小x)

### 第十五部分: EIC唯象学展望
- TMD普适性检验
- Sivers符号反转
- EIC与格点QCD协同

## 主要参考源 (本库路径)

### 核心文档
- `要求.md` — 项目规格说明
- `docs/Note of gluon PDFs.pdf` — 胶子PDF内部笔记 (Eq 20, Eq 25)
- `文档/gluon_pdf_derivation.tex` — 胶子PDF推导
- `文档/gluon_PDF_continuum.tex` — 胶子PDF连续极限

### 补充材料 (LaTeX)
- `补充/格点QCD中的TMD_PDF.tex` — TMD-PDF完整理论
- `补充/格点QCD中的梯度流重整化.tex` — 梯度流完整推导
- `补充/格点QCD中的胶子算符.tex` — 胶子算符分类
- `补充/格点QCD中的场强张量.tex` — Clover项构造
- `补充/格点QCD中的重整化.tex` — 重整化方案
- `补充/格点QCD中的大动量有效理论.tex` — LaMET框架
- `补充/格点QCD中的外推.tex` — 连续极限与手征外推

### 代码
- `code/gluon_pdf_full_workflow.py` — 完整10步计算工作流 (~1900行)
- `examples/Calc_ope_unpol.py` — OPE部分GPU代码
- `examples/2pt_proton_Cg5gmu_L32x64_mom2_xdir_gpu.py` — 质子2pt代码
- `examples/Operator.py` — F_{μν}构造与Wilson线
- `examples/gamma_matrix_cupy_DR.py` — γ矩阵 (DeGrand-Rossi基)

### 外部文献 (在 docs/ 中)
- `Unpolarized gluon PDF of the nucleon from lattice QCD in the continuum limit.pdf` (Chen et al. 2025)
- `Renormalization of TMD parton distribution on the lattice.pdf` (Zhang et al. 2022)
- `Unpolarized TMD Parton Distributions of the Nucleon from Lattice QCD.pdf` (He et al. 2024)
- `Quasi parton distributions and the gradient flow.pdf` (Monahan & Orginos 2017)
- `Accessing gluon parton distributions in large momentum effective theory.pdf` (Zhang et al. 2019)
- `Multiplicative renormalizability of quasi-parton operators.pdf` (Li, Ma, Qiu 2019)
- `Lattice-QCD Calculations of TMD Soft Function Through LaMET.pdf` (LPC 2020)
- `Nonperturbative Determination of CS Kernel from Quasi TMDWFs.pdf` (Chu et al. 2022)
- `Lattice Calculation of the Intrinsic Soft Function and the CS Kernel.pdf` (Chu et al. 2023)
- `Gradient Flow for Parton Distribution Functions First Application to the Pion.pdf` (Francis et al. 2025)

## 关键公式索引

| 公式 | 位置 | 说明 |
|------|------|------|
| ∂_t B_μ = D_ν G_{νμ} | 第一部分 | 梯度流方程 |
| F_{μν} = -i/(8a²)(Q_{μν} - Q_{νμ}) | 第二部分 | Clover场强张量 |
| M = Σ M_{pp,zz,zp,pz,ppzz,gg} | 第三部分 | BMR不变振幅分解 |
| O(z) = M_{tx;tx} + M_{ty;ty} - 2M_{xy;xy} | 第三部分 | 非极化胶子算符 |
| h̃(b_⊥,z,L,P^z) = <N\|FWF\|N> | 第四部分 | 准TMD裸矩阵元 |
| h = h̃ / √Z_E(2L+z,b_⊥) | 第六部分 | Wilson圈减除 |
| S_I = 2N_c F/\|φ(0)\|² | 第七部分 | 内禀软函数提取 |
| K = ln\|φ₁/φ₂\| / ln(P₁/P₂) | 第七部分 | CS核提取 |
| h(a) = h(0) + c₁a² | 第八部分 | 连续极限外推 |
| f̃/√S_I = H·Exp[½ln(ζ_z/ζ)K]·f | 第九部分 | TMD因子化定理 |
| H = 1 + (α_sC_F/2π)[-2+π²/12+ln-½ln²] | 第九部分 | NLO硬匹配核 |
