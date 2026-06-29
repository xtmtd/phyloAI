# PhyloAI 安装指南

[English](installation.md) | [中文](installation.zh.md)

## 目的

本指南说明如何安装 PhyloAI、如何让外部工具在当前 shell 环境中可见，以及如何通过 `phyloai doctor` 验证安装。

PhyloAI 不会自动安装第三方系统发育工具。请通过操作系统包管理器、Conda/Mamba 环境、集群模块系统或上游项目说明自行安装。

## 获取源代码

```bash
git clone https://github.com/xtmtd/phyloAI.git
cd phyloAI
```

## Python 环境

选择其中一种环境管理方式。

### uv

推荐用于本地开发与快速复现。

```bash
uv venv
source .venv/bin/activate
uv pip install -e .
```

### Conda / Mamba

推荐在 Python 包与生信命令行工具需要共存于同一环境时使用。

```bash
mamba create -n phyloai python=3.11
mamba activate phyloai
pip install -e .
```

如果系统只提供 `conda`，则将上面命令中的 `mamba` 替换为 `conda`。

### venv

适用于纯 Python 环境，前提是外部工具已通过其他方式加入 `PATH`。

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .
```

## 验证环境

```bash
phyloai doctor
phyloai doctor --output-format json
```

`phyloai doctor` 会检查当前 shell 环境。运行前请先激活你打算使用的环境。它不会安装工具，也不会修改 `PATH`。

## 一键流水线依赖

| 命令模式 | 主要外部工具 |
|----------|--------------|
| `phyloai run --mode supermatrix --speed normal` | MAFFT、trimAl、Julia（TAPER 用）、IQ-TREE3 |
| `phyloai run --mode supermatrix --speed fast` | MAFFT、trimAl、FastTree |
| `phyloai run --mode supertree --speed normal` | MAFFT、trimAl、Julia（TAPER 用）、IQ-TREE3（基因树）、wASTRAL |
| `phyloai run --mode supertree --speed fast` | MAFFT、trimAl、FastTree、wASTRAL |

## 外部工具

### 核心工作流

| 工具 | 用途 | 安装入口 | 检测名称 | 验证方式 |
|------|------|----------|----------|----------|
| IQ-TREE3 | `tree ml iqtree`、拓扑检验、定年 Hessian、normal 模式 supermatrix `run` | https://github.com/iqtree/iqtree3/releases | `iqtree3` | `phyloai doctor` |
| MAFFT | `pretree align`、`phyloai run` | https://mafft.cbrc.jp/alignment/software/ | `mafft` | `phyloai doctor` |
| trimAl | `pretree trim`、回译、`phyloai run` | https://github.com/inab/trimal | `trimal` | `phyloai doctor` |

请将这些视为外部工具：PhyloAI 仅检查它们是否对当前 shell 可见，不会重新分发它们。

### 树推断与后树分析工具

| 工具 | 用途 | 安装入口 | 检测名称 | 验证方式 |
|------|------|----------|----------|----------|
| FastTree | `tree ml fasttree`、fast 模式 `phyloai run` | http://www.microbesonline.org/fasttree/ | `FastTree` | `phyloai doctor` |
| wASTRAL | `tree msc`、supertree 模式 `phyloai run` | https://github.com/chaoszhang/ASTER | `wastral` | `phyloai doctor` |
| MCMCtree / PAML | `posttree dating mcmc` | https://github.com/abacus-gene/paml/releases | `mcmctree` | `phyloai doctor` |

### 贝叶斯推断

`phyloai tree bi` 需要 PhyloBayes-MPI 工具组。

| 工具 | 用途 | 安装入口 | 检测名称 | 验证方式 |
|------|------|----------|----------|----------|
| pb_mpi | MCMC 采样器 | https://github.com/bayesiancook/pbmpi | `pb_mpi` | `phyloai doctor` |
| bpcomp | 拓扑收敛诊断 | https://github.com/bayesiancook/pbmpi | `bpcomp` | `phyloai doctor` |
| tracecomp | 参数收敛诊断 | https://github.com/bayesiancook/pbmpi | `tracecomp` | `phyloai doctor` |
| mpirun | MPI 启动器 | https://www.open-mpi.org/ | `mpirun` | `phyloai doctor` |
| readpb_mpi | 可选的链文件读取器 | https://github.com/bayesiancook/pbmpi | `readpb_mpi` | `phyloai doctor` |

如果工具安装在 `PATH` 之外，使用 `phyloai tree bi --pb-path /path/to/pbmpi/bin` 指定目录。

### 过滤与修剪扩展工具

| 工具 | 用途 | 安装入口 | 检测名称 | 验证方式 |
|------|------|----------|----------|----------|
| TreeShrink | `pretree filter treeshrink` | https://github.com/uym2/TreeShrink | `run_treeshrink.py` | `phyloai doctor` |
| MAGUS | `pretree align --method magus` | https://github.com/vlasmirnov/MAGUS | `magus` | `phyloai doctor` |
| ClipKIT | `pretree trim --tool clipkit` | https://github.com/JLSteenwyk/ClipKIT | `clipkit` | `phyloai doctor` |
| BMGE | `pretree trim --tool bmge` | https://github.com/BMGE/BMGE 或上游 BMGE 发行版 | `BMGE.jar` | `phyloai doctor` |
| TAPER | `pretree filter taper` | 上游 TAPER 发行版 | `correction_multi.jl` | `phyloai doctor` |

### 运行时依赖

| 工具 | 用途 | 安装入口 | 检测名称 | 验证方式 |
|------|------|----------|----------|----------|
| Java | BMGE 工作流 | https://www.java.com/ | `java` | `phyloai doctor` |
| Julia | TAPER masking | https://julialang.org/downloads/ | `julia` | `phyloai doctor` |

如果 BMGE 或 TAPER 安装在 `PATH` 之外，通过 `--bmge-path /path/to/BMGE.jar` 或 `--taper-path /path/to/correction_multi.jl` 显式传入文件路径。

## 操作系统说明

- macOS：通常 Homebrew 或 Conda/Mamba 是让命令行工具进入 `PATH` 的最简方式。
- Linux：根据环境使用 Conda/Mamba、系统发行版包、集群模块或上游二进制。
- WSL：在 Linux 发行版内部安装工具（不要只装在 Windows 主机上），这样 `phyloai doctor` 才能在 WSL shell 中检测到。