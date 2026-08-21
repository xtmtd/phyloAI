# phyloai tree ml iqtree

[English](tree-ml-iqtree.md) | [中文](tree-ml-iqtree.zh.md)

## 目的

使用 IQ-TREE3 推断最大似然系统发育树。支持同质模型（带或不带 ModelFinder 与分区）、异质 AA 模型（混合与 PMSF）、NT 异质模型（MIX+MF）以及分支支持计算（UFBoot、SH-aLRT）。

## 用法

```bash
# 批量基因树（仅同质工作流）
phyloai tree ml iqtree --msa-dir <dir> [OPTIONS]

# 单一超矩阵（所有工作流）
phyloai tree ml iqtree --matrix <file> [OPTIONS]
```

## 示例

```bash
# 批量：20 个基因树，AA，LG 模型，4 个并行作业（默认）
phyloai tree ml iqtree --msa-dir msas/ --seq-type AA

# 单一矩阵：固定模型，UFBoot + SH-aLRT
phyloai tree ml iqtree --matrix matrix.fa --model LG --boot 1000 --alrt 1000

# 仅 ModelFinder（模型选择，不建树）
phyloai tree ml iqtree --matrix matrix.fa --modelfinder MF --mset LG,WAG

# ModelFinder + 建树 + 分区 + 合并
phyloai tree ml iqtree --matrix matrix.fa --modelfinder MFP --partitions parts.nex

# 禁用分支支持
phyloai tree ml iqtree --matrix matrix.fa --model LG --boot 0

# AA 混合模型（直接使用）
phyloai tree ml iqtree --matrix matrix.fa --model C20

# PMSF AA 混合
phyloai tree ml iqtree --matrix matrix.fa --model C20 --guide-tree guide.nwk

# 自定义交换率矩阵 + 逐位点频率 profile（CAT-PMSF 风格）
phyloai tree ml iqtree \
  --matrix runs/test/matrix.fa \
  --seq-type AA \
  --model runs/test/chain1.exchangeabilities \
  --site-freq-file runs/test/chain1.sitefreq \
  --state-freq none --rate-heterogeneity +R4 \
  --boot 0 --threads 1 --output-dir runs/test/cat-pmsf

# NT 异质
phyloai tree ml iqtree --matrix matrix.fa --seq-type NT --model MIX+MF
```

## 输入

必须提供 `--msa-dir` 或 `--matrix` 之一。批量模式用于同质工作流，单矩阵模式支持完整 IQ-TREE 工作流。

## 参数

### Input (mutually exclusive)

| Flag | Type | Description |
|------|------|-------------|
| `--msa-dir` | Path | 批量基因树推断的 MSA 文件目录 |
| `--matrix` | Path | 单一拼接矩阵，用于超矩阵推断 |

支持的格式：`.fa`、`.fas`、`.fasta`、`.faa`、`.fna`、`.phy`、`.phylip`、`.nex`、`.nxs`、`.nexus`、`.aln`。

### Data Type

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--seq-type` | `AA\|NT\|auto` | `auto` | 分子类型。`NT` → `--seqtype DNA` |

### Model (when `--modelfinder none`)

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--model` | str | `LG` (AA) / `GTR` (NT) | 替换模型，或已有的 IQ-TREE 自定义交换率模型文件 |
| `--state-freq` | `+F\|+FO\|+FQ\|+FU\|none` | `+F` | 状态频率类型 |
| `--rate-heterogeneity` | `+I\|+G4\|+I+G4\|+R4\|+I+R4\|none` | `+R4` | 速率异质性 |

它们组合形成 `-m`（例如 `LG+F+R4`）。当 `--modelfinder` 为 `MF` 或 `MFP` 时被忽略。

### ModelFinder

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--modelfinder` | `MF\|MFP\|none` | `none` | `MF` = 仅模型；`MFP` = 模型 + 建树 |
| `--mset` | str | `LG,WAG` / `GTR,HKY` | 模型搜索空间限制 |
| `--msub` | `nuclear\|mitochondrial\|chloroplast\|viral` | — | AA 模型来源（仅 AA） |

### Partitions (`--matrix` only)

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--partitions` | Path | — | 分区文件（NEXUS 或 RAxML 风格） |
| `--rclusterf` | int (1–100) | `10` | MF/MFP 的合并百分比 |
| `--rcluster-max` | int | — | 最大合并对数（与 `--rclusterf` 互斥） |

### Heterogeneous Models (`--matrix` only)

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--pmsf-base-model` | str | `LG` | 内置 C10–C60 PMSF 的基础 AA 模型 |
| `--guide-tree` | Path | — | 内置 PMSF 的引导树（NEWICK） |
| `--site-freq-file` | Path | — | 自定义 `--model` 的逐位点 AA 频率 profile，映射到 IQ-TREE `-fs`，必须配合 `--state-freq none` |
| `--qmax` | int | `10` | MIX+MF 的速率类别数 |

### Tree Search

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--mode` | `normal\|fast` | `normal` | `fast` → `--fast` |
| `--constraint` | Path | — | 约束树（NEWICK），映射到 IQ-TREE `-g` |

### Branch Support

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--boot` | int (>=0) | `1000` | UFBoot 重复数。`0` = 跳过 |
| `--alrt` | int (>=0) | — | SH-aLRT 重复数。`0` = 参数化 aLRT |
| `--bnni` | flag | `False` | 通过 NNI 优化 UFBoot 树 |

### Output

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--rate` | flag | `False` | 将位点速率写入 `.rate` 文件 |
| `--wslr` | flag | `False` | 将位点对数似然写入 `.sitelh` 文件 |
| `--outgroup` | str | — | 外类群分类单元（逗号分隔） |
| `--prefix` | str | — | 输出前缀（仅 `--matrix`） |

### Execution

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-o`, `--output-dir` | Path | `runs/tree/ml/iqtree` | 输出目录 |
| `--threads`, `-t` | int\|auto | `4` / `auto` | 批量：并行作业。单文件：IQ-TREE 线程 |
| `--overwrite` | flag | `False` | 先移除已有输出目录 |
| `--resume` | flag | `False` | 从 checkpoint 恢复 |
| `--dry-run` | flag | `False` | 打印命令而不执行 |
| `--keep-extra` | flag | `False` | 在 `logs/` 中保留额外 IQ-TREE 文件 |
| `-q`, `--quiet` | flag | `False` | 除错误外不打印输出 |
| `--iqtree-path` | Path | — | 自定义 `iqtree3` 路径 |
| `--tool-args` | str | — | 额外 IQ-TREE 标志 |

## 输出

### 单文件模式（`--matrix`）

```
runs/tree/ml/iqtree/
├── <prefix>.iqtree
├── <prefix>.treefile
├── <prefix>.log
├── ... (其他 IQ-TREE 输出)
├── result.json
```

### 批量模式（`--msa-dir`）

```
runs/tree/ml/iqtree/
├── trees/
│   ├── <gene1>.treefile
│   └── <gene2>.treefile
├── logs/
│   ├── <gene1>.iqtree
│   ├── <gene1>.log
│   └── ... （仅在使用 --keep-extra 时保留额外文件）
├── checkpoint.json
├── result.json
```

## 退出码

| Code | Meaning |
|------|---------|
| `0` | 成功 |
| `1` | 用户输入错误 |
| `2` | 所有 IQ-TREE 运行失败 |
| `3` | 未找到 `iqtree3` |

## 警告与错误

| 条件 | 行为 |
|------|------|
| `--modelfinder MF` 与 `--boot`/`--alrt`/`--bnni` | 警告：仅模型模式下忽略分支支持标志 |
| `--bnni` 不带 `--boot > 0` | 警告：`--bnni` 无效 |
| `--prefix` 在 `--msa-dir` 批量模式下 | 警告：忽略 `--prefix`；使用基因名 |
| `--rclusterf`/`--rcluster-max` 不带 `--partitions` | 警告：这些标志无效 |
| `--qmax` 不带 `--model MIX+MF` | 警告：`--qmax` 仅在 MIX+MF 时生效 |
| 异质模型与 `--msa-dir` | 错误：异质工作流要求 `--matrix` |
| `--partitions` 与 `--msa-dir` | 错误：`--partitions` 要求 `--matrix` |
| `--overwrite` 与 `--resume` 同时使用 | 错误：互斥 |
| 输出目录非空且未加 `--overwrite` | 错误：目录已存在 |
| `--site-freq-file` 未配合 AA 自定义 `--model`、`--matrix`、`--modelfinder none` 或 `--state-freq none` | 错误：无效的自定义 profile 组合 |

## 备注

- `--boot` 默认为 `1000`（启用），与 `phyloai tree ml fasttree` 一致。使用 `--boot 0` 跳过分支支持。
- `--threads` 默认在批量模式下为 `4` 个并行作业，单文件模式下为 `auto`（IQ-TREE 自行决定最优线程数）。
- 单 `--matrix` 模式将 IQ-TREE stdout 流式输出到终端以查看进度。
- 批量 `--msa-dir` 模式默认在 `logs/` 中仅保留 `.iqtree` 与 `.log` 文件。使用 `--keep-extra` 保留所有 IQ-TREE 输出文件。
- 当 `--modelfinder` 为 `MF` 或 `MFP` 时，`--model`、`--state-freq`、`--rate-heterogeneity` 不会传给 IQ-TREE，且在 `result.json` 中记为 `null`。
- IQ-TREE 输入及用户提供的路径参数（`--partitions`、`--guide-tree`、`--constraint`、自定义 `--model`、`--site-freq-file`）在内部解析为绝对路径。
- 自定义模型/profile 文件由 IQ-TREE 直接读取；PhyloAI 不复制或修改它们。`--tool-args "-fs /absolute/profile"` 会覆盖 `--site-freq-file`，仍作为原始 `tool_args` 记录，也必须使用 `--state-freq none`。
- `--matrix` 模式下的 `--resume` 会重新运行 IQ-TREE 命令；IQ-TREE 通过自己的机制（`--redo`）原生处理 checkpoint/resume。`--msa-dir` 批量模式下，PhyloAI 管理 checkpoint 状态以跳过已完成的基因树。
