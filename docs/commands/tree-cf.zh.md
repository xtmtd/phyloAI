# phyloai tree cf

[English](tree-cf.md) | [中文](tree-cf.zh.md)

一致性因子计算 —— 使用基因一致性（gCF）、位点一致性（sCF、sCFl）与四分体一致性（qCF）量化分支支持。

## 目的

一致性因子衡量参考物种树中支持每个 bipartition 的基因树、位点或四分体的比例。它们提供超越标准 bootstrap 的互补分支支持信息。

## 用法

```bash
phyloai tree cf --cf MODE --ref-tree REF_TREE [INPUTS...] [OPTIONS]
```

## 模式

| `--cf`    | Index  | Tool    | Description | Origin |
|-----------|--------|---------|-------------|--------|
| `gcf`     | gCF    | IQ-TREE3 | 基因一致性因子 | Minh et al. 2020 |
| `scf`     | sCF    | IQ-TREE3 | 位点一致性因子（简约） | Minh et al. 2020 |
| `scfl`    | sCFl   | IQ-TREE3 | 位点一致性因子（似然） | Mo et al. 2023 |
| `gcf+scf` | gCF+sCF | IQ-TREE3 | 一次运行同时计算 gCF + sCF | Minh et al. 2020 |
| `qcf`     | qCF    | wASTRAL  | 四分体一致性因子 | Mirarab et al. 2014 |

### CF 索引详情

- **gCF**（基因一致性因子）：包含给定 bipartition 的基因树比例。评估基因树异质性。
- **sCF**（位点一致性因子，简约）：在简约标准下支持 bipartition 的比对位点比例。计算快。
- **sCFl**（位点一致性因子，似然）：使用最大似然支持 bipartition 的位点比例。比 sCF 更准确，但需要模型选择。
- **gCF+sCF**：在单次 IQ-TREE3 调用中同时计算 gCF 和 sCF，节省计算时间。
- **qCF**（四分体一致性因子）：支持每个 bipartition 的四分体（来自基因树）比例。由 wASTRAL 使用其四分体评分引擎计算。

## 各模式输入要求

| Mode     | `--ref-tree` | `--tree`/`--tree-dir` | `--matrix` | `--model-expr`/`--partitions` |
|----------|-------------|----------------------|-----------|-------------------------|
| `gcf`      | Required    | Required             | —         | —                       |
| `scf`      | Required    | —                    | Required  | —                       |
| `scfl`     | Required    | —                    | Required  | Optional (speedup)      |
| `gcf+scf`  | Required    | Required             | Required  | —                       |
| `qcf`      | Required    | Required             | —         | —                       |

## 选项

| Option | Default | Description |
|--------|---------|-------------|
| `--cf MODE` | *required* | 一致性因子类型（gcf/scf/scfl/gcf+scf/qcf） |
| `--ref-tree FILE` | *required* | 参考物种树（NEWICK） |
| `--tree FILE` | — | 单个基因树文件（与 `--tree-dir` 互斥） |
| `--tree-dir DIR` | — | 基因树文件目录 |
| `--matrix FILE` | — | 多序列比对（scf/scfl/gcf+scf 必需） |
| `--model-expr TEXT` | — | scfl 的替换模型（如 `LG+F+R4`） |
| `--partitions FILE` | — | scfl 的分区文件（如 `*.best_model.nex`） |
| `--scf-quartets N` | 100 | sCF/sCFl 的四分体数（推荐 ≥ 100） |
| `--prefix TEXT` | auto | 输出前缀（默认：gCF/sCF/sCFl/gCFsCF/qCF） |
| `-o, --output-dir DIR` | `runs/tree/cf` | 输出目录 |
| `-t, --threads N` | 4 | CPU 线程数 |
| `--iqtree-path PATH` | auto | 显式 iqtree3 可执行路径 |
| `--wastral-path PATH` | auto | 显式 wastral 可执行路径 |
| `--lpp` | off | 同时将局部后验概率附加到 qCF 标签 |
| `--overwrite` | off | 移除已有输出目录 |
| `--dry-run` | off | 显示命令而不执行 |
| `-q, --quiet` | off | 抑制非错误输出 |

## 示例

```bash
# gCF：基因树 + 参考树
phyloai tree cf --cf gcf --ref-tree species.nwk --tree-dir ./genetrees/

# gCF 单文件
phyloai tree cf --cf gcf --ref-tree species.nwk --tree merged.trees

# sCF：比对 + 参考树（理想情况下是 gCF 注解的）
phyloai tree cf --cf scf --ref-tree gCF.cf.tree --matrix msa.fa

# sCFl（似然）带模型加速
phyloai tree cf --cf scfl --ref-tree gCF.cf.tree --matrix msa.fa --model-expr LG+F+R4

# sCFl 使用预计算分区模型
phyloai tree cf --cf scfl --ref-tree gCF.cf.tree --matrix msa.fa \
    --partitions msa.best_model.nex

# gCF + sCF 组合
phyloai tree cf --cf gcf+scf --ref-tree species.nwk --tree-dir ./genetrees/ \
    --matrix msa.fa

# 通过 wASTRAL 计算 qCF
phyloai tree cf --cf qcf --ref-tree species.nwk --tree merged.trees

# qCF 附加局部后验概率（pp1）
phyloai tree cf --cf qcf --ref-tree species.nwk --tree merged.trees --lpp

# 自定义输出前缀与线程
phyloai tree cf --cf gcf --ref-tree species.nwk --tree merged.trees \
    --prefix myCF -t 8
```

## 输出文件

### IQ-TREE3 模式（gcf、scf、scfl、gcf+scf）

| File | Description |
|------|-------------|
| `<prefix>.cf.stat`  | 一致性因子统计表 |
| `<prefix>.cf.branch` | 带分支 ID 的树 |
| `<prefix>.cf.tree`  | 注解了 CF 值的树 |
| `<prefix>.cf.tree.nex` | 供 FigTree 使用的 NEXUS 注解树 |
| `<prefix>.log`       | IQ-TREE3 原生日志（引用为 `data.tool_log`） |
| `result.json`        | PhyloAI 结构化结果（stderr 内联于 `data.tool_stderr`） |
| `merged.trees`       | 合并后的基因树（若使用 `--tree-dir`） |

### qCF 模式

| File | Description |
|------|-------------|
| `<prefix>.cf.tree` | 注解了 qCF（以及可选 pp1）的参考树 |
| `wastral.tre`      | 原始 wASTRAL 输出（中间产物） |
| `result.json`      | PhyloAI 结构化结果（stderr 内联于 `data.tool_stderr`） |
| `merged.trees`     | 合并后的基因树（若使用 `--tree-dir`） |

## qCF 输出格式

qCF 值在 [0,1] 范围内保留为原始小数（不乘 100）。为可读性去除尾随零（例如 `1` 而非 `1.0000`，`0.95` 而非 `0.9500`）。当附加到现有支持值时，格式为：

- 不使用 `--lpp`：`<support>/<q1>`（例如 `100/0.4221`）
- 使用 `--lpp`：`<support>/<q1>/<pp1>`（例如 `100/0.4221/0.95`）

如果不存在现有支持，qCF 值成为唯一标签：`0.75`。

## 警告与错误

- 输入树、比对和参考树的分类单元必须符合所选模式要求。
- 输出目录非空时必须使用 `--overwrite`。

## 备注

- 为获得最佳 sCF/sCFl 结果，请使用 gCF 注解的树作为 `--ref-tree`（例如先运行 `--cf gcf`）。
- `--cf scfl` 不带 `--model-expr` 或 `--partitions` 会自动计算最佳拟合模型 —— 这很慢。提供 `--model-expr` 或 `--partitions` 可以加速。
- `--scf-quartets` 应 ≥ 100 以获得可靠结果。更高的值在运行时成本上提高准确性。
- gCF+sCF 计算在一次 IQ-TREE3 调用中同时运行两种模式，相比两次单独运行节省大量时间。
- qCF 使用 wASTRAL 的校准四分体评分（`-u 2 -C --mode 4`），可处理基因树估计误差。

## 参考文献

- Minh BQ, Hahn MW, Lanfear R (2020) New methods to calculate concordance factors for phylogenomic datasets. *Molecular Biology and Evolution* **37**(5):1530–1534.
- Mo YK, Lanfear R, Hahn MW, Minh BQ (2023) Updated site concordance factors minimize effects of homoplasy and taxon sampling. *Systematic Biology* **72**(3):559–574.
- Mirarab S, Reaz R, Bayzid MS, Zimmermann T, Swenson MS, Warnow T (2014) ASTRAL: genome-scale coalescent-based species tree estimation. *Bioinformatics* **30**(17):i541–i548.
- Zhang C, Rabiee M, Sayyari E, Mirarab S (2018) ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees. *BMC Bioinformatics* **19**(Suppl 6):153.