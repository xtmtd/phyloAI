# phyloai posttree signal

[English](posttree-signal.md) | [中文](posttree-signal.zh.md)


## 目的

通过三个独立子命令进行系统发育信号分布分析：

| 子命令 | 分析内容 | 核心工具 |
|--------|----------|----------|
| `lnl` | 位点级和基因级对数似然分数分布 | IQ-TREE3 `-wslr` |
| `consistent` | 一致基因识别（GLS + GQS，Shen et al. 2021） | IQ-TREE3 + wASTRAL |
| `fclm` | 四簇似然映射 | IQ-TREE3 `-lmap -lmclust` |

这些分析用于检查系统发育信号如何在位点和基因间分布，识别对拓扑有不成比例影响的离群基因，比较支持不同拓扑的基因群之间的指标差异，并通过似然映射评估争议分支的系统发育信号。

## 用法

```bash
# 位点级 lnL 分布（同质模型）
phyloai posttree signal lnl --matrix ./matrix.fa --candidate-trees trees --model-expr LG+F+R4

# 分区模型 + 基因级输出
phyloai posttree signal lnl --matrix ./matrix.fa --candidate-trees trees --partitions partitions.txt

# 基因级 lnL（含基因座范围）+ 离群分析
phyloai posttree signal lnl --matrix ./matrix.fa --candidate-trees trees --model-expr LG+F+R4 --locus-ranges partitions.txt --metrics metrics.csv

# 一致基因识别（GLS + GQS）
phyloai posttree signal consistent --matrix ./matrix.fa --candidate-trees T1.tre,T2.tre --tree-dir ./gene_trees --model-expr LG+F+R4 --locus-ranges partitions.txt

# 四簇似然映射
phyloai posttree signal fclm --matrix ./matrix.fa --taxset-csv taxsets.csv --model-expr LG+C60+F+R4

# 四簇似然映射 + 分区模型
phyloai posttree signal fclm --matrix ./matrix.fa --taxset-csv taxsets.csv --partitions matrix.best_model.nex
```

## signal lnl — 位点级与基因级 lnL 分布

### 目的

使用 IQ-TREE3 的 `-wslr` 计算多棵候选树下位点级与基因级的对数似然分数。依据 Shen et al. (2017) 识别具有不成比例系统发育信号（ΔGLS）的基因。

### 输入

| 输入 | 说明 |
|------|------|
| `--matrix` | 单一超矩阵比对（FASTA、PHYLIP、NEXUS）。必填。对应 IQ-TREE `-s`。 |
| `--candidate-trees` | 树列表文件或以逗号分隔的独立 NEWICK 文件。必填。对应 IQ-TREE `-z`。与 `posttree topology` 格式相同。 |
| `--model-expr` | 完整的 IQ-TREE `-m` 表达式（如 `LG+F+R4`、`C20+F+R4`）。与 `--partitions` 结合使用时，同一模型会独立应用于每个分区。 |
| `--partitions` | 分区文件，按 `--partition-mode` 以 `-p` 或 `-Q` 传给 IQ-TREE。同时提取基因座边界用于基因级计算。与 `--locus-ranges` 互斥。与 `--model-expr` 结合使用时，每个分区独立估计该模型的参数。 |
| `--partition-mode` | `p` = `-p`（边连锁比例模型）；`Q` = `-Q`（边独立模型）。默认 `p`。仅当提供 `--partitions` 时有效。 |
| `--locus-ranges` | 仅用于基因座边界提取的分区文件（不传给 IQ-TREE）。与 `--partitions` 互斥。 |
| `--guide-tree` | PMSF 模型的引导树。对应 IQ-TREE `-ft`。 |
| `--metrics` | `phyloai pretree metrics` 输出的指标 CSV。生成离群 vs 非离群基因对比，以及支持不同候选树的基因群之间的两两指标对比。 |
| `--threads` | IQ-TREE `-T` 值（整数或 `auto`，默认 `auto`）。 |
| `--tool-args` | 额外的 IQ-TREE 参数。被阻止：`-s`、`-z`、`-wslr`、`--prefix`、`-p`、`-Q`。 |
| `--prefix`   | IQ-TREE 输出前缀（默认：`lnl`）。 |
| `--resume`   | 从 IQ-TREE 原生检查点恢复未完成的任务。 |

### 基因级 ΔSLS/ΔGLS 公式

- **2 棵候选树时：** ΔSLS = lnL_T1 − lnL_T2（带符号）；ΔGLS = 基因上求和的结果
- **>2 棵候选树时：** ΔSLS/ΔGLS = 所有成对 |lnL_Ta − lnL_Tb| 的均值

`support_sig` 列（|ΔGLS| ≥ 2）仅在 2 棵树对比的基因级表格中出现。

### 输出

```
runs/posttree/signal/lnl/
├── result.json
├── candidate.trees              # 合并后的树文件（仅当提供了多个树文件时）
├── site_lnl.csv                 # 位点级表格，按 ΔSLS 降序
├── site_support.pdf             # 位点支持分布柱状图
├── support_summary_sites.csv    # 各拓扑支持的位点数量
├── gene_lnl.csv                 # [若有基因座边界]
├── gene_support.pdf             # [若有基因座边界]
├── support_summary_genes.csv    # [若有基因座边界] 各拓扑支持的基因数量
├── outlier_genes.txt            # [若有基因座边界]
├── outlier_comparison.csv       # [若提供了 --metrics]
├── outlier_comparison.pdf       # [若提供了 --metrics]
├── support_comparison.csv       # [若提供了 --metrics 且有 ≥2 个支持组]
├── support_comparison.pdf       # [若提供了 --metrics 且有 ≥2 个支持组]
└── iqtree/
    ├── <prefix>.sitelh          # IQ-TREE 原始位点对数似然
    ├── <prefix>.iqtree          # IQ-TREE 原生报告
    └── <prefix>.log             # IQ-TREE 日志
```

### 示例

```bash
# 仅位点级分析，不产生基因级输出
phyloai posttree signal lnl --matrix matrix.fa --candidate-trees trees --model-expr LG+F+R4

# 带基因座范围 + 基因级输出
phyloai posttree signal lnl --matrix matrix.fa --candidate-trees trees --model-expr LG+F+R4 --locus-ranges partitions.txt

# 离群基因与正常基因，以及不同候选树支持组间的指标对比
phyloai posttree signal lnl --matrix matrix.fa --candidate-trees trees --model-expr LG+F+R4 --locus-ranges partitions.txt --metrics metrics.csv
```

---

## signal consistent — 一致基因识别

### 目的

识别基于似然信号（GLS）和基于四重奏信号（GQS）均一致支持某一候选拓扑的基因。需要恰好 2 棵候选树。GLS 使用 IQ-TREE3 计算，GQS 使用 wASTRAL 计算。参考 Shen et al. (2021)。

### 输入

| 输入 | 说明 |
|------|------|
| `--matrix` | 单一超矩阵比对。必填。 |
| `--candidate-trees` | 恰好 2 棵候选树（文件或逗号分隔）。必填。 |
| `--tree-dir` | 用于 GQS 计算的基因树文件目录。必填。 |
| `--model-expr` | IQ-TREE 模型表达式。与 `--partitions` 结合使用时，同一模型会独立应用于每个分区。 |
| `--partitions` | 分区文件，按 `--partition-mode` 以 `-p` 或 `-Q` 传给 IQ-TREE。同时提取基因座边界。 |
| `--partition-mode` | `p`（边连锁）或 `Q`（边独立）。默认 `p`。仅与 `--partitions` 配合使用。 |
| `--locus-ranges` | 仅用于基因座边界提取的分区文件。与 `--partitions` 互斥。 |
| `--guide-tree` | PMSF 引导树。 |
| `--metrics` | 指标 CSV，用于一致 vs 不一致基因对比。 |
| `--threads` | IQ-TREE `-T`（默认 `auto`）。同时控制 wASTRAL 并行度。 |
| `--tool-args` | 额外 IQ-TREE 参数。被阻止：`-s`、`-z`、`-wslr`、`--prefix`、`-p`、`-Q`。 |
| `--prefix`   | IQ-TREE 输出前缀（默认：`consistent`）。 |
| `--resume`   | 从 IQ-TREE 原生检查点恢复未完成的任务。 |
校验规则：
- 候选树恰好 2 棵（1 棵或超过 2 棵均报硬错误）。
- `--partitions`/`--locus-ranges` 中的所有基因座必须在 `--tree-dir`
  中有对应的基因树文件。`--tree-dir` 中多余的文件会被静默忽略。

### 输出

```
runs/posttree/signal/consistent/
├── result.json
├── candidate.trees              # [若合并了多个树文件]
├── gls.csv                      # 基因级 lnL 对比
├── gqs.csv                      # 基因级 GQS 对比
├── consistent_genes.txt
├── inconsistent_genes.txt
├── gls_support.pdf
├── gqs_support.pdf
├── consistent_comparison.csv    # [若提供了 --metrics]
├── consistent_comparison.pdf    # [若提供了 --metrics]
└── iqtree/
    ├── <prefix>.sitelh
    ├── <prefix>.iqtree
    └── <prefix>.log
```

### 示例

```bash
phyloai posttree signal consistent \
  --matrix matrix.fa \
  --candidate-trees T1.tre,T2.tre \
  --tree-dir gene_trees/ \
  --model-expr LG+F+R4 \
  --locus-ranges partitions.txt
```

---

## signal fclm — 四簇似然映射

### 目的

执行四簇似然映射（FcLM），评估支持四个类群簇间不同关系假设的系统发育信号。使用 IQ-TREE3 的 `-lmap` 和 `-lmclust` 参数。

### 输入

| 输入 | 说明 |
|------|------|
| `--matrix` | 单一超矩阵比对。必填。 |
| `--taxset-csv` | 两列 CSV（`taxon,taxset`），定义簇成员。至少 4 个 taxset。必填。 |
| `--model-expr` | IQ-TREE 模型表达式（如 `LG+C60+F+R4`）。与 `--partitions` 结合使用时，同一模型会独立应用于每个分区。 |
| `--partitions` | 分区文件（如 IQ-TREE 输出的 `.best_model.nex`）。按 `--partition-mode` 以 `-p` 或 `-Q` 传给 IQ-TREE。与 `--model-expr` 结合使用时，每个分区独立估计该模型的参数。 |
| `--partition-mode` | `p` = `-p`（边连锁）；`Q` = `-Q`（边独立）。默认 `p`。仅当提供 `--partitions` 时有效。 |
| `--lmap` | 四重奏采样数量：`ALL` 表示全部，整数为固定数量，不填则为 `50 × 物种数`。对应 IQ-TREE `-lmap`。 |
| `--guide-tree` | PMSF 引导树。 |
| `--threads` | IQ-TREE `-T`（默认 `auto`）。 |
| `--tool-args` | 额外 IQ-TREE 参数。被阻止：`-s`、`-lmap`、`-lmclust`、`-n`、`-p`、`-Q`、`--prefix`。 |
| `--prefix`   | IQ-TREE 输出前缀（默认：`fclm`）。 |
| `--resume`   | 从 IQ-TREE 原生检查点恢复未完成的任务。 |

校验规则：
- CSV 中所有类群名称必须与 `--matrix` 中的名称完全匹配。
- `taxset` 分配必须互斥（每个类群在且仅在一个 taxset 中）。
- 至少需要 4 个 taxset（IQ-TREE FcLM 要求）。

### 输出

```
runs/posttree/signal/fclm/
├── result.json
├── cluster.nexus                # 由 --taxset-csv 自动生成
└── iqtree/
    ├── <prefix>.lmap.eps        # IQ-TREE 似然映射图
    ├── <prefix>.iqtree          # IQ-TREE 原生报告（含所有 lmap 统计量）
    └── <prefix>.log
```

### 示例

```bash
# 同质模型
phyloai posttree signal fclm --matrix matrix.fa --taxset-csv taxsets.csv --model-expr LG+C60+F+R4

# 分区模型
phyloai posttree signal fclm --matrix matrix.fa --taxset-csv taxsets.csv --partitions matrix.best_model.nex
```

---

## 通用说明

- 三个子命令均为单矩阵模式（不支持批量）。
- `--model-expr` 与 `--partitions` 可以组合使用：`--model-expr` 指定替代
  模型公式，`--partitions` 提供分区边界。每个分区独立估计该模型公式的参数。
  `--partition-mode`（默认 `p`）控制 `--partitions` 作为 `-p`（边连锁）
  还是 `-Q`（边独立）传给 IQ-TREE。
- `--partitions` 与 `--locus-ranges`：`--partitions` 传给 IQ-TREE 并用于提取
  基因座边界；`--locus-ranges` 仅用于提取边界。
- IQ-TREE 输出文件（`.sitelh`、`.iqtree`、`.lmap.eps`、`.log`）放置在输出
  目录下的 `iqtree/` 子目录中。IQ-TREE 的标准输出在执行过程中显示在终端上。
- 默认 `--output-dir`：`runs/posttree/signal/lnl`、
  `runs/posttree/signal/consistent`、`runs/posttree/signal/fclm`。
- `--dry-run` 仅打印 IQ-TREE 命令，不实际运行。
- `signal lnl` 为 CCA 准备每个模型对应的 `site_lnl.csv`：对每个模型运行一次，
  并使用相同且顺序一致、包含两棵树的 `--candidate-trees` 输入，使每个 CSV 都含
  `lnL_Tree1` 与 `lnL_Tree2`。CCA 输入准备和解释见
  [系统误差工作流参考](../../skills/phyloai-workflow/references/syserror-workflow.md)。
- 参考文献：Shen et al. (2017) *Nature Ecology & Evolution*；Shen et al. (2021) *Systematic Biology*。

## 退出码

| 代码 | 含义 |
|------|------|
| 0 | 成功 |
| 1 | 用户输入错误（文件缺失、参数无效、输出冲突） |
| 2 | 外部工具执行失败 |
| 3 | 未找到外部工具可执行文件 |
