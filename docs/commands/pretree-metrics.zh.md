# phyloai pretree metrics

[English](pretree-metrics.md) | [中文](pretree-metrics.zh.md)

## 目的

`phyloai pretree metrics` 从 MSA 和/或基因树计算分子标记属性，用于统计探索与下游过滤。它在标准流水线中位于 `pretree trim` 之后、`pretree filter` 之前。

命令组有三个入口：

- `phyloai pretree metrics` 计算度量，写出 `metrics.csv`，生成分布图，创建相关性热图。
- `phyloai pretree metrics plot` 从已有 `metrics.csv` 重新生成单个度量分布图。
- `phyloai pretree metrics correlate` 从已有 `metrics.csv` 重新生成相关性矩阵与热图。

它不运行标记过滤、UMAP 聚类或树推断。

## 用法

从 MSA 文件计算度量：

```bash
phyloai pretree metrics --msa-dir ./runs/pretree/trim/seqs/faa
```

同时计算 MSA 与基因树度量：

```bash
phyloai pretree metrics \
  --msa-dir ./runs/pretree/trim/seqs/faa \
  --tree-dir ./data/gene_trees \
  --output-dir ./runs/pretree/metrics
```

不写文件地预览工作：

```bash
phyloai pretree metrics --msa-dir ./msa --dry-run
```

用自定义样式与 Tukey 过滤重新绘制单个度量：

```bash
phyloai pretree metrics plot \
  --csv ./runs/pretree/metrics/metrics.csv \
  --metric entropy \
  --tukey-k 1.5 \
  --color '#2E86AB' \
  --title 'Entropy distribution'
```

重新生成紧凑的上三角相关性热图：

```bash
phyloai pretree metrics correlate \
  --csv ./runs/pretree/metrics/metrics.csv \
  --triangle upper \
  --label-angle 60 \
  --output-dir ./runs/pretree/metrics/correlate
```

包含所有数值列（含 `freq*` 与 `sd_*`），并以 full 模式绘制聚类矩形：

```bash
phyloai pretree metrics correlate \
  --csv ./runs/pretree/metrics/metrics.csv \
  --metrics all \
  --triangle full \
  --cluster-rectangles 5
```

## 参数

### `phyloai pretree metrics`

| Parameter | Default | Description |
|---|---|---|
| `--msa-dir` | none | 已比对 FASTA 文件目录（`.fa`、`.fasta`、`.fas`、`.fna`、`.faa`、`.aln`）。`--msa-dir` 与 `--tree-dir` 至少需要一个。 |
| `--tree-dir` | none | Newick 树文件目录（`.tre`、`.tree`、`.nwk`、`.newick`、`.treefile`、`.bestTree`、`.contree`）。 |
| `--seq-type` | `auto` | 分子类型：`AA`、`NT`，或每个标记自动检测。 |
| `--outgroup-list` | none | 每行列出一个外类群名用于 DVMC 剪枝；需要 `--tree-dir`。 |
| `--ref-tree` | none | 用于归一化 RF 距离的参考物种树；需要 `--tree-dir`。 |
| `--skip-freq-statistics` | off | 跳过每字符频率列（`freqA`、`freqC` 等）。 |
| `--pseudo-tree-metrics` | off | 计算 FastTree 派生的伪树度量，使用 `_FT` 后缀；需要 `--msa-dir`。 |
| `--fasttree-path` | `FastTree` | 显式 FastTree 可执行路径。 |
| `--skip-pairwise-identity` | off | 跳过 `average_pairwise_identity`；分类单元多的标记建议跳过。 |
| `--round` | `6` | 数值 CSV 值的小数位数；范围 0-12。 |
| `--table-format` | `csv` | 辅助表格输出格式：`csv` 或 `tsv`。所有辅助表（`metrics`、`basic_statistics`、`correlation_matrix`）使用同一格式。 |
| `--output-dir`, `-o` | `runs/pretree/metrics` | 输出目录，存放度量表、图、相关性输出与 `result.json`。 |
| `--threads`, `-t` | `4` | 工作进程数；至少 1。 |
| `--dry-run` | off | 校验输入并显示计划工作，不写文件。 |
| `--overwrite` | off | 删除并重建非空输出目录。 |
| `--quiet`, `-q` | off | 抑制终端进度与汇总输出。 |

### `phyloai pretree metrics plot`

| Parameter | Default | Description |
|---|---|---|
| `--csv` | required | 已有 `metrics.csv`（或 `.tsv`）。 |
| `--input-format` | `auto` | 输入文件格式：`csv`、`tsv` 或 `auto`（通过内容检测 —— tab/comma 计数，失败则回退到扩展名）。 |
| `--metric` | required | 精确的度量列名。 |
| `--bins` | `50` | 直方图分箱数；有效范围 1-500。 |
| `--xmin` | auto | 强制 X 轴下限。 |
| `--xmax` | auto | 强制 X 轴上限。 |
| `--tukey-k` | disabled | Tukey Fences 乘数。设置时，被过滤的位点写入 `<output_dir>/<metric>.tukey_filtered.csv`。 |
| `--title` | `Distribution of <metric>` | 图标题。 |
| `--xlabel` | 度量显示名 | X 轴标签。 |
| `--ylabel` | `Density` | Y 轴标签。 |
| `--color` | `#2E86AB` | 直方图条形填充颜色。 |
| `--fig-width` | `10.0` | 图宽（英寸）。 |
| `--fig-height` | `8.0` | 图高（英寸）。 |
| `--dpi` | `150` | 输出分辨率。 |
| `--font-size` | `12` | 基础字号。 |
| `--output-dir`, `-o` | `<csv_parent>/plot_<metric>/` | PDF 与 `result.json` 的输出目录。 |
| `--overwrite` | off | 替换已有的输出 PDF/目录。 |
| `--quiet`, `-q` | off | 抑制终端输出。 |

### `phyloai pretree metrics correlate`

| Parameter | Default | Description |
|---|---|---|
| `--csv` | required | 已有 `metrics.csv`（或 `.tsv`）。 |
| `--input-format` | `auto` | 输入文件格式：`csv`、`tsv` 或 `auto`（通过内容检测 —— tab/comma 计数，失败则回退到扩展名）。 |
| `--metrics` | core numeric | 逗号分隔的度量列。使用 `all` 包含所有数值列。省略意味着自动选择可读的核心度量。 |
| `--include-freq` | off | 在自动度量选择中包含 `freq*` 列。 |
| `--include-sd` | off | 在自动度量选择中包含 `sd_*` 列。 |
| `--method` | `spearman` | 相关性方法：`spearman` 或 `pearson`。 |
| `--triangle` | `full` | 矩阵显示：`full`、`lower` 或 `upper`。Lower 模式使用左/下标签；upper 模式使用上/右标签。 |
| `--annot` / `--no-annot` | `--no-annot` | 在单元格中显示数值相关性。 |
| `--cluster-rectangles` | none | 仅在 full 矩阵上绘制 N 个聚类矩形。若与 `--triangle lower` 或 `upper` 同时使用，PhyloAI 警告并忽略。 |
| `--cmap` | `RdBu_r` | Matplotlib colormap。 |
| `--fmt` | `.2f` | 标注的数值格式。 |
| `--fig-width` | `12.0` | 图宽（英寸）。 |
| `--fig-height` | `10.0` | 图高（英寸）。 |
| `--dpi` | `150` | 输出分辨率。 |
| `--font-size` | `10` | 度量标签的基础字号。 |
| `--label-angle` | `45.0` | X 轴度量标签的旋转角度（度）。 |
| `--title` | none | 可选图标题。 |
| `--output-dir`, `-o` | `runs/pretree/metrics/correlate` | 存放 `correlation_heatmap.pdf`、`correlation_matrix.csv` 与 `result.json` 的目录。 |
| `--overwrite` | off | 替换已有相关性输出。 |
| `--quiet`, `-q` | off | 抑制终端输出与非关键警告。 |

## 输入

MSA 与树文件按共享的全局匹配策略以逻辑位点配对。

- MSA 逻辑位点：文件名最后一个 `.` 之前的全部内容。
- 树逻辑位点：尝试移除最后一个后缀段，然后移除最后两个后缀段。
- 若恰好有一个树候选与某 MSA 位点匹配，则使用该对。
- 若两个树候选匹配不同的位点，PhyloAI 退出并报歧义错误，而非猜测。

示例：`EOG090X002Z.fas` -> `EOG090X002Z`；`EOG090X002Z.fas.treefile` 尝试 `EOG090X002Z.fas` 与 `EOG090X002Z`。

## 输出

### 主 `metrics` 输出目录

所有辅助表遵循 `--table-format`（默认 `csv`，产出 `.csv` 文件）。若给定 `--table-format tsv`，则改写 `.tsv` 文件。

| 文件或目录 | 描述 |
|---|---|
| `metrics.csv`（或 `.tsv`） | 每个标记一行，含标识符、MSA 度量、树度量、可选频率列与可选伪树度量。 |
| `plots/` | 每个数值度量的密度直方图 PDF。 |
| `metrics.basic_statistics.csv`（或 `.tsv`） | 每个度量的均值、中位数、最小值、最大值、q25、q75、标准差、非 NA 计数与总计数。 |
| `correlate/correlation_heatmap.pdf` | 默认的核心数值度量紧凑相关性热图。 |
| `correlate/correlation_matrix.csv`（或 `.tsv`） | 相关性矩阵，行列标签为度量名。 |
| `result.json` | 结构化的状态、参数、关键计数、警告与数据路径。 |

### 重要的 `metrics.csv` 列

- `loci` 与 `DataType` 标识标记及其分子类型。
- MSA 度量包括 `num_taxa`、`taxa_occupancy`、`num_sites`、`num_patterns`、`proportion_patterns`、`num_parsimony_sites`、`proportion_parsimony`、`num_singletons`、`proportion_singletons`、`proportion_gaps`、`proportion_invariant`、`entropy`、`bollback`、`pattern_entropy`、`rcfv`、`nrcfv`、`average_pairwise_identity` 与 `GC_content`。
- 树度量包括 `average_BS`、`sd_BS`、`total_tree_length`、支长汇总、patristic distance 汇总、`evo_rate`、`treeness`、`dvmc`、`saturation` 与 `RF_distance`。
- 频率度量使用 `freq*` 列名，设置 `--skip-freq-statistics` 时被省略。
- 伪树度量使用 `_FT` 后缀，仅在 `--pseudo-tree-metrics` 下产生。

## 相关性说明

默认相关性图刻意排除标识列、`freq*` 与 `sd_*` 列，以保持 PDF 可读。当需要可视化这些列时，使用 `--include-freq`、`--include-sd`、显式 `--metrics` 或 `--metrics all`。

变量按基于 `1 - |corr|` 的 Ward 聚类排序，但不绘制树状图。Triangle 模式仅绘制可见的那一半，在对角线单元格之外绘制阶梯状三角边框。`--triangle upper` 将度量标签置于上/右，并将 colorbar 移到左侧；`--triangle lower` 将标签置于左/下，colorbar 保持在右侧。Full 矩阵的 colorbar 也在右侧。

## 警告与错误

| 情况 | 行为 |
|---|---|
| 同时缺少 `--msa-dir` 与 `--tree-dir` | 报错退出。 |
| `--pseudo-tree-metrics` 但无 `--msa-dir` | 报错退出。 |
| `--outgroup-list` 或 `--ref-tree` 但无 `--tree-dir` | 报错退出。 |
| 输出目录非空且未加 `--overwrite` | 主 `metrics` 报错退出。 |
| MSA/树分类单元不匹配 | 记录警告；计算继续。 |
| 未配对的 MSA/树文件 | 为每个未匹配 stem 记录警告。 |
| 标记含 >200 个分类单元且启用了 pairwise identity | 警告建议 `--skip-pairwise-identity`。 |
| FastTree 不可用 | 伪树度量为空；核心度量继续。 |
| `--cluster-rectangles` 与 `--triangle lower/upper` | 除 `--quiet` 外打印警告；忽略矩形。 |

## 备注

- `metrics.csv` 是被 `pretree filter` 消费的规范中间产物。
- `rcfv` 是经典的相对组成频率变异性；`nrcfv` 是 Fleming 和 Struck 提出的偏差修正度量。
- `saturation` 是 patristic 距离对未校正序列距离经过原点的斜率。
- `average_pairwise_identity` 是 O(n² × L)；对分类单元多的数据集使用 `--skip-pairwise-identity`。