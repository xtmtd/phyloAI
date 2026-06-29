# phyloai pretree filter

[English](pretree-filter.md) | [中文](pretree-filter.zh.md)

## 目的

`phyloai pretree filter` 是介于树推断与超矩阵拼接之间的质量控制闸门。它提供五种互补的过滤工作流：

- `taper` —— 屏蔽单个位点内的错误位点（位点级）
- `treeshrink` —— 从基因树中剪除长支外类群（分类单元级）
- `metrics` —— 根据度量表上的数值/字符串规则丢弃或保留整个基因座（基因座级）
- `cluster` —— 按度量特征对基因座分组；可选地丢弃离群簇（群体级）

Filter **不**计算标记度量（使用 `phyloai pretree metrics`）、推断基因树，也不拼接保留的 MSA（使用 `phyloai pretree concat`）。它按需读取 `pretree metrics` 的输出表，写出结构化决策文件以及可选的过滤后 MSA/树目录。

模块将度量计算与过滤决策分开，以便在不重计算度量的情况下探索多种阈值组合。

流水线中的位置：

```
phyloai pretree align    →  已比对 MSA
phyloai pretree trim     →  修剪后 MSA
phyloai pretree filter   →  质量控制后的 MSA 与树  ← 你在这里
phyloai pretree concat   →  超矩阵
```

### 推荐工作流

对于典型的系统发育 MSA，五个子命令按顺序使用：

1. **`taper`** —— 在修剪后的 MSA 上屏蔽潜在的位点级错误。这能产出更干净的比对而不丢弃基因座或分类单元。

2. **从屏蔽后的 MSA 构建基因树**（使用外部树推断工具）。这些树反映了修正后的序列。

3. **`treeshrink`** —— 将基因树（以及可选的屏蔽后 MSA）喂给 TreeShrink 以识别并剪除外类群长支。结果是一组剪缩后的树，以及可选的剪缩后 MSA（已剔除问题分类单元）。

4. **（可选）重新推断基因树** 基于剪缩后的 MSA —— TreeShrink 保证了长支分类单元被移除，但使用修剪后的比对可能进一步改善树拓扑。

5. **`metrics`**、**`symtest`** 与/或 **`cluster`** —— 仅在位点屏蔽与分类单元剪除之后才评估每个基因座的质量。这些子命令可选且互补：`metrics` 计算用于规则过滤的标记属性；`symtest` 通过 IQ-TREE3 测试系统发育对称性假设（平稳性/同质性）并按 p 值过滤；`cluster` 按度量特征对基因座分组用于探索与离群移除。根据需要选择其一或多个在清洗后的数据集上运行。

子命令也可以独立使用。例如，如果你已有基因树且只想剪除分类单元，从第 3 步开始。如果只需要应用度量阈值，直接跳到 `filter metrics`。

所有子命令将 `result.json` 写入各自的输出目录。终端输出使用 Rich 表格；用 `--quiet` 抑制。

这些子命令写出的任何 PhyloAI 作者化 FASTA 系列输出均按 60 字符换行。

### 共享选项

| Option | Default | Purpose |
|--------|---------|---------|
| `--output-dir` / `-o` | `runs/pretree/filter/<subcommand>` | 输出目录 |
| `--table-format` | `csv` | 辅助表的分隔符与后缀；不影响 `result.json` |
| `--overwrite` | off | 删除并重建 `--output-dir`（若已存在） |
| `--dry-run` | off | 校验输入并显示计划动作；不写文件（不写 `result.json`） |
| `--quiet` / `-q` | off | 除错误外不打印终端输出 |

### 文件匹配策略

当接受 `--msa-dir` 或 `--tree-dir` 时，所有子命令使用来自 `phyloai/core/file_matching.py` 的、与扩展名无关的逻辑位点匹配：

| 文件 | 逻辑位点 |
|------|----------|
| `gene1.fa` | `gene1` |
| `gene2.v1.ALI` | `gene2.v1` |
| `gene3.treefile` | `gene3` |
| `gene4.fa.treefile` | `gene4.fa`，然后 `gene4` |

每个常规非空文件都会被扫描；格式在解析时校验，而非通过扩展名。有歧义的树匹配（文件名候选已被占用）会报错。`phyloai pretree metrics` 使用同样的助手，因此 `metrics` 与 `filter` 对非标准命名的行为一致。

---

## `filter taper` —— TAPER 错误位点屏蔽

### 目的

运行 TAPER 屏蔽 MSA 内错误的氨基酸或核苷酸位点。TAPER 识别与比对其余部分相比异常分歧的连续残基段，并将其替换为 `X`（AA）或 `N`（NT）。仅统计新引入的屏蔽；忽略输入中原始的歧义字符。

这是位点级质量控制：基因座被保留，但有问题的位置被中和。TAPER **不**移除基因座（使用 `filter metrics`）或移除分类单元（使用 `filter treeshrink`）。

### 用法

```bash
phyloai pretree filter taper \
  --msa-dir <aa_or_nt_msa_dir> \
  [--nt-dir <codon_aligned_nt_msa_dir>] \
  [--seq-type AA|NT|auto] \
  [--cutoff 3] \
  [--taper-path <correction_multi.jl>] \
  [--julia-path <julia>] \
  [--tool-args "..."] \
  [--show-masked-sites] \
  [--output-dir runs/pretree/filter/taper] \
  [--threads 4] [--resume] [--dry-run] [--overwrite]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--msa-dir` | required | 已比对 MSA 文件目录。扫描每个常规非空文件。 |
| `--nt-dir` | — | AA+CDS 模式的密码子比对 NT MSA 目录（NT 长度 == 3 × AA 长度）。 |
| `--seq-type` | `auto` | `AA`、`NT` 或 `auto`（从首个文件检测：`EFILPQWYZ` → AA）。 |
| `--cutoff` | 3 | TAPER 校正 cutoff（≥1）。越低越激进。3 为 TAPER 默认；噪声数据用 1-2，保守屏蔽用 5+。 |
| `--taper-path` | — | 显式 `correction_multi.jl` 路径。默认使用捆绑副本。 |
| `--julia-path` | — | 显式 Julia 可执行路径。通过 PATH 解析；用 `phyloai doctor` 验证。 |
| `--tool-args` | — | 透传其他 TAPER 标志。不能包含 `-m`、`-a`、`-c`、`-l`、输入路径或输出重定向（这些由 PhyloAI 管理；出现则报错）。 |
| `--threads` / `-t` | 4 | 工作进程数，每个工作进程处理一个位点。使用 `ProcessPoolExecutor` 与 checkpoint/resume（与 `pretree align`、`pretree trim` 模式相同）。 |
| `--show-masked-sites` | off | 在 `filter_decisions.csv` 中添加 `masked_taxa_detail` 列（`taxonA:3; taxonB:5`）。 |
| `--table-format` | `csv` | `retained_loci`、`dropped_loci`、`filter_decisions` 的格式。 |
| `--resume` | off | 从 `checkpoint.json` 恢复；参数必须匹配；输出通过校验的已完成位点会被跳过。 |
| `--overwrite` | off | 与 `--resume` 互斥。 |
| `--dry-run` | off | 显示检测到的模式、配对位点、命令模板、输出布局。不写文件。 |
| `--output-dir` / `-o` | `runs/pretree/filter/taper` | 输出目录。 |
| `--quiet` / `-q` | off | 除错误外不打印终端输出。 |

### 输入

三种运行模式：

| 模式 | 输入 | 输出 |
|------|------|------|
| AA-only | AA MSA 文件位于 `--msa-dir` | 屏蔽后 AA → `seqs/` |
| NT-only | NT MSA 文件 + `--seq-type NT` | 屏蔽后 NT → `seqs/` |
| AA+CDS | AA MSA + `--nt-dir` 密码子比对 NT | 屏蔽后 AA → `seqs/faa/`，投影后 CDS → `seqs/fna/` |

AA+CDS 模式要求：NT 记录形成有效密码子 MSA（等长、能被 3 整除）；每个位点 AA 与 NT 分类单元完全匹配；AA 长度 == NT 长度 / 3。

AA+CDS 的投影规则：
- 输入中原始的 `X` → 不变（不算作 TAPER 屏蔽）
- 标准 AA → 被 TAPER 改为 `X` → 对应密码子替换为 `NNN`
- Gap `-` → `X` → 警告；CDS 不变（防御性检查，正常情况下不应出现）

### 输出

```
runs/pretree/filter/taper/
├── seqs/                              （AA+CDS 模式下为 seqs/faa/ + seqs/fna/）
├── retained_loci.csv|tsv
├── dropped_loci.csv|tsv               （locus, reason）
├── filter_decisions.csv|tsv           （locus, status, new_masked_sites, masked_taxa_count,
│                                       使用 --show-masked-sites 时含 masked_taxa_detail）
├── checkpoint.json                    （内部；仅在使用 --resume 时存在）
├── logs/
│   ├── gene1.log
│   └── ...
├── result.json
```

终端输出：两张 Rich 表格 —— Filter Results（输入/保留/丢弃/屏蔽位点/分类单元/位点数）与 Retained MSA Statistics（MSA 数、总/平均/最小/最大比对长度、平均分类单元数）。Julia 版本通过 `julia -v` 自动检测并记录在 `result.json` 中。

### 示例

```bash
# 默认 AA 屏蔽
phyloai pretree filter taper --msa-dir ./trimmed

# 对噪声数据激进屏蔽
phyloai pretree filter taper --msa-dir ./trimmed --cutoff 1 --threads 8

# 保守屏蔽
phyloai pretree filter taper --msa-dir ./trimmed --cutoff 5

# NT-only 模式
phyloai pretree filter taper --msa-dir ./trimmed_nt --seq-type NT

# AA+CDS：屏蔽 AA，投影到密码子比对 NT
phyloai pretree filter taper --msa-dir ./trimmed_aa --nt-dir ./trimmed_fna

# 中断后恢复
phyloai pretree filter taper --msa-dir ./trimmed --resume

# 包含每分类单元屏蔽细节以便检查
phyloai pretree filter taper --msa-dir ./trimmed --show-masked-sites
```

### 警告与错误

| 条件 | 行为 |
|------|------|
| `--nt-dir` 与 `--seq-type NT` 同时使用 | Exit 1 |
| `--threads` < 1 | Exit 1 |
| `--resume` + `--overwrite` | Exit 1 |
| 未找到 Julia | Exit 3 |
| 输出目录非空且未加 `--overwrite` 或 `--resume` | Exit 1 |
| `--msa-dir` 内无有效 MSA 文件 | Exit 1 |
| TAPER 对某位点退出非零 | 该位点被跳过；原因记录在 `dropped_loci.csv` |
| TAPER 输出缺失或 FASTA 校验失败 | 该位点被跳过 |
| 所有位点失败 | Exit 2 |

### 备注

TAPER 始终是第一个过滤步骤，因为屏蔽应发生在树推断之前。屏蔽后，用 `phyloai pretree metrics` 在屏蔽后的 MSA（以及可选的重新推断基因树）上计算每个位点的度量，然后应用 `filter metrics` 或 `filter cluster`。

支持 `--resume` 是因为对大量位点进行屏蔽计算密集。checkpoint 模式与 `pretree align` 和 `pretree trim` 一致。

---

## `filter treeshrink` —— TreeShrink 分类单元剪除

### 目的

运行 TreeShrink 检测并移除基因树中的长支离群分类单元。TreeShrink 使用统计检验在多个树之间联合识别分支异常长的分类单元。当提供 `--msa-dir` 时，匹配的 MSA 也会被剪缩以移除相同的剪除分类单元。

这是分类单元级过滤：从特定基因树中移除分类单元。基因座被保留。TreeShrink **不**移除整个基因座（使用 `filter metrics`）。

### 用法

```bash
phyloai pretree filter treeshrink \
  --tree-dir <gene_tree_dir> \
  [--msa-dir <msa_dir>] \
  [--threshold 0.05] \
  [--treeshrink-mode auto|per-gene|all-genes|per-species] \
  [--treeshrink-path <run_treeshrink.py>] \
  [--tool-args "..."] \
  [--output-dir runs/pretree/filter/treeshrink] \
  [--keep-work-dir] [--dry-run] [--overwrite]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--tree-dir` | required | 基因树文件目录。扫描每个常规非空文件。 |
| `--msa-dir` | — | 可选 MSA 目录。按逻辑位点名配对文件，与树一起剪缩。 |
| `--threshold` | 0.05 | TreeShrink 假阳性阈值（`-q`）。越小移除越多分类单元。0.05 为 TreeShrink 默认。 |
| `--treeshrink-mode` | `auto` | `auto` 省略 `-m`（TreeShrink 默认）。`per-gene`：每基因独立。`all-genes` / `per-species`：跨基因合并。 |
| `--treeshrink-path` | — | 显式 `run_treeshrink.py` 路径。通过 PATH 解析。 |
| `--tool-args` | — | 其他 TreeShrink 标志。不能包含 `-i`、`-t`、`-a`、`-q`、`-m`、`-o`、`-O`。 |
| `--keep-work-dir` | off | 保留 `output_dir/work/` 下的每基因工作目录以便调试。 |
| `--output-dir` / `-o` | `runs/pretree/filter/treeshrink` | 输出目录。 |
| `--table-format` | `csv` | 辅助表的格式。 |
| `--overwrite` | off | 删除并重建输出目录。 |
| `--dry-run` | off | 打印解析后的命令与位点数。 |
| `--quiet` / `-q` | off | 抑制终端输出。 |

### 输入

`--tree-dir` 是必需的。TreeShrink 在整个数据集上仅调用一次（而非按位点），因为其统计模型可以从多个树中合并信息。PhyloAI 在临时目录中创建一个按基因的工作布局：

```
<work_dir>/input/
├── gene1/
│   ├── input.tree
│   └── input.fasta     （仅当提供 --msa-dir 时）
├── gene2/
│   ├── input.tree
│   └── input.fasta
```

### 输出

```
runs/pretree/filter/treeshrink/
├── trees/                              （剪缩后的基因树）
├── seqs/                               （仅当提供 --msa-dir 时）
├── retained_loci.csv|tsv
├── modified_loci.csv|tsv               （剪除 ≥1 个分类单元的位点）
├── dropped_loci.csv|tsv                （输出缺失/无效的位点）
├── removed_taxa.csv|tsv                （每行：locus, taxon）
├── filter_decisions.csv|tsv            （locus, status, removed_count）
├── work/                               （仅当使用 --keep-work-dir 时）
├── logs/treeshrink.log                 （单一共享工具 stderr）
└── result.json
```

决策类别：retained（含未修改）、modified（剪除分类单元）、dropped（输出缺失）。

终端输出：Filter Results 表（输入/保留/修改/丢弃/移除分类单元数）+ Retained MSA Statistics 表（提供 `--msa-dir` 时）。一个上下文提示会提醒用户过滤后的比对可用于重建系统发育树，可能比 TreeShrink 剪除的版本更准确。

### 示例

```bash
# 基础分类单元剪除
phyloai pretree filter treeshrink --tree-dir ./genetrees

# 树 + 匹配 MSA，保守阈值
phyloai pretree filter treeshrink \
  --tree-dir ./genetrees --msa-dir ./trimmed --threshold 0.1

# per-species 模式（跨基因合并）
phyloai pretree filter treeshrink \
  --tree-dir ./genetrees --treeshrink-mode per-species

# 调试输出
phyloai pretree filter treeshrink --tree-dir ./genetrees --keep-work-dir
```

### 警告与错误

| 条件 | 行为 |
|------|------|
| 未找到 `run_treeshrink.py` | Exit 3 |
| `--tree-dir` 内无有效树文件 | Exit 1 |
| 树与 MSA 之间位点匹配有歧义 | Exit 1 并附细节 |
| TreeShrink 退出非零 | 所有位点标记失败 |
| 所有位点失败 | Exit 2 |
| 输出目录非空且未加 `--overwrite` | Exit 1 |

### 备注

不支持 `--resume` 和 `--threads`：TreeShrink 在整个数据集上仅运行一次，按位点并行化会改变其统计模型。TreeShrink 的 `-q` 阈值控制假阳性率；0.05 意味着错误移除分类单元的概率约为 5%。

TreeShrink 之后，你可能希望基于剪缩后的 MSA 重新推断基因树以获得更准确的拓扑，然后计算度量并过滤。

---

## `filter metrics` —— 度量规则过滤

### 目的

按度量 CSV/TSV 表（通常来自 `phyloai pretree metrics`）上的显式数值或字符串条件过滤整个位点。`--keep` 中的所有条件按 AND 组合：一个位点必须满足所有条件才会被保留。

这是基因座级过滤：整个基因被保留或丢弃。`filter metrics` **不**计算度量（使用 `pretree metrics`），也不在位点/分类单元级过滤（使用 `filter taper` 或 `filter treeshrink`）。

### 用法

```bash
phyloai pretree filter metrics \
  --table <metrics.csv|metrics.tsv> \
  --keep "col>=val,col<=val,..." \
  [--input-format auto|csv|tsv] \
  [--loci-column loci] \
  [--msa-dir <msa_dir>] [--tree-dir <tree_dir>] [--copy] \
  [--output-dir runs/pretree/filter/metrics] \
  [--table-format csv|tsv] [--dry-run] [--overwrite]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--table` | required | 度量 CSV/TSV 路径。除非指定 `--input-format`，分隔符自动检测。 |
| `--keep` | required | 逗号分隔的 AND 条件。操作符：`>=`, `>`, `<=`, `<`, `==`, `!=`。字符串列仅允许 `==`/`!=`。 |
| `--input-format` | `auto` | `csv`、`tsv` 或 `auto`。自动检测有歧义时使用。 |
| `--loci-column` | `loci` | 持有位点标识符的列名。 |
| `--msa-dir` | — | 用于保留 MSA 统计的 MSA 目录（终端 + `result.json`）。 |
| `--tree-dir` | — | `--copy` 模式下的树目录。 |
| `--copy` | off | 将保留的 MSA/树复制到输出目录。需要 `--msa-dir` 或 `--tree-dir`。 |
| `--output-dir` / `-o` | `runs/pretree/filter/metrics` | 输出目录。 |
| `--table-format` | `csv` | `retained_loci`、`dropped_loci`、`filter_decisions` 的格式。 |
| `--overwrite` | off | 删除并重建输出目录。 |
| `--dry-run` | off | 解析规则并报告通过/未通过的位点数，不写文件。 |
| `--quiet` / `-q` | off | 抑制终端输出。 |

### 输入

`--table` 文件必须是带表头的 CSV 或 TSV。位点标识符列（默认 `loci`）标识每行。除非指定，否则分隔符通过前 1024 字节检测。空文件会报错。

规则语法：

```
column operator value
```

```bash
# 数值阈值（仅 AND）
--keep "dvmc>=0,dvmc<=0.3,average_BS>=0.8"

# 数值 + 字符串混合
--keep "DataType==AA,num_sites>=300"

# 单一条件
--keep "num_sites>=1000"
```

在字符串列上使用 `>=`/`>`/`<=`/`<` 会报错退出。v1 不支持 OR 逻辑。

### 输出

```
runs/pretree/filter/metrics/
├── retained_loci.csv|tsv
├── dropped_loci.csv|tsv
├── filter_decisions.csv|tsv
├── seqs/                              （仅在使用 --copy --msa-dir 时）
├── trees/                             （仅在使用 --copy --tree-dir 时）
├── result.json
```

终端输出：Filter Results 表（总数/保留/丢弃）+ Retained MSA Statistics 表（提供 `--msa-dir` 时）。

`result.json.key_results.condition_failure_counts` 将每个条件映射到其拒绝的位点数 —— 用它来识别最严格的阈值。

### 示例

```bash
# 基础质量过滤
phyloai pretree filter metrics \
  --table ./metrics/metrics.csv \
  --keep "dvmc<=0.3,average_BS>=0.8"

# 数值 + 字符串混合 + 文件复制
phyloai pretree filter metrics \
  --table ./metrics.csv \
  --keep "DataType==AA,num_sites>=300" \
  --copy --msa-dir ./trimmed

# 先 dry-run 探索阈值再写入
phyloai pretree filter metrics \
  --table ./metrics.csv \
  --keep "num_sites>=500,average_BS>=0.7" \
  --dry-run

# 自定义位点列名
phyloai pretree filter metrics \
  --table ./table.tsv \
  --keep "average_BS>=0.9" \
  --loci-column gene_id

# 比较两种策略（不同输出目录）
phyloai pretree filter metrics \
  --table ./metrics.csv --keep "average_BS>=0.8" \
  -o ./runs/strategy_conservative

phyloai pretree filter metrics \
  --table ./metrics.csv --keep "average_BS>=0.5" \
  -o ./runs/strategy_lenient
```

### 警告与错误

| 条件 | 行为 |
|------|------|
| `--table` 不存在 | Exit 1 |
| `--table` 为空 | Exit 1 |
| `--keep` 语法错误 | Exit 1 并附解析错误详情 |
| `--keep` 引用未知列 | Exit 1 |
| 在字符串列上使用数值操作符（`>=`, `>`, `<=`, `<`） | Exit 1 |
| `--copy` 没有 `--msa-dir` 或 `--tree-dir` | Exit 1 |
| 没有位点满足所有条件 | 结果报告保留数为 0（不是错误） |
| 输出目录非空且未加 `--overwrite` | Exit 1 |

### 备注

`filter metrics` 故意与 `pretree metrics` 计算分离，以便在不重计算度量的情况下探索阈值组合。使用 `--dry-run` 快速迭代，然后用 `--copy` 应用最终阈值生成过滤文件。

没有 `--copy` 时仅写入决策表 —— 对阈值探索既快又省盘。`result.json` 中的 `condition_failure_counts` 准确显示哪条条件丢弃最多位点。

若要实现 OR 逻辑，用不同的 `--keep` 与 `--output-dir` 跑两次该命令。

---

## `filter symtest` —— 对称性检验过滤

### 目的

通过 IQ-TREE3 的对称性检验（Naser-Khdour 等，2019）调用 `--symtest-only` 检测违反系统发育平稳性、同质性或可逆性假设的位点，然后过滤掉 p 值低于可配置阈值的位点。

这是按统计检验的基因座级过滤：未通过对称性的位点被丢弃。这不计算度量（使用 `pretree metrics`），也不屏蔽/剪除单个位点/分类单元（使用 `filter taper` 或 `filter treeshrink`）。

### 用法

```bash
phyloai pretree filter symtest \
  --msa-dir <msa_dir> \
  [--symtest-type MAR|INT] \
  [--symtest-pval 0.05] \
  [--symtest-keep-zero] \
  [--iqtree-path <path>] \
  [--threads 4] \
  [--tree-dir <tree_dir>] \
  [--output-dir runs/pretree/filter/symtest] \
  [--table-format csv|tsv] \
  [--dry-run] [--quiet] [--overwrite]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--msa-dir` | required | 每个位点 MSA 文件目录。扫描每个常规非空文件。 |
| `--symtest-type` | — | `MAR`（marginal/平稳性检验），`INT`（internal/同质性检验）。省略时使用合并的 Sym 检验（SymPval 列）。 |
| `--symtest-pval` | 0.05 | P 值阈值。p ≥ 阈值的位点被保留；低于阈值的被丢弃。必须 > 0 且 ≤ 1。 |
| `--symtest-keep-zero` | off | 将 `--symtest-keep-zero` 传给 IQ-TREE（在检验中保留 NA）。 |
| `--iqtree-path` | — | 显式 iqtree 二进制路径。要求 IQ-TREE3 >= 2.3.0。 |
| `--threads` / `-t` | 4 | IQ-TREE 线程数（`-T`）。 |
| `--tree-dir` | — | 可选基因树目录。与保留位点匹配的树被复制到 `trees/`。 |
| `--output-dir` / `-o` | `runs/pretree/filter/symtest` | 输出目录。 |
| `--table-format` | `csv` | 辅助表的格式。 |
| `--overwrite` | off | 删除并重建输出目录。 |
| `--dry-run` | off | 显示解析后的 IQ-TREE 命令与位点数。 |
| `--quiet` / `-q` | off | 抑制终端输出。 |

### 输入

`--msa-dir` 是唯一必需的输入。命令流程：

1. 从所有 MSA 构建临时超矩阵 + RAxML 风格分区文件
2. 运行 `iqtree -s <matrix> -p <partitions> --symtest-only`
3. 从 `<partitions>.symtest.csv` 解析各分区的 p 值
4. 在所选检验列上应用 p 值阈值

使用的 p 值列取决于 `--symtest-type`：
- （默认） -> `SymPval`（合并的平稳性 + 同质性）
- `MAR` -> `MarPval`（marginal / 平稳性检验）
- `INT` -> `IntPval`（internal / 同质性检验）

临时文件在运行后被清理。

### 输出

```
runs/pretree/filter/symtest/
├── seqs/                              （保留的 MSA）
├── trees/                             （仅当提供 --tree-dir 时）
├── retained_loci.csv|tsv
├── dropped_loci.csv|tsv               （locus, reason）
├── filter_decisions.csv|tsv           （locus, status, p_value, symtest_type,
│                                       sym_pval, mar_pval, int_pval,
│                                       sym_sig, sym_non, mar_sig, mar_non,
│                                       int_sig, int_non）
├── logs/symtest.log                    （单一共享工具 stderr）
└── result.json
```

终端输出：Filter Results 表（输入/保留/丢弃/p 值阈值/symtest 类型）+ Retained MSA Statistics 表 + 可选的 Trees Copied 表。

`result.json.key_results`：`n_input`、`n_retained`、`n_dropped`、`p_value_threshold`、`symtest_type`、`retained_trees_copied`。

### 示例

```bash
# 默认对称性检验（Sym），p < 0.05 被丢弃
phyloai pretree filter symtest --msa-dir ./trimmed

# Marginal 对称性（平稳性）检验
phyloai pretree filter symtest --msa-dir ./trimmed --symtest-type MAR

# 更严格的阈值
phyloai pretree filter symtest --msa-dir ./trimmed --symtest-pval 0.01

# 带树目录：保留匹配的基因树
phyloai pretree filter symtest \
  --msa-dir ./trimmed --tree-dir ./genetrees

# Internal 同质性检验
phyloai pretree filter symtest --msa-dir ./trimmed --symtest-type INT

# Dry-run 检查命令
phyloai pretree filter symtest --msa-dir ./trimmed --dry-run
```

### 备注

对称性检验应在比对、修剪之后、超矩阵拼接之前运行，因为平稳性或同质性的违反可能使系统发育推断产生偏差。`--symtest-type` 默认（合并的 Sym 检验）是最通用、应用最广的。

参考文献：Naser-Khdour 等（2019）"Assessing the Goodness of Fit of Phylogenetic Models..." doi:10.1093/gbe/evz193。

---

## `filter cluster` —— 基于聚类的探索

### 目的

使用降维（PCA 或 UMAP）后接层次聚类按度量特征对位点分组。这主要是一个探索性工具：默认仅写出聚类、诊断图与各聚类的度量汇总，不移除任何位点。

使用 `--drop-outlier-clusters auto` 可选择性地移除表现最差的聚类。`filter cluster` **不**应用基于规则的过滤（使用 `filter metrics`），也不屏蔽/剪除单个位点/分类单元（使用 `filter taper` 或 `filter treeshrink`）。

### 用法

```bash
phyloai pretree filter cluster \
  --table <metrics.csv|metrics.tsv> \
  [--input-format auto|csv|tsv] \
  [--metrics all|col1,col2,...] \
  [--exclude-regex REGEX] [--exclude-regex REGEX] \
  [--reduction pca|umap] \
  [--n-clusters N] [--max-clusters N] \
  [--cluster-linkage ward|average|complete|single] \
  [--cluster-distance euclidean|cosine|manhattan] \
  [--drop-outlier-clusters none|auto] \
  [--outlier-metric average_BS] [--outlier-direction low|high] \
  [--max-drop-fraction 0.2] \
  [--plot-metrics-cols N] [--plot-label-angle 45] \
  [--outlier-boxplot-cols N] \
  [--umap-n-neighbors 15] [--umap-min-dist 0.001] \
  [--umap-replicates 1] [--umap-random-state 42] \
  [--threads 1] \
  [--msa-dir <msa_dir>] [--tree-dir <tree_dir>] [--copy] \
  [--output-dir runs/pretree/filter/cluster] \
  [--dry-run] [--overwrite]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--table` | required | 度量 CSV/TSV。 |
| `--input-format` | `auto` | `csv`、`tsv` 或 `auto`。 |
| `--metrics` | all | 逗号分隔的度量列名。默认：除位点 ID、`DataType`、常数列和 `--exclude-regex` 匹配外的所有数值列。 |
| `--exclude-regex` | — | 可重复。`--exclude-regex '^freq' --exclude-regex '^sd_'` 或组合：`--exclude-regex '^(freq\|sd_)'`。 |
| `--reduction` | `pca` | `pca`（3 个分量，确定性，`scikit-learn`）或 `umap`（3 个分量，随机性，可选 `umap-learn`）。 |
| `--n-clusters` | auto | 固定聚类数。省略时通过多度量投票自动选择。 |
| `--max-clusters` | auto | 自动选择的上界。默认：`min(30, max(6, ceil(sqrt(n_loci)/3)))`。 |
| `--cluster-linkage` | `ward` | `ward`（最小化簇内方差），`average`、`complete`、`single`。 |
| `--cluster-distance` | `euclidean` | `euclidean`、`cosine`、`manhattan`。Ward 要求欧氏距离。 |
| `--drop-outlier-clusters` | `none` | `none`：仅诊断。`auto`：移除最差聚类。 |
| `--outlier-metric` | `average_BS` | 移除排序时使用的度量。 |
| `--outlier-direction` | `low` | `low`：更小值更差。`high`：更大值更差。 |
| `--max-drop-fraction` | 0.2 | 可被移除位点的最大比例（0.0–1.0）。 |
| `--plot-metrics-cols` | 2 | 每个图的箱线图列数（cluster_metric_boxplots）。行数自动计算。 |
| `--plot-label-angle` | 45.0 | 图中 X 轴标签的旋转角度。 |
| `--outlier-boxplot-cols` | 4 | 每个图的箱线图列数（outlier_comparison_boxplots）。行数自动计算。 |
| `--umap-n-neighbors` | 15 | UMAP 局部/全局平衡（PCA：忽略）。 |
| `--umap-min-dist` | 0.001 | UMAP 点分布紧密度（PCA：忽略）。 |
| `--umap-replicates` | 1 | UMAP 运行次数；通过聚类验证秩和评分选择最佳（PCA：忽略）。 |
| `--umap-random-state` | 42 | UMAP 随机种子。仅在 `--umap-replicates 1` 时应用；>1 时不设置种子以便 `--threads` 并行。 |
| `--threads` | 1 | UMAP n_jobs 的 CPU 线程数。仅在 `--reduction umap` 且 `--umap-replicates > 1` 时使用。 |
| `--msa-dir` / `--tree-dir` | — | `--copy` 模式的输入目录。 |
| `--output-dir` / `-o` | `runs/pretree/filter/cluster` | 输出目录。 |
| `--table-format` | `csv` | 所有输出表的格式。 |
| `--overwrite` | off | 删除并重建输出目录。 |
| `--dry-run` | off | 显示所选特征、降维、聚类范围、丢弃计划。 |
| `--quiet` / `-q` | off | 抑制终端输出。 |

### 输入

特征选择：输入表中的所有数值列，排除位点 ID、`DataType`、常数列与 `--exclude-regex` 匹配项。所有特征在降维前进行 z-score 标准化。

PCA 产出 `PC1`/`PC2`/`PC3`（通过 `sklearn.decomposition.PCA` 的 3 个分量）。UMAP 产出 `UMAP1`/`UMAP2`/`UMAP3`（需要 `pip install umap-learn`；缺失依赖会退出并提示安装）。

聚类数选择（未设置 `--n-clusters` 时）：评估 `k=2..max_clusters`。三个内部验证度量投票 —— silhouette（越高越好）、Calinski-Harabasz（越高越好）、Davies-Bouldin（越低越好）。平局时优先 silhouette，再优先更小的 `k`。UMAP 重复选择使用相同三个度量的秩和评分。

离群移除（`--drop-outlier-clusters auto` 时）：按 `--outlier-metric` 的均值排序聚类，丢弃最差的直到累计比例超过 `--max-drop-fraction`。

### 输出

核心（始终写入）：
```
runs/pretree/filter/cluster/
├── 01-input/
│   └── features_used.csv|tsv          （列、是否包含、原因）
├── 02-reduction/
│   ├── reduction.csv|tsv              （每个位点的 PC1/PC2/PC3 或 UMAP1/UMAP2/UMAP3）
│   ├── cluster_selection.csv|tsv      （k、silhouette、calinski_harabasz、davies_bouldin）
│   └── umap_replicates.csv|tsv        （仅当 UMAP replicates > 1 时）
├── 03-clustering/
│   ├── clusters.csv|tsv               （locus, cluster）
│   ├── cluster_summary.csv|tsv        （各聚类大小）
│   └── cluster_loci/cluster_*.csv|tsv  （每个聚类中的位点）
├── 04-diagnostics/
│   ├── cluster_metric_means.csv|tsv   （每个数值度量在每个聚类上的均值）
│   └── plots/
│       ├── cluster_2d.pdf             （前 2 个降维维度）
│       ├── cluster_3d.pdf             （3D 散点）
│       ├── cluster_metric_heatmap.pdf （z-score 热图：度量 × 聚类）
│       └── cluster_metric_boxplots.pdf （按聚类分组的各度量分布）
├── result.json
```

使用 `--drop-outlier-clusters auto` 时额外产出：
```
├── 05-outlier-drop/
│   ├── retained_loci.csv|tsv
│   ├── dropped_loci.csv|tsv
│   ├── filter_decisions.csv|tsv
│   ├── outlier_comparison.csv|tsv     （normal vs outlier：每个度量的均值、中位数、sd、计数）
│   ├── outlier_wilcoxon.csv|tsv       （各度量的 Mann-Whitney U p 值）
│   ├── plots/
│   │   └── outlier_comparison_boxplots.pdf  （所有度量，* p<0.05 ** p<0.01 *** p<0.001）
│   ├── seqs/                          （仅在使用 --copy --msa-dir 时）
│   └── trees/                         （仅在使用 --copy --tree-dir 时）
```

### 示例

```bash
# 探索性：查看位点如何按度量特征聚类
phyloai pretree filter cluster --table ./metrics/metrics.csv

# UMAP 固定聚类数
phyloai pretree filter cluster \
  --table ./metrics.csv --reduction umap --n-clusters 5

# 按 average_BS 丢弃离群聚类，复制存活文件
phyloai pretree filter cluster \
  --table ./metrics.csv \
  --drop-outlier-clusters auto \
  --outlier-metric average_BS \
  --max-drop-fraction 0.15 \
  --copy --msa-dir ./trimmed

# 排除频率与标准差列
phyloai pretree filter cluster \
  --table ./metrics.csv \
  --exclude-regex '^freq' --exclude-regex '^sd_'

# 指定度量子集
phyloai pretree filter cluster \
  --table ./metrics.csv \
  --metrics "average_BS,dvmc,gc_content,num_sites"
```

### 警告与错误

| 条件 | 行为 |
|------|------|
| `--table` 不存在 | Exit 1 |
| `--table` 为空或无数值列 | Exit 1 |
| `--reduction umap` 且未安装 `umap-learn` | Exit 1 并提示 `pip install umap-learn` |
| `--cluster-linkage ward` + 非欧氏 `--cluster-distance` | Exit 1 |
| `--n-clusters` < 2 或 > 位点数 | Exit 1 |
| `--copy` 没有 `--msa-dir` 或 `--tree-dir` | Exit 1 |
| `--drop-outlier-clusters auto` 与 `--copy` 但无位点被丢弃 | Copy 为 no-op（警告） |
| 输出目录非空且未加 `--overwrite` | Exit 1 |

### 备注

无 `--drop-outlier-clusters auto` 时，命令为只读 —— 不移除任何位点，仅写入诊断。这是刻意的：聚类解读应保持由用户引导；自动移除有潜在风险。

`features_used.csv` 是审计痕迹 —— 显示每一列、是否包含在特征集中，以及排除原因。

`plots/cluster_metric_heatmap.pdf` 是一个 z-score 标准化热图，其中：
- **行** = 聚类，**列** = 度量，**颜色** = 聚类均值与全局均值之差（以标准差计），高于为红色，低于为蓝色。
- 深红单元格表示该聚类在该度量上远高于均值；深蓝则远低于均值。
- 每个单元格用精确 z-score 值（如 `+1.53` 或 `−0.87`）标注。
- 用于快速识别：哪些聚类在哪些度量上得分高/低；离群聚类是否在 `average_BS` 等质量代理上系统性地偏低；以及聚类分离是否由少量主导度量驱动。

PCA 是默认降维方法，因为它是确定性的、稳定的，并且仅依赖 `scikit-learn`。UMAP 可用于探索非线性结构，但增加了可选的 `umap-learn` 依赖和随机性。`--umap-replicates 1` 时，固定 `--umap-random-state` 种子保证可重复性。`--umap-replicates > 1` 时不设种子，以便 `--threads` 并行化 UMAP。

终端输出按步骤展示进度（特征选择 → 降维 → 聚类 → 诊断 → 离群检测）并描述每个输出文件。`--drop-outlier-clusters auto` 之后，若提供 `--msa-dir`，屏幕还会显示保留 MSA 统计。

---

## result.json schema

所有子命令遵循相同的结果格式：

```json
{
  "status": "success | error",
  "command": "phyloai pretree filter ...",
  "wall_time": 1.23,
  "tool_versions": {},
  "params": {},
  "key_results": {},
  "error": null,
  "data": {}
}
```

各模式特定的 `key_results` 与 `data`：

| Subcommand | key_results | data |
|------------|-------------|------|
| `taper` | n_input, n_retained, n_dropped, masked_loci, total_masked_taxa, total_masked_aa_sites | retained_msa_stats, dry_run_cmds |
| `treeshrink` | n_input, n_retained, n_modified, n_dropped, n_removed_taxa_total | retained_loci, modified_loci, dropped_loci, removed_taxa, retained_msa_stats |
| `metrics` | n_total, n_retained, n_dropped, condition_failure_counts | copied_msa, copied_tree, retained_msa_stats |
| `cluster` | n_loci, n_valid_loci, n_features, n_clusters, reduction, selected_umap_replicate, n_retained, n_dropped | features, cluster_sizes, drop_clusters, retained_loci, retained_msa_stats, plot_paths, umap_replicates |

外部工具子命令（`taper`、`treeshrink`、`symtest`）的工具 stderr 写入每位点 `logs/<locus>.log`（taper）或共享的 `logs/<tool>.log` 文件（treeshrink、symtest）。纯 Python 子命令（`metrics`、`cluster`）不写外部日志。