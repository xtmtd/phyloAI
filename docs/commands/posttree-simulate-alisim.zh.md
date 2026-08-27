# phyloai posttree simulate alisim

[English](posttree-simulate-alisim.md) | [中文](posttree-simulate-alisim.zh.md)

## 目的（Purpose）

基于 IQ-TREE3 AliSim 的三步序列模拟工作流，用于**保留真实数据集的实证特征**，而不是使用任意模型参数：

1. **params** —— 从已有的 IQ-TREE 报告中按位点提取模拟参数（模型字符串组件、频率、不变位点比例、速率异质性、树），生成 `params.tsv` 表。
2. **iqtree** —— 用 AliSim 模拟比对：既可从显式输入执行单次 IQ-TREE 调用，也可从 `params.tsv` 中采样行进行可恢复的批量模拟。
3. **transfergaps** —— 将原始的按类群 gap 掩码重新引入模拟比对（AliSim 生成无 gap 的 MSA），可处理单条（`--simulated-msa`）或一个目录下的批量（`--simulated-dir`）。

## 用法（Usage）

```bash
# 1. 从之前的 IQ-TREE 运行中提取参数
phyloai posttree simulate alisim params --iqtree-dir runs/tree/ml/iqtree --tree-dir runs/tree/ml/iqtree

# 2. 从提取的表中进行批量模拟
phyloai posttree simulate alisim iqtree --model-params params.tsv --strategy complete --num-simulations 100

# 3. 将原始 gap 模式转移到一条模拟比对
phyloai posttree simulate alisim transfergaps --original-msa original.fa --simulated-msa sim001.fa

# 4. ...或将原始 gap 模式转移到目录下的每一条模拟比对
phyloai posttree simulate alisim transfergaps --original-msa original.fa --simulated-dir MSAs/
```

## alisim params

每个成功解析的 `.iqtree` 报告提取一行参数，并按逻辑位点名与树文件配对。

### 参数

| 选项 | 说明 |
|--------|-------------|
| `--iqtree-dir` | 存放 `.iqtree` 报告文件的目录（任意嵌套深度，以 `**/*.iqtree` 通配）。必填。 |
| `--tree-dir` | 存放树文件的目录。按逻辑位点名（忽略后缀）与 `.iqtree` 文件匹配。必填。 |

### 输出

- `params.tsv` —— 列：`id, seqtype, length, subs_model, subs_rate, freq, prop_inv, rate_heterogeneity, rate_categories, rate_param, tree_path`。多值列（`subs_rate`、`freq`、`rate_param`）使用 IQ-TREE 的逗号分隔符，与 `.iqtree` 报告一致。
- `result.json` —— `n_loci_parsed`、`n_loci_matched`、`n_loci_unmatched`、`seq_types`，以及未匹配位点列表。

报告无匹配树的位点会被标记为 unmatched 并跳过。多个树文件匹配同一逻辑位点名（歧义）是硬错误。

## alisim iqtree

两种互斥模式。

### 单模式（一次 IQ-TREE 调用）

提供参考树、模型字符串或分区文件、序列类型和比对长度：

```bash
phyloai posttree simulate alisim iqtree --ref-tree ref.nwk --model LG+G4 --seq-type AA --length 2000
phyloai posttree simulate alisim iqtree --ref-tree ref.nwk --model-partitions matrix.best_model.nex --seq-type AA
```

#### 参数

| 选项 | 说明 |
|--------|-------------|
| `--ref-tree` | 参考树（Newick）。对应 IQ-TREE `-t`。必填。 |
| `--model` | IQ-TREE 模型字符串（如 `GTR{XXX}+F{XXX}+G4{XXX}`）。对应 `-m`。与 `--model-partitions` 互斥。 |
| `--model-partitions` | NEXUS 分区模型文件。对应 `-p`。与 `--model` 互斥。使用时 `--length` 由分区定义推断，必须省略。 |
| `--seq-type` | `AA` 或 `DNA`（不区分大小写）。对应 IQ-TREE `--seqtype`。必填。 |
| `--length` | 比对长度。对应 `--length`。使用 `--model-partitions` 时可不填。 |
| `--msa-prefix` | 输出 MSA 文件前缀（默认 `sim`）。 |
| `--num-alignments` | 每次 IQ-TREE 调用生成的 MSA 数量（默认 1）。仅单模式。 |
| `--out-format` | 输出 MSA 格式：`fasta`（默认）或 `phy`。对应 `--out-format`。 |
| `--iqtree-threads` | 每次 IQ-TREE 调用的线程数（默认 1）。对应 `-T`。 |
| `--seed` | 随机种子。对应 `--seed`。 |

#### 输出

- `MSAs/<prefix>.*` —— 模拟 MSA（扩展名随 `--out-format`）。
- `logs/<prefix>.log` —— 捕获的 IQ-TREE 控制台输出。
- `result.json` —— `n_msas_generated`、执行的 IQ-TREE 命令及输出文件路径。

AliSim 不产生 `.iqtree` 报告。

### 批量模式（可恢复，每个 MSA 一次 AliSim 调用）

提供 `--model-params` 表（来自 `alisim params`）以及采样策略和目标数量：

```bash
phyloai posttree simulate alisim iqtree --model-params params.tsv --strategy pdf --num-simulations 100 --seed 42
```

#### 参数

| 选项 | 说明 |
|--------|-------------|
| `--model-params` | 来自 `alisim params` 的 TSV 表。激活批量模式。必填。 |
| `--strategy` | `complete`（默认）、`mixed` 或 `pdf`。见[采样策略](#采样策略)。 |
| `--num-simulations` | 要模拟的 MSA 总数。必填。 |
| `--override` | 逗号分隔的 `key=value`，固定应用于所有模拟，如 `length=500,prop_inv=0.1`。对全部策略生效。有效键：`length`、`prop_inv`。 |
| `--noise-scale` | PDF 重采样噪声（0.0-1.0，默认 1.0）：0 = 箱中心，1 = 箱内完全均匀抖动。需 `--strategy pdf`。 |
| `--pdf-params` | 通过密度估计重采样的参数列表（默认 `length,prop_inv,rate_param`）。有效：`length`、`prop_inv`、`rate_param`。需 `--strategy pdf`。 |
| `--seed` | 主种子；每个模拟从主种子生成的独立随机数发生器各取一个随机种子。 |
| `-t, --threads` | 并行模拟任务数（默认 4）。 |
| `--out-format` | 每个模拟的输出 MSA 格式：`fasta`（默认）或 `phy`。 |
| `--iqtree-threads` | 每次 IQ-TREE 调用的线程数（默认 1）。对应 `-T`。 |

#### 采样策略

- **complete** —— 每个模拟比对复刻单个源基因模型的完整参数集（模型核心、速率异质性、比对长度、不变位点比例和树都取自同一行）。
- **mixed** —— 模型核心组、速率异质性组、比对长度、不变位点比例和参考树各自从实证基因模型分布中独立采样，保留单个参数的实证分布及其存在/缺失比例。
- **pdf**（概率密度函数，probability density function）—— 在 mixed 基础上构建；`--pdf-params` 中的参数通过直方图密度估计（Freedman-Diaconis 箱、`--noise-scale` 噪声）对实证概率密度进行重采样，而非直接从实证列抽取。

采样说明：
- 仅 `""` 视为 `prop_inv` 缺失；非空值如 `"0"` 会被保留并重建为 `+I{0}`。
- `rate_param` 仅当采样到的速率组为 `G`（Gamma）时才做密度重采样；FreeRate（`R`）组始终经验采样。
- `--override` 对所有策略生效，同时覆盖表中数值与密度重采样。
- PDF 密度图（经验 vs 模拟的 Gaussian-KDE 曲线，`server.R` 配色 `#2E86AB`/`#A23B72`，仅曲线）仅在 `--strategy pdf` 下生成；`complete`/`mixed` 运行不产生 `plots/` 目录。

#### 输出

- `MSAs/sim001.<ext>, ...` —— 每个模拟一个 MSA。
- `logs/<simulation_id>.log` —— 每个模拟捕获的 IQ-TREE 控制台输出。
- `params_sampled.tsv` —— 实际使用的每一行。`source_id` 列仅 `complete` 策略存在；每行 `seed` 独立随机。
- `plots/*_density.pdf` —— 仅 pdf 策略。
- `checkpoint.json`、`result.json` —— `result.json` 含 `source_loci`、`n_simulations_completed`、`n_simulations_failed`；`complete`/`mixed` 时 `noise_scale` 与 `pdf_params` 为 `null`。

恢复规则：
- `--resume` 从 `checkpoint.json` 恢复批量运行；已完成的模拟被跳过，未完成的被重试。仅批量模式。
- `--overwrite` 与 `--resume` 互斥。
- 非空输出目录需要 `--overwrite`（或 `--resume`）。
- 恢复要求与原始运行相同的参数。

## alisim transfergaps

`--simulated-msa`（单条模式）与 `--simulated-dir`（批量模式）互斥，且必须恰好提供其中一个。

### 参数

| 选项 | 说明 |
|--------|-------------|
| `--original-msa` | 单个原始（带 gap）MSA 文件。必填。 |
| `--simulated-msa` | 来自 `alisim iqtree` 的单个模拟（无 gap）MSA 文件。与 `--simulated-dir` 互斥。 |
| `--simulated-dir` | 来自 `alisim iqtree` 的模拟（无 gap）MSA 文件目录（仅比对格式扩展名）。每个输入输出一个文件，命名为 `<stem>.gaps.fa`。与 `--simulated-msa` 互斥。 |
| `--seq-type` | `AA`、`NT` 或 `auto`（默认，不区分大小写）。决定有效字符集。 |
| `--exclude-ambiguity` | 设置后仅转移真正的 gap 字符（`-`、`.`）；模糊码保留为模拟字符。默认会遮蔽标准字母表以外的每个字符（包括模糊码）。 |

### 输出

- 单条模式：`<original_stem>.gaps.fa`。
- 批量模式：每个输入一个 `<simulated_stem>.gaps.fa`。
- `result.json` —— `n_sequences`、`alignment_length`、`n_positions_masked`、`mean_positions_masked_per_taxon`、`detected_seq_type`、`n_msas`，以及 `data.output_files` 下的转移文件列表。

无论输入格式如何，输出恒为 60 列 FASTA。

## 通用标志

三个命令均支持 `-o, --output-dir`（默认 `runs/posttree/simulate/alisim/<sub>`）、`--overwrite`、`--dry-run`、`-q, --quiet`。`dry-run` 只做校验并打印计划，不写任何文件。

## 警告与错误（Warnings & Errors）

| 情况 | 行为 |
|-----------|----------|
| 输出目录已存在且非空，且未给 `--overwrite`/`--resume` | 硬错误，指明目录 |
| `params`：某报告无匹配树 | 警告；位点列入 `data.unmatched` 并跳过 |
| `params`：多个树匹配同一报告 | 硬错误，列出歧义文件 |
| `transfergaps`：类群集合或长度不匹配、重复类群、无法解析或空输入 | 硬错误（批量模式指明失败文件） |
| `transfergaps`：模拟长度 ≠ 原始长度 | 硬错误，说明 AliSim `--length` 必须等于原始列数 |
| `iqtree`：`--tool-args` 含受管 I/O 标志 | 硬错误，列出该标志 |
| `iqtree`：IQ-TREE 非零退出 | 退出码 2，附 stderr |

## 退出码

| 代码 | 含义 |
|------|---------|
| 0 | 成功 |
| 1 | 用户输入错误（文件缺失、参数无效、输出冲突） |
| 2 | IQ-TREE 执行失败 |
| 3 | 未找到 IQ-TREE 可执行文件 / 环境错误 |

## 备注

- 在 `transfergaps` 后运行 [posttree simulate adequacy](posttree-simulate-adequacy.zh.md)，比较观测与模拟 MSA 的 PPA-DIV、PPA-CONV、PPA-VAR 和 PPA-COMP 统计量。`phybase` 仍为未来占位命令。
- `alisim iqtree` 的 `--tool-args` 透传额外的 IQ-TREE 标志；仅受管的 I/O 标志（`--alisim`、`-t`、`--prefix`、`--out-format`、`-af`，`--flag` 或 `--flag=value` 两种形式均被拦截）被阻止。其它标志（如 `--seqtype`、`--length`、`--num-alignments`、`-T`）可覆盖 PhyloAI 默认值；PhyloAI 会抑制自身已覆盖的受管标志，最终 IQ-TREE 命令中每个标志只出现一次。
- AliSim 不产生 `.iqtree` 报告；`logs/` 存放捕获的 IQ-TREE 控制台输出。
- IQ-TREE3（`iqtree3`）必须在 `PATH` 上（或通过 `--iqtree-path` 指定）；`phyloai doctor` 报告其检测状态。
- 在系统误差工作流中，使用 `readpb --mode ss,rr,r` 输入的 AliSim 是后验均值
  plug-in 参数模拟，不是严格的后验预测模拟。主要的 `readpb --mode ppred` 路径以及
  无 gap/gap-transferred 比较见
  [系统误差工作流参考](../../skills/phyloai-workflow/references/syserror-workflow.md)。
