# phyloai posttree simulate alisim

[English](posttree-simulate-alisim.md) | [中文](posttree-simulate-alisim.zh.md)

基于 IQ-TREE3 AliSim 的三步序列模拟工作流，用于**保留真实数据集的实证特征**，而不是使用任意模型参数。

## 工作流

1. **params** —— 从已有的 IQ-TREE 报告中按位点提取模拟参数（模型字符串组件、频率、不变位点比例、速率异质性、树），生成 `params.tsv` 表。
2. **iqtree** —— 用 AliSim 模拟比对：既可从显式输入执行单次 IQ-TREE 调用，也可从 `params.tsv` 中采样行进行可恢复的批量模拟。
3. **transfergaps** —— 将原始的按类群 gap 掩码重新引入模拟比对（AliSim 生成无 gap 的 MSA）。

```bash
# 1. 从之前的 IQ-TREE 运行中提取参数
phyloai posttree simulate alisim params --iqtree-dir runs/tree/ml/iqtree --tree-dir runs/tree/ml/iqtree

# 2. 从提取的表中进行批量模拟
phyloai posttree simulate alisim iqtree --model-params params.tsv --strategy complete --num-simulations 100

# 3. 将原始 gap 模式转移到一条模拟比对
phyloai posttree simulate alisim transfergaps --original-msa original.fa --simulated-msa sim001.fa
```

## `alisim params`

| 选项 | 说明 |
|--------|-------------|
| `--iqtree-dir` | 存放 `.iqtree` 报告文件的目录（任意嵌套深度，以 `**/*.iqtree` 通配）。必填。 |
| `--tree-dir` | 存放树文件的目录。按逻辑位点名（忽略后缀）与 `.iqtree` 文件匹配。必填。 |

输出 `params.tsv`（列：`id, seqtype, length, subs_model, subs_rate, freq, prop_inv, rate_heterogeneity, rate_categories, rate_param, tree_path`）以及 `result.json`（含 `n_loci_parsed`、`n_loci_matched`、`n_loci_unmatched`、`seq_types`）。报告无匹配树的位点会被标记为 unmatched 并跳过。

## `alisim iqtree`

两种互斥模式。

### 单模式（一次 IQ-TREE 调用）

提供参考树、模型字符串或分区文件、序列类型和比对长度：

```bash
phyloai posttree simulate alisim iqtree --ref-tree ref.nwk --model LG+G4 --seq-type AA --length 2000
phyloai posttree simulate alisim iqtree --ref-tree ref.nwk --model-partitions matrix.best_model.nex --length 2000
```

| 选项 | 说明 |
|--------|-------------|
| `--ref-tree` | 参考树（Newick）。对应 IQ-TREE `-t`。 |
| `--model` | IQ-TREE 模型字符串（如 `GTR{XXX}+F{XXX}+G4{XXX}`）。对应 `-m`。与 `--model-partitions` 互斥。 |
| `--model-partitions` | NEXUS 分区模型文件。对应 `-p`。与 `--model` 互斥。 |
| `--seq-type` | `AA` 或 `DNA`。对应 `--seqtype`。 |
| `--length` | 比对长度。对应 `--length`。 |
| `--num-alignments` | 每次 IQ-TREE 调用生成的 MSA 数量（默认 1）。仅单模式。 |
| `--msa-prefix` | 输出 MSA 文件前缀（默认 `sim`）。 |
| `--out-format` | 输出 MSA 格式：`fasta`（默认）或 `phy`。对应 `--out-format`。 |
| `--iqtree-threads` | 每次 IQ-TREE 调用的线程数（默认 1）。对应 `-T`。 |
| `--seed` | 随机种子。对应 `--seed`。 |

输出：`MSAs/<prefix>.*`、`logs/<prefix>.iqtree`、`logs/<prefix>.log`、`result.json`。

### 批量模式（可恢复，每个 MSA 一次 AliSim 调用）

提供 `--model-params` 表（来自 `alisim params`）以及采样策略和目标数量：

```bash
phyloai posttree simulate alisim iqtree --model-params params.tsv --strategy pdf --num-simulations 100 --seed 42
```

| 选项 | 说明 |
|--------|-------------|
| `--model-params` | 来自 `alisim params` 的 TSV 表。激活批量模式。 |
| `--strategy` | `complete`（均匀采样完整行）、`mixed`（随机化模型类别 + 长度）或 `pdf`（直方图密度重采样）。批量模式必填。 |
| `--num-simulations` | 要模拟的 MSA 总数。批量模式必填。 |
| `--override` | 逗号分隔的 `key=value`，固定应用于所有模拟，如 `length=500,prop_inv=0.1`。有效键：`length`、`prop_inv`。 |
| `--noise-scale` | PDF 重采样噪声（0.0-1.0，默认 1.0）：0 = 箱中心，1 = 箱内完全均匀抖动。需 `--strategy pdf`。 |
| `--pdf-params` | 通过密度重采样采样的参数列表（默认 `length,prop_inv,rate_param`）。有效：`length`、`prop_inv`、`rate_param`。需 `--strategy pdf`。 |
| `--seed` | 主种子；每个任务的种子 = 主种子 + 任务序号。 |
| `-t, --threads` | 并行模拟任务数（默认 4）。 |

采样说明：
- 完整行作为原子单元采样：模型核心与速率组保持一致，`+I` 的存在性在重采样其值之前决定。
- 仅 `""` 视为 `prop_inv` 缺失；非空值如 `"0"` 会被保留并重建为 `+I{0}`。
- PDF 密度图仅为选定的 PDF 参数（未被 override）生成；`length` 使用 Freedman-Diaconis 箱，`prop_inv` 限制在 `[0, 1)`。

输出：`MSAs/sim001.fa, ...`、`logs/<simulation_id>.log`、`params_sampled.tsv`（实际使用的每一行）、`plots/*_density.pdf`（pdf 策略）、`checkpoint.json`、`result.json`（含 `source_loci`、`n_simulations_completed`、`n_simulations_failed`）。

恢复规则：
- `--resume` 从 `checkpoint.json` 恢复批量运行；已完成的模拟被跳过，未完成的被重试。仅批量模式。
- `--overwrite` 与 `--resume` 互斥。
- 非空输出目录需要 `--overwrite`（或 `--resume`）。
- 恢复要求与原始运行相同的参数。

## `alisim transfergaps`

| 选项 | 说明 |
|--------|-------------|
| `--original-msa` | 单个原始（带 gap）MSA 文件。必填。 |
| `--simulated-msa` | 来自 `alisim iqtree` 的单个模拟（无 gap）MSA 文件。必填。 |
| `--seq-type` | `AA`、`NT` 或 `auto`（默认）。决定有效字符集。 |
| `--exclude-ambiguity` | 设置后仅转移真正的 gap 字符（`-`、`.`）；模糊码保留为模拟字符。默认会遮蔽标准字母表以外的每个字符（包括模糊码）。 |

校验：输入可解析且非空、无重复类群 ID、类群集合必须匹配、原始/模拟比对长度必须相等。掩码位点被替换（从不插入）；输出顺序遵循原始 MSA。

输出：`<original_stem>_transferred.<ext>`（60 列 FASTA）和 `result.json`（含 `n_sequences`、`alignment_length`、`n_positions_masked`、`mean_positions_masked_per_taxon`、`detected_seq_type`）。

## 通用标志

三个命令均支持 `-o, --output-dir`（默认 `runs/posttree/simulate/alisim/<sub>`）、`--overwrite`、`--dry-run`、`-q, --quiet`。`dry-run` 只做校验并打印计划，不写任何文件。

## 退出码

| 代码 | 含义 |
|------|---------|
| 0 | 成功 |
| 1 | 用户输入错误（文件缺失、参数无效、输出冲突） |
| 2 | IQ-TREE 执行失败 |
| 3 | 未找到 IQ-TREE 可执行文件 / 环境错误 |

## 备注

- `adequacy` 与 `phybase` 子命令保留给未来工作，目前返回未实现消息。
- `alisim iqtree` 的 `--tool-args` 透传额外的 IQ-TREE 标志；受管标志（`--alisim`、`-t`、`-m`、`-p`、`-q`、`-Q`、`--seqtype`、`--length`、`--out-format`、`-af`、`--num-alignments`、`-T`、`--seed`、`--prefix`）被阻止。
- IQ-TREE3（`iqtree3`）必须在 `PATH` 上（或通过 `--iqtree-path` 指定）；`phyloai doctor` 报告其检测状态。
