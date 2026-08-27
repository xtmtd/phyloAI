# phyloai posttree syserror taxcomp

[English](posttree-syserror-taxcomp.md) | [中文](posttree-syserror-taxcomp.zh.md)

## 用途

在单个核苷酸或氨基酸比对中筛查**跨类群（taxa）**的组成异质性。`taxcomp` 构建一个"类群 × 状态"计数表，并从同一张表提供两个互补视图：

1. **Pearson 共同组成卡方检验**：一个总体同质性统计量（P4 式 omnibus），加上每个类群的行贡献，以及名义 per-taxon p 值（即 IQ-TREE/P4 式的 per-taxon 筛查值）和 Holm 校正后的探索性 p 值。
2. **`phyloai posttree simulate adequacy` 已实现的观测 PPA-COMP 描述统计**：每个类群的平方组成距离、`comp_max` 和 `comp_mean`。

`taxcomp` 是纯 Python 诊断命令。它不使用树、替换模型或外部可执行程序，因此不需要 `phyloai doctor`、断点续跑（checkpoint）或 resume 支持。

**非目标。** 本命令不删除类群、不重编码数据、不选择模型或拓扑、不推荐阈值、不进行后验预测模拟。其 p 值只是名义上的探索性结果，不是经系统发育校准的模型充分性检验。低 p 值或大组成距离是检查的信号，而非结论。

## 用法

```bash
phyloai posttree syserror taxcomp \
  --matrix matrix.aa.fa \
  [--seq-type AA|NT|auto] \
  [--table-format csv|tsv] \
  [-o runs/posttree/syserror/taxcomp] \
  [--overwrite] [--dry-run] [-q]
```

## 示例

对氨基酸比对进行显式序列类型筛查：

```bash
phyloai posttree syserror taxcomp --matrix matrix.aa.fa --seq-type AA
```

对核苷酸比对输出 TSV 汇总并写入自定义目录：

```bash
phyloai posttree syserror taxcomp --matrix matrix.nt.fa --seq-type NT --table-format tsv -o runs/posttree/syserror/taxcomp-nt
```

## 输入

| 选项 | 必需 | 默认 | 说明 |
|---|---|---|---|
| `--matrix` | 是 | -- | 一个已比对的 FASTA、PHYLIP、PHYLIP-PAML 或 Nexus 比对。格式用现有 PhyloAI 读取器自动检测；不支持 Clustal。不需要树或模型输入。 |
| `--seq-type` | 否 | `auto` | `AA`、`NT` 或自动检测。AA 只统计标准氨基酸，NT 只统计 `ACGT`；其余所有字符（gap、未知、简并码、终止符）都视为缺失，不产生小数计数。 |
| `--table-format` | 否 | `csv` | 两个汇总表的定界符和后缀（`csv` 或 `tsv`）。 |
| `-o`, `--output-dir` | 否 | `runs/posttree/syserror/taxcomp` | 输出目录。应用标准的非空目录冲突策略。 |
| `--overwrite` | 否 | false | 仅在验证成功后删除并重建非空输出目录。 |
| `--dry-run` | 否 | false | 校验、解析并计算所有汇总，但不写入任何文件。 |
| `-q`, `--quiet` | 否 | false | 除错误外抑制终端输出。 |

比对必须包含至少两个名称唯一的类群和至少两个全局观测到的标准状态；所有序列必须比对到相同长度。重复类群标识符、空或不可读的比对、零有效字符的类群、序列长度不等都是硬错误。

## 输出

```
<output_dir>/
├── overall_summary.csv|tsv
├── taxon_summary.csv|tsv
└── result.json
```

### overall_summary

一行，列为：

| 列 | 含义 |
|---|---|
| `n_taxa` | 类群数 |
| `n_states` | 有效 Pearson `K`：进入期望计数和自由度的全局观测标准状态数 |
| `x2` | 总体 Pearson 卡方（各类群行贡献之和） |
| `df` | `(n_taxa - 1) * (n_states - 1)` |
| `p_nominal` | 总体检验的名义卡方生存概率 |
| `sparse_count_check` | `triggered` 或 `not_triggered` |
| `expected_cells_total` | 期望单元格总数（`n_taxa * n_states`） |
| `expected_cells_below_1` | `< 1` 的期望单元格数 |
| `expected_cells_below_5` | `< 5` 的期望单元格数 |
| `expected_cells_below_5_fraction` | `< 5` 的期望单元格占比 |
| `comp_max` | 每个类群平方组成距离的最大值 |
| `comp_mean` | 每个类群平方组成距离的均值 |

### taxon_summary

每类群一行，按输入顺序：

| 列 | 含义 |
|---|---|
| `taxon` | 类群标识符，原样保留 |
| `x2_contribution` | 该类群对总体 X2 的行贡献 |
| `df` | `n_states - 1`（per-taxon 筛查自由度） |
| `p_nominal` | 名义 per-taxon 卡方 p 值（IQ-TREE/P4 式筛查值） |
| `p_holm` | Holm 逐步校正的名义 p 值，还原到输入顺序 |
| `squared_composition_distance` | 与等类群平均组成的无量纲平方欧氏频率偏差；即每个类群的观测 PPA-COMP 值（即 PhyloBayes `chain1.comp` 的 `obs comp`，以及 `posttree simulate adequacy` 的 `adequacy_taxon_comp.csv` 中 `obs` 列） |

## 解读

- **总体 X2**：在常规卡方筛查下，反对"所有类群共享一个组成"的合并证据。它是 PhyloAI/P4 式扩展，不是 IQ-TREE 输出。
- **类群 X2 贡献**：每行对总体 X2 的贡献。其 p 值是筛查值，**不是**独立的"一类群对剩余类群"列联检验。
- **`p_holm`**：对名义 per-taxon p 值应用 Holm 多重性校正。当边际 p 值有效时，Holm 在任意依赖下控制族系错误率，但这里它继承了名义卡方校准的局限。
- **`sparse_count_check`**：报告常规稀疏单元格规则是否触发（任一期望单元格 `< 1`，或超过 20% 的期望单元格 `< 5`）。`not_triggered` **不**意味着独立、**不**验证系统发育零模型，也**不**把筛查 p 值变成后验预测检验；`triggered` 则警告名义 p 值尤其不可靠。
- **平方组成距离**：与等类群平均组成的无量纲平方欧氏频率偏差。它**不是**进化距离，也**没有**通用阈值；`comp_max` 和 `comp_mean` 描述最大和平均的观测偏差。不同字母表、类群集合和缺失模式之间的数值未必可直接比较。

即使稀疏单元格规则未触发，系统发育依赖也限制了常规卡方 p 值的解释。报告模板会重复这些注意事项，绝不会把类群归类为 significant、failed、outlier、biased 或 removable。

## 敏感性后续分析

大距离或低 p 值是检查注释、覆盖度、污染、谱系组成和模型充分性的提示。重编码是独立、需用户批准、使用现有 concat 流程的敏感性分析，必须与原始分析并列报告：

- AA：`phyloai pretree concat --recoding Dayhoff-6`
- NT：`phyloai pretree concat --recoding RY-nucleotide`

删除类群是人工数据整理决策，需要独立证据，并在适用时与经模型校准的 `phyloai posttree simulate adequacy` 结果对照。本命令从不产生删除清单。

## 警告 / 错误

- 非空输出目录在未给 `--overwrite` 时被拒绝。
- 目录被占用（claim）之前发生的验证失败不写任何文件。
- `--overwrite` 仅在验证成功后删除并重建目录；出错时保留现有文件，可能只替换根 `result.json`。
- `--dry-run` 不写任何文件，并把计算出的汇总以 JSON 打印。
- `result.json` 遵循标准 PhyloAI schema：成功时 `error_category: null`，出错时用标准类别。

## 备注

- `taxcomp` 只是筛查：它从不删除类群、不重编码数据、不选模型或拓扑、不推荐阈值。
- 所有 p 值都是名义/探索性的；即使稀疏单元格规则未触发，系统发育依赖也限制常规卡方 p 值的解释。
- `squared_composition_distance` 即每个类群的观测 PPA-COMP 值（`obs comp`）；经模型校准的解释需用
  `posttree simulate adequacy`，它提供基于模拟的 z 分数与后验预测 pp 值。
- 不同字母表、类群集合和缺失模式之间的数值未必可直接比较。
- 关于理论、证据边界、重编码/类群采样敏感性选择，以及与非平稳组成模型的界限，见
  [系统误差工作流参考](../../skills/phyloai-workflow/references/syserror-workflow.md)。
