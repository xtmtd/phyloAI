# phyloai tree ml fasttree

[English](tree-ml-fasttree.md) | [中文](tree-ml-fasttree.zh.md)

## 目的

使用 FastTree 推断最大似然系统发育树。

## 用法

```bash
# 从 MSA 目录批量推断基因树（并行）
phyloai tree ml fasttree --msa-dir ./trimmed/seqs \
    --seq-type AA --model lg --mode normal --boot 1000 \
    --cat 20 --gamma --threads 8 -o runs/tree/ml/fasttree

# 单一超矩阵树
phyloai tree ml fasttree --matrix ./concat/matrix.fa \
    --seq-type NT --model gtr --mode slow --boot 1000 \
    -o runs/tree/ml/fasttree

# 禁用 bootstrap（无节点支持）
phyloai tree ml fasttree --msa-dir ./trimmed --boot 0

# Fast 模式，JTT 模型（AA 默认）
phyloai tree ml fasttree --msa-dir ./trimmed --mode fastest --model jtt
```

## 输入

必须提供 `--msa-dir` 或 `--matrix` 之一，二者互斥。

## 示例

```bash
phyloai tree ml fasttree --msa-dir ./aligned --seq-type AA
```

## 参数

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--msa-dir` | — | MSA 文件目录。与 `--matrix` 互斥。 |
| `--matrix` | — | 单一拼接矩阵文件。与 `--msa-dir` 互斥。 |
| `--seq-type` | auto | AA、NT 或 auto（从输入检测）。 |
| `--model` | lg (AA) / gtr (NT) | 替换模型。AA：jtt、lg、wag。NT：jc、gtr。 |
| `--mode` | normal | 速度/精度：normal、fastest、slow。 |
| `--boot` | 1000 | Bootstrap 重复数。0 = 无支持（-nosupport）。 |
| `--cat` | 20 | 速率类别数。 |
| `--gamma` | on | Gamma 分布的速率异质性。默认始终启用；使用 --tool-args 禁用。 |
| `--output-dir` / `-o` | runs/tree/ml/fasttree | 输出目录。 |
| `--threads` / `-t` | 4 | 并行工作数（仅 `--msa-dir`）。 |
| `--fasttree-path` | — | 显式 FastTree 路径。 |
| `--tool-args` | — | FastTree 的额外策略标志。 |
| `--overwrite` | — | 覆盖已有输出目录。 |
| `--resume` | — | 从 checkpoint 恢复（仅 `--msa-dir`）。 |
| `--dry-run` | — | 打印命令而不执行。 |
| `--quiet` / `-q` | — | 除错误外不打印输出。 |

## 输出

- `result.json`：结构化结果（树、失败、跳过）
- `trees/`：Newick 树文件（每个输入一个，`--msa-dir` 模式）
- `logs/`：每基因 FastTree 日志（`--msa-dir` 模式）
- `checkpoint.json`：恢复状态（`--msa-dir` 模式）
- 单个 `.tre` 文件（`--matrix` 模式）

## 支持的格式

FastTree 原生读取 FASTA（`.fa`、`.fas`、`.fasta`、`.faa`、`.fna`）和 phylip-relaxed（`.phy`、`.phylip`）。

不支持 NEXUS 文件（`.nex`、`.nxs`、`.nexus`）。请先转换：
```bash
phyloai pretree convert --input data.nex --to fasta
```

## 警告与错误

| 条件 | 行为 |
|------|------|
| `--msa-dir` 与 `--matrix` 都提供或都不提供 | 错误：恰好需要一个 |
| `--overwrite` 与 `--resume` 同时使用 | 错误：互斥 |
| `--resume` 在 `--matrix` 模式下 | 错误：仅批量模式支持 resume |
| `--threads` 在 `--matrix` 模式下 | 警告：`--threads` 在单文件模式下无效 |
| `--msa-dir` 不存在 | 错误：未找到目录 |
| `--msa-dir` 中无有效输入 | 错误：无有效输入文件 |

## 备注

- `--threads` 仅控制 `--msa-dir` 模式下的并行基因树推断。FastTree 本身是单线程的。
- `--resume` 仅在 `--msa-dir` 批量模式下可用。
- 模型默认：AA 为 LG，NT 为 GTR。