# phyloai posttree syserror rate

[English](posttree-syserror-rate.md) | [中文](posttree-syserror-rate.zh.md)

## 目的

对来自 IQ-TREE3 或 PhyloBayes 的逐位点速率进行排序，并可从已有比对中提取最慢或
最快的位点比例。这是位点速率排序/提取的敏感性分析工具：应比较不同保留比例的下游
分析；它不证明某个拓扑、不自动校正系统误差，也不自动选择比例。

该命令为本地纯 Python 工具，不调用外部可执行程序，因此无需 `phyloai doctor` 检查，
也不支持 resume。

## 用法

```bash
phyloai posttree syserror rate \
  (--iqtree-rate matrix.rate | --pb-rate chain.meansiterates) \
  [--matrix raw.fa] [--subset slow|fast] [--fraction 0.25,0.5,0.75] [options]
```

## 输入

| 选项 | 必填 | 默认值 | 说明 |
|---|---|---|---|
| `--iqtree-rate` | 二选一 | -- | 含 `Site`、`Rate` 列的 IQ-TREE3 `--rate` 表。 |
| `--pb-rate` | 二选一 | -- | PhyloBayes `readpb -r` 的无表头 `<site> <rate>` 表。 |
| `--matrix` | 否 | -- | 原始 MSA，提供后才能提取子集。支持 FASTA、relaxed PHYLIP、PAML PHYLIP、NEXUS。 |
| `--subset` | 有 matrix 时 | `slow` | 提取 `slow` 或 `fast` 位点；无 `--matrix` 时无效。 |
| `--fraction` | 有 matrix 时 | -- | 一个或多个 `(0, 1]` 的逗号分隔比例；有 matrix 时必填，无 matrix 时无效。 |
| `-o`, `--output-dir` | 否 | `runs/posttree/syserror/rate` | 输出目录。 |
| `--overwrite` | 否 | false | 删除并重建非空输出目录。 |
| `--dry-run` | 否 | false | 校验输入、打印已校验的 payload，但不写文件。 |
| `-q`, `--quiet` | 否 | false | 除错误外抑制终端输出。 |

两个速率来源必须且只能提供一个。只有提供 `--matrix` 时，未指定的选择方向才默认
为 `slow`。

## 速率输入与索引

IQ-TREE 位点标识必须是严格连续的 1 基索引。PhyloBayes 位点标识必须是严格连续的
0 基索引，并会加一归一化。因此两种来源都会产生严格的 `1..N` 归一化位点。空输入、
格式错误、重复、非整数、非有限值、负值或不连续索引都会被拒绝。

`rates.csv` 使用 1 基归一化位点，按 `(rate, site)` 从慢到快确定性排序：

```csv
site,rate
21,0.19145
```

## 输出

未提供 `--matrix` 时：

```text
runs/posttree/syserror/rate/
├── rates.csv
└── result.json
```

提供 `--matrix --subset slow --fraction 0.25,0.5` 时：

```text
runs/posttree/syserror/rate/
├── rates.csv
├── slow25/
│   ├── positions.txt
│   └── matrix.fa
├── slow50/
│   ├── positions.txt
│   └── matrix.fa
└── result.json
```

每个比例保留 `ceil(N * fraction)` 个位点，不会因边界速率并列而扩展。`positions.txt`
每行一个 1 基原始位点位置，按原始比对顺序排列。生成的子比对为 FASTA，每行序列
长度为 60。`result.json` 记录输入、解析后的设置、汇总统计、子集数量和所有生成文件。

## 示例

```bash
# 仅对 IQ-TREE 速率排序
phyloai posttree syserror rate --iqtree-rate matrix.rate

# 慢位点敏感性分析
phyloai posttree syserror rate --iqtree-rate matrix.rate --matrix raw.fa \
  --subset slow --fraction 0.25,0.5,0.75 -o runs/posttree/syserror/rate

# 从 PhyloBayes 速率提取快位点
phyloai posttree syserror rate --pb-rate chain.meansiterates --matrix raw.phy \
  --subset fast --fraction 0.1
```

## 警告 / 错误

- `--iqtree-rate` 与 `--pb-rate` 必须且只能提供一个。
- 速率来源路径必须存在，且为可读取的普通文件。
- `--subset`、`--fraction` 需要 `--matrix`；`--matrix` 需要 `--fraction`。
- MSA 必须可解析，包含非空且 ID 唯一的记录，并且所有序列长度相同；其长度必须等于归一化后的速率数。
- 比例必须有效、唯一，并产生唯一的目录标签。
- 非空输出目录需要 `--overwrite`；不提供 resume/checkpoint。

## 说明

- 建议使用多个可解释比例（常用 `0.25,0.5,0.75`）做敏感性分析，而不要将任意一个
  阈值视为具有特殊生物学意义。
- 慢位点子集降低快速演化位点的贡献；快位点子集隔离该贡献。两者均不能确定真实拓扑。
- 对提取矩阵进行树推断和比较永不自动执行。
- 关于速率模型比较、慢/快位点子集敏感性解释，以及它与 heterotachy 的关系，见
  [系统误差工作流参考](../../skills/phyloai-workflow/references/syserror-workflow.md)。
