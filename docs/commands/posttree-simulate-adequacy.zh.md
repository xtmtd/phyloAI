# phyloai posttree simulate adequacy

[English](posttree-simulate-adequacy.md) | [中文](posttree-simulate-adequacy.zh.md)

## 用途

将观测 MSA 的 PPA-DIV、PPA-CONV、PPA-VAR 和 PPA-COMP 统计量与模拟 MSA 的经验零分布比较。该命令仅在本地运行，不需要外部可执行程序。

## 用法

```bash
phyloai posttree simulate adequacy \
  --original-msa matrix.fa \
  --simulated-dir runs/sim/MSAs \
  --threads 4 \
  --table-format csv \
  --output-dir runs/adequacy
```

## 参数

| 选项 | 说明 |
|--------|-------------|
| `--original-msa` | 观测比对，必需。自动检测 FASTA、PHYLIP-relaxed、PHYLIP-PAML 和 NEXUS；分类单元 ID 必须唯一。 |
| `--simulated-dir` | 模拟比对目录，必需。每个非空常规文件独立自动检测，且必须与观测 MSA 有相同的分类单元集合和比对长度。 |
| `--seq-type` | `AA`、`NT` 或 `auto`（默认，大小写不敏感）。`auto` 仅根据观测 MSA 判定类型，并在工作进程前完成解析。 |
| `--threads` | 模拟 MSA 统计的并行工作进程数（默认 `4`）。 |
| `--table-format` | 三个输出表的分隔符和后缀：`csv`（默认）或 `tsv`。 |
| `-o, --output-dir` | 输出目录（默认 `runs/posttree/simulate/adequacy`）。 |
| `--overwrite` | 删除并重建非空输出目录；不能与 `--resume` 同时使用。 |
| `--resume` | 从 `checkpoint.json` 恢复。观测 MSA 和恢复所需参数必须与原运行一致；已替换的模拟文件会重算。 |
| `--dry-run` | 验证观测 MSA、判定序列类型、统计模拟文件，但不写任何输出。 |
| `-q, --quiet` | 除错误外不输出终端信息。 |

## 示例

```bash
# 氨基酸数据
phyloai posttree simulate adequacy --original-msa concat.aa.fa \
  --simulated-dir runs/sim/MSAs --seq-type AA --threads 4

# 恢复写入 TSV 表格的运行
phyloai posttree simulate adequacy --original-msa concat.aa.fa \
  --simulated-dir runs/sim/MSAs --table-format tsv -o runs/adequacy --resume
```

## 输出

```
runs/posttree/simulate/adequacy/
├── adequacy_summary.csv      # 标量 PPA-DIV/CONV/VAR/COMP 结果
├── adequacy_taxon_comp.csv   # 逐分类单元 PPA-COMP 结果
├── per_simulation_stats.csv  # 每个有效模拟的原始标量统计量
├── checkpoint.json           # 可恢复的逐模拟统计量
└── result.json               # 机器可读结果
```

使用 `--table-format tsv` 时，三个表会使用 `.tsv` 后缀和制表符分隔；`result.json` 在 `data.output_files` 中记录其解析后的路径。

`adequacy_summary` 为 `div`、`siteconvprob`、`sitecomp`、`comp_max` 和 `comp_mean` 报告观测值、模拟均值和总体 SD、经验 95% 区间、z-score、后验预测 p 值及有效模拟数。`adequacy_taxon_comp` 对每个观测分类单元给出相同比较。

## 警告与错误

| 条件 | 行为 |
|-----------|----------|
| 少于 10 个有效模拟 | 处理后硬错误；不生成适当性汇总。 |
| 观测 MSA 有重复分类单元 ID、序列长度不等、分类单元全缺失或无信息位点 | 硬输入错误。 |
| 模拟 MSA 有重复 ID、分类单元/长度不匹配、分类单元全缺失或无法解析 | 跳过该文件，记录为失败并写入警告。 |
| 非空输出目录未使用 `--overwrite` 或 `--resume` | 预检错误；不改变任何输出。 |
| `--resume` 的 checkpoint 缺失、不兼容、缺少观测 MSA 指纹或观测 MSA 已改变 | 预检错误；使用 `--overwrite` 开始新运行。 |

## 退出码

| 代码 | 含义 |
|------|---------|
| `0` | 成功。 |
| `1` | 用户输入、验证、输出冲突或恢复错误。 |

## 说明

- `div` 的 pp 为 `P(sim <= obs)`；其他 pp 为 `P(sim > obs)`。低 pp（`< 0.05`）或 `|z| > 2` 表示潜在模型适当性不足。
- 观测和模拟 MSA 可以混用上述支持格式。
- 当观测 MSA 缺失数据较多而模拟序列无 gap 时，应在 adequacy 前运行 `phyloai posttree simulate alisim transfergaps`。
- 所有统计量排除非标准字符：AA 为 `ACDEFGHIKLMNPQRSTVWY`，NT 为 `ACGT`。模拟分布 SD 为零时，JSON 中 pp 为 `null`，表格 pp 单元格为空。
