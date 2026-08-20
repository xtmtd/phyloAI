# phyloai posttree topology

[English](posttree-topology.md) | [中文](posttree-topology.zh.md)

## 目的

使用 IQ-TREE 对一组候选树与超矩阵比对进行拓扑检验（AU / KH / SH / WKH / WSH / c-ELW）。该命令检验备选拓扑是否显著差于最佳候选——**不会**推断新树。

## 用法

```bash
# 同质模型
phyloai posttree topology --matrix matrix.fa --candidate-trees candidates.trees --model-expr LG+F+R4

# PMSF 模型 + 引导树
phyloai posttree topology --matrix matrix.fa --candidate-trees candidates.trees --model-expr LG+C20+F+R4 --guide-tree guide.nwk

# 使用先前优化的分区模型
phyloai posttree topology --matrix matrix.fa --candidate-trees candidates.trees --partitions matrix.best_model.nex

# 多个独立树文件（逗号分隔，PhyloAI 合并）
phyloai posttree topology --matrix matrix.fa --candidate-trees h1.nwk,h2.nwk,h3.nwk --model-expr LG+F+R4

# 自定义交换矩阵 + 位点频率，通过 --tool-args
phyloai posttree topology --matrix matrix.fa --candidate-trees trees --model-expr custom.exchangeabilities+R4 --tool-args "-fs custom.sitefreq" -t 30

# 异质模型
phyloai posttree topology --matrix matrix.fa --candidate-trees trees --model-expr C20+F+R4
```

## 示例

```bash
phyloai posttree topology --matrix matrix.fa --candidate-trees candidates.trees --model-expr LG+F+R4
```

## 输入

| Input | Description |
|-------|-------------|
| `--matrix` | 单一超矩阵比对（FASTA、PHYLIP、NEXUS）。对应 IQ-TREE `-s`。 |
| `--candidate-trees` | 一个树列表文件（每行一棵 NEWICK 树）或多个以逗号分隔的独立 NEWICK 文件（如 `h1.nwk,h2.nwk`）。多个文件由 PhyloAI 按顺序合并为 `candidate.trees`。对应 IQ-TREE `-z`。 |
| `--input-format` | PhyloAI 侧的矩阵格式提示（`auto\|fasta\|phylip-relaxed\|nexus`，默认 `auto`）。不传给 IQ-TREE。 |

## 模型来源

请提供恰好一种模型来源。PhyloAI **不会**重跑 ModelFinder——模型选择请使用 `phyloai tree ml iqtree`。

| Option | Description |
|--------|-------------|
| `--model-expr` | 完整的 IQ-TREE `-m` 表达式。示例：`LG+F+R4`、`C20+F+R4`、`LG+C20+F+R4`、`custom.exchangeabilities+R4`。 |
| `--partitions` | 先前优化的分区文件（如 IQ-TREE 的 `.best_model.nex`）。对应 IQ-TREE `-p`。 |

`--guide-tree` 与 PMSF 模型（如 `LG+C20+F+R4`）配合使用。对应 IQ-TREE `-ft`。

## 默认检验

PhyloAI 生成标准的拓扑检验标志：

```
-n 0 -zb <replicates> -zw -au
```

| Test | Description |
|------|-------------|
| bp-RELL | Bootstrap proportion (RELL) |
| KH | Kishino-Hasegawa test |
| SH | Shimodaira-Hasegawa test |
| WKH | Weighted KH test |
| WSH | Weighted SH test |
| c-ELW | Expected likelihood weight |
| AU | Approximately unbiased test |

KH、SH、WKH、WSH、AU 是 **p 值**。p < 0.05 的树在该检验下被拒绝。bp-RELL 与 c-ELW 是 **权重**，不是 p 值。推荐：AU、WSH、WKH。

## 高级 IQ-TREE 参数

| Flag | Description |
|------|-------------|
| `--tool-args` | 额外的 IQ-TREE 策略参数。**被阻止的：** `-s`（矩阵）、`-z`（候选树）。Shell I/O 重定向（`<`、`>`、`|`）被拒绝。 |
| `--iqtree-path` | 显式 `iqtree3` 可执行文件路径。 |
| `--prefix` | IQ-TREE 输出前缀（默认：矩阵文件 stem）。 |

PhyloAI 构建的标志在 `--tool-args` 中出现同名标志时被抑制（suppress-if-present）。可被覆盖的标志：`-m`、`-p`、`-ft`、`-n`、`-zb`、`-zw`、`-au`、`-T`、`--prefix`。

IQ-TREE stdout 实时流式输出到终端，便于查看进度。结果从 `.iqtree` 报告中解析并以格式化表格呈现。

## 输出

IQ-TREE 原生文件：
- `<prefix>.iqtree` —— 完整的 IQ-TREE 报告，含拓扑检验表
- `<prefix>.log` —— IQ-TREE 日志
- `<prefix>.treels.trees` —— IQ-TREE 优化后的候选树（后缀可能变化）

PhyloAI 文件：
- `result.json` —— 结构化结果，含解析后的检验表
- `candidate.trees` —— 合并后的树文件（仅当提供了多个独立文件时）

## 退出码

| Code | Meaning |
|------|---------|
| 0 | 成功 |
| 1 | 用户输入错误（文件缺失、参数无效、输出冲突） |
| 2 | IQ-TREE 执行失败 |
| 3 | 未找到 IQ-TREE 可执行文件 |

## 警告与错误

- `--overwrite` 与 `--resume` 互斥。
- 输出目录非空时必须使用 `--overwrite` 或 `--resume`。

## 备注

- 该命令仅支持单矩阵（不支持批量模式）。
- `--replicates` 默认为 10000。过大的值会让 RELL 重采样变慢。
- 通过 `--resume` 支持 IQ-TREE 原生 resume（`.ckp.gz`）。
- 所有文件路径（矩阵、分区、引导树）在传入 IQ-TREE 前都会被解析为绝对路径。