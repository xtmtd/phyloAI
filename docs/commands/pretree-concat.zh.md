# phyloai pretree concat

[English](pretree-concat.md) | [中文](pretree-concat.zh.md)

## 目的

`phyloai pretree concat` 将多个 MSA 文件拼接为一个超矩阵，用于下游系统发育推断。支持占用率过滤、字符重编码、密码子相关变体（翻译、排除第三位密码子）、外类群重排以及多格式输出。

输入必须是已比对 MSA 文件（FASTA 或 Phylip）。请先运行 `phyloai pretree align`；若需删除对齐较差的列，先运行 `phyloai pretree trim`。

## 用法

最简调用：
```bash
phyloai pretree concat --msa-dir ./aligned
```

带重编码与 CODON 变体：
```bash
phyloai pretree concat \
  --msa-dir ./runs/pretree/trim/seqs/fna \
  --seq-type CODON \
  --recoding RY-nucleotide \
  --translate-codon \
  --exclude-codon3 \
  --outgroup Homo_sapiens \
  --taxa-occupancy 0.7 \
  --to fasta \
  --output-dir ./runs/pretree/concat
```

## 参数

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--msa-dir` | required | 输入 MSA 文件目录 |
| `--output-dir` / `-o` | `runs/pretree/concat` | 输出目录 |
| `--prefix` | `matrix` | 输出文件名前缀 |
| `--seq-type` | `auto` | `AA`、`NT`、`CODON` 或 `auto`；`CODON` 永远不会被自动检测 |
| `--taxa-occupancy` | `0.5` | MSA 包含的最低分类单元占比（0.0–1.0） |
| `--recoding` | — | `RY-nucleotide`（仅 NT；A/G→R，C/T/U→Y），`Dayhoff-6/9/12/15/18`、`SandR-6`、`KGB-6`（仅 AA） |
| `--outgroup` | — | 移到第一位的单个分类单元名 |
| `--to` | `fasta` | 输出格式：`fasta`、`phylip-relaxed`、`phylip-paml`、`nexus` |
| `--translate-codon` | off | 同时生成 CDS→AA 翻译矩阵（仅 CODON） |
| `--exclude-codon3` | off | 同时生成 codon1+2 矩阵（仅 CODON） |
| `--dry-run` | off | 校验输入并报告计划动作，不写、不删、不替换任何文件 |
| `--quiet` / `-q` | off | 除错误外不打印终端输出 |
| `--overwrite` | off | 删除并重建非空输出目录 |

## 输入

扫描 `--msa-dir` 一层深度，匹配扩展名：`.fa`、`.fas`、`.fasta`、`.faa`、`.fna`、`.phy`。跳过子目录、空文件以及无法识别的扩展名。

序列在拼接前通过 `core/sequence_normalization.py` 按基因归一化。某基因中缺失的分类单元用 `?` 填充。

## 输出

```
runs/pretree/concat/
├── matrix.fa                   # （或 .phy/.nex）
├── matrix.partitions           # RAxML 风格分区文件
├── matrix.recoded.fa           # 使用 --recoding 时
├── matrix.recoded.partitions   # 使用 --recoding 时
├── matrix.translated.fa        # 使用 --translate-codon 时
├── matrix.translated.partitions # 使用 --translate-codon 时
├── matrix.cds12.fa             # 使用 --exclude-codon3 时
├── matrix.cds12.partitions     # 使用 --exclude-codon3 时
├── dropped_alignments.csv      # 有 MSA 被丢弃时
├── result.json
```

### 变体

| 变体 | 条件 | `seq_type` |
|------|------|-----------|
| Original | 总是 | 从输入解析 |
| Recoded | `--recoding` | `other` |
| Translated | `--translate-codon` | `AA` |
| Codon1+2 | `--exclude-codon3` | `NT` |

### 分区文件

每个矩阵都附带一个 RAxML 风格的 `.partitions` 文件，描述超矩阵中各基因的边界，供 IQ-TREE 等工具做分区分析（`-p` 选项）。

文件中每行格式为：

```
TYPE, gene_name = start-end
```

**前缀规则：**

| 矩阵变体 | `TYPE` |
|---|---|
| Original（NT / CODON） | `DNA` |
| Original（AA） | `LG` |
| Recoded（任意） | `AUTO` |
| Translated（CODON→AA） | `LG` |
| Codon1+2（CODON→NT） | `DNA` |

基因名使用输入文件 basename（不含扩展名）。位置为 1-based 的闭区间。对于 translated/cds12 变体，位置会重算以匹配变体矩阵的长度。`--dry-run` 下不写分区文件。

示例：
```
DNA, COI = 1-654
DNA, 16S = 655-1203
```

### result.json

已生成与计划生成的变体输出以完整路径形式记录在 `key_results.variants_produced` 与 `data.variants[].path` 中。本命令产出的所有 PhyloAI 作者化 FASTA 系列输出均按 60 字符换行。

## 屏幕显示（Rich）

三个面板：

1. **Overview** —— prefix、to_format、n_taxa、n_msa_* 计数、taxon_occupancy_threshold、recoding、outgroup、生成的变体文件。

2. **Character Summary** —— 按变体的表格：seq_type、total_length、gap_ratio、ambiguous_ratio、gap_ambiguous_ratio、standard_ratio。Recoded 变体（`seq_type = "other"`）显示有意义的 gap_ratio，不确定指标显示 `—`。

3. **Site Patterns** —— 按变体的表格：alignment_length、distinct_patterns（计数 + 比例）、constant_sites、parsimony_informative、singleton_sites。比例保留 4 位小数。Distinct-pattern 计数把所有非标准字符归并为 gap 符号，与 IQ-TREE 约定一致。

## 示例

```bash
# 基础 NT 拼接
phyloai pretree concat --msa-dir ./aligned_nt --seq-type NT --to fasta

# CODON 带所有变体
phyloai pretree concat --msa-dir ./aligned_codon --seq-type CODON \
  --translate-codon --exclude-codon3

# 重编码 + 外类群
phyloai pretree concat --msa-dir ./aligned_aa --seq-type AA \
  --recoding Dayhoff-6 --outgroup Sp_A

# 严格占用率，先 dry-run
phyloai pretree concat --msa-dir ./aligned --taxa-occupancy 1.0 --dry-run
```

## 警告与错误

| 条件 | 行为 |
|------|------|
| 缺少 `--msa-dir` 或未找到 MSA 文件 | Exit 1 |
| 输出目录非空且未加 `--overwrite` | Exit 1 |
| 在非 CODON seq_type 上使用 `--translate-codon` / `--exclude-codon3` | Exit 1 |
| 在 NT 输入上使用仅 AA 的重编码方案（或反之） | Exit 1 |
| `--outgroup` 指定的分类单元在矩阵中找不到 | Exit 1 |
| 没有 MSA 通过 `--taxa-occupancy` 过滤 | Exit 1 |
| `--taxa-occupancy` 超出 0.0–1.0 | Exit 1 |

## 备注

- 两遍处理以节约内存：第一遍仅扫描表头（FASTA）收集分类单元集合并按占用率过滤；丢弃的文件不会被完整读取。第二遍进行流式拼接。
- 翻译按基因进行（先于拼接）以保持基因边界处的密码子相位。
- 原始变体的统计在格式转换前于内存中计算，避免 Phylip-PAML 输出时的名称截断问题。
- `result.json` 的 `variant_stats` 包含各变体的 character summary 与 site patterns。
- `.partitions` 文件随每个矩阵一起生成，供分区系统发育分析使用（如 `iqtree -s matrix.fa -p matrix.partitions`）。
- `--dry-run --overwrite` 仍保留任何已有输出目录不动；`--overwrite` 仅在真实运行时删除文件。