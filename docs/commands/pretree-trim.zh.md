# phyloai pretree trim

[English](pretree-trim.md) | [中文](pretree-trim.zh.md)

## 目的

`phyloai pretree trim` 是 PhyloAI 的批量 MSA 修剪命令，用于已比对 FASTA 文件。它使用选定的后端工具（`trimAl`、`BMGE` 或 `ClipKIT`）一次修剪一个基因比对，并支持 AA-only、NT-only、CODON 以及 AA+NT 双输出工作流。

该命令不进行序列比对或格式转换。输入必须已经是已比对 FASTA 文件。请先运行 `phyloai pretree align`；当源文件不是 FASTA，或希望对 FASTA 输入多做一次归一化与校验时，先运行 `phyloai pretree convert --to fasta`。

## 用法

最简调用：
```bash
phyloai pretree trim --msa-dir ./aligned
```

AA + NT 带 trimAl 手动阈值：
```bash
phyloai pretree trim \
  --msa-dir ./runs/pretree/align/seqs/faa \
  --nt-dir ./runs/pretree/align/seqs/fna \
  --tool trimal \
  --tool-args "-gt 0.9 -cons 60" \
  --output-dir ./runs/pretree/trim
```

## 参数

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--msa-dir` | required | 输入已比对 MSA 文件目录 |
| `--output-dir` / `-o` | `runs/pretree/trim` | 输出目录 |
| `--tool` | `trimal` | `trimal`、`bmge` 或 `clipkit` |
| `--seq-type` | `auto` | `AA`、`NT`、`CODON` 或 `auto`；`auto` 仅区分 AA 与 NT |
| `--nt-dir` | — | 仅 AA+NT 模式；trimAl 接受原始 CDS 或带 gap 的密码子比对 NT，因为 PhyloAI 在回译前会先去除 NT gap；BMGE 与 ClipKIT 要求密码子比对 NT MSA |
| `--trimal-method` | `automated1` | trimAl 自动预设；当 `--tool-args` 包含 `-gt`、`-cons` 等手动阈值时被忽略 |
| `--bmge-matrix` | dynamic | AA/CODON 默认 `BLOSUM62`；NT 默认 `DNAPAM100:2` |
| `--bmge-entropy` | `0.5` | 越低越严格 |
| `--clipkit-method` | `smart-gap` | ClipKIT 模式 |
| `--trimal-path` | — | 显式 trimAl 可执行路径 |
| `--bmge-path` | — | 显式 BMGE.jar 路径 |
| `--clipkit-path` | — | 显式 clipkit 可执行路径 |
| `--threads` / `-t` | 4 | 并行修剪任务数 |
| `--tool-args` | — | 仅工具策略参数；PhyloAI 管理输入/输出/日志/密码子/线程标志 |
| `--resume` | off | 从 `checkpoint.json` 恢复 |
| `--overwrite` | off | 删除并重建非空输出目录 |
| `--dry-run` | off | 仅打印命令，不创建文件 |
| `--quiet` / `-q` | off | 除错误外不打印终端输出 |

## 输入

扫描 `--msa-dir` 一层深度，匹配扩展名：`.fa`、`.fas`、`.fasta`、`.faa`、`.fna`。跳过子目录、空文件以及无法识别的扩展名。

对于 AA+NT 模式：

- `trimAl`：`--msa-dir` 是已比对 AA MSA；`--nt-dir` 按 stem 提供匹配的 NT 文件。NT 输入可以是原始 CDS 或带 gap 的密码子比对 NT。PhyloAI 在调用 trimAl `-backtrans` 前会先去除 NT gap。
- `BMGE`：`--msa-dir` 是已比对 AA MSA；`--nt-dir` 是匹配的密码子比对 NT MSA。BMGE 修剪 AA，然后 PhyloAI 将保留的列投影到 NT。
- `ClipKIT`：`--msa-dir` 是已比对 AA MSA；`--nt-dir` 是匹配的密码子比对 NT MSA。ClipKIT 修剪 AA，然后 PhyloAI 将保留的列投影到 NT。

## 输出

**AA-only 或 NT-only：**
```
runs/pretree/trim/
├── seqs/
│   ├── gene1.fa
│   └── ...
├── logs/
│   ├── gene1.log
│   └── ...
├── checkpoint.json
└── result.json
```

**CODON 或 AA+NT：**
```
runs/pretree/trim/
├── seqs/
│   ├── faa/
│   │   └── gene1.fa
│   └── fna/
│       └── gene1.fa
├── logs/
│   ├── gene1.log
│   └── ...
├── checkpoint.json
└── result.json
```

`result.json` 报告修剪/跳过计数、跳过原因、修剪前/后的比对长度、解析后的参数与警告。

本命令产出的所有 PhyloAI 作者化 FASTA 系列输出均按 60 字符换行。

## 示例

```bash
# 默认 trimAl 自动修剪
phyloai pretree trim --msa-dir ./aligned_aa --tool trimal

# 通过 --tool-args 使用 trimAl 手动阈值
phyloai pretree trim --msa-dir ./aligned_aa --tool trimal \
  --tool-args "-gt 0.9 -cons 60"

# BMGE 处理密码子比对 NT 输入
phyloai pretree trim --msa-dir ./aligned_codon --seq-type CODON --tool bmge

# ClipKIT AA+NT 投影模式
phyloai pretree trim --msa-dir ./aligned_aa --nt-dir ./aligned_nt \
  --tool clipkit --clipkit-method gappy
```

## 警告与错误

| 条件 | 行为 |
|------|------|
| `--seq-type CODON` 与 `--nt-dir` 同时使用 | Exit 1 |
| 缺失外部工具 | Exit 3 |
| 输出目录非空且未加 `--overwrite` 或 `--resume` | Exit 1 |
| 无有效输入文件 | Exit 1 |
| AA+NT 模式缺失 NT 配对 | 跳过该基因 |
| trimAl backtrans 收到带 gap 的 NT 输入 | PhyloAI 先去除 gap 再继续 |
| 工具退出非零 | 跳过该基因，原因来自 stderr |
| 所有基因被跳过 | Exit 2 |

## 备注

- `--tool-args` 仅用于策略参数。不要传入工具管理的标志，如 trimAl `-in/-out/-backtrans`、BMGE `-i/-of/-t`、ClipKIT `-o/--codon`。
- 对于 trimAl，若 `--tool-args` 包含 `-gt`、`-cons` 等手动阈值，PhyloAI 不再追加默认自动预设。
- 每个基因的工具 stderr 写入 `logs/<locus>.log`。`result.json` 的 `data.files[]` 条目通过 `log_file` 引用这些日志。