# phyloai pretree convert

[English](pretree-convert.md) | [中文](pretree-convert.zh.md)

## 目的

`phyloai pretree convert` 归一化序列字符并转换 PhyloAI 工作流支持的序列/比对格式。仅适用于 FASTA、Phylip-relaxed、Phylip-PAML、Nexus 文件；不是通用格式转换工具。

## 用法

```bash
phyloai pretree convert --input ./raw --output-dir ./runs/pretree/convert --to fasta
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | 输入目录或单个文件 |
| `--output-dir` | `runs/pretree/convert` | 转换后文件与 result.json 所在目录 |
| `--to` | `fasta` | 目标格式：`fasta`、`phylip-relaxed`、`phylip-paml`、`nexus` |
| `--input-format` | `auto` | 覆盖输入格式检测 |
| `--seq-type` | `auto` | 覆盖分子类型检测 |
| `--aa-special` | `x` | 将 `B/Z/J/X/U/O` 转为 `X`，或用 `keep` 保留 |
| `--threads` | `4` | 目录模式工作进程数 |
| `--quiet` | false | 除错误外不打印 Rich 终端输出 |
| `--overwrite` | false | 删除并重建非空输出目录 |

## 终端输出

默认情况下，命令在终端显示 Rich 进度条与汇总表格。完整 JSON 结果写入 `--output-dir` 内的 `result.json`。使用 `--quiet` 抑制终端输出。

## 输入

`--input` 可以是目录或单个文件。目录模式仅扫描一层，跳过子目录、空文件、非序列文件以及无法解析的文件。

## 输出

转换后的序列文件写入 `--output-dir` 内的 `seqs/` 子目录。目标扩展名为 `.fa`、`.phy`、`.paml.phy`、`.nex`。

JSON 结果写入 `--output-dir` 内的 `result.json`。payload 在 `data` 下包含 `summary`、`files`、`skipped`、`warnings`。`key_results` 为空，因为 `convert` 是工具型命令。

对于 `phylip-paml`，输出记录使用 PAML 顺序式 header（含 `S`），将归一化后的分类单元名写入 30 字符字段，名称字段与序列之间至少两个空格。

`convert` 命令产出的所有 PhyloAI 作者化 FASTA 系列输出均按 60 字符换行。

示例输出结构：
```
runs/pretree/convert/
├── seqs/
│   ├── gene1.fa
│   ├── gene2.fa
│   └── ...
└── result.json
```

## 示例

```bash
phyloai pretree convert --input ./raw
phyloai pretree stats --seq-dir ./runs/pretree/convert/seqs
phyloai pretree convert --input ./gene.phy --output-dir ./converted --to fasta --seq-type NT
phyloai pretree convert --input ./aligned --to phylip-paml --overwrite
```

## 警告与错误

若部分文件无效，它们会被跳过并列在输出中。若所有输入都失败或被跳过，命令以退出码 1 退出。若输出目录已存在且非空，使用 `--overwrite` 覆盖。

## 备注

当原始输入文件可能包含混合格式或非标准字符时，在 `pretree stats` 之前使用 `pretree convert`。

即使源文件已经是 FASTA，在数据来源不确定时，下游步骤前仍建议运行 `pretree convert --to fasta`。它会重新检查解析、归一化序列字符，并能尽早发现畸形记录、意外符号或其他输入问题。