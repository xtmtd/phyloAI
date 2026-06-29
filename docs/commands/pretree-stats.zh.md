# phyloai pretree stats

[English](pretree-stats.md) | [中文](pretree-stats.zh.md)

## 目的

`phyloai pretree stats` 是用于序列与比对文件的只读 QC 与检查命令。它汇总格式、序列类型、比对状态、分类单元数、长度分布、gap 比例、歧义字符比例，以及对齐文件的位点模式统计。

它不归一化、转换、比对、修剪或修改输入文件。当原始输入可能包含混合格式或应在校验前标准化的字符时，先使用 `phyloai pretree convert`。

## 用法

```bash
phyloai pretree stats [OPTIONS]
```

必须且只能指定 `--seq` 或 `--seq-dir` 之一。

## 参数

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--seq-dir DIRECTORY` | none | 目录模式。扫描一个文件夹中支持的序列/比对文件并计算数据集级汇总。 |
| `--seq FILE` | none | 单文件模式。详细检查一个序列或比对文件。 |
| `--unaligned` | `False` | 仅目录模式。将输入视为未比对序列。每基因 CSV 排除 `alignment_length` 与位点模式列，包含 `seq_length_*` 列。 |
| `--per-gene` | `False` | 仅目录模式。将每基因结果写入输出目录的 CSV/TSV 文件。不打印终端表格；用此文件做详细的逐位点检查。 |
| `--table-format csv\|tsv` | `csv` | 仅目录模式。配合 `--per-gene` 写入的每基因表格式。 |
| `--output-dir DIRECTORY` | `runs/pretree/stats` | 写入 `result.json` 与每基因文件的目录。 |
| `--input-format fasta\|phylip-relaxed\|nexus` | auto | 覆盖自动格式检测。 |
| `--seq-type AA\|NT` | auto | 覆盖自动分子类型检测。 |
| `--threads INTEGER`, `-t INTEGER` | `4` | 仅目录模式。工作进程数。至少 `1`。 |
| `--quiet`, `-q` | `False` | 除错误外不打印终端输出。 |
| `--overwrite` | `False` | 若输出目录已存在且非空，则删除并重建。 |

## 输入

单文件模式用 `--seq` 读一个文件。目录模式用 `--seq-dir` 扫描一个目录，处理支持的序列/比对扩展名。

支持的输入格式为 FASTA、Phylip-relaxed、Nexus。经典 `phylip` 不作为独立的 `stats` 选项暴露，以免与 Phylip-relaxed 混淆。

每个文件的比对状态自动检测。当文件包含多于一条序列且所有序列等长时视为已比对。单序列文件视为未比对。

## 输出

除 `--quiet` 外，终端输出使用 Rich 表格与面板。

单文件模式显示概览、字符汇总、每个分类单元统计，以及文件已比对时的位点模式统计。

目录模式显示汇总表。`--per-gene` 将每基因表写入输出目录（仅文件，不打印到终端）。

结果始终写入输出目录：

```
runs/pretree/stats/
├── result.json           # JSON 结果（始终写入）
└── per-gene.csv          # 每基因表（使用 --per-gene 时）
```

每基因文件使用 `--table-format` 指定的格式（默认 csv）。

## 示例

检查一个已比对氨基酸文件：

```bash
phyloai pretree stats --seq ref/phylogenomics_examples/test/EOG090X0971.faa
```

检查一个未比对核苷酸文件：

```bash
phyloai pretree stats --seq ref/phylogenomics_examples/2-loci_filter/fna/EOG090X0971.fna
```

汇总一个目录：

```bash
phyloai pretree stats --seq-dir ref/phylogenomics_examples/3-align/faa
```

保存目录汇总与每基因 CSV：

```bash
phyloai pretree stats \
  --seq-dir ref/phylogenomics_examples/2-loci_filter/fna \
  --per-gene \
  --output-dir runs/pretree/stats
```

汇总未比对序列并使用适当的每基因列：

```bash
phyloai pretree stats \
  --seq-dir ref/phylogenomics_examples/2-loci_filter/fna \
  --per-gene \
  --unaligned \
  --output-dir runs/pretree/stats
```

当设置了 `--unaligned` 时，每基因 CSV 包含 `seq_length_*` 列，不含 `alignment_length` 与位点模式列。默认（无 `--unaligned`）产出相反的列集，适合已比对 MSA。

将每基因表保存为 TSV：

```bash
phyloai pretree stats \
  --seq-dir ./data \
  --per-gene \
  --table-format tsv \
  --output-dir runs/pretree/stats
```

原始输入归一化后的推荐顺序：

```bash
phyloai pretree convert --input ./raw --output-dir ./runs/pretree/convert --to fasta
phyloai pretree stats --seq-dir ./runs/pretree/convert/seqs
```

## 警告与错误

`--seq` 与 `--seq-dir` 互斥。同时不指定或同时指定两者都是输入错误。

`--per-gene` 仅用于目录模式。将其与 `--seq` 一起使用是输入错误。

`--threads` 必须至少为 `1`。

若输出目录已存在且非空，命令报错退出。使用 `--overwrite` 替换。

若序列类型检测有歧义，命令默认为 `AA` 并发出警告。

若任意序列中出现 `*`，命令发出警告，因为它可能表示终止密码子或上游处理问题。

**混合状态警告：** 若每个文件的自动检测与声明的比对模式（默认已比对或 `--unaligned`）不一致，命令会发出警告。例如，在大多数文件序列长度不一致的目录上不指定 `--unaligned` 运行会产生：
> `--unaligned is NOT set, but 1056 of 1066 files were detected as unaligned ... Use --unaligned to write unaligned-specific columns to the per-gene table.`

反过来，在已比对数据上使用 `--unaligned` 会警告去掉该标志。

## 备注

字符类别被报告为标准、gap/missing 或歧义。位点模式统计在确定简约信息位点与 singleton 位点时排除 gap/missing/歧义字符。