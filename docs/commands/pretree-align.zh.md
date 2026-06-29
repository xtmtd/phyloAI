# phyloai pretree align

[English](pretree-align.md) | [中文](pretree-align.zh.md)

## 目的

使用 MAFFT 或 MAGUS 对一个目录中的未比对序列文件进行多重比对。每个输入基因产出一个已比对 FASTA 文件。可通过 trimAl 回译额外产出密码子级 NT 比对。

该命令不进行格式转换。输入必须是 FASTA。若输入为 PHYLIP、Nexus 等其他格式，请先运行 `phyloai pretree convert --to fasta`。

## 用法

最简调用：
```bash
phyloai pretree align --seq-dir ./raw_aa
```

完整调用：
```bash
phyloai pretree align \
  --seq-dir ./raw_aa \
  --method linsi \
  --seq-type AA \
  --output-dir ./runs/pretree/align \
  --threads 4
```

带回译：
```bash
phyloai pretree align \
  --seq-dir ./raw_aa \
  --method linsi \
  --seq-type AA \
  --backtrans \
  --nt-dir ./raw_nt \
  --output-dir ./runs/pretree/align \
  --threads 4
```

## 参数

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--seq-dir` | required | 未比对序列文件目录 |
| `--method` | `linsi` | fftns1, fftns2, auto, linsi, einsi, ginsi, magus |
| `--seq-type` | `auto` | AA、NT 或 auto（从前几个基因自动检测） |
| `--backtrans` | off | 生成 NT 密码子比对；需要 `--nt-dir` |
| `--nt-dir` | — | 用于回译的未比对 CDS 目录 |
| `--output-dir` / `-o` | `runs/pretree/align` | 输出目录 |
| `--threads` / `-t` | 4 | 并行比对任务数（每个任务单线程） |
| `--tool-args` | — | 仅 MAGUS 策略参数；MAFFT 方法下被忽略 |
| `--mafft-path` | — | MAFFT 方法的显式可执行路径 |
| `--magus-path` | — | `--method magus` 的显式可执行路径 |
| `--trimal-path` | — | `--backtrans` 的显式 trimAl 可执行路径 |
| `--resume` | off | 从 `checkpoint.json` 恢复；要求完全相同的解析参数 |
| `--overwrite` | off | 删除并重建非空输出目录 |
| `--dry-run` | off | 仅打印命令，不创建文件 |
| `--quiet` / `-q` | off | 除错误外不打印终端输出 |

## 输入

扫描 `--seq-dir` 一层深度，匹配扩展名：`.fa`、`.fas`、`.fasta`、`.faa`、`.fna`。跳过子目录、空文件以及无法识别的扩展名。

## 输出

**AA 或 NT 模式：**
```
runs/pretree/align/
├── seqs/
│   ├── gene1.fa
│   └── ...
├── logs/
│   ├── gene1.log
│   └── ...
├── checkpoint.json
└── result.json
```

**AA + 回译模式：**
```
runs/pretree/align/
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

PhyloAI 作者化的 FASTA 输出使用统一的 FASTA 写入器，序列按 60 字符换行。

`result.json` 的 `key_results` 包含 `n_aligned`、`method`、`mean_alignment_length` 和 `mean_n_taxa`，便于整合到报告中。

每个基因的工具 stderr 写入 `logs/<locus>.log`。`result.json` 的 `data.files[]` 条目通过 `log_file` 引用这些日志。MAFFT 比对的 stdout 作为 FASTA 文件保存到 `seqs/`，不会重复记录到日志。

## 示例

```bash
# 大数据集的快速比对
phyloai pretree align --seq-dir ./raw_aa --method fftns2 --threads 8

# 高精度蛋白比对 + 密码子 NT 比对
phyloai pretree align --seq-dir ./raw_aa --seq-type AA \
  --backtrans --nt-dir ./raw_nt --method linsi --threads 4

# 直接 NT 比对
phyloai pretree align --seq-dir ./raw_nt --seq-type NT --method linsi

# MAGUS 带工具策略选项
phyloai pretree align --seq-dir ./raw_aa --method magus \
  --tool-args "--maxsubsetsize 200" --threads 4

# 预览命令而不执行
phyloai pretree align --seq-dir ./raw_aa --method linsi --dry-run

# 恢复中断的运行
phyloai pretree align --seq-dir ./raw_aa --method linsi --seq-type AA \
  --output-dir ./runs/pretree/align --resume
```

## 恢复行为

`pretree align` 支持 `--resume`，可以从中断、停电或外部工具失败中恢复，已完成的工作无需重做。

- 输出目录必须已包含 `checkpoint.json`。
- 当前调用的解析参数必须与 checkpoint 完全一致。包括分析参数以及 `--threads`、`--quiet` 等运行控制参数。
- 状态为 `success` 且输出文件仍有效的任务会被跳过。
- 状态为 `failed`、`pending`、`running`，或 `success` 但输出缺失/无效的任务会被重跑。
- 进度条仅在 checkpoint 校验后统计剩余可运行任务。已校验通过的任务会被汇总，不会重新走进度条。
- `--resume` 与 `--overwrite` 互斥。
- Resume 会以分隔符追加到每个位点的 `logs/<locus>.log`，并在完成时重写 `result.json`。

## 警告与错误

| 条件 | 行为 |
|------|------|
| `--backtrans` 未带 `--nt-dir` | Exit 1 |
| `--seq-type NT` 与 `--backtrans` 同时使用 | Exit 1 |
| `--seq-type auto` 检测到 NT 且使用 `--backtrans` | Exit 1 |
| 未找到 `mafft` 或 `magus` | Exit 3 |
| 使用 `--backtrans` 但未找到 `trimal` | Exit 3 |
| 在非 Linux 上使用 `--method magus` | Exit 1（MAGUS 捆绑二进制仅限 Linux） |
| 输出目录非空 | Exit 1（使用 `--overwrite`） |
| `--resume` 但无 checkpoint | Exit 1 |
| Resume 参数不匹配 | Exit 1 |
| CDS 长度不是 3 的倍数 | 该基因跳过回译，在 result.json 中产生警告 |
| CDS 中含有内部终止密码子 | 该基因跳过回译，在 result.json 中产生警告 |
| trimAl 退出非零 | 该基因跳过回译，stderr 捕获为警告 |
| 生成的 MSA 为空、不可解析或序列长度不一致 | 跳过该基因，原因记录在 result.json |
| 所有基因失败 | Exit 1 |
| MAFFT 方法下使用 `--tool-args` | 打印警告，参数被忽略 |

## 备注

- 下游：将 `phyloai pretree trim` 的 `--msa-dir` 指向 `seqs/`（模式 1/2）或 `seqs/faa/`（模式 3 AA）。
- `result.json` 的 `key_results` 供 Methods 段落使用："X 个基因使用 MAFFT L-INS-i 比对；平均比对长度 Y aa。"
- 运行 `phyloai doctor` 验证 MAFFT、MAGUS、trimAl 是否被检测到。