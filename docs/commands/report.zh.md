# phyloai report

[English](report.md) | [中文](report.zh.md)

## 目的

从一个 PhyloAI 运行目录生成可复现、可审计的分析报告。产出两个文件：

- **`report.json`** —— 每一步的机器可读记录（`params`、`key_results`、`tool_versions`、`methods_text`、输出文件路径）。作为 AI/MCP 诊断与可复现性审计的入口。
- **`report.html`** —— 人类可读报告，含嵌入图表、可排序表格、内联 CSV 数据，以及带一键复制的期刊级 Methods 段落草稿。

单次调用即可覆盖整个运行目录；无需子命令。

## 示例

```bash
phyloai report --run-dir <运行目录> [选项]
```

## 输入

`--run-dir` 必须包含一个或多个由 PhyloAI 生成的 `result.json` 文件。

## 输出

命令会在报告输出目录中写入 `report.json` 和 `report.html`。

## 用法

```bash
phyloai report --run-dir <运行目录> [选项]
```

## 快速开始

```bash
# 报告单条流水线运行（phyloai run 的输出）
phyloai report --run-dir ./runs/run/faa

# 报告模块级运行（如所有 pretree 步骤）
phyloai report --run-dir ./runs/pretree

# 覆盖已有报告
phyloai report --run-dir ./runs/run/faa --overwrite

# 自定义输出位置
phyloai report --run-dir ./runs/pretree -o ./my-report
```

## 输入要求

`--run-dir` 必须是一个包含一个或多个由 PhyloAI 命令生成的 `result.json` 文件的目录。报告自动识别两种结构：

| 结构 | 检测 | 典型用途 |
|------|------|----------|
| **pipeline** | 顶层 `run-dir/result.json` 由 `phyloai run` 产生，且子目录也包含 `result.json` | `phyloai run` 的输出 |
| **module** | 无顶层 `result.json`（子目录包含 `result.json`），或顶层 `result.json` 来自非 `run` 命令（可带或不带辅助子结果，如 `brlen label_nodes`） | `phyloai pretree`、`phyloai tree`、`phyloai posttree syserror brlen` 等 |

步骤发现完全基于文件系统 —— 目录扫描会排除 `report/`、`logs/` 与隐藏目录。

## 报告内容

HTML 报告包含五个面板：

### Panel A —— 运行汇总

汇总卡片，显示步骤计数（成功/失败）、总耗时，以及流水线运行的进度条。失败的步骤按名称列出。

### Panel B —— Methods

一段期刊级 Methods 段落，每个分析步骤对应一组句子。每个步骤都有可点击的 `[step_id]` 徽章，链接到其 Step Detail 卡片。"复制到剪贴板"按钮复制纯文本以便用于稿件。

### Panel C —— 步骤详情

每个步骤一张可折叠卡片，显示：
- 状态指示与耗时
- `↑ Methods` 返回链接
- 科学参数表
- 关键结果表
- 嵌入 CSV 表（≤200 行，≤500 KB），列可排序
- 警告与错误消息
- 完整 CLI 命令

### Panel D —— 图表

所有由分析命令生成的 PDF/PNG 图表，本地嵌入。矢量 PDF 保留质量。按分析阶段编号：`Fig-3.x`（pretree）、`Fig-4.x`（tree）、`Fig-5.x`（posttree）。

### Panel E —— 输出文件索引

可排序的表格，列出所有步骤的输出文件。大型表（>20 行）可折叠。每个条目链接到实际文件。

## 不完整运行

报告始终会成功，即使某些步骤已经失败：

- 失败的步骤会被包含，并显示完整错误详情，默认展开
- 失败步骤的 `methods_text` 为空，从 Methods 段落中排除
- 部分步骤失败时，流水线状态为 `"partial"`；全部失败时为 `"failed"`

## 选项

| Flag | Default | Description |
|------|---------|-------------|
| `--run-dir PATH` | *(required)* | 要报告的运行目录。 |
| `-o, --output-dir PATH` | `<run-dir>/report` | 报告文件输出目录。 |
| `--overwrite` | off | 覆盖已有报告文件。 |
| `-q, --quiet` | off | 除错误外不打印终端输出。 |

## 输出结构

```
<run-dir>/report/
├── report.json    # 机器可读的真相之源
└── report.html    # 自包含 HTML 报告（无外部依赖）
```

## report.json Schema

```json
{
  "phyloai_version": "0.1.0",
  "generated_at": "2026-06-27T14:23:00Z",
  "run_dir": "/abs/path/runs/pretree",
  "run_mode": "module",
  "status": "complete",
  "pipeline_summary": {
    "status": "complete",
    "n_steps_total": 5,
    "n_steps_success": 5,
    "n_steps_failed": 0,
    "n_steps_skipped": 0,
    "total_wall_time": 142.3
  },
  "steps": [
    {
      "step_id": "pretree.align",
      "command": "phyloai pretree align --seq-dir ./raw --method linsi",
      "status": "success",
      "wall_time": 31.4,
      "tool_versions": {"mafft": "7.526"},
      "params": {"method": "linsi", "seq_type": "AA", "threads": 8},
      "key_results": {"n_aligned": 100, "n_skipped": 0},
      "methods_text": "Multiple sequence alignments were performed...",
      "output_files": {},
      "warnings": [],
      "error": null
    }
  ],
  "methods_paragraph": "Multiple sequence alignments were performed...",
  "methods_blocks": [
    {"step_id": "pretree.align", "text": "Multiple sequence alignments...", "step_index": 0}
  ],
  "figures_index": [
    {
      "figure_id": "Fig-3.1",
      "step_id": "pretree.metrics",
      "label": "correlation_heatmap",
      "caption": "Correlation heatmap",
      "path": "/abs/path/correlation_heatmap.pdf",
      "type": "pdf"
    }
  ],
  "tables_index": [
    {
      "table_id": "Table-3.1",
      "step_id": "pretree.metrics",
      "label": "metrics_table",
      "caption": "Phylogenetic informativeness metrics per locus",
      "path": "/abs/path/metrics.csv",
      "type": "csv"
    }
  ]
}
```

## 关键结果补全

报告会从某些命令放在标准 `key_results` 字段之外的数据中自动补全 `key_results`：

- `data.summary` 中的标量（int、float、str）合并入 `key_results`
- `data.*` 顶层标量被合并（如 stats 单文件模式）
- 嵌套数值字典被展平：`{length_before: {mean: 10}}` → `length_before_mean`
- concat 特定度量（`gap_ratio`、`pi_ratio`）从 `data.variant_stats[0]` 提取

这确保模板无论哪个模块产出 `result.json`，总能获得完整数据。

## 退出码

| Code | Meaning |
|------|---------|
| `0` | 报告生成成功（即使运行中有失败步骤）。 |
| `1` | 用户输入错误 —— `--run-dir` 无效、未找到 `result.json`，或未加 `--overwrite` 时报告文件已存在。 |

## 警告与错误

| 条件 | 行为 |
|------|------|
| `--run-dir` 不是有效的运行目录 | Exit 1；"No result.json found." |
| 报告文件已存在 | Exit 1；使用 `--overwrite`。 |
| 步骤 `result.json` 含 `status: "error"` | 该步骤包含错误详情；Methods 段落排除之。 |
| `result.json` 损坏或不可读 | 该步骤记为错误；报告继续。 |

## 示例

```bash
# 报告已完成的流水线运行
phyloai report --run-dir ./runs/run/faa

# 报告 pretree 模块（多步）
phyloai report --run-dir ./runs/pretree

# 报告单个命令（如仅 1-convert）
phyloai report --run-dir ./runs/pretree/1-convert

# 用更新后的代码重新生成报告
phyloai report --run-dir ./runs/run/faa --overwrite

# 静默模式（仅错误输出到 stderr）
phyloai report --run-dir ./runs/pretree -q

# 自定义输出路径
phyloai report --run-dir ./runs/tree -o ./documents/methods
```

## 备注

- `report.html` 完全由 `report.json` 派生，可随时重新渲染而无需再次扫描运行目录。
- Methods 文本是确定性的（Python 模板函数，无 LLM）。所有科学上有意义的参数都会被描述；技术参数（线程、路径、标志）会被省略。
- PDF 图表通过 `<object>` 标签嵌入，保留矢量质量。不需要外部依赖（字体、CDN、JavaScript 库）。
- 报告被设计为 AI/MCP 诊断的主要入口 —— `report.json` 将所有步骤记录、参数、关键结果和图表路径聚合到一个可查询文档中。
- 小型 CSV/TSV 表（≤200 行，≤500 KB）以内联方式嵌入到 Step Detail 卡片中，作为可排序的 HTML 表。