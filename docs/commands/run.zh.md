# phyloai run

[English](run.md) | [中文](run.zh.md)

## 目的

从原始序列文件到物种树的一键式系统发育流水线。使用合理默认编排所有预处理与推断步骤。若需对任何单个步骤进行细粒度控制，请直接使用相应子命令。

## 用法

```bash
phyloai run --seq-dir <目录> [选项]
```

## 输入

`--seq-dir` 是包含原始序列文件的必需输入目录。

## 输出

命令在 `--output-dir` 下写入各步骤输出和顶层 `result.json`。

## 快速开始

```bash
# 典型用法：目录中的 AA 标记 → 物种树（supermatrix，normal 速度）
phyloai run --seq-dir ./markers

# 快速探索性运行
phyloai run --seq-dir ./markers --speed fast

# Supertree（基因树 → 通过 wASTRAL 的物种树）
phyloai run --seq-dir ./markers --mode supertree

# 恢复中断的运行
phyloai run --seq-dir ./markers --output-dir ./runs/run --resume

# 预览步骤而不执行
phyloai run --seq-dir ./markers --mode supertree --dry-run
```

## 输入要求

`--seq-dir` 必须是一个包含序列文件的目录。接受任何格式 —— 所有文件都会先转换为归一化的 FASTA。

支持的输入格式：`.fa`、`.fas`、`.fasta`、`.faa`、`.fna`、`.phy`、`.phylip`、`.nex`、`.nxs`、`.nexus`、`.aln`。

序列类型（AA 或 NT）自动检测。对于混合目录，请直接用 `phyloai pretree convert` 预先转换。

## 流水线模式

### `--mode supermatrix`（默认）

将每个基因比对、修剪、可选过滤，拼接为单一矩阵，并从该矩阵推断物种树。适用于基因共享共同历史（无基因树/物种树不一致）的数据集。

```
convert → align → trim → [filter] → concat → species tree
```

### `--mode supertree`

将每个基因比对、修剪、可选过滤，推断单个基因树，然后使用 wASTRAL 从基因树构建物种树。适用于预期存在基因树异质性（ILS、HGT、重复）时。

```
convert → align → trim → [filter] → gene trees → species tree
```

### 模式对比

| | Supermatrix | Supertree |
|---|---|---|
| 物种树方法 | IQ-TREE3 对拼接矩阵 | 基于基因树的 wASTRAL |
| 处理不一致 | 否 | 是 |
| 计算成本 | 较低（一棵树） | 较高（N 个基因树 + 溯祖） |
| 适用场景 | 基因树一致；基因共享历史 | 预期有 ILS、HGT 或基因重复 |
| Normal 速度树工具 | IQ-TREE3（不分区） | 每个基因用 IQ-TREE3 + wASTRAL |
| Fast 速度树工具 | 对矩阵用 FastTree | 每个基因用 FastTree + wASTRAL |

## 速度模式

| 步骤 | `--speed normal` | `--speed fast` |
|------|-----------------|---------------|
| Align | MAFFT `linsi`（最高精度） | MAFFT `auto`（启发式选择） |
| Trim | trimAl `-automated1` | trimAl `-automated1` |
| Filter | TAPER 错误位点屏蔽 | **跳过** |
| Gene trees | IQ-TREE3 带 ModelFinder | FastTree |
| Species tree（supermatrix） | IQ-TREE3（ModelFinder，不分区） | FastTree |
| Species tree（supertree） | wASTRAL mode 1 | wASTRAL mode 1 |

`--speed fast` 完全跳过 TAPER 过滤，并将 IQ-TREE3 换为 FastTree。用于快速探索性分析；用于发表质量结果请使用 `--speed normal`。

## 流水线步骤详解

### Step 1：Convert（`1-convert/`）

将所有输入序列文件转换为归一化 FASTA。处理格式检测、删除空行、多行序列合并、非法字符去除。输出文件写入 `1-convert/seqs/*.fa`。

使用 `phyloai pretree convert`。

### Step 2：Align（`2-align/`）

对每个 FASTA 文件运行多重序列比对。

- `normal`：MAFFT 使用 `--linsi`（迭代细化，最高精度）
- `fast`：MAFFT 使用 `--auto`（启发式方法选择）

输出：`2-align/seqs/*.fa`。

使用 `phyloai pretree align`。

### Step 3：Trim（`3-trim/`）

使用 trimAl `-automated1`（启发式选择 gap/相似度阈值）修剪对齐较差的列。减少比对不确定性带来的噪声。

输出：`3-trim/seqs/*.fa`。

使用 `phyloai pretree trim`。

### Step 4：Filter（`4-filter/`，仅 normal 速度）

使用 TAPER（位点级错误概率估计）屏蔽易出错位点。高错误概率的位点被移除，提升下游系统发育信号。

在 `--speed fast` 下，此步骤被跳过，不创建 `4-filter/`。

输出：`4-filter/seqs/*.fa`。

使用 `phyloai pretree filter`。

### Step 5：拼接或基因树

**Supermatrix（`5-concat/`）：** 将所有过滤/修剪后的比对拼接为一个超矩阵（`matrix.fa`）。输出一个指示基因边界的分区文件。

使用 `phyloai pretree concat`。

**Supertree（`5-genetrees/`）：** 为每个过滤/修剪后的比对推断一棵最大似然基因树。

- `normal`：每个基因用 IQ-TREE3 带 ModelFinder
- `fast`：每个基因用 FastTree

输出：`5-genetrees/trees/*.treefile`。

使用 `phyloai tree ml iqtree --msa-dir` 或 `phyloai tree ml fasttree --msa-dir`。

### Step 6：Species Tree（`6-tree/`）

**Supermatrix：** 从拼接矩阵推断物种树。

- `normal`：IQ-TREE3 带自动 ModelFinder（不分区；分区分析需要直接运行 `phyloai tree ml iqtree`）
- `fast`：FastTree

**Supertree：** 使用 wASTRAL（mode 1，无根）从基因树推断物种树。

使用 `phyloai tree ml iqtree`、`phyloai tree ml fasttree` 或 `phyloai tree msc`。

## 选项

| Flag | Default | Description |
|------|---------|-------------|
| `--seq-dir PATH` | *(required)* | 输入序列目录。任何格式；总是先转换。 |
| `--mode supermatrix\|supertree` | `supermatrix` | 流水线模式。 |
| `--speed normal\|fast` | `normal` | 速度/精度权衡。`fast` 跳过 TAPER 并使用 FastTree。 |
| `-o, --output-dir PATH` | `runs/run` | 所有流水线步骤的根输出目录。 |
| `-t, --threads INT` | `4` | 传递给所有步骤的线程数。 |
| `--resume` | off | 从 `run_checkpoint.json` 恢复。 |
| `--overwrite` | off | 删除并重建输出目录。与 `--resume` 互斥。 |
| `--dry-run` | off | 打印步骤列表而不执行。 |
| `-q, --quiet` | off | 除错误外不打印输出。 |

## 输出结构

```
runs/run/
├── run_checkpoint.json        # 流水线级 checkpoint（支持 resume）
├── result.json                # 整体流水线结果
├── 1-convert/
│   ├── result.json
│   └── seqs/                  # 归一化 FASTA 文件
├── 2-align/
│   ├── result.json
│   └── seqs/                  # 已比对 FASTA 文件
├── 3-trim/
│   ├── result.json
│   └── seqs/                  # 修剪后比对
├── 4-filter/                  # 仅 --speed normal
│   ├── result.json
│   └── seqs/                  # TAPER 过滤后比对
├── 5-concat/                  # supermatrix 模式
│   ├── result.json
│   └── matrix.fa              # 拼接后的超矩阵
├── 5-genetrees/               # supertree 模式
│   ├── result.json
│   └── trees/                 # 基因树文件
├── 6-tree/
│   └── result.json            # 物种树结果
```

每个步骤子目录包含自己的 `result.json`，含详细结果（`tool_versions`、`key_results`、`data`）。

## Checkpoint 与 Resume

`--resume` 从输出目录加载 `run_checkpoint.json`。checkpoint 记录每步的状态（`pending`、`running`、`success`、`failed`）、输出目录和参数哈希。

### Resume 行为

1. 校验 checkpoint 的 `params_hash` 是否与当前参数匹配。若参数已变，Exit 1（使用 `--overwrite` 进行全新运行）。
2. 跳过标记为 `success` 且其 `result.json` 包含 `"status": "success"` 的步骤。
3. 重跑标记为 `running` 或 `interrupted` 的步骤（从头开始，而非工具级 resume）。
4. 失败的步骤被标记为 `failed` 并使流水线中止。

注意：大多数单个步骤有自己的 checkpoint/resume 机制。流水线级 resume 较为粗粒度 —— 它重新运行被中断的步骤，而不是在工具原生 checkpoint 中恢复。

### Checkpoint JSON 结构

```json
{
  "schema_version": 1,
  "step": "run",
  "command": "phyloai run --seq-dir ./markers --mode supermatrix",
  "status": "success",
  "params_hash": "abc123...",
  "params": { "seq_dir": "/abs/path/to/markers", "mode": "supermatrix", ... },
  "started_at": "2025-01-01T00:00:00Z",
  "updated_at": "2025-01-01T01:00:00Z",
  "completed_at": "2025-01-01T01:00:00Z",
  "steps": [
    { "name": "convert",      "status": "success", "output_dir": "/abs/path/1-convert" },
    { "name": "align",        "status": "success", "output_dir": "/abs/path/2-align" },
    ...
  ]
}
```

## Result JSON

顶层 `result.json` 提供整个流水线的概览：

```json
{
  "status": "success",
  "command": "phyloai run --seq-dir ./markers ...",
  "wall_time": 1234.5,
  "tool_versions": { "mafft": "7.520", "iqtree3": "2.3.6", ... },
  "params": { ... },
  "key_results": {
    "n_input_genes": 50,
    "n_genes_after_filter": 47,
    "final_tree": "/abs/path/6-tree/matrix.fa.treefile",
    "matrix_length": 35000,
    "matrix_taxa": 100
  },
  "data": {
    "mode": "supermatrix",
    "speed": "normal",
    "steps": [
      { "name": "convert", "status": "success", "output_dir": "...", "result_json": "..." },
      ...
    ]
  }
}
```

`key_results` 字段因模式而异 —— `matrix_length` 和 `matrix_taxa` 仅出现在 supermatrix 模式。

## 必需工具

| Tool | Used By | Check With |
|------|---------|------------|
| MAFFT | Step 2（align） | `mafft --version` |
| trimAl | Step 3（trim） | `trimal --version` |
| TAPER | Step 4（filter，仅 normal 速度） | `taper --version` 或 `taper -h` |
| IQ-TREE3 | Steps 5/6（基因树、物种树、normal 速度） | `iqtree3 --version` 或 `iqtree2 --version` |
| FastTree | Steps 5/6（基因树、物种树、fast 速度） | `fasttree -h` |
| wASTRAL | Step 6（物种树、supertree 模式） | `wastral -h` |

运行 `phyloai doctor` 检查所有必需工具的环境。

## 退出码

| Code | Meaning |
|------|---------|
| `0` | 成功 —— 流水线完成，最终树写入 `6-tree/`。 |
| `1` | 输入错误 —— 缺少 `--seq-dir`、resume 时参数不匹配、`--resume` + `--overwrite` 同时使用、未加 `--overwrite` 时输出目录非空。 |
| `2` | 步骤失败 —— 一个流水线步骤（外部工具错误或数据错误）失败。检查失败步骤的 `result.json` 与日志。 |
| `3` | 环境错误 —— 必需工具未安装。运行 `phyloai doctor`。 |

## 警告 / 错误

| 条件 | 行为 |
|------|------|
| `--seq-dir` 不存在或为空 | 在任何步骤启动前 Exit 1。 |
| `--resume` 但无已有 checkpoint | Exit 1；使用 `--overwrite` 进行全新运行。 |
| `--resume` + `--overwrite` 同时使用 | Exit 1；互斥。 |
| resume 时参数不匹配 | Exit 1；checkpoint 哈希与当前参数不匹配。 |
| 输出目录非空且未加 `--overwrite` | Exit 1；使用 `--overwrite` 或 `--resume`。 |
| 步骤工具返回非零 | Exit 2；步骤在 checkpoint 中标记为 `failed`。 |
| 必需工具未找到（`FileNotFoundError`） | Exit 3；例如 `mafft`、`iqtree3`、`trimal`。 |
| 树步骤后未找到最终树文件 | Exit 2；树步骤声明成功但未产出输出。 |
| 中间 result.json 缺失或在 `--resume` 时不为 `"success"` | 重跑该步骤（忽略 checkpoint 状态）。 |

## 示例

```bash
# 默认：supermatrix，normal 速度
phyloai run --seq-dir ./markers

# 带显式选项的完整命令
phyloai run --seq-dir ./markers --mode supermatrix --speed normal --threads 8

# Supertree + fast 速度 + 16 线程
phyloai run --seq-dir ./markers --mode supertree --speed fast --threads 16

# 自定义输出目录
phyloai run --seq-dir ./markers -o ./runs/my_analysis

# 恢复之前中断的运行
phyloai run --seq-dir ./markers --output-dir ./runs/run --resume

# 预览步骤而不运行
phyloai run --seq-dir ./markers --mode supertree --dry-run

# 覆盖已有运行
phyloai run --seq-dir ./markers --overwrite
```

## 备注

- `phyloai run` 使用每步的默认参数。对于非默认设置（如分区 IQ-TREE、自定义 TAPER cutoff、特定的 MAFFT 方法），通过相应子命令单独运行各步骤。
- 在 supermatrix normal 模式下，IQ-TREE3 运行自动 ModelFinder 而不带分区文件。这是一遍不分区结果。分区分析需要直接运行 `phyloai tree ml iqtree`。
- 在 supertree normal 模式下，基因树推断使用 `phyloai tree ml iqtree --msa-dir`（带 ModelFinder 的批量 IQ-TREE3，每个基因一次）。
- `--threads` 传递给所有步骤。每个步骤可能以不同方式解读它：MAFFT 用作线程数；IQ-TREE 批量模式用作并行作业数；trimAl 用于并行化。
- TAPER（`4-filter/`）仅在 `--speed normal` 下运行。在 `--speed fast` 下，跳过该目录，修剪后的比对直接进入拼接或基因树推断。
- 每一步的工具版本都会聚合并写入顶层 `result.json` 的 `tool_versions`。
- 流水线在第一个失败步骤处停止。Checkpoint 状态在每步开始前和完成后保存，便于从失败步骤 `--resume`。
- 所有 convert/align/trim/filter 步骤在内部使用 `overwrite=True`（强制在流水线运行中重新生成输出），因为流水线通过其 checkpoint 管理目录状态。
- 在 `--speed fast` 下，基因树推断（supertree）和物种树推断（supermatrix）都使用 FastTree，速度更快但只产出近似树。FastTree 对 AA 使用 JTT+CAT 模型，对 NT 使用 GTR+CAT。