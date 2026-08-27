# phyloai posttree syserror brlen

[English](posttree-syserror-brlen.md) | [中文](posttree-syserror-brlen.zh.md)


## 目的

从系统发育树中提取枝长统计量，用于诊断**跨类群速率异质性（枝长异质性）**——系统发育组学中系统误差的一个主要来源。通过在**不同替换模型**下推断的树之间比较枝长，评估长枝吸引（LBA）与模型依赖的枝长估计。本命令提供原子化的枝长提取测量；它本身不判定 LBA 因果性。

纯 Python 实现（Bio.Phylo）——不调用任何外部可执行程序，因此无需 `phyloai doctor` 检查。

## 用法

```bash
phyloai posttree syserror brlen [options]
phyloai posttree syserror brlen label-nodes --tree <tree.nwk> [options]
```

## 输入

| 选项 | 类型 | 必填 | 默认值 | 说明 |
|--------|------|----------|---------|-------------|
| `--tree` | Path | 互斥 | — | 单个树文件（Newick）。与 `--tree-dir` 互斥。 |
| `--tree-dir` | Path | 互斥 | — | 树文件目录。与 `--tree` 互斥。 |
| `--mode` | str | 是 | — | 逗号分隔的模式列表（见“模式”）。 |
| `--map` | Path | 否 | — | 节点-物种映射文件，用于节点识别。 |
| `--node1` | str | 否 | — | 第一个节点名（端点模式）。 |
| `--node2` | str | 否 | — | 第二个节点名（node-to-node 模式）。 |
| `--tip1` | str | 否 | — | 第一个末梢物种名。 |
| `--tip2` | str | 否 | — | 第二个末梢物种名（tip-to-tip 模式）。 |
| `-o`, `--output-dir` | Path | 否 | `runs/posttree/syserror/brlen` | 输出目录。 |
| `--table-format` | csv\|tsv | 否 | csv | 输出表格分隔符。 |
| `-t`, `--threads` | int | 否 | 4 | 批处理模式的并行工作进程数。 |
| `--max-rows` | int | 否 | 5000000 | patristic 输出行数安全上限。0 = 不限。 |
| `--overwrite` | flag | 否 | False | 删除并重建输出目录。 |
| `--dry-run` | flag | 否 | False | 仅校验输入（含端点解析），不写入文件。 |
| `-q`, `--quiet` | flag | 否 | False | 除错误外抑制终端输出。 |

## 模式

### 批量模式（可逗号组合，无需端点参数）

| 模式 | 说明 | 输出列 |
|------|-------------|----------------|
| `total` | 每棵树所有枝长之和 | tree_file, total_branch_length |
| `terminal` | 每个末梢（tip）枝 | tree_file, taxon, branch_length |
| `internal` | 每个非根内部枝 | tree_file, representation, edge_taxa, branch_length |
| `patristic` | 所有两两末梢间距离 | tree_file, tip1, tip2, distance |
| `all` | = total + terminal + internal + patristic | 四张 CSV |

`internal.csv` 对根/非根混合批次使用统一 schema：`representation` 为
`rooted` 或 `unrooted`，`edge_taxa` 对有根边列出后代类群，对无根边列出
**规范分割侧**（较小的叶集合，平分时按字典序取第一侧）。根被排除，因为它没有入边。

### 端点模式（一次一个，不可与批量模式组合）

| 模式 | 必填参数 | 输出列 |
|------|--------------------|--------------------|
| `tip-to-tip` | `--tip1`、`--tip2` | tree_file, tip1, tip2, distance |
| `node-to-node` | `--node1`、`--node2`、（`--map` 或带标签的树） | tree_file, node1, node2, node1_type, node2_type, distance |
| `node-to-tip` | `--node1`、（`--map` 或带标签的树；`--tip1` 可选） | tree_file, node, node_type, tip, distance |

`node_type` 为 `internal` 或 `tip`（单物种映射项解析为该末梢；此列使区分显式化）。
端点模式与批量模式互斥，端点模式之间也互斥。

## 节点识别

基于节点的模式（`node-to-node`、`node-to-tip`）按以下顺序识别内部节点：

1. 若提供了 `--map` 文件，则优先使用（始终覆盖标签）。
2. 否则使用内部节点标签（如来自 `label-nodes` 的 `N1`）。

映射文件格式（冒号分隔节点名与逗号分隔的物种列表；空白会被去除；空行和
无冒号的行被跳过）：

```
NodeName:sp1,sp2,sp3
Outgroup:spA,spB
```

`--map` 使用每棵树中存在的类群。当存在的子集恰好构成一个有根支序或无根分割
时可解析；空交集或不兼容的组会发出警告并跳过该树的端点计算。`Nxx` 标签仅
适用于生成它们的参考拓扑。

- **有根树**（根恰好 2 个子节点）：映射组必须等于 MRCA 的精确后代叶集合（单系支序）。
- **无根树**（根 3 个以上子节点）：映射组必须等于某内部分割的一侧。两侧都会被测试；
  若仅补集（非直观）侧匹配，则记录一条警告。

`node-to-tip` 省略 `--tip1` 时：
- 有 `--map`：对（映射组 ∩ 树末梢）中的每个类群输出一行。
- 无 `--map`、有根带标签树：对标签节点的每个后代输出一行。
- 无 `--map`、无根树：报错——后代推断有歧义。

## 有根与无根表示

表示方式依据 Newick 根子节点数结构性地判定：2 个子节点 = 有根，3 个以上 =
无根。这遵循 IQ-TREE、RAxML、FastTree、wASTRAL 和 gotree 的惯例。这是**结构性
启发式，并非生物学有根性的证明**——同一无根树可以写成二叉根，单根或根处多歧
也被视为无根。需要根语义（有根 `edge_taxa`、有根标记、基于后代的 `node-to-tip`）
的用户必须提供二叉有根 Newick 表示。

## 输出

```
runs/posttree/syserror/brlen/
├── result.json
└── tables/
    ├── total.csv
    ├── terminal.csv
    ├── internal.csv
    ├── patristic.csv
    ├── tip_to_tip.csv
    ├── node_to_node.csv
    └── node_to_tip.csv
```

仅为所请求的模式创建表格；扩展名随 `--table-format`（`.csv` 或 `.tsv`）。
所有成功的非 dry-run 运行只写一个根 `result.json`，包含状态、解析后的参数、
关键结果（树数量、模式、各模式的 `n_values`/均值/总体标准差/最小/最大）、
警告和 `data.output_files`。`data.warnings` 记录被跳过的树、patristic 行估计
和端点跳过原因。`tool_versions` 为 `{}`（纯 Python）。不支持 `--resume` 和
检查点：这是一次性工具。

`terminal` 报告每棵树的**全部**末梢枝；如需单个类群请过滤其表格，而不是请求
单末梢终端模式。

## 示例

```bash
# 批量：从后验树目录提取所有枝长（支持多树文件）
phyloai posttree syserror brlen --tree-dir ./posterior_trees --mode all --max-rows 0

# 单个：末端和内部枝长
phyloai posttree syserror brlen --tree LG.tre --mode terminal,internal

# 单个：两个类群之间的末梢距离
phyloai posttree syserror brlen --tree LG.tre --mode tip-to-tip --tip1 Neelus_murinus --tip2 Folsomia_candida

# 批量：使用 map 跨模型树计算节点间距离
phyloai posttree syserror brlen --tree-dir ./model_trees --mode node-to-node \
    --map nodes.map.txt --node1 Collembola --node2 Outgroup

# 批量：跨后验树计算节点到末梢（map 定义类群）
phyloai posttree syserror brlen --tree-dir ./posterior_trees --mode node-to-tip \
    --map nodes.map.txt --node1 Collembola

# 有根带标签树上指定末梢的节点到末梢（无需 map）
phyloai posttree syserror brlen --tree species.labeled.nwk --mode node-to-tip \
    --node1 N5 --tip1 Folsomia_candida

# 有根带标签树节点到全部后代的节点到末梢（无需 map）
phyloai posttree syserror brlen --tree species.labeled.nwk --mode node-to-tip --node1 N5

# 为参考树生成带标签的树和 map 模板
phyloai posttree syserror brlen label-nodes --tree species.nwk
```

### label-nodes

为单个树的内部节点编号 `N1..Nxx`（前序），并写出：

| 文件 | 说明 |
|------|-------------|
| `<stem>.labeled.nwk` | 带内部节点标签的 Newick |
| `<stem>.map.txt` | 供主命令使用的节点-物种映射模板 |
| `<stem>.labeled.pdf` | 树可视化（matplotlib Agg） |
| `result.json` | 标准 PhyloAI 结果 |

有根树标记包括根在内的所有内部节点；无根树排除人为根。默认输出目录：
`runs/posttree/syserror/brlen/label_nodes`。

标签不补零（`N1`、`N2`、...、`Nxx`）。标记会清除每个节点的数值 support
（`confidence` 字段），因此 `labeled.nwk` 只含干净的 `Nxx:length` 记号，
不会混入标签与 support 文本，PDF 也只显示标签。枝长以无损方式写出
（`format_branch_length="%r"`）——保留解析出的浮点值，不做默认的 5 位小数截断。

## 警告 / 错误

- 退出码 **0**：成功。**1**：输入校验错误（无有效树、模式组合非法、缺少必填
  参数、patristic 行数超限、单模式下端点无法解析）。不调用外部工具，因此退出
  码 2/3 不适用。
- 单个 `--tree` 输入：解析失败为退出码 1 错误。`--tree-dir` 模式下，无效文件/
  树被跳过并告警；仅在无有效树剩余时失败。
- 单端点模式：无法解析的 tip/node/map 端点为退出码 1 错误。批量端点模式：
  仅警告并跳过该棵树，并计入 `n_trees_skipped`。
- 少于两个末梢的树被跳过并告警。枝长全部缺失的树按 0.0 处理并告警。
- patristic 输出为每棵树 O(n²)。写入前先按 `--max-rows`（默认 5,000,000；
  0 禁用限制）检查估计行数，并记录一条估计警告。
- `Bio.Phylo` 可能把任意文本解析为退化的单末梢树；此类输入发出“少于两个末梢”
  警告而非解析失败警告。

## 说明

- 枝长单位为每位点替代数（substitutions per site）；它们不能独立区分经过的
  时间与替换速率。
- 多树文件以 `filename:index`（零基）标识；单树文件保留裸文件名。
- 批量处理会显示 transient Rich 进度条（"Processing trees"），除非
  `--quiet`/`--dry-run`。快速运行在可见前即完成；大型后验批（尤其 patristic）
  会稳定显示进度。`--threads` 并行化 `--tree-dir`（非 patristic）处理；单个
  `--tree` 文件串行执行。
- `result.json` 中记录的 `command` 为可复现调用：恒含必填输入与 `-o`，外加
  显式提供的选项与非默认值；默认值标志省略。`params` 始终携带全部解析值。
- 本命令仅提取诊断测量值，不判定哪个模型更优。请跨模型运行比较分布。
- 系统误差的解释、模型/类群敏感性选择及可选的后验预测模拟见
  [系统误差工作流参考](../../skills/phyloai-workflow/references/syserror-workflow.md)。
  使用后验树分布时，应先检查收敛并通常准备去除适当 burn-in 后的树输入；
  `brlen` 不选择 burn-in，也不对 treelist 做 thinning/filtering。
