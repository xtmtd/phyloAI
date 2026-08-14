# phyloai posttree syserror cca

[English](posttree-syserror-cca.md) | [中文](posttree-syserror-cca.zh.md)

## 目的

执行组成约束分析（Compositional Constraint Analysis, CCA）：这是一个本地组成约束诊断，比较两个模型分析中位点对候选拓扑的偏好。CCA 对每个位点由 20 个氨基酸的位点频率计算有效氨基酸数（Keff），并对每个模型独立计算 Tree2 减 Tree1 的位点对数似然差。

CCA 仅作诊断。它不生成位点频率或似然值、不推断树、不筛选位点、不选择模型，也不能确定真实拓扑。通常先用 IQ-TREE PMSF 或 `phyloai tree bi readpb --mode ss` 准备位点频率，并用两次 `phyloai posttree signal lnl` 运行准备似然表。

该命令为纯 Python，不调用外部可执行程序；无需 `phyloai doctor`、checkpoint 或 resume。

## 用法

```bash
phyloai posttree syserror cca \
  --site-freq chain1.sitefreq \
  --site-lnl1 lnl_LG/site_lnl.csv \
  --site-lnl2 lnl_C20/site_lnl.csv \
  [--model1-name model1] [--model2-name model2] \
  [--title TEXT] \
  [--xlabel "Effective number of amino acids"] \
  [--ylabel "Log-likelihood difference"] \
  [--fig-width 10] [--fig-height 6] [--dpi 300] [--font-size 16] \
  [-o runs/posttree/syserror/cca] [--overwrite] [--dry-run] [-q]
```

## 输入

| 选项 | 必填 | 默认值 | 说明 |
|---|---|---|---|
| `--site-freq` | 是 | -- | IQ-TREE PMSF 或 PhyloAI 转换的 `readpb --mode ss` `.sitefreq` 表；不接受原始 PhyloBayes `.siteprofiles`。 |
| `--site-lnl1` | 是 | -- | 模型分析 1 的 `site_lnl.csv`；必须含 `site`、`lnL_Tree1`、`lnL_Tree2`。 |
| `--site-lnl2` | 是 | -- | 模型分析 2 的 `site_lnl.csv`；必须含 `site`、`lnL_Tree1`、`lnL_Tree2`。 |
| `--model1-name` | 否 | `model1` | 模型 1 的非空且不同于模型 2 的标签，用于 CSV、图例、JSON 和报告。 |
| `--model2-name` | 否 | `model2` | 模型 2 的非空且不同于模型 1 的标签，用于 CSV、图例、JSON 和报告。 |
| `--title` | 否 | 空 | 可选 PDF 标题。 |
| `--xlabel` | 否 | `Effective number of amino acids` | X 轴标签。 |
| `--ylabel` | 否 | `Log-likelihood difference` | Y 轴标签。 |
| `--fig-width` | 否 | `10` | 正的 PDF 宽度（英寸）。 |
| `--fig-height` | 否 | `6` | 正的 PDF 高度（英寸）。 |
| `--dpi` | 否 | `300` | 正的 PDF 栅格化 DPI 元数据。 |
| `--font-size` | 否 | `16` | 图例文字大小（pt）。 |
| `-o`, `--output-dir` | 否 | `runs/posttree/syserror/cca` | 输出目录。 |
| `--overwrite` | 否 | false | 删除并重建非空输出目录。 |
| `--dry-run` | 否 | false | 校验输入并计算结果载荷，但不写文件。 |
| `-q`, `--quiet` | 否 | false | 除错误外抑制终端输出。 |

`.sitefreq` 的每个数据行必须含 1 基位点 ID 和恰好 20 个有限、非负、和为 1（容差 `1e-6`）的频率。位点 ID 必须唯一且严格连续为 `1..N`。两个 LNL CSV 都必须有字面量表头 `site`、`lnL_Tree1`、`lnL_Tree2`；`ΔSLS`、`support` 等额外列会被忽略。行的顺序可以不同，但三个输入的完整 1 基连续位点集合必须完全一致。

## 计算与 CSV

对每个位点计算：

```text
Keff = 1 / sum(p_i^2), i = 1..20
delta_lnl_tree2_tree1 = lnl_tree2 - lnl_tree1
```

CCA 会重新计算第二个表达式，不会使用 `site_lnl.csv:ΔSLS`；后者在 signal LNL 中的符号是 Tree1 减 Tree2。CCA 值为正表示支持 Tree2，负表示支持 Tree1。

`cca.csv` 固定使用以下 ASCII snake_case 列，并按位点升序、每个位点模型 1 后模型 2 排序：

```csv
model,site,keff,lnl_tree1,lnl_tree2,delta_lnl_tree2_tree1
LG,1,11.974845235298696,-14.2296,-14.3580,-0.1284
C20,1,11.974845235298696,-13.8521,-13.9077,-0.0556
```

CSV 中 `lnl_tree1`、`lnl_tree2` 和 `delta_lnl_tree2_tree1` 固定保留 4 位小数；`keff` 保留完整浮点精度。下例是历史培训 `cca.txt` 参考；当前随附 `site_lnl.csv` 输入的 site 1 delta 分别为 LG `-0.0999`、C20 `-0.0436`。

它与培训 `cca.txt` 字段一一对应：`keff` → `Keff`、`lnl_tree1` → `LnL_T1`、`lnl_tree2` → `LnL_T2`、`delta_lnl_tree2_tree1` → `deltaLnL_T2_T1`。

## 图形

CCA 按 `floor(Keff)` 为每个模型分箱，汇总每箱 CCA 差值并输出 PDF。分箱固定为 1 至 20；缺失的模型/箱组合补零；Keff 为 20 时属于第 20 箱。

图形遵循培训 ggplot 语义：成对柱中心为 `bin + 0.5`，width/dodge 为 1。绘图区延至 21、但 ticks 仍为 1–20，从而完整显示有效的 Keff=20 分箱。与 ggplot 默认离散填充色标一致，模型标签按字母排序后依次分配 `#F8766D`、`#00BFC4`：对培训中的 `C20`/`LG`，C20 为橙红色，LG 为蓝绿色。X 轴 break 为 1–20、无扩展，并关闭 x 主网格；保留 y 主网格。1–20 的竖向分箱边界为灰色实线 0.1 pt，零基线为黑色。Y 轴上下限为 `min(0, 1.1 * minimum_bin_sum)` 和 `max(0, 1.1 * maximum_bin_sum)`；负值区域为 `#ffdab9`、正值区域为浅蓝色，均为 50% 透明度。右上图例锚定于 `(0.99, 0.9)`、无标题、白色半透明背景、黑色 0.5 pt 边框，并使用设置的图例字号；其余文字遵循约 11 pt 的 `theme_bw()` 基线。

## 输出

```text
runs/posttree/syserror/cca/
├── cca.csv
├── cca.pdf
└── result.json
```

不会写出 PNG 或重复的汇总表。`result.json` 记录所有解析后的参数、Keff 与各模型 delta 汇总、20 个作图分箱汇总，以及 CSV/PDF 的绝对路径和描述。

## 示例

```bash
# 用共享位点频率比较 LG 与 C20 似然分析
phyloai posttree syserror cca \
  --site-freq chain1.sitefreq \
  --site-lnl1 lnl_LG/site_lnl.csv --site-lnl2 lnl_C20/site_lnl.csv \
  --model1-name LG --model2-name C20

# 校验输入和计算，但不创建输出目录
phyloai posttree syserror cca --site-freq chain1.sitefreq \
  --site-lnl1 lnl1/site_lnl.csv --site-lnl2 lnl2/site_lnl.csv --dry-run
```

## 警告 / 错误

- 所有输入路径必须是可读取的普通文件。
- 每个 site frequency 行必须有恰好 20 个有效氨基酸频率；原始 `.siteprofiles` 和 0 基位点 ID 会因格式校验失败。
- 每个 LNL 表必须使用精确的必需表头；自定义拓扑标签的列名不被接受。
- 空输入、格式错误、重复、非整数、非有限值或非连续位点/数值，以及任一输入间的位点集合不匹配，都会报错。
- 模型标签必须非空且彼此不同；图尺寸、DPI、图例字号必须为正。
- 非空输出目录需要 `--overwrite`。dry-run 不写文件；验证在 overwrite 删除前完成。若带 `--overwrite` 的验证失败，保留已有文件，但将根 `result.json` 替换为输入错误记录；没有 resume/checkpoint。

## 说明

- Keff 是 20 个频率的 inverse homozygosity，不是熵；CCA 不实现熵版本。
- CCA 可显示拓扑偏好是否随组成约束而改变，但不能单独证明任一模型或拓扑在生物学上正确。
