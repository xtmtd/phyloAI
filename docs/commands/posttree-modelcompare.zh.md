# phyloai posttree modelcompare

[English](posttree-modelcompare.md) | [中文](posttree-modelcompare.zh.md)


## 目的

通过两个独立子命令进行相对替代模型比较与选择：

| 子命令 | 分析内容 | 核心工具 |
|--------|----------|----------|
| `iqtree` | ModelFinder BIC/AIC/AICc 模型比较 | IQ-TREE3 `-m MF` |
| `pb` | 留一交叉验证（LOO-CV）/ wAIC 模型比较 | PhyloBayes `.sitelogl`（纯 Python） |

`iqtree` 使用 ModelFinder（Kalyaanamoorthy et al. 2017）对单条比对进行替代模型比较，并报告 BIC/AIC/AICc 分数、权重及 95% 置信集归属。`pb` 使用 LOO-CV 和广泛适用信息准则（wAIC，Lartillot 2023）根据位点对数似然文件比较拟合模型的预测准确性。

## 用法

```bash
# 在同质模型搜索空间上进行 ModelFinder 比较（BIC/AIC/AICc）
phyloai posttree modelcompare iqtree --matrix ./matrix.aa.fa --homogeneous-model LG,WAG --mrate E,G,R

# 加入异质混合模型（仅限 AA）
phyloai posttree modelcompare iqtree --matrix ./matrix.aa.fa --homogeneous-model LG --heterogeneous-model C10,C20 --het-mrate G,R

# 对两个拟合模型进行 LOO-CV / wAIC 比较（每个目录代表一个模型）
phyloai posttree modelcompare pb --sitelogl-dir ./cat_sitelogl,./gtr_sitelogl

# 显式指定链分组（重复 --sitelogl；每次出现代表一个模型）
phyloai posttree modelcompare pb --sitelogl model1/c1.sitelogl,model1/c2.sitelogl --sitelogl model2/c1.sitelogl,model2/c2.sitelogl

# 单模型拟合报告（一个目录）
phyloai posttree modelcompare pb --sitelogl-dir ./gtr_sitelogl --model-names GTR
```

## modelcompare iqtree — ModelFinder 模型比较

### 目的

运行 IQ-TREE3 `-m MF`，指定同质模型搜索空间（`-mset`）与速率异质性类型（`-mrate`），并可选择通过 `-madd` 展开异质混合模型。从 `.iqtree` 报告的 “List of models sorted by BIC scores:” 部分解析出比较表。

### 输入

| 输入 | 说明 |
|------|------|
| `--matrix` | 单一超矩阵比对（FASTA、PHYLIP、NEXUS）。必填。对应 IQ-TREE `-s`。 |
| `--homogeneous-model` | 同质搜索空间的逗号分隔标准模型。必填。对应 IQ-TREE `-mset`。 |
| `--mrate` | 同质模型的速率异质性类型。有效值：`E`、`G`、`R` 的任意非空子集（逗号分隔）。默认 `E,G`。对应 IQ-TREE `-mrate`。 |
| `--heterogeneous-model` | 逗号分隔的 AA 混合模型（`C10`–`C60`、`EX*`、`EHO`、`UL*`、`EX_EHO`、`LG4M`、`LG4X`），通过 `-madd` 评估。仅限 AA；对 NT 数据会报错。 |
| `--het-mrate` | 异质模型展开的速率异质性。每个 token 选择一个变体族：`E` = 基础模型（`C10, C10+F`）、`G` = `+G4`、`R` = `+R4`。有效值：`E`、`G`、`R` 的任意子集（逗号分隔）。默认 `E,G`。仅限 AA。 |
| `--seq-type` | `AA`、`NT` 或 `auto`（默认 `auto`）。为 `auto` 时，会在模型校验前读取比对以检测 AA 或 NT。 |
| `--prefix` | IQ-TREE 输出前缀（默认：`modelcompare`）。必须是单个文件名——路径分隔符、`..` 和绝对路径均被拒绝，以确保输出留在运行目录内。 |
| `--threads` | IQ-TREE `-T` 值（整数或 `auto`，默认 `auto`）。 |
| `--iqtree-path` | iqtree3 可执行文件的显式路径。 |
| `--tool-args` | 额外的 IQ-TREE 参数。被阻止：`-s`、`--prefix`。若包含 PhyloAI 也管理的参数（`-m`、`-mset`、`-mrate`、`-madd`、`-cmin`、`-cmax`、`-T`），`--tool-args` 中的值将覆盖 PhyloAI 生成的参数（命令中不出现重复）。 |
| `--output-dir` | 输出目录（默认：`runs/posttree/modelcompare/iqtree`）。 |
| `--overwrite` | 运行前删除并重建输出目录。 |
| `--resume` | 从 IQ-TREE 原生检查点恢复未完成的任务。 |
| `--dry-run` | 打印 IQ-TREE 命令而不执行。 |
| `--quiet` | 除错误外抑制终端输出。 |

### 异质模型展开算法

`--het-mrate` 中的每个 token 为 `--heterogeneous-model` 中的每个模型 M 选择一个变体族，与 `--mrate` 语义一致：

- `E` → `M, M+F`（经验状态频率，无速率类别）
- `G` → `M+G4, M+F+G4`
- `R` → `M+R4, M+F+R4`

仅产生被请求的变体族。例如，`--heterogeneous-model C10 --het-mrate E,G` 产生：

```
C10, C10+F, C10+G4, C10+F+G4
```

而单独 `--het-mrate G` 只产生 `C10+G4, C10+F+G4`（不含基础 `C10`）。

展开后的列表以逗号连接，传给 `-madd`。

### 输出

```
runs/posttree/modelcompare/iqtree/
├── result.json
├── model_fit.csv                 # Rank,Model,LogL,AIC,w_AIC,In_AIC_95,AICc,w_AICc,In_AICc_95,BIC,w_BIC,In_BIC_95
└── iqtree/
    ├── modelcompare.iqtree       # IQ-TREE 原生报告（BIC/AIC/AICc 表）
    ├── modelcompare.log
    ├── modelcompare.model.gz
    └── modelcompare.treefile
```

`model_fit.csv` 按 BIC 排序。`In_AIC_95` / `In_AICc_95` / `In_BIC_95` 表示各准则 95% 置信集归属（IQ-TREE 输出中的 `+`）。

### 示例

```bash
phyloai posttree modelcompare iqtree --matrix concat.aa.fa --homogeneous-model LG,WAG --mrate E,G,R --heterogeneous-model C10,C20 --het-mrate G,R
```

---

## modelcompare pb — LOO-CV / wAIC 模型比较

### 目的

使用留一交叉验证（LOO-CV）和广泛适用信息准则（wAIC，Watanabe 2009）比较拟合模型的预测准确性，依据 Lartillot (2023)。完全在 Python 中由 PhyloBayes `.sitelogl` 位点对数似然文件计算（无需外部工具）。

### 输入

| 输入 | 说明 |
|------|------|
| `--sitelogl-dir` | 逗号分隔的目录；每个目录代表一个模型，扫描 `*.sitelogl`（每目录 ≥2 个文件）。与 `--sitelogl` 互斥。支持 shell 路径补全。 |
| `--sitelogl` | 可重复选项；每次出现是单个模型的一组逗号分隔 `.sitelogl` 文件路径（每模型 ≥2 个文件）。与 `--sitelogl-dir` 互斥。支持 shell 路径补全。 |
| `--model-names` | 与模型组数量一致的逗号分隔模型标签。省略时模型命名为 `model_1`、`model_2` 等。每个标签必须是单个路径组件（不含 `/`、`..`），因为它作为 `sitelogl/` 下的输出子目录名。 |
| `--output-dir` | 输出目录（默认：`runs/posttree/modelcompare/pb`）。 |
| `--overwrite` | 运行前删除并重建输出目录。 |
| `--quiet` | 除错误外抑制终端输出。 |

校验规则：
- 必须且只能提供 `--sitelogl-dir` 与 `--sitelogl` 之一。
- 至少需要 1 个模型组；每个模型组至少 2 个 `.sitelogl` 文件。
- 同一模型组内所有 `.sitelogl` 文件的数据行数必须相同。
- **跨模型位点校验：** 当提供 ≥2 个模型组时，所有组必须具有相同的位点数量与完全一致的按序 `site` 标识。来自不同比对或位点顺序的分数不可比较；不一致会硬报错。
- 提供 `--model-names` 时，标签数量必须与模型组数量一致，标签必须唯一，且每个标签必须是单个路径组件（不含 `/`、`..`）。
- 组内重复的文件名会追加数字后缀消歧（`chain1_1.sitelogl`、`chain1_2.sitelogl`）。

### 计算方法

每个模型（Lartillot 2023）：
- **LOO-CV：** 先对每个 run 求位点 `logcpo` 均值，再对 run 间取均值；去偏为 `分数 − 0.5 × 位点间 logcpo 方差（跨 run）的均值`。
- **wAIC：** 每个 run 为 `mean(logpostmeanl) − mean(var)`，再对 run 间取均值；去偏为 `分数 + 0.5 × 位点间 logpostmeanl 方差（跨 run）的均值`。
- **ESS 质量：** `%(ess<10)`（ESS < 10 的位点占比）与 `f(ess<10)`（此类位点贡献的分数占比）。质量为 `good`（max < 0.1）、`ok`（max < 0.3）或 `no`（max ≥ 0.3）。
- **置信区间：** 使用 `n_runs − 1` 自由度的 Student's t 临界值（df 1–30 用精确表值；30 以上线性插值逼近 1.96）。
- **Δ 值（≥2 个模型）：** 每个指标独立选出最优模型（去偏分数最高）；Δ = 模型分数 − 最优分数，因此最优模型 Δ = 0，其余 ≤ 0。

### 输出

```
runs/posttree/modelcompare/pb/
├── result.json
├── model_fit.csv
└── sitelogl/
    ├── model_1/
    │   ├── chain1.sitelogl
    │   └── chain2.sitelogl
    └── model_2/
        ├── chain1.sitelogl
        └── chain2.sitelogl
```

单模型时 `model_fit.csv` 使用按指标格式（Metric, Score, Bias, StDev, CI95_min, CI95_max, ESS, Pct_ESS_lt10, Frac_ESS_lt10, Quality）；多模型时使用按模型宽格式，含 Delta_LOOCV 与 Delta_wAIC。

`result.json` 的 `key_results` 同时暴露 `best_loocv_quality` 与 `best_waic_quality`（各为 `good` / `ok` / `no`）。

### 示例

```bash
# 两个模型目录，各含 >= 2 个链文件
phyloai posttree modelcompare pb --sitelogl-dir model1,model2 --model-names CAT,GTR

# 显式链分组
phyloai posttree modelcompare pb \
  --sitelogl model1/c1.sitelogl,model1/c2.sitelogl \
  --sitelogl model2/c1.sitelogl,model2/c2.sitelogl
```

---

## 通用说明

- `--seq-type auto` 会在模型校验*之前*读取比对以检测 AA 或 NT；`--heterogeneous-model` 对 NT 数据会报错。显式指定 `--seq-type AA|NT` 时会与实际检测类型交叉校验，不一致则报错。
- `iqtree` 的 IQ-TREE 输出文件放在 `iqtree/` 子目录下；执行时 IQ-TREE stdout 实时输出到终端。
- `--dry-run` 打印 IQ-TREE 命令并校验输入，但不运行外部工具。
- 输出目录必须为空（或使用 `--overwrite`）；`pb` 在未指定 `--overwrite` 时会拒绝非空目录。
- 参考文献：Kalyaanamoorthy et al. (2017) *Nature Methods*；Lartillot (2023) *Systematic Biology* 72(3):616–638；Watanabe (2009) *JMLR*。

## 退出码

| 代码 | 含义 |
|------|------|
| 0 | 成功 |
| 1 | 用户输入错误（文件缺失、参数无效、输出冲突） |
| 2 | 外部工具执行失败 |
| 3 | 未找到外部工具可执行文件 |
