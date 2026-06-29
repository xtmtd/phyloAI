# phyloai posttree dating

[English](posttree-dating.md) | [中文](posttree-dating.zh.md)

## 目的

基于 MCMCtree 的两步法贝叶斯分子定年，使用近似似然（`usedata=2`）：

1. **`hessian`** —— 用 IQ-TREE3 `--dating mcmctree` 计算梯度与 Hessian。
2. **`mcmc`** —— 使用 IQ-TREE 输出的 MCMCtree，并行运行独立的 posterior + prior 链，提供实时进度与诊断图。

仅实现近似似然（`usedata=2`）。精确似然（`usedata=1`，使用序列数据）与 `usedata=3`（输出 `in.BV`）在 PhyloAI 中未实现。

## 快速开始

```bash
# Step 1：计算 Hessian（AA，自动模型 LG+F+G4）
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre

# Step 2：运行 MCMCtree（2 条独立链，独立速率分子钟）
phyloai posttree dating mcmc --hessian-dir runs/posttree/dating/hessian

# Step 3：查看 runs/posttree/dating/mcmc/diagnostics/ 下的诊断结果
```

## 定年树格式

`--rooted-tree` 必须是带 MCMCtree 风格定年注解、并带有根节点年龄约束的有根 NEWICK 树。定年单位为 100 Mya（例如 `'>3.1<3.8'` 表示 310–380 Mya）。

```
# 单下界约束（C-D 支为 310 Mya）
(A,((B,(C,D)'>3.1<3.8'),(E,F)'>2.9<3.6'))'<4.2';

# 点定年 + 软约束（均匀分布）
(A,((B,(C,D)'>0.5<0.7'),(E,F)'>0.3<0.5'))'<1.0';

# Tip dating：节点标签标注现存分类单元
(A_tip,((B,(C,D_tip)'>3.1<3.8'),E_tip)'<4.2')'<5.0';
```

约束使用 `>`、`<` 或两者组合：`'>L<U'` 定义下界 `L` 与上界 `U`。根节点（最外层节点）必须带约束——PhyloAI 会拒绝缺少根约束的树。

## Hessian 步骤

计算 `usedata=2`（近似似然）所需的梯度与 Hessian 矩阵，并产出三个供 `mcmc` 步骤消费的文件。

### 用法

```bash
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre [OPTIONS]
```

### 输入

| Flag | Default | Description |
|------|---------|-------------|
| `--matrix PATH` | *(required)* | 比对文件：FASTA（`.fa`/`.fas`/`.fasta`/`.faa`/`.fna`/`.aln`）、PHYLIP 或 NEXUS。`--seq-type auto` 通过 PhyloAI 的统一格式检测器识别三种格式。 |
| `--rooted-tree PATH` | *(required)* | MCMCtree 定年树（见上）。对应 IQ-TREE `-te`。 |
| `--seq-type AA\|NT\|auto` | `auto` | 序列类型。`AA` → `LG+F+G4`；`NT` → `GTR+G4`。`auto` 根据比对内容自动检测。 |
| `--model-expr STR` | — | IQ-TREE 模型字符串（如 `C10+F+G4`）。与 `--partitions` 互斥。 |
| `--partitions PATH` | — | 分区文件（RAxML 风格、NEXUS `.best_model.nex` 或 cluster 文件）。分区数 < 10 时直接运行；≥ 10 时自动用 `--merge --rclusterf 10` 合并。对应 `-Q`。 |
| `-o, --output-dir PATH` | `runs/posttree/dating/hessian` | 输出目录。 |
| `-t, --threads INT` | `4` | IQ-TREE 线程数（`-T`）。 |
| `--iqtree-path PATH` | — | 覆盖自动检测到的 `iqtree3` 二进制路径。 |
| `--tool-args STR` | — | 额外的 IQ-TREE 参数。被阻止的：`-s`、`--dating`、`-te`、`--prefix`。 |
| `--overwrite` | off | 删除并重建输出目录。 |
| `--resume` | off | 恢复中断的 IQ-TREE 运行（IQ-TREE 原生 checkpoint）。 |
| `--dry-run` | off | 仅打印 IQ-TREE 命令而不执行。 |
| `-q, --quiet` | off | 除错误外不打印终端输出。 |

### 输出

```
runs/posttree/dating/hessian/
├── result.json
├── iqtree.dummy.phy          → mcmc seqfile
├── iqtree.rooted.nwk         → mcmc treefile
└── iqtree.mcmctree.hessian   → 在每个 run 目录中被重命名为 in.BV
```

`iqtree.*` 输出文件使用固定前缀——PhyloAI 在 `--tool-args` 中阻止 `--prefix`，因为 `mcmc` 步骤依赖这些确切的文件名。

### 模型选择逻辑

| 条件 | 模型 |
|------|------|
| 设置了 `--model-expr` | 直接使用该表达式（`-m <expr>`） |
| 设置了 `--partitions`，AA | `-m MF -Q <file> --mset LG -mfreq F -mrate G` |
| 设置了 `--partitions`，NT | `-m MF -Q <file> --mset GTR -mrate G` |
| 不分区，AA | `-m LG+F+G4` |
| 不分区，NT | `-m GTR+G4` |

当 `--partitions` 与 ModelFinder（`-m MF`）联用时，模型搜索空间被约束到定年场景下最简单合适的模型族——不含自由速率模型（`+R`）或会显著增加计算时间的复杂混合模型。

## MCMC 步骤

使用 hessian 输出运行 MCMCtree 贝叶斯定年。并行启动独立 posterior 链，每条与一条 prior-predictive 链（`usedata=0`）配对，使用相同随机种子。

### 用法

```bash
phyloai posttree dating mcmc --hessian-dir runs/posttree/dating/hessian [OPTIONS]
```

### 输入

| Flag | Default | Description |
|------|---------|-------------|
| `--hessian-dir PATH` | *(required)* | 包含 `iqtree.dummy.phy`、`iqtree.rooted.nwk`、`iqtree.mcmctree.hessian` 的目录。 |
| `-o, --output-dir PATH` | `runs/posttree/dating/mcmc` | 输出目录。 |
| `--runs INT` | `2` | 独立 posterior 链数（每条与一条 prior 链配对）。 |
| `--clock 1\|2\|3` | `2` | 分子钟模型：`1` 全局；`2` 独立速率；`3` 相关速率。提供 `--ctl` 时被忽略。 |
| `--burnin INT` | `100000` | MCMC burn-in 迭代数。提供 `--ctl` 时被忽略。 |
| `--sample-freq INT` | `10` | 采样频率（每 N 次迭代采一个样本）。提供 `--ctl` 时被忽略。 |
| `--nsamples INT` | `10000` | burn-in 后保留的样本数。提供 `--ctl` 时被忽略。 |
| `--ctl PATH` | — | 预配置的 `mcmctree.ctl`。与非默认 `--clock`/`--burnin`/`--sample-freq`/`--nsamples` 互斥。 |
| `--mcmctree-path PATH` | — | 覆盖自动检测到的 `mcmctree` 二进制路径。 |
| `--overwrite` | off | 删除并重建输出目录。 |
| `--dry-run` | off | 仅生成 ctl 并退出，不运行。 |
| `-q, --quiet` | off | 不输出 MCMC 日志尾部与进度条。 |

**总迭代数** = `--burnin` + (`--sample-freq` × `--nsamples`)。默认 200,000。

### 分子钟模型指南

| Clock | Model | 适用场景 |
|-------|-------|---------|
| `1` | Global | 所有支系速率相同（罕见成立）。最快、最简单。 |
| `2` | Independent（推荐） | 每条支有独立速率，从对数正态分布采样。多数数据集的良好默认。 |
| `3` | Correlated | 相邻支速率相关。更具生物学真实性但计算更重。 |

### 输出

```
runs/posttree/dating/mcmc/
├── result.json
├── mcmctree.ctl                   # 生成的模板（或 --ctl 的副本）
├── run1/
│   ├── mcmctree.ctl               # posterior ctl，注入种子
│   ├── mcmc.txt                   # MCMC 参数轨迹
│   ├── mcmctree.out               # 分歧时间汇总
│   ├── mcmctree.log
│   ├── FigTree.tre                # 带注解的定年树
│   ├── FigTree.node.tre           # 节点标签树（从 mcmctree.out 解析）
│   ├── SeedUsed                   # MCMCtree 写入的种子
│   ├── iqtree.dummy.phy -> <hessian-dir>/
│   ├── iqtree.rooted.nwk -> <hessian-dir>/
│   ├── in.BV -> <hessian-dir>/iqtree.mcmctree.hessian
│   └── prior/
│       ├── mcmctree.ctl           # usedata=0，与 posterior 同种子
│       ├── mcmc.txt
│       ├── mcmctree.out
│       └── FigTree.node.tre
├── run2/                          # 结构相同
│   └── ...
└── diagnostics/
    ├── traces/                    # 各参数的 MCMC 轨迹 PDF
    ├── convergence/
    │   ├── posterior_times.csv    # 所有链合并
    │   ├── prior_times.csv
    │   └── convergence_*_runX_vs_runY.pdf  # 散点图 + 拟合线
    ├── infinite_sites/            # 平均年龄 vs 95% CI 宽度图
    ├── posterior_vs_prior/        # 各节点 posterior vs prior 平均年龄
    └── spearman_correlations.csv
```

### 诊断解读

| 诊断 | 看什么 |
|------|--------|
| **Traces** (`traces/`) | 参数轨迹应混合良好（"毛毛虫"形态），无趋势或停滞。早期收敛、保持稳定区间。 |
| **Convergence** (`convergence/`) | 独立链的平均分歧时间应沿 y=x 排列。斜率接近 1.0、Spearman ρ 接近 1.0、RMSE 低。系统性偏离提示 burnin 或链长不足。 |
| **Infinite-sites** (`infinite_sites/`) | CI 宽度不应随节点年龄显著增加。强正斜率意味着老节点约束差——检查化石定年。 |
| **Posterior vs Prior** | posterior 强烈偏离 prior 意味着数据有信息；若多数节点 posterior ≈ prior，数据携带的时间信号弱。 |
| **Spearman correlations** (`spearman_correlations.csv`) | 报告每对比较的 ρ、p 值、线性拟合（斜率、截距、RMSE）。用作快速收敛检查。 |

## 示例

```bash
# ── Hessian ───────────────────────────────────────

# 1. 不分区 AA，默认模型（LG+F+G4）
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre

# 2. 自定义混合模型
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre --model-expr C10+F+G4

# 3. 显式指定 NT 序列类型
phyloai posttree dating hessian --matrix concat.nt.fa --rooted-tree calib.tre --seq-type NT

# 4. 分区分析（< 10 个分区）
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre \
    --partitions partitions.nex -o runs/dating/hessian

# 5. 自定义输出目录
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre \
    -o runs/dating/hessian

# 6. 恢复中断的运行（IQ-TREE 原生 checkpoint）
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre --resume

# ── MCMC ──────────────────────────────────────────

# 7. 默认：2 条链，独立速率分子钟
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian

# 8. 三条链，相关分子钟，更长链长
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian \
    --runs 3 --clock 3 --burnin 200000 --nsamples 20000

# 9. Dry-run：运行前查看生成的 ctl
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian --dry-run

# 10. 使用预配置的 mcmctree.ctl 完全控制
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian --ctl my_run.ctl

# 11. 单链（不生成收敛诊断）
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian --runs 1

# 12. 覆盖已有结果
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian --overwrite
```

## 退出码

| Code | Meaning |
|------|---------|
| `0` | 成功 |
| `1` | 输入错误（文件缺失/为空、参数冲突、输出目录非空） |
| `2` | 工具失败（IQ-TREE 非零退出、MCMCtree posterior 失败） |
| `3` | 环境错误（未找到 `iqtree3` 或 `mcmctree`） |

## 警告 / 错误

| 条件 | 行为 |
|------|------|
| `hessian-dir` 为空或缺失 | Exit 1，提示缺失的文件 |
| `--matrix` 或 `--rooted-tree` 不存在或为空 | 在 IQ-TREE 启动前 Exit 1 |
| `--rooted-tree` 缺少根节点年龄约束 | Exit 1；最外层节点必须有 `'>L<U'` 或 `'<U'` |
| `--model-expr` 与 `--partitions` 同时使用 | Exit 1；互斥 |
| `--ctl` 与非默认 `--clock`/`--burnin`/`--sample-freq`/`--nsamples` | Exit 1；这些参数仅在从头生成 ctl 时生效 |
| 输出目录非空且未加 `--overwrite` | Exit 1；用 `--overwrite` 覆盖 |
| IQ-TREE 返回非零 | `result.json` 中 status 为 `error`；stderr 记录在 `data.warnings` |
| IQ-TREE 产出空文件 | 记为警告；可能是写入过程中崩溃 |
| IQ-TREE 报告缺少 "Total CPU time used" | 警告；IQ-TREE 可能被中断 |
| 未找到 `mcmctree` 二进制 | Exit 3（环境错误）。运行 `phyloai doctor`。 |
| posterior 退出码 ≠ 0 | 顶层 `status` 变为 `"error"` |
| prior 退出码 ≠ 0 或文件缺失 | 记录在 `data.warnings`；posterior 结果与诊断仍生成 |
| `mcmc.txt` 或 `mcmctree.out` 为空 | 记录在 `data.warnings`；诊断优雅降级 |
| `--runs=1` | 跳过收敛图；其他诊断仍生成 |

## 备注

- **仅 usedata=2。** PhyloAI 未实现 `usedata=1`（基于序列数据的精确似然）和 `usedata=3`（输出 `in.BV`）。
- **种子注入。** PhyloAI 为每条链生成唯一随机种子（`random.randint(1, 2³¹-1)`），注入 ctl 后启动。同一 run 内的 posterior 与 prior 共享种子，便于比较。
- **ndata。** 数据块数始终直接从 `iqtree.dummy.phy` 计数——而非 `result.json`。这保证了 IQ-TREE 在合并分区（`--merge --rclusterf 10`）时的正确性。
- **seqtype。** 从 `hessian/result.json`（`params.seq_type`）读取；若 result.json 缺失或损坏，则从 dummy.phy 内容自动检测。
- **版本检测。** `mcmctree` 二进制通过 `shutil.which` 查找；版本通过无参运行二进制并匹配 `paml version (\d+(?:\.\d+)+)` 提取。运行完成后回退解析 `mcmctree.log`。
- **前缀锁定。** PhyloAI 在 hessian 步骤硬编码 `--prefix iqtree`。这样 mcmc 步骤无需用户追踪命名就能找到 `iqtree.dummy.phy`、`iqtree.rooted.nwk`、`iqtree.mcmctree.hessian`。`--prefix` 在 `--tool-args` 中被阻止。
- **自定义 ctl (`--ctl`)。** PhyloAI 始终将标准 hessian 文件（`iqtree.dummy.phy`、`iqtree.rooted.nwk`、`in.BV`）软链到每个 `runN/` 目录。ctl 中自定义的相对 `seqfile`/`treefile` 路径相对 ctl 所在目录解析，并软链到每个 `runN/`。
- **诊断。** 收敛图要求 `--runs >= 2`。所有散点图使用一条虚线拟合线带方程；节点标签使用 `nXX` 格式。posterior-vs-prior 图使用等长坐标轴便于直接对比。
- **`result.json`** 遵循 PhyloAI JSON Output Standard。
- **Prior-predictive 链。** Prior `mcmctree.ctl` 由 posterior ctl 通过将 `usedata = 2` 改为 `usedata = 0` 派生，保留相同种子。其他参数完全一致，确保 prior 与 posterior 可直接比较。
- **OMP_NUM_THREADS。** 对 MCMCtree 子进程始终设为 `1`。MCMC 采样器不能从多线程中受益，`OMP_NUM_THREADS > 1` 时甚至会产生错误结果。