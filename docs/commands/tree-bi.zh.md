# phyloai tree bi

[English](tree-bi.md) | [中文](tree-bi.zh.md)

使用 [PhyloBayes-MPI](https://github.com/bayesiancook/pbmpi) 进行贝叶斯系统发育推断。

## 概述

`phyloai tree bi` 是一个包含四个子命令的 Click 组：

| 子命令 | 用途 |
|---|---|
| `tree bi pb` | 使用 `pb_mpi` 运行 MCMC 链推断 |
| `tree bi bpcomp` | 使用 `bpcomp` 进行拓扑收敛分析 |
| `tree bi tracecomp` | 使用 `tracecomp` 进行参数收敛分析 |
| `tree bi readpb` | 使用 `readpb_mpi` 进行后验分析和预测检验 |

默认输出根目录为 `runs/tree/bi/`。

## 要求

每个子命令仅从 `PATH`（或 `--pb-path` 指定的目录）解析其所需的 PhyloBayes-MPI 工具：

| 工具 | 所需子命令 | 用途 |
|------|-----------|------|
| `pb_mpi` | `bi pb` | MCMC 采样器 |
| `bpcomp` | `bi pb`, `bi bpcomp` | 拓扑收敛 |
| `tracecomp` | `bi pb`, `bi tracecomp` | 参数收敛 |
| `mpirun` | `bi pb`, `bi readpb` | Open MPI 启动器 |
| `readpb_mpi` | `bi readpb` | 后验分析 |

运行 `phyloai doctor` 确认安装。

---

# phyloai tree bi pb

并行运行 N 条独立的 MCMC 链，实时监控收敛，生成一致树。

## 用法

```bash
phyloai tree bi pb --matrix <比对文件> [OPTIONS]
```

## 示例

```bash
# 默认：3 条链，CAT-GTR 模型，无限运行
phyloai tree bi pb --matrix concat/matrix.phy

# 同质 LG+G4，10000 个 MCMC 总循环后停止
phyloai tree bi pb --matrix concat/matrix.phy --model lg --mixture 1 --nsamples 10000

# 向现有运行添加两条新链
phyloai tree bi pb --matrix concat/matrix.phy --chain-names chain4,chain5 -o runs/tree/bi

# 恢复所有链
phyloai tree bi pb -o runs/tree/bi --resume

# 恢复指定链
phyloai tree bi pb -o runs/tree/bi --resume chain1,chain3

# 自定义 PhyloBayes 工具目录
phyloai tree bi pb --matrix concat/matrix.phy --pb-path /opt/pbmpi/bin

# 干运行
phyloai tree bi pb --matrix concat/matrix.phy --dry-run
```

## 参数

### 输入 / 输出

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--matrix / -m` | Path | 必需 | 输入比对（PHYLIP 或 FASTA）。使用 `--resume` 时不需要。 |
| `--output-dir / -o` | Path | `runs/tree/bi` | 输出目录。 |
| `--overwrite` | flag | False | 删除并重建输出目录。与 `--resume` 互斥。 |

### 模型

| 参数 | 可选值 | 默认值 | pb_mpi 标志 | 说明 |
|------|--------|--------|-------------|------|
| `--model` | `gtr`, `poisson`, `lg`, `wag`, `jtt`, `mtrev`, `mtzoa`, `mtart` | `gtr` | `-gtr`, `-poisson`, … | 速率矩阵。 |
| `--mixture` | `auto`、`1` 或整数 N | `auto` | `-cat` / 无混合标志 / `-ncat N` | `auto` = CAT Dirichlet 过程；`1` = 单矩阵同质模型；整数 N = 固定 N 组分混合。 |
| `--gamma-cats` | int ≥ 1 | 4 | `-dgam N` | 离散 Gamma 速率类别。 |
| `--start-tree` | Path | None | `-t <file>` | 起始树。与 `--fix-tree` 互斥。 |
| `--fix-tree` | Path | None | `-T <file>` | 固定拓扑。与 `--start-tree` 互斥。 |

### 链与并行

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--chains` | int ≥ 1 | 3 | 独立链数。 |
| `--chain-prefix` | str | `chain` | 自动命名链的前缀。 |
| `--chain-names` | str | None | 逗号分隔的名称。覆盖 `--chains`/`--chain-prefix`。 |
| `--threads / -t` | int ≥ 2 | 4 | 每条链的 MPI 进程数（`mpirun -np`）。 |

### 采样

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--sample-freq` | int ≥ 1 | 1 | 每 N 个循环保存一个点。 |
| `--nsamples` | int | `-1` | 每条链在 N 个 MCMC 总循环后停止。`-1` = 无限运行。使用 `--sample-freq N` 时，保存点数 = 循环数 / N。 |

### 收敛监控

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--monitor-freq` | int ≥ 1 | 100 | 每 N 个新样本运行收敛检查。 |
| `--burnin-frac` | float `[0.0, 1.0)` | 0.5 | 监控中丢弃的样本比例。 |
| `--poll-interval` | int ≥ 1 | 60 | .trace 文件轮询间隔（秒）。 |

### 恢复

| 参数 | 类型 | 说明 |
|------|------|------|
| `--resume [CHAINS]` | 可选 str | 从 `run_state.json` 恢复。与 `--overwrite` 互斥。 |

### 工具

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--pb-path` | Path | None | 包含 `pb_mpi`、`bpcomp`、`tracecomp` 和 `mpirun` 的目录。 |
| `--dry-run` | flag | False | 打印命令不执行。 |
| `-q, --quiet` | flag | False | 抑制终端输出。 |

## 工作原理

### 启动

1. 验证输入（矩阵存在、互斥检查、参数范围）。
2. 检测工具（`pb_mpi`, `bpcomp`, `tracecomp`, `mpirun`）。缺失 → 退出码 3。
3. 准备输出目录。写入 `run_state.json`。
4. 如需将 FASTA 自动转换为 PHYLIP。
5. 启动所有链进程（`subprocess.Popen`，工作目录 `chains/`）。
6. 进入监控循环。

### 监控循环

- **每 `--poll-interval` 秒：** 轮询 `.trace` 文件，更新进度显示，检查 `--nsamples` 目标。
- **每 `--monitor-freq` 新样本：** 触发收敛检查（bpcomp + tracecomp）。
- **Ctrl+C：** 软停止所有链，最终收敛检查，写入 `result.json`。

## 收敛阈值

| 指标 | 良好 | 可接受 | 未收敛 |
|------|------|--------|--------|
| `bpcomp maxdiff` | < 0.1 | < 0.3 | ≥ 0.3 |
| `tracecomp min effsize` | > 300 | > 50 | ≤ 50 |
| `tracecomp max rel_diff` | < 0.1 | < 0.3 | ≥ 0.3 |

## 恢复语义

`run_state.json` 是恢复的权威来源。使用 `--chain-names` 添加链时会验证模型参数一致性。`--resume` 读取存储状态，可选择覆盖 `--nsamples`。

## 退出码

| 码 | 含义 |
|----|------|
| 0 | 成功（链完成或被软停止） |
| 1 | 输入验证错误 |
| 2 | pb_mpi 链非零退出 |
| 3 | 必需工具未找到 |

---

# phyloai tree bi bpcomp

使用用户指定的整数 burn-in 运行一次 `bpcomp` 进行最终拓扑收敛分析。

## 用法

```bash
phyloai tree bi bpcomp --chain-dir <链目录> [OPTIONS]
```

## 示例

```bash
phyloai tree bi bpcomp --chain-dir runs/tree/bi/chains --burnin 1000
phyloai tree bi bpcomp --chain-dir runs/tree/bi/chains --burnin 5000 --sample-freq 10
phyloai tree bi bpcomp --chain-dir runs/tree/bi/chains --chain-names chain1,chain2 --burnin 5000
```

## 参数

### 输入 / 输出

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--chain-dir` | Path | **必需** | 包含 `.chain` 文件的目录。 |
| `--chain-names` | str | `all` | 逗号分隔名称。`all` = 自动发现。 |
| `--output-dir / -o` | Path | `runs/tree/bi/bpcomp` | 输出目录。 |
| `--overwrite` | flag | False | 删除并重建输出目录。 |

### 分析

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--burnin` | int ≥ 0 | 0 | 丢弃的保存样本数。 |
| `--sample-freq` | int ≥ 1 | 1 | burn-in 后子采样频率。 |
| `--until` | str | `all` | 停止样本索引。`all` = 整条链。 |
| `--cutoff` | float (0,1) | 0.5 | 多数规则一致树阈值。 |

### 工具

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--pb-path` | Path | None | 包含 `bpcomp` 的目录。 |
| `--dry-run` | flag | False | 打印命令不执行。 |
| `-q, --quiet` | flag | False | 抑制终端输出。 |

## 退出码

| 码 | 含义 |
|----|------|
| 0 | 成功 |
| 1 | 输入验证错误 |
| 2 | bpcomp 非零退出 |
| 3 | `bpcomp` 未找到 |

---

# phyloai tree bi tracecomp

使用用户指定的整数 burn-in 运行一次 `tracecomp` 进行最终参数收敛分析。

## 用法

```bash
phyloai tree bi tracecomp --chain-dir <链目录> [OPTIONS]
```

## 示例

```bash
phyloai tree bi tracecomp --chain-dir runs/tree/bi/chains --burnin 5000
```

## 参数

### 输入 / 输出

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--chain-dir` | Path | **必需** | 包含 `.trace` 文件的目录。 |
| `--chain-names` | str | `all` | 逗号分隔名称。`all` = 自动发现。 |
| `--output-dir / -o` | Path | `runs/tree/bi/tracecomp` | 输出目录。 |
| `--overwrite` | flag | False | 删除并重建输出目录。 |

### 分析

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--burnin` | int ≥ 0 | 0 | 丢弃的保存样本数。 |

### 工具

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--pb-path` | Path | None | 包含 `tracecomp` 的目录。 |
| `--dry-run` | flag | False | 打印命令不执行。 |
| `-q, --quiet` | flag | False | 抑制终端输出。 |

tracecomp 输出被捕获并以逐行 `[good]`/`[ok]`/`[no]` 标注打印。

## 退出码

| 码 | 含义 |
|----|------|
| 0 | 成功 |
| 1 | 输入验证错误 |
| 2 | tracecomp 非零退出 |
| 3 | `tracecomp` 未找到 |

---

# phyloai tree bi readpb

在单条链上运行 `readpb_mpi` 进行后验分析，支持多种分析模式并自动进行格式转换。

## 用法

```bash
phyloai tree bi readpb --chain <链路径> --mode <模式> [OPTIONS]
```

## 示例

```bash
phyloai tree bi readpb --chain chains/chain1 --mode ss,rr --burnin 5000
phyloai tree bi readpb --chain chains/chain1 --mode ss,rr,r --burnin 5000
phyloai tree bi readpb --chain chains/chain1 --mode r --burnin 1000
phyloai tree bi readpb --chain chains/chain1 --mode allppred --burnin 2000
```

## 参数

### 输入 / 输出

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--chain` | Path | **必需** | 无扩展名的链文件路径。 |
| `--mode` | str | **必需** | 逗号分隔的分析模式。 |
| `--output-dir / -o` | Path | `runs/tree/bi/readpb` | readpb 输出和 `result.json` 的输出目录。 |
| `--overwrite` | flag | False | 删除并重建 `--output-dir`。 |

### 分析

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--burnin` | int ≥ 0 | 0 | 丢弃的保存样本数。 |
| `--sample-freq` | int ≥ 1 | 1 | burn-in 后子采样频率。 |
| `--until` | str | `all` | 停止样本索引。 |
| `--threads / -t` | int ≥ 2 | 4 | MPI 进程数。 |

### 工具

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--pb-path` | Path | None | 包含 `readpb_mpi` 和 `mpirun` 的目录。 |
| `--dry-run` | flag | False | 打印命令不执行。 |
| `-q, --quiet` | flag | False | 抑制终端输出。 |

## `--mode` 值

| 值 | readpb_mpi 标志 | 输出 | 说明 |
|---|-----------------|------|------|
| `rr` | `-rr` | `.meanrr` → `.exchangeabilities` | 后验平均交换率（PAML 格式）。 |
| `ss` | `-ss` | `.siteprofiles` → `.sitefreq` | 位点特异频率（IQ-TREE 格式）。 |
| `r` | `-r` | `.meansiterates` | 后验平均位点速率。 |
| `sitelogl` | `-sitelogl` | `.sitelogl`, `.cpo` | 位点边缘对数似然及交叉验证。 |
| `ppred` | `-ppred` | `.ppred` | 后验预测分布的 MSA 模拟。 |
| `div` | `-div` | `.div` | 多样性检验（PPA-DIV）。 |
| `sitecomp` | `-sitecomp` | `.sitecomp` | 组成异质性（PPA-VAR）。 |
| `siteconvprob` | `-siteconvprob` | `.siteconvprob` | 收敛概率（PPA-CONV）。 |
| `comp` | `-comp` | `.comp` | 组成同质性检验。 |
| `allppred` | `-allppred` | `.ppred` | 联合后验预测检验。 |

`allppred` 与 `div`, `sitecomp`, `siteconvprob`, `comp` 互斥。

## 后处理

### `rr` → exchangeabilities

将 `<chain>.meanrr` 转换为 PAML 下三角格式（`<chain>.exchangeabilities`），使用 IQ-TREE 兼容的氨基酸顺序和均匀先验频率。

### `ss` → sitefreq

将 `<chain>.siteprofiles` 转换为 IQ-TREE `-fs` 格式（`<chain>.sitefreq`），从 PhyloBayes AA 顺序重排序到 IQ-TREE 顺序，`1e-8` 下限并重新归一化。

### `ss,rr,r` → PMSF 模拟分区

`r` 输出是无 header、从 0 开始的 `site rate` 后验平均位点速率；链的 trace 根据指定 burn-in/子采样窗口提供后验平均 alpha；链 log 提供离散 Gamma 类别数。它们与每个位点频率 profile 及同步生成的 `<chain>.exchangeabilities` 模型组合，写出含 `+Gk{alpha}` 的 `partition.PMSF.nex`。示例中的 `-p` 表示 edge-proportional 分区；仅在需要 edge-equal 分支长度时改用 `-q`。

```bash
iqtree3 --alisim simulated.phy -t tree.nwk -p runs/tree/bi/readpb/partition.PMSF.nex
```

## 退出码

| 码 | 含义 |
|----|------|
| 0 | 成功 |
| 1 | 输入验证错误 |
| 2 | readpb_mpi 非零退出 |
| 3 | `readpb_mpi` 或 `mpirun` 未找到 |
