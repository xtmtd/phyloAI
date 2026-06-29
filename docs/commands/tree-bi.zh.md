# phyloai tree bi

[English](tree-bi.md) | [中文](tree-bi.zh.md)

使用 [PhyloBayes-MPI](https://github.com/bayesiancook/pbmpi)（`pb_mpi`）进行贝叶斯系统发育推断。

## 目的

`phyloai tree bi` 并行运行 N 条独立的 MCMC 链（采用 PhyloBayes 替换模型），使用 `bpcomp` 和 `tracecomp` 实时监控链收敛，并在链停止后生成一致树。

与 `tree ml` 和 `tree msc` 不同，`tree bi` 是一个长时间运行的交互式命令。链可能运行数小时或数天。命令在整个过程中保持活动，显示实时进度展示和周期性收敛统计，直到用户使用 Ctrl+C（软停止）终止链或达到链目标。

`tree bi` 没有子命令层 —— 它以 `phyloai tree bi [OPTIONS]` 形式调用。

默认输出目录为 `runs/tree/bi/`。

## 要求

PhyloBayes-MPI 工具必须安装且可在 `PATH` 上发现（或通过 `--pb-path`）：

| Tool | Purpose |
|------|---------|
| `pb_mpi` | MCMC 采样器（始终必需） |
| `bpcomp` | 拓扑收敛（始终必需） |
| `tracecomp` | 参数收敛（始终必需） |
| `mpirun` | Open MPI 启动器（始终必需） |
| `readpb_mpi` | 读取链文件（可选） |

运行 `phyloai doctor` 确认安装。

## 用法

```bash
phyloai tree bi --matrix <alignment> [OPTIONS]
```

输入比对可以是 PHYLIP 或 FASTA。pb_mpi 要求 PHYLIP；FASTA 会在链启动前自动转换。

## 示例

```bash
# 默认：3 条链，CAT-GTR 模型，永久运行
phyloai tree bi --matrix concat/matrix.phy

# 同质 LG+G4，保存 10000 个采样点后停止
phyloai tree bi --matrix concat/matrix.phy --model lg --mixture 1 --nsamples 10000

# WAG+C20 混合模型
phyloai tree bi --matrix concat/matrix.phy --model wag --mixture 20

# 向已有运行添加两条额外链
phyloai tree bi --matrix concat/matrix.phy --chain-names chain4,chain5 -o runs/tree/bi

# 从之前状态恢复所有链
phyloai tree bi -o runs/tree/bi --resume

# 仅恢复 chain1 与 chain3
phyloai tree bi -o runs/tree/bi --resume chain1,chain3

# 恢复并扩展到新的目标
phyloai tree bi -o runs/tree/bi --resume --nsamples 10000

# 恢复并永久运行（覆盖之前目标）
phyloai tree bi -o runs/tree/bi --resume --nsamples -1

# 自定义 PhyloBayes 工具目录
phyloai tree bi --matrix concat/matrix.phy --pb-path /opt/pbmpi/bin

# 打印将运行的命令然后退出
phyloai tree bi --matrix concat/matrix.phy --dry-run
```

## 参数

### Input / Output

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--matrix / -m` | Path | required | 输入比对（PHYLIP 或 FASTA）。FASTA 在链启动前转换为 PHYLIP。使用 `--resume` 时不需要。 |
| `--output-dir / -o` | Path | `runs/tree/bi` | 输出目录。包含 `chains/`、`convergence/`、`run_state.json` 与 `result.json`。 |
| `--overwrite` | flag | False | 启动前删除并重建输出目录。与 `--resume` 互斥。 |

### Model

| Flag | Choice | Default | pb_mpi flag | Description |
|------|--------|---------|-------------|-------------|
| `--model` | `gtr`, `poisson`, `lg`, `wag`, `jtt`, `mtrev`, `mtzoa`, `mtart` | `gtr` | `-gtr`, `-poisson`, … | 速率矩阵。 |
| `--mixture` | str | `auto` | `-cat` / `-ncat N` | `auto` = CAT Dirichlet 过程；`1` = 同质（如 LG+G4）；整数 N > 1 = 固定 N 分量混合。 |
| `--gamma-cats` | int ≥ 1 | 4 | `-dgam N` | 离散 Gamma 速率类别数。 |
| `--start-tree` | Path | None | `-t <file>` | 起始树（Newick）。拓扑可在 MCMC 中改变。与 `--fix-tree` 互斥。 |
| `--fix-tree` | Path | None | `-T <file>` | 固定拓扑（Newick）。仅采样支长与其他参数。必须二叉。与 `--start-tree` 互斥。 |

**简写：**

| PhyloAI 调用 | pb_mpi 等价 | IQ-TREE 类比 |
|---|---|---|
| (defaults) | `-cat -gtr -dgam 4` | CAT-GTR |
| `--model lg --mixture 1` | `-lg -ncat 1 -dgam 4` | LG+G4 |
| `--model poisson` | `-cat -poisson -dgam 4` | CAT-Poisson |
| `--model wag --mixture 20` | `-wag -ncat 20 -dgam 4` | WAG+C20 |

### Chains & Parallelism

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--chains` | int ≥ 1 | 3 | 独立链数；自动命名为 `<prefix>1`、`<prefix>2`、… |
| `--chain-prefix` | str | `chain` | 自动命名链的前缀。 |
| `--chain-names` | str | None | 逗号分隔的名字（如 `chain4,chain5`）。覆盖 `--chains` 与 `--chain-prefix`。用于向已有运行添加链。 |
| `--threads / -t` | int ≥ 2 | 4 | 每条链的 MPI 进程数（`mpirun -np`）。最少 2（1 master + N-1 slaves）。 |

生效的名字列表为：若提供 `--chain-names`，使用它；否则生成 `[prefix+str(i) for i in range(1, chains+1)]`。

### Sampling

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--sample-freq` | int ≥ 1 | 1 | 每 N 个 cycle 保存一个 MCMC 点（pb_mpi `-x <every>`）。 |
| `--nsamples` | int | `-1` | 每条链运行 N 个 MCMC cycle 后停止（pb_mpi `-x <until>`）。`-1` = 永久运行；用 Ctrl+C 停止。保存点数 = `nsamples / sample-freq`。 |

要停止永久运行的链：使用 Ctrl+C，或 `echo 0 > chains/<chainname>.run`。

### Convergence Monitoring

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--monitor-freq` | int ≥ 1 | 100 | 每 N 个新采样点运行一次 `bpcomp` + `tracecomp`（按最小链长增量计）。 |
| `--burnin-frac` | float `[0.0, 1.0)` | 0.5 | 丢弃的保存样本比例作为 burn-in，**仅在收敛监控期间**。动态应用：`burnin = floor(min_chain_length × burnin_frac)`。最小 burn-in 为 10 个样本；若链太短则跳过检查并警告。不传递给 pb_mpi。 |
| `--poll-interval` | int ≥ 1 | 60 | `.trace` 文件轮询间隔（秒），用于进度展示与收敛触发。 |

### Resume

| Flag | Type | Description |
|------|------|-------------|
| `--resume [CHAINS]` | optional str | 从 `run_state.json` 恢复。无值 = 恢复所有链。逗号分隔的名字 = 仅恢复这些链（如 `--resume chain1,chain3`）。与 `--overwrite` 互斥。 |

Resume 使用原生 pb_mpi 机制：`mpirun -np <threads> pb_mpi <chainname>`（无 `-d` 或模型标志；pb_mpi 读取已有 `.chain` 文件）。当 `--resume` 与 `--nsamples` 同时使用时，新值覆盖 `run_state.json` 中存储的目标，链从当前状态继续。已达到解析目标的链会被跳过。运行中的链在达到目标时收到软停止。

Click 实现：`@click.option('--resume', default=None, is_flag=False, flag_value='__ALL__', help='...')` —— 不出现 = `None`，裸 `--resume` = `'__ALL__'`，`--resume chain1,chain2` = `'chain1,chain2'`。

### Tool

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--pb-path` | Path | None | 包含 `pb_mpi`、`bpcomp`、`tracecomp` 的目录。覆盖 PATH 查找。若此处存在 `readpb_mpi` 也会被检测，但不是必需的。 |
| `--dry-run` | flag | False | 打印所有命令而不执行。 |
| `-q, --quiet` | flag | False | 抑制终端输出（进度展示与收敛统计）。 |

## 工作原理

### 启动

1. 校验输入（矩阵存在、`--start-tree`/`--fix-tree` 互斥、`--mixture` 值、`--threads ≥ 2`）。
2. 通过 `ToolEnv` 检测工具（`pb_mpi`、`bpcomp`、`tracecomp`、`mpirun`）。缺失必需工具 → 退出码 3。
3. 准备输出目录：使用 `--overwrite` 时删除并重建；非空目录在无 `--overwrite` 时拒绝，除非 `--resume` 或 `--chain-names` 添加新链。
4. 创建 `chains/` 与 `convergence/` 子目录。
5. 写入 `run_state.json`（仅在全新运行；新增链时更新而非替换）。
6. 构建并启动所有链进程（`subprocess.Popen`），工作目录 `chains/`。
7. 进入监控循环。

### 链命令

全新运行：

```
mpirun -np <threads> pb_mpi -d <abs_matrix_path> [model_flags] -x <sample_freq> <nsamples> <chainname>
```

恢复：

```
mpirun -np <threads> pb_mpi <chainname>
```

无 `-d` 或模型标志；pb_mpi 从已有 `.chain` 文件读取所有设置。

### 监控循环

单个主线程轮询 `.trace` 文件并在链子进程运行时触发收敛检查：

- **每 `--poll-interval` 秒**：轮询每条链的 `.trace` 文件，更新进度展示，并检查是否达到 `--nsamples`（对已达到的链向 `chains/<chain>.run` 写入 `0`）。
- **每 `--monitor-freq` 个新采样点**：触发 `bpcomp` + `tracecomp` 检查（"收敛监控"一节）。收敛输出通过 `live_display.stop() → print() → live_display.start()` 渲染，因此进度条得以保留。
- **Ctrl+C 时**：向所有 `chains/<chainname>.run` 文件写入 `0`，等待子进程完成当前 cycle，然后运行最终收敛检查并写入 `result.json`。
- **子进程非零退出时**：停止剩余链并写入 `status: "error"` 的 `result.json`。

进度条任务预先初始化为各链的现有 trace 采样数，因此恢复运行立即显示正确值（而非 0）。

## 收敛监控

### 触发与 Burn-in

当 `min(chain_lengths) - last_check_min ≥ monitor-freq` 时触发检查。动态 burn-in 为 `floor(min_chain_length × burnin_frac)`。若 `burnin < 10`，跳过检查并发出警告 `"Skipping convergence check: chains too short (burnin < 10)"`。在 `--resume` 时，新的 `--nsamples`（若提供）成为链目标；已经达到或超过目标的链被跳过。

### bpcomp / tracecomp 调用

对于 N 条链，PhyloAI 运行：

- **所有链**：`bpcomp -x <burnin> -o bpcomp_all ../chains/chain1 ../chains/chain2 …`（cwd `convergence/`）
- **所有两两组合**：`bpcomp -x <burnin> -o bpcomp_chain1_chain2 ../chains/chain1 ../chains/chain2`，等等。
- **tracecomp** 结构相同，输入 `.trace` 文件，输出 `.contdiff`。

两个命令的工作目录都是 `convergence/`，因此链文件被引用为 `../chains/<chain>`，trace 文件为 `../chains/<chain>.trace`。

### 阈值

| Metric | Good | Acceptable | Not converged |
|--------|------|------------|---------------|
| `bpcomp maxdiff` | < 0.1 | < 0.3 | ≥ 0.3 |
| `tracecomp min effsize` | > 300 | > 50 | ≤ 50 |
| `tracecomp max rel_diff` | < 0.1 | < 0.3 | ≥ 0.3 |

每列状态：`bpcomp` 仅使用 `maxdiff`；`tracecomp` 使用 `min effsize` 与 `max rel_diff` 中较差者。整体状态 = 所有三个度量中最差者。该命令**不**自动停止链。

### 终端展示

```
  All chains
  bpcomp    maxdiff  0.081   meandiff  0.006   [good]
  tracecomp  min effsize  312   max rel_diff  0.094   [good]

  Pairwise
    pair              maxdiff  min effsize  max rel_diff  bpcomp  tracecomp
    chain1 x chain2   0.073       340           0.094     good       good
    chain1 x chain3   0.432        76           0.210       no         ok
    chain2 x chain3   0.065       355           0.072     good       good
```

表格下方的分层通知：

- 所有对 `good` → `"*** All convergence criteria met (all pairs good). You may stop chains with Ctrl+C when ready. ***"`
- 所有对至少 `ok` → `"Convergence acceptable across all chain pairs (N good, M ok). Consider stopping when ready."`
- 一些对已收敛 → `"Some chain pairs agree (N good, M ok, K not converged)."`
- 都未收敛 → 无通知（仅表格本身足够）。

### Trace 图

在每次收敛检查时，PhyloAI 使用 matplotlib 重新生成 `convergence/trace_plots.pdf`。每个 trace 参数列一页（除 `iter` 与 `time` 之外的所有列）；每条链一条线。垂直虚线标记当前 burn-in 位置。若 matplotlib 未安装，则静默跳过 PDF 生成并一次性打印 `"matplotlib not available; trace plots disabled."`。

## Resume 语义

`run_state.json` 是 resume 的真相之源。它在全新启动时创建，在通过 `--chain-names` 添加新链时更新（而非替换）。在 `--resume` 与新 `--nsamples` 一起使用时再次更新。该文件在全新启动前永远不存在。

Schema：

```json
{
  "chain_names": ["chain1", "chain2", "chain3"],
  "matrix": "/abs/path/matrix.phy",
  "model_flags": ["-cat", "-gtr", "-dgam", "4"],
  "sample_freq": 1,
  "nsamples": 10000,
  "threads": 4
}
```

**添加链（对已有目录使用 `--chain-names chain4,chain5`）：** PhyloAI 校验新调用的模型参数（`model_flags`、`sample_freq`、`nsamples`、`threads`）与存储值匹配。若任意不同，退出码 1 并提示 `"Model parameters conflict with existing run_state.json. Use --resume to continue existing chains or choose a different --output-dir."`。它还校验新名字都不已存在。

**恢复：**

1. 读取 `run_state.json` 获取存储的 `nsamples` 与 `chain_names`。
2. 若用户提供 `--nsamples` 且与存储值不同，使用新值并更新 `run_state.json`（这允许扩展已完成的运行）。
3. 解析要恢复的链：`'__ALL__'` → 所有名字；逗号分隔 → 仅这些名字。
4. 对每条要恢复的链：读取其 `.trace`。若 `nsamples != -1` 且 `current_length ≥ nsamples`，跳过该链。
5. 以恢复命令启动剩余链。监控循环在各链达到目标时发出软停止。

之前的 `result.json` 在写入新结果前会自动备份并加时间戳（如 `result_20260624_134500.json`）。多次 resume 各自产生独立备份。

## 安全停止

使用 **Ctrl+C**。PhyloAI 向每条 `chains/<chain>.run` 文件写入 `0`，并等待 pb_mpi 完成当前 cycle。直接中断 pb_mpi（如 `kill -9`）可能使 trace 文件中残留不完整样本。若链尚未达到目标，`result.json` 中的 `data.interrupted` 设为 `true`。

## 输出

```
runs/tree/bi/
├── chains/
│   ├── chain1.trace        # MCMC 轨迹（TSV：iter time topo loglik ...）
│   ├── chain1.treelist     # 采样树（Newick，每行一棵）
│   ├── chain1.chain        # 完整参数状态（二进制；供 readpb_mpi 使用）
│   ├── chain1.param        # 当前参数快照（文本）
│   ├── chain1.monitor      # 混合统计（文本）
│   ├── chain1.run          # 运行标志：1=运行中，0=停止
│   ├── chain1.log          # 合并 stdout+stderr（PhyloAI 写入）
│   ├── chain2.{trace,...}
│   └── chain3.{trace,...}
│
├── convergence/
│   ├── trace_plots.pdf              # 所有 trace 参数，每列一页
│   ├── bpcomp_all.bpdiff            # 全链 bpcomp 汇总（已解析）
│   ├── bpcomp_all.bplist            # 全链 bipartition 列表
│   ├── bpcomp_all.con.tre           # 全链一致树
│   ├── bpcomp_chain1_chain2.bpdiff  # 两两 bpcomp 输出
│   ├── bpcomp_chain1_chain2.con.tre
│   ├── ...
│   ├── tracecomp_all.contdiff             # 全链 tracecomp 输出
│   ├── tracecomp_chain1_chain2.contdiff
│   └── ...
│
├── run_state.json           # Resume 元数据
└── result.json              # PhyloAI 结构化结果
```

`result.json` 遵循标准 schema（status / command / wall_time / tool_versions / params / key_results / error / data）。BI 命令使用**多链扩展**的单模式：

- `data.chain_cmds` 是一个以链名为键的字典；每个值是实际执行的 argv 列表（全新或恢复格式）。
- `data.tool_stderr` 是一个以链名为键的字典；每个值是 `chains/<chain>.log` 中合并的 stdout+stderr 内容。
- `data.tool_logs` 引用每条链的日志文件路径。
- 若运行在达到 `--nsamples` 前通过 Ctrl+C 结束，`data.interrupted` 为 `true`。
- 正常退出与 Ctrl+C 软停止时 `status` 都是 `"success"`；仅当 pb_mpi 链以非零返回码退出时为 `"error"`。
- `key_results.consensus_tree` 仅当 `bpcomp` 运行实际产出该文件时（即收敛检查运行并成功）才设置为 `convergence/bpcomp_all.con.tre`。在检查被跳过、失败或链太短时为 `null`。
- `tool_versions.pb_mpi`、`bpcomp`、`tracecomp` 在有文件名（如 `pb_mpiManual1.9.pdf`、`VERSION`）时从文件名检测，否则为 `null`。`mpirun` 从 `--version` 检测。

## 退出码

| Code | Meaning |
|------|---------|
| 0 | 成功（链已完成或软停止） |
| 1 | 输入校验错误（矩阵缺失、链名冲突、参数组合无效、输出目录冲突） |
| 2 | 工具执行失败（一条 pb_mpi 链以非零返回码退出） |
| 3 | 必需工具未找到（pb_mpi、bpcomp、tracecomp 或 mpirun） |

校验、工具检测与输入冲突失败始终会写入 `status: "error"` 与错误消息的 `result.json`。

## 提示与警告

- **让链预热。** 仅当 `burnin = floor(min_chain_length × burnin_frac) ≥ 10` 时才会触发收敛检查。短链会跳过检查并警告。
- **根据模型调整 `--burnin-frac`。** 较大的值对慢混合模型更安全，但会延迟首次检查。
- **使用 `--threads ≥ 2`。** pb_mpi 要求 1 master + N-1 slaves；CLI 拒绝 `--threads 1`。
- **添加链而无需重跑整组。** `--chain-names chain4,chain5`（使用相同的输出目录）将新链追加到现有的 `run_state.json` 并以全新链方式启动它们。模型参数必须与存储值匹配。
- **扩展已完成的运行。** 使用 `--resume --nsamples <new_target>` 从链停止处继续。已经达到新目标的链会被静默跳过。
- **可视化混合。** `convergence/trace_plots.pdf` 在每次检查时重新生成。在任何 PDF 阅读器中打开以检查链混合情况。
- **慢速运行。** 较大的 `--poll-interval` 在网络文件系统上减少 I/O；较大的 `--monitor-freq` 减少 `bpcomp`/`tracecomp` 的开销。
- **建立 screen / tmux 会话。** 长时间运行的链应从会话管理器启动，以便 Ctrl+C 和 resume 正常工作。