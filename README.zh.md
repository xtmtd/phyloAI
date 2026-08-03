# PhyloAI

[English](README.md) | [中文](README.zh.md)

AI-native 模块化系统发育分析平台。

PhyloAI 将常见的 pre-tree、tree 与 post-tree 任务连接为有据可查的命令行工作流，同时通过结构化的 `result.json` 输出与可选的 HTML 报告让每一步都可被检查。它旨在把标记序列文件夹转化为可追溯的系统发育学证据：每个主要步骤都会记录运行内容、使用的工具、保留的位点或分类单元，以及支持最终树的诊断信息。

借助 MCP server 与引导式工作流 Skill，PhyloAI 同样可以对话式地驱动：AI 助手可以检查环境、选择正确的命令、审阅参数、查看运行状态、读取 `result.json`、诊断失败，并帮助解释结果。完成的分析可以汇总为可读的 HTML 报告，其中嵌入图表、表格、出处与方法草稿。

## 为什么选择 PhyloAI

现代系统发育学已不再只是推断一棵最佳树。实际的分析需要数据清洗、比对、标记诊断、基因树与物种树推断、拓扑检验、定年、一致性分析以及可发表的报告。现有工具通常在某一阶段表现出色，但需要用户手工拼接参数、日志、文件格式与诊断。

PhyloAI 专注于这个集成层并使其更强：

- **一套框架应对多种研究设计：** 从同一标记目录运行 supermatrix 或 supertree 分析，无需手工重建工作流即可延伸至拓扑检验、一致性因子、贝叶斯推断、分子定年与报告。
- **透明而非黑盒：** 每个命令都将参数、检测到的工具版本、日志、决策表、输出路径与摘要写到可预测的输出目录。
- **以质量控制为中心：** 标记统计、对称性检验、TAPER 屏蔽、TreeShrink 剪枝、基于聚类的标记探索以及相关性图，在最终推断前让数据问题变得可见。
- **无锁定地集成最佳实践后端：** PhyloAI 编排 MAFFT、trimAl、BMGE、ClipKIT、TAPER、IQ-TREE3、FastTree、wASTRAL、PhyloBayes-MPI、TreeShrink、MCMCtree 等成熟工具，同时保留其原生行为与引用。
- **可恢复且可审计的运行：** 结构化的 `result.json` 文件使中断的分析、失败的位点、保留/丢弃决策、工具版本差异易于检查或重跑。
- **AI-native 设计：** MCP server 与引导式工作流 Skill 让 AI 助手检查 schema、查看运行状态、读取结果、解释分析、支持对话式系统发育学，同时不隐藏命令行细节。
- **可读的报告，不只是文件：** `phyloai report` 将运行目录转为自包含的 HTML 报告，其中嵌入图表、可排序表格、出处与方法草稿，便于稿件准备。

输入边界是刻意的：PhyloAI 不组装 reads、不调用靶标、不推断同源组，也不提取 BUSCO/UCE 风格的标记。这些上游步骤应在使用 PhyloAI 前完成。一旦标记文件就绪，PhyloAI 就提供一个轻量、可脚本化的框架，用于比对、过滤、supermatrix 或 supertree 推断、诊断与报告。

## 工作流概览

PhyloAI 命令遵循系统发育研究的常见流程：

1. **准备标记序列：** 转换格式、检查序列统计、比对未比对位点、修剪 MSA，并移除有问题的位点或分类单元。
2. **评估标记：** 计算占用率、熵、成对一致性、组成偏差、饱和度、树距离与相关性摘要。
3. **构建矩阵或基因树：** 拼接保留位点用于 supermatrix 分析，或为溯祖工作流推断每位点基因树。
4. **推断物种关系：** 运行 ML supermatrix 树、贝叶斯分析或 wASTRAL 物种树推断。
5. **诊断冲突与稳健性：** 计算一致性因子、运行拓扑检验、比较标记聚类，并检查定年分析中的 posterior/prior 行为。
6. **报告运行：** 将参数、工具版本、表格、图表与方法草稿收集为可复现的 JSON 与 HTML 报告。

每个命令写入自己的输出目录。重要的文件通常包括 `result.json`、日志、决策表与最终的序列/树/报告输出。这让失败或部分的运行更易于检查、恢复或借助 AI 助手解释。

## 安装

```bash
git clone https://github.com/xtmtd/phyloAI.git
cd phyloAI
pip install -e .
```

PhyloAI 不捆绑第三方系统发育学可执行文件。安装工作流所需的外部工具，然后用以下命令验证：

```bash
phyloai doctor
```

参见 [docs/commands/installation.md](docs/commands/installation.md) 了解 Python 环境选项、外部工具组与操作系统相关说明。

## 快速开始

```bash
phyloai doctor
```

从原始序列到物种树的一键系统发育流水线：

```bash
phyloai run --seq-dir ./markers
phyloai run --seq-dir ./markers --mode supertree --speed fast --threads 16
```

三步 AliSim 模拟工作流：从 IQ-TREE 报告提取每个位点的实证参数，
批量模拟比对，再重新应用原始 gap 掩码：

```bash
phyloai posttree simulate alisim params --iqtree-dir reports --tree-dir trees -o runs/params
phyloai posttree simulate alisim iqtree --model-params runs/params/params.tsv --strategy pdf --num-simulations 100 -o runs/sim
phyloai posttree simulate alisim transfergaps --original-msa markers/concat.aa.fa --simulated-dir runs/sim/MSAs -o runs/transfer
```

显示所有可用命令：

```bash
phyloai --help
```

## Shell 补全

PhyloAI 可以为 Bash、Zsh、Fish 生成静态的 shell 补全脚本：

```bash
phyloai completion bash
phyloai completion zsh
phyloai completion fish
```

生成脚本一次并配置你的 shell 加载保存的文件。参见 [docs/commands/completion.md](docs/commands/completion.md) 了解 Bash、Zsh、Fish 的配置示例。

## AI 集成

PhyloAI 包含一个 MCP server 与一个引导式工作流 Skill 以支持对话式分析。参见 [docs/commands/ai-integration.md](docs/commands/ai-integration.md) 了解与 OpenCode、Claude Code 或 Codex 的配置与使用。

## 许可证

PhyloAI 作者编写的代码在学术、教育、非商业研究目的下可自由使用、复制、修改与分发。商业用途、商业再分发、再许可、销售或集成到商业产品或服务需要版权所有者的事先书面许可。详见 [LICENSE](LICENSE)。

本仓库同时与第三方软件互操作，每个第三方组件保留自己的许可证。[docs/commands/installation.md](docs/commands/installation.md) 列出的工具是外部依赖，必须在其上游许可证下安装与使用。本节是项目级许可证声明，不能替代那些第三方许可证。

## 作者与联系方式

张峰  
南京农业大学  
Email: <xtmtd.zf@gmail.com>

## 命令

`phyloai doctor` 是唯一支持 `--output-format text|json` 的命令。其他命令将结构化结果写入各自输出目录的 `result.json`，除非设置 `--quiet`，否则使用 Rich 终端输出。

| 命令 | 用途 | 文档 |
|---------|---------|---------------|
| `phyloai doctor` | 检查外部工具可用性。 | [docs/commands/doctor.md](docs/commands/doctor.md) |
| 安装 | 设置 Python 环境与外部工具，然后用 `phyloai doctor` 验证。 | [docs/commands/installation.md](docs/commands/installation.md) |
| `phyloai completion` | 生成静态的 Bash、Zsh 或 Fish shell 补全脚本。 | [docs/commands/completion.md](docs/commands/completion.md) |
| `phyloai run`     | 从原始序列到物种树的一键式系统发育流水线。 | [docs/commands/run.md](docs/commands/run.md) |
| `phyloai pretree convert` | 在下游分析前归一化与转换序列文件。 | [docs/commands/pretree-convert.md](docs/commands/pretree-convert.md) |
| `phyloai pretree stats`   | 检查一个序列/比对文件，或汇总一个目录的文件。 | [docs/commands/pretree-stats.md](docs/commands/pretree-stats.md)     |
| `phyloai pretree align`   | 使用 MAFFT 或 MAGUS 比对序列。 | [docs/commands/pretree-align.md](docs/commands/pretree-align.md)     |
| `phyloai pretree trim`    | 使用 trimAl、BMGE 或 ClipKIT 后端批量修剪已比对 MSA。 | [docs/commands/pretree-trim.md](docs/commands/pretree-trim.md)       |
| `phyloai pretree metrics` | 计算 MSA/树度量、生成分布图与紧凑的相关性热图，用于标记评估。 | [docs/commands/pretree-metrics.md](docs/commands/pretree-metrics.md) |
| `phyloai pretree filter`  | 标记级过滤：TAPER 错误位点屏蔽、TreeShrink 分类单元剪枝、度量规则过滤、对称性检验过滤、基于聚类的探索。 | [docs/commands/pretree-filter.md](docs/commands/pretree-filter.md) |
| `phyloai pretree concat`  | 将多个 MSA 拼接为带占用率过滤、重编码、密码子变体与外类群重排的超矩阵。 | [docs/commands/pretree-concat.md](docs/commands/pretree-concat.md) |
| `phyloai tree ml fasttree` | 使用 FastTree 推断 ML 基因树或超矩阵树。 | [docs/commands/tree-ml-fasttree.md](docs/commands/tree-ml-fasttree.md) |
| `phyloai tree ml iqtree`   | 使用 IQ-TREE3 推断 ML 树：同质、异质、分区、ModelFinder，以及自定义交换率/位点频率 profile 工作流。 | [docs/commands/tree-ml-iqtree.md](docs/commands/tree-ml-iqtree.md) |
| `phyloai tree bi pb`    | 使用 PhyloBayes-MPI 进行 MCMC 链推断：多链并行、实时收敛监控、轨迹图与 resume。 | [docs/commands/tree-bi.md](docs/commands/tree-bi.md) |
| `phyloai tree bi bpcomp` | 使用 bpcomp 进行最终拓扑收敛分析（用户指定 burn-in）。 | [docs/commands/tree-bi.md](docs/commands/tree-bi.md) |
| `phyloai tree bi tracecomp` | 使用 tracecomp 进行最终参数收敛分析（用户指定 burn-in）。 | [docs/commands/tree-bi.md](docs/commands/tree-bi.md) |
| `phyloai tree bi readpb` | 使用 readpb_mpi 进行后验分析与预测检验；`--mode ss,rr,r` 还会生成使用后验参数的 PMSF 模拟分区。 | [docs/commands/tree-bi.md](docs/commands/tree-bi.md) |
| `phyloai tree msc`   | 使用 wASTRAL 进行多物种溯祖物种树推断。 | [docs/commands/tree-msc.md](docs/commands/tree-msc.md) |
| `phyloai tree cf`    | 一致性因子计算：gCF、sCF、sCFl（IQ-TREE3）和 qCF（wASTRAL）。 | [docs/commands/tree-cf.md](docs/commands/tree-cf.md) |
| `phyloai posttree topology` | 树拓扑检验（AU / KH / SH / WKH / WSH / c-ELW），将候选树与超矩阵进行比较。 | [docs/commands/posttree-topology.md](docs/commands/posttree-topology.md) |
| `phyloai posttree dating`  | 使用 MCMCtree 进行贝叶斯分子定年：IQ-TREE Hessian 计算 + 带诊断的 MCMC 分歧时间估计。 | [docs/commands/posttree-dating.md](docs/commands/posttree-dating.md) |
| `phyloai posttree signal` | 系统发育信号分布分析：位点/基因 lnL 打分、一致基因识别、四簇似然映射。命令：`signal lnl`、`signal consistent`、`signal fclm`。 | [docs/commands/posttree-signal.md](docs/commands/posttree-signal.md) |
| `phyloai posttree modelcompare iqtree` | 使用 IQ-TREE3 ModelFinder 进行相对模型比较（BIC/AIC/AICc），支持通过 `-madd` 展开异质混合模型。 | [docs/commands/posttree-modelcompare.md](docs/commands/posttree-modelcompare.md) |
| `phyloai posttree modelcompare pb` | 使用 PhyloBayes `.sitelogl` 位点对数似然文件进行 LOO-CV / wAIC 相对模型比较（Lartillot 2023），纯 Python 实现。 | [docs/commands/posttree-modelcompare.md](docs/commands/posttree-modelcompare.md) |
| `phyloai posttree simulate alisim` | 基于 IQ-TREE3 AliSim 的序列模拟，保留数据集实证特征：`params` 从 IQ-TREE 报告中提取逐位点参数，`iqtree` 模拟单个或批量 MSA（complete/mixed/pdf 策略，可恢复），`transfergaps` 将原始 gap 掩码重新引入一条或多条模拟比对。 | [docs/commands/posttree-simulate-alisim.md](docs/commands/posttree-simulate-alisim.md) |
| `phyloai report`   | 生成可复现的分析报告（JSON + 自包含的 HTML，含嵌入图表、可排序表格与方法段落草稿）。自动生成的方法文本在发表前应仔细核对。 | [docs/commands/report.md](docs/commands/report.md) |
