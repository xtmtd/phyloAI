# Parameter Annotations

Runtime schema from `get_command_schema` is authoritative for names, types, defaults, and choices.
These notes add Chinese scientific context and recommended values.
Parameters that are purely operational (`--output-dir`, `--threads`, `--overwrite`, `--resume`,
`--dry-run`, `--quiet`, `-h/--help`) or visual rendering (`--dpi`, `--fig-width`, `--color`, etc.)
are self-explanatory and omitted below.

## Common Parameters

### --msa-dir
输入多序列比对（MSA）目录。包含 `.fa`/`.fasta`/`.faa`/`.fna`/`.aln` 格式文件。
每一轮分析使用一个独立的 `msa-dir`，方便追踪每步进出。

### --tree-dir
输入基因树目录。包含 `.tre`/`.tree`/`.nwk`/`.newick`/`.treefile`/`.bestTree`/`.contree` 格式文件。
通过逻辑基因名（去掉后缀）与 MSA 自动配对。

### --seq-type
分子类型。`AA` 蛋白质，`NT` 核苷酸，`CODON` 编码序列。默认为 `auto`，自动检测但不保证正确；
建议对编码序列显式指定 `CODON`。

### --tool-args
传递给底层工具（MAFFT/IQ-TREE/trimAl/TAPER 等）的策略参数。
**禁止放输入/输出路径**（这些由 PhyloAI 显式管理，会触发 blocked flag 错误）。
典型用法：`--tool-args "--maxsubsetsize 50"` 或 `--tool-args "-g 0.8"`。

### --input-format
输入文件格式提示。`csv`/`tsv` 用于表格文件；`fasta`/`phylip-relaxed`/`nexus` 用于比对文件。
默认为 `auto`（根据内容和文件扩展名检测）；不确定时显式指定。

### --table-format
辅助表格输出格式。默认 `csv`，可选 `tsv`。

### --seq-dir
输入未比对序列目录，用于 `pretree align`（原始序列）和 `pretree stats`。
文件匹配后缀：`.fa`/`.fasta`/`.fas`/`.fna`/`.faa`/`.fq`/`.fastq`。

---

## Doctor

### --output-format
输出格式。`text` 为人眼可读终端表格（默认），`json` 为机器可读 JSON，供 MCP 使用。

---

## Run (One-Click Pipeline)

### --seq-dir
输入未比对序列目录，整个管线的起点。

### --mode
管线模式。`supermatrix`：比对→修剪→拼接→ML 树；`supertree`：比对→修剪→基因树→wASTRAL 物种树。

### --speed
速度策略。`normal`（默认）：MAFFT linsi + TAPER + IQ-TREE；`fast`：MAFFT auto + 跳过 TAPER + FastTree。

---

## Pretree

### pretree convert

#### --input-path, --output-dir, --seq-type
见 Common Parameters。

#### --target-format
目标格式。支持 `fasta`/`phylip`/`nexus`/`phylip-paml`。
PHYLIP-PAML 为 PAML 系列工具的专用双空格 PHYLIP 格式。

#### --aa-special
为蛋白质序列启用特殊字符检查（B/Z/J/X 等非标准氨基酸），对 auto 类型也生效。

---

### pretree stats

#### --seq-dir, --seq-type, --table-format
见 Common Parameters。

#### --seq
同时输出每个序列的长度等逐条统计；不加则仅输出每个文件级别的汇总统计。

#### --unaligned
输入为未比对序列。在此模式下跳过 MSA 特有指标（间隙率等），仅计算序列级统计。

#### --per-gene
为每个基因单独生成一个统计文件（`<gene>.stats.csv`），存放在 `stats/` 子目录中。

---

### pretree align

#### --seq-dir, --seq-type, --tool-args
见 Common Parameters。

#### --method
比对算法/工具。
- `linsi`：MAFFT 迭代精炼，高精度，适合 <500 序列 [推荐：AA 小数据集]
- `einsi`：MAFFT 迭代精炼 + E-INS-i 算法，适合含大内部间隙的序列
- `ginsi`：MAFFT 全局比对，序列长度差异大时避免使用
- `auto`：MAFFT 自动选择策略，适合大型数据集 [推荐：NT 大数据集]
- `fftns1`/`fftns2`：MAFFT 快速渐进法，速度最快但精度最低
- `magus`：MAGUS 引导树分解比对，精度最高但最慢，仅 Linux 可用 [推荐：难比对数据集]

#### --backtrans
启用 trimAl backtranslation：将蛋白质比对结果反翻译为密码子级别的核苷酸比对。
需要 `--nt-dir` 提供原始 CDS 序列。

#### --nt-dir
`--backtrans` 所需的未比对 CDS 核苷酸序列目录。

#### --mafft-path, --magus-path, --trimal-path
显式指定外部工具（MAFFT/MAGUS/trimAl）可执行文件路径。默认通过 PATH 查找。

---

### pretree trim

#### --msa-dir, --seq-type, --tool-args
见 Common Parameters。

#### --tool
修剪工具。
- `trimal`：trimAl，参数灵活，支持 gap/consistency 策略 [推荐：通用]
- `bmge`：BMGE，使用熵和 BLOSUM 矩阵选择区块，偏保守 [推荐：信息位点保留]
- `clipkit`：ClipKIT，四种内置修剪模式，零参数 [推荐：快速标准修剪]

#### --nt-dir
用于 backtrans 模式（AA 比对时保留对应的 NT 比对）。配合 `pretree align --backtrans` 使用。

#### --trimal-method
trimAl 模式。`automated1`（默认，自动选择）、`gappyout`（按 gap 比例裁剪）、`strict`（严格保守）、`strictplus`（加强版 strict）。

#### --bmge-matrix
BMGE 替换矩阵。可选 `BLOSUM30`/`BLOSUM45`/`BLOSUM62`/`BLOSUM90`/`PAM120`/`PAM250` 等。
数值越低越允许远缘序列；BLOSUM90 适合近缘物种。

#### --bmge-entropy
BMGE 熵阈值（0-1）。低于此值的区块被移除，阈值越高越保守。

#### --clipkit-method
ClipKIT 修剪模式。`smart-gap`（默认，间隙率 + 简约信息位点）、`gappy`（仅间隙率）、
`kpic`（简约信息位点）、`kpic-smart-gap`、`kpic-gappy`、`kpi`（信息位点）、`crop-only`。

#### --trimal-path, --bmge-path, --clipkit-path
显式指定工具可执行文件路径。

---

### pretree metrics

#### --msa-dir, --tree-dir, --seq-type, --table-format
见 Common Parameters。

#### --outgroup-list
包含外群分类单元名的文件，每行一个。用于计算 DVMC（方向性向量组分变异），需要 `--tree-dir`。

#### --ref-tree
参考物种树（NEWICK），用于计算每个基因树的标准化 RF 距离。需要 `--tree-dir`。

#### --skip-freq-statistics
跳过逐字符频率列（`freqA`/`freqC`/`freqG`/`freqT` 及 `freq*` 簇）。数据量大时启用可显著减少 CSV 体积和计算时间。

#### --pseudo-tree-metrics
启用 FastTree 推导伪树（约化精度）计算树相关指标。仅需要 `--msa-dir`。指标名以 `_FT` 后缀标识。

#### --fasttree-path
显式 FastTree 可执行文件路径，用于 `--pseudo-tree-metrics`。

#### --skip-pairwise-identity
跳过 `average_pairwise_identity`（O(n² × L) 复杂度）。当标记含许多分类单元时建议启用。

#### --decimal-places
CSV 中小数位数，0–12。默认 6。降低对小数据集可减少文件大小。

---

### pretree metrics plot

#### --csv-path
现有 `metrics.csv` 路径。从 `pretree metrics` 输出中获取。

#### --metric
要绘制的指标列名（必须存在于 CSV 中）。例如 `entropy`、`treeness`、`evo_rate`。

#### --tukey-k
Tukey's Fences 异常值检测乘数。设置后将过滤后的基因列表写入 `<metric>.tukey_filtered.csv`。

---

### pretree metrics correlate

#### --csv-path
现有 `metrics.csv` 路径。

#### --metrics
要计算相关性的指标列名，逗号分隔。设为 `all` 包含全部数值列。省略则自动选择核心可读指标。

#### --include-freq
在自动选择中纳入 `freq*` 列。

#### --include-sd
在自动选择中纳入 `sd_*`（标准差）列。

#### --method
相关性方法。`spearman`（秩相关，默认，稳健）或 `pearson`（线性相关，要求正态性）。

#### --triangle
矩阵显示方式。`full`（完整）、`lower`（下三角）、`upper`（上三角，用于紧凑版）。

#### --cluster-rectangles
在完整矩阵 Ward 聚类图上绘制 N 个分类矩形。仅 `--triangle full` 有效；`lower`/`upper` 模式会被警告忽略。

---

### pretree filter taper

#### --msa-dir, --seq-type, --tool-args, --table-format
见 Common Parameters。

#### --nt-dir
与 AA MSA 对应的核苷酸比对目录。提供后将同时生成蛋白质和密码子水平的错误位点掩码。

#### --cutoff
TAPER 错误位点评分阈值（0-1）。高于此值的位点被标记为错误。默认值随 TAPER 版本而定；越低越激进（屏蔽更多位点）。

#### --taper-path
显式 TAPER（Julia）可执行文件路径。

#### --julia-path
显式 Julia 运行时路径。

#### --show-masked-sites
在输出 FASTA 中用 `X`（AA）或 `N`（NT）标记被屏蔽位点，而非删除。

---

### pretree filter treeshrink

#### --tree-dir, --msa-dir, --tool-args, --table-format
见 Common Parameters。MSA 为可选项：提供后会将剪枝同步应用于其 MSA。

#### --threshold
TreeShrink 异常分支长度阈值（α）。默认 0.05。值越小越激进（更多分类单元被移除）。

#### --treeshrink-mode
TreeShrink 运行模式。`per-locus`（单基因模式）或 `all-gene`（利用多基因联合信息）。

#### --treeshrink-path
显式 TreeShrink（R）可执行文件路径。

#### --keep-work-dir
保留 TreeShrink 中间文件（默认清理）。调试时有用。

---

### pretree filter metrics

#### --table-path
`metrics.csv` 路径（来自 `pretree metrics`）。

#### --keep
保留条件表达式。支持 `column operator value`（如 `entropy > 0.5`）及 `and`/`or` 组合。
例：`num_taxa >= 5 and proportion_gaps < 0.3`。

#### --input-format, --table-format
见 Common Parameters。

#### --loci-column
基因名列名（默认 `loci`）。仅当 CSV 列名不同时需要设置。

#### --msa-dir, --tree-dir
可选。提供后将同步复制筛选通过的 MSA/Tree 文件到输出 `seqs/`/`trees/` 子目录。

#### --copy
输出时可同时复制关联的 MSA/树文件，需配合 `--msa-dir`/`--tree-dir`。

---

### pretree filter symtest

#### --msa-dir, --tree-dir, --iqtree-path, --table-format
见 Common Parameters。`--tree-dir` 为可选，IQ-TREE 也可以从 MSA 推导对称测试树。

#### --symtest-type
对称性检验类型。`overall`（全局检验，默认，速度快）；`pairwise`（逐对检验，精度高但 O(n²)）。

#### --symtest-pval
显著性阈值（p-value）。默认 0.05。p-value 低于此值的基因被标记为组成异质性显著（建议过滤）。

#### --symtest-keep-zero
保留 p-value 恰好为 0 的基因（默认移除）。DJ 检验常在非常小的比对中产生零值。

---

### pretree filter cluster

#### --table-path, --input-format, --table-format
见 Common Parameters。

#### --metrics
用于聚类的指标列名（逗号分隔）。默认使用核心数值指标。

#### --exclude-regex
排除匹配正则表达式的基因名。例如 `^mitochondrial_.*` 排除线粒体基因。

#### --reduction
降维方法。`umap`（默认，保留局部/全局结构）或 `pca`（主成分分析）。

#### --n-clusters
固定聚类数目。省略则自动确定最优 k。

#### --max-clusters
允许的最大聚类数（自动选择时）。默认 20。

#### --cluster-linkage
聚类连接方法。`ward`（默认，最小方差）、`average`、`complete`、`single`。

#### --cluster-distance
聚类距离度量。`euclidean`（默认）、`correlation` 等。

#### --drop-outlier-clusters
自动移除异常聚类（离群基因群）。配合 `--max-drop-fraction` 控制上限。

#### --outlier-metric
用于异常值筛选的指标名。

#### --outlier-direction
异常值方向。`high`（大于阈值）、`low`（小于阈值）、`both`。

#### --max-drop-fraction
最多删除的基因比例（0-1，默认 0.3）。防止过度筛选。

#### --plot-metrics-cols
替代图中使用的指标列（默认使用 `--metrics` 指定的列）。

#### --plot-label-angle
图中标签旋转角度。

#### --outlier-boxplot-cols
箱线图中显示的列。

#### --umap-n-neighbors, --umap-min-dist, --umap-replicates, --umap-random-state
UMAP 参数。`n_neighbors`（局部邻域大小，默认 15）、`min_dist`（最小点间距，默认 0.1）、
`replicates`（重复运行次数）、`random_state`（随机种子）。

#### --msa-dir, --tree-dir, --copy
可选。同 `filter metrics`，提供后将同步复制筛选通过的 MSA/Tree 文件。

---

### pretree concat

#### --msa-dir, --seq-type
见 Common Parameters。

#### --prefix
输出文件前缀。默认 `matrix`，生成 `matrix.fa`、`matrix.partitions`。

#### --taxa-occupancy
物种占有率和基因筛选。值为 0–1。如设为 `0.75`，仅保留在 >=75% 物种中存在的基因。
阈值越高，矩阵越完整但基因越少。

#### --recoding
氨基酸重编码方案，减少组成偏差。选项：`Dayhoff4`、`Dayhoff6`、`SR4`、`SR6`、`KGB6` 等。

#### --outgroup
外群参考序列。确保外群序列排在矩阵第一位。

#### --to
输出序列类型。在拼接后自动转换为 `AA` 或 `NT`（不支持 `CODON`）。

#### --translate-codon
将 NT 比对翻译为 AA 再拼接。

#### --exclude-codon3
拼接 NT 时排除每个密码子的第三个位点（高度饱和）。

---

## Tree

### tree ml fasttree

#### --msa-dir, --matrix, --seq-type, --tool-args
见 Common Parameters。`--msa-dir`（批量为每个基因推树）和 `--matrix`（单超级矩阵推树）二选一。

#### --model
替换模型名。对蛋白质常用 `LG`、`WAG`；对核苷酸为 `GTR`。默认为 `LG`。
FastTree 不支持复杂模型（无 C10/C20/C60 等 mixture 模型）。

#### --mode
运行模式。`default`（单树）、`gene-trees`（批量基因树）、`pseudo-tree`（指标用途快速推导）。

#### --boot
Bootstrap 伪复制次数。FastTree 使用 SH-like 局部支持值，非标准 bootstrap。默认无。

#### --cat
速率类别数（CAT 近似）。默认 20。越大越接近 GAMMA 连续分布。

#### --gamma
启用 Gamma 速率异质性（更精确但更慢）。关闭则使用 CAT 近似。

#### --fasttree-path
显式 FastTree 可执行文件路径。

---

### tree ml iqtree

#### --msa-dir, --matrix, --seq-type, --tool-args
见 Common Parameters。`--msa-dir`（批量为每个基因推树）和 `--matrix`（单超级矩阵推树）二选一。

#### --model
替换模型名（`LG+F+R4`、`C20+F+G4`、`GTR+F+R6` 等）。
模型选择需结合数据类型和先验知识。使用 `--modelfinder` 可自动搜索最优模型。

#### --state-freq
状态频率类型。`F`（经验频率，推荐）、`FO`（优化频率）、`FC`（计数频率）。

#### --rate-heterogeneity
速率异质性模型。`G4`（Gamma 4 类）、`R4`（FreeRate 4 类）、`I+G4`（不变位点 + Gamma）等。

#### --modelfinder
启用 ModelFinder 自动模型选择。配合 `--mset`/`--msub` 限制搜索空间。

#### --mset
ModelFinder 的候选模型集合。例如 `LG,WAG,JTT` 仅在这三者间选择。

#### --msub
ModelFinder 的候选混合模型设置。`all` 或 `none`。限制 `--msub none` 可跳过 C10/C20 等混合模型搜索。

#### --mode
运行模式。`default`（单树）、`gene-trees`（批量基因树）、`tree-search` 等。

#### --boot
超快 Bootstrap（UFBoot）伪复制次数。标准值 1000，小数据集可用 5000-10000。
仅当 `--bnni` 同时设置时才执行树优化。

#### --alrt
SH-aLRT 检验伪复制次数（默认 1000）。UFBoot >= 95% 且 SH-aLRT >= 80% 的分枝被视为强支持。

#### --bnni
Bootstrap 后执行 NNI 树搜索优化（减少 UFBoot 假阳性）。
建议与 `--boot` 和 `--alrt` 共同使用以获得标准分支支持值。

#### --partitions
分区文件（RAxML 格式或 NEXUS `.best_model.nex`）。映射 IQ-TREE `-p` 和 `-Q`。

#### --rclusterf, --rcluster-max
分区合并参数。IQ-TREE 会合并比例 < `--rclusterf`（默认 0.1）的分区，
直到达到 `--rcluster-max` 个分区。用于防止过度分区。

#### --pmsf-base-model
C10/C20/C60 等 PMSF 混合模型的基础模型（如 `LG`、`WAG`）。
重参数模型（如 C20+F+R4）中的经验频率需匹配此基础模型。

#### --guide-tree
引导树（NEWICK 文件），映射 IQ-TREE `-ft`。用于 PMSF 模型推导。

#### --qmax
混合物类别的搜索上限（默认 10）。设为 100 可搜索全部类别。

#### --constraint
拓扑约束文件（NEWICK 格式）。强制某些节点必须（或不得）出现在结果树中。
格式：`(A,B,(C,D));` 约束 (A,B) 和 (C,D) 必须为姐妹群。

#### --rate
位点速率权重的进化速率值。用于多基因加权分析。

#### --wslr
启用加权似然比检验（支持值稳健性诊断）。

#### --outgroup
外群分类单元名。固定该序列在树的根部（但不会强制外群为最早分化的类群）。

#### --prefix
IQ-TREE 输出前缀。默认与输入文件同名。

#### --iqtree-path
显式 IQ-TREE3 可执行文件路径。

#### --keep-extra
保留 IQ-TREE 生成的全部中间文件（默认仅保留 `.treefile`、`.log`、`.iqtree`）。

---

### tree bi pb

#### --matrix, --threads
见 Common Parameters。

#### --model
替换率矩阵（如 `LG`、`GTR`）。

#### --mixture
位点混合模型：`auto` = CAT Dirichlet 过程，`1` = 单矩阵同质模型，整数 N = 固定 N 组分混合。

#### --gamma-cats
位点速率类别数（4 为推荐默认值）。

#### --start-tree
起始树文件。若无，PhyloBayes 随机生成起始树。

#### --fix-tree
固定拓扑不更新，仅估计分枝长度（最快，适合分歧时间估计的前置分析）。

#### --chains
独立 MCMC 链数量。默认 3，建议至少 2 条以交叉验证收敛性。

#### --chain-prefix, --chain-names
链输出前缀和自定义名称。默认以 `chain` 为前缀。

#### --sample-freq
MCMC 采样频率（每隔 N 次迭代记录一个样本）。

#### --nsamples
每条链的 MCMC 总循环数。`-1` 表示无限运行；使用 `--sample-freq N` 时，保存点数为总循环数 / N。

#### --monitor-freq
监控器输出频率（控制屏幕/日志中的进度更新频率）。

#### --burnin-frac
burnin 比例（0-1）。例如 0.2 表示丢弃前 20% 样本（用于排除 MCMC 尚未收敛的区域）。

#### --poll-interval
轮询间隔（秒）。PhyloAI 每隔 N 秒检查一次 PhyloBayes 运行状态。

#### --pb-path
PhyloBayes-MPI 工具目录。对 `tree bi pb`，目录应包含 `pb_mpi`、`bpcomp`、`tracecomp` 和 `mpirun`。

---

### tree bi bpcomp

#### --chain-dir
包含 `.chain` 文件的目录（必需）。

#### --chain-names
逗号分隔的链名。`all` = 自动发现 `--chain-dir` 中所有 `.chain` 文件。

#### --burnin
丢弃的 saved-sample 数量（整数 ≥ 0）。`0` = 无 burn-in。

#### --sample-freq
burn-in 后的子采样频率（每 N 个样本取 1 个）。

#### --until
停止的样本索引。`all` = 到链末尾；整数 = 到指定 saved-sample 索引。

#### --cutoff
多数规则一致树阈值（0-1，默认 0.5）。后验概率低于此值的节点被折叠。

#### --pb-path
PhyloBayes 工具目录（包含 bpcomp）。

---

### tree bi tracecomp

#### --chain-dir
包含 `.trace` 文件的目录（必需）。

#### --chain-names
逗号分隔的链名。`all` = 自动发现 `--chain-dir` 中所有 `.trace` 文件。

#### --burnin
丢弃的 saved-sample 数量（整数 ≥ 0）。

#### --pb-path
PhyloBayes 工具目录（包含 tracecomp）。

---

### tree bi readpb

#### --chain
无扩展名的链文件路径（必需）。如 `runs/tree/bi/chains/chain1`。

#### --mode
逗号分隔的分析模式（必需）。可选: `rr`, `ss`, `r`, `sitelogl`, `ppred`, `div`, `sitecomp`, `siteconvprob`, `comp`, `allppred`。
`allppred` 与 `div`/`sitecomp`/`siteconvprob`/`comp` 互斥。

#### --burnin
丢弃的 saved-sample 数量（整数 ≥ 0）。

#### --sample-freq
burn-in 后的子采样频率。

#### --until
停止的样本索引。`all` = 到链末尾；整数 = 到指定索引。

#### --threads
MPI 进程数（≥ 2，默认 4）。

#### --pb-path
PhyloBayes 工具目录（包含 readpb_mpi 和 mpirun）。

#### --output-dir
readpb 输出目录。每个模式完成后，PhyloAI 将其生成的文件移入此目录；`allppred` 的 `<chain>.ppred` 直接位于输出目录根部。仅 `ppred` 模式的 `<chain>_ppred*.ali` 位于 `ppred/` 子目录。

#### --mode ss,rr,r
当三个模式同时指定时，PhyloAI 自动生成 `partition.PMSF.nex`。`r` 的 `.meansiterates` 提供从 0 开始的后验平均位点速率；`.trace` 的 `alpha` 列按 `--burnin`、`--sample-freq`、`--until` 取样求均值；`.log` 提供离散 Gamma 类别数。因此标准 G4 链写出 `+G4{alpha}`，并与 `<chain>.exchangeabilities` 和 site-specific frequencies 组合，可用于 `iqtree3 --alisim`。

---

### tree msc

#### --tree, --tree-dir, --mode, --outgroup, --tool-args
见 Common Parameters。`--tree` 提供单基因树文件列表；`--tree-dir` 提供基因树目录。
至少提供一个。

#### --boot
Bootstrap 伪复制次数。wASTRAL 每次复制从每个基因随机抽取一棵树进行种树推断。

#### --extra-rounds
额外轮数。wASTRAL 在标准搜索后追加的额外搜索轮数（提高精度，但更耗时）。

#### --tree-boot-type
Bootstrap 类型。`gene-only`（默认，每次从每个基因树文件中随机选择一棵树）、
`gene-site`（先从每个基因 bootstrap 再建立基因树）、
`site-only`（从每个基因比对重新采样位点）。

#### --tree-boot-min, --tree-boot-max
用于 Bootstrap 的每基因最小/最大树数。

#### --wastral-path
显式 wASTRAL 可执行文件路径。

---

### tree cf

#### --cf
一致性因子类型。`gcf`（基因一致性因子）、`scf`（位点一致性因子）、
`scfl`（每位点对数似然值）、`qcf`（四重奏 CF）。
gCF + sCF 联合使用可区分基因树不完整谱系分选（gCF 低）和信号冲突（sCF 低）。

#### --ref-tree
参考物种树（NEWICK 格式）。所有 CF 值以该树为基准计算。

#### --tree, --tree-dir
与 `tree msc` 相同：基因树输入，至少提供一个。

#### --matrix
超级矩阵比对文件。仅 `sCF` 和 `sCFl` 需要；`gCF` 省略。

#### --partitions
分区文件。仅 `sCF` 和 `sCFl` 需要，与 `--matrix` 配对使用。

#### --model
进化模型。仅 `sCF` 和 `sCFl` 需要。

#### --scf-quartets
sCF 中采样的四重奏数量。默认随数据大小调整；值越大，sCF 越精确但越慢。

#### --lpp
计算位点对数似然值（sCFl 模式标志）。

#### --prefix
输出前缀。默认根据参考物种树名设置。

#### --iqtree-path, --wastral-path
显式 IQ-TREE3 和 wASTRAL 可执行文件路径。

---

## Posttree

### posttree topology

#### --matrix, --guide-tree, --partitions, --iqtree-path, --tool-args, --input-format
见 Common Parameters。

#### --candidate-trees
候选树输入。支持：
1. 单个树列表文件（每行一个 NEWICK 树）
2. 逗号分隔的多个独立树文件（如 `h1.nwk,h2.nwk,h3.nwk`）
PhyloAI 自动合并并去重名称。

#### --model-expr
完�的 IQ-TREE 模型表达式。例如 `LG+F+R4`（同质模型）或 `C20+F+R4`（混合物模型）。
分区模型中省略，使用 `--partitions` 替代。

#### --replicates
RELL Bootstrap 复制次数（默认 10000）。映射 IQ-TREE `-zb`。

#### --prefix
IQ-TREE 输出前缀。默认使用矩阵文件名主干。

---

### posttree dating hessian

#### --matrix, --rooted-tree, --seq-type, --partitions, --iqtree-path, --tool-args
见 Common Parameters。
`--rooted-tree` 必须包含 MCMCtree 校准格式的节点年龄约束（单位：100 Mya），
且必须有根年龄约束。例：`(A,(B,C)'>3.1<3.8')'<4.2';`

#### --model-expr
自定义 IQ-TREE 模型表达式（与 `--partitions` 互斥）。默认：AA = `LG+F+G4`，NT = `GTR+G4`。

---

### posttree dating mcmc

#### --hessian-dir
`phyloai posttree dating hessian` 输出目录（包含 `in.BV` 和校准树文件）。
命令从此处读取 IQ-TREE 的梯度和海森矩阵以及校准信息。

#### --ctl
使用已有 `mcmctree.ctl` 文件（不自动生成）。此模式下 `--clock`/`--burnin`/`--sample-freq`/`--nsamples` 被忽略。

#### --clock
分子钟模型。`1`（全局钟）、`2`（独立速率，推荐）、`3`（相关速率）。

#### --burnin
MCMC burnin 迭代数。默认 100,000。

#### --sample-freq
MCMC 采样间隔（每隔 N 次迭代采样一次）。默认 10。

#### --nsamples
MCMC 保留样本数。默认 10,000。总迭代 = burnin + (sample_freq × nsamples) = 200,000。

#### --runs
独立后验 MCMC 链数量（默认 2，每条链配一条先验链）。--runs=1 跳过收敛诊断。

#### --mcmctree-path
显式 mcmctree（PAML 套件）可执行文件路径。

---

## Report

### --run-dir
要生成报告的运行目录。支持两种结构：
- `pipeline`：`phyloai run` 产生的双层输出（含顶层 `result.json` 和各步骤子目录）
- `module`：单体模块输出（一个或多个步骤子目录）
report 自动检测结构。

### --output-dir, --overwrite, --quiet
见 Common Parameters。报告默认写入 `<run_dir>/report/`。
