# Parameter Annotations

Runtime schema from `get_command_schema` is authoritative. These notes add Chinese scientific context only.

## Common Parameters

### --output-dir
输出目录。建议每一步使用独立目录，避免覆盖已有结果。

### --threads
并行线程数。推荐设为可用 CPU 的一半到全部；共享服务器上避免占满资源。

### --overwrite
覆盖已有输出。只在确认旧结果不需要时使用。

### --resume
从检查点恢复长任务。不要和 `--overwrite` 同时使用。

### --tool-args
传递给底层工具的策略参数。不要放输入/输出路径；这些由 PhyloAI 管理。

## Frequently Used Commands

### pretree align --method
多序列比对算法。`linsi` 精确但慢，适合小型数据；`auto` 适合大型数据。

### pretree trim --tool
修剪工具。BMGE 偏保守，trimAl/ClipKIT 可用于不同缺失和保守性策略。

### pretree concat --taxa-occupancy
保留达到指定物种占有率的基因。阈值越高，矩阵更完整但基因更少。

### tree ml iqtree --model
最大似然模型。模型选择应结合数据类型、分区和前期分析结果。

### tree ml fasttree
快速树推断。适合探索或大批量基因树，不替代最终高精度 ML 分析。

### posttree dating hessian
为 MCMCtree 近似似然计算梯度和 Hessian。需要带校准的有根树。

### posttree dating mcmc
执行贝叶斯分歧时间估计。重点检查校准、链长、采样频率和收敛诊断。

### report
生成 `report.json` 和 `report.html`，用于恢复会话、归档和 Methods 初稿。
