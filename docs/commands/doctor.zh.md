# phyloai doctor

[English](doctor.md) | [中文](doctor.zh.md)

## 目的

`phyloai doctor` 检查本地 PhyloAI 环境所依赖的外部工具是否能被找到、位于何处、以及是否能解析出版本字符串。

它不会安装工具、修改环境，也不会校验输入数据集。

## 用法

```bash
phyloai doctor [--output-format text|json]
```

## 参数

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--output-format text\|json` | `text` | 选择人类可读终端输出或机器可读 JSON。`doctor` 是唯一默认输出 `text` 的命令。 |

## 输入

`doctor` 不读取序列或树输入文件。它读取当前 shell 环境（含 `PATH`），并检查已知外部工具是否在该环境中可见。

所有第三方工具都依赖环境，应由用户根据其操作系统、包管理器、集群或 Conda 环境自行安装。BMGE 检测名称为 `BMGE.jar`；TAPER 检测名称为 `correction_multi.jl`。

实际安装指南见 [installation.md](installation.md)。它列出了 Python 环境选项、外部工具组、`phyloai run` 依赖模式以及操作系统相关注意事项。

## 输出

文本输出是一个 Rich 表格，展示每个工具、状态、检测到的版本、路径，以及相关的安装提示。

JSON 输出是从工具名到结构化状态字段的映射，如 `status`、`path`、`version` 和 `note`。

## 示例

人类可读的环境检查：

```bash
phyloai doctor
```

供脚本或 CI 使用的 JSON：

```bash
phyloai doctor --output-format json
```

## 警告与错误

缺失某个可选工具会被报告，但不意味着整个安装不可用。后续依赖该工具的工作流步骤可能会失败或不可用。

如果 `doctor` 中某个工具缺失，说明 PhyloAI 看到的环境与后续 CLI 命令所看到的一致。请先激活目标 Conda 或虚拟环境，再重新运行该命令。

如果 `TAPER`（`correction_multi.jl`）或 `BMGE`（`BMGE.jar`）报告缺失，请自行安装或在使用到它们的命令中显式传入路径。

## 备注

当前注册表包含必需工具：`iqtree3`、`mafft`、`trimal`；以及可选工具：`wastral`、`pb_mpi`、`mcmctree`、`run_treeshrink.py`、`magus`、`clipkit`、`BMGE.jar`、`correction_multi.jl`、`java`、`julia`、`FastTree`。