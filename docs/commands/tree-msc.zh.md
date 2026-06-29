# phyloai tree msc

[English](tree-msc.md) | [中文](tree-msc.zh.md)

使用 [wASTRAL](https://github.com/chaoszhang/ASTER)（ASTER）进行多物种溯祖物种树推断。

## 目的

`phyloai tree msc` 消费基因树并使用 wASTRAL 产出带局部后验概率分支支持的物种树。wASTRAL 是 ASTRAL 的重新实现，用于多物种溯祖模型下的物种树推断。

wASTRAL 计算是一次性完成的 —— 不支持 `--resume`。

## 用法

```bash
# 单个基因树文件输入
phyloai tree msc --tree gene_trees.trees -o runs/tree/msc

# 基因树文件目录（自动合并）
phyloai tree msc --tree-dir ./genetrees/

# 传统未加权 Astral 与穷举搜索
phyloai tree msc --tree-dir ./genetrees/ --mode 4 -R

# 使用自定义范围的 Bootstrap 输入支持
phyloai tree msc --tree-dir ./genetrees/ \
    --mode 1 --boot 2 -R \
    --tree-boot-type bootstrap --tree-boot-min 10 --tree-boot-max 95 \
    -t 8 -o runs/tree/msc

# 通过 --tool-args 覆盖
phyloai tree msc --tree input.trees --tool-args "-r 32 -s 32"

# 使用外类群物种作为根
phyloai tree msc --tree-dir ./genetrees/ --outgroup Oryza_sativa
```

## 输入

| Option | Description |
|--------|-------------|
| `--tree` | 单个基因树文件（newick，每行一棵树）。与 `--tree-dir` 互斥。 |
| `--tree-dir` | 基因树文件目录。扫描 `.nwk`、`.tre`、`.tree`、`.nw`、`.trees`、`.newick` 扩展名，合并为一个输入。与 `--tree` 互斥。 |

## 参数

| Option | Default | Description |
|--------|---------|-------------|
| `--mode` | 1 | 1=hybrid，2=分支支持加权，3=分支长度加权，4=传统未加权 |
| `--boot` | 1 | wastral -u。0=仅拓扑，1=局部后验概率，2=四分体+local-PP，3=2 + freqQuad.csv |
| `-R` / `--extra-rounds` | off | 启用穷举搜索（wastral -R）。 |
| `--tree-boot-type` | auto | 基因树分支支持类型预设：`auto`（检测），`likelihood`（wastral -L/--lrt），`abayes`（wastral -B/--bayes），`bootstrap`（wastral -S/--bootstrap）。设置预设的 -x/-n 值。 |
| `--tree-boot-min` | -- | 最小支持值（wastral -n）。覆盖预设默认值。 |
| `--tree-boot-max` | -- | 最大支持值（wastral -x）。覆盖预设默认值。 |
| `--outgroup` | -- | 用于定根的外类群物种名（wastral --root）。 |

## 输出

```
runs/tree/msc/
├── result.json            # PhyloAI 结构化结果（stderr 内联于 data.tool_stderr）
├── wastral.tre            # 物种树输出（newick）
├── merged.trees           # 合并后的输入（仅 `--tree-dir` 模式）
└── freqQuad.csv           # 四分体频率数据（仅 `--boot 3`）
```

## 退出码

| Code | Meaning |
|------|---------|
| 0 | 成功 |
| 1 | 用户输入错误 |
| 2 | wastral 执行失败 |
| 3 | 未找到 wastral |

## 备注

- wASTRAL 必须安装并在 PATH 上（或使用 `--wastral-path`）。
- wASTRAL 的 stderr 在 `result.json` 中作为 `data.tool_stderr` 内联。不写单独的日志文件。
- `--boot 2` 计算四分体支持 + local-PP，并将值嵌入输出树；不写单独的数据文件。使用 `--boot 3` 获取 `freqQuad.csv`。
- `--tree-dir` 模式将所有有效基因树文件合并为一个输入文件，保存为 `merged.trees`。
- `--tool-args` 透传额外标志给 wastral。`-i` 与 `-o` 被阻止。策略标志覆盖 phyloAI 默认值。
- `--outgroup` 指定单个物种名对树进行定根（wastral `--root`）。