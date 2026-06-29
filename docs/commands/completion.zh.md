# phyloai completion

[English](completion.md) | [中文](completion.zh.md)

## 目的

`phyloai completion <shell>` 为 Bash、Zsh、Fish 生成静态的 shell 补全脚本。

请在已安装 `phyloai` 的环境中生成一次脚本，将其保存为持久文件，并配置你的 shell 加载该脚本。

不要在 `.bashrc`、`.zshrc`、`config.fish` 等 shell 启动文件中动态运行 `phyloai completion ...`。

## 用法

```bash
phyloai completion bash
phyloai completion zsh
phyloai completion fish
```

## Bash

```bash
mkdir -p ~/.config/phyloai/completion
phyloai completion bash > ~/.config/phyloai/completion/phyloai.bash
```

把这一行加入 `~/.bashrc`：

```bash
source ~/.config/phyloai/completion/phyloai.bash
```

如果只在当前终端手动执行 `source`，补全只对该 shell 会话生效。

## Zsh

```bash
mkdir -p ~/.config/phyloai/completion
phyloai completion zsh > ~/.config/phyloai/completion/phyloai.zsh
```

把这一行加入 `~/.zshrc`：

```bash
source ~/.config/phyloai/completion/phyloai.zsh
```

如果只在当前终端手动执行 `source`，补全只对该 shell 会话生效。

## Fish

```bash
mkdir -p ~/.config/fish/completions
phyloai completion fish > ~/.config/fish/completions/phyloai.fish
```

Fish 会在新 shell 中自动从 `~/.config/fish/completions/` 加载补全文件，无需额外的 `source` 行。