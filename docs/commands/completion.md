# phyloai completion

[English](completion.md) | [中文](completion.zh.md)


## Purpose

`phyloai completion <shell>` generates static shell completion scripts for Bash, Zsh, and Fish.

Generate the script once from an environment where `phyloai` is installed, save it to a persistent file, and configure your shell to load that saved script.

Do not run `phyloai completion ...` dynamically from `.bashrc`, `.zshrc`, `config.fish`, or other shell startup files.

## Usage

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

Add this line to `~/.bashrc`:

```bash
source ~/.config/phyloai/completion/phyloai.bash
```

If you only run the `source` command manually in the current terminal, completion only works for that shell session.

## Zsh

```bash
mkdir -p ~/.config/phyloai/completion
phyloai completion zsh > ~/.config/phyloai/completion/phyloai.zsh
```

Add this line to `~/.zshrc`:

```bash
source ~/.config/phyloai/completion/phyloai.zsh
```

If you only run the `source` command manually in the current terminal, completion only works for that shell session.

## Fish

```bash
mkdir -p ~/.config/fish/completions
phyloai completion fish > ~/.config/fish/completions/phyloai.fish
```

Fish loads completion files from `~/.config/fish/completions/` automatically in new shells, so no extra `source` line is required.
