# PhyloAI AI 集成

[English](ai-integration.md) | [中文](ai-integration.zh.md)

PhyloAI 通过 **MCP Server** + **Skill** 提供 AI 交互层，让你可以用对话方式完成系统发育分析，无需记忆 CLI 参数。

```
用户  ←→  AI 客户端（Skill 驱动对话）  ←→  MCP Server（执行桥）  ←→  phyloai CLI  ←→  文件系统
```

| 层 | 角色 |
|----|------|
| **MCP Server** | 将分析、报告、doctor、run 等 CLI 命令以可调用工具形式暴露，工具 schema 自动生成。处理任务启动、状态跟踪、结果读取。 |
| **Skill** | 引导对话：参数确认卡片、结果解读、下一步建议、会话恢复、错误诊断。 |

## 快速开始

### 1. 安装 phyloai

```bash
pip install -e .
```

验证 MCP server 可正常启动：

```bash
phyloai mcp-server --help
```

### 2. 配置 AI 客户端

将 MCP server 添加到客户端配置，然后让 Skill 可被发现。

不要只启用 MCP server 而跳过 Skill。MCP 暴露的是执行工具，而 `phyloai-workflow` 是引导层——它负责运行 `doctor`、渲染参数卡片、等待确认、解读结果、诊断缺失工具。

`AGENTS.md` 是位于项目根目录的 AI 智能体约定文件。打开本仓库作为工作区的 AI 客户端（OpenCode、Claude Code、Codex）会自动加载它并遵循其中的规则——无需用户手动操作。当工作目录不在本仓库内时它不会生效。

如需在仓库外也启用 Skill-first 行为：
- **OpenCode：** 把 `AGENTS.md` 中的规则（或整个文件）复制到 `~/.config/opencode/AGENTS.md`，使其在所有会话中全局生效。
- **Claude Code：** 复制到 `~/.claude/CLAUDE.md`。
- 或者在 AI 客户端的 skills 目录中安装 `skills/phyloai-workflow`（见下方各小节）。

---

## OpenCode（首选）

### MCP Server

添加到 `~/.config/opencode/opencode.jsonc`：

```jsonc
"mcp": {
  // ... existing servers ...
  "phyloai": {
    "type": "local",
    "command": ["phyloai", "mcp-server"],
    "enabled": true
  }
}
```

### Skill

将 phyloai-workflow skill 复制（或软链）到 OpenCode 的 skills 目录：

```bash
cp -r skills/phyloai-workflow ~/.config/opencode/skills/phyloai-workflow
```

或建立软链（与仓库更新保持同步）：

```bash
ln -s "$(pwd)/skills/phyloai-workflow" ~/.config/opencode/skills/phyloai-workflow
```

如果你始终在 phyloai 仓库内工作，可以跳过复制——OpenCode 在工作区为 phyloai 项目根时会自动发现本地的 `skills/` 目录。

### 开始对话

重启 OpenCode 后，试试下面任一指令：

| 你说 | 会发生什么 |
|------|------------|
| "我有蛋白序列，做系统发育分析" | Skill 检测到新任务 → 运行 `doctor` → 询问 `seq-dir` → 引导 pretree 流程 |
| "我的 iqtree 跑完了，下一步？" | Skill 调用 `read_report` 或 `read_result` → 总结已完成步骤 → 建议下一步逻辑命令 |
| "在 ./raw 上跑 pretree stats" | Skill 调用 `get_command_schema("pretree_stats")` → 渲染带中文注解的参数卡片 → 等待你确认 |
| "帮我做系统发育分析" | 同上，对话跟随你的语言 |
| "先给我看个 demo" | Skill 使用内置 demo 数据（20 个基因 × 6 个物种）→ 一步步运行到用户/自动生成的 run 目录 |

关键行为：
- **每个命令都需要显式确认** —— Skill 展示参数卡片并等待你确认后才执行。
- **长时任务是 fire-and-forget** —— 你拿到一个 `output_dir` 句柄；随时调用 `check_status` 查看进度。
- **会话恢复** —— 新建对话后说"恢复我在 runs/run/ 的分析"，Skill 读取 `report.json` 重建完整上下文。

---

## Claude Code

### MCP Server

添加到 `~/.claude/claude_desktop_config.json`：

```json
{
  "mcpServers": {
    "phyloai": {
      "command": "phyloai",
      "args": ["mcp-server"]
    }
  }
}
```

或者在项目根目录的 `.mcp.json` 中添加以做项目级 server。

### Skill

复制到 Claude Code 的 skills 目录：

```bash
cp -r skills/phyloai-workflow ~/.claude/skills/phyloai-workflow
```

### 开始对话

与 OpenCode 相同的交互模式——Skill 定义工作流，与使用哪个客户端无关。

---

## Codex（OpenAI Codex CLI）

### MCP Server

添加到 `~/.codex/config.toml`：

```toml
[mcp_servers.phyloai]
command = "phyloai"
args = ["mcp-server"]
```

### Skill

复制到 Codex 的 skills 目录：

```bash
cp -r skills/phyloai-workflow ~/.codex/skills/phyloai-workflow
```

---

## 内置工具

MCP server 为每个 CLI 子命令注册一个工具，外加四个工具型工具。所有工具 schema 都从 Click 命令树动态生成——无需手动同步。

### CLI 工具

| Tool | CLI equivalent | Execution |
|------|---------------|
| `doctor` | `phyloai doctor --output-format json` | synchronous |
| `pretree_convert` | `phyloai pretree convert` | fire-and-forget |
| `pretree_stats` | `phyloai pretree stats` | fire-and-forget |
| `pretree_align` | `phyloai pretree align` | fire-and-forget |
| `pretree_trim` | `phyloai pretree trim` | fire-and-forget |
| `pretree_metrics` | `phyloai pretree metrics` | fire-and-forget |
| `pretree_filter_taper` | `phyloai pretree filter taper` | fire-and-forget |
| `pretree_filter_treeshrink` | `phyloai pretree filter treeshrink` | fire-and-forget |
| `pretree_filter_metrics` | `phyloai pretree filter metrics` | fire-and-forget |
| `pretree_filter_symtest` | `phyloai pretree filter symtest` | fire-and-forget |
| `pretree_filter_cluster` | `phyloai pretree filter cluster` | fire-and-forget |
| `pretree_concat` | `phyloai pretree concat` | fire-and-forget |
| `tree_ml_fasttree` | `phyloai tree ml fasttree` | fire-and-forget |
| `tree_ml_iqtree` | `phyloai tree ml iqtree` | fire-and-forget |
| `tree_bi_pb` | `phyloai tree bi pb` | fire-and-forget |
| `tree_bi_bpcomp` | `phyloai tree bi bpcomp` | fire-and-forget |
| `tree_bi_tracecomp` | `phyloai tree bi tracecomp` | fire-and-forget |
| `tree_bi_readpb` | `phyloai tree bi readpb` | fire-and-forget |
| `tree_msc` | `phyloai tree msc` | fire-and-forget |
| `tree_cf` | `phyloai tree cf` | fire-and-forget |
| `posttree_topology` | `phyloai posttree topology` | fire-and-forget |
| `posttree_dating_hessian` | `phyloai posttree dating hessian` | fire-and-forget |
| `posttree_dating_mcmc` | `phyloai posttree dating mcmc` | fire-and-forget |
| `posttree_syserror_brlen` | `phyloai posttree syserror brlen` | fire-and-forget |
| `posttree_syserror_brlen_label_nodes` | `phyloai posttree syserror brlen label-nodes` | fire-and-forget |
| `posttree_syserror_cca` | `phyloai posttree syserror cca` | fire-and-forget |
| `report` | `phyloai report` | synchronous |
| `run` | `phyloai run` | fire-and-forget |

占位工具（返回"not yet available"）：`posttree_simulate`、`posttree_syserror_sites`。

### 工具型工具（只读、同步）

| Tool | Description |
|------|-------------|
| `check_status` | 按 `output_dir` 检查任务状态。返回 `not_started` / `running` / `success` / `error` / `unknown`。 |
| `read_result` | 读取某步骤输出目录中的 `result.json`。 |
| `read_report` | 读取 run 目录下的 `report.json`（在 `<run_dir>/report/` 下查找）。 |
| `get_command_schema` | 获取任意工具的运行时参数 schema（名称、类型、默认值、可选值、说明）。 |

## Skill 参考文件

`skills/phyloai-workflow/` 中的 Skill 包含：

| File | Content |
|------|---------|
| `SKILL.md` | 核心规则、入口模式、工作流阶段、错误处理策略 |
| `references/parameter-annotations.md` | 常用参数的中文科学注解；运行时 schema（来自 `get_command_schema`）是权威来源；省略运维与可视化参数 |
| `references/error-catalog.md` | Exit 1/3 的已知错误模式及修复卡片 |
| `references/dialog-templates.md` | 参数卡片、结果卡片、恢复卡片、demo 提示词 |
| `references/demo-data.md` | 内置 demo 数据集路径及各步骤入口 |
| `references/workflow.md` | 各阶段的执行参考：输入/输出/检查/下一步指南 |

## 架构说明

- **传输方式**：stdio（本地进程）。无 HTTP server、无端口、无鉴权。
- **Schema 来源**：所有 MCP 工具的参数 schema 都在 server 启动时通过内省 Click 命令树生成。CLI 参数变更会自动反映——零维护。
- **任务句柄**：所有 fire-and-forget 命令都返回 `output_dir`（绝对路径）作为持久任务句柄。跨会话用 `check_status` 跟踪。
- **Result.json**：每个 CLI 命令都写结构化 JSON 结果。Skill 用它做结果解读；MCP server 用它做工具响应。
- **Skill 版本**：Skill 位于仓库内，与 CLI 版本绑定。参数注解更新与 CLI 变更一起发布。