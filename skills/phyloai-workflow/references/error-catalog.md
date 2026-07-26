# Error Catalog

Use this for exit code 1 and 3 errors. Exit code 2 is external-tool stderr diagnosis.

## Exit 1: Missing Required Parameter
Pattern: `Missing option`, `is required`
Fix: Re-render the parameter card and request the missing field.

## Exit 1: Output Directory Exists
Pattern: `Output directory exists and is not empty`
Fix: Prefer a new `--output-dir` or `--resume` when available. Use `--overwrite` only after separate explicit confirmation.

Overwrite confirmation template:

```text
`--overwrite` 会删除或替换已有输出目录:
<output-dir>

请明确确认是否覆盖这个目录。未确认前不要执行。
```

## Exit 1: Invalid Enum
Pattern: `Invalid value for`, `Choose from`
Fix: Show allowed enum values from runtime schema.

## Exit 1: File Not Found
Pattern: `does not exist`, `not found`
Fix: Ask user to confirm the path; prefer absolute paths for recovery.

## Exit 1: Empty File
Pattern: `empty`
Fix: Exclude the file or regenerate upstream input.

## Exit 1: Unsupported Format
Pattern: `Unsupported`, `unrecognized format`
Fix: Run `pretree convert` or set `--input-format` if appropriate.

## Exit 1: Overwrite Resume Conflict
Pattern: `--overwrite and --resume are mutually exclusive`
Fix: Choose one: overwrite for fresh run, resume for interrupted run.

## Exit 1: Blocked Tool Args
Pattern: `blocked flag`, `managed by PhyloAI`
Fix: Remove path/control flags from `--tool-args`; use PhyloAI parameters instead.

## Exit 3: Required Tool Missing
Pattern: `not installed`, `not detectable`, `not found`, `Missing required tool`
Fix: Show a concise missing-tool fix card. Include the missing tool name, whether it is required or optional for the requested command, the affected PhyloAI command, why the tool is needed, and a link to `docs/commands/installation.md`.

Template:

```text
缺少工具: <tool>
状态: required | optional
影响命令: <phyloai command>
为什么需要: <short purpose>
下一步:
  1. 查看 docs/commands/installation.md#<tool-or-group-anchor>
  2. 安装后运行 phyloai doctor
  3. 如果已安装但未检测到，检查 PATH 或使用该命令的显式工具路径参数
```

Notes:
- Missing required tools block the requested command and may make `phyloai doctor` exit 3.
- Missing optional tools only make the dependent module unavailable; `phyloai doctor` can still exit 0 with warnings.
- For PhyloBayes-MPI tools, required tools vary by subcommand:
  - `phyloai tree bi pb`: requires `pb_mpi`, `bpcomp`, `tracecomp`, `mpirun`; `readpb_mpi` is optional.
  - `phyloai tree bi bpcomp`: requires `bpcomp`.
  - `phyloai tree bi tracecomp`: requires `tracecomp`.
  - `phyloai tree bi readpb`: requires `readpb_mpi` and `mpirun`.

## Exit 3: Runtime Missing
Pattern: `Java`, `Julia`, `MPI`
Fix: Show the same missing-tool fix card. Link Java and Julia to the Runtime Dependencies section of `docs/commands/installation.md`; link MPI-related failures to the Bayesian Inference section.
