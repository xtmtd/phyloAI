# Dialog Templates

## Parameter Card

```text
命令: phyloai <command>
目的: <one-line purpose>

参数:
  --paramA      <value>    (默认: <default>)  <中文说明> [推荐: <value>]
  --paramB      <value>    (默认: <default>)  <中文说明 or --help text>
  --paramC      <value>    (默认: <default>)  <中文说明 or --help text>
  ...

Schema source: runtime CLI via get_command_schema

确认执行？还是需要调整参数？
```

Rules:
- List **every** parameter from `get_command_schema`. Do not omit any.
- Show the schema default in parentheses for every parameter.
- If the schema marks a parameter as required, it MUST have a value before approval — do not launch with it empty.
- Parameters with annotations in `references/parameter-annotations.md` get Chinese descriptions.
- Parameters without annotations use their CLI `--help` text verbatim.
- If `--overwrite true` is present, add a separate line after the parameter block: `WARNING: --overwrite 将删除 <output-dir>。请单独确认覆盖。`

## Result Card

```text
状态: <success/error>
耗时: <wall_time>
关键结果: <key_results summary>
警告: <warnings if any>
建议下一步: <options>
```

## Recovery Card

```text
已完成: <steps>
运行中: <steps>
失败: <steps/errors>
建议下一步: <next action>
```

## Demo Prompt

```text
如果需要，我可以先用内置示例数据演示这个步骤，然后再切回您的数据。
```
