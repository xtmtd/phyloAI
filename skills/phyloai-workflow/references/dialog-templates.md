# Dialog Templates

## Parameter Card

```text
命令: phyloai <command>
目的: <one-line purpose>

参数:
  --name    <value>    <中文说明> [推荐: <value if applicable>]

Schema source: runtime CLI via get_command_schema

确认执行？还是需要调整参数？
```

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
