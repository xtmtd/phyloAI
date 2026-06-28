"""Generate MCP tool schemas from the PhyloAI Click command tree."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import click

_EXCLUDED_TOOL_NAMES = {"mcp-server", "completion_bash", "completion_zsh", "completion_fish"}


def walk_click_tree(root: click.Group) -> list[dict[str, Any]]:
    """Return leaf Click commands as MCP command descriptors."""
    descriptors: list[dict[str, Any]] = []

    def walk(group: click.Group, parts: list[str]) -> None:
        for name in group.list_commands(None):
            command = group.get_command(None, name)
            if command is None:
                continue
            next_parts = parts + [name]
            if isinstance(command, click.Group):
                if command.invoke_without_command and command.callback is not None:
                    descriptor = _descriptor(next_parts, command)
                    if descriptor is not None:
                        descriptors.append(descriptor)
                walk(command, next_parts)
            else:
                descriptor = _descriptor(next_parts, command)
                if descriptor is not None:
                    descriptors.append(descriptor)

    walk(root, [])
    return descriptors


def _descriptor(parts: list[str], command: click.Command) -> dict[str, Any] | None:
    tool_name = "_".join(parts)
    if tool_name in _EXCLUDED_TOOL_NAMES:
        return None
    return {
        "tool_name": tool_name,
        "command_path": parts,
        "click_command": command,
        "help": command.help or "",
    }


def click_param_to_json_schema(param: click.Parameter) -> dict[str, Any]:
    """Convert a Click parameter into a JSON-schema property."""
    schema: dict[str, Any] = {"description": getattr(param, "help", None) or ""}

    if isinstance(param, click.Option) and param.is_flag:
        schema["type"] = "boolean"
    elif isinstance(param.type, (click.IntRange, click.types.IntParamType)):
        schema["type"] = "integer"
    elif isinstance(param.type, (click.FloatRange, click.types.FloatParamType)):
        schema["type"] = "number"
    elif isinstance(param.type, click.Choice):
        schema["type"] = "string"
        schema["enum"] = list(param.type.choices)
    elif isinstance(param.type, (click.Path, click.File)):
        schema["type"] = "string"
        schema["format"] = "path"
    else:
        schema["type"] = "string"

    if getattr(param, "default", None) is not None and param.default is not ...:
        from click._utils import Sentinel

        if isinstance(param.default, Sentinel):
            pass
        else:
            default = param.default
            if isinstance(default, Path):
                default = str(default)
            schema["default"] = default
    return schema


def build_mcp_tool(descriptor: dict[str, Any]) -> dict[str, Any]:
    """Build one MCP tool definition from a Click command descriptor."""
    command: click.Command = descriptor["click_command"]
    properties: dict[str, dict[str, Any]] = {}
    required: list[str] = []
    for param in command.params:
        if getattr(param, "hidden", False):
            continue
        properties[param.name] = click_param_to_json_schema(param)
        if param.required:
            required.append(param.name)
    return {
        "name": descriptor["tool_name"],
        "description": descriptor["help"] or command.help or "",
        "inputSchema": {"type": "object", "properties": properties, "required": required},
    }
