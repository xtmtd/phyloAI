from __future__ import annotations

from phyloai.mcp.server import create_server


def test_server_creates_without_error() -> None:
    server = create_server()

    assert server is not None
