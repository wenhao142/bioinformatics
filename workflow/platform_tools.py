from pathlib import Path


def resolve_tool_command(tool_name: str, configured: str | None = None) -> str:
    """
    Resolve a tool command from optional config input.

    Supported config forms:
    - None/"": fallback to bare executable name
    - directory path: <dir>/<tool_name> or <dir>/<tool_name>.exe
    - executable path: use as-is
    """
    raw = (configured or "").strip()
    if not raw:
        return tool_name

    base = Path(raw)
    if base.is_file():
        return str(base)

    if base.is_dir():
        native = base / tool_name
        win = base / f"{tool_name}.exe"
        if native.exists():
            return str(native)
        if win.exists():
            return str(win)

    # Allow users to pass command names directly (e.g., "prefetch", "prefetch.exe")
    return raw
