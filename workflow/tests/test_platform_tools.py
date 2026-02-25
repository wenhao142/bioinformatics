from pathlib import Path

from workflow.platform_tools import resolve_tool_command


def test_resolve_tool_command_defaults_to_tool_name():
    assert resolve_tool_command("prefetch", "") == "prefetch"
    assert resolve_tool_command("fasterq-dump", None) == "fasterq-dump"


def test_resolve_tool_command_accepts_explicit_binary_path(tmp_path: Path):
    exe = tmp_path / "prefetch"
    exe.write_text("", encoding="utf-8")
    assert resolve_tool_command("prefetch", str(exe)) == str(exe)


def test_resolve_tool_command_resolves_directory_binary(tmp_path: Path):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    cmd = bin_dir / "fasterq-dump"
    cmd.write_text("", encoding="utf-8")
    assert resolve_tool_command("fasterq-dump", str(bin_dir)) == str(cmd)


def test_resolve_tool_command_resolves_directory_windows_suffix(tmp_path: Path):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    cmd = bin_dir / "prefetch.exe"
    cmd.write_text("", encoding="utf-8")
    assert resolve_tool_command("prefetch", str(bin_dir)) == str(cmd)
