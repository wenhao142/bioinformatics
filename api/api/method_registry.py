import copy
import re
import time
from typing import Any

_SEMVER_RE = re.compile(r"^\d+\.\d+\.\d+(?:[-+][A-Za-z0-9.-]+)?$")

# plugin_id -> version -> manifest_record
_REGISTRY: dict[str, dict[str, dict[str, Any]]] = {}


def _semver_key(version: str) -> tuple[int, int, int, int, str]:
    # Semver precedence ignores build metadata (`+...`) and treats release > prerelease.
    base = version.split("+", 1)[0]
    release = base
    prerelease = ""
    if "-" in base:
        release, prerelease = base.split("-", 1)
    major_s, minor_s, patch_s = release.split(".")
    release_flag = 1 if not prerelease else 0
    return (int(major_s), int(minor_s), int(patch_s), release_flag, prerelease)


def _set_latest_flags(plugin_id: str) -> None:
    versions = _REGISTRY.get(plugin_id, {})
    if not versions:
        return
    latest = list_versions(plugin_id)[0]
    for version, record in versions.items():
        record["is_latest"] = version == latest


def register_method(record: dict[str, Any]) -> dict[str, Any] | None:
    plugin_id = str(record.get("plugin_id", "")).strip()
    version = str(record.get("version", "")).strip()
    if not plugin_id or not version or not _SEMVER_RE.match(version):
        raise ValueError("Invalid plugin_id or version")

    versions = _REGISTRY.setdefault(plugin_id, {})
    if version in versions:
        return None

    stored = copy.deepcopy(record)
    stored["registered_at"] = time.time()
    versions[version] = stored
    _set_latest_flags(plugin_id)
    return copy.deepcopy(versions[version])


def get_method(plugin_id: str, version: str | None = None) -> dict[str, Any] | None:
    versions = _REGISTRY.get(plugin_id)
    if not versions:
        return None
    if version:
        record = versions.get(version)
        return copy.deepcopy(record) if record else None
    latest_version = list_versions(plugin_id)[0]
    return copy.deepcopy(versions[latest_version])


def list_versions(plugin_id: str) -> list[str]:
    versions = _REGISTRY.get(plugin_id, {})
    return sorted(versions.keys(), key=_semver_key, reverse=True)


def list_methods(*, latest_only: bool = True) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for plugin_id in sorted(_REGISTRY.keys()):
        versions = _REGISTRY[plugin_id]
        if latest_only:
            latest_version = list_versions(plugin_id)[0]
            rows.append(copy.deepcopy(versions[latest_version]))
            continue
        for version in list_versions(plugin_id):
            rows.append(copy.deepcopy(versions[version]))
    return rows


def set_enabled(plugin_id: str, enabled: bool, version: str | None = None) -> dict[str, Any] | None:
    versions = _REGISTRY.get(plugin_id)
    if not versions:
        return None
    target_version = version or list_versions(plugin_id)[0]
    record = versions.get(target_version)
    if record is None:
        return None
    record["enabled"] = enabled
    return copy.deepcopy(record)


def clear_registry() -> None:
    _REGISTRY.clear()
