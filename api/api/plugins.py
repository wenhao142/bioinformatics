from typing import Any

from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel, Field

from api.auth import current_user

router = APIRouter(prefix="/plugins", tags=["plugins"])

_MANIFESTS: dict[str, dict[str, Any]] = {}


class PluginResources(BaseModel):
    cpu_millicores: int = Field(default=1000, ge=100, le=64000)
    memory_mb: int = Field(default=1024, ge=128, le=262144)
    gpu: bool = False


class PluginManifest(BaseModel):
    plugin_id: str = Field(pattern=r"^[a-z0-9][a-z0-9._-]{1,63}$")
    name: str = Field(min_length=1, max_length=120)
    version: str = Field(pattern=r"^\d+\.\d+\.\d+(?:[-+][A-Za-z0-9.-]+)?$")
    image: str = Field(min_length=3, max_length=512)
    description: str | None = Field(default=None, max_length=2000)
    input_schema: dict[str, Any]
    output_schema: dict[str, Any]
    resources: PluginResources = PluginResources()
    enabled: bool = True
    tags: list[str] = Field(default_factory=list, max_length=30)


def _validate_json_schema(name: str, schema: dict[str, Any]) -> None:
    if not isinstance(schema, dict):
        raise HTTPException(status_code=400, detail=f"{name} must be a JSON object")
    if schema.get("type") != "object":
        raise HTTPException(status_code=400, detail=f"{name}.type must be 'object'")

    properties = schema.get("properties")
    if not isinstance(properties, dict) or not properties:
        raise HTTPException(status_code=400, detail=f"{name}.properties must be a non-empty object")

    allowed_types = {"string", "number", "integer", "boolean", "object", "array"}
    for field_name, field_spec in properties.items():
        if not isinstance(field_spec, dict):
            raise HTTPException(status_code=400, detail=f"{name}.properties.{field_name} must be an object")
        field_type = field_spec.get("type")
        if field_type not in allowed_types:
            raise HTTPException(
                status_code=400,
                detail=f"{name}.properties.{field_name}.type must be one of {sorted(allowed_types)}",
            )

    required = schema.get("required", [])
    if required is not None:
        if not isinstance(required, list) or any(not isinstance(item, str) for item in required):
            raise HTTPException(status_code=400, detail=f"{name}.required must be a list of strings")
        missing = [item for item in required if item not in properties]
        if missing:
            raise HTTPException(status_code=400, detail=f"{name}.required contains unknown fields: {missing}")

    additional = schema.get("additionalProperties", True)
    if not isinstance(additional, (bool, dict)):
        raise HTTPException(status_code=400, detail=f"{name}.additionalProperties must be bool or object")


def _admin_only(user: dict[str, Any]) -> None:
    if user["role"] != "admin":
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Admin only")


@router.post("/register", status_code=201)
def register_plugin(manifest: PluginManifest, user=Depends(current_user)):
    _admin_only(user)
    if manifest.plugin_id in _MANIFESTS:
        raise HTTPException(status_code=409, detail="Plugin already registered")

    _validate_json_schema("input_schema", manifest.input_schema)
    _validate_json_schema("output_schema", manifest.output_schema)

    record = manifest.model_dump()
    record["registered_by"] = user["email"]
    _MANIFESTS[manifest.plugin_id] = record
    return {"plugin": record}


@router.get("")
def list_plugins(user=Depends(current_user)):
    return {"plugins": list(_MANIFESTS.values())}


@router.get("/{plugin_id}")
def get_plugin(plugin_id: str, user=Depends(current_user)):
    plugin = _MANIFESTS.get(plugin_id)
    if plugin is None:
        raise HTTPException(status_code=404, detail="Plugin not found")
    return {"plugin": plugin}
