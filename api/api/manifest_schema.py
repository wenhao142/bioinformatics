from typing import Any

from jsonschema import Draft202012Validator

PLUGIN_MANIFEST_SCHEMA: dict[str, Any] = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "type": "object",
    "additionalProperties": False,
    "required": [
        "plugin_id",
        "name",
        "version",
        "image",
        "input_schema",
        "output_schema",
    ],
    "properties": {
        "plugin_id": {"type": "string", "pattern": "^[a-z0-9][a-z0-9._-]{1,63}$"},
        "name": {"type": "string", "minLength": 1, "maxLength": 120},
        "version": {"type": "string", "pattern": r"^\d+\.\d+\.\d+(?:[-+][A-Za-z0-9.-]+)?$"},
        "image": {"type": "string", "minLength": 3, "maxLength": 512},
        "description": {"type": ["string", "null"], "maxLength": 2000},
        "input_schema": {"$ref": "#/$defs/io_schema"},
        "output_schema": {"$ref": "#/$defs/io_schema"},
        "resources": {
            "type": "object",
            "additionalProperties": False,
            "properties": {
                "cpu_millicores": {"type": "integer", "minimum": 100, "maximum": 64000},
                "memory_mb": {"type": "integer", "minimum": 128, "maximum": 262144},
                "gpu": {"type": "boolean"},
            },
        },
        "enabled": {"type": "boolean"},
        "deprecated": {"type": "boolean"},
        "deprecation_message": {"type": ["string", "null"], "maxLength": 400},
        "superseded_by": {"type": ["string", "null"], "maxLength": 64},
        "adapter": {
            "type": ["object", "null"],
            "additionalProperties": False,
            "required": ["from_type", "to_type"],
            "properties": {
                "from_type": {"type": "string", "minLength": 1, "maxLength": 120},
                "to_type": {"type": "string", "minLength": 1, "maxLength": 120},
                "strategy": {"type": "string", "minLength": 1, "maxLength": 80},
            },
        },
        "tags": {
            "type": "array",
            "maxItems": 30,
            "items": {"type": "string", "maxLength": 64},
        },
    },
    "$defs": {
        "io_schema": {
            "type": "object",
            "required": ["type", "properties"],
            "properties": {
                "type": {"const": "object"},
                "properties": {
                    "type": "object",
                    "minProperties": 1,
                    "additionalProperties": {"type": "object"},
                },
                "required": {
                    "type": "array",
                    "items": {"type": "string"},
                    "uniqueItems": True,
                },
                "additionalProperties": {
                    "anyOf": [
                        {"type": "boolean"},
                        {"type": "object"},
                    ]
                },
            },
            "additionalProperties": True,
        }
    },
}

_VALIDATOR = Draft202012Validator(PLUGIN_MANIFEST_SCHEMA)


def validate_plugin_manifest_or_raise(payload: dict[str, Any]) -> None:
    errors = sorted(_VALIDATOR.iter_errors(payload), key=lambda e: list(e.path))
    if not errors:
        return
    err = errors[0]
    if err.path:
        path = ".".join(str(p) for p in err.path)
        raise ValueError(f"manifest.{path}: {err.message}")
    raise ValueError(f"manifest: {err.message}")
