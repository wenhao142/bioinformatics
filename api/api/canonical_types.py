from typing import Any

CANONICAL_TYPE_REGISTRY_VERSION = "1.0.0"

_CANONICAL_TYPES: list[dict[str, Any]] = [
    {
        "id": "reads.fastq.gz",
        "description": "Compressed FASTQ reads",
        "extensions": [".fastq.gz", ".fq.gz"],
    },
    {
        "id": "align.bam",
        "description": "Coordinate-sorted BAM alignment file",
        "extensions": [".bam"],
        "requires_sidecars": [".bai"],
    },
    {
        "id": "variants.vcf.gz",
        "description": "Compressed VCF variant file",
        "extensions": [".vcf.gz"],
        "requires_sidecars": [".tbi"],
    },
    {
        "id": "expression.counts.tsv",
        "description": "Raw expression counts table",
        "extensions": [".tsv"],
    },
    {
        "id": "expression.diff_table.tsv",
        "description": "Differential expression results table",
        "extensions": [".tsv"],
    },
    {
        "id": "report.html",
        "description": "HTML report artifact",
        "extensions": [".html"],
    },
]

_CANONICAL_IDS = {entry["id"] for entry in _CANONICAL_TYPES}


def list_canonical_types() -> list[dict[str, Any]]:
    return [dict(entry) for entry in _CANONICAL_TYPES]


def is_known_canonical_type(type_id: str) -> bool:
    return type_id in _CANONICAL_IDS


def validate_schema_canonical_types(schema_name: str, schema: dict[str, Any]) -> None:
    properties = schema.get("properties", {})
    if not isinstance(properties, dict):
        return
    for field_name, field_spec in properties.items():
        if not isinstance(field_spec, dict):
            continue
        canonical_type = field_spec.get("canonical_type")
        if canonical_type is None:
            continue
        if not isinstance(canonical_type, str) or not canonical_type.strip():
            raise ValueError(f"{schema_name}.properties.{field_name}.canonical_type must be a non-empty string")
        if not is_known_canonical_type(canonical_type):
            raise ValueError(
                f"{schema_name}.properties.{field_name}.canonical_type is unknown: {canonical_type}. "
                f"Known types: {sorted(_CANONICAL_IDS)}"
            )
