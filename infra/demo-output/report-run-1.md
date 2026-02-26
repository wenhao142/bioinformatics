# Analysis Report

## Summary
- Project: `demo-project`
- Run ID: `run-1`
- Kind: `evidence_rank`
- Created by: `admin@example.com`
- Created at (UTC): `2026-02-26T01:46:27.527792+00:00`

## Top Genes
| Rank | Gene | Score |
| --- | --- | ---: |
| 1 | GENE1 | 7.5031 |
| 2 | GENE2 | 5.8709 |

## Top Loci Evidence
| Item | Chr | Start/Pos | End/Alt | Score |
| --- | --- | --- | --- | ---: |
| GENE1 | chr1 | 80 | 140 | 7.5031 |
| GENE2 | chr1 | 700 | 760 | 5.8709 |

## Research Direction (Inference)
- Disease: `Alzheimer disease`
- Inference label: `inference`
- Summary: Offline template summary for Alzheimer disease: candidate genes are GENE1, GENE2. Interpret as inference pending literature and lab validation.

### Hypotheses
- H1 (GENE1): GENE1 may modulate Alzheimer disease risk through transcriptomic and variant-linked mechanisms.
- H2 (GENE2): GENE2 may modulate Alzheimer disease risk through transcriptomic and variant-linked mechanisms.

### Experiment Plans
- E1 (GENE1): Targeted validation plan for GENE1
- E2 (GENE2): Targeted validation plan for GENE2

### Evidence Citations
| ID | Type | Label |
| --- | --- | --- |
| C1 | input-gene | Requested gene GENE1 |
| C2 | input-gene | Requested gene GENE2 |

## Run Metadata
```json
{
  "counts": {
    "omics_rows": 4,
    "variants_in_region": 5
  },
  "input_hashes": {
    "omics_table_sha256": "aa88c368a69613b37028d76f2d11cf7445970d21c5f6cc7b342e664519aba7d0",
    "variants_region_sha256": "17e75bf844a2ef69cbf57dde22800967cb25fb8a89005f9aa5470fa9d64c6072"
  },
  "params": {
    "chr": "chr1",
    "end": 1000,
    "start": 1,
    "top_n": 5
  },
  "rerun_of": null,
  "result_hash": "bc6037a8339a03b9fd8985ffdb0804381cc3785ddbaaae0c04b973b0e6a465b6",
  "selected_method": "evidence-rank",
  "signature_hash": "eee334fc8cdf77fa7a9dd2420690a36149e1d1b2f7f9c760fd69672456797b50",
  "stable_with_previous": false,
  "tool_versions": {
    "fastapi": "0.110.0",
    "pydantic": "2.6.1",
    "pyjwt": "2.8.0",
    "python": "3.11.14",
    "scoring": "evidence-v1"
  }
}
```
