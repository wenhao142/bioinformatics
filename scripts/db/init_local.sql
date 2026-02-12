-- Minimal local schema for offline genomics experiments.
-- Safe for free-tier/local use; adjust for production (roles, tablespaces, retention).

CREATE EXTENSION IF NOT EXISTS btree_gist;
CREATE EXTENSION IF NOT EXISTS vector;

-- Samples/metadata
CREATE TABLE IF NOT EXISTS samples (
    sample_id      TEXT PRIMARY KEY,
    project        TEXT,
    sex            TEXT,
    ancestry       TEXT,
    platform       TEXT,
    created_at     TIMESTAMPTZ DEFAULT now()
);

-- Variants (per-sample optional)
CREATE TABLE IF NOT EXISTS variants (
    variant_id     BIGSERIAL PRIMARY KEY,
    chr            TEXT NOT NULL,
    pos            INTEGER NOT NULL,
    ref            TEXT NOT NULL,
    alt            TEXT NOT NULL,
    qual           DOUBLE PRECISION,
    filter_status  TEXT,
    info           JSONB DEFAULT '{}'::JSONB,
    sample_id      TEXT REFERENCES samples(sample_id),
    af             DOUBLE PRECISION,
    ac             INTEGER,
    dp             INTEGER,
    gt             TEXT,
    embedding      VECTOR(768), -- optional for later ML; can stay NULL
    created_at     TIMESTAMPTZ DEFAULT now(),
    updated_at     TIMESTAMPTZ DEFAULT now()
);

CREATE INDEX IF NOT EXISTS variants_chr_pos_idx ON variants (chr, pos);
CREATE INDEX IF NOT EXISTS variants_pos_brin    ON variants USING brin (pos);
CREATE INDEX IF NOT EXISTS variants_sample_idx  ON variants (sample_id);

-- Phenotypes / traits
CREATE TABLE IF NOT EXISTS phenotypes (
    phenotype_id   BIGSERIAL PRIMARY KEY,
    sample_id      TEXT REFERENCES samples(sample_id),
    trait          TEXT NOT NULL,
    value          DOUBLE PRECISION,
    unit           TEXT,
    case_control   BOOLEAN,
    created_at     TIMESTAMPTZ DEFAULT now()
);

CREATE INDEX IF NOT EXISTS phenotypes_trait_idx ON phenotypes (trait);

-- Gene and functional annotations
CREATE TABLE IF NOT EXISTS gene_annotations (
    gene_id        TEXT PRIMARY KEY,
    symbol         TEXT,
    chr            TEXT,
    start_pos      INTEGER,
    end_pos        INTEGER,
    transcript     TEXT,
    consequence    TEXT,
    source         TEXT,
    created_at     TIMESTAMPTZ DEFAULT now()
);

CREATE INDEX IF NOT EXISTS gene_annotations_range_gist
    ON gene_annotations USING gist (chr, int4range(start_pos, end_pos, '[]'));

-- File manifest for traceability
CREATE TABLE IF NOT EXISTS file_manifest (
    file_id        BIGSERIAL PRIMARY KEY,
    path           TEXT UNIQUE NOT NULL,
    md5            CHAR(32),
    size_bytes     BIGINT,
    source         TEXT,
    license        TEXT,
    version        TEXT,
    created_at     TIMESTAMPTZ DEFAULT now()
);

-- Simple run tracking (ETL / analysis)
CREATE TABLE IF NOT EXISTS runs (
    run_id         BIGSERIAL PRIMARY KEY,
    name           TEXT UNIQUE NOT NULL,
    config         JSONB DEFAULT '{}'::JSONB,
    status         TEXT,
    started_at     TIMESTAMPTZ DEFAULT now(),
    finished_at    TIMESTAMPTZ
);
