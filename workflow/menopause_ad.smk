import os

configfile: "config/menopause_ad.yaml"

GSE_IDS = config["gse_ids"]
OUT_BASE = "data/geo"

rule all:
    input:
        expand("results/de/{gse}.de.tsv", gse=GSE_IDS),
        expand("results/signatures/{gse}.up.txt", gse=config["menopause_sets"]),
        "results/ssgsea/ad_signature_scores.tsv"

rule fetch_geo:
    output:
        expr=f"{OUT_BASE}/{{gse}}/expression.tsv.gz",
        pheno=f"{OUT_BASE}/{{gse}}/pheno.tsv"
    params:
        outdir=lambda wildcards: f"{OUT_BASE}/{wildcards.gse}"
    threads: 2
    shell:
        "mkdir -p {params.outdir} && "
        "/Users/yuewenhao/miniforge3/bin/python scripts/py_fetch_geo.py {wildcards.gse} {params.outdir}"

rule de_ttest:
    input:
        expr=f"{OUT_BASE}/{{gse}}/expression.tsv.gz",
        pheno=f"{OUT_BASE}/{{gse}}/pheno.tsv"
    output:
        tsv="results/de/{gse}.de.tsv",
        volcano="results/de/{gse}.volcano.tsv"
    params:
        group=lambda wc: config["designs"][wc.gse]["group_col"],
        case=lambda wc: config["designs"][wc.gse]["case"],
        ctrl=lambda wc: config["designs"][wc.gse]["ctrl"],
        prefix=lambda wc: f"results/de/{wc.gse}"
    threads: 2
    shell:
        "mkdir -p results/de && "
        "/Users/yuewenhao/miniforge3/bin/python scripts/py_de_ttest.py {input.expr} {input.pheno} "
        "\"{params.group}\" \"{params.case}\" \"{params.ctrl}\" {params.prefix}"

rule build_signature:
    input:
        de="results/de/{gse}.de.tsv"
    output:
        up="results/signatures/{gse}.up.txt",
        down="results/signatures/{gse}.down.txt"
    params:
        logfc=lambda wc: config["signature"]["logfc"],
        fdr=lambda wc: config["signature"]["fdr"],
        prefix=lambda wc: f"results/signatures/{wc.gse}"
    threads: 1
    shell:
        "mkdir -p results/signatures && "
        "/Users/yuewenhao/miniforge3/bin/python scripts/py_build_signature.py {input.de} {params.logfc} {params.fdr} {params.prefix}"

rule score_ad:
    input:
        expr=f"{OUT_BASE}/{config['ad_gse']}/expression.tsv.gz",
        up=lambda wc: f"results/signatures/{config['menopause_sets'][0]}.up.txt",
        down=lambda wc: f"results/signatures/{config['menopause_sets'][0]}.down.txt"
    output:
        "results/ssgsea/ad_signature_scores.tsv"
    params:
        out="results/ssgsea/ad_signature_scores.tsv"
    threads: 2
    shell:
        "mkdir -p results/ssgsea && "
        "/Users/yuewenhao/miniforge3/bin/python scripts/py_score_signature.py {input.expr} {input.up} {input.down} {output}"
