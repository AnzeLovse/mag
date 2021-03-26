rule deseq2_preprocess:
    input:
        counts="counts/merged.tsv",
    output:
        "deseq2/dds.rds",
    params:
        samples=config["samples"],
        units=config["replicates"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/preprocess.log",
    script:
        "../scripts/deseq2_preprocess.R"

rule deseq2:
    input:
        "deseq2/dds.rds",
    output:
        table="deseq2/{contrast}.diffexp.tsv",
        ma_plot="deseq2/{contrast}.ma-plot.svg",
    params:
        contrast=get_contrast,
    log:
        "logs/deseq2/{contrast}.diffexp.log",
    conda:
        "../envs/deseq2.yaml"
    threads: config["diffexp"]["threads"]
    script:
        "../scripts/deseq2.R"