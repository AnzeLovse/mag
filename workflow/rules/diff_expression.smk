rule deseq2:
    input:
        counts="counts/merged.tsv",
    output:
        rds="deseq2/{contrast}_dds.rds",
        table="deseq2/{contrast}.diffexp.tsv",
        ma_plot="deseq2/{contrast}.ma-plot.svg",
    params:
        samples=config["samples"],
        units=config["replicates"],
        contrast=get_contrast,
    log:
        "logs/deseq2/{contrast}.diffexp.log",
    threads: config["diffexp"]["threads"]
    script:
        "../scripts/deseq2.R"