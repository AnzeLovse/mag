rule deseq2:
    input:
        counts="counts/merged.tsv",
    output:
        rds="deseq2/dds.rds",
        table="deseq2/diffexp.tsv",
        pca_plot="deseq2/pca_plot.svg",
    params:
        samples=config["samples"],
        units=config["replicates"],
        model=config["diffexp"]["model"],
        reduced_model=config["diffexp"]["reduced_model"],
        alpha=config["diffexp"]["alpha"],
    log:
        "logs/deseq2/diffexp.log",
    script:
        "../scripts/deseq2.R"