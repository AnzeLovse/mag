rule feature_counts:
    input:
        bam="sorted_reads/{sample}-{rep}.bam",
        bai="sorted_reads/{sample}-{rep}.bam.bai",
    output:
        counts="counts/{sample}-{rep}.txt",
        summary="qc/feature_counts/{sample}-{rep}.txt"
    params:
        annotation=config["feature_counts"]["annotation"],
        extra=feature_counts_params
    threads:
        config["feature_counts"]["threads"]
    log:
        "logs/feature_counts/{sample}-{rep}.txt"
    shell:
        "featureCounts "
        "{params.extra} "
        "-a {params.annotation} "
        "-o {output.counts} "
        "-T {threads} "
        "{input.bam} > {log} 2>&1; "
        "mv \"{output.counts}.summary\" {output.summary}"


rule merge_counts:
    input:
        [f"counts/{sample}-{rep}.txt" for sample, rep in zip(replicates["sample"], replicates["replicate"])]

    output:
        "counts/merged.tsv",
        "counts/merged_raw.tsv"
    run:
        indices = ["Geneid", "Chr", "Start", "End", "Strand","Length"]
        frames = [
            pd.read_csv(fp, sep="\t", skiprows=1, index_col=indices)
            for fp in input
        ]
        merged = pd.concat(frames, axis=1)

        merged = merged.rename(columns=lambda c: Path(c).stem)
        merged.to_csv(output[0], sep="\t", index=True)

        # Save counts compatible with rnanorm
        merged = merged.reset_index(level="Geneid")
        merged = merged.rename(columns={"Geneid":"FEATURE_ID"})
        merged.to_csv(output[1], sep="\t", index=False)
rule normalize_counts:
    input:
        counts="counts/merged_raw.tsv",
        annotation=config["feature_counts"]["annotation"]
    output:
        cpm="counts/merged_cpm.tsv",
        tpm="counts/merged_tpm.tsv"
    shell:
        "rnanorm {input.counts} --cpm-output={output.cpm} --tpm-output={output.tpm} --annotation={input.annotation}"
