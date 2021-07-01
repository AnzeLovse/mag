def prepare_annotation(annotation, strand, color):
    annotation = annotation.loc[annotation["strand"] == strand]
    annotation = annotation.drop(columns=["type", "strand"])
    annotation["color"] = color

    return json.loads(annotation.to_json(orient="records"))


rule gtf_to_bed:
    input:
        annotation=config["feature_counts"]["annotation"]
    output:
        complete="visualisations/gene_annotation.bed"
    threads: 1
    shell:
        "cat {input.annotation} | "
        "awk 'OFS=\"\\t\" {{if ($3==\"gene\") {{print $1,$4-1,$5,$10,$16,$7}}}}' | "
        "tr -d '\";' > {output.complete}"


rule create_annotation:
    input:
        annotation="visualisations/gene_annotation.bed"
    output:
        json="visualisations/annotation.json",
    threads: 1
    run:
        annotation = pd.read_csv(
            input.annotation,
            sep="\t",
            header=None,
            names=["chr", "start", "end", "name", "type", "strand"]
        )
        annotation = annotation.rename(columns = {"chr":"block_id"})

        pos = prepare_annotation(annotation, "+", "rgb(255,51,51)")
        neg = prepare_annotation(annotation, "-", "rgb(255,51,51)")
        joined = {"positive": pos, "negative": neg}
        with open(output.json, "w") as handle:
            json.dump(joined, handle, indent=4)


rule calculate_coverage:
    input:
        bam="sorted_reads/{sample}-t{time}-{rep}.bam"
    output:
        coverage = "bedgraphs/{sample}-t{time}-{rep}_coverage_{strand}.bedgraph"
    threads: 1
    params:
        paired = "-pc" if is_paired else ""
    log:
        "logs/bedgraphs/{sample}-t{time}-{rep}_coverage_{strand}.log"
    shell:
        "genomeCoverageBed "
        "-ibam {input.bam} "
        "{params.paired} "
        "-bga "
        "-strand {wildcards.strand} "
        "> {output.coverage}"


rule index_fasta:
    input:
        config["genome"]
    output:
        f"{config['genome']}.fai"
    threads: 1
    shell:
        "samtools faidx {input}"


rule coverage_to_json:
    input:
        coverage_pos = "bedgraphs/{sample}-t{time}-{rep}_coverage_+.bedgraph",
        coverage_neg = "bedgraphs/{sample}-t{time}-{rep}_coverage_-.bedgraph",
        fai = f"{config['genome']}.fai"
    output:
        json="visualisations/{sample}-t{time}-{rep}_coverage.json"
    run:
        fai = pd.read_csv(
            input.fai,
            sep="\t",
            header=None,
            names=["chr", "len"],
            usecols=[0,1],
            index_col=0
        )
        lengths = dict(zip(fai.index, fai.len))

        coverage = {"length":lengths, "positive":{}, "negative":{}}
        coverage_pos = pd.read_csv(input.coverage_pos, sep="\t", header=None, names=["block_id", "start", "end", "value"])
        coverage_neg = pd.read_csv(input.coverage_neg, sep="\t", header=None, names=["block_id", "start", "end", "value"])

        for chrom in lengths.keys():
            pos = coverage_pos.loc[coverage_pos["block_id"] == chrom]
            neg = coverage_neg.loc[coverage_neg["block_id"] == chrom]

            coverage["positive"][chrom] = json.loads(pos.to_json(orient="records"))
            coverage["negative"][chrom] = json.loads(neg.to_json(orient="records"))

        with open(output.json, "w") as handle:
            json.dump(coverage, handle, indent=4)
