rule fastqc:
    input:
        get_single_fastq
    output:
        html="qc/fastqc/{sample}-t{time}-{rep}_R{pair}.html",
        zip="qc/fastqc/{sample}-t{time}-{rep}_R{pair}_fastqc.zip"
    log:
        "logs/fastqc/{sample}-t{time}-{rep}_R{pair}.log"
    threads: config["fastqc"]["threads"]
    wrapper:
        "0.72.0/bio/fastqc" # TODO update version to correctly remove .fq.gz ending

rule fastqc_post_trim:
    input:
        "trimmed/{sample}-t{time}-{rep}_{pair}.fastq.gz",
    output:
        html="qc/fastqc_posttrim/{sample}-t{time}-{rep}_R{pair}.html",
        zip="qc/fastqc_posttrim/{sample}-t{time}-{rep}_R{pair}_fastqc.zip"
    log:
        "logs/fastqc_posttrim/{sample}-t{time}-{rep}_R{pair}.log"
    threads: config["fastqc"]["threads"]
    wrapper:
        "0.72.0/bio/fastqc" # TODO update version to correctly remove .fq.gz ending

rule samtools_stats:
    input:
        bam="sorted_reads/{sample}-t{time}-{rep}.bam",
        idx="sorted_reads/{sample}-t{time}-{rep}.bam.bai"
    output:
        "qc/samtools_stats/{sample}-t{time}-{rep}.txt"
    log:
        "logs/samtools_stats/{sample}-t{time}-{rep}.log"
    shell:
        "samtools stats {input.bam} 1> {output} 2> {log}"

rule multiqc:
    input:
        get_fastqc_outputs,
        get_qc_ouptuts,
        # expand("qc/samtools_stats/{sample}-t{time}-{rep}.txt", sample=replicates["sample"], rep=replicates["replicate"], time=replicates["time"]),
        # expand("logs/bowtie2/{sample}-t{time}-{rep}.log", sample=replicates["sample"], rep=replicates["replicate"], time=replicates["time"]),
        # expand("qc/feature_counts/{sample}-t{time}-{rep}.summary", sample=replicates["sample"], rep=replicates["replicate"], time=replicates["time"]),
        # expand("trimmed/{sample}-t{time}-{rep}.qc.txt", sample=replicates["sample"], rep=replicates["replicate"], time=replicates["time"])
    output:
        "qc/multiqc.html"
    params:
        "--dirs"
    log:
        "logs/multiqc.log"
    wrapper:
        "0.72.0/bio/multiqc"