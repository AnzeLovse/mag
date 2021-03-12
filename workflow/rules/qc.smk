rule fastqc:
    input:
        get_single_fastq
    output:
        html="qc/fastqc/{sample}-{rep}_R{pair}.html",
        zip="qc/fastqc/{sample}-{rep}_R{pair}_fastqc.zip"
    log:
        "logs/fastqc/{sample}-{rep}_R{pair}.log"
    threads: 1
    wrapper:
        "0.72.0/bio/fastqc" # TODO update version to correctly remove .fq.gz ending

rule fastqc_post_trim:    
    input:
        "trimmed/{sample}-{rep}_{pair}.fastq.gz",
    output:
        html="qc/fastqc_posttrim/{sample}-{rep}_R{pair}.html",
        zip="qc/fastqc_posttrim/{sample}-{rep}_R{pair}_fastqc.zip"
    log:
        "logs/fastqc_posttrim/{sample}-{rep}_R{pair}.log"
    threads: 1
    wrapper:
        "0.72.0/bio/fastqc" # TODO update version to correctly remove .fq.gz ending

rule samtools_stats:
    input:
        bam="sorted_reads/{sample}-{rep}.bam",
        idx="sorted_reads/{sample}-{rep}.bam.bai"
    output:
        "qc/samtools_stats/{sample}-{rep}.txt"
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools_stats/{sample}-{rep}.log"
    shell:
        "samtools stats {input.bam} 1> {output} 2> {log}"
