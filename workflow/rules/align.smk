rule bowtie2:
    input:
        sample=get_fastq,
        index=rules.bowtie2_index.output,
    output:
        "mapped_reads/{sample}-{rep}.bam",
    log:
        "logs/bowtie2/{sample}-{rep}.log",
    params:
        index="index/genome",
        extra="",
    threads: 8
    wrapper:
        "0.72.0/bio/bowtie2/align"


rule samtools_sort:
    input:
        "mapped_reads/{sample}-{rep}.bam"
    output:
        "sorted_reads/{sample}-{rep}.bam"
    params:
        threads=config["samtools"]["threads"]
    shell:
        "samtools sort -@ {params.threads} -T sorted_reads/{wildcards.sample}-{wildcards.rep} -O bam -o {output}  {input}"


rule samtools_index:
    input:
        "sorted_reads/{sample}-{rep}.bam"
    output:
        "sorted_reads/{sample}-{rep}.bam.bai"
    params:
        threads=config["samtools"]["threads"] - 1
    shell:
        "samtools index -@ {params.threads} {input}"