if is_paired:
    rule cutadapt:
        input:
            get_raw_fastq
        output:
            fastq1="trimmed/{sample}-{rep}_1.fastq.gz",
            fastq2="trimmed/{sample}-{rep}_2.fastq.gz",
            qc="trimmed/{sample}-{rep}.qc.txt",
        params:
            adapters="",
            extra=get_trim_params
        log:
            "logs/cutadapt/{sample}-{rep}.log"
        threads: 4
        wrapper:
            "0.72.0/bio/cutadapt/pe"
else:
    rule cutadapt:
        input:
            get_raw_fastq
        output:
            fastq="trimmed/{sample}-{rep}_1.fastq.gz",
            qc="trimmed/{sample}-{rep}.qc.txt",
        params:
            adapters="",
            extra=get_trim_params
        log:
            "logs/cutadapt/{sample}-{rep}.log"
        threads: 4
        wrapper:
            "0.72.0/bio/cutadapt/se"