def get_single_fastq(wildcards):
    """Get file name of reads file."""
    return replicates.loc[(wildcards.sample, wildcards.rep, wildcards.time), [f"fq{wildcards.pair}"]].dropna()

def get_raw_fastq(wildcards):
    """Get file names of raw reads."""
    return replicates.loc[(wildcards.sample, wildcards.rep, wildcards.time), ["fq1", "fq2"]].dropna()

def get_fastq(wildcards):
    """Get file names of reads for FASTQC."""
    if config["trimming"]["skip"]:
        return replicates.loc[(wildcards.sample, wildcards.rep, wildcards.time), ["fq1", "fq2"]].dropna()
    else:
        if is_paired:
            return expand(
                    "trimmed/{sample}-t{time}-{rep}_{pair}.fastq.gz", pair=[1, 2], **wildcards
                )
        else:
            return ["trimmed/{sample}-t{time}-{rep}_1.fastq.gz"]

def get_fastqc_outputs(wildcards):
    """Get all FASTQC outputs for a sample."""
    if config["trimming"]["skip"]:
        if is_paired:
            return expand(
                    "qc/fastqc/{sample}-t{time}-{rep}_R{pair}_fastqc.zip", pair=[1, 2],
                    sample=replicates["sample"], rep=replicates["replicate"], time=replicates["time"],
                )
        else:
            return expand(
                "qc/fastqc/{sample}-t{time}-{rep}_R1_fastqc.zip",
                sample=replicates["sample"], rep=replicates["replicate"], time=replicates["time"],)
    else:
        if is_paired:
            samples = [
                f"qc/{{step}}/{sample}-t{time}-{rep}_R{{pair}}_fastqc.zip"
                for sample, rep, time in replicates.index.to_list()
            ]
            return expand(
                    samples, step=["fastqc", "fastqc_posttrim"], pair=[1, 2]
                )
        else:
            return expand(
                    "qc/{step}/{sample}-t{time}-{rep}_R1_fastqc.zip", step=["fastqc", "fastqc_posttrim"],
                     sample=replicates["sample"], rep=replicates["replicate"], time=replicates["time"],
                )

def get_qc_ouptuts(wildcards):
    """Get all QC outputs for a sample."""
    samtools = [
        f"qc/samtools_stats/{sample}-t{time}-{rep}.txt"
        for sample, rep, time in replicates.index.to_list()
    ]
    bowtie = [
        f"logs/bowtie2/{sample}-t{time}-{rep}.log"
        for sample, rep, time in replicates.index.to_list()
    ]
    feature_counts = [
        f"qc/feature_counts/{sample}-t{time}-{rep}.summary"
        for sample, rep, time in replicates.index.to_list()
    ]
    trimming = [
        f"trimmed/{sample}-t{time}-{rep}.qc.txt"
        for sample, rep, time in replicates.index.to_list()
    ]
    return samtools + bowtie + feature_counts + trimming

def get_visualisation_ouptuts(wildcards):
    """Get all visualisation outputs for a sample."""
    return [
        f"visualisations/{sample}-t{time}-{rep}_coverage.json"
        for sample, rep, time in replicates.index.to_list()
    ]

def get_trim_params(wildcards):
    """Get a string with Cutadapt parameters."""
    params=[
        "--quality-cutoff", config["trimming"]["min_qual"],
        "--minimum-length", config["trimming"]["min_len"]]
    if config["trimming"]["adapter_1"]:
        params.extend(["-a", config["trimming"]["adapter_1"]])
    if config["trimming"]["adapter_2"]:
        params.extend(["-A", config["trimming"]["adapter_2"]])
    if config["trimming"]["additional"]:
        params.append(config["trimming"]["additional"])

    return " ".join(params)

def feature_counts_params(wildcards):
    """Generate a string with featureCounts parameters."""
    params = config["feature_counts"]["additional"]
    if is_paired:
        params += " -p"
    return params
