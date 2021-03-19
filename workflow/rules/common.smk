def get_single_fastq(wildcards):
    return replicates.loc[(wildcards.sample, wildcards.rep), [f"fq{wildcards.pair}"]].dropna()

def get_raw_fastq(wildcards):
    return replicates.loc[(wildcards.sample, wildcards.rep), ["fq1", "fq2"]].dropna()

def get_fastq(wildcards):
    if config["trimming"]["skip"]:
        return replicates.loc[(wildcards.sample, wildcards.rep), ["fq1", "fq2"]].dropna()
    else:
        if is_paired:
            return expand(
                    "trimmed/{sample}-{rep}_{pair}.fastq.gz", pair=[1, 2], **wildcards
                )
        else:
            return ["trimmed/{sample}-{rep}_1.fastq.gz"]

def get_fastqc_outputs(wildcards):
    if config["trimming"]["skip"]:
        if is_paired:
            return expand(
                    "qc/fastqc/{sample}-{rep}_R{pair}_fastqc.zip", pair=[1, 2], **wildcards
                )
        else:
            return "qc/fastqc/{sample}-{rep}_R1_fastqc.zip"
    else:
        if is_paired:
            return expand(
                    "qc/{step}/{sample}-{rep}_R{pair}_fastqc.zip", step=["fastqc", "fastqc_posttrim"],
                    pair=[1, 2], **wildcards
                )
        else:
            return expand(
                    "qc/{step}/{sample}-{rep}_R1_fastqc.zip", step=["fastqc", "fastqc_posttrim"],
                     **wildcards
                )


def get_trim_params(wildcards):
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
    params = config["feature_counts"]["additional"]
    if is_paired:
        params += " -p"
    return params
