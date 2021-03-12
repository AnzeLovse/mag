rule bowtie2_index:
    input:
        reference=config["genome"]
    output:
        multiext(
            "index/genome",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    log:
        "logs/bowtie2_build/build.log"
    params:
        extra=""
    threads: 8
    wrapper:
        "0.72.0-5-gd0d054c/bio/bowtie2/build"