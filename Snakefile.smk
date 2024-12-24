# List of samples and their NCBI accessions
samples = {
    "sample1": "GCF_000009865.1",
    "sample2": "GCF_000195755.2"
}

rule all:
    input:
        expand("data/results/{sample}/antismash.json", sample=samples.keys())

rule fetch_genomes:
    output:
        "data/raw/{sample}.fasta"
    params:
        accession=lambda wildcards: samples[wildcards.sample]
    script:
        "scripts/fetch_genomes.py"

rule run_antismash:
    input:
        "data/raw/{sample}.fasta"
    output:
        "data/results/{sample}/antismash.json"
    params:
        outdir="data/results/{sample}"
    conda:
        "envs/bgc_env.yaml"
    shell:
        """
        antismash {input} --output-dir {params.outdir} --json
        """

rule parse_results:
    input:
        "data/results/{sample}/antismash.json"
    output:
        "data/processed/{sample}_parsed.tsv"
    script:
        "scripts/parse_antismash.py"
