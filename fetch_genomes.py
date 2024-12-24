import os
from Bio import Entrez

Entrez.email = "your_email@example.com"

def fetch_genome(accession, output_path):
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    with open(output_path, "w") as out_file:
        out_file.write(handle.read())

# Access parameters directly from Snakemake
accession = snakemake.params.accession
output_path = snakemake.output[0]

fetch_genome(accession, output_path)

