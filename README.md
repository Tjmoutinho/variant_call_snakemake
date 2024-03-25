# variant_call_snakemake
Snakemake pipeline to call variants using illumina fastq data

# Create the docker image:
docker build vcf_sm:latest

# Execute using a command similar to this:
docker run -v /path/to/repo/variant_call_snakemake:/usr/src/vcf vcf:latest snakemake --cores all

You will need to replace '/path/to/repo/' with your actual path
This command will create a directory called 'results' with all of the output files.

# Run pytest using:
docker run -v /path/to/repo/variant_call_snakemake:/usr/src/vcf vcf:latest pytest test.py

# Results:
results_final/ Contains the results of this pipeline when run with the data and config file provided. 

# General results for each gene by name:
Summary of the results found in /path/to/repo/vcf_nextflow_test/results_final/reports/sample1_snpEff_report.genes.txt
E: 0 MODERATE impact
M: 1 MODERATE impact, 1 missense variant
N: 4 MODERATE impact, 4 missense variants
ORF10: 0 MODERATE impact
ORF3a: 2 MODERATE impact, 2 missense variants
ORF6: 0 MODERATE impact
ORF7a: 2 MODERATE impact, 2 missense variants
ORF8: 1 MODERATE impact
S: 7 MODERATE impact, 6 missense variants, 1 disruptive inframe deletion
orf1ab: 10 MODERATE impact, 10 missense variants

# Reports:
There are several reports available. The MultiQC report is the most informative for sequencing QC and some basic information about the VCF results. There is also a snpEff Report that has more detailed information about the VCF. Generally speaking, the fastq files look good enough for reliable results. 