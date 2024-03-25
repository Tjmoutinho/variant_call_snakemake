configfile: "config.yaml"

def get_fastq_paths(sample, pair):
    return config['samples'][sample][pair]

rule all:
    input:
        expand("results/marked_reads/{sample}.bam", sample=config["samples"]),
        expand("results/marked_reads/{sample}.bam.bai", sample=config["samples"]),
        expand("results/calls/{sample}.vcf", sample=config["samples"]),
        expand("results/calls/hq_{sample}.vcf", sample=config["samples"]),
        expand("results/annotated_calls/hq_{sample}.annotated.vcf", sample=config["samples"]),
        expand("results/reports/{sample}_snpEff_report.csv", sample=config["samples"]),
        "results/multiqc_report.html"

rule unzip_fastq:
    input:
        lambda wildcards: get_fastq_paths(wildcards.sample, wildcards.pair)
    output:
        temp("data/{sample}_{pair}.fastq")
    threads: 8
    shell:
        "gunzip -c {input} > {output}"

rule fastqc_before:
    input:
        fastq="data/{sample}_{pair}.fastq"
    output:
        html="results/qc/{sample}_{pair}_before_fastqc.html",
        zip="results/qc/{sample}_{pair}_before_fastqc.zip"
    threads: 8
    params:
        outdir="results/qc/"
    shell:
        """
        fastqc {input.fastq} --outdir {params.outdir}
        mv {params.outdir}{wildcards.sample}_{wildcards.pair}_fastqc.html {output.html}
        mv {params.outdir}{wildcards.sample}_{wildcards.pair}_fastqc.zip {output.zip}
        """

rule trimmomatic:
    input:
        fq1="data/{sample}_R1.fastq",
        fq2="data/{sample}_R2.fastq"
    output:
        trimmed_fq1="results/trimmed_data/{sample}_R1.fastq",
        trimmed_fq2="results/trimmed_data/{sample}_R2.fastq"
    threads: 8
    shell:
        "trimmomatic PE {input.fq1} {input.fq2} {output.trimmed_fq1} /dev/null {output.trimmed_fq2} /dev/null SLIDINGWINDOW:4:20 MINLEN:36"

rule downsample:
    input:
        fq1="results/trimmed_data/{sample}_R1.fastq",
        fq2="results/trimmed_data/{sample}_R2.fastq"
    output:
        downsampled_fq1="results/downsampled_data/{sample}_R1.fastq",
        downsampled_fq2="results/downsampled_data/{sample}_R2.fastq"
    threads: 8
    shell:
        """
        seqtk sample -s100 {input.fq1} 16000 > {output.downsampled_fq1}
        seqtk sample -s100 {input.fq2} 16000 > {output.downsampled_fq2}
        """

rule fastqc_after:
    input:
        fastq="results/downsampled_data/{sample}_{pair}.fastq"
    output:
        html="results/qc/{sample}_{pair}_after_fastqc.html",
        zip="results/qc/{sample}_{pair}_after_fastqc.zip"
    threads: 8
    params:
        outdir="results/qc/"
    shell:
        """
        fastqc {input.fastq} --outdir {params.outdir}
        mv {params.outdir}/{wildcards.sample}_{wildcards.pair}_fastqc.html {output.html}
        mv {params.outdir}/{wildcards.sample}_{wildcards.pair}_fastqc.zip {output.zip}
        """

rule bwa_index:
    input:
        fa=config["reference"]
    output:
        amb=temp(config["reference"] + ".amb"),
        ann=temp(config["reference"] + ".ann"),
        bwt=temp(config["reference"] + ".bwt"),
        pac=temp(config["reference"] + ".pac"),
        sa=temp(config["reference"] + ".sa")
    threads: 8
    shell:
        "bwa index {input.fa}"

rule faidx:
    input:
        fa=config["reference"]
    output:
        fai=temp(config["reference"] + ".fai")
    threads: 8
    shell:
        "samtools faidx {input} > {output}"

rule bwa:
    input:
        fa=config["reference"],
        fq1="results/downsampled_data/{sample}_R1.fastq",
        fq2="results/downsampled_data/{sample}_R2.fastq",
        idx=lambda wc: expand("{reference}.{ext}", reference=config["reference"], ext=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        temp("results/mapped_reads/{sample}.bam")
    params:
        rg="@RG\\tID:{sample}\\tSM:{sample}"
    log:
        "results/logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input.fa} {input.fq1} {input.fq2} | samtools view -Sb - > {output}) 2> {log}" 

rule sambamba_sort:
    input:
        bam="results/mapped_reads/{sample}.bam"
    output:
        sorted_bam=temp("results/mapped_reads/{sample}.sorted.bam")
    threads: 8
    shell:
        """
        sambamba sort -t {threads} -m 5GB --tmpdir /tmp/ -o {output.sorted_bam} {input.bam}
        """

rule sambamba_markdup:
    input:
        sorted_bam="results/mapped_reads/{sample}.sorted.bam"
    output:
        marked_bam="results/marked_reads/{sample}.bam"
    threads: 8
    shell:
        """
        sambamba markdup -t {threads} {input.sorted_bam} {output.marked_bam}
        """

rule sambamba_index:
    input:
        "results/marked_reads/{sample}.bam"
    output:
        "results/marked_reads/{sample}.bam.bai"
    threads: 8
    shell:
        """
        sambamba index -t {threads} {input} {output}
        """

rule freebayes_call:
    input:
        fa=config["reference"],
        fai=expand("{reference}.fai", reference=config["reference"]),
        bam=expand("results/marked_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("results/marked_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        vcf=expand("results/calls/{sample}.vcf", sample=config["samples"])
    threads: 8
    log:
        expand("results/logs/free_bayes/{sample}.log", sample=config["samples"])
    shell:
        "(freebayes -f {input.fa} --ploidy 1 {input.bam} > {output.vcf}) 2> {log}"
        
rule vcffilter:
    input:
        vcf=expand("results/calls/{sample}.vcf", sample=config["samples"])
    output:
        hq_vcf=expand("results/calls/hq_{sample}.vcf", sample=config["samples"])
    threads: 8
    params:
        qs = config["vcf_quality_threshold"]
    shell:
        "vcffilter -f 'QUAL > {params.qs}' {input.vcf} > {output}"

rule snpeff_annotation:
    input:
        vcf="results/calls/hq_{sample}.vcf"
    output:
        annotated_vcf="results/annotated_calls/hq_{sample}.annotated.vcf",
        csv_report="results/reports/{sample}_snpEff_report.csv",
        txt_report="results/reports/{sample}_snpEff_report.genes.txt",
        html_report="results/reports/{sample}_snpEff_report.html"
    threads: 8
    params:
        snpeff_db= "MN908947.3",
    log:
        "results/logs/snpeff/hq_{sample}.log"
    shell:
        """
        snpEff download -v {params.snpeff_db}
        snpEff ann -csvStats {output.csv_report} -Stats {output.html_report} {params.snpeff_db} {input.vcf} > {output.annotated_vcf} 2> {log}
        """

rule multiqc:
    input:
        expand("results/qc/{sample}_{pair}_after_fastqc.html", sample=config["samples"], pair=['R1', 'R2']),
        snpeff_report=expand("results/reports/{sample}_snpEff_report.csv", sample=config["samples"])
    output:
        "results/multiqc_report.html"
    threads: 8
    shell:
        "multiqc results/qc/ {input.snpeff_report} --filename {output} --verbose"