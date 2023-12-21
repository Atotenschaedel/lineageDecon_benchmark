# -*- coding: utf-8 -*-

# snakemake workflow
# generate simulated wastewater sequencing data for testing
# snakemake -s simulation.smk -c1 --use-conda 
# 1. call gernateAbundance script --> can't be called in same run
# 2. call swampy tool

# get parameters from config file
import os, yaml
with open("workflow_config.yaml", "r") as cfile:
    config = yaml.safe_load(cfile)

#get parameters from config file
EXPERIMENT_NAME = config["EXPERIMENT_NAME"]

# get newly generated file names
def get_abundances_file_names(wildcards):
    check_output=checkpoints.generateAbundances.get(**wildcards).output[0]
    s_name, = glob_wildcards(os.path.join(check_output, "{sample}.tsv"))
    return expand(os.path.join(check_output, "{SAMPLE}.tsv"), SAMPLE=s_name)

def get_sequences_file_names(wildcards):
    seq_output=checkpoints.generateAbundances.get(**wildcards).output[0]
    s_name, = glob_wildcards(os.path.join(seq_output, "{sample}.tsv"))
    return expand(os.path.join("experiments/"+EXPERIMENT_NAME+"/data", "{SM}_{ext}"), SM=s_name, ext=["R1.fastq","R2.fastq"])


def get_genomes():
    genome_list = ["reference/pango-sequences/data/pango-consensus-sequences_genome-nuc.fasta",
                    "reference/consensus_sequences/others.fasta",
                    "reference/consensus_sequences/B.1.1.529.fasta"]

    return genome_list


## check qc for simulated files - adapter sequences ect
rule all:
    input:
        expand("experiments/{exp}/simulation/QualityControl/{exp}_multiqc_report.html", exp=EXPERIMENT_NAME)

checkpoint generateAbundances:
    output: 
        dir=directory(expand("experiments/{exp}/simulation/abundances/", exp=EXPERIMENT_NAME))
    conda: 
        "envs/python3.yaml"
    params:
        exp=EXPERIMENT_NAME,
        config_file="workflow_config.yaml"
    log:
        expand("logs/{exp}/simulation/genSamples.out.log", exp=EXPERIMENT_NAME)
    shell: 
        "python3 bin/custom_scripts/generateSamples.py --output_dir {output.dir} "
        "--parameter_file {params.config_file} --experiment_name {params.exp} &>{log}"

rule generate_cons_genome:
    input:
        genomes = get_genomes()
    output:
        temp("experiments/{exp}/simulation/genomes.fasta")
    shell:
        # conc all necessary genomes. make sure that reads do not contain gaps
        "cat {input}| sed  '/^>/!s/-//g' >{output} "

rule init_swampy:
    input:
        tsv= get_abundances_file_names,
        genome=expand("experiments/{exp}/simulation/genomes.fasta", exp=EXPERIMENT_NAME)
    output: 
        r1="experiments/{exp}/data/{exp}_{sample}_R1.fastq",
        r2="experiments/{exp}/data/{exp}_{sample}_R2.fastq",
        tsv="experiments/{exp}/simulation/abundances/{exp}_{sample}_amplicon_abundances_summary.tsv",
        vcf=temp("experiments/{exp}/data/{exp}_{sample}_PCR_errors.vcf"),
        l=temp("experiments/{exp}/data/{exp}_{sample}.log")
    conda: 
        "envs/swampy.yaml"
    params:
        exp=EXPERIMENT_NAME,
        primer_set=config["SIMULATION"]["PRIMER_SET"],
        n_reads=config["SIMULATION"]["N_READS"],
        seqSys=config["SIMULATION"]["SEQ_SYS"],
        read_length=config["SIMULATION"]["READ_LENGTH"],
        amplicon_distribution=config["SIMULATION"]["AMPLICON_DISTRIBUTION"],
        amplicon_pseudocounts=config["SIMULATION"]["AMPLICON_PSEUDOCOUNTS"]
    log:    
        "logs/{exp}/simulation/SWAMPy_{sample}.out.log"
    shell:
        "python bin/SWAMPy/src/simulate_metagenome.py --genomes_file {input.genome} "
        "--genome_abundances experiments/{wildcards.exp}/simulation/abundances/{wildcards.exp}_{wildcards.sample}.tsv " 
        "--output_folder experiments/{params.exp}/data --output_filename_prefix {wildcards.exp}_{wildcards.sample} "
        "--n_reads {params.n_reads} --seqSys {params.seqSys} --read_length {params.read_length} --primer_set {params.primer_set} "
        "--amplicon_distribution {params.amplicon_distribution} --amplicon_pseudocounts {params.amplicon_pseudocounts} "
        "--temp_folder bin/SWAMPy/temp/{wildcards.sample} --autoremove &>{log} "
        "&& mv experiments/{params.exp}/data/{params.exp}_{wildcards.sample}_amplicon_abundances_summary.tsv experiments/{params.exp}/simulation/abundances"

rule qualityControl:
    input:
        get_sequences_file_names
    output:
        "experiments/{exp}/simulation/QualityControl/{exp}_multiqc_report.html"
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc {input} --quiet -o experiments/{wildcards.exp}/simulation/QualityControl && "
        "multiqc experiments/{wildcards.exp}/simulation/QualityControl/. --outdir experiments/{wildcards.exp}/simulation/QualityControl/ "
        "--quiet --no-data-dir --filename {wildcards.exp}_multiqc_report.html"
