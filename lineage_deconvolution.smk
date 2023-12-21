# snakemake workflow
"""
Created on Tue Jun 27 11:25:22 2023
This script is doing ....
@author: aschedl
"""

# get parameters from config file
import yaml
from pathlib import Path
with open("workflow_config.yaml", "r") as cfile:
    config = yaml.safe_load(cfile)

# get parameters from config file
EXPERIMENT_NAME = config["EXPERIMENT_NAME"]
REF_GENOME_COMB = config["REFERENCES"]["REF_GENOME_COMB"][0]
REF_GENOME = config["REFERENCES"]["REF_GENOME"][0]
PRIMER_SET = config["REFERENCES"]["PRIMER_SET"][0]

if config["APP"]["DEBUG"]:
    SAMPLES = config["APP"]["DEBUG_SAMPLES"]
else:
    p = Path.cwd() / "experiments" / EXPERIMENT_NAME / "data"
    SAMPLES = list(set([x.stem[:-3] for x in p.glob("*") if x.is_file()]))

# get sample metadata
def get_sample_date(EXPERIMENT_NAME, sample):
    with open("experiments/"+ EXPERIMENT_NAME + "/simulation/" + EXPERIMENT_NAME + "_metadata.tsv") as f:
        lines = f.readlines()
        for l in lines:
            if l.split("\t")[0] == sample:
                return l.split("\t")[1].strip()

rule all:
    input: 
        expand("experiments/{exp}/results/postPrediction/data.csv", exp=EXPERIMENT_NAME)


rule postPrediction:
    input:
        expand("experiments/{exp}/results/vaquero/globalFullData.csv.gz", exp=EXPERIMENT_NAME),
        expand("experiments/{exp}/results/freyja/{sample}.tsv", exp=EXPERIMENT_NAME, sample=SAMPLES),
        expand("experiments/{exp}/results/VLQ/{sample}_predictions.tsv", exp=EXPERIMENT_NAME, sample=SAMPLES),
        expand("experiments/{exp}/results/lollipop/deconvoluted.tsv", exp=EXPERIMENT_NAME)
    output:
        expand("experiments/{exp}/results/postPrediction/data.csv", exp=EXPERIMENT_NAME)
    params:
        exp = EXPERIMENT_NAME
    conda:
        "envs/python3.yaml"
    log:
        expand("logs/{exp}/postPrediction/postPrediction.out", exp=EXPERIMENT_NAME)
    shell:
        "python bin/custom_scripts/PostPred_dataFilter.py -e {params.exp} >>{log} 2>&1"

rule vaquero_generateInput_af: #generate tidy table from vcf files 
    input: "experiments/{exp}/results/variantCall/12_phasedLofreqVcf/{sample}_samp.lofreq_filtered_phased.vcf.gz"
    output: temp("experiments/{exp}/results/vaquero/{sample}_af.tsv")
    params: 
        sample_suffix = config["VAQUERO"]["SAMPLE_SUFFIX"]
    conda: 
        "envs/pysam.yaml"
    shell:
        "python3 bin/VaQuERo/scripts/vcf2tsv_long.py -i {input} -m 0.1 | "
        "sed 's/{wildcards.sample}/{wildcards.sample}{params.sample_suffix}/g' >{output} "

rule vaquero_merge_afs: # merge tidy tables togehter, remove all but first header
    input: expand("experiments/{exp}/results/vaquero/{sample}_af.tsv", exp=EXPERIMENT_NAME, sample=SAMPLES)
    output: expand("experiments/{exp}/results/vaquero/Input_afs.tsv", exp=EXPERIMENT_NAME)
    shell: "cat {input} | sed -e '2,${{ /^SAMPLEID/d }}' >{output}"

rule vaquero_generateInput_metaData: # generate metadata for vaquero based on generic metadata file
    input: "experiments/{exp}/simulation/{exp}_metadata.tsv"
    output: "experiments/{exp}/results/vaquero/Input_metadata.tsv"
    conda: "envs/python3.yaml"
    shell: "python3 bin/custom_scripts/generateToolsMetadata.py -i {input} -t vaquero -o {output}"

rule vaquero_main:
    input:
        afs="experiments/{exp}/results/vaquero/Input_afs.tsv",
        meta="experiments/{exp}/results/vaquero/Input_metadata.tsv"
    output: "experiments/{exp}/results/vaquero/globalFullData.csv.gz"
    params:
        markerDate = config["VAQUERO"]["MARKER_DATE"],
        markerLocation = config["VAQUERO"]["MARKER_LOCATION"]
    conda: 
        "envs/vaquero.yaml"
    log: 
        "logs/{exp}/vaquero/vaquero.out"
    shell:
        "Rscript bin/custom_scripts/install_devtools.R &>{log} && "
        "Rscript bin/VaQuERo/scripts/VaQuERo_v2.r --metadata {input.meta} --data {input.afs} --inputformat tidy "
        "--marker bin/VaQuERo/resources/mutations_list_grouped_pango_codonPhased_{params.markerDate}_{params.markerLocation}.csv "
        "--smarker bin/VaQuERo/resources/mutations_special_2022-12-21.csv "
        "--pmarker bin/VaQuERo/resources/mutations_problematic_vss1_v3.csv "
        "--dir experiments/{wildcards.exp}/results/vaquero/ "
        "--zero=0.01 --depth=50 --minuniqmark=1 --minuniqmarkfrac=0.5 "
        "--minqmark=5 --minmarkfrac=0.5 --smoothingsamples=2 --smoothingtime=22 "
        "--recent 9999 --voi=XBB.1.16,XBB.2.3,FE.1 --colorBase=XBB,BA.1,BA.2,BA.5 >>{log} 2>&1" 

rule freyja_generateInput:
    input:
        gz="experiments/{exp}/results/variantCall/11_lofreqrefilt/{sample}_samp.lofreq_filtered.vcf.gz",
        bam="experiments/{exp}/results/variantCall/07_clipOverlap/{sample}_ds_clip_viral.bam",
        ref="reference/RefSeq_sequence_Wuhan-Hu1.fa"
    output:
        vcf="experiments/{exp}/results/variantCall/11_lofreqrefilt/{sample}_samp.lofreq_filtered.vcf",
        depth="experiments/{exp}/results/freyja/{sample}.depth"
    conda: 
        "envs/freyja.yaml"
    log:
        "logs/{exp}/freyja/{sample}.input.out"
    shell:  # generate depth file, unzip vcf, variants and demix
        "samtools mpileup -aa -A -d 600000 -Q 20 -q 0 -B -f {input.ref} {input.bam} | cut -f1-4 >{output.depth} &>{log} && "
        "gunzip -k {input.gz} "

rule freyja_main:
    input:
        vcf="experiments/{exp}/results/variantCall/11_lofreqrefilt/{sample}_samp.lofreq_filtered.vcf",
        bam="experiments/{exp}/results/variantCall/07_clipOverlap/{sample}_ds_clip_viral.bam",
        depth="experiments/{exp}/results/freyja/{sample}.depth",
        ref="reference/RefSeq_sequence_Wuhan-Hu1.fa"
    output:
        tsv="experiments/{exp}/results/freyja/{sample}.tsv"
    conda: 
        "envs/freyja.yaml"
    log:
        "logs/{exp}/freyja/{sample}.demix.out"
    shell:
        "freyja variants --ref {input.ref} --variants {input.vcf} --depths {input.depth} {input.bam} &>{log} && " 
        "freyja demix {input.vcf} {input.depth} --output {output.tsv} >>{log} 2>&1"


if config["VLQ"]["CREATE_REFSET"]:
    rule VLQ_generateRefset:
        input:
            metadata="reference/GISAID/metadata_tsv_2023_07_01.tar.xz",
            fasta="reference/GISAID/sequences_fasta_2023_07_01.tar.xz",
            ref="reference/RefSeq_sequence_Wuhan-Hu1.fa"
        output:
            gisaid_seq=temp("reference/GISAID/sequences.fasta"),
            gisaid_meta=temp("reference/GISAID/metadata.tsv"),
            metadata="reference/GISAID/VLQ_reference_set/metadata.tsv",
            kallisto_idx="reference/GISAID/VLQ_reference_set/sequences.kallisto_idx"
        params:
            seed = config["APP"]["SEED"],
            continent="Europe",
            ref_set="reference/GISAID/VLQ_reference_set"
        conda:
            "envs/VLQ.yaml"
        log:
            expand("logs/{exp}/VLQ/generateRefSet.out", exp=EXPERIMENT_NAME)
        shell: 
            """
            # build a reference set --results in sequences.fasta and metadata.tsv
            unxz {input.metadata} metadata.tsv {input.fasta} sequences.fasta >{log} &&
            cat {output.gisaid_meta} | perl -F"\t" -lane 'if($.==1 || $F[5] =~ m/20\d\d-\d\d-\d\d/) {{ print }}' > {output.gisaid_meta} &&

            # preprocessing reference sequencing data
            python bin/VLQ/wastewater_analysis/pipeline/preprocess_references.py -m {output.gisaid_meta} 
            -f {output.gisaid_seq} -k 1000 --seed {params.seed} --continent {params.continent} -o {params.ref_set}/ >>{log} && 

            # call variants compared to orginal sars cov-2 and compute allele frequencies per lineage
            # bash script does not recognize symbolic links, need orginal ABSOLUTE path
            bash pipeline/call_variants.sh {params.ref_set}/ {input.ref} >>{log} && 

            # select sequences per lineage with all mutation with al freq >50% were captured at least once
            python bin/VLQ/wastewater_analysis/pipeline/select_samples.py -m {output.gisaid_meta} -f {output.gisaid_seq} -o {params.ref_set}
            --vcf {params.ref_set}/*_merged.vcf.gz --freq {params.ref_set}/*_merged.frq >>{log} &&

            # index reference set
            kallisto index -i {output.kallisto_idx} {params.ref_set}/sequences.fasta >>{log} &&

            # clean-up
            rm {params.ref_set}/*_merged.vcf.gz {params.ref_set}/*_merged.frq {params.ref_set}/*_merged.frq {params.ref_set}/*_merged.log {params.ref_set}/*_merged.sites.pi &&
            rm -rf {params.ref_set}/*
            """
 
else:
    rule VLQ_dummy:
        output:
            metadata="reference/GISAID/VLQ_reference_set/metadata.tsv",
            kallisto_idx="reference/GISAID/VLQ_reference_set/sequences.kallisto_idx"

rule VLQ_main:
    input:
        kallisto_idx="reference/GISAID/VLQ_reference_set/sequences.kallisto_idx",
        metadata="reference/GISAID/VLQ_reference_set/metadata.tsv",
        fastq1="experiments/{exp}/results/variantCall/02_correctOverlapPE/{sample}_ecco_1.fq.gz",
        fastq2="experiments/{exp}/results/variantCall/02_correctOverlapPE/{sample}_ecco_2.fq.gz"
    output:
        abundance=temp("experiments/{exp}/results/VLQ/{sample}/abundance.tsv"),
        prediction="experiments/{exp}/results/VLQ/{sample}_predictions.tsv"
    params:
        threads = 20
    conda:
        "envs/VLQ.yaml"
    log:
        "logs/{exp}/VLQ/{sample}.out"
    shell:  # calculate abundance per reference sequence and obtain abundances
        "kallisto quant -t {params.threads} -i {input.kallisto_idx} -o experiments/{wildcards.exp}/results/VLQ/{wildcards.sample} -t {params.threads} {input.fastq1} {input.fastq2} &>{log} && "
        "python bin/VLQ/wastewater_analysis/pipeline/output_abundances.py -o {output.prediction} --metadata reference/GISAID/VLQ_reference_set/metadata.tsv {output.abundance} >>{log} 2>&1"

if config["LOLLIPOP"]["CREATE_REFSET"]:
    rule lollipop_mutationList:
        input:
            gff3="reference/gffs/Genes_NC_045512.2.GFF3"
        output: 
            mutlist="reference/GISAID/phe-genomics/phe2cojac/mutlist.tsv",
            pangovars="reference/GISAID/phe-genomics/phe2cojac/variants_pangolin.yaml"
        params:
            ref_dir="reference/GISAID/phe-genomics"
        shell:
            # fetch the repository of standardised variant definitions
            "mkdir {params.ref_dir} && cd {params.ref_dir} && git clone https://github.com/phe-genomics/variant_definitions.git && cd ../.. "

            # generate a YAML for lineages using the corresponding standardised variant definitions in cojac format

            # cojac phe2cojac {params.ref_dir}/variant_definitions/variant_yaml/*.yml -y {params.ref_dir}/phe2cojac/*.yml

            ## phe curated list now with new category (probable), not included in cojac
            # cojac sig-generate -u https://lapis.cov-spectrum.org/open/v1 --var B.1.324.1 > {params.ref_dir}/phe2cojac/coha.yml
            # cojac sig-generate -u https://lapis.cov-spectrum.org/open/v1 --var BA.5 > {params.ref_dir}/phe2cojac/dibu.yml
            # cojac sig-generate -u https://lapis.cov-spectrum.org/open/v1 --var XBB > {params.ref_dir}/phe2cojac/edta.yml
            # cojac sig-generate -u https://lapis.cov-spectrum.org/open/v1 --var B.1.617.2 > {params.ref_dir}/phe2cojac/emse.yml
            # cojac sig-generate -u https://lapis.cov-spectrum.org/open/v1 --var P.2 > {params.ref_dir}/phe2cojac/mado.yml
            "lollipop generate-mutlist --genes {input} --output {output.mutlist} --out-pangovars {output.pangovars} -- {params.ref_dir}/phe2cojac/*.yml"

else:
    rule lollipop_dummy:
        output: 
            mutlist="reference/GISAID/phe-genomics/phe2cojac/mutlist.tsv",
            pangovars="reference/GISAID/phe-genomics/phe2cojac/variants_pangolin.yaml"


rule lollipop_main_singleSample:
    input:
        bam="experiments/{exp}/results/variantCall/05_softclipPrimer/{sample}_sorted_viral_rh_trim.bam",
        mutlist="reference/GISAID/phe-genomics/phe2cojac/mutlist.tsv",
    output:
        basecnt=temp("experiments/{exp}/results/lollipop/{sample}.basecnt.tsv.gz"),
        cov=temp("experiments/{exp}/results/lollipop/{sample}.coverage.tsv.gz"),
        mut="experiments/{exp}/results/lollipop/{sample}.mut.tsv"
    conda:
        "envs/lollipop.yaml"
    params:
        sample_date=lambda wildcards: get_sample_date(EXPERIMENT_NAME, f"{wildcards.sample}"),
        location=config["LOLLIPOP"]["LOCATION_NAME"]
    log:
        "logs/{exp}/lollipop/{sample}.out"
    shell:
        # search mutation in a single sample - create basecount table
        "aln2basecnt --basecnt {output.basecnt} --first 1 --coverage {output.cov} --name \"{wildcards.sample}\" {input.bam} &>{log} && "
        "lollipop getmutations from-basecount --based 1 --output {output.mut} "
        "--location \"{params.location}\" --date \"{params.sample_date}\" -m {input.mutlist} -- {output.basecnt} >>{log} 2>&1"

rule lollipop_main_combinedSeries:
    input:
        pangovars="reference/GISAID/phe-genomics/phe2cojac/variants_pangolin.yaml",
        decon_config="bin/lollipop/" + config["LOLLIPOP"]["KERNEL_DECONVOLUTION_CONFIG"],
        s_mut= lambda wildcards: expand("experiments/{exp}/results/lollipop/{sample}.mut.tsv", exp=EXPERIMENT_NAME, sample=SAMPLES)
    output:
        mutgz=temp("experiments/{exp}/results/lollipop/tallymut.tsv.gz"),
        decon="experiments/{exp}/results/lollipop/deconvoluted.tsv"
    params: 
        seed=config["APP"]["SEED"]
    conda:
        "envs/lollipop.yaml"
    log:
        "logs/{exp}/lollipop/deconvolution.out"
    shell:
        # combine time series and run deconvolution
        "xsv cat rows {input.s_mut} | xsv fmt --out-delimiter '\\t' | gzip >{output.mutgz} && "
        "lollipop deconvolute --output {output.decon} --var {input.pangovars} --seed={params.seed} --dec {input.decon_config} -- {output.mutgz} &>{log}"