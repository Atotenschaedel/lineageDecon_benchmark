# lineageDecon_benchmark
Pipeline for assessment of deconvolution tools for virus lineage abundance estimation from wastewater sequencing


# How to run
Exact parameter settings are to be set in workflow_config.yaml

Workflow consists out of three separate snakemake workflows to be exectuded in series

1) simulation.smk  - simulate sequencing dataset
2) variantCalling.smk - variant calling pipeline
3) lineage deconvolution and post-prediction analysis

R markdown file PostPredict_plots.Rmd used to further post-prediciton analysis and plotting of results 
