# Configs README

Here are the configuration files required to run the included snakemake workflows.

## `pangenome_analysis.yaml`

used for the following snakemake workflows:

- `../annotate_all_assemblies_v3.smk`
- `../kraken2_QC.smk`
- `../pan_genome_analysis_v5.smk`

To run these workflows, ensure you have the following:

- All relevant containers are downloaded and converted to singularity/apptainers. (Alternatively, update the snakefiles to use docker.)
- ensure paths are changed relative to their location in your working directory.
- runtime resources are changed to reflect the machine this workflow is being executed on.

## `singularity_config.yaml`

used for the following snakemake workflow:
- `../get_mlst_type.smk`.

To run this workflows, ensure the following:

- all gff3 files are in a single directory at the same level and this path is specified.
- all paths are changed relative to your locations in your working directory.
- resources are updated to reflect the machine this workflow is being executed on.
