# Complex exchanges among plasmids and clonal expansion of lineages shape the population structure and virulence of Borrelia burgdorferi

This repository contains all code and sequences used in the above manuscript.

Ensure you have [git lfs](https://git-lfs.com/) installed to access all files and data.

**Total repository size: 23G**


## Directory Structure

- `assemblies/` - all assemblies and annotations used in our analyses are housed here. The subdirectory name is `Strain_ID` from table_s1.
    - Each assembly has the following contents:
        - `*.embl`
        - `*.faa`
        - `*.ffn`
        - `*.fna`
        - `*.gbff`
        - `*.gff3`
        - `*.hypotheticals.faa`
        - `*.hypotheticals.tsv`
        - `*.inference.tsv`
        - `*.json`
        - `*.log`
        - `*.tsv`
        - `*.txt`
- `containers/` - all containers used in these analyses are defined here. Contains submodule: [mjf-containers](https://github.com/mjfos2r/containers.git)
- `group2BB/` - Code used to map roary pangene groups back to B31 genes.
- `metadb/` - Code used to setup and query an sqlite3 database containing all assemblies and their annotations.
- `notebooks/` - Jupyter notebooks used throughout development. Some notebooks may be superceded by scripts present in `scripts/`
- `output/` - Contains all output used in the analysis and figure generation.
    - `alignments/` - contains alignments for all vs all, all vs B31.
    - `genotyping/` - contains OspC, MLST, RST typing results and plasmid calls.
    - `homology_networks/` - contains homology network graphs and network.json used to render graph in Blender.
    - `reports/` - contains AGAT stats, kraken2, quast, and multiqc reports.
    - `results/` - contains roary pangenome output for split and non-split paralogs.
- `ref/` - contains database of assemblies, B31 reference genome, and replicon references used in classification.
- `scripts/` - contains scripts and commands used in various analyses.
- `snakefiles/` - contains snakefiles and their configs.

README in progress.
