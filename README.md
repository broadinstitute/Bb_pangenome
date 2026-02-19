# Long read assemblies reveal Borrelia burgdorferi population variation structured by plasmids, plasmid subtypes, and genetic networks

Preprint: https://doi.org/10.1101/2025.01.29.635312

This repository contains all code and sequences used in the above manuscript.

Ensure you have [git lfs](https://git-lfs.com/) installed to access all files and data. (Instructions below)

**Total repository size: 27G**

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
    - `genotyping/` - contains OspC and MLST typing results as well as plasmid calls.
    - `homology_networks/` - contains homology network graphs and network.json used to render graph in Blender.
    - `reports/` - contains AGAT stats, kraken2, quast, and multiqc reports.
    - `results/` - contains roary and panaroo pangenome output for split and non-split paralogs.
- `ref/` - contains database of assemblies, B31 reference genome, and replicon references used in classification.
- `scripts/` - contains scripts and commands used in various analyses.
- `VAE/` - All code and data for variable autoencoder analyses
- `scaffolding/` - contains data for scaffolding analyses.
- `figures/` - contains notebooks, data, scripts, and relevant files to reproduce published figures. Paths in notebooks are relative to repository root.
- `interactive/` - interactive RShiny app for Figure 2 and Figure 6. Requires R-Studio.

***

In situations where code or data are unable to be located, refer to our [working repository](https://github.com/mjfos2r/longread_pangenome) or open an issue!

*** 

## Using LFS to access files and data
To get started, clone the repository and setup LFS:
```bash
git clone https://github.com/broadinstitute/Bb_pangenome.git
cd Bb_pangenome
git lfs install 
```

To retrieve LFS tracked files individually:
```bash
git lfs fetch path/to/lfs-tracked-file.ext
git lfs checkout
```

To retrieve LFS tracked directories:
```bash
git lfs fetch path/to/lfs-tracked-directory/
git lfs checkout
```

To retrieve all files and data:
```bash
git lfs fetch --all
git lfs checkout
```
