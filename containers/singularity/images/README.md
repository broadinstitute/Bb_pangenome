# Singularity Images

To preserve space in the repository, the images used have been listed in `singularity_images.txt`.
The images that were defined can be built using the definitions present in `containers/singularity/defs`

To perform the conversion of docker images into singularity images, you can execute the following:

```bash
singularity build <image-name>.sif docker://dockeruser/image-name:tag
```

For example, to download the roary container and convert to singularity:

```bash
singularity build roary_latest.sif docker://sangerpathogens/roary:latest
```

## Images Used

agat: `quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0`
bakta: `https://hub.docker.com/r/oschwengers/bakta:1.10.1`
fasttree2: `containers/singularity/defs/fasttree2.def`
genometools: `containers/singularity/defs/genometools.def`
kraken2: `containers/singularity/defs/kraken2.def`
mlst-check: `https://hub.docker.com/r/sangerpathogens/mlst_check:latest`
multiqc: `https://hub.docker.com/r/staphb/multiqc:latest`
quast: `https://hub.docker.com/r/staphb/quast:latest`
raxml: `https://hub.docker.com/r/pegi3s/raxml:latest`
roary: `https://hub.docker.com/r/sangerpathogens/roary:latest`
