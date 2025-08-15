#!/bin/bash
# Michael J. Foster
# github.com/mjfos2r
# This is a short helper script to automate the command for getting MLST sequence classifications using 
# docker.

assemblies='assemblies/**/*.fna'
outdir='mlst/results/'

# first thing we need to do is pull all of the fasta files together in the same directory. 
# we'll just do this temporarily and wipe it after execution.

mkdir -p tmp_fastas
find assemblies -type f -name "*.fna" -exec cp {} tmp_fastas/ \;
# and since it's given me trouble on non fasta extensions, let's rename everything from .fna to .fasta

for i in tmp_fastas/*.fna; do
    mv "${i}" "${i%.*}.fasta"
done

# Container is sangerpathogens/mlst_check:sha256:f328196956a0cf9bdeffaaf7fb8d1f1fd24f45659ab3d1b979b52b4456be9b6b
# output path needs /data/ appended so that it can be located within the container.
docker run --rm -v $(pwd):/data -it sangerpathogens/mlst_check \
    get_sequence_type -c \
    -s 'Borrelia' \
    -d $(nproc) \
    -o "${OUTDIR}" \
    tmp_fastas/*.fasta
