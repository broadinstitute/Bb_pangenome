import pandas
import glob

configfile: local("configs/singularity_config.yaml")
analysis_dir = config['analysis_dir']
input_dir = config['input_dir']

rule all:
    input:
        f"{analysis_dir}/genotyping/mlst/mlst_results.allele.csv"

##
# MLST Type each assembly!
## 

rule getMLST_type:
    input: glob.glob(f"{input_dir}/*/*.fna")
    output:
        outdir = directory(f"{analysis_dir}/genotyping/mlst/"),
        final_out = f"{analysis_dir}/genotyping/mlst/mlst_results.allele.csv"
    singularity: "singularity/images/mlst_check_latest.sif"
    log: f"{analysis_dir}/logs/mlst_classification.log"
    threads: 360
    message: "running MLST Check for the following files: {input}"
    shell:
        """
        (get_sequence_type -d {threads} -s 'Borrelia' -c -o {output.outdir} {input}) 2>&1 >{log}
        """

