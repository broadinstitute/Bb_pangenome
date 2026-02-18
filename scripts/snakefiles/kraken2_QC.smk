import pandas
import glob

configfile: local("configs/pangenome_analysis.yaml")
analysis_dir = config['analysis_dir']
version = config['version']

ASSEMBLIES = [os.path.basename(os.path.dirname(f)) for f in glob.glob(f"{analysis_dir}/assemblies/expanded_v2/*/*.fna")]
print(ASSEMBLIES)
rule all:
    input:
        expand(f'{analysis_dir}/kraken2/reports/{{assembly}}.report', assembly=ASSEMBLIES),
        expand(f'{analysis_dir}/kraken2/reports/{{assembly}}.output', assembly=ASSEMBLIES),
        f'{analysis_dir}/kraken2/kraken2_report.krona',
        f'{analysis_dir}/kraken2/kraken2_report.krona.html'

rule run_Kraken2:
    input: f'{analysis_dir}/assemblies/expanded_v2/{{assembly}}/{{assembly}}.fna'
    output:
        output = f'{analysis_dir}/kraken2/reports/{{assembly}}.output',
        report = f'{analysis_dir}/kraken2/reports/{{assembly}}.report'
        #outdir = directory(f'{analysis_dir}/kraken2/reports/{{assembly}}/')
    singularity: 'singularity/images/kraken2_latest.sif'
    log: f'{analysis_dir}/logs/kraken2/{{assembly}}_kraken2.log'
    params:
        db = config['kraken2_db']
    threads: 32
    shell:
        """
        kraken2 --use-names --memory-mapping --db /kraken2_dbs \
            --threads {threads} \
            --output {output.output} \
            --report {output.report} \
            {input}
        """

rule krona_plot:
    input:
        report = f'{analysis_dir}/kraken2/reports/{{assembly}}.report'
    output:
        krona = f'{analysis_dir}/kraken2/krona/{{assembly}}.krona',
    singularity: 'singularity/images/kraken2_latest.sif'
    threads: 8
    shell:
        """
        kreport2krona.py -r {input} -o {output.krona}
        """

rule html_plot:
    input:
        krona = expand(f'{analysis_dir}/kraken2/krona/{{assembly}}.krona', assembly=ASSEMBLIES)
    output:
        html = f'{analysis_dir}/kraken2/kraken2_report.krona.html',
        report = f'{analysis_dir}/kraken2/kraken2_report.krona'
    singularity: 'singularity/images/kraken2_latest.sif'
    threads: 8
    params:
        r_files = lambda wildcards, input: " ".join([""] + input.krona)
    shell:
        """
        ktImportText {params.r_files} -o {output.html}
        """