import glob

configfile: 'configs/pangenome_analysis.yaml'

analysis_dir = config['analysis_dir']
gffs_dir = config['gffs_dir']
version = config['version']
type = config['type']
total_threads = config['threads']

rule all:
    input:
        #f'{analysis_dir}/output/results/{version}/{version}_{type}/roary_{version}_{type}/core_gene_alignment.aln',
        #f'{analysis_dir}/output/results/{version}/{version}_{type}/fasttree/fasttree_roary_pangenome_{version}_{type}.newick',
        #f'{analysis_dir}/output/results/{version}/{version}_{type}/tree/RAxML_bestTree.roary_{version}_{type}_gtrgam',
        f'{analysis_dir}/output/results/{version}/{version}_{type}/bakta_{version}_{type}/roary_pangenome_{version}_{type}.gff3'

rule Roary:
    input:
        glob.glob(f'{gffs_dir}/*.gff3')
    output:
        outfile = f'{analysis_dir}/output/results/{version}/{version}_{type}/roary_{version}_{type}/core_gene_alignment.aln',
        outasm = f'{analysis_dir}/output/results/{version}/{version}_{type}/roary_{version}_{type}/pan_genome_reference.fa',
    params:
        opts = config['roary']['opts'],
        anal = f'{analysis_dir}/output/results/{version}',
        outdir = f'{analysis_dir}/output/results/{version}/{version}_{type}/roary_{version}_{type}',
        temp_dir = f'{analysis_dir}/output/results/{version}/{version}_{type}/roary_{version}_{type}_temp'
    threads: total_threads
    message: "running roary on all of our annotations! : \n roary -p {threads} -f {params.outdir} {input} 2>&1 >{log}"
    log: f'{analysis_dir}/output/results/{version}/{version}_{type}/logs/roary_log.txt'
    singularity: 'containers/singularity/images/roary_latest.sif'
    shell: # -s to not split paralogs
        """
        roary {params.opts} -p {threads} -f {params.temp_dir} {input} >{log} 2>&1
        # move results from the roary_temp dir to the snakemake dir.
        # I hate snakes
        ls {params.anal}
        mv {params.temp_dir}*/* {params.outdir}/
        rmdir {params.temp_dir}*
        """
rule FastTree:
    input:
        f'{analysis_dir}/output/results/{version}/{version}_{type}/roary_{version}_{type}/core_gene_alignment.aln'
    output:
        f'{analysis_dir}/output/results/{version}/{version}_{type}/fasttree/fasttree_roary_pangenome_{version}_{type}.newick'
    params:
        log_dir = f'{analysis_dir}/output/results/{version}/{version}_{type}/fasttree/logs'
    threads: total_threads
    singularity: 'containers/singularity/images/fasttree2.sif'
    message: "running fasttree on the roary output ==>> {input}"
    log: f'{analysis_dir}/output/results/{version}/{version}_{type}/logs/fasttree.txt'
    shell:
        """
        export OMP_NUM_THREADS={threads}
        fasttreeMP -gamma -nt -gtr -log {params.log_dir} {input} > {output} 2>{log}
        """

rule RAxML:
    input:
        f'{analysis_dir}/output/results/{version}/{version}_{type}/roary_{version}_{type}/core_gene_alignment.aln'
    output:
        outdir = directory(f'{analysis_dir}/output/results/{version}/{version}_{type}/tree/'),
        outfile = f'{analysis_dir}/output/results/{version}/{version}_{type}/tree/RAxML_bestTree.roary_{version}_{type}_gtrgam'
    params:
        prefix = f'roary_{version}_gtrgam'
    threads: total_threads
    singularity: 'containers/singularity/images/raxml_latest.sif'
    message: "running RAxML on the roary output ==>> {input}"
    log: f'{analysis_dir}/output/results/{version}/{version}_{type}/logs/RAxML_log.txt'
    shell:
        "raxmlHPC -m GTRGAMMA -p 71235 -s {input} -n {params.prefix} -w {output.outdir} 2>&1 >{log}"

rule bakta:
    input: f'{analysis_dir}/output/results/{version}/{version}_{type}/roary_{version}_{type}/pan_genome_reference.fa'
    output:
        outdir = directory(f'{analysis_dir}/output/results/{version}/{version}_{type}/bakta_{version}_{type}/'),
        outgff = f'{analysis_dir}/output/results/{version}/{version}_{type}/bakta_{version}_{type}/roary_pangenome_{version}_{type}.gff3',
    params:
        gram = config['bakta']['gram'],
        opts = config['bakta']['opts'],
        prefix = f'roary_pangenome_{version}_{type}',
    threads: total_threads
    log: f'{analysis_dir}/output/results/{version}/{version}_{type}/logs/bakta.log'
    singularity: "containers/singularity/images/bakta_latest.sif"
    message:
        "Running bakta for pangenome!"
    shell:
        "(bakta --db /db {params.gram} --threads {threads} "
        "{params.opts} --output {output.outdir} "
        "--prefix {params.prefix} --genus Borrelia --species burgdorferi --strain sensu-stricto --keep-contig-headers -v --skip-plot {input}) 2>&1 >{log} &&"
        "sed -i -E '/^##feature-ontology/s|http[s]?://[^ ]+|https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/refs/tags/v3.1/so.obo|' {output.outgff}"
