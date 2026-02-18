import os
import subprocess
from datetime import datetime
import pickle
import argparse

# General code for annotating
def call_bakta(outpt, prefix, inpt):
    result = subprocess.run([
        "docker", "run", "--rm", "-t", # "-it",
        "-v", "/home/mf019/db/bakta:/db",
        "-v", "/home/kotzen/scaffolding_analysis:/data",
        "mjfos2r/bakta:latest",
        "--db", "/db",
        "--gram", "-",
        "--threads", "12",
        "--verbose",
        "--force",
        "--output", outpt,
        "--prefix", prefix,
        "--keep-contig-headers",
        "--skip-plot",
        inpt
    ], capture_output=True, text=True)
    return result

# Annotate genes
def annotate_genes(input_fasta_dir, file_suffix, trouble_suffix):
    trouble_asms = {}
    
    input_fastas = os.listdir(input_fasta_dir)
    for file in input_fastas:
        iso = file.replace(file_suffix, '')
        print(iso)
        print("Start time:", datetime.now().strftime("%H:%M:%S"))
        prefix = iso
        outpt = f"/data/{input_fasta_dir.split('/')[0]}/bakta/{iso}"
        inpt = f"/data/{input_fasta_dir}/{file}"
        result = call_bakta(outpt, prefix, inpt)

        if result.stdout != None or result.stderr != None:
            trouble_asms[iso] = result

    if len(trouble_asms)>0:
        with open(f'trouble_{trouble_suffix}.pkl', 'wb') as f:
            pickle.dump(trouble_asms, f)
        
        raise ValueError(trouble_asms.keys())

def main(args):
    if args.SRA:
        annotate_genes('SRAs/fasta', '.fasta', 'SRA')

    if args.scaf:
        annotate_genes('scaffolded_asms/fasta', '_scaffold.fasta', 'scaf')

    if args.LRA:
        annotate_genes('LRAs/fasta', '.fna', 'LRA')

    else:
        annotate_genes('SRAs/fasta', '.fasta', 'SRA')
        annotate_genes('scaffolded_asms/fasta', '_scaffold.fasta', 'scaf')
        annotate_genes('LRAs/fasta', '.fna', 'LRA')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Handle SRA, scaf, and LRA flags.")

    parser.add_argument('--SRA', action='store_true', help='Enable SRA processing')
    parser.add_argument('--scaf', action='store_true', help='Enable scaffold processing')
    parser.add_argument('--LRA', action='store_true', help='Enable LRA processing')

    args = parser.parse_args()

    main(args)
    
    
        


