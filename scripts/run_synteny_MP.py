import os
import sys
import glob
import math
import argparse
import subprocess
from tqdm import tqdm
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

def get_fasta_inputs(in_dir):
    genomes = glob.glob(f'{in_dir}/*.fasta')
    return genomes

def get_gb_inputs(in_dir):
    genomes = glob.glob(f'{in_dir}/*.gbff')
    if len(genomes) == 0:
        genomes = glob.glob(f'{in_dir}/**/*.gbff', recursive=True)
    return genomes

def set_output(in_file, out_dir):
    isolate_name = os.path.basename(in_file).replace('.gbff','')
    return os.path.join(out_dir, isolate_name)

def set_multi_output(file1, file2, out_dir):
    plasmid1_id = os.path.basename(file1).replace('.gbff','')
    plasmid2_id = os.path.basename(file2).replace('.gbff','')
    filename = f'{plasmid1_id}_vs_{plasmid2_id}'
    return os.path.join(out_dir, filename)

def get_alignment_command(file1, file2, output_path, prog):
    # Construct the alignment command based on the chosen program
    if prog == "mauve":
        return get_mauve_command(file2, file1, output_path)
    elif prog == "mummer":
        return get_mummer_command(file2, file1, output_path)
    elif prog == "mmseq":
        return get_mmseq_command(file2, file1, output_path)

def get_mauve_command(file1, file2, output_path):
    # Set up the command to run one genome vs the B31 reference.
    isolate_name = file1.split('/')[-1].replace('.gbff','')
    out_file = f'{output_path}/{isolate_name}'
    progressive_mauve_path = "progressiveMauve"
    command = [f'{progressive_mauve_path}', f'--output={out_file}/{isolate_name}.xmfa',
                 f'--backbone-output={out_file}.backbone',
                 f'--output-guide-tree={out_file}.gtree',
                 f'{file2}', f'{file1}']
    return command

def get_mummer_command(file1, file2, output_path):
    # Set up the command to run one genome vs another.
    # Changed params for length_thr from 1000 to 250
    # Changed params to not output plot and png 
    command = ['pgv-mummer', '-q', '--threads=2', f'{file2}', f'{file1}', '-o', output_path, '--formats=png',]
    #'--show_scale_bar', '--curve',] # '--length_thr=250', '--identity_thr=75', apparently both default to 0 lol.
    return command

def get_mmseq_command(file1, file2, output_path):
    # Set up the command to run one genome vs another.
    command = ['pgv-mmseqs', '-q', '--threads=1', f'{file2}', f'{file1}', '-o', output_path, '--formats=png',]# '--show_scale_bar'] 
               #'--length_thr=250', '--identity_thr=75', '--show_scale_bar',]
    return command

def run_command(command):
    try:
        # Run the command
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        return f'Successfully completed: {' '.join(command)}' if result.returncode == 0 else result.stderr
    
    except subprocess.CalledProcessError as e:
        return f"Command '{command}' failed with error: {e.stderr}"
    
    except Exception as e:
        return f"An unexpected error occurred: {str(e)}"

def parse_results(results):
    # Parse out the useful results?
    # visualize? ## ITS LITERALLY ALREADY DOING THIS :)
    # This will actually probably be the updated multiparser stuff I've written in the other file!
    pass
    
#def run_alignment(prog, genome, output_path, reference):
#    if prog == "mauve":
#        command = get_mauve_command(genome, output_path, reference)
#        run_mauve(command)
#        return f"progMauve against {os.path.basename(reference)} for {genome} finished!"
#    elif prog == "mummer":
#        command = get_mummer_command(genome, output_path, reference)
#        run_pgv_mummer(command)os.path.basename(
#   )gv-mummer command for {genome} against {os.path.basename(reference)} finished!"

def run_alignment(command):
    # Execute the alignment command
    try:
        run_command(command)
        return f"Finished alignment of {' '.join(command)}"
    except Exception as e:
        return f"Error running {' '.join(command)}"

def progress_bar(futures):
    with tqdm(total=len(futures)) as pbar:
        for future in as_completed(futures):
            try:
                result = future.result()
                tqdm.write(result)
            except Exception as e:
                tqdm.write(f"Error: {e}")
            pbar.update(1)

def main():
    # Create the arg parser
    parser = argparse.ArgumentParser(description="A script to run mauve on a directory of assemblies against the B31 reference genome.")
    
    # Add the arguments
    parser.add_argument('--mode',      type=str, help='Please specify which mode to run in: all_v_all or all_v_one [--ava, --av1]', required=True)
    parser.add_argument('--prog',      type=str, help='type of alignment to run [mauve, mummer, mmseq]', required=True)
    parser.add_argument('--input_dir', type=str, help='The directory of inputs', required=True)
    parser.add_argument('--reference', type=str, help='path to reference genome for alignment (if running av1 mode)')
    parser.add_argument('--output',    type=str, help='The directory for outputs', required=True)
    parser.add_argument('--cores',     type=int, help='how many cores we boggin?', required=True)
    
    # Parse the arguments
    args = parser.parse_args()
    
    print(f"Parsed arguments: {args}", file=sys.stderr) # Log parsed arguments for debugging
    if args.mode == 'ava':
        print("running in all_vs_all mode!")
    elif args.mode == 'av1':
        print("running in all_vs_one mode!")
    else:
        parser.error("ERROR: selected mode is not allowed. try again with --ava or --av1!")
        
    if args.mode == 'av1' and not args.reference:
        parser.error("--reference is required when running in all_v_one (av1) mode")
                        
    # figure out how many workers to spin up based on 4 threads per instance of mmseqs2. round num workers down. 
    num_workers = math.floor(args.cores/2) # commented out but this is preferable. Ugh. now that mmseqs is on path we should be groovy.
    
    # Print the arguments
    print(f"alignment mode {args.mode}")
    print(f"alignment type {args.prog}")
    print(f"Input dir {args.input_dir}")
    print(f"Output dir: {args.output}")
    print(f"Cores: {args.cores}")
    print(f"Workers: {num_workers}")
    output_dir = args.output.strip()

    #if args.prog == "mauve":
    #    genomes = get_fasta_inputs(args.input_dir)
    #elif args.prog == "mummer" or args.prog == "mmseq":
    #    if args.input_dir.endswith('.gbff'):
    #        genomes = [args.input_dir]
    #    else:
    #        genomes = get_gb_inputs(args.input_dir)
    #else:
    #    raise ValueError("Unsupported alignment type. Use 'mauve' 'mmseq' or 'mummer'.")

    if args.mode == 'ava':
        if args.reference:
            print("ERROR: single reference should not be provided in all_vs_all mode")
            return 0
        genomes = get_gb_inputs(args.input_dir)
        print(f"Number of genomes being processed: {len(genomes)}")
        # okay let's get this pool party started.
        with ProcessPoolExecutor(max_workers=num_workers) as executor: # execute the pool party
            futures = []
            #job_count = 0
            for i, genome1 in enumerate(genomes): # big loop through all genomes
                for j, genome2 in enumerate(genomes): # in
                    if i < j: # process each pair only once, do not compare against self.
                        output_path = set_multi_output(genome1, genome2, output_dir)
                        command = get_alignment_command(genome1, genome2, output_path, args.prog)
                        futures.append(executor.submit(run_alignment, command))
                        #print(f"job:{job_count} added!")
                        #job_count += 1
            progress_bar(futures)
            
    elif args.mode == 'av1':
        genomes = get_gb_inputs(args.input_dir)
        print(f"Number of genomes being processed: {len(genomes)}")
        reference = args.reference
        print(f"Reference: {reference}")
                                    
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = []
            for genome in genomes:
                output_path = set_output(genome, output_dir)
                command = get_alignment_command(genome, reference, output_path, args.prog)
                futures.append(executor.submit(run_alignment, command))
            progress_bar(futures)
            
        
    ##print(f"Cores for Mauve: {args.cores_for_mauve}")
    #match args.prog:
    #    case "mauve":
    #        genomes = get_fasta_inputs(args.input_dir)
    #        for genome in genomes:
    #            output_path = set_output(genome, args.output_dir)
    #            command = get_mauve_command(genome, output_path, args.reference)
    #            run_mauve(command)
    #            print(f"progMauve against {ref_name} for {genome} finished!")
    #    case "mummer":
    #        if args.input_dir.endswith('.gbff'):
    #            genomes = args.input_dir
    #        else:
    #            genomes = get_gb_inputs(args.input_dir)
    #        for genome in genomes:
    #            output_path = set_output(genome, args.output_dir)
    #            command = get_mummer_command(genome, output_path, args.reference)
    #            print(' '.join(command))
    #            print(f"running pgv-mummer for {genome} against {ref_name}!")
    #            run_pgv_mummer(command)
    #            print(f"pgv-mummer command finished!\n")
#
if __name__ == "__main__":
    main()
