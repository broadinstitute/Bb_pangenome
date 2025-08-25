import os
import sys
import csv
import glob
import re
import subprocess
import tempfile
from Bio import SeqIO
from collections import defaultdict
import argparse
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from queue import Queue
from threading import Thread
from typing import Dict, List, Tuple
import json
import fcntl
import select
import time

def setup_dialog():
    """Setup named pipe and dialog for progress display"""
    pipe_path = "/tmp/homology_progress"
    if os.path.exists(pipe_path):
        os.unlink(pipe_path)
    os.mkfifo(pipe_path)
    
    # Start dialog in background
    dialog_cmd = f"dialog --gauge 'Processing files...' 10 70 0 < {pipe_path}"
    dialog_process = subprocess.Popen(dialog_cmd, shell=True)
    return pipe_path, dialog_process

def update_progress(pipe_path, percent, message):
    """Update dialog progress"""
    try:
        with open(pipe_path, 'w') as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            f.write(f"XXX\n{percent}\n{message}\nXXX\n")
            f.flush()
            fcntl.flock(f, fcntl.LOCK_UN)
    except IOError:
        pass

def parse_genbank(path: str):
    try:
        records = list(SeqIO.parse(path, 'genbank'))
        return records[0] if records else None
    except Exception as e:
        print(f"Error parsing GenBank file {path}: {str(e)}", file=sys.stderr)
        return None

def process_alignment(alignment_data: str, genbank_paths: Dict) -> Dict:
    """Process a single alignment entry from the pipeline"""
    try:
        # Parse the alignment path and data
        alignment_file, content = alignment_data.split('\0', 1)
        asm1_id, asm2_id = os.path.basename(os.path.dirname(alignment_file)).split('_vs_')
        
        asm1_record = parse_genbank(genbank_paths[asm1_id])
        asm2_record = parse_genbank(genbank_paths[asm2_id])
        
        if not content.strip():
            return None
            
        # Parse alignment data
        lines = content.splitlines()
        if not lines:
            return None
            
        keys = lines[0].strip().split('\t')
        alignments = [dict(zip(keys, line.strip().split('\t'))) for line in lines[1:]]
        results = []
        
        for aln in alignments:
            try:
                results.append({
                    'asm1_id': asm1_id,
                    'asm1_coords': f"{aln['QUERY_START']}-{aln['QUERY_END']}",
                    'asm2_id': asm2_id,
                    'asm2_coords': f"{aln['REF_START']}-{aln['REF_END']}",
                    'identity': float(aln['IDENTITY'])
                })
            except (KeyError, ValueError):
                continue
                
        return results
    except Exception as e:
        print(f"Error processing alignment: {e}", file=sys.stderr)
        return None

def process_stream(pipe, genbank_paths: Dict, progress_pipe: str, total_files: int):
    """Process streaming input from find | pv pipeline"""
    results = []
    processed = 0
    current_file = []
    
    while True:
        ready, _, _ = select.select([pipe], [], [], 1.0)
        if not ready:
            if pipe.closed:
                break
            continue
            
        line = pipe.readline()
        if not line:
            break
            
        current_file.append(line)
        if line.strip() == "":  # Empty line indicates end of file
            if current_file:
                file_content = "".join(current_file)
                result = process_alignment(file_content, genbank_paths)
                if result:
                    results.extend(result)
                
                processed += 1
                percent = min(int((processed / total_files) * 100), 100)
                update_progress(progress_pipe, percent, 
                              f"Processed {processed}/{total_files} files ({percent}%)")
                
            current_file = []
    
    return results

def main():
    parser = argparse.ArgumentParser(description="Homology Analysis with Progress Tracking")
    parser.add_argument('--annotations_dir', required=True)
    parser.add_argument('--alignments_dir', required=True)
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--cores', type=int, required=True)
    args = parser.parse_args()

    # Setup progress tracking
    progress_pipe, dialog_process = setup_dialog()
    
    try:
        # Count total files with progress
        count_cmd = f"find {args.alignments_dir} -name 'align_coords.tsv' | pv -l -N 'Counting files' | wc -l"
        total_files = int(subprocess.check_output(count_cmd, shell=True, stderr=subprocess.PIPE).decode().strip())
        
        update_progress(progress_pipe, 0, f"Found {total_files} files to process")
        
        # Index GenBank files
        genbank_paths = {}
        for path in glob.glob(f'{args.annotations_dir}/*.gbff'):
            base_name = os.path.basename(path).split('.')[0]
            genbank_paths[base_name] = path
        
        # Create processing pipeline
        find_cmd = f"find {args.alignments_dir} -name 'align_coords.tsv' -type f"
        xargs_cmd = f"xargs -I % sh -c 'echo -n \"%\\0\"; cat \"%\"; echo'"
        pv_cmd = f"pv --line-mode --size {total_files} --name 'Processing:'"
        
        pipeline = f"{find_cmd} | {pv_cmd} | {xargs_cmd}"
        
        # Start pipeline and process results
        process = subprocess.Popen(pipeline, shell=True, stdout=subprocess.PIPE, 
                                 stderr=sys.stderr, text=True, bufsize=1)
        
        results = process_stream(process.stdout, genbank_paths, progress_pipe, total_files)
        
        # Write results
        os.makedirs(args.output_dir, exist_ok=True)
        output_path = os.path.join(args.output_dir, 'homology_results.tsv')
        
        update_progress(progress_pipe, 100, "Writing results...")
        
        # Use pv to show write progress
        df = pd.DataFrame(results)
        temp_path = tempfile.mktemp()
        df.to_csv(temp_path, sep='\t', index=False)
        
        file_size = os.path.getsize(temp_path)
        subprocess.run(f"pv -s {file_size} {temp_path} > {output_path}", shell=True)
        os.unlink(temp_path)
        
        update_progress(progress_pipe, 100, "Analysis complete!")
        time.sleep(1)  # Give dialog time to show final message
        
    finally:
        # Cleanup
        dialog_process.terminate()
        if os.path.exists(progress_pipe):
            os.unlink(progress_pipe)

if __name__ == "__main__":
    main()