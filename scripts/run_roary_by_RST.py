import os
import glob
import subprocess

# set the base directory for the assemblies
base_datadir = '/home/mf019/assemblies/longread'

# get CWD
cwd = os.getcwd()

# Define the folder names
rst_types = ["RST1", "RST2", "RST3"]


# Define the path to the Roary Apptainer image
roary_image = "roary_latest.sif"

# Function to run Roary using Apptainer
def run_roary(input_folder, output_folder):
    gffs = glob.glob(f'{input_folder}/*.gff3')
    command = [
        "apptainer", "exec",
        "--bind", f"{input_folder}:/data/input:r", # bind input to container, give explicit read perms
        "--bind", f"{output_folder}:/data/output:rw", # bind output to container, give explicit rw perms
        roary_image, "roary", # call container image and then run roary inside.
        "-p", "32", "-z", # 32 threads, keep int files
        "-e", "--mafft", # prank aln, mafft
        "-f", "/data/output", # Output folder
    ] + [f"/data/input/{os.path.basename(file)}" for file in gffs]
    subprocess.run(command, check=True)

# Loop through the folders and run Roary
for rst in rst_types:
    input_folder = os.path.join(base_datadir, rst)
    output_folder = os.path.join(cwd, f'roary_{rst}')

    # Create the output directory if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    os.chmod(output_folder, 0o777) # ugh change perms so roary can write to output
    print(f'{input_folder}\n{output_folder}')

    # Run Roary
    run_roary(input_folder, output_folder)

print("Roary pangenome analysis completed for all folders.")
