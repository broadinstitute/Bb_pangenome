# ENA Assembly Submission Toolkit

Scripts for submitting genome assemblies to EMBL-ENA via webin-cli.

## Prerequisites

- Java 11+
- Python 3.10+
- webin-cli.jar ([download latest](https://github.com/enaDataSubmission/webin-cli/releases))
- ENA Webin account credentials

```bash
# Download webin-cli
curl -L -o webin-cli.jar \
  https://github.com/enaDataSubmission/webin-cli/releases/latest/download/webin-cli.jar
```

## Input Files

1. **Metadata TSV** — genome-info table with columns:
   `assembly_name, biosample_accession, genome_coverage, assembly_methods, sequencing_technologies, filename`

2. **Accession mapping** (optional) — links assembly_name to bioproject/biosample/SRA accessions.
   If not provided, biosample is pulled from the metadata TSV.

3. **Assembly FASTAs** — one per genome, in a directory. 

## Workflow

### Step 1: Generate Manifests

```bash
python 01_generate_manifests.py \
    --metadata genome-info.tsv \
    --fasta-dir /path/to/assemblies/ \
    --output-dir manifests/ \
    --study PRJNA_XXXXXXX \
    --dry-run  # preview first
```

```bash
python 01_generate_manifests_annotated.py \
    --metadata genome-info.tsv \
    --accessions accession_map.tsv \
    --fasta-dir /path/to/assemblies/ \
    --output-dir manifests/ \
    --study PRJNA_XXXXXXX
```

### Step 2: Validate

```bash
chmod +x 02_validate.sh
./02_validate.sh manifests/ Webin-XXXXX
```

### Step 3: Test Submit (3 assemblies to test server)

```bash
chmod +x 03_test_submit.sh
./03_test_submit.sh manifests/ Webin-XXXXX 3
```

Verify at: https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/

### Step 4: Production Submit

```bash
chmod +x 04_batch_submit.sh
./04_batch_submit.sh manifests/ Webin-XXXXX
```

### Step 5: Collect Accessions

```bash
# Immediately after submission (ERZ receipts only)
python 05_collect_accessions.py \
    --log-dir manifests/logs/submit/ \
    --output accessions.tsv

# after GCA accessioning.
python 05_collect_accessions.py \
    --log-dir manifests/logs/submit/ \
    --output accessions.tsv \
    --query-ena
```

## Column Name Configuration


```python
# Columns in the genome-info metadata TSV
META_ASSEMBLY_NAME = "assembly_name"
META_BIOSAMPLE = "biosample_accession"
META_COVERAGE = "genome_coverage"
META_ASSEMBLER = "assembly_methods"
META_PLATFORM = "sequencing_technologies"
META_FILENAME = "filename"

# Columns in the accession mapping file
ACC_ASSEMBLY_NAME = "assembly_name"
ACC_BIOPROJECT = "bioproject"
ACC_BIOSAMPLE = "biosample"
```

## Directory Structure After Running

```
manifests/
├── manifest_sample001.txt
├── manifest_sample002.txt
├── ...
├── manifest_summary.csv
├── fastas/
│   ├── sample001.fasta.gz
│   └── ...
├── logs/
│   ├── validate/
│   ├── test_submit/
│   └── submit/
├── submission_receipts_YYYYMMDD_HHMMSS.tsv
└── accessions.tsv
```
