#!/bin/bash
# Michael J. Foster
# 2025-JAN-15
# https://github.com/mjfos2r

set -e  # Exit on error
set -u  # Exit on undefined variable

# Default values
INPUT_FILE=""
OUTPUT_DIR=""
PREFIX=""
VERSION=""
ID_MAPPING=""
MINLEN=1000
MINWEIGHT=0.5

function usage() {
    cat << EOF
Usage: $0 [-h] -i <input> -o <output> -p <prefix> -v <version> -m <mapping>

Pipeline runner for sequence homology network analysis.

Required arguments:
    -i <input>    Input multifasta file path
    -o <output>   Output directory path
    -p <prefix>   Analysis prefix
    -v <version>  Analysis version
    -m <mapping>  ID mapping file path
    -l <int>      minimum length for edge filtering
    -w <float>    minimum weight for edge filtering

Optional arguments:
    -h           Show this help message

Example:
    $0 -i sequences.fasta -o results -p myanalysis -v v1
EOF
}

function check_dependencies() {
    # Check for required tools
    local missing_deps=()
    
    if ! command -v docker &> /dev/null; then
        missing_deps+=("docker")
    fi
    
    if ! command -v python &> /dev/null; then
        missing_deps+=("python")
    fi
    
    if ((${#missing_deps[@]} > 0)); then
        echo "Error: Missing required dependencies: ${missing_deps[*]}" >&2
        exit 1
    fi
}

function validate_inputs() {
    # Check if input file exists
    if [[ ! -f "$INPUT_FILE" ]]; then
        echo "Error: Input file does not exist: $INPUT_FILE" >&2
        exit 1
    fi
    
    # Check if input file is readable
    if [[ ! -r "$INPUT_FILE" ]]; then
        echo "Error: Input file is not readable: $INPUT_FILE" >&2
        exit 1
    fi
    
    # Create output directory if it doesn't exist
    if [[ ! -d "$OUTPUT_DIR" ]]; then
        mkdir -p "$OUTPUT_DIR" || {
            echo "Error: Could not create output directory: $OUTPUT_DIR" >&2
            exit 1
        }
    fi
    
    # Check if output directory is writable
    if [[ ! -w "$OUTPUT_DIR" ]]; then
        echo "Error: Output directory is not writable: $OUTPUT_DIR" >&2
        exit 1
    fi
}

function setup_directories() {
    # Create directory structure
    ALIGNMENTS_DIR="${OUTPUT_DIR}/alignments/ava/${VERSION}"
    NETWORK_DIR="${OUTPUT_DIR}/homology_networks/${VERSION}"
    PLOTS_DIR="${NETWORK_DIR}/plots"
    
    mkdir -p "$ALIGNMENTS_DIR" "$NETWORK_DIR" "$PLOTS_DIR" || {
        echo "Error: Failed to create directory structure" >&2
        exit 1
    }
}

function run_alignments() {
    echo "Running all-vs-all alignment with MUMmer4..."
    docker run -v "$(pwd):/data" -it mjfos2r/mummer4 \
        -run \
        "/data/${INPUT_FILE}" \
        "/data/${ALIGNMENTS_DIR}/${PREFIX}" || {
        echo "Error: MUMmer4 alignment failed" >&2
        exit 1
    }
}

function parse_alignments() {
    echo "Parsing MUMmer output..."
    python scripts/parse_mummer_output.py \
        "${ALIGNMENTS_DIR}/${PREFIX}.coords.tab" \
        "${NETWORK_DIR}" || {
        echo "Error: Failed to parse alignments" >&2
        exit 1
    }
}

function generate_matrix() {
    echo "Generating homology matrix..."
    python scripts/make_homology_matrix.py \
        "${NETWORK_DIR}/${PREFIX}_homology_network.tsv" \
        "${INPUT_FILE}" \
        "${NETWORK_DIR}" || {
        echo "Error: Failed to generate homology matrix" >&2
        exit 1
    }
}

function create_visualizations() {
    echo "Generating 3D visualizations..."
    python scripts/draw_3D_graph.py \
        --matrix "${NETWORK_DIR}/${PREFIX}_homology_network_alignment_matrix.tsv" \
        --edge-matrix "${NETWORK_DIR}/${PREFIX}_matrix_edges.csv" \
        --id-mapping "$ID_MAPPING" \
        --output-dir "${PLOTS_DIR}" \
        --layout "kk3d" \
        --min-contig-length $MINLEN \
        --min-weight $MINWEIGHT || {
        echo "Error: Failed to generate visualizations" >&2
        exit 1
    }
}

# Parse command line arguments
while getopts "hi:o:p:v:m:l:w:" opt; do
    case $opt in
        h)
            usage
            exit 0
            ;;
        i)
            INPUT_FILE="$OPTARG"
            ;;
        o)
            OUTPUT_DIR="$OPTARG"
            ;;
        p)
            PREFIX="$OPTARG"
            ;;
        v)
            VERSION="$OPTARG"
            ;;
        m)
            ID_MAPPING="$OPTARG"
            ;;
        l)
            MINLEN="$OPTARG"
            ;;
        w)
            MINWEIGHT="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument" >&2
            usage
            exit 1
            ;;
    esac
done

# Check if all required arguments are provided
if [[ -z "$INPUT_FILE" || -z "$OUTPUT_DIR" || -z "$PREFIX" || -z "$VERSION" || -z "$ID_MAPPING" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
    exit 1
fi

# Main pipeline execution
echo "Starting pipeline with:"
echo "  Input file: $INPUT_FILE"
echo "  Output directory: $OUTPUT_DIR"
echo "  Prefix: $PREFIX"
echo "  Version: $VERSION"
echo "  ID mapping file: $ID_MAPPING"
echo "  Minimum Length: $MINLEN"
echo "  Minimum Weight: $MINWEIGHT"

check_dependencies
validate_inputs
setup_directories

# Run pipeline steps
echo "Running pipeline steps..."
run_alignments
parse_alignments
generate_matrix
create_visualizations

echo "Pipeline completed successfully!"
echo "Results can be found in: $OUTPUT_DIR"