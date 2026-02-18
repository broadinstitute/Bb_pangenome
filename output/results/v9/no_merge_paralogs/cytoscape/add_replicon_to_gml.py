#!/usr/bin/env python3
"""
Add replicon identity to a Panaroo pangenome GML file.

Uses gene_data.csv to build a mapping from (isolateIdx, scaffoldIdx) to
replicon name, then annotates each GML node by parsing its geneIDs.

Supports a best_hits lookup table to map NCBI accession numbers to
standardized replicon names (e.g., cp019844.1 -> chromosome).

Usage:
    python add_replicon_to_gml.py \
        --gene-data gene_data.csv \
        --gml final_graph.gml \
        --output final_graph_replicon.gml \
        --threshold 0.9 \
        --best-hits best_hits_1000bp_v10_2.csv
"""

import argparse
import csv
import logging
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime

import networkx as nx


def setup_logging(log_file: str = None) -> logging.Logger:
    """Setup logging to both console and file."""
    logger = logging.getLogger("replicon_annotator")
    logger.setLevel(logging.DEBUG)
    
    # Console handler (INFO level)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_format = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)
    
    # File handler (DEBUG level - more detailed)
    if log_file:
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setLevel(logging.DEBUG)
        file_format = logging.Formatter('%(asctime)s | %(levelname)s | %(message)s')
        file_handler.setFormatter(file_format)
        logger.addHandler(file_handler)
    
    return logger


def build_accession_lookup(best_hits_path: str, logger: logging.Logger) -> dict:
    """
    Build a lookup from NCBI accession numbers to standardized replicon names.
    
    best_hits CSV has columns: assembly_id, contig_id, ..., plasmid_name
    contig_id format: "NZ_CP019844.1 Borreliella burgdorferi strain PAli chromosome..."
    
    We extract the accession (CP019844.1), lowercase it, and map to plasmid_name.
    """
    accession_lookup = {}
    rows_parsed = 0
    
    with open(best_hits_path, "r") as f:
        reader = csv.DictReader(f)
        
        if "contig_id" not in reader.fieldnames or "plasmid_name" not in reader.fieldnames:
            logger.warning(f"best_hits file missing required columns (contig_id, plasmid_name)")
            return {}
        
        for row in reader:
            rows_parsed += 1
            contig_id = row.get("contig_id", "")
            plasmid_name = row.get("plasmid_name", "")
            
            if not contig_id or not plasmid_name:
                continue
            
            # Extract accession from contig_id
            # Format: "NZ_CP019844.1 Borreliella burgdorferi..." or just the accession
            # Also handles: "contig_1 [gcode=11]..." for non-NCBI
            
            # Try to find accession pattern (letters + numbers + .version)
            match = re.search(r'([A-Z]{2}_)?([A-Z]{2}\d+\.\d+)', contig_id)
            if match:
                # Get just the core accession (e.g., CP019844.1), lowercase
                accession = match.group(2).lower()
                if accession not in accession_lookup:
                    accession_lookup[accession] = plasmid_name.lower()
                    logger.debug(f"Accession mapping: {accession} -> {plasmid_name.lower()}")
    
    logger.info(f"  Parsed {rows_parsed:,} rows from best_hits")
    logger.info(f"  Built {len(accession_lookup):,} accession -> replicon mappings")
    
    if accession_lookup:
        logger.info(f"\n  Example accession mappings:")
        for i, (acc, rep) in enumerate(accession_lookup.items()):
            if i >= 5:
                break
            logger.info(f"    {acc} -> {rep}")
    
    return accession_lookup


def parse_replicon_from_scaffold(scaffold_name: str, accession_lookup: dict = None) -> str:
    """
    Extract replicon identity from scaffold name.

    Expected: B331P_chromosome_contig_1 -> chromosome
              B331P_lp54_contig_1      -> lp54
              B331P_cp32-3_contig_1    -> cp32-3
              GCF_xxx_cp019844.1_contig_1 -> lookup in accession table -> chromosome
    
    If the replicon looks like an accession (cpXXXXXX.1 pattern) and we have
    a lookup table, resolve it to the standardized name.
    """
    if accession_lookup is None:
        accession_lookup = {}
    
    parts = scaffold_name.split("_")
    contig_indices = [i for i, p in enumerate(parts) if p.lower() == "contig"]

    if contig_indices:
        contig_idx = contig_indices[0]
        if contig_idx > 1:
            replicon = "_".join(parts[1:contig_idx]).lower()
        else:
            return "unknown"
    else:
        if len(parts) > 1:
            replicon = "_".join(parts[1:]).lower()
        else:
            return "unknown"
    
    # Check if this looks like an NCBI accession (e.g., cp019844.1)
    # Pattern: 2 letters + 6 digits + . + version number
    if re.match(r'^[a-z]{2}\d{6}\.\d+$', replicon):
        if replicon in accession_lookup:
            return accession_lookup[replicon]
        # If not in lookup, keep the accession as-is (will show up in output)
    
    return replicon


def build_scaffold_lookup(gene_data_path: str, accession_lookup: dict = None, logger: logging.Logger = None) -> dict:
    """
    Parse gene_data.csv and build a lookup from (isolateIdx, scaffoldIdx)
    to replicon name.

    clustering_id format: isolateIdx_scaffoldIdx_geneIdx
    scaffold_name format: ISOLATE_REPLICON_contig_N
    """
    if accession_lookup is None:
        accession_lookup = {}
    
    scaffold_lookup = {}
    rows_parsed = 0
    skipped = 0
    accession_resolved = 0

    with open(gene_data_path, "r") as f:
        header_line = f.readline()
        delimiter = "\t" if "\t" in header_line else ","

    with open(gene_data_path, "r") as f:
        reader = csv.DictReader(f, delimiter=delimiter)

        expected = {"scaffold_name", "clustering_id"}
        if not expected.issubset(set(reader.fieldnames or [])):
            alt_delim = "," if delimiter == "\t" else "\t"
            f.seek(0)
            reader = csv.DictReader(f, delimiter=alt_delim)
            if not expected.issubset(set(reader.fieldnames or [])):
                logger.error(f"gene_data.csv must contain columns: {expected}")
                logger.error(f"Found columns: {reader.fieldnames}")
                sys.exit(1)

        for row in reader:
            rows_parsed += 1
            scaffold = row.get("scaffold_name", "")
            clustering_id = row.get("clustering_id", "")

            if not scaffold or not clustering_id:
                skipped += 1
                continue

            parts = clustering_id.split("_")
            if len(parts) >= 3:
                isolate_idx = parts[0]
                scaffold_idx = parts[1]
                key = (isolate_idx, scaffold_idx)
                if key not in scaffold_lookup:
                    replicon = parse_replicon_from_scaffold(scaffold, accession_lookup)
                    scaffold_lookup[key] = replicon
                    logger.debug(f"Scaffold mapping: {key} ({scaffold}) -> {replicon}")
                    # Track if this was resolved via accession lookup
                    if re.match(r'^[a-z]{2}\d{6}\.\d+$', "_".join(scaffold.split("_")[1:]).lower().split("_contig")[0] if "_contig" in scaffold.lower() else ""):
                        accession_resolved += 1

    logger.info(f"  Parsed {rows_parsed:,} rows from gene_data.csv")
    logger.info(f"  Skipped {skipped:,} rows")
    logger.info(f"  Built {len(scaffold_lookup):,} unique (isolate, scaffold) -> replicon mappings")

    logger.info(f"\n  Example mappings:")
    for i, (key, val) in enumerate(scaffold_lookup.items()):
        if i >= 8:
            break
        logger.info(f"    ({key[0]}, {key[1]}) -> {val}")

    return scaffold_lookup


def assign_consensus_replicon(replicon_list: list, threshold: float = 0.9) -> dict:
    """Determine consensus replicon for a gene cluster."""
    counts = Counter(replicon_list)
    total = len(replicon_list)

    if total == 0:
        return {
            "consensus_replicon": "unknown",
            "top_replicon": "unknown",
            "consensus_fraction": 0.0,
            "n_isolates": 0,
            "all_replicons": "unknown",
        }

    top_replicon, top_count = counts.most_common(1)[0]
    fraction = top_count / total
    consensus = top_replicon if fraction >= threshold else "multi-replicon"
    all_replicons = ",".join(f"{rep}({cnt})" for rep, cnt in counts.most_common())

    return {
        "consensus_replicon": consensus,
        "top_replicon": top_replicon,
        "consensus_fraction": round(fraction, 3),
        "n_isolates": total,
        "all_replicons": all_replicons,
    }


def classify_replicon_type(replicon: str) -> str:
    """Classify into broad categories for coarse coloring."""
    if replicon == "multi-replicon":
        return "multi-replicon"
    elif replicon == "chromosome":
        return "chromosome"
    elif replicon.startswith("cp"):
        return "circular_plasmid"
    elif replicon.startswith("lp"):
        return "linear_plasmid"
    elif replicon in ("unknown", "unmatched"):
        return replicon
    else:
        return "other"


def classify_replicon_family(replicon: str) -> str:
    """
    Classify replicon into family groups.
    
    lp28-1, lp28-2, lp28-3, etc. -> lp28
    cp32-1, cp32-3, cp32-1+5, etc. -> cp32
    chromosome -> chromosome
    """
    if replicon in ("unknown", "unmatched", "multi-replicon", "unknownreplicon"):
        return replicon
    
    # Handle chromosome
    if replicon == "chromosome":
        return "chromosome"
    
    # Handle fusion plasmids - classify by first component
    if "+" in replicon or "-cp" in replicon or "-lp" in replicon:
        # e.g., cp32-1+5 -> cp32, lp21-cp9 -> lp21
        base = replicon.split("+")[0].split("-cp")[0].split("-lp")[0]
        # Now extract family from base
        match = re.match(r'^([a-z]+\d+)', base)
        if match:
            return match.group(1)
    
    # Standard pattern: lp28-1 -> lp28, cp32-3 -> cp32, cp26 -> cp26
    match = re.match(r'^([a-z]+\d+)', replicon)
    if match:
        return match.group(1)
    
    # Fallback - return as-is
    return replicon


def get_replicon_color(replicon: str) -> str:
    """Return hex color for a replicon (distinct per-plasmid scheme)."""
    color_map = {
        # Chromosome
        "chromosome": "#bceddb",
        # cp26
        "cp26": "#d60000",
        # cp32 family - each distinct
        "cp32-1": "#ff28fd",
        "cp32-2": "#f2cdff",
        "cp32-3": "#d3008c",
        "cp32-4": "#c86e66",
        "cp32-5": "#ff6200",
        "cp32-6": "#93ac83",
        "cp32-7": "#ff8ec8",
        "cp32-8": "#b8ba01",
        "cp32-9": "#afa5ff",
        "cp32-10": "#953f1f",
        "cp32-11": "#9a6900",
        "cp32-12": "#9ee2ff",
        "cp32-13": "#a0e491",
        "cp32-1+5": "#fdf490",
        "cp32-3+10": "#72b8ff",
        "cp32-9-4": "#bc9157",
        # cp9 variants
        "cp9": "#56642a",
        "cp9-3": "#ae083f",
        # lp54
        "lp54": "#018700",
        # lp28 family - each distinct
        "lp28-1": "#366962",
        "lp28-2": "#f4bfb1",
        "lp28-3": "#97ff00",
        "lp28-4": "#79525e",
        "lp28-5": "#009e7c",
        "lp28-6": "#a877ac",
        "lp28-7": "#90318e",
        "lp28-8": "#ff3464",
        "lp28-9": "#e48eff",
        "lp28-11": "#c6a5c1",
        # Other linear plasmids
        "lp17": "#b500ff",
        "lp25": "#ffa52f",
        "lp36": "#05acc6",
        "lp38": "#00fdcf",
        "lp56": "#00c846",
        "lp5": "#8c9ab1",
        "lp21": "#829026",
        "lp21-cp9": "#77c6ba",
        "lp32-3": "#ff9070",
        # Unknown/Unclassified
        "unknown": "#d3c37c",
        "unmatched": "#FDFEFE",
        "unknownreplicon": "#d3c37c",
    }
    return color_map.get(replicon, "#d3c37c")  # Default to Unclassified color


def get_family_color(family: str) -> str:
    """Return hex color for a replicon family."""
    family_color_map = {
        "chromosome": "#48C9B0",
        "cp26": "#E67E22",
        "cp32": "#E74C3C",  # Red for cp32 family
        "cp9": "#F1C40F",
        "lp54": "#1E8449",
        "lp28": "#8E44AD",  # Purple for lp28 family
        "lp17": "#2ECC71",
        "lp25": "#27AE60",
        "lp36": "#16A085",
        "lp38": "#1ABC9C",
        "lp56": "#8B4513",
        "lp5": "#82E0AA",
        "lp21": "#58D68D",
        "lp32": "#73C6B6",
    }
    return family_color_map.get(family, "#95A5A6")


def generate_pie_chart_string(replicon_detail: str) -> str:
    """
    Generate Enhanced Graphics pie chart string from replicon_detail.
    
    Input: "cp32-3(49),cp32-10(5),lp56(2)"
    Output: "piechart: attributelist="cp32-3,cp32-10,lp56" colorlist="#F1948A,#E53935,#8B4513" valuelist="49,5,2""
    
    This format is used by Cytoscape's Enhanced Graphics plugin.
    """
    if not replicon_detail or replicon_detail == "unknown":
        return ""
    
    # Parse the detail string
    parts = replicon_detail.split(",")
    replicons = []
    values = []
    colors = []
    
    for part in parts:
        # Parse "replicon(count)" format
        if "(" in part and ")" in part:
            replicon = part.split("(")[0].strip()
            count = part.split("(")[1].replace(")", "").strip()
            replicons.append(replicon)
            values.append(count)
            colors.append(get_replicon_color(replicon))
    
    if not replicons:
        return ""
    
    # Build Enhanced Graphics pie chart string
    pie_str = f'piechart: attributelist="{",".join(replicons)}" colorlist="{",".join(colors)}" valuelist="{",".join(values)}" showlabels="false"'
    
    return pie_str


def analyze_family_diversity(replicon_list: list) -> dict:
    """
    Analyze family-level diversity for a gene cluster.
    
    Returns metrics about cross-family vs within-family variation.
    """
    from collections import Counter
    
    if not replicon_list:
        return {
            "families": [],
            "family_counts": {},
            "n_families": 0,
            "top_family": "unknown",
            "top_family_count": 0,
            "family_consensus_frac": 0.0,
            "is_single_family": True,
            "cross_family_score": 0.0,
        }
    
    # Get family for each replicon
    families = [classify_replicon_family(r) for r in replicon_list]
    family_counts = Counter(families)
    n_families = len(family_counts)
    total = len(families)
    
    top_family, top_count = family_counts.most_common(1)[0]
    family_consensus_frac = top_count / total if total > 0 else 0
    
    is_single_family = n_families == 1
    
    # Cross-family score: 0 = single family, 1 = maximally diverse
    # Score increases with number of families and evenness of distribution
    if n_families <= 1:
        cross_family_score = 0.0
    else:
        # Shannon entropy normalized by max possible entropy
        import math
        entropy = -sum((c/total) * math.log2(c/total) for c in family_counts.values() if c > 0)
        max_entropy = math.log2(n_families) if n_families > 1 else 1
        cross_family_score = entropy / max_entropy if max_entropy > 0 else 0
    
    return {
        "families": list(family_counts.keys()),
        "family_counts": dict(family_counts),
        "n_families": n_families,
        "top_family": top_family,
        "top_family_count": top_count,
        "family_consensus_frac": round(family_consensus_frac, 3),
        "is_single_family": is_single_family,
        "cross_family_score": round(cross_family_score, 3),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Add replicon identity to Panaroo pangenome GML"
    )
    parser.add_argument("--gene-data", required=True, help="Panaroo gene_data.csv")
    parser.add_argument("--gml", required=True, help="Input GML file")
    parser.add_argument("--output", default="final_graph_replicon.gml", help="Output GML")
    parser.add_argument("--threshold", type=float, default=0.9, help="Consensus threshold")
    parser.add_argument("--summary", default="replicon_summary.tsv", help="Summary TSV")
    parser.add_argument("--best-hits", default=None, help="Best hits CSV for accession->replicon lookup")
    parser.add_argument("--log", default="replicon_annotation.log", help="Log file path")
    args = parser.parse_args()

    # Setup logging
    logger = setup_logging(args.log)
    logger.info(f"=" * 60)
    logger.info(f"Replicon Annotation Tool")
    logger.info(f"Run started: {datetime.now().isoformat()}")
    logger.info(f"=" * 60)
    logger.info(f"\nParameters:")
    logger.info(f"  gene-data:  {args.gene_data}")
    logger.info(f"  gml:        {args.gml}")
    logger.info(f"  output:     {args.output}")
    logger.info(f"  threshold:  {args.threshold}")
    logger.info(f"  summary:    {args.summary}")
    logger.info(f"  best-hits:  {args.best_hits}")
    logger.info(f"  log:        {args.log}")

    # Step 0: Build accession lookup if provided
    accession_lookup = {}
    if args.best_hits:
        logger.info(f"\nBuilding accession lookup from: {args.best_hits}")
        accession_lookup = build_accession_lookup(args.best_hits, logger)

    # Step 1: Build scaffold lookup
    logger.info(f"\nBuilding scaffold lookup from: {args.gene_data}")
    scaffold_lookup = build_scaffold_lookup(args.gene_data, accession_lookup, logger)

    unique_replicons = sorted(set(scaffold_lookup.values()))
    logger.info(f"\n  Unique replicons found ({len(unique_replicons)}):")
    for r in unique_replicons:
        logger.info(f"    {r}")

    # Step 2: Load GML
    logger.info(f"\nReading GML from: {args.gml}")
    G = nx.read_gml(args.gml)
    logger.info(f"  {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Step 3: Annotate nodes via geneIDs
    logger.info(f"\nAnnotating nodes (threshold={args.threshold})...")
    matched = 0
    unmatched = 0
    partially_matched = 0
    refound_count = 0

    for node_id, attrs in G.nodes(data=True):
        gene_ids_str = attrs.get("geneIDs", "")
        if not gene_ids_str:
            G.nodes[node_id]["replicon"] = "no_geneIDs"
            G.nodes[node_id]["top_replicon"] = "no_geneIDs"
            G.nodes[node_id]["replicon_type"] = "unknown"
            G.nodes[node_id]["replicon_detail"] = ""
            G.nodes[node_id]["replicon_consensus_frac"] = 0.0
            G.nodes[node_id]["replicon_n_isolates"] = 0
            unmatched += 1
            logger.debug(f"Node {node_id}: no geneIDs found")
            continue

        gene_ids = gene_ids_str.split(";")
        replicons = []
        refound_in_node = 0

        for gid in gene_ids:
            parts = gid.strip().split("_")
            if len(parts) >= 3:
                isolate_idx = parts[0]
                scaffold_idx = parts[1]

                if scaffold_idx == "refound":
                    refound_in_node += 1
                    refound_count += 1
                    continue

                key = (isolate_idx, scaffold_idx)
                if key in scaffold_lookup:
                    replicons.append(scaffold_lookup[key])

        if replicons:
            result = assign_consensus_replicon(replicons, args.threshold)
            family_result = analyze_family_diversity(replicons)
            
            G.nodes[node_id]["replicon"] = result["consensus_replicon"]
            G.nodes[node_id]["top_replicon"] = result["top_replicon"]
            G.nodes[node_id]["replicon_type"] = classify_replicon_type(
                result["consensus_replicon"]
            )
            G.nodes[node_id]["top_replicon_type"] = classify_replicon_type(
                result["top_replicon"]
            )
            G.nodes[node_id]["replicon_detail"] = result["all_replicons"]
            G.nodes[node_id]["replicon_consensus_frac"] = result["consensus_fraction"]
            G.nodes[node_id]["replicon_n_isolates"] = result["n_isolates"]
            
            # Family-level metrics
            G.nodes[node_id]["top_family"] = family_result["top_family"]
            G.nodes[node_id]["n_families"] = family_result["n_families"]
            G.nodes[node_id]["family_consensus_frac"] = family_result["family_consensus_frac"]
            G.nodes[node_id]["is_single_family"] = 1 if family_result["is_single_family"] else 0
            G.nodes[node_id]["cross_family_score"] = family_result["cross_family_score"]
            G.nodes[node_id]["family_detail"] = ",".join(f"{f}({c})" for f, c in family_result["family_counts"].items())
            
            # Generate pie chart string for Enhanced Graphics
            pie_str = generate_pie_chart_string(result["all_replicons"])
            G.nodes[node_id]["pie_chart"] = pie_str
            
            # Store colors for styling
            G.nodes[node_id]["top_replicon_color"] = get_replicon_color(result["top_replicon"])
            G.nodes[node_id]["top_family_color"] = get_family_color(family_result["top_family"])

            gene_name = attrs.get("name", attrs.get("label", str(node_id)))
            logger.debug(f"Node {node_id} ({gene_name}): {result['consensus_replicon']} "
                        f"(frac={result['consensus_fraction']}, families={family_result['n_families']}, "
                        f"cross_family={family_result['cross_family_score']}, detail={result['all_replicons']})")

            if len(replicons) < len(gene_ids) - refound_in_node:
                partially_matched += 1
            matched += 1
        else:
            G.nodes[node_id]["replicon"] = "unmatched"
            G.nodes[node_id]["top_replicon"] = "unmatched"
            G.nodes[node_id]["replicon_type"] = "unmatched"
            G.nodes[node_id]["top_replicon_type"] = "unmatched"
            G.nodes[node_id]["replicon_detail"] = ""
            G.nodes[node_id]["replicon_consensus_frac"] = 0.0
            G.nodes[node_id]["replicon_n_isolates"] = 0
            G.nodes[node_id]["top_family"] = "unmatched"
            G.nodes[node_id]["n_families"] = 0
            G.nodes[node_id]["family_consensus_frac"] = 0.0
            G.nodes[node_id]["is_single_family"] = 1
            G.nodes[node_id]["cross_family_score"] = 0.0
            G.nodes[node_id]["family_detail"] = ""
            G.nodes[node_id]["pie_chart"] = ""
            G.nodes[node_id]["top_replicon_color"] = "#FDFEFE"
            G.nodes[node_id]["top_family_color"] = "#FDFEFE"
            unmatched += 1
            logger.debug(f"Node {node_id}: no replicons matched")

    logger.info(f"\n  Results:")
    logger.info(f"    Matched:           {matched}")
    logger.info(f"    Unmatched:         {unmatched}")
    logger.info(f"    Partially matched: {partially_matched}")
    logger.info(f"    Refound genes:     {refound_count} (skipped for replicon assignment)")

    # Replicon distribution
    replicon_counts = Counter(
        attrs.get("replicon", "unknown") for _, attrs in G.nodes(data=True)
    )
    logger.info(f"\n  Replicon distribution ({len(replicon_counts)} categories):")
    for rep, count in sorted(replicon_counts.items(), key=lambda x: -x[1]):
        logger.info(f"    {rep:25s} {count:6d} nodes")

    # Replicon type distribution
    type_counts = Counter(
        attrs.get("replicon_type", "unknown") for _, attrs in G.nodes(data=True)
    )
    logger.info(f"\n  Replicon type distribution:")
    for rtype, count in sorted(type_counts.items(), key=lambda x: -x[1]):
        logger.info(f"    {rtype:25s} {count:6d} nodes")

    # Multi-replicon detail
    multi_nodes = [
        (nid, attrs) for nid, attrs in G.nodes(data=True)
        if attrs.get("replicon") == "multi-replicon"
    ]
    if multi_nodes:
        logger.info(f"\n  Multi-replicon nodes ({len(multi_nodes)}):")
        for nid, attrs in multi_nodes[:20]:
            name = attrs.get("name", attrs.get("label", str(nid)))
            detail = attrs.get("replicon_detail", "")
            logger.info(f"    {name:30s} {detail}")
        if len(multi_nodes) > 20:
            logger.info(f"    ... and {len(multi_nodes) - 20} more")
        
        # Log ALL multi-replicon nodes to the log file (DEBUG level)
        logger.debug(f"\n  FULL multi-replicon node list:")
        for nid, attrs in multi_nodes:
            name = attrs.get("name", attrs.get("label", str(nid)))
            detail = attrs.get("replicon_detail", "")
            logger.debug(f"    {name:30s} {detail}")

    # Family diversity analysis
    single_family_count = sum(1 for _, attrs in G.nodes(data=True) if attrs.get("is_single_family") == 1)
    cross_family_count = sum(1 for _, attrs in G.nodes(data=True) if attrs.get("n_families", 0) > 1)
    
    logger.info(f"\n  Family diversity analysis:")
    logger.info(f"    Single-family nodes:     {single_family_count} (within-family variation only)")
    logger.info(f"    Cross-family nodes:      {cross_family_count} (true inter-replicon movement)")
    
    # Distribution of cross-family scores
    high_cross_family = [(nid, attrs) for nid, attrs in G.nodes(data=True) 
                         if attrs.get("cross_family_score", 0) > 0.5]
    if high_cross_family:
        logger.info(f"\n  High cross-family score nodes (>0.5): {len(high_cross_family)}")
        for nid, attrs in sorted(high_cross_family, key=lambda x: -x[1].get("cross_family_score", 0))[:15]:
            name = attrs.get("name", attrs.get("label", str(nid)))
            score = attrs.get("cross_family_score", 0)
            family_detail = attrs.get("family_detail", "")
            logger.info(f"    {name:30s} score={score:.2f}  {family_detail}")
        if len(high_cross_family) > 15:
            logger.info(f"    ... and {len(high_cross_family) - 15} more")

    # Step 4: Write output
    logger.info(f"\nWriting annotated GML to: {args.output}")
    nx.write_gml(G, args.output)

    logger.info(f"Writing summary to: {args.summary}")
    with open(args.summary, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "node_id", "gene_name", "label", "size", "replicon", "top_replicon",
            "replicon_type", "top_replicon_type", "consensus_fraction", "n_isolates", 
            "replicon_detail", "pie_chart", "top_replicon_color"
        ])
        for node_id, attrs in G.nodes(data=True):
            writer.writerow([
                node_id,
                attrs.get("name", ""),
                attrs.get("label", ""),
                attrs.get("size", ""),
                attrs.get("replicon", ""),
                attrs.get("top_replicon", ""),
                attrs.get("replicon_type", ""),
                attrs.get("top_replicon_type", ""),
                attrs.get("replicon_consensus_frac", ""),
                attrs.get("replicon_n_isolates", ""),
                attrs.get("replicon_detail", ""),
                attrs.get("pie_chart", ""),
                attrs.get("top_replicon_color", ""),
            ])

    logger.info(f"\nLog file written to: {args.log}")
    logger.info("\n" + "=" * 60)
    logger.info("Done! In Cytoscape:")
    logger.info("  For pie charts: Image/Chart -> 'pie_chart' column (Enhanced Graphics)")
    logger.info("  Shape -> 'top_replicon_type' (triangle=linear, circle=circular, diamond=chromosome)")
    logger.info("  Tooltip -> 'replicon_detail' (full breakdown per node)")
    logger.info("=" * 60)
    logger.info(f"\nRun completed: {datetime.now().isoformat()}")


if __name__ == "__main__":
    main()
