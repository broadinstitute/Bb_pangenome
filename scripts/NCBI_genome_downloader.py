#!/home/mf019/bin/miniconda3/bin/python
from datetime import datetime
from Bio import Entrez
from pathlib import Path
import time
import os
import logging
from typing import List, Dict, Optional
import argparse

class NCBIGenomeFetcher:
    def __init__(self, email: str, output_dir, api_key: Optional[str] = None):
        """
        Initialize the NCBI genome fetcher with your email and optional API key.

        Args:
            email (str): Your email address (required by NCBI)
            api_key (str, optional): Your NCBI API key for higher rate limits
        """
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

        # Setup logging
        timestamp = datetime.now().strftime('__%Y%m%d__%H%M%S__')
        logfile = f'NCBI_genome_downloader_{timestamp}.log'
        logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s - %(levelname)s - %(message)s',
            filename=logfile
        )
        self.logger = logging.getLogger(__name__)
        self.logger.debug("creating output directory")
        # Create output directory
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.logger.debug("output directory created!")

    def search_assemblies(self, taxonomy_id: int = 139) -> List[Dict]:
        """
        Search for B. burgdorferi assemblies with non-Illumina sequencing.

        Args:
            taxonomy_id (int): Taxonomy ID for the organism

        Returns:
            List[Dict]: List of assembly information dictionaries
        """
        query = f"""
            txid{taxonomy_id}[Organism:exp] AND
            (("PacBio SMRT"[Properties] OR "Oxford Nanopore"[Properties])) AND
            (latest[filter] AND ("chromosome level"[filter] OR "complete genome"[filter]))
        """
        self.logger.debug("Using Query: %s", query)
        try:
            # Search assemblies
            self.logger.info("Searching for matching assemblies...")
            handle = Entrez.esearch(db="assembly", term=query, retmax=1000)
            record = Entrez.read(handle)
            handle.close()

            if not record["IdList"]:
                self.logger.warning("No matching assemblies found")
                return []

            # Get assembly details
            assembly_list = []
            for assembly_id in record["IdList"]:
                time.sleep(0.5)  # Rate limiting
                handle = Entrez.esummary(db="assembly", id=assembly_id)
                summary = Entrez.read(handle)
                handle.close()

                assembly_list.append({
                    'accession': summary['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession'],
                    'organism': summary['DocumentSummarySet']['DocumentSummary'][0]['Organism'],
                    'ftp_path': summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq'] or
                               summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
                })

            return assembly_list

        except Exception as e:
            self.logger.error(f"Error searching assemblies: {str(e)}")
            return []

    def download_genome(self, assembly_info: Dict, **kwargs) -> None:
        """
        Download genome sequence and annotations for a given assembly.

        Args:
            assembly_info (Dict): Assembly information dictionary
        """
        accession = assembly_info['accession']
        assembly_dir = self.output_dir / accession
        assembly_dir.mkdir(exist_ok=True)

        try:
            # Try RefSeq first
            query = f"{accession}[Assembly] AND refseq[filter]"
            handle = Entrez.esearch(db="nucleotide", term=query)
            record = Entrez.read(handle)
            handle.close()

            if record["Count"] == "0":
                # If no RefSeq, try GenBank
                query = f"{accession}[Assembly]"
                handle = Entrez.esearch(db="nucleotide", term=query)
                record = Entrez.read(handle)
                handle.close()

            if record["Count"] == "0":
                self.logger.warning(f"No sequences found for {accession}")
                return

            # Download sequences and annotations
            for seq_id in record["IdList"]:
                time.sleep(0.5)  # Rate limiting

                # Download FASTA
                self.logger.info(f"Downloading FASTA for {accession}")
                if kwargs.get('dryrun') is True:
                    handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
                    with open(assembly_dir / f"{accession}_genomic.fasta", "w") as f:
                        f.write(handle.read())
                    handle.close()

                    # Download GenBank
                    time.sleep(0.5)  # Rate limiting
                    self.logger.info(f"Downloading GenBank file for {accession}")
                    handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
                    with open(assembly_dir / f"{accession}_annotation.gbk", "w") as f:
                        f.write(handle.read())
                    handle.close()
                else:
                    print("Dryrun enabled! Iterating to next accession!")

        except Exception as e:
            self.logger.error(f"Error downloading {accession}: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description='a simple tool to query NCBI and download genomes')
    parser.add_argument('--dryrun', required=False, action='store_true', help='perform a dry run before downloading?')
    parser.add_argument('--output_dir', '--o', required=True, help='specify the directory to output downloaded genomes')
    
    args = parser.parse_args()

    # Initialize fetcher with your email
    fetcher = NCBIGenomeFetcher(email="SANITIZED FOR REPO", output_dir=args.output_dir)  # Replace with your email

    # Search for assemblies
    assemblies = fetcher.search_assemblies()

    # Download each assembly
    for assembly in assemblies:
        fetcher.download_genome(assembly, dryrun=args.dryrun)

if __name__ == "__main__":
    main()
