
exec > >(cat >> $outfile)
exec 1> >(tee -a $outfile >&1)
exec 2> >(tee -a $errfile >&2)
log()
{
    echo "$1"
}

cd $path_data

LOG_DIR="/Users/mh057/Documents/BbSensuStricto/fastas"
cd $LOG_DIR
for genome in $(ls *.fna | cut -f1 -d "."); 
do
    echo -e " *** formatting "$genome" using anvio ..." ;
    anvi-script-reformat-fasta -o "$genome"_clean.fna \
            --prefix "$genome" \
            --simplify-names \
            -r "$genome"_report.txt \
			--seq-type NT\
            "$genome".fna; 
    echo " ###############################################################"
    echo -e " *** Generating contig db for "$genome" ..." ;
    anvi-gen-contigs-database -f "$genome"_clean.fna \
                              -o "$genome".db \
                              --num-threads 4 \
                              -n "$genome";
    echo " ################## done for "$genome" ##########################"
done; 
# anvi-setup-ncbi-cogs --num-threads $thread  --cog-data-dir /Users/mh057/anvio/COG20

for genome in $(ls *.db | cut -f1 -d "."); 
do
	echo -e " ****** (1) running HMMs annotation for $genome ";
	anvi-run-hmms -c $genome.db \
		--num-threads $thread \
		-I Bacteria_71 ;
	echo -e " ****** (2) running NCBI COG 2020v annotation for $genome ";
	anvi-run-ncbi-cogs -c $genome.db \
	--cog-data-dir $db_cog20 \
	--num-threads $thread;
	echo -e " ****** (3) running PFAM annotation for $genome ";
	anvi-run-pfams -c $genome.db \
		--num-threads $thread \
		--pfam-data-dir $db_pfam;
    echo -e " ****** (4) running KEGG KOfams annotation for $genome ";
    anvi-run-kegg-kofams -c $genome.db \
        --num-threads $thread \
        -H 0.6 -E 1e-06;
	echo -e " ######################### done for "$genome" ############################ "
done

	anvi-script-gen-genomes-file --input-dir . -o external-genomes.txt;
	
#	anvi-display-contigs-stats *db --report-as-text -o anvi_report.txt;

#	anvi-estimate-genome-completeness -e external-genomes.txt > anvi_genome_completeness.txt

	anvi-gen-genomes-storage -e external-genomes.txt -o LD83-GENOMES.db
   

	anvi-pan-genome -g LD83-GENOMES.db --project-name LD83_pangenome\
			 --num-threads 5\
			 --minbit 0.7\
			 --mcl-inflation 10\
             --description text_LD94.txt\
             --enforce-hierarchical-clustering \
             --output-dir LD83_run; 
	
	anvi-get-sequences-for-gene-clusters -p LD83_run/LD83_pangenome-PAN.db\
                                     -g LD83-GENOMES.db\
                                     --min-num-genomes-gene-cluster-occurs 83 \
                                     --max-num-genes-from-each-genome 1 \
                                     --concatenate-gene-clusters \
									 --report-DNA-sequences\
                                     --output-file SGC_LD83.aln;
	
	anvi-compute-genome-similarity -e external-genomes.txt \
					-p LD83_run/LD83_pangenome-PAN.db \
					-T 4 \
					-o pyANI_LD83_run \
					--program pyANI \
                    --min-alignment-fraction 0.1 \
                    --significant-alignment-length 2500;

	anvi-display-pan -p LD83_run/LD83_pangenome-PAN.db -g LD83-GENOMES.db;

	anvi-script-add-default-collection -p LD83_run/LD83_pangenome-PAN.db;

	anvi-summarize -p LD83_run/LD83_pangenome-PAN.db\
					-g LD83-GENOMES.db\
                    -C DEFAULT\
					--report-DNA-sequences\
                	-o SUMMARY;
	
	anvi-get-sequences-for-gene-clusters -p LD83_run/LD83_pangenome-PAN.db\
                                     -g LD83-GENOMES.db\
                                     --min-num-genomes-gene-cluster-occurs 82 \
                                     --max-num-genes-from-each-genome 1 \
                                     --concatenate-gene-clusters \
                                     --output-file SGCs_LD83.faa;

	trimal -in SGCs_LD83.faa -gt 0.5 -out SGCs_LD83_gt50.faa  
	
	iqtree -s SGCs_LD83_gt25.faa \
            -nt 4 -m WAG -bb 1000 --seed 123\
            -asr --boot-trees --runs 3 -T 5; 

# generate file for the import of phylogenetic tree to anvio
echo -e "item_name\tdata_type\tdata_value" > SGCS_LD83_PhyloTree_layer.txt
# add the newick tree as an order in the file
echo -e "SCGs_Bayesian_Tree_Rooted\tnewick\t`cat SCGs_LD83_newick.new`" >> SGCS_LD83_PhyloTree_layer.txt
echo -e "SCGs_Bayesian_Tree\tnewick\t`cat SCGs_LD83_nogaps.aln.contree`" >> SGCS_LD83_PhyloTree_layer.txt
# import the Phylogenetic file as layers order
anvi-import-misc-data -p LD83_run/LD83_pangenome-PAN.db \
                      -t layer_orders SGCS_LD90_PhyloTree_layer.txt;   

anvi-display-pan -p LD83_run/LD83_pangenome-PAN.db -g LD83-GENOMES.db; 

echo -e "########## analysis done ########## ";

