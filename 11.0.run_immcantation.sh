# bind docker image
#singularity exec -B /data/NCBR/projects/NCBR-323/:/NCBR-323,/data/zhangy68:/zhangy68 -c /data/OpenOmics/SIFs/immcantation_suite_4.4.0.sif /bin/bash

# change directory to Immcantation

path=../multi-3rd-pass_finalMultiUsingHTODemuxAssignment/download_for_downstream_analysis/VDJ_B_files

for i in TRIAD-18 TRIAD-7 TRIAD-2 TRIAD-3 TRIAD-4 TRIAD-1 HV-1 HV-3 HV-2; \
do AssignGenes.py igblast -s ${path}/${i}_vdj_b_filtered_contig.fasta -b /usr/local/share/igblast --organism human --loci ig --format blast --outdir . --outname ${i}_BCR_data_sequences; \
MakeDb.py igblast -i ${i}_BCR_data_sequences_igblast.fmt7 -s ${path}/${i}_vdj_b_filtered_contig.fasta -r /usr/local/share/germlines/imgt/human/vdj/ --10x ${path}/filtered_contig_annotations/${i}.csv --extended --outname ${i}_BCR_data_sequences_igblast_db-pass.tsv;
done
