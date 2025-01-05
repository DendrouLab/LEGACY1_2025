# 01_AssignGenes
singularity exec --no-home -B /legacy/:/legacy/  \
/packages/immcantation/immcantation_suite-4.4.0.sif \
    AssignGenes.py igblast -s $DATA_DIR/all_contig.fasta -b /usr/local/share/igblast\
   --organism human --loci $VDJ_TYPE --format blast --outdir $OUT_DIR --outname $OUT_NAME --nproc 2

# 02_MakeDB
singularity exec --no-home -B /well/legacy/:/well/legacy/  \
/packages/immcantation/immcantation_suite-4.4.0.sif \
    MakeDb.py igblast -i $FMT7_DIR/$FILE_NAME -s $DATA_DIR/all_contig.fasta \
   -r /usr/local/share/germlines/imgt/human/vdj/ \
   --10x $DATA_DIR/all_contig_annotations.csv --extended --outdir $OUT_DIR --outname $OUT_NAME

   