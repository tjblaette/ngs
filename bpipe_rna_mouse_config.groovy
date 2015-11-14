//file to contain pipeline configuration of RNA-seq pipeline for mouse


// VARIANT FILTER SCRIPT
//FILTER="filter_exome_wFlankingSeq.sh"

// KNOWN SEQUENCES
REF="${NGS}/refgenome/mm10/mm10_reordered.fa"
REFDIR="${NGS}/refgenome/mm10/starIndexFixedGTF/" //ref indexed by STAR with annotation
GTF="${NGS}/known_sites/mm10/gencode.vM7.chr_patch_hapl_scaff.annotation_UCSCcontigs.gtf"

// TARGET SEQUENCES FOR COVERAGE CALC
MOUSE="${NGS}/refgenome/mm10/mouse_wes_padded_mm9lifted2mm10_galaxy_sorted.bed"
EXON_TARGET=MOUSE
