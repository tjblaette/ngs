//file to contain pipeline configuration of MOUSE pipeline


// VARIANT FILTER SCRIPT
FILTER="filter_mouse_wFlankingSeq.sh"

// KNOWN SEQUENCES
REF="${NGS}/refgenome/mm10/mm10_reordered.fa"

// TARGET SEQUENCES FOR COVERAGE CALC
MOUSE="${NGS}/refgenome/mm10/mouse_wes_padded_mm9lifted2mm10_galaxy_sorted.bed"
EXON_TARGET=MOUSE

