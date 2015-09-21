//file to contain pipeline configuration of WES pipeline


// VARIANT FILTER SCRIPT
FILTER="filter_exome_wFlankingSeq.sh"

// KNOWN SEQUENCES
REF="${NGS}/varsim_run/hs37d5.fa" //also includes unplaced/unlocalized contigs and alternative haplotypes
GOLD_STANDARD_1000G_INDELS="${NGS}/known_sites/b37/Mills_and_1000G_gold_standard.indels.b37.vcf"
PHASE1_1000G_INDELS="${NGS}/known_sites/b37/1000G_phase1.indels.b37.vcf"

// TARGET SEQUENCES FOR COVERAGE CALC

