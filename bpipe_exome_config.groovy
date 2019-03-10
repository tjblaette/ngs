//file to contain pipeline configuration of WES pipeline


// VARIANT FILTER SCRIPT
FILTER="filter_exome.sh"

// KNOWN SEQUENCES
REF="${NGS}/refgenome/GATK/ucsc.hg19.fasta" //also includes unplaced/unlocalized contigs and alternative haplotypes
CANDIDATES="${NGS}/known_sites/candidate_genes_aml.txt"
DBSNP="${NGS}/known_sites/hg19/dbsnp_138.hg19.vcf"
GOLD_STANDARD_1000G_INDELS="${NGS}/known_sites/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
PHASE1_1000G_INDELS="${NGS}/known_sites/hg19/1000G_phase1.indels.hg19.sites.vcf"

// TARGET SEQUENCES FOR COVERAGE CALC
TRUSEQ="${NGS}/known_sites/hg19/TruSeq-Exome-Targeted-Regions-BED-file_sorted.bed" 
NEXTERA="${NGS}/known_sites/hg19/nexterarapidcapture_exome_targetedregions_v1.2_sorted.bed"
NEXTERA_EXPANDED="${NGS}/known_sites/hg19/nexterarapidcapture_expandedexome_targetedregions_sorted.bed"
EXON_TARGET=TRUSEQ

