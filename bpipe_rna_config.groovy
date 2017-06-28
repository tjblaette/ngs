//file to contain pipeline configuration of WES pipeline


// VARIANT FILTER SCRIPT
//FILTER="filter_exome_wFlankingSeq.sh"

// KNOWN SEQUENCES
REF="${NGS}/refgenome/GATK/ucsc.hg19_noAltHaps.fasta" //also includes unplaced/unlocalized contigs but no alternative haplotypes
REFDIR="${NGS}/refgenome/GATK/starIndexFixedGTFnoAltHaps/" //ref indexed by STAR with annotation
DBSNP="${NGS}/known_sites/hg19/dbsnp_138.hg19.vcf"
GOLD_STANDARD_1000G_INDELS="${NGS}/known_sites/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
PHASE1_1000G_INDELS="${NGS}/known_sites/hg19/1000G_phase1.indels.hg19.sites.vcf"
GTF="${NGS}/known_sites/hg19/gencode.v19.chr_patch_hapl_scaff.annotation_UCSCcontigs_noAltHaps.gtf"

// TARGET SEQUENCES FOR COVERAGE CALC
TRUSEQ="${NGS}/known_sites/hg19/TruSeq-Exome-Targeted-Regions-BED-file.bed" 
NEXTERA="${NGS}/known_sites/hg19/nexterarapidcapture_exome_targetedregions_v1.2.bed"
NEXTERA_EXPANDED="${NGS}/known_sites/hg19/nexterarapidcapture_expandedexome_targetedregions.bed"
EXON_TARGET=TRUSEQ

