# NGS pipelines &  Co

<div align="justify">
This repository contains the collection of next-generation (NGS) analysis scripts that have accumulated during my PhD at Ulm University. They include pipelines for mutation calling and gene expression analysis, mostly for Illumina data, and several downstream / helper scripts. All analyses were developed to study acute myeloid leukemia (AML) but can also be applied to other samples.  They have, in parts, been used in the context of the following publications:


* Cocciardi et al. "Clonal evolution patterns in acute myeloid leukemia with *NPM1* mutation". In: Nature Communications 10.1 (2019)
* Scheffold et al. "IGF1R as druggable target mediating PI3K-δ inhibitor resistance in a murine model of chronic lymphocytic leukemia". In: Blood (2019)
* Rücker et al. "Chromothripsis is linked to TP53 alteration, cell cycle impairment, and dismal outcome in acute myeloid leukemia with complex karyotype". In: haematologica 103.1 (2018)
* L'Abbate et al. "MYC-containing amplicons in acute myeloid leukemia: genomic structures, evolution, and transcriptional consequences". In: Leukemia (2018)
* Hirsch et al. "Circular RNAs of the nucleophosmin (NPM1) gene in acute myeloid leukemia". In: haematologica (2017)
* Thiel et al. "Heterodimerization of AML1/ETO with CBFβ is required for leukemogenesis but not for myeloproliferation". In: Leukemia (2017)


## Author / Support
T.J. Blätte \
tamara.blaette@charite.de


## Overview

#### Pipelines
Actual pipeline scripts are prefixed by `bpipe_` and implemented using
the same. They include three main types of files:

* `*stages` files, which define the individual pipeline steps
* `*config` files, which define some default parameters
* `*pipeline` files, which piece together configs and stages to define a specific pipeline

#### Requirements
Bpipe and all other third-party tools must be installed and correctly assigned to the
respective variables defined within the `*stages` file.
Correctly formatted / pre-processed references must be provided for alignment; these are defined within the `*config` files. Depending on the pipeline to be used, additional database files are required for filtering and annotation.
Required inputs to all pipelines are paired FASTQ files of one or more samples. Paired FASTQ files are recognized based on the substring `_R?` - this should therefore be the only capital `R` contained in the FASTQ name and while additional characters may follow, different samples should be distinguished by whatever comes before.

To run a given pipeline, provide these to the pipeline file:
```
bpipe run your-choice.pipeline \
    your-input_R1.fastq \
    your-input_R2.fastq
```


#### Optional inputs and arguments
Several parameters can be controlled via optional command line arguments. To use these, follow the `bpipe run` command with `-p ` and the argument assignment.
Common to all pipelines is for example `nkern` for the number of cores to use for parallelization. Note that these are the cores used **per sample** and must therefore be adjusted when multiple samples are to be processed in parallel:

```
bpipe run -p nkern=42 your-choice.pipeline \
    your-input_R1.fastq \
    your-input_R2.fastq
```


---


## Mutation calling

Four pipelines are available to call mutations from Illumina sequencing data:

* `bpipe_amplicon.pipeline` for naive amplicons
* `bpipe_haloplexHS.pipeline` for molecularly barcoded amplicons
* `bpipe_exome.pipeline` for targeted enrichment protocols
* `bpipe_rna_callVariants.pipeline` for RNA sequencing data

Primary outputs of interest for all of these are the CSV files containing the filtered lists of mutations. Intermediate alignments and statistics are saved as well, at nearly every step of the analysis.

### Amplicon pipeline
The amplicon pipeline performs the following steps: Quality assessment with the NGS QC Toolkit, adapter trimming via cutadapt, alignment to the reference via BWA-MEM, optical duplicate removal via Picard, local realignment via GATK, coverage calculation with BEDTools, variant calling using VarScan2, annotation via ANNOVAR and filtering using our in-house developed `filter_*` scripts.


Expected inputs are paired FASTQ files and a BED file defining the target intervals of the reference covered by the design. Variant calling and on-target coverage calculation, but not alignment, will be restricted to these intervals. FASTQ files must be passed as positional arguments on the command line; the BED file may be defined as a default in the respective config or passed as an optional argument:

```
bpipe run -p target_region=your-targets.BED bpipe_amplicon.pipeline \
    your-input_R1.fastq \
    your-input_R2.fastq
```

To simultaneously change the default number of cores to use for the analysis, run:
```
bpipe run -p target_region=your-targets.BED -p nkern=42 bpipe_amplicon.pipeline \
    your-input_R1.fastq \
    your-input_R2.fastq
```


To analyze multiple samples in parallel using the default BED defined in the config, do:
```
bpipe run bpipe_amplicon.pipeline \
    your-input_R1.fastq \
    your-input_R2.fastq \
    your-other-input_R1.fastq \
    your-other-input_R2.fastq \
    your-yet-another-input_R1.fastq \
    your-yet-another-input_R2.fastq
```

### HaloplexHS pipeline
The haloplexHS pipeline processes molecularly barcoded amplicons and performs the following steps: Quality assessment with the NGS QC Toolkit, alignment to the reference via BWA-MEM, optical and barcode duplicate removal via Picard, local realignment via GATK, coverage calculation with BEDTools, variant calling using VarScan2, annotation via ANNOVAR and filtering using our in-house developed `filter_*` scripts.


Expected inputs are paired FASTQ files containing the sequencing reads, a third FASTQ named `*_R3*.fastq` containing the molecular barcode of each read pair and a BED file defining the target amplicons of the design. Variant calling and on-target coverage calculation, but not alignment, will be restricted to these intervals. FASTQ files must be passed as positional arguments on the command line; the BED file may be defined as a default in the respective config or passed as an optional argument:

```
bpipe run -p target_region=your-targets.BED bpipe_haloplexHS.pipeline \
    your-input_R1.fastq \
    your-input_R2.fastq \
    your-inputs-molecular-barcodes_R3.fastq
```


### Exome / targeted enrichment pipeline
The exome / targeted enrichment pipeline processes Illumina NGS data from targeted enrichment based protocols, including whole-exome sequencing (WES) but also designs to enrich for other target regions or subsets of genes of interest. It performs the following steps: Quality assessment with the NGS QC Toolkit, alignment to the reference via BWA-MEM, optical and sequence duplicate removal via Picard, local realignment via GATK, coverage calculation with BEDTools, variant calling using VarScan2, annotation via ANNOVAR and filtering using our in-house developed `filter_*` scripts.


Expected inputs are paired FASTQ files of two matched samples, the second being a healthy / control sample used to distinguish germline, somatic and LOH mutations, and a BED file defining the regions targeted by the design. Variant calling and on-target coverage calculation, but not alignment, will be restricted to these intervals. FASTQ files must be passed as positional arguments on the command line; the BED file may be defined as a default in the respective config or passed as an optional argument:

```
bpipe run bpipe_exome.pipeline \
        your-tumor_R1.fastq  \
        your-tumor_R2.fastq  \
        your-tumors-healthy-control_R1.fastq  \
        your-tumors-healthy-control_R2.fastq
```


### RNA calling
The RNA mutation calling pipeline processes Illumina RNA sequencing but is currently still somewhat experimental, as we have so far only used it to check for aberrations in specific candidate genes or in combination with matched DNA sequencing data. It mimics the exome / targeted enrichment pipeline and performs quality assessment with the NGS QC Toolkit, alignment to the reference via STAR, optical and sequence duplicate removal via Picard, local realignment via GATK, coverage calculation with BEDTools, variant calling using VarScan2, annotation via ANNOVAR and filtering using our in-house developed `filter_*` scripts:

```
bpipe run bpipe_rna_callVariants.pipeline \
        your-input_R1.fastq \
        your-input_R2.fastq
```

---


## Gene expression analysis
There is one bpipe pipeline available for gene expression analysis, which performs quality assessment with the NGS QC Toolkit, reference alignment and gene-wise quantification via STAR, optical duplicate removal via Picard and requantification with HTSeq.
With the help of different downstream analyses, this can be used to study both linear and circular transcripts.

Expected inputs are the paired FASTQ files of a single sample:
```
bpipe run bpipe_rna.pipeline \
        your-input_R1.fastq \
        your-input_R2.fastq
```

Primary outputs of interest are the STAR and HTSeq read count files. Note that STAR writes hidden files by default (filename starts with `.`).

### Count format conversion
HTSeq counts can be directly input into DESeq2 for downstream analysis of differential expression.
Counts output by STAR have to be reformatted first.
The main difference is that while HTSeq generates a different counts file for each of the possible stranded / unstranded library preparation protocols, STAR writes all variants to one file.
To extract the desired counts into a 2-column file that can be processed further, use `star_counts_to_htseq_format.sh`.

The script takes as input the counts file to be converted and the column to extract. When no column is given, it defaults to column `2`, for unstranded designs. Columns 3 and 4 contain counts for the respective stranded protocols:

```
star_counts_to_htseq_format.sh \
        your-star-counts.tab \
        2
```


### Differential gene expression analysis

Heavily based on the DESeq2 manual, the `deseq2*` scripts test for differentially expressed genes.

##### Inputs
Inputs to the wrapper script, which then itself calls the respective R script, are:

* the design table
* the design formula to test for
* the reference level
* alpha
* the minimum log fold change to test for
* a database file for ensembl gene ID to gene symbol conversion

The design table is a TSV file specifying at least the columns `fileName`, `sampleName` and the condition to test on for differential expression. Each line describes one sample. Filenames must refer to the respective read count files, which are expected to reside within the `intermediate_files` subfolder, in the directory from which the script is run. These count files are those output by HTSeq or STAR, following the respective formatting in the latter case. All of the information provided in this table will be used to annotate the respective output plots.

The design formula is something like `~ condition` for an unpaired analysis or `~ patient + condition` for a paired analysis. Provided factors must be defined as columns in the design table.

The reference level defines relative to what fold changes are calculated. Typically, this is something like `untreated`, `control` or similar. The term must be specified in the tested condition column of the design table.

Alpha corresponds to the false discovery rate (FDR) cutoff that defines differentially expressed genes. Default is `0.1`.

The minimum log fold change to test for between groups defaults to `0.6` (1.5 on a linear / non-log scale).
Set this to `0` (1 on a linear / non-log scale) to test for **any** difference between groups.

The database file defaults to `/NGS/known_sites/hg19/gencode.v19.chr_patch_hapl_scaff.annotation_UCSCcontigs.gtf`.


In addition, if a file named `candidates.txt` is present in the same folder as the input design table, containing gene IDs, one per line, present in the analyzed count files, all files generated for the differentially expressed genes will additionally be generated for the set of genes in this file.
This includes the PC and clustering plots as well as the count and test statistic tables.

If a file named `replace_sizeFactors.txt` is present, containing the size factors calculated by DESeq2 in a previous analysis to normalize for sequencing depth, these are used in the current analysis instead of whatever would be recalculated by DESeq2.



##### Usage
To compare treated versus untreated samples of the same patients using all of the provided defaults, run:
```
deseq2_wrapper.sh \
        your-design-table.tsv \
        "~ patient + treatment" \
        "untreated"
```
The provided TSV file must then contain a `patient` column, containing patient IDs of each sample, and a `treatment` column containing entries `treated` and `untreated`.

To change the provided defaults, the values with which to override the defaults must be added to command.
Because arguments are assigned based on their position, all arguments up to the last default that shall be changed must be provided. For example, to change the minimum log fold change to `1` (2 on the linear scale), run:

```
deseq2_wrapper.sh \
        your-design-table.tsv \
        "~ patient + treatment" \
        "untreated"  \
        0.1  \
        1
```

##### Outputs
Primary outputs of interest created by the script are:

* tables containing the normalized and log-transformed read counts of all samples
* a table containing the test statistics, fold changes, p- and FDR- estimates of each gene
* PDF plots of principal components 1-6 from PCA and heatmaps of different hierarchical clusterings, for differentially expressed genes and increasing subsets of genes with a high coefficient of variation (CV), as a sort of unsupervised exploratory comparison

### Gene Set Enrichment Analysis (GSEA)
This script takes DESeq2 results, ranks genes based on their log2 fold change estimates and calls on Broad's GSEA for gene set enrichment analysis.
It requires a database file for ensembl ID to HUGO gene symbol conversion, defined within the script, the GSEA JAR file and one or more gene set collections to test. The latter can be downloaded from MSigDB or custom-built.

To process the paired DESeq2 analysis from above, run:
```
gsea.sh \
        your-design-table_DESeq2~patient+treatment_all_woutNA.txt \
        your-favorite-gene-set-collection.gmt \
        your-second-favorite-gene-set-collection.gmt
```

A dated output folder is created with a subfolder for each individual analysis which in turn contains all of the generated output files.

---

## Circular RNAs

Several scripts are available to detect and quantify circular RNAs from our gene expression analysis pipelines:

* `get_circs.sh`
* `count_circs.sh`
* `count_linear-alternatives.sh`

### Identifying circular transcripts
The `get_circs.sh` script identifies circular RNAs supported by backsplice junctions, where downstream splice donors are joined to upstream splice acceptors to form circular RNAs.
It takes STAR output files containing the discovered chimeric junctions and chimeric reads, the BAM file containing all of the alignments, and the genome file of the reference that was used for this alignment:

```
get_circs.sh \
        your-input_Chimeric.out.junction \
        your-input_Chimeric.out.sam \
        your-input_alignSTAR.bam \
        stars_chrNameLength.txt
```

Two output files are of importance: 
`your-input_Chimeric.out.junction_circsAnnotatedFinal.txt` contains a table listing the coordinates, genes and supporting reads of each backsplice junction that was discovered. 
`your-input_Chimeric.out.junction_circsAnnotatedFinal_linearJunctionsCounts.tsv` contains in addition the number of reads supporting the backsplice junctions' linear alternatives, that is linear junctions involving the backsplice donor or acceptor. The supporting reads of these alternative linear junction can then be used in a downstream analysis to compare circular to linear transcript expression.


### Generating circular RNA / backsplice count files
The `count_circs.sh` script takes the `*linearJunctionCounts.tsv` files created by `get_circs.sh` as input and generates read count files, based on HTSeq's output format, which can then be further processed by DESeq2.
Two of these count files are created per sample: 
`*_circs.tsv` contains the supporting reads **of each circular RNA** / backsplice junction; `*_circs-per-gene.counts` contains the supporting reads of the circular RNAs / backsplice junctions **of each gene**.
Thus, for the latter, counts are summed when multiple circular RNAs / backsplice junctions are discovered within the same gene.

Specifically, the script takes as input a single output file prefix and all of the `*_linearJunctionCounts.tsv` files generated by `get_circs.sh`, which are to be processed together by DESeq2.
It is important that all of these samples are preprocessed together for count generation to ensure that the same features are counted for all samples, even when distinct sets of circRNAs were discovered.

```
count_circs.sh \
        your-output-prefix \
        your-input_*linearJunctionCounts.tsv \
        your-other-input_*linearJunctionCounts.tsv \
        your-yet-another-input_*linearJunctionCounts.tsv
```


The generated count files are all saved to the `intermediate_files` subfolder, which is also created by the script, and can then be input into DESeq2 to test for differential expression, as described [above](#differential-gene-expression-analysis).
A dummy DESeq2 design table listing the respective files is also created: `*_circs.tsv` lists per-junction read count files and `*_circs-per-gene.tsv` lists the per-gene files.
Results from DESeq2 can then be further studied using GSEA, just like conventional counts from linear transcripts.


### Comparing circular with linear transcript expression
A comparison of backsplice vs all non-backsplice reads is biased and does not properly reflect circular to linear transcript expression.
Instead, we compare backsplice-supporting read counts to those supporting alternative linear junctions involving either the backsplice donor or acceptor.
To quantify support (= expression) of these linear junctions, `count_linear-alternatives.sh` was implemented.

It works analogously to `count_circs.sh` and also takes as input an output file prefix and the `*_linearJunctionCounts.tsv` files which are to be processed further:

```
count_linear-alternatives.sh \
        your-output-prefix \
        your-input_*linearJunctionCounts.tsv \
        your-other-input_*linearJunctionCounts.tsv \
        your-yet-another-input_*linearJunctionCounts.tsv
```

Count files `*_linear-alternatives.counts` and `*_linear-alternatives-per-gene.counts` are again saved to the `intermediate_files` subfolder and dummy design tables for DESeq2 are prepared as `*_linear-alternatives.tsv` and `linear-alternatives-per-gene.tsv`.



---

## Additional helper scripts

### Sub-sampling FASTQ files
To extract a random set of `n` reads, use the `sample_fastqs.sh` script.
It takes as input the two paired FASTQ files to subsample and the number of reads to extract.
Thus, the command to extract 100 reads would be:

```
sample_fastqs.sh \
        your-input_R1.fastq \
        your-input_R2.fastq \
        100
```

Two output files will be written, one for each input file, containing the respectively extracted reads. Their filename prefix will be `your-input`, the suffix will contain the number of reads extracted together with a timestamp, dedicating the day, hour, minute and second that the script was run at. This serves to distinguish individual subsamples created from the same pair of original FASTQ files.

### Sorting & Padding BED files
All BED files must be sorted before they can be passed to our pipelines by chromosome contig and coordinate.
This can be done using the `sort_bed.sh` script, which takes as input the BED file to be sorted and the reference genome index file (.fai), which defines the proper order of reference contigs:

```
sort_bed.sh \
        your-input.bed \
        your-reference-genome-index.fai
```

In addition to a sorted BED file, `sort_bed.sh` also generates a padded BED file, which is also sorted, whose entries are extended by 5 bp each in both the 3' and 5' direction.
These padded BED files are required by targeted-enrichment pipeline to limit variant calling to the actual targets but include directly adjacent flanking regions as well, whose coverage is often also sufficent for analysis.



### Finding shared differentially expressed genes (DEGs)
The `get_shared_degs.sh` scripts takes as input an FDR cutoff and one or more DESeq2 result tables.
When one table is given, it returns all those genes which are differentially expressed with an FDR at or below the given cutoff.
When multiple files are given, the script returns only those DEGs shared by **all** of the provided files at or below the given FDR cutoff:

```
get_shared_degs.sh \
        your-fdr-cutoff \
        your-input_DESeq2*woutNA.txt \
        your-other-input_DESeq2*woutNA.txt \
        your-yet-another-input_DESeq2*woutNA.txt
```


### Finding shared GSEA gene sets
Two helper scripts were implemented:
* `get_shared_gsea_sets.sh`
* `get_shared_gsea_sets_for_all_combinations.sh`

The first of these, `get_shared_gsea_sets.sh` takes as input an FDR cutoff and one or more GSEA result tables.
It works analogously to `get_shared_degs.sh`:

```
get_shared_gsea_sets.sh \
        your-fdr-cutoff \
        your-gsea-result-tables-to-extract-shared-significant-gene-sets-from
```

The second script, `get_shared_gsea_sets_for_all_combinations.sh`, is essentially a wrapper for the first: Provided with an FDR cutoff, a ranking metric and one or more GSEA output file prefixes, shared significant gene sets per combination of ranking metric, gene set collection and direction of enrichment (enrichment vs depletion) are returned.
The current version of the `gsea.sh` script generates only analyses for gene lists ranked by log fold change, for which `lfc` has to be passed as the ranking metric.
The prefix has to uniquely identify the respective GSEA analyses:

```
get_shared_gsea_sets_for_all_combinations.sh \
        your-fdr-cutoff \
        lfc \
        your-design-table_DESeq2~patient+treatment_all_woutNA.txt \
        your-other-design-table_DESeq2~patient+treatment_all_woutNA.txt \
        your-yet-another-design-table_DESeq2~patient+treatment_all_woutNA.txt
```


### Testing for independence of a gene's expression and / or clinical variable(s)
The `test_association.R` script applies fisher's exact tests to test
for the independence of two variables.
These may be any attributes provided in the DESeq2 design table and / or
the expression level of a given gene. For genes, samples are split based on
the median expression and high vs low expression is the attribute that goes
into testing.

Required inputs are the DESeq2 design table and the name of the two attributes to test:
```
Rscript test_association.R \
        your-design-table \
        one-column-in-your-design-table \
        another-column-in-your-design-table
```

When gene symbols are given, the normalized counts table from DESeq2 is required as a fourth argument:
```
Rscript test_association.R \
        your-design-table \
        one-column-in-your-design-table \
        your-gene-of-interest \
        your-design-table_DESeq2~patient+treatment_all_countsNormalized.txt
```
(Except when the gene symbol is actually also a column in the design table. In this case, the design table will be used instead.)


### Creating bpipe stat summaries
Three scripts exist to merge the per-sample stat files created by bpipe into a single file describing all the samples of a given project:

* `summarize_coverage.sh`
* `summarize_duplicates.sh`
* `summarize_filtering.sh`

All of these scripts work analogously: They take as input a folder from which they collect the per-sample stat files to merge and then create a single output table with data of all of these.
When no input folder is given, it defaults to the current working directory from which the script is run.


##### Coverage stats
`summarize_coverage.sh` collects the `.summary` files from the given folder and merges them to a single file, `coverage_summary.tsv`.
The original `.summary` files contain the coverage statistics written by bpipe for a single given sample.

```
summarize_coverage.sh \
        folder-with-individual-coverage-stat-files
```

##### Read duplicate stats
`summarize_duplicates.sh` collects Picard's `.metrics.txt` files from the given folder and merges them to a single one, `bamStat_summary.tsv`.
These files contain the number and proportion of sequence and optical duplicates detected by Picard in the respective samples.

```
summarize_duplicates.sh  \
        folder-with-individual-duplicate-stat-files
```


##### Variant filtering stats
`summarize_filtering.sh` collects the `filter_statistic.txt` files from the given folder and merges them into one file, `filter_summary.tsv`.
The original `.summary` files contain the coverage statistics written by bpipe for a single given sample.

```
summarize_coverage.sh \
        folder-with-individual-coverage-stat-files
```

### Per-gene coverage
The `get_coverage_per_gene.sh` script calculates the average per-gene coverage for one or more samples processed by one of the mutation calling pipelines.
It takes as input the BED file which was also provided to the respective bpipe pipeline followed by all of the `*coverBED.txt` files generated by bpipe, for which the per-gene coverage is to be calculated. Column 4 of this BED file must contain the gene IDs / symbols of interest.
Output is a single file, `average_coverage_summary.tsv`, analogously to the summary scripts above.

```
get_coverage_per_gene.sh \
        your-targets.BED \
        your-input_*coverBED.txt \
        your-other-input_*coverBED.txt \
        your-yet-another-input_*coverBED.txt
```

### Filtering (for) specific mutations
The `filter_specific_mutations.sh` script filters all of the provided input files for a specified list of mutations. It was originally developed to filter recurrent mutations *not* of interest in a certain project. For each input file, mutations matching and not matching the provided list are written to distinct output files.

Note that the search is based on exact string matching.
Therefore, the list of mutations to filter for may contain any string or substring present (or not) in the CSV files.
For example, exact mutation coordinates, chromosomes or genes may be supplied, but all must be exactly formatted as expected in the CSV, with one search key per line.

```
filter_specific_mutations.sh \
        your-list-of-muts-to-filter-for.txt \
        your-inputs-mutations.csv \
        your-other-inputs-mutations.csv \
        your-yet-another-inputs-mutations.csv
```


### Annotating mutations' flanking sequences
The `get_flanking_sequence*.sh` scripts annotate mutations output by our pipelines with flanking reference sequences of a given length. Originally, `get_flanking_sequence.sh` was used, which is slow but works offline:

```
get_flanking_sequence.sh \
        your-inputs-mutations.csv \
        the-reference-genome.fasta \
        number-of-flanking-bases-to-annotate
```

Currently, our pipelines use `get_flanking_sequence_online.sh` for annotation of all of the detected mutations.
It requires the input CSV file and the number of flanking bases to annotate as well as a running internet connection:

```
get_flanking_sequence_online.sh \
        your-inputs-mutations.csv \
        number-of-flanking-bases-to-annotate
```


### Extracting molecular barcodes from within main-reads
For one of our protocols, molecular barcodes were sequenced not as separate index reads but as part of the main sequence reads.
To extract these into separate files, `get_barcodes.py` was implemented.

It takes as input the two paired FASTQ files from which the barcodes are to be extracted, an output file prefix for the extracted barcodes and the forward and reverse reads' primer sequences that signal the start / end of the respective barcode - these primers are then matched with the reads' sequences and preceding / succeeding sequence is extracted:

```
python3 get_barcodes.py \
        your-input_R1.fastq \
        your-input_R2.fastq \
        your-output-prefix \
        your-primer-sequence-in-your-input_R1 \
        your-primer-sequence-in-your-input_R2
```

Five output files are generated:
* `your-output-prefix_indexed_R1.fastq`: your-input_R1.fastq minus the extracted barcodes
* `your-output-prefix_indexed_R2.fastq`: your-input_R2.fastq minus the extracted barcodes
* `your-output-prefix_indexed_I1.fastq`: Barcodes extracted from your-input_R1.fastq
* `your-output-prefix_indexed_I2.fastq`: Barcodes extracted from your-input_R2.fastq
* `your-output-prefix_indexed_I.fastq`: Concatenation of the extracted barcodes into one per read pair

The file containing the concatenated indices can then, for example, be renamed / symlinked to `_R3.fastq` and passed as the index file of our HaloplexHS pipeline.


### Oxford Nanopore FAST5 to FASTQ conversion + stats
The `poretools_wrapper.sh` takes as input a folder containing Oxford Nanopore FAST5 files and converts these to a single FASTQ file.
Several statistics on run and data quality are also generated.
All output files are prefixed with the input folder's basename.

```
poretools_wrapper.sh \
        your-input-folder
```



### JAR file wrapper
The `run_jar.sh` script is a simple wrapper script that calls on a single JAR file, residing in the same folder as the wrapper script, via `java -jar ` and passes all arguments on to the JAR file. I implemented this to make JAR files callable from `$PATH`.




---

