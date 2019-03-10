//file to contain pipeline stages of all DNA-seq pipeline: exome, amplicon, mouse, pindel

//// COMMON VARIABLE DEFINITIONS
// BASE DIRECTORY WITH PROGRAMS
NGS="/NGS" //symlink to respective folder

// PROGRAM DIRECTORIES
ANNOVAR="/NGS/links/annovar"
BEDTOOLS="/NGS/links/bedtools"
BWA="bwa"
BWA_LONG="${NGS}/bwa-long/bwa/bwa"
CUTADAPT="cutadapt"
GATK="gatk"
QC="ngsqc"
PICARD="picard"
PINDEL="/NGS/links/pindel"
PLATYPUS="platypus"
SAMTOOLS="samtools"
VARSCAN="varscan"
NGSBITS="/NGS/links/ngs-bits"

getVersions = {
    exec "ls -l ${NGS}/links"
}

qc = {
    var nkern : 24
    // subsample FASTQ files for QC -> take every 25th read only -> (NR = read fraction ^-1 * 4 = 25 * 4 = 100)
    exec "awk 'NR % 100 > 0 && NR % 100 < 5' $input1.fastq > $output1.subsample"
    exec "awk 'NR % 100 > 0 && NR % 100 < 5' $input2.fastq > $output2.subsample"
    exec "$QC -pe $output1.subsample $output2.subsample 2 A -c $nkern -onlyStat -o FASTQ_QC && touch $output.success"
}

trim = { 
    //trim reads in fastq by triml bases on the left and trimr bases on the right (also adjust bqs!)
    var triml : 0
    var trimr : 0
    exec """sed '2~4s/^.\\{${triml}\\}\\(.*\\).\\{${trimr}\\}\$/\\1/g' $input | sed '4~4s/^.\\{${triml}\\}\\(.*\\).\\{${trimr}\\}\$/\\1/g' > $output.fastq"""
}

trim_haloplexC = {
    exec "${NGSBITS}/FastqTrim -in $input1.fastq -out $output1.fastq -end 1"
    exec "${NGSBITS}/FastqTrim -in $input2.fastq -out $output2.fastq -start 1"
    exec "cut -f1 -d ' ' $input3.fastq > $output3.fastq"
}

cutadapt = {
    var nkern : 24
    exec "$CUTADAPT --cores=$nkern -m 50 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o $output1.fastq -p $output2.fastq $input1.fastq $input2.fastq"
}


alignMEM = {
    //bpipe causes additional tabs in the samfile-header which is not good as it is supposed to be tab-deliminated for field-tags -> this sed reduces these additional tabs to spaces
    var nkern : 24 
    output.dir="intermediate_files"
    exec """$BWA mem 
	-M
        -t $nkern 
	-R "@RG\tID:$input1.prefix\tSM:$input1.prefix\tPL:illumina\tCN:exome" 
        $REF
        $input1.fastq
        $input2.fastq | sed -e '/^@PG/s/\t/ /5' 
                		-e '/^@PG/s/\t/ /5' 
               			-e '/^@PG/s/\t/ /5' 
                		-e '/^@PG/s/\t/ /5' | $PICARD SortSam
                                                        INPUT=/dev/stdin
                                                        OUTPUT=$output.bam
                                                        CREATE_INDEX=true
                                                        SORT_ORDER=coordinate"""
}

// forward 3rd fastq = index file -> that is the only diff to alignMEM above!
alignMEMhaloplex = {
    //bpipe causes additional tabs in the samfile-header which is not good as it is supposed to be tab-deliminated for field-tags -> this sed reduces these additional tabs to spaces
    var nkern : 24
    output.dir="intermediate_files"
    exec """$BWA mem 
        -M
        -t $nkern 
        -R "@RG\tID:$input1.prefix\tSM:$input1.prefix\tPL:illumina\tCN:exome" 
        $REF
        $input1.fastq
        $input2.fastq | sed -e '/^@PG/s/\t/ /5' 
                                -e '/^@PG/s/\t/ /5' 
                                -e '/^@PG/s/\t/ /5' 
                                -e '/^@PG/s/\t/ /5' | $PICARD SortSam
                                                        INPUT=/dev/stdin
                                                        OUTPUT=$output.bam
                                                        CREATE_INDEX=true
                                                        SORT_ORDER=coordinate"""
      forward(output.bam, input3.fastq)
}


alignMEMlong = {
    var nkern : 48 
    output.dir="intermediate_files"
    exec """$BWA_LONG mem 
	-M
        -t $nkern 
	-R "@RG\tID:$input1.prefix\tSM:$input1.prefix\tPL:illumina\tCN:exome" 
	-x ont2d
        $REF
        $input.fastq | sed -e '/^@PG/s/\t/ /5' 
                		-e '/^@PG/s/\t/ /5' 
               			-e '/^@PG/s/\t/ /5' 
                		-e '/^@PG/s/\t/ /5' | $PICARD SortSam
                                                        INPUT=/dev/stdin
                                                        OUTPUT=$output.bam
                                                        CREATE_INDEX=true
                                                        SORT_ORDER=coordinate"""
}

sortSAM = {
    output.dir="intermediate_files"
    exec "$SAMTOOLS view -bSh $input.sam > ${input.prefix}.bam"
    exec "$SAMTOOLS sort -o ${input.prefix}.bam $output.prefix > $output.bam"
}


indexSAM = {
    output.dir="intermediate_files"
    exec  "$SAMTOOLS index $input.bam"
}


idxstatSAM = {
    output.dir="intermediate_files"
    exec  "$SAMTOOLS idxstats $input.bam > $output.txt"
}


dedupSAM = {
    output.dir="intermediate_files"
    exec "$SAMTOOLS rmdup $input.bam $output.bam"
}


sortPIC = {
    output.dir="intermediate_files"
    exec """  
	$PICARD SortSam
        INPUT=$input.bam
        OUTPUT=$output.bam
	CREATE_INDEX=true
        SORT_ORDER=coordinate"""
}


indexPIC = {
    output.dir="intermediate_files"
    transform("bam") to ("bai") {
        exec """$PICARD BuildBamIndex
            INPUT=$input.bam"""
    }
}


idxstatPIC = {
    output.dir="intermediate_files"
    exec """$PICARD BamIndexStats
        INPUT=$input.bam > $output.tsv"""
    forward input.bam
}


dedupPIC = {
    output.dir="intermediate_files"
    exec """$PICARD MarkDuplicates
	INPUT=$input.bam
	OUTPUT=$output.bam
	REMOVE_DUPLICATES=true
	CREATE_INDEX=true
	METRICS_FILE=${output.bam.prefix}_metrics.txt"""
}

dedupOptPIC = {
    output.dir="intermediate_files"
    exec """$PICARD MarkDuplicates
	INPUT=$input.bam
	OUTPUT=$output.bam
	REMOVE_SEQUENCING_DUPLICATES=true
	CREATE_INDEX=true
	METRICS_FILE=${output.bam.prefix}_metrics.txt"""
}

unmarkDupsPIC = {
    output.dir="intermediate_files"
    exec """$PICARD RevertSam
        INPUT=$input.bam
        OUTPUT=$output.bam
        CREATE_INDEX=true
        SORT_ORDER=coordinate
        REMOVE_DUPLICATE_INFORMATION=true
        REMOVE_ALIGNMENT_INFORMATION=false"""
}

dedupBarcode = {
    output.dir="intermediate_files"
    var exon_cover : EXON_TARGET
    exec "${NGSBITS}/BamDeduplicateByBarcode -bam $input.bam -index $input.fastq -out $output.bam -dist 1 -min_group 1 -hs_file $exon_cover"
}

coverBED = {
    output.dir="intermediate_files"
    var exon_cover : EXON_TARGET
    exec "${BEDTOOLS}/coverageBed -b $input.bam -a $exon_cover -hist -sorted -g ${REF}.genomeFile > $output.txt"
    exec "grep '^all' $output.txt | sort -k2,2nr | awk -v OFS='\t' -v CUMSUM=0 -v CUMFREQ=0.0 -v TOTALONTARGET=0 '{CUMSUM=CUMSUM+\$3; CUMFREQ=(100 * CUMSUM/\$4); TOTALONTARGET=TOTALONTARGET + (\$2*\$3); print \$0,CUMSUM,CUMFREQ,TOTALONTARGET;}' | sort -k2,2n > $output.nice"
    exec "formatCoverage.sh $output.nice $input1.fastq $input2.fastq 0 1 10 15 50 100 120 200 500 1000 1500 2000 2500 > $output.summary"
    forward input.bam
}

realignGATK = {
    output.dir="intermediate_files"
    var nkern : 24
    exec """$GATK
        -T RealignerTargetCreator
        -nt $nkern
        -R $REF
        -known $GOLD_STANDARD_1000G_INDELS
        -known $PHASE1_1000G_INDELS
        -o $output.list"""
    exec """$GATK
        -T IndelRealigner
        -R $REF
        -targetIntervals $output.list
        -I $input.bam
        -known $GOLD_STANDARD_1000G_INDELS
        -known $PHASE1_1000G_INDELS
        -o $output.bam"""
}


realignGATKwoutKnown = {
    output.dir="intermediate_files"
    var nkern : 24
    exec """$GATK
        -T RealignerTargetCreator
        -nt $nkern
        -R $REF
	-I $input.bam
        -o $output.list"""
    exec """$GATK
        -T IndelRealigner
        -R $REF
        -targetIntervals $output.list
        -I $input.bam
        -o $output.bam"""
}


bqsrGATK = {
    var nkern : 24
    output.dir="intermediate_files"
    exec """$GATK
        -T BaseRecalibrator
        -nct $nkern
        -R $REF
        -l INFO
        -knownSites $DBSNP
        -knownSites $GOLD_STANDARD_1000G_INDELS
        -knownSites $PHASE1_1000G_INDELS
        -cov ReadGroupCovariate
        -cov QualityScoreCovariate
        -cov CycleCovariate
        -cov ContextCovariate
        -I $input.bam
        --out $output.table"""
    exec """$GATK
        -T PrintReads
        -R $REF
        -BQSR $output.table
        -l INFO
        -I $input.bam
        -o $output.bam"""
}


mpileupSAMpad = {
    var exon_cover : EXON_TARGET
    var bqs : 25
    exec """$SAMTOOLS mpileup -f $REF -q 1 -Q $bqs -B -l ${exon_cover}_padded -d 1000000 $input.bam > $output.pileup"""
}

//additional stage to use non-padded BED for Amplicons
mpileupSAMexact = {
    var exon_cover : EXON_TARGET
    var bqs : 25
    exec "$SAMTOOLS mpileup -f $REF -q 1 -Q $bqs -B -l $exon_cover -d 1000000 $input.bam > $output.pileup"
}

//additional stage without BED for variant calling from RNA-seq
mpileupSAM_noBed = {
    var bqs : 25
    exec "$SAMTOOLS mpileup -f $REF -q 1 -Q $bqs -B -d 1000000 $input.bam > $output.pileup"
}

//additional stage without BED and quality filters for VarSim
mpileupSAM_varsim = {
    exec "$SAMTOOLS mpileup -f $REF -d 1000000 $input.bam > $output.pileup"
}

processSAM = segment {
	idxstatSAM +
        dedupSAM + indexSAM + idxstatSAM
}


processPICARD = segment {
	idxstatPIC +
        dedupPIC + idxstatPIC
}


runPINDEL = {
    var nkern : 24 
    var region : 'chr13:28,577,689-28,675,147'
    output.dir="intermediate_files"
    exec "$SAMTOOLS view -hb $input.bam $region > $output.bam"
    exec "$SAMTOOLS index $output.bam && touch $output.index_dummy"
    exec "echo $output.bam 250 ${output.bam.prefix} > $output.config"
    exec "${PINDEL}/pindel -f $REF -i $output.config -c ALL -T $nkern -x 5 -r -t -l -k -s -c $region -o $output.prefix"
    exec """${PINDEL}/pindel2vcf
       -r $REF
       -R ucsc.hg19
       -d 20130526
       -P $output.prefix
       -v $output.vcf
       -G"""
    exec "touch $output"
    forward output.vcf
}


processBED = {
    output.dir="intermediate_files"
    exec """${BEDTOOLS}/intersectBed -a $input2.vcf -b $input1.vcf > ${input2.prefix}.intersect.vcf"""
    exec """${BEDTOOLS}/subtractBed -a $input2.vcf -b $input1.vcf > ${input2.prefix}.subtract.vcf"""
    exec """${BEDTOOLS}/subtractBed -b $input2.vcf -a $input1.vcf > ${input1.prefix}.subtract.vcf"""
    forward(glob("intermediate_files/*ct.vcf"))
}


conv2ANNO = {
    output.dir="intermediate_files"
    exec "${ANNOVAR}/convert2annovar.pl -format vcf4 $input.vcf -includeinfo --outfile ${input.prefix}.avinput"
    forward(glob("intermediate_files/*.avinput"))
}


unpairedPLATYPUS = {
    var nkern : 24
    output.dir="intermediate_files"
    exec "$PLATYPUS callVariants --refFile=$REF --bamFiles=$input.bam --assemble=1 --assembleBrokenPairs=1 --nCPU=$nkern --output=$output.vcf"
}

pairedPLATYPUS = {
    var nkern : 48
    output.dir="intermediate_files"
    exec "$PLATYPUS callVariants --refFile=$REF --bamFiles=$input1.bam,$input2.bam --assemble=1 --assembleBrokenPairs=1 --nCPU=$nkern --output=$output.vcf"
}

haplocGATK = {
    var nkern : 48
    output.dir="intermediate_files"
    exec """$GATK
        -T HaplotypeCaller
        -nct $nkern
        -R $REF
        -I $input1.bam
        -I $input2.bam
        -o $output.vcf"""
}


runEXOME_VARSCAN = segment {
	alignMEM +
	processPICARD + 
	realignGATK +
        [ coverBED,  mpileupSAMpad ]
}


runEXOME_VARSCAN_mouse = segment {
        alignMEM +
        processPICARD +
        realignGATKwoutKnown +
        [ coverBED, mpileupSAMpad ]
}


runEXOME_PINDEL = segment {
	alignMEM +
	processPICARD +
	realignGATK +
	indexSAM +
	runPINDEL 
}


runEXOME_PLATYPUS = segment {
        alignMEM +
        processPICARD +
	realignGATK +
        coverBED +
        unpairedPLATYPUS
}


somVARSC = {
    output.dir="results_varscan"
    exec """$VARSCAN somatic
        $input2.pileup
        $input1.pileup
        --output-snp $output.snp
        --output-indel $output.indel
        --min-coverage 1
        --min-var-freq 0.01
        --min-freq-for-hom 0.75
        --normal-purity 1.0
        --tumor-purity 1.0
        --p-value 0.99
        --somatic-p-value 0.05"""
}


somVARSCvcf = {
    output.dir="results_varscan"
    exec """$VARSCAN somatic
        $input1.pileup
        $input2.pileup
        results_varscan/${input2.prefix}_somVARSC
        --output-snp results_varscan/${input2.prefix}.somVARSC.snp
        --output-indel results_varscan/${input2.prefix}.somVARSC.indel
        --min-coverage 1
        --min-var-freq 0.01
        --min-freq-for-hom 0.75
        --normal-purity 1.0
        --tumor-purity 1.0
        --p-value 0.99
        --somatic-p-value 0.05
        --output-vcf 1"""
    forward(glob("results_varscan/*.somVARSC.*"))
}

somVARSCunpaired = {
    output.dir="results_varscan"
    exec """$VARSCAN somatic
        $input1.pileup
        $input1.pileup
        --output-snp $output.snp
        --output-indel $output.indel
        --min-coverage 1
        --min-var-freq 0.01
        --min-freq-for-hom 0.75
        --normal-purity 1.0
        --tumor-purity 1.0
        --p-value 0.99 
        --somatic-p-value 0.05"""
}  

amplicon = segment {
        alignMEM +
	dedupOptPIC + unmarkDupsPIC +
        realignGATK + 
        [ coverBED, mpileupSAMexact ] +
        somVARSCunpaired
}

haloplex = segment {
	trim_haloplexC +
        alignMEMhaloplex +
	dedupBarcode +
	sortPIC +
	idxstatPIC +
        realignGATK + 
        [ coverBED, mpileupSAMexact ] +
        somVARSCunpaired
}

processVARSC = {
    produce(input + ".Germline", input + ".LOH", input + ".Somatic") {
        exec """$VARSCAN processSomatic $input
        &&
        rm -f results_varscan/*.hc"""
    }
}


fixFormat = {
    output.dir="intermediate_files"
    produce(input + ".fixFormat.pre_rwrt", input + ".fixFormat.rwrt", input + ".fixFormat.fixed") {
        // delete header and fix indel ref and alt alleles to match annovar input specifications
        exec """sed
            -e '1d'
            -e 's/^\\([^\t]*\t[^\t]*\t\\).[^\t]*\t+\\([^\t]*\\)/\\1-\t\\2/'
            -e '0,/chr.*\$/{s/chr.*\$/&\\n&/}'
                $input  > $output.pre_rwrt"""
        exec "${NGS}/snp_rwrt/snp_rwrt $output.pre_rwrt $output.rwrt"

        // fix del coords:
        exec """
            tail -n +2 $output.rwrt |
            awk -v OFS="\t" 'substr(\$5,1,1) == "-" {\$2=\$2+1; \$3=\$3+length(\$5)-1; \$4=substr(\$5,2,length(\$5)); \$5="-"; print \$0} substr(\$5,1,1) != "-" {print \$0}' > $output.fixed"""
    }
    forward """$output.fixed"""
}


tableANNOVAR = {
    produce(input + ".hg19_multianno.csv"){
        exec """${ANNOVAR}/table_annovar.pl $input
            ${ANNOVAR}/humandb/
            -buildver hg19
            -out $input
            -remove
            -protocol refGene,genomicSuperDups,esp6500_all,1000g2014sep_all,snp138,cosmic70,ljb23_pp2hdiv,ljb23_sift
            -operation g,r,f,f,f,f,f,f
            -nastring '"."'
            -csvout
            -otherinfo"""
    }
}


tableANNOVARmm10 = {
    produce(input + ".mm10_multianno.csv") {
        exec """${ANNOVAR}/table_annovar.pl $input
            ${ANNOVAR}/mm10db/
            -buildver mm10
            -out $input
            -remove
            -protocol refGene,genomicSuperDups,snp138
            -operation g,r,f
            -nastring '"."'
            -csvout
            -otherinfo"""
    }
}


merged = {
    output.dir="results_csv"
    from("*multianno.csv") produce(input1.fastq.prefix + "_merged.csv") {
        exec """
            head -n 1 --quiet $inputs |
            uniq |
            sed -e 's/\t"\$//'
                -e 's/,"/\t/g'
                -e 's/"//g'
                -e 's/,/\t/4'
                -e 's/,/\t/3'
                -e 's/,/\t/2'
                -e 's/,/\t/1'
                -e 's/,/;/g'
                -e 's/\t/","/g'
                -e 's/^/"/'
                -e 's/\$/"/'
                -e 's/Otherinfo/normal_reads1","normal_reads2","normal_var_freq","normal_gt","tumor_reads1","tumor_reads2","tumor_var_freq","tumor_gt","somatic_status","variant_p_value","somatic_p_value","tumor_reads1_plus","tumor_reads1_minus","tumor_reads2_plus","tumor_reads2_minus","normal_reads1_plus","normal_reads1_minus","normal_reads2_plus","normal_reads2_minus/'
                -e '1s/;/","/g' > $output
            &&
            tail -n +2 --quiet $inputs |
            sed -e 's/\t"\$//'
                -e 's/,"/\t/g'
                -e 's/"//g'
                -e 's/,/\t/4'
                -e 's/,/\t/3'
                -e 's/,/\t/2'
                -e 's/,/\t/1'
                -e 's/,/;/g'
                -e 's/\t/","/g'
                -e 's/^/"/'
                -e 's/\$/"/' >> $output"""
    }
}

mergedAmplicon = {
    output.dir="results_csv"
    produce(input.split("/")[-1].split("R")[0] + "merged.csv") {
        exec """
            head -n 1 --quiet $inputs |
            uniq |
            sed -e 's/\t"\$//'
                -e 's/,"/\t/g'
                -e 's/"//g'
                -e 's/,/\t/4'
                -e 's/,/\t/3'
                -e 's/,/\t/2'
                -e 's/,/\t/1'
                -e 's/,/;/g'
                -e 's/\t/","/g'
                -e 's/^/"/'
                -e 's/\$/"/'
                -e 's/Otherinfo/normal_reads1","normal_reads2","normal_var_freq","normal_gt","tumor_reads1","tumor_reads2","tumor_var_freq","tumor_gt","somatic_status","variant_p_value","somatic_p_value","tumor_reads1_plus","tumor_reads1_minus","tumor_reads2_plus","tumor_reads2_minus","normal_reads1_plus","normal_reads1_minus","normal_reads2_plus","normal_reads2_minus/'
                -e '1s/;/","/g' > $output
            &&
            tail -n +2 --quiet $inputs |
            sed -e 's/\t"\$//'
                -e 's/,"/\t/g'
                -e 's/"//g'
                -e 's/,/\t/4'
                -e 's/,/\t/3'
                -e 's/,/\t/2'
                -e 's/,/\t/1'
                -e 's/,/;/g'
                -e 's/\t/","/g'
                -e 's/^/"/'
                -e 's/\$/"/' >> $output"""
    }
}


filterOutput = {
    var candidates : CANDIDATES
    produce(input.prefix + "_filter_statistic.txt", input.prefix + "*.csv") {
        exec "echo \$(date) running $FILTER"
        exec "$FILTER $input $input.prefix $candidates $REF 15"
    }
}


cleanUp = {
    exec "bpipe_finish.sh \$(basename $input)"
}



finalSedPINDEL = {
    exec "cut -f1-5 -d',' $input > intermediate_files/\$(basename ${input}_not_otherinfo1.txt)"
    exec """cut -f6- -d',' $input | sed -e '2,\$ s/","/"___re___"/g' -e 's/""/"/g' -e '2,\$ s/,/;/g' -e 's/___re___/,/g'  > intermediate_files/\$(basename ${input}_not_otherinfo2.txt)"""
    exec "cut -f1-12 -d',' intermediate_files/\$(basename ${input}_not_otherinfo2.txt) > intermediate_files/\$(basename ${input}_not_otherinfo3.txt)"
    exec """cut -f13- -d',' intermediate_files/\$(basename ${input}_not_otherinfo2.txt) | sed -e 's/    /","/g' > intermediate_files/\$(basename ${input}_otherinfo.txt)"""
    exec "paste -d',' intermediate_files/\$(basename ${input}_not_otherinfo1.txt) intermediate_files/\$(basename ${input}_not_otherinfo3.txt) intermediate_files/\$(basename ${input}_otherinfo.txt) > $input"
    exec """sed  -i -e 's/      /,/g' -e 's/,"\$//g' $input"""
}

finalSedPLATYPUS = {
    exec """paste  -d ',' <(cut -d',' -f1-5 $input.csv)  <(cut -d',' -f6- $input.csv | sed -e 's/,"/\t/g' -e 's/"//g' -e 's/,/;/g' -e 's/\t/,/g' | cut -d',' -f1-12) <(cut -f1 $input.csv | sed 's/^.*,\\([^,]*\\)\$/\\1/') <(cut -f2-7 $input.csv | sed -e 's/\t/,/g')  <(cut -f8 $input.csv | sed -e 's/;/\t/g' -e 's/,/;/g' -e 's/\t/,/g') <(cut -f9 $input.csv) <(cut -f10 $input.csv | sed -e 's/,/;/g' -e 's/:/,/g') | sed -e 's/^/"/' -e 's/,/","/g' -e 's/""/"/g' | sed '1 s/^.*\$/"Chr","Start","End","Ref","Alt","Func.refGene","Gene.refGene","GeneDetail.refGene","ExonicFunc.refGene","AAChange.refGene","genomicSuperDups","esp6500_all","1000g2014sep_all","snp138","cosmic70","ljb23_pp2hdiv","ljb23_sift","chr1","pos","NA","ref","alt",,"PASS","Fraction of reads around this variant that failed filters","Estimated population frequency of variant","Homopolymer run length around variant locus","Haplotype score measuring the number of haplotypes the variant is segregating into in a window","Worst goodness-of-fit value reported across all samples","Median minimum base quality for bases around variant","Root mean square of mapping qualities of reads at the variant position","tumor_reads_2plus|Total number of forward reads containing this variant","Total number of reverse reads tumor_reads_2_minuscontaining this variant","Posterior probability (phred scaled) that this variant segregates","Variant-quality|read-depth for this variant","Variants fail sequence-context filter. Surrounding sequence is low-complexity","Binomial P-value for strand bias test","Sourse_of_variants","Total coverage at this locus","Total forward strand coverage at this locus","Total reverse strand coverage at this locus","Total number of reads containing this variant","End position of calling window","Starting position of calling window","description","Genotype","Genotype log10-likelihoods for AA;AB and BB genotypes; where A = ref and B = variant","Variant fails goodness-of-fit test","Genotype quality as phred score","Number of reads covering variant location in this sample","Number of reads containing variant in this sample"/' > $output.csv"""
}
