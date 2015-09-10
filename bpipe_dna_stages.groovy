//file to contain pipeline stages of all DNA-seq pipeline: exome, amplicon, mouse, pindel

//// COMMON VARIABLE DEFINITIONS
// BASE DIRECTORY WITH PROGRAMS
NGS="/NGS" //symlink to respective folder

// PROGRAM DIRECTORIES
ANNOVAR="/NGS/links/annovar"
BEDTOOLS="/NGS/links/bedtools"
BWA="bwa"
CUTADAPT="cutadapt"
GATK="gatk"
PICARD="picard"
PINDEL="/NGS/links/pindel"
PLATYPUS="platypus"
SAMTOOLS="samtools"
VARSCAN="varscan"

getVersions = {
    exec "ls -l ${NGS}/links"
}

trim = { 
    //trim reads in fastq by triml bases on the left and trimr bases on the right (also adjust bqs!)
    var triml : 0
    var trimr : 0
    exec """sed '2~4s/^.\\{${triml}\\}\\(.*\\).\\{${trimr}\\}\$/\\1/g' $input | sed '4~4s/^.\\{${triml}\\}\\(.*\\).\\{${trimr}\\}\$/\\1/g' > $output.fastq"""
}


cutadapt = {
    exec "$CUTADAPT -m 50 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o $output1.fastq -p $output2.fastq   $input1.fastq   $input2.fastq"
}


alignSAMPE = { //it is not recommended to use this stage -> use alignMEM instead
    var nkern : 24 
    output.dir="intermediate_files"
    exec "$BWA aln -t $nkern $REF $input1.fastq > $output1.sai"
    exec "$BWA aln -t $nkern $REF $input2.fastq > $output2.sai"
    exec """$BWA sampe -P 
            -r "@RG\tID:Sample\tSM:Sample\tPL:illumina\tCN:exome" 
            $REF 
            $output1.sai 
            $output2.sai 
            $input1.fastq 
            $input2.fastq > 
            $output.sam"""
}


alignMEM = {
    var nkern : 24 
    output.dir="intermediate_files"
    exec """$BWA mem 
	-M
        -t $nkern 
	-R '@RG\tID:Sample\tSM:Sample\tPL:illumina\tCN:exome' 
        $REF
        $input1.fastq
        $input2.fastq > $output.sam"""
}


sed = { 
    //bpipe causes additional tabs in the samfile-header which is not good as it is supposed to be tab-deliminated for field-tags -> this sed reduces these additional tabs to spaces
    exec """sed -i -e '/^@PG/s/\t/ /5' 
		-e '/^@PG/s/\t/ /5' 
		-e '/^@PG/s/\t/ /5' 
		-e '/^@PG/s/\t/ /5' 
		$input.sam"""
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
        INPUT=$input.sam
        OUTPUT=$output.bam
	CREATE_INDEX=true
        SORT_ORDER=coordinate"""
}


indexPIC = {
    output.dir="intermediate_files"
    exec """$PICARD BuildBamIndex 
        INPUT=$input.bam"""
}


idxstatPIC = {
    output.dir="intermediate_files"
    exec """$PICARD BamIndexStats
        INPUT=$input.bam"""
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


coverBED = {
    output.dir="intermediate_files"
    var exon_cover : EXON_TARGET
    exec "${BEDTOOLS}/coverageBed -b $input.bam -a $exon_cover -d -sorted -g ${REF}.genomeFile > ${output.txt.prefix}_exon.txt"
    exec "${NGS}/calc_cvg/calc_cvg ${output.txt.prefix}_exon.txt $input1.fastq $input2.fastq > $output.txt"
    exec "rm ${output.txt.prefix}_exon.txt" //irgendwie loescht er das nicht?? 
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


baseRecalGATK = {
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
        -R $ref
        -BQSR $output.table
        -l INFO
        -I $input.bam
        -o $output.bam"""
}


mpileupSAM = {
    exec "$SAMTOOLS mpileup -f $REF -q 1 $input.bam > $output.pileup"
}


processSAM = segment {
	sortSAM +
        indexSAM + idxstatSAM +
        dedupSAM + indexSAM + idxstatSAM
}


processPICARD = segment {
	sortPIC + idxstatPIC +
        dedupPIC + idxstatPIC
}


PINDEL = {
    var nkern : 24 
    output.dir="intermediate_files"
    exec "echo $input.bam 250 ${input.bam.prefix}.PINDEL > $output.pindel_cfg.txt"
    exec "${PINDEL}/pindel -f $REF -i $output.pindel_cfg.txt -c ALL -T $nkern -x 5 -r -t -l -k -s -o $output.prefix"
    exec """${PINDEL}/pindel2vcf
       -r $REF
       -R ucsc.hg19
       -d 20130526
       -P $output.prefix
       -v $output.vcf
       -G"""
    exec "touch $output"
    exec """sed -e 's/chr//' $output.vcf | awk '{OFS="\t"; if (!/^#/){print \$1,\$2-1,\$2,\$4"/"\$5,"+"}}' > $output.bed"""
    forward output.vcf
}


processBED = {
    output.dir="intermediate_files"
    exec """${BEDTOOLS}/intersectBed -a $input2.vcf -b $input1.vcf > ${input2.prefix}.intersect.vcf"""
    exec "touch ${input2}_is_tumor.txt"
    exec """${BEDTOOLS}/subtractBed -a $input2.vcf -b $input1.vcf > ${input2.prefix}.subtract.vcf"""
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


haplocGATK = {
    var nkern : 24
    output.dir="intermediate_files"
    exec """$GATK
        -T HaplotypeCaller
        -nct $nkern
        -R $REF
        -I $input.bam
        -o $output.vcf"""
}


runEXOME_VARSCAN = segment {
	alignMEM +
        sed +
	processPICARD + 
	realignGATK +
        coverBED +
        mpileupSAM
}


runEXOME_VARSCAN_mouse = segment {
        alignMEM +
        sed +
        processPICARD +
        realignGATKwoutKnown +
        coverBED +
        mpileupSAM
}


runEXOME_PINDEL = segment {
	alignMEM +
	sed +
	processPICARD +
	realignGATK +
	indexSAM +
	PINDEL 
}


runEXOME_PLATYPUS = segment {
        alignMEM +
        sed +
        processPICARD +
	realignGATK +
        coverBED +
        unpairedPLATYPUS
}


somVARSC = {
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
        --somatic-p-value 0.05"""
    forward(glob("results_varscan/*.somVARSC.*"))
}


somVARSCunpaired = {
    output.dir="results_varscan"
    exec """$VARSCAN somatic
        $input1.pileup
        $input1.pileup
        results_varscan/${input1.prefix}_somVARSC
        --output-snp results_varscan/${input1.prefix}.somVARSCunpaired.snp
        --output-indel results_varscan/${input1.prefix}.somVARSCunpaired.indel
        --min-coverage 1
        --min-var-freq 0.01
        --min-freq-for-hom 0.75
        --normal-purity 1.0
        --tumor-purity 1.0
        --p-value 0.99 
        --somatic-p-value 0.05"""
    forward(glob("results_varscan/*.somVARSCunpaired.*"))
}  


amplicon = segment {
        alignMEM +
        sed +
        sortPIC +
        indexPIC + idxstatPIC +
        realignGATK + 
        coverBED +
        mpileupSAM +
        somVARSCunpaired
}


sed_rfmt = {
    exec "${NGS}/snp_filter/snp_filter $input intermediate_files/\$(basename ${input}.pass) intermediate_files/\$(basename ${input}.dumped) 0 0 0 0 0 0"
    forward(glob("intermediate_files/\$(basename ${input}.pass)"))
}


processVARSC = {
    exec "$VARSCAN processSomatic intermediate_files/\$(basename ${input}.pass)"
    forward(glob("intermediate_files/*.pass.*"))
}


doublePos = {
    exec "sed -e '0,/chr.*\$/{s/chr.*\$/&\\n&/}' $input > intermediate_files/\$(basename ${input}.pre_rwrt)"
    exec "${NGS}/snp_rwrt/snp_rwrt intermediate_files/\$(basename ${input}.pre_rwrt) intermediate_files/\$(basename ${input}.rwrt)"
    exec "sed -i -e '1d' intermediate_files/\$(basename ${input}.rwrt)"
    forward """${input}.rwrt"""
}


tableANNOVAR = {
    exec """${ANNOVAR}/table_annovar.pl $input ${ANNOVAR}/humandb/ -buildver hg19 -out $input -remove -protocol refGene,genomicSuperDups,esp6500_all,1000g2014sep_all,snp138,cosmic70,ljb23_pp2hdiv,ljb23_sift -operation g,r,f,f,f,f,f,f -nastring '"."' -csvout -otherinfo"""
    forward(glob("intermediate_files/*.csv"))
}


tableANNOVARmm10 = {
    exec """${ANNOVAR}/table_annovar.pl $input ${ANNOVAR}/mm10db/ -buildver mm10 -out $input -remove -protocol refGene,genomicSuperDups,snp138 -operation g,r,f -nastring '"."' -csvout -otherinfo"""
    forward(glob("intermediate_files/*.csv"))
}


merged = {
    output.dir="results_csv"
    exec "head -n 1 intermediate_files/*snp*Somatic.hc*.csv > results_csv/${input3.fastq.prefix}_merged.csv"
    exec "tail -n +2 intermediate_files/*Germline.rwrt*.csv | cat >> results_csv/${input3.fastq.prefix}_merged.csv"
    exec "tail -n +2 intermediate_files/*LOH.rwrt*.csv | cat >> results_csv/${input3.fastq.prefix}_merged.csv"
    exec "tail -n +2 intermediate_files/*Somatic.rwrt*.csv | cat >> results_csv/${input3.fastq.prefix}_merged.csv" 
    exec """sed -i -e '/^\$/d' results_csv/${input3.fastq.prefix}_merged.csv"""
    exec """sed -i -e '/^==>/d' results_csv/${input3.fastq.prefix}_merged.csv"""
    forward(glob("results_csv/*merged.csv"))
}


mergedAmplicon = {
    output.dir="results_csv"
    exec "head -n 1 intermediate_files/${input1.fastq.prefix}*snp*Germline.rwrt*.csv > results_csv/${input1.fastq.prefix}_merged.csv"
    exec "tail -n +2 intermediate_files/${input1.fastq.prefix}*Germline.rwrt*.csv | cat >> results_csv/${input1.fastq.prefix}_merged.csv"
    exec """sed -i -e '/^\$/d' results_csv/${input1.fastq.prefix}_merged.csv"""
    exec """sed -i -e '/^==>/d' results_csv/${input1.fastq.prefix}_merged.csv"""
    forward(glob("results_csv/${input1.fastq.prefix}*merged.csv"))
}


final_sed = {
    exec "cut -f1-5 -d',' $input > intermediate_files/\$(basename ${input}_not_otherinfo1.txt)"
    exec """cut -f6- -d',' $input | sed -e '2,\$ s/","/"___re___"/g' -e 's/""/"/g' -e '2,\$ s/,/;/g' -e 's/___re___/,/g'  > intermediate_files/\$(basename ${input}_not_otherinfo2.txt)"""
    exec "cut -f1-12 -d',' intermediate_files/\$(basename ${input}_not_otherinfo2.txt) > intermediate_files/\$(basename ${input}_not_otherinfo3.txt)"
    exec """cut -f13- -d',' intermediate_files/\$(basename ${input}_not_otherinfo2.txt) | sed -e 's/^"//g' > intermediate_files/\$(basename ${input}_otherinfo.txt)"""
    exec "paste intermediate_files/\$(basename ${input}_not_otherinfo1.txt) intermediate_files/\$(basename ${input}_not_otherinfo3.txt) intermediate_files/\$(basename ${input}_otherinfo.txt) > $input"
    exec """sed  -i -e 's/	/,/g' -e 's/,"\$//g' -e 's/Otherinfo/Otherinfo=(normal_reads1,normal_reads2,normal_var_freq,normal_gt,tumor_reads1,tumor_reads2,tumor_var_freq,tumor_gt,somatic_status,variant_p_value,somatic_p_value,tumor_reads1_plus,tumor_reads1_minus,tumor_reads2_plus,tumor_reads2_minus,normal_reads1_plus,normal_reads1_minus,normal_reads2_plus,normal_reads2_minus)?/' $input"""
}

final_sed_mm10 = {
    exec "cut -f1-5 -d',' $input > intermediate_files/\$(basename ${input}_not_otherinfo1.txt)"
    exec """cut -f6- -d',' $input | sed -e '2,\$ s/","/"___re___"/g' -e 's/""/"/g' -e '2,\$ s/,/;/g' -e 's/___re___/,/g'  > intermediate_files/\$(basename ${input}_not_otherinfo2.txt)"""
    exec "cut -f1-7 -d',' intermediate_files/\$(basename ${input}_not_otherinfo2.txt) > intermediate_files/\$(basename ${input}_not_otherinfo3.txt)"
    exec """cut -f8- -d',' intermediate_files/\$(basename ${input}_not_otherinfo2.txt) | sed -e 's/^"//g' > intermediate_files/\$(basename ${input}_otherinfo.txt)"""
    exec "paste intermediate_files/\$(basename ${input}_not_otherinfo1.txt) intermediate_files/\$(basename ${input}_not_otherinfo3.txt) intermediate_files/\$(basename ${input}_otherinfo.txt) > $input"
    exec """sed  -i -e 's/      /,/g' -e 's/,"\$//g' -e 's/Otherinfo/Otherinfo=(normal_reads1,normal_reads2,normal_var_freq,normal_gt,tumor_reads1,tumor_reads2,tumor_var_freq,tumor_gt,somatic_status,va
riant_p_value,somatic_p_value,tumor_reads1_plus,tumor_reads1_minus,tumor_reads2_plus,tumor_reads2_minus,normal_reads1_plus,normal_reads1_minus,normal_reads2_plus,normal_reads2_minus)?/' $input"""
}


finalSedPINDEL = {
    exec "cut -f1-5 -d',' $input > intermediate_files/\$(basename ${input}_not_otherinfo1.txt)"
    exec """cut -f6- -d',' $input | sed -e '2,\$ s/","/"___re___"/g' -e 's/""/"/g' -e '2,\$ s/,/;/g' -e 's/___re___/,/g'  > intermediate_files/\$(basename ${input}_not_otherinfo2.txt)"""
    exec "cut -f1-12 -d',' intermediate_files/\$(basename ${input}_not_otherinfo2.txt) > intermediate_files/\$(basename ${input}_not_otherinfo3.txt)"
    exec """cut -f13- -d',' intermediate_files/\$(basename ${input}_not_otherinfo2.txt) | sed -e 's/    /","/g' > intermediate_files/\$(basename ${input}_otherinfo.txt)"""
    exec "paste -d',' intermediate_files/\$(basename ${input}_not_otherinfo1.txt) intermediate_files/\$(basename ${input}_not_otherinfo3.txt) intermediate_files/\$(basename ${input}_otherinfo.txt) > $input"
    exec """sed  -i -e 's/      /,/g' -e 's/,"\$//g' $input"""
}


filterOutput = {
      var candidates : "${NGS}/candidate_genes_aml.txt"
      exec "echo \$(date) running $FILTER"
      exec "$FILTER $input.csv $input.csv.prefix $candidates $REF 15"
}



