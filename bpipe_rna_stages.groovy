//file to contain pipeline stages of all RNA-seq pipelines

//// COMMON VARIABLE DEFINITIONS
// PROGRAM DIRECTORIES
STAR="star"

alignSTAR = {
        var nkern : 48 
        output.dir="intermediate_files"
        exec """$STAR
                --runThreadN $nkern
                --genomeDir $REFDIR 
                --twopassMode Basic
                --outStd BAM_SortedByCoordinate
                --outSAMtype BAM SortedByCoordinate
                --outFileNamePrefix ./.$input1.fastq.prefix
                --outSAMunmapped Within
                --outSAMmapqUnique 60 
                --outSAMattrRGline ID:0 PL:ILLUMINA LB:rna SM:$input1.fastq.prefix
                --outSAMstrandField intronMotif 
                --outSAMattributes NH HI AS nM NM
		--chimSegmentMin 15   
    		--chimJunctionOverhangMin 15
                --readFilesIn  $input1.fastq $input2.fastq > $output.bam"""
}

addrgPIC = {
        output.dir="intermediate_files"
        exec """$picard AddOrReplaceReadGroups
                INPUT=$input.bam
                OUTPUT=$output.bam
                RGID=id
                RGLB=library
                RGPL=Illumina
                RGPU=machine
                RGSM=$input1.fastq.prefix"""
}

splitNtrimGATK = {
        output.dir="intermediate_files"
        exec """$GATK
                -T SplitNCigarReads
                -U ALLOW_N_CIGAR_READS
                -R $REF
                -I $input.bam
                -o $output.bam"""
}

countHTSeq = {
        output.dir="intermediate_files"
        exec "python -m HTSeq.scripts.count $input.bam $GTF --stranded=no --format=bam --order=pos > $output.txt"
}

callVariants_RNA = segment {
        alignSTAR +
        dedupPIC +
        splitNtrimGATK + 
        realignGATK +
	mpileupSAM + 
	somVARSCunpaired 
}
