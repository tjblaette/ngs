//file to contain pipeline stages of all RNA-seq pipelines

//// COMMON VARIABLE DEFINITIONS
// PROGRAM DIRECTORIES
STAR="star"
STAR_LONG="/NGS/star_long/source/STARlong"


qc = {
    var nkern : 24
    // subsample FASTQ files for QC -> take every 25th read only -> (NR = read fraction ^-1 * 4 = 25 * 4 = 100)
    exec "awk 'NR % 100 > 0 && NR % 100 < 5' $input1.fastq > $output1.subsample"
    exec "awk 'NR % 100 > 0 && NR % 100 < 5' $input2.fastq > $output2.subsample"
    exec "$QC -pe $output1.subsample $output2.subsample N A -c $nkern -onlyStat -o FASTQ_QC && touch $output.success"
}


alignSTAR = {
        var nkern : 48 
        output.dir="intermediate_files"
        exec """$STAR
                --runThreadN $nkern
                --genomeDir $REFDIR 
                --twopassMode Basic
		--quantMode GeneCounts
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


alignSTARlong = {
        var nkern : 48 
        output.dir="intermediate_files"
	//	--quantMode GeneCounts
	//	--chimSegmentMin 15   
    	//	--chimJunctionOverhangMin 15
        //        --twopassMode Basic
        exec """$STAR_LONG
                --runThreadN $nkern
                --genomeDir $REFDIR 
                --outStd BAM_SortedByCoordinate
                --outSAMtype BAM SortedByCoordinate
                --outFileNamePrefix ./.$input.fastq.prefix
                --outSAMunmapped Within
                --outSAMmapqUnique 60 
                --outSAMattrRGline ID:0 PL:ILLUMINA LB:rna SM:$input.fastq.prefix
                --outSAMstrandField intronMotif 
                --outSAMattributes NH HI AS nM NM
		--outFilterMultimapScoreRange 20   
		--outFilterScoreMinOverLread 0  
		--outFilterMatchNminOverLread 0.66   
		--outFilterMismatchNmax 1000   
		--winAnchorMultimapNmax 200   
		--seedSearchLmax 12   
		--seedSearchStartLmax 12   
		--seedPerReadNmax 100000   
		--seedPerWindowNmax 100   
		--alignTranscriptsPerReadNmax 100000   
		--alignTranscriptsPerWindowNmax 10000
		--alignIntronMax 1
                --readFilesIn  $input.fastq > $output.bam"""
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

countHTSeq_unstranded = {
        output.dir="intermediate_files"
        exec "python -m HTSeq.scripts.count $input.bam $GTF --stranded=no --format=bam --order=pos > $output.counts"
}

countHTSeq_sameStranded = {
        output.dir="intermediate_files"
        exec "python -m HTSeq.scripts.count $input.bam $GTF --stranded=yes --format=bam --order=pos > $output.counts"
}

countHTSeq_otherStranded = {
        output.dir="intermediate_files"
        exec "python -m HTSeq.scripts.count $input.bam $GTF --stranded=reverse --format=bam --order=pos > $output.counts"
}
