#!/bin/R

# this script takes as input a tab-separated file containing mutational data for multiple samples
# table should contain one column for sample names followed by one additional column for each gene described
# mutations should be marked with "1", WT genotypes with "0" (but without quotation marks!)

# example input file:
#sampleName	ASXL1	BCOR	EZH2	MLLptd	RUNX1	SF3B1	U2AF1
#2805	0	0	0	0	1	0	0
#2813	0	0	0	0	0	0	0
#2814	0	0	0	0	0	0	0
#2817	0	0	1	0	0	0	0
#2820	0	0	0	0	0	0	0
#2821	1	0	0	0	1	0	1
#2830	0	0	0	1	0	0	0
#2843	0	0	0	0	0	0	1


library(ggplot2)
library(reshape2)
library(ggthemes)

args <- commandArgs(TRUE)
print(args)
input_file <- args[1]
input_file_prefix <- tools::file_path_sans_ext(input_file)

mm <- as.matrix(read.table(input_file, header=TRUE, row.names=1, na.string=0))


mutcounts <- apply(mm,2,sum,na.rm=T)
mms <- mm[,(order(mutcounts))]

mms <- transform(data.frame(mms))
mms <- mms[ do.call(order,rev(mms)), ]



mmsd <- data.frame(mms)
mmsdpat <- data.frame(mms,patients=rownames(mms))
mmsdpat$patients <- factor(mmsdpat$patients, levels=mmsdpat$patients)

melted <- melt(mmsdpat)
pdf(paste(input_file_prefix,"_mutationHeatmmap.pdf", sep=""))
ggplot(data = melted, aes(patients,variable, fill=value*as.numeric(variable))) + geom_tile(color="grey60", show.legend=FALSE) + scale_fill_gradient(na.value = "white", low = "#CCCCFF", high="#333399") + theme_tufte(base_family="Helvetica") + ylab(NULL) + xlab("Samples") + ggtitle("Mutational landscape of TCGA AML subset") + theme(axis.ticks=element_blank(), axis.text.x=element_blank())
dev.off()



