#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


# requires as input first the *_Covered.bed file that defines genomic intervals covered by the design
# all additional inputs must be *coverBED.txt files of samples
# example call: Rscript $PATH/getCoveragePerGene.R /NGS/known_sites/hg19/haloplex_montse_27156-1457076787_Covered.bed *coverBED.txt




#map <- read.table("/NGS/known_sites/hg19/haloplex_montse_27156-1457076787_Covered.bed")
map <- read.table(args[1])
colnames(map) <- c("chr","start","end","gene")

genes <- unique(map$gene)
gene_def <- matrix(0,ncol=4,nrow=length(genes))
colnames(gene_def) <- c("gene","chr","start","end")


for(gene_ind in 1:length(genes))
{
  this_gene <- as.character(genes)[gene_ind]
  gene_def[gene_ind,1] <- this_gene
  gene_def[gene_ind,2] <- as.character(unique(map[map$gene == this_gene, 1]))
  gene_def[gene_ind,3] <- as.numeric(min(map[map$gene == this_gene,c(2,3)])) - 5
  gene_def[gene_ind,4] <- as.numeric(max(map[map$gene == this_gene,c(2,3)])) + 5
}

gene_def <- as.data.frame(gene_def)
gene_def$chr <- as.character(gene_def$chr)
gene_def$gene <- as.character(gene_def$gene)
gene_def$start <- as.numeric(as.character(gene_def$start))
gene_def$end <- as.numeric(as.character(gene_def$end))

summary <- matrix(0,ncol=length(args)-1,nrow=length(genes))
colnames(summary) <- args[2:length(args)]
rownames(summary) <- genes

for (file_index in 2:length(args))
{
input <- args[file_index]

#sample <- read.delim("example.txt", header=FALSE)
sample <- read.delim(input, header=FALSE)
colnames(sample) <- c("chr","start","end","amplicon","score","strand","cov","bases","lengthTotal","fraction")
sample <- sample[sample$chr != "all",]


sum <- 0

averageCov <- c()
dd <- c();

for(i in 1:(dim(gene_def)[1]))
{
  d <- sample[sample$chr == gene_def$chr[i] & gene_def$start[i] <= sample$start & gene_def$end[i] >= sample$end,]
  sum <- sum + dim(d)[1]
  
  totalBases <- sum(d$cov * d$bases)
  totalLength <- sum(unique(d[,c(4,9)])[2]) 

  averageCov <- c(averageCov, totalBases/totalLength)
}   

test <- sum == dim(sample)[1]

if(test)
{
  res <- data.frame(genes=genes,averageCoverage=averageCov)
  #write.table(res,file=paste("averageCoverage_",input,sep=""), sep='\t', row.names=FALSE)

  summary[,file_index-1] <- res[,2];
}

}

average <- apply(summary,1,mean)
summary <- cbind(summary,average)
#colnames(summary) <- c(colnames(summary),"average")
write.table(summary,file="averageCoverage_summary.txt", sep='\t')


