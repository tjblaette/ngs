args <- commandArgs(TRUE)
print(args)
input_file <- args[1]
my_alpha <- args[2]
my_lfc <- args[3]

library("DESeq2")

# DOES NOT WORK WITH ARGUMENT PASSED!
my_lfc <- 0.6

fromFile <- function(input) {
  myTable <- read.table(input, header=TRUE)
  myDir <- "intermediate_files"

  myddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=myTable, directory=myDir, design=~ condition)
  mydds <- DESeq(myddsHTSeq)
 # myresults <- results(mydds)

  return(mydds)
}


mydds <- fromFile(input_file)
myresults <- results(mydds, alpha=my_alpha, altHypothesis="greaterAbs", lfcThreshold=my_lfc)


#sort according to adjusted p-value
myresultsOrdered <- myresults[order(myresults$padj),]

#to convert ENSEMBL gene ids to gene symbols, I have to remove the decimal
rownames(myresultsOrdered) <- unlist(strsplit(rownames(myresultsOrdered), split='\\.'))[2*(1:length(rownames(myresultsOrdered)))-1]


#system(paste("head -1", f," > blub_sorted"))
#system(paste("tail -n +2", f," | sort >> blub_sorted"))


##This annotation gives different results from ENSEMBLE's biomart or browser queries so I am not using it!
#library("AnnotationDbi")
#library("org.Hs.eg.db")
##add gene symbol annotation
#myresultsOrdered$geneSymbol <- mapIds(org.Hs.eg.db, keys=row.names(myresultsOrdered), column="SYMBOL", keytype="ENSEMBL",multiVals="first")

#print summary of deg analysis
sink(paste(input_file,"_DESeq2results_summary.txt",sep=""))
summary(myresultsOrdered, alpha=my_alpha)
sink()

#save normalized gene counts to text file
countsTable <- c();

#plot normalized gene counts to pdf
pdf(paste(input_file,"_DESeq2results_geneCountPlots.pdf",sep=""))
  for (i in 1:(dim(mydds)[1]))
  {
    plotCounts(mydds, gene=i, intgroup="condition")
    countsTable <- rbind(countsTable,plotCounts(mydds, gene=i, intgroup="condition", returnData=TRUE)[,1]) 
    colnames(countsTable)[i] <- rownames(mydds)[i]
  }
dev.off()

#write deg analysis results to file
write.table(myresultsOrdered, file=paste(input_file,"_DESeq2results.txt",sep=""))
write.table(countsTable, file=paste(input_file,"_DESeq2results_CountsTable.txt",sep=""))




