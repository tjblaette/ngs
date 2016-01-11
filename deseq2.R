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

###########################################################################################
###########################################################################################
# PREPARE FOR EXPLORATORY ANALYSIS

# should we filter out genes with 0 expression?
# mydds <- mydds[ rowSums(counts(mydds)) > 0, ]

# get normalized read counts
counts <- counts(mydds,normalized=TRUE)

# calculate coefficient of variation
standev <- apply(counts,1,sd)
avrg <- apply(counts,1,mean)
cv <- standev/avrg

# select genes with highest cv
select <- order(cv, decreasing=TRUE)[1:min(length(cv),30000)]

# transform counts and apply cv-based selection
trans <- rlog(mydds)
transCounts <- assay(trans)[select,]


# PREPARE COUNT CLUSTERING
# extract sample labels
df <- as.data.frame(colData(mydds)[,"condition"])
rownames(df) <- colnames(mydds)
colnames(df) <- "condition"


# PREPARE INTERSAMPLE DISTANCE CLUSTERING
sampleDists <- dist(t(transCounts))
sampleDistMatrix <- as.matrix(sampleDists)


# PREPARE PCA
pca <- plotPCA(trans, intgroup="condition")

# PLOT
#1: complete linkage clustering based on Euclidean distance of transformed read counts -> based on top X genes with max CV
#2: complete linkage clustering based on Euclidean distance of Euclidean intersample distances -> based on top X genes with max CV
#3: PCA -> based on whole data set
library(pheatmap)

pdf(paste(input_file,"_DESeq2results_exploratory.pdf",sep=""))
print(pca)
pheatmap(transCounts, show_rownames=FALSE, treeheight_row=0, annotation_col=df)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists)
dev.off()


###################################################################################################
# PROCESS RESULTS OF DGE ANALYSIS 
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
  }
dev.off()

colnames(countsTable) <- colnames(mydds)
rownames(countsTable) <- rownames(mydds)

#write deg analysis results to file
write.table(myresultsOrdered, file=paste(input_file,"_DESeq2results.txt",sep=""),sep="\t")
write.table(countsTable, file=paste(input_file,"_DESeq2results_CountsTable.txt",sep=""),sep="\t")

