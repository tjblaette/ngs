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
  # set the reference/base level for differential testing to a specific group within a column in the design table: (otherwise the first group in alphabetical order is chosen)
  # myddsHTSeq$condition <- relevel(myddsHTSeq$condition, "control")

  # should we filter out genes with basically no expression? YES
  myddsHTSeq <- myddsHTSeq[ rowSums(counts(myddsHTSeq)) > 1, ]

  # differential expression testing  
  mydds <- DESeq(myddsHTSeq)
  return(mydds)
}

mydds <- fromFile(input_file)


###########################################################################################
###########################################################################################
# PREPARE FOR EXPLORATORY ANALYSIS

# get normalized read counts
counts <- counts(mydds,normalized=TRUE)
write.table(counts, sep="\t", file=paste(input_file,"_DESeq2results_CountsNormalized.txt", sep=""))


# calculate coefficient of variation
standev <- apply(counts,1,sd)
avrg <- apply(counts,1,mean)
cv <- standev/avrg

# select genes with highest cv
select <- order(cv, decreasing=TRUE)[1:min(length(cv),30000)]

# transform counts and apply cv-based selection
trans <- rlog(mydds, blind=FALSE)
transCounts <- assay(trans)[select,]
write.table(assay(trans), sep="\t",file=paste(input_file,"_DESeq2results_CountsNormalizedTransformed.txt", sep=""))


# PREPARE COUNT CLUSTERING
# extract sample labels -> for columns get all factors for each sample
all_cols <- names(colData(mydds))
cols <- all_cols[!all_cols %in% c("sampleName", "fileName", "sizeFactor","replaceable")]
df <- as.data.frame(colData(mydds)[,cols])
rownames(df) <- rownames(colData(mydds))
colnames(df) <- cols

# PREPARE INTERSAMPLE DISTANCE CLUSTERING
# euclidean distance
sampleDists <- dist(t(transCounts))
sampleDistMatrix <- as.matrix(sampleDists)

# pearson correlation distance (taken from pheatmap source code)
sampleDists_corr <- as.dist(1 - cor(transCounts))

# to append sample info to sampleName for row labels:
#   rownames(sampleDistMatrix) <- paste(rownames(sampleDistMatrix), apply( df[ , cols ] , 1 , paste, collapse = "-" ), sep="-")



# PREPARE PCA
pca <- plotPCA(trans, intgroup="condition")
pca_full <- plotPCA(trans, intgroup=cols)


# PREPARE HEATMAP ANNOTATION COLORS (otherwise they do not match colors from PCA) -> hard-coded "condition" column!
# function to retrieve colors (evenly spaced on color ring)
ggplotColours <- function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

# retrieve colors for required number of groups to differentiate within "condition"
colors <- ggplotColours(length(unique(df$condition)))

# label colors with corresponding group
names(colors) <- sort(unique(df$condition))

# create a list to pass to pheatmap command
anno_colors <- list(condition=colors)




# PLOT
#1: complete linkage clustering based on Euclidean distance of transformed read counts -> based on top X genes with max CV
#2: complete linkage clustering based on Euclidean distance of Euclidean intersample distances -> based on top X genes with max CV
#3: PCA -> based on whole data set
library(pheatmap)

pdf(paste(input_file,"_DESeq2results_exploratory.pdf",sep=""))
print(pca)
if (length(cols) > 1)
{
  print(pca_full)
}
# to scale by row after clustering samples, provide distances calculated above -> rows are still clustered after scaling
pheatmap(transCounts, show_rownames=FALSE, treeheight_row=0, annotation_col=df, fontsize=7, scale="row", annotation_color=anno_colors, main="Clustered by Euclidean distance", clustering_distance_cols=sampleDists)
pheatmap(transCounts, show_rownames=FALSE, treeheight_row=0, annotation_col=df, fontsize=7, scale="row", annotation_color=anno_colors, main="Clustered by Pearson correlation", clustering_distance_cols=sampleDists_corr)
pheatmap(sampleDistMatrix, annotation_col=df, fontsize=7, annotation_color=anno_colors, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists)
dev.off()


###################################################################################################
# PROCESS RESULTS OF DGE ANALYSIS 
myresults <- results(mydds, alpha=my_alpha, altHypothesis="greaterAbs", lfcThreshold=my_lfc)


#sort according to adjusted p-value and write results to file
myresultsOrdered <- myresults[order(myresults$padj),]
write.table(myresultsOrdered, file=paste(input_file,"_DESeq2results.txt",sep=""),sep="\t")

#print summary of deg analysis
sink(paste(input_file,"_DESeq2results_summary.txt",sep=""))
summary(myresultsOrdered, alpha=my_alpha)
sink()

####################################################################################################
# ANNOTATE WITH GENE SYMBOL IN ADDITION TO ENSEMBL ID

#to convert ENSEMBL gene ids to gene symbols, I have to remove the decimal
rownames(myresultsOrdered) <- unlist(strsplit(rownames(myresultsOrdered), split='\\.'))[2*(1:length(rownames(myresultsOrdered)))-1]

#plot normalized gene counts to pdf
pdf(paste(input_file,"_DESeq2results_geneCountPlots.pdf",sep=""))
for (i in 1:(dim(mydds)[1]))
  {
    plotCounts(mydds, gene=i, intgroup="condition")
  }
dev.off()

#write deg analysis results to file
write.table(myresultsOrdered, file=paste(input_file,"_DESeq2results.txt",sep=""),sep="\t")

