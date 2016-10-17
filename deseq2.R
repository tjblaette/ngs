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
save(mydds, file=paste(input_file,"_DESeq2results.RData", sep=""))
#load(paste(input_file,"_DESeq2results.RData", sep=""))

###########################################################################################
###########################################################################################
# PREPARE FOR EXPLORATORY ANALYSIS

# get normalized read counts
counts <- counts(mydds,normalized=TRUE)
write.table(counts, sep="\t", file=paste(input_file,"_DESeq2results_CountsNormalized.txt", sep=""))

# get transformed read counts
trans <- rlog(mydds)
save(trans, file=paste(input_file,"_DESeq2results_trans.RData", sep=""))
#load(file=paste(input_file,"_DESeq2results_trans.RData", sep=""))
write.table(assay(trans), sep="\t", file=paste(input_file,"_DESeq2results_CountsNormalizedTransformed.txt", sep=""))


# calculate coefficient of variation
standev <- apply(counts,1,sd)
avrg <- apply(counts,1,mean)
cv <- standev/avrg



# PREPARE PLOTTING
library(pheatmap)
pdf(paste(input_file,"_DESeq2results_exploratory.pdf",sep=""), height=10)

# extract sample labels -> for columns get all factors for each sample
all_cols <- names(colData(mydds))
cols <- all_cols[!all_cols %in% c("sampleName", "fileName", "sizeFactor","replaceable")]
df <- as.data.frame(colData(mydds)[,cols])
rownames(df) <- rownames(colData(mydds))
colnames(df) <- cols


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
color_genes <- c("forestgreen","red3")
names(color_genes) <- c("WT","mut")
color_CKs <- color_genes
names(color_CKs) <- c("nCK","CK")

# label colors with corresponding group
names(colors) <- sort(unique(df$condition))

# create a list to pass to pheatmap command
anno_colors <- list(condition=colors)
anno_colors <- list(condition=colors,ASXL1=color_genes,BCOR=color_genes,EZH2=color_genes,MLLptd=color_genes,RUNX1=color_genes,SF3B1=color_genes, STAG2=color_genes,U2AF1=color_genes,TP53=color_genes, CK=color_CKs)

# plot PCA -> if there is more than one annotation column, print a second PCA with all of that information
print(pca)
if (length(cols) > 1)
{
  print(pca_full)
}



for(maxGenes in c(50,100,500,1000,5000,10000,20000,30000))
{
 # select genes with highest cv
 select <- order(cv, decreasing=TRUE)[1:min(length(cv),maxGenes)]

 # transform counts and apply cv-based selection
 transCounts <- assay(trans)[select,]
 ##write.table(transCounts, sep="\t",file=paste(input_file,"_DESeq2results_CountsNormalizedTransformed_maxGenes_",maxGenes,".txt", sep=""))


 # PREPARE INTERSAMPLE DISTANCE CLUSTERING
 # euclidean distance
 sampleDists <- dist(t(transCounts))
 sampleDistMatrix <- as.matrix(sampleDists)

 # pearson correlation distance (taken from pheatmap source code)
 sampleDists_corr <- as.dist(1 - cor(transCounts))

 # to append sample info to sampleName for row labels:
 #   rownames(sampleDistMatrix) <- paste(rownames(sampleDistMatrix), apply( df[ , cols ] , 1 , paste, collapse = "-" ), sep="-")



 #1: complete linkage clustering based on Euclidean distance of transformed read counts -> based on top X genes with max CV, scaled by row
 #2: complete linkage clustering based on Pearson correlation of transformed read counts -> based on top X genes with max CV, scaled by row
 #3: complete linkage clustering based on Euclidean distance of Euclidean intersample distances -> based on top X genes with max CV, scaled by row
 # to scale by row after clustering samples, provide distances calculated above -> rows are still clustered after scaling
 pheatmap(transCounts, show_rownames=FALSE, treeheight_row=0, annotation_col=df, fontsize=5, scale="row", annotation_color=anno_colors, main=paste("Clustered by Euclidean distance - ",maxGenes,sep=""), clustering_distance_cols=sampleDists)
 pheatmap(transCounts, show_rownames=FALSE, treeheight_row=0, annotation_col=df, fontsize=5, scale="row", annotation_color=anno_colors, main=paste("Clustered by Pearson correlation - ",maxGenes,sep=""), clustering_distance_cols=sampleDists_corr)
 pheatmap(sampleDistMatrix, annotation_col=df, fontsize=5, annotation_color=anno_colors, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists)

} # end plotting

dev.off()


###################################################################################################
# PROCESS RESULTS OF DGE ANALYSIS 
myresults <- results(mydds, alpha=my_alpha, altHypothesis="greaterAbs", lfcThreshold=my_lfc)

# plot heatmap with DEGs only
sig <- which(myresults$padj < my_alpha)

# check if there are significant DEGs to plot
if(length(sig) > 1)
{
  sigCounts <- assay(trans)[sig, ]
print("should write degs counts") 
print(paste(input_file,"_DESeq2results_CountsNormalized_degs.txt", sep=""))
  write.table(sigCounts, sep="\t", file=paste(input_file,"_DESeq2results_CountsNormalized_degs.txt", sep=""))


  # euclidean distances
  sig_sampleDists <- dist(t(sigCounts))
  sig_sampleDistMatrix <- as.matrix(sig_sampleDists)

  # pearson correlation distances (taken from pheatmap source code)
  sig_sampleDists_corr <- as.dist(1 - cor(sigCounts))


  # prepare PCA
  sig_pca <- plotPCA(trans[sig,], intgroup="condition")
  sig_pca_full <- plotPCA(trans[sig,], intgroup=cols)

  # PLOT
  pdf(paste(input_file,"_DESeq2results_degs.pdf",sep=""), height=10)
  print(sig_pca) 
  if (length(cols) > 1)
  {
    print(sig_pca_full)
  }

  pheatmap(sigCounts, show_rownames=FALSE, treeheight_row=0, annotation_col=df, fontsize=5, scale="row", annotation_color=anno_colors, main="DEGs clustered by Euclidean distance", clustering_distance_cols=sig_sampleDists)
  pheatmap(sigCounts, show_rownames=FALSE, treeheight_row=0, annotation_col=df, fontsize=5, scale="row", annotation_color=anno_colors, main="DEGs clustered by Pearson correlation", clustering_distance_cols=sig_sampleDists_corr)
  pheatmap(sig_sampleDistMatrix, annotation_col=df, fontsize=5, annotation_color=anno_colors, clustering_distance_rows=sig_sampleDists, clustering_distance_cols=sig_sampleDists)
  dev.off()
}


#sort according to adjusted p-value and write results to file
myresultsOrdered <- myresults[order(myresults$padj),]
write.table(myresultsOrdered, file=paste(input_file,"_DESeq2results.txt",sep=""),sep="\t")

#print summary of deg analysis
sink(paste(input_file,"_DESeq2results_summary.txt",sep=""))
summary(myresultsOrdered, alpha=my_alpha)
sink()

####################################################################################################
# ANNOTATE WITH GENE SYMBOL IN ADDITION TO ENSEMBL ID

#plot normalized gene counts to pdf
pdf(paste(input_file,"_DESeq2results_geneCountPlots.pdf",sep=""))
for (i in 1:(dim(mydds)[1]))
  {
    plotCounts(mydds, gene=i, intgroup="condition")
  }
dev.off()

## PRINT SESSION INFO
sink(paste(input_file,"_sessionInfo.txt",sep=""))
print(sessionInfo())
sink()

#to convert ENSEMBL gene ids to gene symbols later, I have to remove the decimal now
rownames(myresultsOrdered) <- unlist(strsplit(rownames(myresultsOrdered), split='\\.'))[2*(1:length(rownames(myresultsOrdered)))-1]

#write deg analysis results to file -> overwrite if decimal removal was successful (will fail for Susi's circRNAs) and keep existing version otherwise
write.table(myresultsOrdered, file=paste(input_file,"_DESeq2results.txt",sep=""),sep="\t")

