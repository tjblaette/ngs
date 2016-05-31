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


mydds_all <- fromFile(input_file)

# open pdf for plotting
pdf(paste(input_file,"_DESeq2results_exploratory.pdf",sep=""))

for (blinded in c(TRUE))#,TRUE))
{
print(blinded)
 # should we filter out genes with basically no expression? YES
 for (minRowSum in c(100))#1,2,5,10,50,100))
 {
  mydds <- mydds_all[ rowSums(counts(mydds_all)) >= minRowSum, ]
print(minRowSum)

  ###########################################################################################
  ###########################################################################################
  # PREPARE FOR EXPLORATORY ANALYSIS

  # get normalized read counts
  counts <- counts(mydds,normalized=TRUE)
  #write.table(counts, file="normlized_counts.txt")

  # calculate coefficient of variation
  standev <- apply(counts,1,sd)
  avrg <- apply(counts,1,mean)
  cv <- standev/avrg

  plot(x=1:length(cv), y=rowSums(counts(mydds))[order(cv)], main="Rowsums ranked acc. to CV")
  plot(x=1:length(cv), y=rowSums(counts(mydds))[order(standev*standev)], main="Rowsums ranked acc. to variance")

  for (max_genes in c(500))#,1000,5000,10000,15000,20000,25000,30000))
  { 
print(max_genes)
   # select genes with highest cv
   num_genes <- min(length(cv),max_genes)
   select <- order(cv, decreasing=TRUE)[1:num_genes]

   # transform counts and apply cv-based selection
   trans <- rlog(mydds, blind=blinded)
   transCounts <- assay(trans)[select,]
   #write.table(transCounts, file="normlized_counts_transformed.txt")


   # PREPARE COUNT CLUSTERING
   # extract sample labels
   #df <- as.data.frame(colData(mydds)[,c("condition","sex")])
   all_cols <- names(colData(mydds))
   cols <- all_cols[!all_cols %in% c("sampleName", "fileName", "sizeFactor","replaceable")]
   df <- as.data.frame(colData(mydds)[, cols])


   # PREPARE INTERSAMPLE DISTANCE CLUSTERING
   sampleDists <- dist(t(transCounts))
   sampleDistMatrix <- as.matrix(sampleDists)
   # to append sample info to sampleName for row labels:
#   rownames(sampleDistMatrix) <- paste(rownames(sampleDistMatrix), apply( df[ , cols ] , 1 , paste, collapse = "-" ), sep="-")

   # PREPARE PCA
   pca <- plotPCA(trans, intgroup="condition")
   pca_full <- plotPCA(trans, intgroup=cols)

   # PLOT
   #1: complete linkage clustering based on Euclidean distance of transformed read counts -> based on top X genes with max CV
   #2: complete linkage clustering based on Euclidean distance of Euclidean intersample distances -> based on top X genes with max CV
   #3: PCA -> based on whole data set except for minRowSum filter
   library(pheatmap)

   print(pca)
   print(pca_full)
   pheatmap(transCounts, show_rownames=FALSE, treeheight_row=0, annotation_col=df, fontsize=7, main=paste("num_genes: ", num_genes, ", minRowSum: ", minRowSum, ", blind: ", blinded, ", distance= euclidean", sep=''))
   pheatmap(transCounts, show_rownames=FALSE, treeheight_row=0, annotation_col=df, fontsize=7, clustering_distance_cols="correlation", main=paste("num_genes: ", num_genes, ", minRowSum: ", minRowSum, ", blind: ", blinded, ", distance= pearson", sep=''))
   pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, annotation_col=df, fontsize=7, main=paste("num_genes: ", num_genes, ", minRowSum: ", minRowSum, ", blind: ", blinded, ", distance= euclidean", sep=''))
  }
 }
} 

# close pdf
dev.off()

