
args <- commandArgs(TRUE)
IN <- args[1] # DESeq2results_annotated_woutNA.txt
OUT <- args[2] # output file prefix
ORGANISM <- args[3] # Hs for human 

# IN <- "/media/LBseq/tamara/countCircs/forDESeq2.tsv_DESeq2results.txt"
if(!exists("ORGANISM"))
{
  ORGANISM="Hs"
}
print(ORGANISM)

library("topGO")
DATA <- read.table(IN, header=TRUE)

# topGO requires a named vector -> values are p-values, names are gene IDs
geneList <- DATA$padj
names(geneList) <- DATA$geneID
#names(geneList) <- rownames(DATA)


# function to select genes of interest for count-based tests (like classical fisher's exact test)
topDiffGenes <- function(vec)
{
  return(vec <= 0.05)
}

# define algorithms and statistical methods to be used for enrichment analysis -> match based on position on respective vector!
algos <- c("classic", "weight01")
stats <- c("fisher", "ks")


pdf(paste(OUT,"_topGO.pdf", sep=""))
for(method in 1:(length(algos)))
{
  for(onto in c("CC","MF","BP"))
  {
    cat("\n", "\n", "\n")
    cat("======================================================================", "\n")
    cat("======================================================================", "\n", "\n", "\n")
    cat("ENRICHMENT of ", onto, " GO terms based on ", algos[method], stats[method], " test", "\n")

    # create a topGO data object for analysis
    goData <- new("topGOdata", 
              ontology=onto,
              allGenes=geneList,
              geneSel = topDiffGenes,  # requires max p-value of 0.05
              nodeSize=10,
              annot=annFUN.org, mapping = paste("org.",ORGANISM,".eg.db", sep=""), ID="ensembl")
              #annot=annFUN.org, mapping = "org.Hs.eg.db", ID="ensembl")


    # perform the enrichment test
    result <- runTest(goData, algorithm=algos[method], statistic=stats[method])

    # collect some statistics on nr of significant terms 
    cat("\n")
    cat("Number of significant terms prior to adjusting for mutliple testing (p < 0.05): ", sum(score(result) < 0.05), "\n")
    candidates <- sum(score(result) < 0.1) # nr of most sign terms to save to disk later

    # adjust for multiple testing -> this is not recommended by authors of topGO, especialy for weight and elim algos!
    #score(result) <- p.adjust(score(result), method="BH")
    #sigN <- sum(score(result) < 0.05)
    #cat("Number of significant terms after adjusting for mutliple testing (p < 0.05, BH): ", sigN, "\n", "\n")

    # save all candidates to a file...
    longSumm <- GenTable(goData, pValue=result, orderBy="pValue", topNodes=candidates, numChar=500)
    write.table(longSumm, sep="\t", row.names=FALSE, file=paste(c(paste(c(OUT,"topGO",onto,algos[method], stats[method]),collapse="_"),".txt"), collapse=""))

    #...and output the most interesting ones to stdout
    summ <- GenTable(goData, pValue=result, orderBy="pValue", topNodes=min(max(c(10,candidates)),50))
    cat("\n","Summary of top ", min(max(c(10,candidates)),50), " terms", "\n", "\n")
    print(summ)

    # plot the most interesting terms in the GO hierarchy
    showSigOfNodes(goData, score(result), firstSigNodes=min(max(c(5,candidates)),7), useInfo="all"); title(main=paste(c("Top ", min(max(c(5,candidates)),7) , " ", onto, " terms based on ", algos[method], " ", stats[method], " test"), collapse=""), cex.main=0.65);
  }
}
dev.off()

