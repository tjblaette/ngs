
####
# T.J.Blätte
# 2016
####
#
# Performs GO term enrichment analysis using
#       the topGO R package.
#
# Args:
#   IN: Name of / Path to the DESeq2 results
#       file, with NA-containing rows filtered
#       (*woutNA.txt)
#   OUT: Output file prefix
#   ORGANISM: Abbreviation for the organisms
#       from which the data of IN was derived.
#       Default is "Hs" for human.
#
# Output:
#   ${OUT}_topGO.pdf
#
# Requires:
#   topGO, stringr, Rgraphviz
#
####

args <- commandArgs(TRUE)
IN <- args[1]
OUT <- args[2]
ORGANISM <- if (length(args) < 3) "Hs" else args[3]

library("topGO")
library(stringr)
library(Rgraphviz)

data <- read.table(IN, header=TRUE)
if ('geneID_geneSymbol' %in% colnames(data))
{
    data$geneID <- str_split_fixed(data$geneID_geneSymbol, "_", 2)[,1]
    data$geneSymbol <- str_split_fixed(data$geneID_geneSymbol, "_", 2)[,2]
}


# topGO requires a named vector -> values are p-values, names are gene IDs
geneList <- data$padj
names(geneList) <- data$geneID


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
                      annot=annFUN.org,
                      mapping = paste("org.",ORGANISM,".eg.db", sep=""), ID="ensembl")

        # perform the enrichment test
        result <- runTest(goData,
                          algorithm=algos[method],
                          statistic=stats[method])

        # collect some statistics on nr of significant terms
        cat("\n",
            "Number of significant terms prior to adjusting for mutliple testing (p < 0.05): ",
            sum(score(result) < 0.05),
            "\n")
        candidates <- sum(score(result) < 0.1) # nr of most sign terms to save to disk later

        # save all candidates to a file...
        longSumm <- GenTable(goData,
                             pValue=result,
                             orderBy="pValue",
                             topNodes=candidates,
                             numChar=500)
        write.table(longSumm,
                    sep="\t",
                    row.names=FALSE,
                    file=paste(c(paste(c(OUT,"topGO",onto,algos[method], stats[method]),collapse="_"),".txt"), collapse=""))

        #...and output the most interesting ones to stdout
        summ <- GenTable(goData,
                         pValue=result,
                         orderBy="pValue",
                         topNodes=min(max(c(10,candidates)),50))
        cat("\n",
            "Summary of top ",
            min(max(c(10,candidates)),50),
            " terms",
            "\n",
            "\n")
        print(summ)

        # plot the most interesting terms in the GO hierarchy
        showSigOfNodes(goData,
                       score(result),
                       firstSigNodes=min(max(c(5,candidates)),7),
                       useInfo="all")
        title(main=paste(c("Top ", min(max(c(5,candidates)),7), " ", onto, " terms based on ", algos[method], " ", stats[method], " test"), collapse=""),
              cex.main=0.65)
    }
}
dev.off()
