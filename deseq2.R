args <- commandArgs(TRUE)
print(args)
input_file <- args[1]
my_design <- args[2]
my_reference_level <- args[3] # factor level to use as the reference for calculated fold changes
my_alpha <- as.numeric(args[4])
my_lfc <- as.numeric(args[5])
output_prefix <- args[6]

library(DESeq2)
library(pheatmap)
library(ggplot2)


design_str_to_vector_of_str_terms <- function(design_str) {
    vector_of_terms <- attr(terms(formula(design_str)), "term.labels")
    return (vector_of_terms)
}

remove_interaction_terms <- function(design_terms) {
    design_terms <- design_terms[!grepl(":", design_terms)]
    return (design_terms)
}


fromFile <- function(input) {
    myDir <- "intermediate_files"
    myTable <- read.table(input, header=TRUE)
    condition_to_test <- rev(remove_interaction_terms(design_str_to_vector_of_str_terms(my_design)))[1]

    # patient must be a factor, even if ID is numeric
    if ("patient" %in% colnames(myTable))
    {
        myTable[["patient"]] = factor(myTable[["patient"]])
    }

    myddsHTSeq <- DESeqDataSetFromHTSeqCount(
            sampleTable=myTable,
            directory=myDir,
            design= formula(my_design))
    # explicitely set the reference/base level for differential testing
    # --> otherwise the first group in alphabetical order is chosen
    myddsHTSeq[[condition_to_test]] <- relevel(myddsHTSeq[[condition_to_test]], my_reference_level)

    # filter genes basically not expressed
    # --> will be filtered by DESeq2 anyway - doing it here will speed up analysis
    myddsHTSeq <- myddsHTSeq[ rowSums(counts(myddsHTSeq)) > 1, ]

    # differential expression testing
    mydds <- DESeq(myddsHTSeq)
    return(mydds)
}

mydds <- fromFile(input_file)
my_terms_of_interest <- remove_interaction_terms(design_str_to_vector_of_str_terms(my_design))
save(mydds, file=paste(output_prefix,".RData", sep=""))
#load(paste(output_prefix,".RData", sep=""))

###########################################################################################
###########################################################################################
# PREPARE FOR EXPLORATORY ANALYSIS

# plot sparsity
# plot of the concentration of counts in a single sample over the sum of counts per gene.
# --> useful diagnostic for datasets which might not fit a negative binomial assumption:
#     genes with many zeros and individual very large counts are difficult to model with
#     the negative binomial distribution.
pdf(paste(output_prefix,"_sparsity.pdf",sep=""), height=10)
plotSparsity(mydds)
dev.off()


# get normalized read counts
counts <- counts(mydds,normalized=TRUE)
write.table(
        counts,
        sep="\t",
        file=paste(output_prefix,"_CountsNormalized.txt", sep=""))

# get transformed read counts
trans <- rlog(mydds)
save(trans, file=paste(output_prefix,"_trans.RData", sep=""))
#load(file=paste(output_prefix,"_trans.RData", sep=""))
write.table(
        assay(trans),
        sep="\t",
        file=paste(output_prefix,"_CountsNormalizedTransformed.txt", sep=""))


# calculate coefficient of variation
standev <- apply(counts,1,sd)
avrg <- apply(counts,1,mean)
cv <- standev/avrg



# PREPARE PLOTTING
pdf(paste(output_prefix,"_exploratory.pdf",sep=""), height=10)

# extract sample labels -> for columns get all factors for each sample
all_cols <- names(colData(mydds))
cols <- all_cols[!all_cols %in% c("sampleName", "fileName", "sizeFactor","replaceable")]
df <- as.data.frame(colData(mydds)[,cols])
rownames(df) <- rownames(colData(mydds))
colnames(df) <- cols


# PREPARE PCA
# PREPARE HEATMAP ANNOTATION COLORS (otherwise they do not match colors from PCA)
# function to retrieve colors (evenly spaced on color ring)
ggplotColours <- function(n=6, h=c(0, 360) +15){
    if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

# function to prepare data for PCA plotting (based on DESeq2 code)
my_prepPCA <- function (object, xPC, yPC, intgroup = "condition", ntop = 500, returnData = FALSE, output_prefix=output_prefix)
{
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    write.table(pca["rotation"], file=paste(output_prefix, "_pca.txt", sep=""), sep="\t")
    percentVar <- pca$sdev^2/sum(pca$sdev^2)

    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }

    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
    if (length(intgroup) > 1) {
        group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
    } else {
        group <- colData(object)[[intgroup]]
    }

    d <- data.frame(PC1=pca$x[, xPC], PC2=pca$x[, yPC], group=group, intgroup.df, colnames(object))
    if (returnData) {
        attr(d, "percentVar") <- percentVar[xPC:yPC]
        return(d)
    }
}

# function to plot PCA (based on DESeq2 code)
my_plotPCA <- function(data, intgroup, xPC, yPC, xPC_label, yPC_label, pass_outputprefix=output_prefix)
{
    pcaData <- my_prepPCA(data, xPC, yPC, intgroup=intgroup, returnData=TRUE, output_prefix=pass_outputprefix)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes_(as.name("PC1"), as.name("PC2"), color = as.name(intgroup))) +
            geom_point(size = 3) +
            xlab(paste0(xPC_label, percentVar[1], "% variance")) +
            ylab(paste0(yPC_label, percentVar[2], "% variance")) +
            labs(color = intgroup) +
            coord_fixed()
}

# function to pretty print gene counts, based on DESeq2 vignette code
#   --> possibly very slow, try to speed things up
my_plotCounts <- function(mydds, gene, design, terms_of_interest)
{
    plot_data <- plotCounts(
            mydds,
            gene=i,
            xlab=design,
            intgroup=terms_of_interest,
            replaced=("replaceCounts" %in% names(assays(mydds))),
            returnData=TRUE)

    if (length(terms_of_interest) == 1)
    {
        print(ggplot(plot_data, aes_string(x=terms_of_interest[1], y="count")) +
                geom_jitter(size=1.5, position = position_jitter(width=.15)) +
                stat_summary(fun.y=mean, geom="line", colour="red", size=0.8, aes(group=terms_of_interest[1])) +
                xlab(terms_of_interest[1]) +
                ylab("normalized counts") +
                ggtitle(rownames(mydds[i])))
    } else {
        if (length(terms_of_interest) == 2)
        {
            ggplot_data <- ggplot(plot_data, aes_string(x=terms_of_interest[2], y="count", group=terms_of_interest[1])) +
                    facet_wrap(as.formula(paste("~", terms_of_interest[1]))) +
                    stat_summary(fun.y=mean, geom="line", colour="red", size=0.8) +
                    xlab(terms_of_interest[1]) +
                    ylab("normalized counts") +
                    ggtitle(rownames(mydds[i]))
            # if this is a paired analysis, i.e. each facet contains 1 data point per condition, do not jitter
            # --> or the other way around: jitter only, if this is NOT a paired analysis,
            #           i.e. if there is more than data point in plot_data per combination of terms_of_interest
            if (length(unique(plot_data[[terms_of_interest[1]]])) * length(unique(plot_data[[terms_of_interest[2]]]))  < dim(plot_data)[1])
            {
                ggplot_data <- ggplot_data +
                        geom_jitter(size=1.5, position = position_jitter(width=.15))
            }
            print(ggplot_data)
        } else {
            # for designs with more than 2 terms, use default plots
            plotCounts(
                    mydds,
                    gene=i,
                    xlab=design,
                    intgroup=terms_of_interest,
                    replaced=("replaceCounts" %in% names(assays(mydds))))
        }
    }
}

# actually plot PCA, color-coding each group in annotation-df separately
# collect colors used to pass to pheatmap, to match colors between heatmap and PCA
anno_colors <- list()
for (col in cols)
{
    # retrieve colors for required number of groups to differentiate within col
    # colors <- ggplotColours(length(unique(na.omit(df[[col]]))))
    # retrieve colors used by PCA (PCA command copied from below!)
    colors <- unique(ggplot_build(my_plotPCA(trans, col, 1, 2, "PC1: ", "PC2: "))$data[[1]][["colour"]])

    # label colors with corresponding group
    names(colors) <- unique(df[[col]])
    #names(colors) <- unique(na.omit(df[[col]]))
    colors <- colors[!sapply(names(colors), is.na)]

    # sort both
    colors <- colors[order(names(colors))]
    names(colors) <- sort(names(colors))

    # add to list of colors later passed to pheatmap
    anno_colors[[col]] <- colors

    # plot PCA -> if there is more than one annotation column, print a second PCA with all of that information
    print(my_plotPCA(trans, col, 1, 2, "PC1: ", "PC2: "))
    print(my_plotPCA(trans, col, 3, 4, "PC3: ", "PC4: "))
    print(my_plotPCA(trans, col, 5, 6, "PC5: ", "PC6: "))
}

# this is not working! printing one of cols only!
#if (length(cols) > 1)
#{
#  print(my_plotPCA(trans, cols, 1, 2, "PC1: ", "PC2: "))
#  print(my_plotPCA(trans, cols, 3, 4, "PC3: ", "PC4: "))
#  print(my_plotPCA(trans, cols, 5, 6, "PC5: ", "PC6: "))
#}



for(maxGenes in c(50,100,500,1000,5000,10000,20000,30000))
{
    # select genes with highest cv
    select <- order(cv, decreasing=TRUE)[1:min(length(cv),maxGenes)]

    # transform counts and apply cv-based selection
    transCounts <- assay(trans)[select,]

    # PREPARE INTERSAMPLE DISTANCE CLUSTERING
    # euclidean distance
    sampleDists <- dist(t(transCounts))
    sampleDistMatrix <- as.matrix(sampleDists)

    # pearson correlation distance (taken from pheatmap source code)
    sampleDists_corr <- as.dist(1 - cor(transCounts))


    #1: complete linkage clustering based on Euclidean distance of transformed read counts -> based on top X genes with max CV, scaled by row
    #2: complete linkage clustering based on Pearson correlation of transformed read counts -> based on top X genes with max CV, scaled by row
    #3: complete linkage clustering based on Euclidean distance of Euclidean intersample distances -> based on top X genes with max CV, scaled by row
    # to scale by row after clustering samples, provide distances calculated above -> rows are still clustered after scaling
    pheatmap(
            transCounts,
            annotation_col=df,
            annotation_color=anno_colors,
            clustering_distance_cols=sampleDists,
            main=paste("Clustered by Euclidean distance of ",maxGenes," top CV genes", sep=""),
            scale="row",
            show_rownames=FALSE,
            treeheight_row=0,
            fontsize=5,
            border_color=NA)
    pheatmap(
            transCounts,
            annotation_col=df,
            annotation_color=anno_colors,
            clustering_distance_cols=sampleDists_corr,
            main=paste("Clustered by Pearson correlation of ",maxGenes," top CV genes", sep=""),
            scale="row",
            show_rownames=FALSE,
            treeheight_row=0,
            fontsize=5,
            border_color=NA)
    pheatmap(
            sampleDistMatrix,
            annotation_col=df,
            annotation_color=anno_colors,
            clustering_distance_rows=sampleDists,
            clustering_distance_cols=sampleDists,
            main=paste("Euclidean inter-sample distance based on ",maxGenes," top CV genes", sep=""),
            fontsize=5,
            border_color=NA)
}

dev.off()




###########################################################################################
###########################################################################################
# PROCESS RESULTS OF DGE ANALYSIS


myresults <- results(mydds, alpha=my_alpha, altHypothesis="greaterAbs", lfcThreshold=my_lfc)

# plot heatmap with DEGs only
sig <- which(myresults$padj < my_alpha)

# check if there are significant DEGs to plot
if(length(sig) > 1)
{
    sigCounts <- assay(trans)[sig, ]
    write.table(
            sigCounts,
            sep="\t",
            file=paste(output_prefix,"_CountsNormalizedTransformed_degs.txt", sep=""))

    # euclidean distances
    sig_sampleDists <- dist(t(sigCounts))
    sig_sampleDistMatrix <- as.matrix(sig_sampleDists)

    # pearson correlation distances (taken from pheatmap source code)
    sig_sampleDists_corr <- as.dist(1 - cor(sigCounts))

    # plot PCAs, again color-coding each of the annotations separately
    pdf(paste(output_prefix,"_degs.pdf",sep=""), height=10)
    for (col in cols)
    {
        print(my_plotPCA(trans[sig, ], col, 1, 2, "PC1: ", "PC2: "))
        print(my_plotPCA(trans[sig, ], col, 3, 4, "PC3: ", "PC4: "))
        print(my_plotPCA(trans[sig, ], col, 5, 6, "PC5: ", "PC6: "))
    }

    #  if (length(cols) > 1)
    #  {
    #    #sig_pca_full <- plotPCA(trans[sig,], intgroup=cols)
    #    #print(sig_pca_full)
    #    print(my_plotPCA(trans[sig, ], cols, 1, 2, "PC1: ", "PC2: "))
    #    print(my_plotPCA(trans[sig, ], cols, 3, 4, "PC3: ", "PC4: "))
    #    print(my_plotPCA(trans[sig, ], cols, 5, 6, "PC5: ", "PC6: "))
    #  }

    # plot the same heatmaps as above, now for DEGs only
    pheatmap(
            sigCounts,
            annotation_col=df,
            annotation_color=anno_colors,
            clustering_distance_cols=sig_sampleDists,
            main="Clustered by Euclidean distance of DEGs",
            scale="row",
            show_rownames=length(sig)<=50,
            treeheight_row=0,
            fontsize=5,
            border_color=NA)
    pheatmap(
            sigCounts,
            annotation_col=df,
            annotation_color=anno_colors,
            clustering_distance_cols=sig_sampleDists_corr,
            main="Clustered by Pearson correlation of DEGs",
            scale="row",
            show_rownames=length(sig)<=50,
            treeheight_row=0,
            fontsize=5,
            border_color=NA)
    pheatmap(
            sig_sampleDistMatrix,
            annotation_col=df,
            annotation_color=anno_colors,
            clustering_distance_rows=sig_sampleDists,
            clustering_distance_cols=sig_sampleDists,
            main="Euclidean inter-sample distance based on DEGs",
            fontsize=5,
            border_color=NA)
    dev.off()

    # plot DEG counts
    pdf(paste(output_prefix,"_geneCountPlots_degs.pdf",sep=""))
    for (i in sig)
    {
        my_plotCounts(mydds, gene=i, design=my_design, terms_of_interest=my_terms_of_interest)
    }
    dev.off()
}


#sort according to adjusted p-value and write results to file
myresultsOrdered <- myresults[order(myresults$padj),]
write.table(
        myresultsOrdered,
        sep="\t",
        file=paste(output_prefix,".txt",sep=""))

#print summary of deg analysis
sink(paste(output_prefix,"_summary.txt",sep=""))
summary(myresults)
sink()



###########################################################################################
###########################################################################################
# ANNOTATE WITH GENE SYMBOL IN ADDITION TO ENSEMBL ID

pdf(paste(output_prefix,"_geneCountPlots.pdf",sep=""))
for (i in 1:nrow(mydds))
{
    plotCounts(
            mydds,
            gene=i,
            xlab=my_design,
            intgroup=my_terms_of_interest,
            replaced=("replaceCounts" %in% names(assays(mydds))))
}
dev.off()



###########################################################################################
###########################################################################################
## PRINT SESSION INFO

sink(paste(output_prefix,"_sessionInfo.txt",sep=""))
print(sessionInfo())
sink()


###########################################################################################
###########################################################################################
## FIX ANNO

# to convert ENSEMBL gene ids to gene symbols later, I have to remove the decimal now
rownames(myresultsOrdered) <- unlist(strsplit(rownames(myresultsOrdered), split='\\.'))[2*(1:length(rownames(myresultsOrdered)))-1]

# rewrite deg analysis results to file
# -> overwrite if decimal removal was successful (will fail for Susi's circRNAs) and keep existing version otherwise
write.table(
        myresultsOrdered,
        sep="\t",
        file=paste(output_prefix,".txt",sep=""))

