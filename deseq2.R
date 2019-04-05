args <- commandArgs(TRUE)
input_file <- args[1]
my_design <- args[2]
my_reference_level <- args[3] # factor level to use as the reference for calculated fold changes
my_alpha <- as.numeric(args[4])
my_lfc <- as.numeric(args[5])
output_prefix <- args[6]

# save script arguments to file
sink(paste(output_prefix, "_args.txt", sep=""), split=TRUE)
cat("Input file: ", input_file)
cat("\nDesign formula: ", my_design)
cat("\nReference level: ", my_reference_level)
cat("\nFDR cutoff alpha: ", my_alpha)
cat("\nMinimum log fold change: ", my_lfc, "\n")
sink()

# load required libraries
suppressMessages(library(DESeq2))
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
    myddsHTSeq <- myddsHTSeq[ rowSums(counts(myddsHTSeq)) > 2, ] # require min 3 total counts
    #myddsHTSeq <- myddsHTSeq[ rowSums(counts(myddsHTSeq) > 0) > 3, ] # require min 4 samples expressing feature at all

    # differential expression testing
    mydds <- DESeq(myddsHTSeq)
    return(mydds)
}

cat("\nRunning DESeq2...\n")
mydds <- fromFile(input_file)
my_terms_of_interest <- remove_interaction_terms(design_str_to_vector_of_str_terms(my_design))
save(mydds, file=paste(output_prefix,".RData", sep=""))
#load(paste(output_prefix,".RData", sep=""))

###########################################################################################
###########################################################################################
# PREPARE FOR EXPLORATORY ANALYSIS

# save normalization factors to file
# --> normalized counts = raw counts / normlization factors
write.table(
        sizeFactors(mydds),
        col.names=c("sizeFactor"),
        sep="\t",
        file=paste(output_prefix,"_sizeFactors.txt", sep=""))

# if "replace_sizeFactors.txt" is present in the current folder, load and overwrite current ones
#   --> format of "replace_sizeFactors.txt" must be the same as that saved above
sizeFactor_file <- file.path(dirname(input_file), "replace_sizeFactors.txt")
if (file.exists(sizeFactor_file)) {
    tmp_df <- read.table(sizeFactor_file)
    my_sizeFactors <- tmp_df$sizeFactor
    names(my_sizeFactors) <- rownames(tmp_df)

    sizeFactors(mydds) <- my_sizeFactors
}


# plot sparsity
# plot of the concentration of counts in a single sample over the sum of counts per gene.
# --> useful diagnostic for datasets which might not fit a negative binomial assumption:
#     genes with many zeros and individual very large counts are difficult to model with
#     the negative binomial distribution.
cat("\nPrinting sparsity plot...\n")
pdf(paste(output_prefix,"_sparsity.pdf",sep=""), height=10)
plotSparsity(mydds)
cat("done\n")
invisible(dev.off())


# get normalized read counts
counts <- counts(mydds,normalized=TRUE)
write.table(
        counts,
        sep="\t",
        quote=FALSE,
        file=paste(output_prefix,"_all_countsNormalized.txt", sep=""))

# get transformed read counts
trans <- rlog(mydds)
#save(trans, file=paste(output_prefix,"_trans.RData", sep=""))
#load(file=paste(output_prefix,"_trans.RData", sep=""))
write.table(
        assay(trans),
        sep="\t",
        quote=FALSE,
        file=paste(output_prefix,"_all_countsNormalizedTransformed.txt", sep=""))


# calculate coefficient of variation
standev <- apply(counts,1,sd)
avrg <- apply(counts,1,mean)
cv <- standev/avrg



# PREPARE PLOTTING
cat("\nPrinting exploratory PCA and clusters...\n")
pdf(paste(output_prefix,"_all_exploratory.pdf",sep=""), height=10)

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
    write.table(
            pca["rotation"],
            sep="\t",
            quote=FALSE,
            file=paste(output_prefix, "_pca.txt", sep=""))
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
my_plotPCA <- function(data, intgroup, xPC, yPC, xPC_label, yPC_label, ntop=500, pass_outputprefix=output_prefix)
{
    pcaData <- my_prepPCA(data, xPC, yPC, intgroup=intgroup, ntop=ntop, returnData=TRUE, output_prefix=pass_outputprefix)
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
    colors <- unique(ggplot_build(my_plotPCA(trans, col, 1, 2, "PC1: ", "PC2: ", pass_outputprefix=paste(output_prefix, "_all_exploratory", sep="")))$data[[1]][["colour"]])

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
    print(my_plotPCA(trans, col, 1, 2, "PC1: ", "PC2: ", pass_outputprefix=paste(output_prefix, "_all_exploratory", sep="")))
    print(my_plotPCA(trans, col, 3, 4, "PC3: ", "PC4: ", pass_outputprefix=paste(output_prefix, "_all_exploratory", sep="")))
    print(my_plotPCA(trans, col, 5, 6, "PC5: ", "PC6: ", pass_outputprefix=paste(output_prefix, "_all_exploratory", sep="")))
}


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

    # pearson correlation distance (taken from pheatmap source code)
    # --> pearson's r = cov(X,Y) / sd(X)*sd(Y)
    # ==> undefined if sd() is 0 for any sample! --> in that case, omit plot
    if (sum(apply(transCounts, 2, sd) == 0) > 0) {
        cat("Clustering based on Pearson correlation could not be performed\n--> standard deviation of at least one sample was 0\n")
    } else {
        sampleDists_corr <- as.dist(1 - cor(transCounts))
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
    }

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

cat("done\n")
invisible(dev.off())




###########################################################################################
###########################################################################################
# PROCESS RESULTS OF DGE ANALYSIS


myresults <- results(mydds, alpha=my_alpha, altHypothesis="greaterAbs", lfcThreshold=my_lfc)

# sort according to adjusted p-value and write results to file
write.table(
        myresults[order(myresults$padj),],
        sep="\t",
        quote=FALSE,
        file=paste(output_prefix,"_all.txt",sep=""))

# print summary of deg analysis
sink(paste(output_prefix,"_summary.txt",sep=""))
summary(myresults)
sink()



###########################################################################################
###########################################################################################
# PROCESS DEGs

get_candidate_genes <- function() {
    candidates_file <- file.path(dirname(input_file), "candidates.txt")
    if (file.exists(candidates_file)) {
        candidates <- read.table(candidates_file, header=FALSE)

        return(which(rownames(mydds) %in% candidates[,1]))
    } else {
        return(c())
    }
}


degs <- which(myresults$padj < my_alpha)
candidates <- get_candidate_genes()

for (gois in c("degs", "candidates")) {
    goi_type <- gois
    gois <- eval(parse(text=gois))

    # check if there are genes of interest (gois) to plot
    if(length(gois) >= 1)
    {
        goiCounts <- assay(trans)[gois, ]
        write.table(
                goiCounts,
                sep="\t",
                quote=FALSE,
                file=paste(output_prefix, "_", goi_type, "_countsNormalizedTransformed.txt", sep=""))

        if (length(gois) > 1) {
            # euclidean distances
            goi_sampleDists <- dist(t(goiCounts))
            goi_sampleDistMatrix <- as.matrix(goi_sampleDists)

            # plot PCAs, again color-coding each of the annotations separately
            cat(paste("\nPrinting ", goi_type, " PCA and clusters...\n", sep=""))
            pdf(paste(output_prefix, "_", goi_type,".pdf",sep=""), height=10)
            for (col in cols)
            {
                if (length(gois) >= 2) {
                    print(my_plotPCA(trans[gois, ], col, 1, 2, "PC1: ", "PC2: ", ntop=length(gois), pass_outputprefix=paste(output_prefix, "_", goi_type, sep="")))
                }
                if (length(gois) >= 4) {
                    print(my_plotPCA(trans[gois, ], col, 3, 4, "PC3: ", "PC4: ", ntop=length(gois), pass_outputprefix=paste(output_prefix, "_", goi_type, sep="")))
                }
                if (length(gois) >= 6) {
                    print(my_plotPCA(trans[gois, ], col, 5, 6, "PC5: ", "PC6: ", ntop=length(gois), pass_outputprefix=paste(output_prefix, "_", goi_type, sep="")))
                }
            }

            # plot the same heatmaps as above, now for GOIs only
            pheatmap(
                    goiCounts,
                    annotation_col=df,
                    annotation_color=anno_colors,
                    clustering_distance_cols=goi_sampleDists,
                    main=paste("Clustered by Euclidean distance of ", goi_type, sep=""),
                    scale="row",
                    show_rownames=length(gois)<=50,
                    treeheight_row=0,
                    fontsize=5,
                    border_color=NA)

            # pearson correlation distances (taken from pheatmap source code)
            # --> pearson's r = cov(X,Y) / sd(X)*sd(Y)
            # ==> undefined if sd() is 0 for any sample! --> in that case, omit plot
            if (sum(apply(goiCounts, 2, sd) == 0) > 0) {
                cat("Clustering based on Pearson correlation could not be performed\n--> standard deviation of at least one sample was 0\n")
            } else {
                goi_sampleDists_corr <- as.dist(1 - cor(goiCounts))
                pheatmap(
                        goiCounts,
                        annotation_col=df,
                        annotation_color=anno_colors,
                        clustering_distance_cols=goi_sampleDists_corr,
                        main=paste("Clustered by Pearson correlation of ", goi_type, sep=""),
                        scale="row",
                        show_rownames=length(gois)<=50,
                        treeheight_row=0,
                        fontsize=5,
                        border_color=NA)
            }

            pheatmap(
                    goi_sampleDistMatrix,
                    annotation_col=df,
                    annotation_color=anno_colors,
                    clustering_distance_rows=goi_sampleDists,
                    clustering_distance_cols=goi_sampleDists,
                    main=paste("Euclidean inter-sample distance based on ", goi_type, sep=""),
                    fontsize=5,
                    border_color=NA)
            invisible(dev.off())
        }

        # plot GOI counts
        cat(paste("Printing ", goi_type, " counts...\n", sep=""))
        pdf(paste(output_prefix, "_", goi_type, "_geneCountPlots.pdf",sep=""))
        for (i in gois)
        {
            if (length(my_terms_of_interest) > 1) {
                # simple plot
                my_plotCounts(
                        mydds,
                        gene=i,
                        design=paste("~", rev(unlist(strsplit(my_design, split=" ")))[1]),
                        terms_of_interest=rev(my_terms_of_interest)[1])
            }
            # pretty plot
            my_plotCounts(
                    mydds,
                    gene=i,
                    design=my_design,
                    terms_of_interest=my_terms_of_interest)
        }
        cat("done\n")
        invisible(dev.off())
    } else {
        cat(paste("\nNo ", goi_type, " to process!\n", sep=""))
    }
}



###########################################################################################
###########################################################################################
# PRINT ALL GENE COUNTS (simple only because pretty will take too long for all)

# calc max number of factor level combinations
max_number_of_factor_combinations <- function(mydds) {
    n_combinations <- 1
    for (this_factor in my_terms_of_interest) {
        n_levels <- length(levels(colData(mydds)[[this_factor]]))
        n_combinations <- n_combinations * n_levels
    }
    return(n_combinations)
}

plottable_terms <- my_terms_of_interest
plottable_design <- my_design
if (max_number_of_factor_combinations(mydds) > 4) {
    plottable_terms <- rev(my_terms_of_interest)[1]
    plottable_design <- paste("~", rev(unlist(strsplit(my_design, split=" ")))[1])
}


cat("\nPrinting all counts...\n")
pdf(paste(output_prefix,"_all_geneCountPlots.pdf",sep=""))
for (i in 1:nrow(mydds))
{
    plotCounts(
            mydds,
            gene=i,
            xlab=plottable_design,
            intgroup=plottable_terms,
            replaced=("replaceCounts" %in% names(assays(mydds))))
}
cat("done\n")
invisible(dev.off())



###########################################################################################
###########################################################################################
## PRINT SESSION INFO

sink(paste(output_prefix,"_sessionInfo.txt",sep=""))
print(sessionInfo())
sink()
