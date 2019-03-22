#!/usr/bin/Rscript

####
# T.J.Blätte
# March 2017
####
#
# Draw a waterfall plot / mutation heatmap
#       of the data provided.
#
# Args:
#   input_file: Name of / Path to the tsv file
#       specifying matrix to visualize as heatmap.
#       Script was intended for mutation heatmaps
#       of matrices listing which genes were mutated
#       in which sample. Genes are listed in rows,
#       samples are listed in columns, and the
#       first row and column contain gene and
#       sample names respectively. Use '0' to
#       indicate a sample was not mutated and
#       '1' to indicate it was - do not quote
#       values. Non-integer values may be provided
#       as well, to plot for example VAFs of mutated genes.
#   heatmap_title: Title to print above the
#       plot. Leave empty to omit.
#
# Output:
#   *_heatmap.pdf: Heatmap of input_file
#
####

# example input file:
#
# sampleID    ASXL1   BCOR    EZH2    MLLptd  RUNX1   SF3B1   U2AF1
# 2805          0       0       0       0       1       0       0
# 2813          0       0       0       0       0       0       0
# 2814          0       0       0       0       0       0       0
# 2817          0       0       1       0       0       0       0
# 2820          0       0       0       0       0       0       0
# 2821          1       0       0       0       1       0       1
# 2830          0       0       0       1       0       0       0
# 2843          0       0       0       0       0       0       1


library(ggplot2)
library(reshape2)
library(ggthemes)

args <- commandArgs(TRUE)
input_file <- args[1]
input_file_prefix <- tools::file_path_sans_ext(input_file)
output_file <- paste(input_file_prefix, "_heatmap.pdf", sep="")
heatmap_title <- args[2]

# set heatmap title to empty string if none is given
if (is.na(heatmap_title)) {
    heatmap_title <- ""
}

# print input and output parameters
cat("Input file: ", input_file, "\n", sep="")
cat("Plot title: ", heatmap_title, "\n", sep="")
cat("Output file: ", output_file, "\n", sep="")


cat("\nPreparing heatmap...")
mutations <- as.matrix(
        read.table(
                input_file,
                header=TRUE,
                check.names=FALSE,
                row.names=1,
                na.string=0))

mutations_freq <- apply(mutations,2,sum,na.rm=T)

mutations_sorted <- mutations[,(order(mutations_freq))]
mutations_sorted <- transform(data.frame(mutations_sorted))
mutations_sorted <- mutations_sorted[ do.call(order,rev(mutations_sorted)), ]


mutations_df <- data.frame(
        mutations_sorted,
        patients=factor(rownames(mutations_sorted), levels=rownames(mutations_sorted)))

mutations_df_melted <- melt(mutations_df, id.vars="patients")


cat("\nPlotting...\n")
pdf(paste(
        input_file_prefix,
        "_heatmap.pdf",
        sep=""))

plot <- ggplot(
            data=mutations_df_melted,
            aes(patients, variable, fill=value)) +
            #aes(patients, variable, fill=value*as.numeric(variable))) +
            geom_tile(color="grey60", show.legend=TRUE) +
            scale_fill_gradient(na.value="white", low="#CCCCFF", high="#333399") +
            theme_tufte(base_family="Helvetica") +
            ylab(NULL) +
            xlab(NULL) +
            ggtitle(heatmap_title) +
            theme(axis.ticks=element_blank(),  legend.title = element_blank())

# do not print x-axis labels if there are so many they'd overlap
if (dim(mutations)[1] > 15) {
    plot <- plot +
            theme(axis.text.x=element_blank())
}

print(plot)
invisible(dev.off())
cat("done\n")
