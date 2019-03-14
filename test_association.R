args <- commandArgs(TRUE)                                                                  
input_file <- args[1]


normalized_counts_filename <- args[1] #"/media/data/tabl/pollack_ckaml/expression/p53Aberration_vs_wt/pollack_classLabels_p53mut_all_wBatch.tsv_DESeq2results_countsNormalized.txt"
annotation_filename <- args[2] #"/media/data/tabl/pollack_ckaml/expression/p53Aberration_vs_wt/pollack_classLabels_p53mut_all_wBatch.tsv"

gene <- args[3] #"CDKN1C"
attribute <- args[4] #"TP53del"                                                                                                             


counts <- read.table(
        normalized_counts_filename, 
        header=1, 
        check.names=FALSE,
        row.names=1)

gene_symbol <- counts$geneSymbol
counts <- counts[,-1]

anno <- read.table(
        annotation_filename, 
        header=1,
        check.names=FALSE,
        row.names=1)

counts <- counts[,rownames(anno)]

gene_median <- median(t(counts[which(gene_symbol == gene),]))

samples_high <- colnames(counts)[counts[which(gene_symbol == gene),] > gene_median]
samples_low <- colnames(counts)[counts[which(gene_symbol == gene),] <= gene_median]

attribute_in_samples_high <- table(anno[attribute][samples_high,])
attribute_in_samples_low <- table(anno[attribute][samples_low,])



contingency_table <- matrix(
        c(attribute_in_samples_high, attribute_in_samples_low),
        nrow=2, 
        ncol=2, 
        dimnames=list(
                attribute = names(attribute_in_samples_high),
                expression = c(paste(gene,"high", sep="_"), paste(gene,"low", sep="_"))))


#fisher.test(anno[[attribute]], as.factor(counts[which(gene_symbol == gene),] > gene_median))
fisher <- fisher.test(contingency_table)


print(fisher)
cat("Inverse of Odds Ratio: ", 1 / fisher$estimate,"\n\n\n")

cat("Contingency table:\n\n")
print(contingency_table)

cat("\n(The Odds Ratio is the ratio of the odds of a sample in group \"", rownames(contingency_table)[1], "\" being of group \"", colnames(contingency_table)[1], "\" vs \"", colnames(contingency_table)[2],"\")", sep="")
cat("\n(The inverse Odds Ratio is the ratio of the odds of a sample in group \"", rownames(contingency_table)[1], "\" being of group \"", colnames(contingency_table)[2], "\" vs \"", colnames(contingency_table)[1],"\")", sep="")
cat("\nThe median expression used to define high vs low expression was: ", gene_median, " normalized read counts\n\n", sep="")

