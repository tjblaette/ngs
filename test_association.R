args <- commandArgs(TRUE)                                                                  


annotation_filename <- args[1] #"/media/data/tabl/pollack_ckaml/expression/p53Aberration_vs_wt/pollack_classLabels_p53mut_all_wBatch.tsv"
anno <- read.table(
        annotation_filename,
        header=1,
        check.names=FALSE,
        row.names=1)

attr1 <- args[2] #"CDKN1C"
attr2 <- args[3] #"TP53del"



contingency_table <- ""
if (attr1 %in% colnames(anno) && attr2 %in% colnames(anno)) {
    contingency_table <- table(anno[c(attr1, attr2)]) #matrix(
} else {
    normalized_counts_filename <- args[4] #"/media/data/tabl/pollack_ckaml/expression/p53Aberration_vs_wt/pollack_classLabels_p53mut_all_wBatch.tsv_DESeq2results_countsNormalized.txt"
    counts <- read.table(
            normalized_counts_filename,
            header=1,
            check.names=FALSE,
            row.names=1)

    gene_symbol <- counts$geneSymbol
    counts <- counts[,-1]
    counts <- counts[,rownames(anno)]

    if ((attr1 %in% colnames(anno) && attr2 %in% gene_symbol) || (attr1 %in% gene_symbol && attr2 %in% colnames(anno))) {
        if (attr1 %in% colnames(anno) && attr2 %in% gene_symbol) {
            attribute <- attr1
            gene <- attr2
        } else {
            if (attr2 %in% colnames(anno) && attr1 %in% gene_symbol) {
            attribute <- attr2
            gene <- attr1
            }
        }

        gene_median <- median(t(counts[which(gene_symbol == gene),]))

        samples_high <- colnames(counts)[counts[which(gene_symbol == gene),] > gene_median]
        samples_low <- colnames(counts)[counts[which(gene_symbol == gene),] <= gene_median]

        attribute_in_samples_high <- table(anno[attribute][samples_high,])
        attribute_in_samples_low <- table(anno[attribute][samples_low,])

        #print(attribute_in_samples_low)
        #print(attribute_in_samples_high)


        contingency_table <- matrix(
                c(attribute_in_samples_high, attribute_in_samples_low),
                nrow=2,
                ncol=2,
                dimnames=list(
                        attribute = names(attribute_in_samples_high),
                        expression = c(paste(gene,"high", sep="_"), paste(gene,"low", sep="_"))))
    } else {
        cat("\nAttribute not found in input files!\n\n")
        quit()
    }
}


if (ncol(contingency_table) >= 2 && nrow(contingency_table) >= 2) {
    #fisher.test(anno[[attribute]], as.factor(counts[which(gene_symbol == gene),] > gene_median))
    fisher <- fisher.test(contingency_table)


    print(fisher)
    cat("Inverse of Odds Ratio: ", 1 / fisher$estimate,"\n\n\n")

    cat("Contingency table:\n\n")
    print(contingency_table)

    cat("\n(The Odds Ratio is the ratio of the odds of a sample in", dimnames(contingency_table)[[1]], "group \"", rownames(contingency_table)[1], "\" being of group \"", colnames(contingency_table)[1], "\" vs \"", colnames(contingency_table)[2],"\")", sep="")
    cat("\n(The inverse Odds Ratio is the ratio of the odds of a sample in group \"", rownames(contingency_table)[1], "\" being of group \"", colnames(contingency_table)[2], "\" vs \"", colnames(contingency_table)[1],"\")\n", sep="")

    if (exists("gene_median")) {
        cat("\nThe median expression used to define high vs low expression was: ", gene_median, " normalized read counts\n\n", sep="")
    }
} else {
    cat("\nCannot perform Fisher's exact test because the samples all share the same form of one of the attributes!\n")
    cat("p-value = NA\n")
    cat("Contingency table:\n\n")
    print(contingency_table)
}
