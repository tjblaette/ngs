args <- commandArgs(TRUE)
print(args)
IN <- args[1] #*_DESeq2results_annotated_woutNA_wEntrezIDs.txt
ORG <- args[2] # organism abbreviation from KEGG: hsa for human, mmu for mouse
OUT <- args[3] # output file prefix


d <- read.table(IN)

# vector containing Entrez gene IDs of all genes in the experiment
all <- unique(d[,9])

# vector with logFCs of DE genes with names as Entrez gene IDs
de <- d[which(d[,8] < 0.05),4]
names(de) <- d[which(d[,8] < 0.05),9]

de_unique <- de[!duplicated(names(de))]


library("SPIA")

# run SPIA analysis and save result as tsv table
res <- spia(de=de_unique,all=all,organism=ORG,data.dir=paste("/NGS/KEGG/",ORG,"/",sep=""),nB=2000,plots=TRUE,beta=NULL,combine="fisher",verbose=TRUE)
write.table(res, file=paste(OUT,"_spia.txt", sep=""), sep="\t", row.names=FALSE)

res_sig <- res[res$pGFdr <= 0.05,]
write.table(res_sig, file=paste(OUT,"_spia_sig.txt", sep=""), sep="\t", row.names=FALSE)
 

