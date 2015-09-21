args <- commandArgs(TRUE)
#print(args)
INPUT <- args[1]


db <- read.table(INPUT);
bint <- c()

for (i in 1:(dim(db)[1]))
{
	cat("new row: ",i,"\n")
	row <- db[i,]
	#read in exon described in the current row
	if(row[4] == 1 )
	{	
		bint <- rep(TRUE,row[9])
	} else
	{
		bint <- c(bint,rep(!bint[length(bint)],row[9]))
	}
	#test if this was the transcript's last exon and if it is, count exon-exon junctions, otherwise continue
	if(row[10] == 0)
	{	
		#generate mask to sample reads from the transcript
		mask <- c(rep(TRUE,100),rep(FALSE,row[8]-100))
		read <- bint[mask]
		junctions <- sum(xor(read[1:(length(read)-1)],read[2:length(read)]))
		cat(junctions,"\n")

		#generate all possible reads across the transcript
		for(j in 1:(as.matrix(row[8]-100)))
		{	
			#shift mask
			mask <- c(FALSE,mask[1:(length(mask)-1)])				
			read <- bint[mask]
			junctions <- sum(xor(read[1:(length(read)-1)],read[2:length(read)]))
			cat(junctions,"\n")
		}
	}
}
