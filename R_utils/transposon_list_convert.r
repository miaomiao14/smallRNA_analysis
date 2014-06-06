#!/usr/bin/env Rscript
library(reshape)
library(DESeq)
library(edgeR)


plot.MA <- function (input,nftype,type)
{
	NFcountsTable=input
	#add M and A columns by transform
	normalizedCountsTable=transform(NFcountsTable,M=log2(NFcountsTable[,1]/NFcountsTable[,2]),A=log2(NFcountsTable[,1]*NFcountsTable[,2])/2,feature=rownames(NFcountsTable),dt=rep(type,nrow(NFcountsTable)))

	#keep the intermediate results
	outfilename<- paste(dir,"/",filename,".",nftype,".normalizedcounts.txt",sep="")
	write.table(normalizedCountsTable, file = outfilename, append = FALSE, quote = FALSE, sep = "\t",
			eol = "\n", na = "NA", dec = ".", row.names = FALSE,
			col.names = TRUE, qmethod = c("escape", "double"),
			fileEncoding = "")
	
	MA <- new("MAList")
	MA$M = normalizedCountsTable$M
	MA$A = normalizedCountsTable$A
	status1 = ifelse(MA$M>1,1,0)
	status2 = ifelse(MA$M<(-1),-1,0)
	status = status1 + status2
	pdfname <- paste(dir,"/",filename,".",nftype,".MAplot.pdf",sep="")
	pdf(pdfname)
	m=paste("MA-Plot ",colnames(NFcountsTable)[1],"/",colnames(NFcountsTable)[2],sep="")
	plotMA(MA,main=m,status=status,values=c("1","0","-1"),col=c("red","black","green"),cex=rep(1.5,numofRows),legend=FALSE,ylim=c(-max(abs(MA$M)),max(abs(MA$M))),xlim=c(0,max(MA$A)))
	abline(h=c(-1,1),col="darkgrey",lty=2)
	dev.off()
	

	
}

args <- commandArgs (TRUE);
file <- args[1];
filename=basename(file)
dir <- args[2];
type <- args[3]; #RSQ or DEG

group=read.table("/home/wangw1/pipeline/common/Zamore.group",F)
colnames(group)=c("transposon","g")

htv <- read.table(file,F);
colnames(htv) <- c("gt","feature","count");
htvReshape=cast(htv,feature~gt) #cast from reshape
countsTable <- htvReshape[-1] #remove the feature column
genotypes=factor(colnames(htvReshape)[-1]) #create the conditions, type: factor
featureName=unlist(htvReshape[1]) #get the row names from the first column
rownames(countsTable)=featureName 
countsTable=as.data.frame(countsTable)
#remove 0,0 
##Go through each row and determine if a value is zero
row_sub = apply(countsTable, 1, function(row) all(row !=0 ))
##Subset as usual
countsTable=countsTable[row_sub,]

numofRows=nrow(countsTable)
plot.MA(countsTable,"nncNF",type)


###use the sizeFactors from DESeq

## cds<-newCountDataSet (countsTable,genotypes)
## cds<-estimateSizeFactors(cds)
## lib.size=sizeFactors(cds)
## normalizedCountsTable1=sweep( countsTable, 2, lib.size, '/' )
## plot.MA(normalizedCountsTable1,"DESeqNF",type)


