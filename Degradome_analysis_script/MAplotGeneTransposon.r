#!/usr/bin/env Rscript
library(reshape)
library(DESeq)
library(edgeR)


plot.MA <- function (input,nftype,trnames,type)
{
	NFcountsTable=input
	#add M and A columns by transform
	normalizedCountsTable=transform(NFcountsTable,M=log2(NFcountsTable[,1]/NFcountsTable[,2]),A=log2(NFcountsTable[,1]*NFcountsTable[,2])/2)
	normalizedCountsTable=transform(normalizedCountsTable,ann=ifelse(rownames(normalizedCountsTable) %in% trnames$transposon,"transposon","gene"),feature=rownames(normalizedCountsTable),dt=rep(type,nrow(normalizedCountsTable)))
	#plot gene and transposons
	MA <- new("MAList")
	MA$M = normalizedCountsTable$M
	MA$A = normalizedCountsTable$A

	numofRows=nrow(MA)
	status=normalizedCountsTable$ann
	pdfname <- paste(dir,"/",filename,".",nftype,".MAplot.genetransposon.pdf",sep="")
	pdf(pdfname)
	m=paste("MA-Plot ",colnames(NFcountsTable)[1],"/",colnames(NFcountsTable)[2],sep="")
	plotMA(MA,main=m,status=status,values=c("gene","transposon"),col=c("lightgrey","red"),cex=rep(1.5,numofRows),legend=TRUE,ylim=c(-max(abs(MA$M)),max(abs(MA$M))),xlim=c(0,max(MA$A)))
	abline(h=c(-1,1),col="darkgrey",lty=2)
	dev.off()
	
	#keep the intermediate results
	outfilename<- paste(dir,"/",filename,".",nftype,".genetransposon.normalizedcounts.txt",sep="")
	write.table(normalizedCountsTable, file = outfilename, append = FALSE, quote = FALSE, sep = "\t",
			eol = "\n", na = "NA", dec = ".", row.names = FALSE,
			col.names = TRUE, qmethod = c("escape", "double"),
			fileEncoding = "")
	
	
	#plot transposon only
	normalizedCountsTableTRN=as.data.frame(lapply(subset(normalizedCountsTable,as.character(ann)=="transposon"),'[',drop=TRUE))
	
	TRNMA2<- new("MAList")
	TRNMA2$M = normalizedCountsTableTRN$M
	TRNMA2$A = normalizedCountsTableTRN$A

	status1 = ifelse(TRNMA2$M>1,1,0)
	status2 = ifelse(TRNMA2$M<(-1),-1,0)
	status = status1 + status2
	numofRows=nrow(TRNMA2)
	pdfname <- paste(dir,"/",filename,".",nftype,".MAplot.transposon.pdf",sep="")
	pdf(pdfname)
	m=paste("MA-Plot ",colnames(NFcountsTable)[1],"/",colnames(NFcountsTable)[2],sep="")
	plotMA(TRNMA2,main=m,status=status,values=c("1","0","-1"),col=c("red","black","green"),cex=rep(1.5,numofRows),legend=FALSE,ylim=c(-max(abs(TRNMA2$M)),max(abs(TRNMA2$M))),xlim=c(0,max(TRNMA2$A)))
	abline(h=c(-1,1),col="darkgrey",lty=2)
	dev.off()
	
	#keep the intermediate results
	outfilename<- paste(dir,"/",filename,".",nftype,".transposon.normalizedcounts.txt",sep="")
	write.table(normalizedCountsTableTRN, file = outfilename, append = FALSE, quote = FALSE, sep = "\t",
			eol = "\n", na = "NA", dec = ".", row.names = FALSE,
			col.names = TRUE, qmethod = c("escape", "double"),
			fileEncoding = "")
	

}

args <- commandArgs (TRUE);
file <- args[1];
filename=basename(file)
seqDepFile <- args[2]
dir <- args[3];
type <- args[4]; #RSQ or DEG


group=read.table("/home/wangw1/pipeline/common/Zamore.group",F)
colnames(group)=c("transposon","g")

trnames=read.table("/home/wangw1/pipeline_dm3/common/trn.list",F)
colnames(trnames)=c("transposon")

htv <- read.table(file,F);
colnames(htv) <- c("gt","feature","count");
htv$feature=sub("FBgn0000004_17","FBgn0000004_17.6",htv$feature)
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
###use the sizeFactors from DESeq

cds<-newCountDataSet (countsTable,genotypes)
cds<-estimateSizeFactors(cds)
lib.size=sizeFactors(cds)
normalizedCountsTable1=sweep( countsTable, 2, lib.size, '/' )
plot.MA(normalizedCountsTable1,"DESeqNF",trnames,type)

###use edgeR's own NF
cdsedgeR<-DGEList (counts=countsTable,group=genotypes)
cdsedgeR<-calcNormFactors(cdsedgeR)
normalizedCountsTable2=sweep( countsTable, 2, cdsedgeR$samples$norm.factors, '/' )
#add M and A columns by transform
plot.MA(normalizedCountsTable2,"EdgeRNF",trnames,type)


###use seqDep as sizeFactor
#seqDep<-read.delim(seqDepFile,sep =",",F)
#libSizes<-c( assign(as.character(seqDep$V1[1]),seqDep$V2[1]),assign(as.character(seqDep$V1[2]),seqDep$V2[2]) )/1000000
#names(libSizes)=seqDep$V1
#normalizedCountsTable3=sweep( countsTable, 2, libSizes, '/' )
#plot.MA(normalizedCountsTable3,"seqDepNF",trnames,type)


#use upper quantile as sizeFactor

#uq1=quantile(countsTable[,1])[4];
#uq2=quantile(countsTable[,2])[4];	
#uq = c(uq1,uq2)
#write.table(df,paste(dir,"nf.txt",sep="/"),col.names=TRUE)


