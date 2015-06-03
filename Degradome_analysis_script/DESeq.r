library(reshape)
library(DESeq)
library(edgeR)

args <- commandArgs (TRUE);
file <- args[1];
filename=basename(file)
seqDepFile <- args[2]
dir <- args[3];

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

cds<-newCountDataSet (countsTable,genotypes)
cds<-estimateSizeFactors(cds)

#cds<-estimateDispersions (cds) #it does not work here as we don't have replicates
#res<- nbinomTest(cds,genotyps[1],genotypes[2])

#uq1=quantile(htv1$count)[4];
#uq2=quantile(htv2$count)[4];	
#df <-data.frame(quantile(htv1$count),quantile(htv2$count))
#write.table(df,paste(dir,"nf.txt",sep="/"),col.names=TRUE)


###USE edgeR instead,use the sizeFactors from DESeq

cdsedgeR<-DGEList (counts=countsTable,group=genotypes,lib.size=1e7*sizeFactors(cds))
cdsedgeR$genes=rownames(countsTable)
cdsedgeR<-estimateCommonDisp(cdsedgeR) #set to 0
cdsedgeR<-estimateTagwiseDisp(cdsedgeR)
RescdsedgeR <- exactTest (cdsedgeR)
RescdsedgeRPadj <- p.adjust(RescdsedgeR$table$p.value,metho="BH")
statsTable<-topTags(RescdsedgeR,n=numofRows)

outFileName=paste(dir,"/",filename,".DEseqNF.edgeR.stat",sep="")
write.table(statsTable,outFileName,col.names=TRUE,sep="\t")

deStatsSummary <- summary( de <- decideTestsDGE(RescdsedgeR,p=0.05,adjust="BH"))
outFileName=paste(dir,"/",filename,".DEseqNF.edgeR.stat.summary",sep="")
write.table(deStatsSummary,outFileName,col.names=TRUE,sep="\t")

detags <- rownames(cdsedgeR)[as.logical(de)]
pdfname <- paste(dir,"/",filename,".DEseqNF.edgeR.stat.MAplot.pdf",sep="")
pdf(pdfname)
plotSmear(RescdsedgeR,de.tags=detags)
abline(h=c(-1,1),col="blue",lty=2)
dev.off()

###use edgeR's own NF
cdsedgeR<-DGEList (counts=countsTable,group=genotypes)
cdsedgeR<-calcNormFactors(cdsedgeR)
cdsedgeR<-estimateCommonDisp(cdsedgeR) #set to 0
cdsedgeR<-estimateTagwiseDisp(cdsedgeR)
RescdsedgeR <- exactTest (cdsedgeR)
RescdsedgeRPadj <- p.adjust(RescdsedgeR$table$p.value,metho="BH")
statsTable<-topTags(RescdsedgeR,n=numofRows)

outFileName=paste(dir,"/",filename,".EdgeRNF.edgeR.stat",sep="")
write.table(statsTable,outFileName,col.names=TRUE,sep="\t")

deStatsSummary <- summary( de <- decideTestsDGE(RescdsedgeR,p=0.05,adjust="BH"))
outFileName=paste(dir,"/",filename,".EdgeRNF.edgeR.stat.summary",sep="")
write.table(deStatsSummary,outFileName,col.names=TRUE,sep="\t")

detags <- rownames(cdsedgeR)[as.logical(de)]
pdfname <- paste(dir,"/",filename,".EdgeRNF.edgeR.stat.MAplot.pdf",sep="")
pdf(pdfname)
plotSmear(RescdsedgeR,de.tags=detags)
abline(h=c(-1,1),col="blue",lty=2)
dev.off()




###use seqDep as sizeFactor
seqDep<-read.delim(seqDepFile,sep =",",F)
libSizes<-c( assign(as.character(seqDep$V1[1]),seqDep$V2[1]),assign(as.character(seqDep$V1[2]),seqDep$V2[2]) )
names(libSizes)=seqDep$V1
cdsedgeR<-DGEList (counts=countsTable,group=genotypes,lib.size=libSizes)
cdsedgeR<-estimateCommonDisp(cdsedgeR) #set to 0
cdsedgeR<-estimateTagwiseDisp(cdsedgeR)
RescdsedgeR <- exactTest (cdsedgeR)
RescdsedgeRPadj <- p.adjust(RescdsedgeR$table$p.value,metho="BH")
statsTable<-topTags(RescdsedgeR,n=numofRows)
outFileName=paste(dir,"/",filename,".seqDepNF.edgeR.stat",sep="")
write.table(statsTable,outFileName,col.names=TRUE,sep="\t")

deStatsSummary <- summary( de <- decideTestsDGE(RescdsedgeR,p=0.05,adjust="BH"))
outFileName=paste(dir,"/",filename,".seqDepNF.edgeR.stat.summary",sep="")
write.table(deStatsSummary,outFileName,col.names=TRUE,sep="\t")

detags <- rownames(cdsedgeR)[as.logical(de)]
pdfname <- paste(dir,"/",filename,".seqDepNF.edgeR.stat.MAplot.pdf",sep="")
pdf(pdfname)
plotSmear(RescdsedgeR,de.tags=detags)
abline(h=c(-1,1),col="blue",lty=2)
dev.off()



