
argv = commandArgs (TRUE)
input=argv[1]

g = read.table ("/home/wangw1/pipeline/common/Zamore.group", F)
plotPolarHisto=function(filename)
{
	list=read.delim(filename,T)
	group=g
	colnames(group)=c("transposon","g")
	cd=merge(list,group,by="transposon")
	cd= cd[cd[,10]!=0,]
	cd_p=subset(cd,select=c("g","transposon","piRNA_antisense","piRNA_sense"))
	colnames(cd_p)=c("family","item","piRNA_antisense","piRNA_sense")
	cd_p<-melt(cd_p,c("family","item"),variable_name="score")
	
	p<-polarHistogram(cd_p,familyLabel=TRUE)
	print(p)
	pdfname=paste(filename,".pdf",sep="")
	ggsave(pdfname,width=12,height=12)
	dev.off()
}

	plotPolarHisto(input);





