library (ggplot2)
library (reshape)
argv = commandArgs (TRUE)

filename=argv[1]
outdir=argv[2]

g = read.table ("/home/wangw1/pipeline/common/Zamore.group", F)
plot_PolarHisto_senseFraction=function(filename,outdir)
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
	file=basename(filename)
	pdfname=paste(outdir/file,"senseFraction.polarHisto.pdf",sep="")
	ggsave(pdfname,width=12,height=12)
	dev.off()
}

#plot_PolarHisto_senseFraction(input);





