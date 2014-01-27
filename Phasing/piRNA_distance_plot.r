library(gridExtra)
library(ggplot2)


argv = commandArgs (TRUE)
input=argv[1]


plot_distribution <- function (input) {
	file=input
	pp=read.table(file,F)	
	colnames(pp)=c("strand","chr","distance","reads")
	for ( i in 1:length(levels(pp$strand)) )
	{
		ppp=as.data.frame(lapply(subset(pp,strand==levels(pp$strand)[i]),'[',drop=TRUE))
		pdfname= paste(file,"_",levels(pp$strand)[i],"_", 'piRNA_55_distance_distribution.pdf', sep='')
		pdf(pdfname,width=25,height=20)
		grid.newpage()
		pushViewport(viewport(layout=grid.layout(2,2)))
		a=c(1,1,2,2)
		b=c(1,2,1,2)
		n=1
		for (k in 1:length(levels(ppp$chr)))
		{
			pppp=as.data.frame(lapply(subset(ppp,chr==levels(ppp$chr)[k]),'[',drop=TRUE))	
			ppppp=pppp[order(pppp$distance),]
			t=paste(levels(pp$strand)[i],'_',levels(ppp$chr)[k],sep="")
			p<-ggplot(ppppp,aes(distance,reads))+labs(title=t)
			p1<-p+geom_line(size=2,colour="red")+scale_x_continuous(limits = c(0, 100),breaks=seq(0,100,2),labels=seq(0,100,2))+geom_vline(xintercept = c(27,54),col="black",linetype=2)
			print(p1,vp=viewport(layout.pos.row = a[n], layout.pos.col = b[n]))
			n=n+1

		}
		dev.off()
	}
}
	
	

plot_distribution_summary <- function (input) {
		file=input
		filename=basename(input)
		pp=read.table(file,F)
		pdfname= paste(file, '_piRNA_distance_distribution.pdf', sep='')
		pdf(pdfname,width=9,height=5)
		colnames(pp)=c("distance","reads")
		pp=pp[order(pp$distance),]
		p<-ggplot(pp,aes(distance,reads))+labs(title=filename)
		p1<-p+geom_line(size=1,colour="red")+xlab("5'-5'end distance on the same strand")+scale_x_continuous(limits = c(0, 100),breaks=seq(0,100,4),labels=seq(0,100,4))+geom_vline(xintercept = c(27,54),col="black",linetype=2)
		print(p1)
		
		dev.off()
	}

	
fineplot_distribution_summary <- function (input) {
		file=input
		filename=basename(input)
		
		pp=read.table(file,F)
		colnames(pp)=c("distance","reads")
		pp=pp[order(pp$distance),]
		ppf=transform(pp,fr=reads/colSums(pp$reads))
		
		pdfname= paste(file, '_piRNA_distance_distribution.fraction.pdf', sep='')
		pdf(pdfname,width=9,height=5)
		p<-ggplot(ppf,aes(distance,fr))+labs(title=filename)
		p1<-p+geom_line(size=1,colour="red")+xlab("5'-5'end distance on the same strand")+scale_x_continuous(limits = c(0, 100),breaks=seq(0,100,4),labels=seq(0,100,4))+geom_vline(xintercept = c(27,54),col="black",linetype=2)
		print(p1)
		
		dev.off()
	}
	
	summary_masterTable <- function (input,outdir) {
		file=input
		filename=basename(input)
		pp=read.table(file,T)
		colnames(pp)[1]<-c("distance")
		pp=pp[order(pp$distance),]
		pp=pp[,-1]
		pp=as.matrix(pp)
		ppProb=prop.table(pp, margin=2)*100
		
		ppProb=transform(ppProb,distance=rownames(ppProb))
		ppProb$distance=as.numeric(as.character(ppProb$distance))
		outfilename=paste(outdir,"/","allpiRNAs.allgt.5-5.distance.min.distribution.summary.fraction.mastertable.txt",sep="")
		write.table(ppProb,outfilename,row.names = FALSE,sep = "\t",quote=FALSE)
}