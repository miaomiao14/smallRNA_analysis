library(gridExtra)
files=list.files(pattern="*.distribution$")

for (j in (1:length(files)))
{
  pp=read.table(files[j],F)
	colnames(pp)=c("strand","chr","distance","reads")
	for ( i in 1:length(levels(pp$strand)) )
	{
		ppp=as.data.frame(lapply(subset(pp,strand==levels(pp$strand)[i]),'[',drop=TRUE))
		pdfname= paste(files[j],"_",levels(pp$strand)[i],"_", 'piRNA_distance_distribution.pdf', sep='')
		#pdfname= paste(files[j], 'wt_unox_piRNA_distance_distribution.pdf', sep='')
		
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
			p<-ggplot(ppppp,aes(distance,reads))+opts(title=t)
			#p1<-p+geom_line(size=2,colour="red")+scale_x_continuous(limits = c(-50, 50),breaks=seq(-50,50,2),labels=seq(-50,50,2))+geom_vline(xintercept = c(-22,1),col="black",linetype=2)
			p1<-p+geom_line(size=2,colour="red")+scale_x_continuous(limits = c(0, 100),breaks=seq(0,100,2),labels=seq(0,100,2))+geom_vline(xintercept = c(27,54),col="black",linetype=2)
			print(p1,vp=viewport(layout.pos.row = a[n], layout.pos.col = b[n]))
			n=n+1
			#dev.off()
		}
		dev.off()
	}
}
