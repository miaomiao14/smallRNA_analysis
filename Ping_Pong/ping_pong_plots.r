library(ggplot2)
library(calibrate)
library(reshape)
library(plyr)

groupingfun <- function (x) {if (x == "AU" || x == "CG" || x == "GC" || x=="UA") {y=1} else { y=0}}

plot_ua_va <- function (input,gt,outdir) {

	pdfname=paste(outdir,"/",gt,"_UA_VA_pair_counts",".pdf",sep="")
	pdf(pdfname,height=12,width=15)
	layout(matrix(1:4,2,2,byrow=TRUE))
	
	pairg=data.frame(pair=c("AU","UA","GC","CG","BU","VA","HC","DG"),group=c(1,1,1,1,0,0,0,0))
	
	pp=read.table(input,F)
	colnames(pp)=c("genotype","pair","NofPairs","raw");
	for (j in (1:length(levels(pp$genotype))))
	{
		#df=as.data.frame(lapply(subset(L,as.character(rank)==r),'[',drop=TRUE))
		piwipair=levels(pp$genotype)[j]
		f=as.data.frame(lapply( subset(pp,as.character(genotype)==piwipair),'[',drop=TRUE))
		m=paste(piwipair,"PP pair",sep=" ")
		f_order=f[order(f$pair),]
		f=f_order
		
		ff=merge(f,pairg,by="pair")
		ff$group=as.factor(ff$group)
		f1=as.data.frame(lapply( subset(ff,as.character(group)==1),'[',drop=TRUE))
		f2=as.data.frame(lapply( subset(ff,as.character(group)==0),'[',drop=TRUE))
		#f1=subset(ff,group==1)
		#f2=subset(ff,group==0)
		#############################################
		#b<-barplot(f1$NofPairs,names.arg=f1$pair,main=m,axes=F,ylim=c(0,1.2*max(f1$NofPairs)), col=rep("darkblue",4),border="white",cex.names=1,ylab="Paired Number of pairs")
		#text(x=b,y=f1$NofPairs+2.5,label=f1$raw,cex=1)
		#grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")
		#axis(2,ylim=c(0,1.2*max(f1$NofPairs)))
		#par(new=TRUE)
		
		#b<-barplot(f2$NofPairs,names.arg=f2$pair,,axes=F, col=rep("darkgrey",4),border="white",cex.names=1,ylab="Unpaired Number of pairs")
		#text(x=b,y=f2$NofPairs+2.5,label=f2$raw,cex=1)
		#grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")
		#axis(4,ylim=c(0,1.2*max(f2$NofPairs)))
		#############################################
			
		# data a
		a=seq(1,8,2)
	# data b, different scale
		b=seq(1,16,4)
	# width of your bar
		bar_width = .5
	# distance
		bar_dis = 0.5
	# making new plot
		plot.new ()
	# making window for the scale of data a
		plot.window (xlim=c(0,1.2*max(a)), ylim=c(0, 1.2*max(f1$NofPairs)),main=m)
	# draw "bar", which is actually rectangle
		rect (xleft = a-bar_width/2, xright=a+bar_width/2, ybottom=0, ytop=f1$NofPairs, ylim=c(0, 1.2*max(f1$NofPairs)), col="darkblue", border=F)
	# draw X axis
		axis (1,tick=0.01,las=1)
	# draw left Y axis
		axis (2, ylim = c(0, 1.2*max(f1$NofPairs)),  at = seq (0, 1.2*max(f1$NofPairs), as.integer ( 1.2*max(f1$NofPairs)/5 ) ) )
	# draw data b for its scale
		plot.window (xlim=c(0,1.2*max(a)), ylim=c(0, 1.2*max(f2$NofPairs)))
	# draw data b, distance between two differents data sesÊ
		rect (xleft=a+bar_dis-bar_width/2, xright=a+bar_dis+bar_width/2, ybottom=0, ytop=f2$NofPairs, ylim = c(0, 1.2*max(f2$NofPairs)), col="darkgrey", border=F )
		axis (4, ylim = c(0, 1.2*max(f2$NofPairs)), at = seq (0, 1.2*max(f2$NofPairs), as.integer (1.2*max(f2$NofPairs)/5)))
	
	}
	dev.off()	
}