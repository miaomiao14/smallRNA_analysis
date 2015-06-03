

args = commandArgs (T)
file = args[1]
filename = basename(file)

plot (log2(te3[,1]), log2(te3[,2]), 
		xlim=c(-20,20), 
		ylim=c(-20,20), 
		xlab="", 
		ylab="", 
		pch=21, 
		col="white", 
		bg="red", 
		xaxt='n', 
		yaxt='n', 
		frame=F, 
		cex=2.0, 
		main="transposon abundance by RNASeq", 
		ex.lab=1.5, 
		cex.main=2
) + 
	abline (0,1, lty=2, lwd=2) + 
	axis (1, tck=0.01, lwd=3, at=seq(-20,20,5), cex.axis=1.5) + 
	axis (2, tck=0.01, lwd=3, at=seq(-20,20,5), cex.axis=1.5)

identify (log2(te3[,1]), log2(te3[,2]), labels=rownames (te3))
