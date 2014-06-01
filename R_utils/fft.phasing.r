library(plyr)
library(doBy)
#####Download from http://stackoverflow.com/questions/3485456/useful-little-functions-in-r

#main function
#example
x=c(1:1000)
y=sin(2*pi*x*0.01)+3*sin(2*pi*x*0.04)
pdf("fft.transform.example.pdf")
par(mfrow=c(2,1))
plot(x,y,ty='l')
test_ff=plotFFT(x,y,1)
dev.off()

#real data
#real data
plotFFTaxes <- function(x, y, samplingFreq, shadeNyq=FALSE, showPeriod = TRUE)
{
	Nyq.Freq <- samplingFreq/2
	FFTFreqs <- getFFTFreqs(Nyq.Freq, y)
	
	FFT <- fft(y)
	modFFT <- Mod(FFT)
	FFTdata <- cbind(FFTFreqs, modFFT)
	plot(FFTdata[1:nrow(FFTdata)/2,], t="l", pch=20, lwd=2, cex=0.8, axes=F, main="",
			xlab="Frequency (Hz)", ylab="Power",ylim=c(0,6))
	axis(2,tck=0.01, at=seq(0,6,by=2),las=2) ##modify
	axis(1,tck=0.01, at=seq(0,max(FFTFreqs),by=0.02),las=2) ##modify
	if (showPeriod == TRUE)
	{
		# Period axis on top        
		a <- axis(3, lty=0, labels=FALSE)
		axis(3, cex.axis=0.6, labels=format(1/a, digits=2), at=a)
	}
	if (shadeNyq == TRUE)
	{
		# Gray out lower frequencies
		rect(0, 0, 2/max(x), max(FFTdata[,2])*2, col="gray", density=30)
	}
	
	ret <- list("freq"=FFTFreqs, "FFT"=FFT, "modFFT"=modFFT)
	return (ret)
}
distance=x
distancePairFre=test
par(mfrow=c(2,1))
#plot(x,y,ty='l')
#test_fft=plotFFT(distance,distancePairFre,1)
plot(distance,distancePairFre,ty='l')
test_ff=plotFFTaxes(distance,distancePairFre,1)




# Gets the frequencies returned by the FFT function
getFFTFreqs <- function(Nyq.Freq, data)
{
	if ((length(data) %% 2) == 1) # Odd number of samples
	{
		FFTFreqs <- c(seq(0, Nyq.Freq, length.out=(length(data)+1)/2), 
				seq(-Nyq.Freq, 0, length.out=(length(data)-1)/2))
	}
	else # Even number
	{
		FFTFreqs <- c(seq(0, Nyq.Freq, length.out=length(data)/2), 
				seq(-Nyq.Freq, 0, length.out=length(data)/2))
	}
	
	return (FFTFreqs)
}

# FFT plot
# Params:
# x,y -> the data for which we want to plot the FFT 
# samplingFreq -> the sampling frequency
# shadeNyq -> if true the region in [0;Nyquist frequency] will be shaded
# showPeriod -> if true the period will be shown on the top
# Returns a list with:
# freq -> the frequencies
# FFT -> the FFT values
# modFFT -> the modulus of the FFT
plotFFT <- function(x, y, samplingFreq, shadeNyq=FALSE, showPeriod = TRUE)
{
	Nyq.Freq <- samplingFreq/2
	FFTFreqs <- getFFTFreqs(Nyq.Freq, y)
	
	FFT <- fft(y)
	modFFT <- Mod(FFT)
	FFTdata <- cbind(FFTFreqs, modFFT)
	plot(FFTdata[1:nrow(FFTdata)/2,], t="l", pch=20, lwd=2, cex=0.8,  main="",
			xlab="Frequency (Hz)", ylab="Power")
	
	if (showPeriod == TRUE)
	{
		# Period axis on top        
		a <- axis(3, lty=0, labels=FALSE)
		axis(3, cex.axis=0.6, labels=format(1/a, digits=2), at=a)
	}
	if (shadeNyq == TRUE)
	{
		# Gray out lower frequencies
		rect(0, 0, 2/max(x), max(FFTdata[,2])*2, col="gray", density=30)
	}
	
	ret <- list("freq"=FFTFreqs, "FFT"=FFT, "modFFT"=modFFT)
	return (ret)
}
