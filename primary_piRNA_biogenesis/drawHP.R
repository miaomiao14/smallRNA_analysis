library ("RColorBrewer")
argv  = commandArgs (TRUE)
data =read.table(argv[1], F)
pdf (paste(argv[1], ".pdf", sep=''))
heatmap(as.matrix(data), scale="row", Rowv=NA, Colv=NA,col=brewer.pal(9,"Blues"))
g = dev.off ()

