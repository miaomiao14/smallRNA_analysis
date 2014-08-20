
library ("ggplot2")
library ("ggthemes")
library ("RColorBrewer")

draw_agg = function (t1, name) {
	plots = read.table(t1, F, sep="\t")
	colnames(plots) = c('Position','Signal')
	ggplot(plots, aes(x=Position,y=Signal)) +
	theme_few () +
	scale_colour_few() +
	theme( panel.border = element_blank () ,
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		plot.title=element_text(family="Helvetica", lineheight=.8) ) +
	ggtitle(name) +
	geom_line( size=1, alpha=0.75 ) +
	geom_vline(xintercept=26, slope=0,colour = "red") +
	# geom_vline(xintercept=31, slope=0,colour = "green") +
	geom_vline(xintercept=52, slope=0,colour = "red") +
	# geom_vline(xintercept=56, slope=0,colour = "green") +
	geom_vline(xintercept=78, slope=0,colour = "red") +
	xlab ("Position (bp)") +
	ylab("signal")
}
argv  = commandArgs (TRUE)
pdf (paste(argv[2], ".pdf", sep=''))
draw_agg (argv[1], argv[3])
g = dev.off ()
