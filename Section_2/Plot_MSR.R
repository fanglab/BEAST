#!/usr/bin/env Rscript

args <- commandArgs()
data1 <- data.frame(read.table(args[6], header = TRUE))
chrsize <- as.numeric(args[7])

library(ggplot2)
library(RColorBrewer)

pdf(file = "MSR_Plot.pdf", width = 11, height = 8)

ggplot(data1, aes(x=Position, y=MSR_Scale, z=PValue)) + geom_tile(binwidth = c(100,1), stat="summary_2d", na.rm = TRUE) +  
scale_fill_gradient2("SFC", limits=c(-0.5, 1), low = "darkblue", mid = "white", high = "red", na.value = "white") +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = "white"), panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position = "bottom") + scale_y_reverse() +
labs(x = "Chromosome position (bp)", y = "MSR Scale") + scale_x_continuous(limits=c(0, chrsize))

dev.off()
