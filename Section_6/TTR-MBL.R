#!/usr/bin/env Rscript
library(PopGenome)
library(ape)

args <- commandArgs()

#Computing Transition-Transversion ratio of alignment
GENOME.class <- readData(args[6])
TTR <- get.sum.data(GENOME.class)
write.csv(TTR, file = "TTR.csv")

#Computing Mean Branch Length
mytree <- read.tree(args[7])
MBL <- mean(mytree$edge.length)
write.csv(MBL, file = "MBL.csv")