#!/usr/bin/env Rscript

library("Biostrings")
library("knitr")

args <- commandArgs()

## Importing TFBSs in fasta format ##
fastaFile <- readDNAStringSet(args[6])
seq_name = names(fastaFile)
sequences = paste(fastaFile)
df <- print(as.data.frame(sequences))

## Building an alignment matrix with one row per site and one column per aligned position ##
alignment <- t(as.data.frame(strsplit(sequences, "")))
colnames(alignment) <- 1:ncol(alignment)
row.names(alignment) <- paste(sep="_", "site", 1:nrow(alignment))
kable(alignment, align="c")

## Computing a count matrix from the alignment table ##
cts <- data.frame()
residues <- c("A", "C", "G", "T")
for (r in residues) {
  cts <- rbind(cts, apply(alignment==r, 2, sum))
}
colnames(cts) <- 1:ncol(cts)
row.names(cts) <- residues

kable(cts, align="c")

## Raw frequency matrix
freq <- cts / apply(cts, 2, sum)
kable(freq, align="c")

## Defining nucleotide-specific residue priors ##
unequal.priors <- c("A"=0.3547, "C"=0.1453, "G"=0.1453, "T"=0.3547)
kable(data.frame(unequal.priors))

## Defining an arbitrary value for the pseudo-weight ##
k <- 4

## Computing corrected frequencies ##
corrected.freq <- (cts + k*unequal.priors) / apply(cts + k*unequal.priors, 2, sum)
kable(corrected.freq, align="c", digits=2, 
      caption=paste("Frequency matrix with unequal priors and pseudo-weight=",k))
weights <- log(corrected.freq / unequal.priors)

sink('PSSM.txt')
kable(weights, align="c", digits=10, caption = paste("Weight matrix with unequal priors and pseudo-weight =",k))
sink()

