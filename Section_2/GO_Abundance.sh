#!/bin/csh -f

#########################################################################
# File Name: GO_Abundance.sh
# Author(s): Pedro H. Oliveira
# Institution: Mount Sinai School of Medicine, NY, USA
# Mail: pcphco@gmail.com
# Created Time: Thu 12 Sep 2019 16:23:10 PM (EDT)
#########################################################################


if ("$1" == "-h" || "$1" == "-help") then
	echo ""
	echo "Computes number of motifs from a FASTA file and their over/under-representation."
	echo ""
	echo "Requires: RMES (Schbath and Hoebeke, 2011. In Advances in genomic sequence analysis and pattern discovery), and SeqKit (Shen et al., 2016)"
	echo ""
	echo "Usage: GO_Abundance.sh <*.fasta> <motif> <word_length> <markov_model_order> <approximation method>"
	echo ""
	echo "word_length: length l of the word to search"
	echo "markov_model_order: From 1 to a maximum of l-2"
	echo "approximation_method: gauss (gaussian), compoundpoisson (Poisson)"
	echo ""
    exit 0
endif

set fasta_file = "$1"
set motif = "$2"
set motif_up = `echo "$motif" | awk '{print toupper($0)}'`
set RMES_l = "$3"
set RMES_w = "$4"
set apx_method = "$5"
set rev_com = `echo "$motif" | rev | tr ACGT TGCA`
set b = $fasta_file:r

#########################################################################


echo "Mapping motif positions in FASTA file..."

awk '{if($0~">") print $0; else print toupper($0)}' $fasta_file |\
seqkit locate -p "$motif_up" |\
awk -v var1=$motif_up -v var2=$rev_com 'NR>1{if($4=="+") print $5,var1; else print $5,var2}' |\
sort | uniq | sort -k1,1n > $b.motif.txt

tput setaf 2; echo "Done!"; tput sgr0


#########################################################################


echo "Running R'MES..."

rmes --$apx_method -s $fasta_file -l $RMES_l -m $set RMES_w -o temp1
rmes.format --tmax 0 --tmin 0 < temp1.0 > $b.RMES.txt
rm temp*

tput setaf 2; echo "Done!"; tput sgr0


#########################################################################

