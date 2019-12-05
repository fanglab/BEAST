#!/bin/bash


#########################################################################
# File Name: GO_Prophages_Integrons.sh
# Author(s): Pedro H. Oliveira
# Institution: Mount Sinai School of Medicine, NY, USA
# Mail: pcphco@gmail.com
# Created Time: Mon 30 Sep 2019 16:52:14 PM (EDT)
#########################################################################


if [[ "$1" == "-h" || "$1" == "-help" ]];
then
	echo "#############################################################################################################################"
	echo "Searches Prophages using Phage Finder (Fouts et al, 2006) and Integrons using IntegronFinder (Cury et al., 2016)"
	echo ""
	echo "Requires: PhageFinder, Entrez Direct (EDirect), IntegronFinder"
	echo ""
	echo "Usage: GO_Prophages_Integrons.sh <genome_accession_number> <prefix>"
	echo "#############################################################################################################################"
   exit 0
fi


ACCESSION=$1
PREFIX=$2


#########################################################################

echo "Extracting genome and proteome. Please wait..."
esearch -db nucleotide -query "$1" | efetch -format fasta | awk -v v1=$PREFIX '{if($0~">") print ">"v1; else print $0}' > $PREFIX.con
esearch -db nucleotide -query "$1" | efetch -format fasta_cds_na > temp1
elink -db nucleotide -id "$1" -target protein | efetch -format fasta_cds_aa > temp2
tput setaf 2; echo "Done!"; tput sgr0

#########################################################################

echo "Mapping prophages. Please wait..."
chr_length=`awk '/^>/ {if (seqlen) print seqlen;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' "$PREFIX.con" | uniq`
echo $chr_length
awk -v v1=$PREFIX -v v2=$chr_length '{if($0~">" && $0!~"pseudogene" && $0~"complement" && $0!~"join") print v1, v2, v1":"$(NF-1), $(NF-1), "-"; else if($0~">" && $0!~"pseudogene" && $0!~"complement" && $0!~"join") print v1, v2, v1":"$(NF-1), $(NF-1), "+"}' temp1 |\
sed 's/\[location=complement//g' | sed 's/\]//g' | sed 's/\.\./ /g' | sed 's/)//g' |\
sed 's/(//g' | sed 's/\[location=//g' | awk -v OFS='\t' '{print $1,$2,$3"-"$4,$5,$6,$7}' > phage_finder_info.txt

awk -v v1=$PREFIX '{if($0~">" && $0!~"join") print ">"v1":"$(NF-1); else if($0!~">") print $0}' temp2 |\
sed 's/\[location=complement(//g'  | sed 's/\[location=//g' | sed 's/\]//g' |\
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' |\
sed 's/:/ /g' | sed 's/)/ /g' | sort -k2,2n | sed 's/\.\./-/g' | awk '{print $1":"$2,$3}' |\
awk '{print $1"\n"$2}' | fold -w 70 > $PREFIX.pep

./phage_finder_v2.1.sh $PREFIX
rm temp*
tput setaf 2; echo "Done!"; tput sgr0

#########################################################################

echo "Running IntegronFinder..."

integron_finder $fasta_file --local_max --func_annot

tput setaf 2; echo "Done!"; tput sgr0

#########################################################################
