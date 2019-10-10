#!/bin/bash

#########################################################################
# File Name: GO_TFBS.sh
# Author(s): Pedro H. Oliveira
# Institution: Mount Sinai School of Medicine, NY, USA
# Mail: pcphco@gmail.com
# Created Time: Thu 21 Sep 2019 17:49:11 PM (EDT)
#########################################################################


if [[ "$1" == "-h" || "$1" == "-help" ]];
then
	echo "#############################################################################################################################"
	echo "Takes TFBS multifasta, computes PSSM, and TFBS hits in a given genome."
	echo ""
	echo "Requires: GO_Build_PWM.R, MAST (MEME Suite), Bioconductor, R"
	echo ""
	echo "Usage: GO_TFBS.sh <TFBS_multifasta> <*.fasta> <TF_name>"
	echo "#############################################################################################################################"
   exit 0
fi


TFBS=$1 #TFBS multifasta
fasta_file=$2 #Genome in fasta format
tf_name=$3 #TF name
rootdir=`dirname -- "$0"`
wrk_dir=`cd $rootdir && pwd`


#########################################################################


echo "Converting TFBS multifasta to a PSSM. Please wait..."
echo ""
Rscript GO_PSSM.R $TFBS
echo ""
tput setaf 2; echo "Done!"; tput sgr0


#########################################################################


echo "Computing TFBS hits. Please wait..."
echo ""
tmp1=`awk 'NR>4{print $0}' PSSM.txt | sed 's/|//g' | sed 's/://g' |\
awk '{$1=""; print $0}' | awk '{for (i=1; i<=NF; i++)  {a[NR,i] = $i}} NF>p { p = NF } END {for(j=1; j<=p; j++) {str=a[1,j]; for(i=2; i<=NR; i++){str=str" "a[i,j];} print str}}'`
cnt=`echo "$tmp1" | wc -l | awk '{print $1}'`
echo "$tmp1" | awk -v var="$cnt" 'BEGIN{print "log-odds matrix: alength= 4 w="" "var}; {print};' |\
awk -v var="$tf_name" 'BEGIN{print "MOTIF", var}1' | awk 'BEGIN{print "ALPHABET= ACGT"}1' |\
awk 'BEGIN{print "MEME version 4"}1' > $tf_name.mast

for f in "$fasta_file"
do
	for g in "$tf_name.mast"
	do
		mast $g $f -hit_list > $tf_name"_TFBS".txt
	done
done

rm PSSM.txt

tput setaf 2; echo "Done!"; tput sgr0


#########################################################################

