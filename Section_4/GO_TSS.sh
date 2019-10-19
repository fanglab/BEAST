#!/bin/bash

#########################################################################
# File Name: GO_TSS.sh
# Author(s): Pedro H. Oliveira
# Institution: Mount Sinai School of Medicine, NY, USA
# Mail: pcphco@gmail.com
# Created Time: Thu 23 Sep 2019 12:10:27 PM (EDT)
#########################################################################


if [[ "$1" == "-h" || "$1" == "-help" ]];
then
	echo "#############################################################################################################################"
	echo "Runs the Parseq program (Mirauta et al., 2014) for reconstruction of the transcriptional landscape from RNA-Seq."
	echo ""
	echo "Requires: Samtools, IGVtools, Parseq, GSL"
	echo ""
	echo "Usage: GO_TSS.sh <bam_file> <chrom_size_file> <*.fas> <counts_folder_path> <results_folder_path> <parseq_folder_path>"
	echo "#############################################################################################################################"
   exit 0
fi



filein1="$1" # *.bam file
filein2="$2" #*.chrom.sizes file in IGV format
filein3="$3" #fasta file
counts_folder="$4"
results_folder="$5"
parseq_folder="$6"
rootdir=`dirname -- "$0"`
wrk_dir=`cd $rootdir && pwd`


#########################################################################


echo "Sorting and indexing bam file. Please wait..."
samtools sort "$filein1" -o "$filein1".sorted.bam
samtools index "$filein1".sorted.bam
echo ""
tput setaf 2; echo "Done!"; tput sgr0


echo "Computing read counts of 5' ends (wig file). Please wait..."
igvtools count -w 1 --postExtFactor 1 --strands read "$filein1".sorted.bam counts.wig "$filein2"
echo ""
tput setaf 2; echo "Done!"; tput sgr0


echo "Running Parseq - Transcription profile reconstruction. Please wait..."
tmp1=`awk '{print $1}' "$filein2"`
mkdir "$4"
mkdir "$5"
cp "${filein3}" *.wig *.chrom.sizes "$4"
Parseq_pmcmc fast "$tmp1" all genome "$wrk_dir" "$5" "$6"
echo ""
tput setaf 2; echo "Done!"; tput sgr0


echo "Running Parseq - Calling transcripts and post-Parseq. Please wait..."
Parseq_particle2proba 0.1 "$tmp1" all genome "$wrk_dir" "$5"
echo ""
tput setaf 2; echo "Done!"; tput sgr0
rm *.sorted bam *.bai

#########################################################################

