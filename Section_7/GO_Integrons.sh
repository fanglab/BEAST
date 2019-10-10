#!/bin/csh -f

#########################################################################
# File Name: GO_Integron.sh
# Author(s): Pedro H. Oliveira
# Institution: Mount Sinai School of Medicine, NY, USA
# Mail: pcphco@gmail.com
# Created Time: Thu 12 Sep 2019 16:19:10 PM (EDT)
#########################################################################


#Computes Integrons from a FASTA file using Integron Finder.
#Requires: Integron Finder (Cury et al, 2016. Nucleic Acids Res. 44(10): 4539–4550)
#Usage: GO_Integron.sh <*.fasta>
#Example (default): GO_Integron.sh CD630.fasta CRT1.2-CLI.jar 3 19 38 19 48 9


if ("$1" == "-h" || "$1" == "-help") then
	echo ""
	echo "Usage: GO_CRISPRs.sh <*.fasta> <CRT_filename.jar> <minNR> <minRL> <maxRL> <minSL> <maxSL> <searchWL>"
	echo ""
	echo "-minNR: minimum number of repeats a CRISPR must contain; default 3"
	echo "-minRL: minimum length of a CRISPR's repeated region; default 19"
	echo "-maxRL: maximum length of a CRISPR's repeated region; default 38"
	echo "-minSL: minimum length of a CRISPR's non-repeated region (or spacer region); default 19"
	echo "-maxSL: maximum length of a CRISPR's non-repeated region (or spacer region); default 48"
	echo "-searchWL: length of search window used to discover CRISPRs; (range: 6-9)"
	echo ""
    exit 0
endif

set fasta_file = "$1"
set CRT_filename = "$2"
set minNR = "$3"
set minRL = "$4"
set maxRL = "$5"
set minSL = "$6"
set maxSL = "$7"
set searchWL = "$8"
set rootdir = `dirname $0`
set wrk_dir = `cd $rootdir && pwd`


#########################################################################

echo $fasta_file
set b = $fasta_file:r
set c = $CRT_filename:r
if (! -e $b.crispr_raw) then
	touch $b.crispr_raw
    	echo "search crispr in $b"
	java -cp $c.jar crt -minNR $3 -minRL $4 -maxRL $5 -minSL $6 -maxSL $7 -searchWL $8 $b.fasta $b.crispr_raw
endif

echo "Parsing CRT output"

awk '{print FILENAME,$0}' $b.crispr_raw | \
awk '{if($2=="CRISPR") print $1,$2"_"$3,$5,$7; else if($2=="Repeats:") print $1,$3,$6,$9}' |\
awk 'NR%2{printf "%s ",$0;next;}1' | sed 's/.crispr_raw//g' | awk '{print $1,$2,$3,$4,$6,$7,$8}' | sort -k1,1 -k3,3n |\
awk 'BEGIN{print "Genome_ID", "CRISPR_ID", "CRISPR_Start", "CRISPR_END", "Number_repeats", "Av_repeat_size", "Av_spacer_size"}1' > $b.crispr_parsed


tput setaf 2; echo "Done!"; tput sgr0
