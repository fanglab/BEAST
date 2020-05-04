#!/bin/csh -f

#########################################################################
# File Name: GO_CRISPRs.sh
# Author(s): Pedro H. Oliveira
# Institution: Mount Sinai School of Medicine, NY, USA
# Mail: pcphco@gmail.com
# Created Time: Thu 12 Sep 2019 09:31:13 PM (EDT)
#########################################################################


#Computes CRISPRs from a FASTA file using the CRISPR Recognition Tool (CRT) and parses the CRT output (crispr_raw) into a TAB format (crispr_parsed).
#Requires: CRT (Bland et al, 2007. BMC Bioinformatics, 8: 209)
#Usage: GO_CRISPRs.sh <*.fasta> <CRT_filename.jar> <minNR> <minRL> <maxRL> <minSL> <maxSL> <searchWL>
#Example (default): GO_CRISPRs.sh CD630.fasta CRT1.2-CLI.jar 3 19 38 19 48 9


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


#########################################################################

set b = $fasta_file:r
if (! -e $b.crispr_raw) then
	touch $b.crispr_raw
    	echo "search crispr in $b"

	java -cp $CRT_filename crt -minNR $minNR -minRL $minRL -maxRL $maxRL -minSL $minSL -maxSL $maxSL -searchWL $searchWL $fasta_file $b.crispr_raw
endif

echo "Parsing CRT output"

awk '{print FILENAME,$0}' $b.crispr_raw |\
awk '{if($2=="CRISPR") print $1,$2"_"$3,$5,$7; else if($2=="Repeats:") print $1,$3,$6,$9}' |\
awk 'NR%2{printf "%s ",$0;next;}1' | sed 's/.crispr_raw//g' | awk '{print $1,$2,$3,$4,$6,$7,$8}' | sort -k1,1 -k3,3n |\
awk 'BEGIN{print "Genome_ID", "CRISPR_ID", "CRISPR_Start", "CRISPR_END", "Number_repeats", "Av_repeat_size", "Av_spacer_size"}1' > $b.crispr_parsed


tput setaf 2; echo "Done!"; tput sgr0
