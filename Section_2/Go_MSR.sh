#!/bin/csh -f


#########################################################################
# File Name: GO_MSR.sh
# Author(s): Pedro H. Oliveira
# Institution: Mount Sinai School of Medicine, NY, USA
# Mail: pcphco@gmail.com
# Created Time: Thu 19 Sep 2019 10:27:32 AM (EDT)
#########################################################################



if ("$1" == "-h" || "$1" == "-help") then
	echo "#############################################################################################################################"
	echo "Runs the Multi-Scale Representation (MSR) pipeline (Knijnenburg et al., 2014) on a signal file listing all motif positions."
	echo ""
	echo "Requires: SeqKit (Shen et al., 2016), run_msr_runtime_SIGNAL.sh, MATLAB Compiler Runtime (MCR) version 8.1 (R2013a)"
	echo ""
	echo "Usage: GO_MSR.sh <*.fasta> <motif> <mcr_directory> <parameter_file> <MSR_output_filename (same as in parameter file)>"
	echo "#############################################################################################################################"
   exit 0
endif


set fasta_file = "$1"
set motif = "$2"
set motif_up = `echo "$motif" | awk '{print toupper($0)}'`
set pth = "$3"
set parameter = "$4"
set MSR_out = "$5"
set rev_com = `echo "$motif_up" | rev | tr ACGT TGCA`
set chr_size = `awk '/^>/ {if (seqlen) print seqlen;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' $fasta_file`
set motif_size = `expr "${motif}" : '.*'`
set b = $fasta_file:r

#########################################################################

echo "Mapping motif positions in FASTA file..."

awk '{if($0~">") print $0; else print toupper($0)}' $fasta_file |\
seqkit locate -p "$motif_up" |\
awk -v var1=$motif_up -v var2=$rev_com 'NR>1{if($4=="+") print $5,var1; else print $5,var2}' |\
sort | uniq | sort -k1,1n > $b.motif.txt

tput setaf 2; echo "Done!"; tput sgr0


#########################################################################

echo "Building signal file. Please wait..."

awk -v var="$chr_size" 'BEGIN {for(i=1;i<=var;i++) printf "%d ",i; print}' |\
tr ' ' '\n' | awk 'NF' > tmp1.txt

if ("$motif" == "$rev_com") then

	awk -v size=$motif_size '{print $1,$1+size-1}' $b.motif.txt |\
	awk '{for (i=$1; i<=$NF ; i++) print i}' > tmp2.txt
	bash -c 'comm -13 <(sort tmp2.txt) <(sort tmp1.txt)' | awk '{print $1,"0"}' > tmp3.txt
	awk '{print $1,"1"}' tmp2.txt > tmp4.txt
	cat tmp3.txt tmp4.txt | sort -n | awk '{print $2}' > $motif"_"signal.txt
	rm tmp*.txt $b.motif.txt
	
else
	
	awk -v var1="$motif" '{if($2==var1) print $0}' $b.motif.txt > $b.motif"_forward".txt
	awk -v var2="$rev_com" '{if($2==var2) print $0}' $b.motif.txt > $b.motif"_reverse".txt
	
	#Handling forward motif
	awk -v size=$motif_size '{print $1,$1+size-1}' $b.motif"_forward".txt |\
	awk '{for (i=$1; i<=$NF ; i++) print i}' > tmp2.txt
	bash -c 'comm -13 <(sort tmp2.txt) <(sort tmp1.txt)' | awk '{print $1,"0"}' > tmp3.txt
	awk '{print $1,"1"}' tmp2.txt > tmp4.txt
	cat tmp3.txt tmp4.txt | sort -n | awk '{print $2}' > $motif"_forward_"signal.txt

	#Handling reverse motif
	awk -v size=$motif_size '{print $1,$1+size-1}' $b.motif"_reverse".txt |\
	awk '{for (i=$1; i<=$NF ; i++) print i}' > tmp5.txt
	bash -c 'comm -13 <(sort tmp5.txt) <(sort tmp1.txt)' | awk '{print $1,"0"}' > tmp6.txt
	awk '{print $1,"1"}' tmp5.txt > tmp7.txt
	cat tmp6.txt tmp7.txt | sort -n | awk '{print $2}' > $motif"_reverse_"signal.txt
	
	rm tmp*.txt $b.motif.txt $b.motif"_forward".txt $b.motif"_reverse".txt

	echo ""
	tput setaf 2; echo "Done!"; tput sgr0
endif


#########################################################################


echo "Running MSR algorithm. Please wait..."

./run_msr_runtime_SIGNAL.sh $pth $parameter
echo ""
tput setaf 2; echo "Done!"; tput sgr0
echo ""
echo "Parsing MSR output. Please wait..."

awk '{print $1,$2/1,$3/1,$4}' $MSR_out |\
awk '{if($4!="0") print $2,$3}' | awk '{for (i=$1; i<=$NF ; i++) print i}' > tmp1.txt
awk '{print $1,$2/1,$3/1,$4}' $MSR_out |\
awk '{if($4!="0") print $3-$2+1,$1,$4}' | awk '{for(i=1;i<=$1;i++)print $2,$3}' > tmp2.txt
paste tmp1.txt tmp2.txt | awk '{print $2,$1,$3}' |\
awk 'BEGIN{print "MSR_Scale","Position","SFC"}1' > "Pruned_Final_"$motif_up.txt

rm tmp*.txt
echo ""
tput setaf 2; echo "Done!"; tput sgr0


#########################################################################


echo "Plotting MSR graph. Please wait..."
Rscript Plot_MSR.R "Pruned_Final_"$motif_up.txt $chr_size
echo ""
tput setaf 2; echo "Done!"; tput sgr0


#########################################################################

