#!/bin/bash


#########################################################################
# File Name: GO_ConsVar.sh
# Author(s): Pedro H. Oliveira
# Institution: Mount Sinai School of Medicine, NY, USA
# Mail: pcphco@gmail.com
# Created Time: Fri 27 Sep 2019 10:58:41 AM (EDT)
#########################################################################


if [[ "$1" == "-h" || "$1" == "-help" ]];
then
	echo "#############################################################################################################################"
	echo "Performs multiple whole-genome alignment and looks at the conservation of methylation motifs across genomes."
	echo ""
	echo "Requires: ProgressiveMauve, convertAlignment.pl (https://github.com/lskatz), Bedtools, jvarkit (msa2vcf), vcftools, SeqKit"
	echo ""
	echo "Usage: GO_ConsVar.sh <minimal length of LCB> <number of genomes to align> <species_prefix> <MAUVE_DIR> <motif>"
	echo "#############################################################################################################################"
   exit 0
fi

LCB=$1
NUMBERGEN=$2
PREFIX=$3
MAUVE_DIR=$4
motif=$5
motif_up=`echo "$motif" | awk '{print toupper($0)}'`
rev_com=`echo "$motif" | rev | tr ACGT TGCA`

count=`ls -1 *.{fas,fasta} 2>/dev/null | wc -l`
if [ $count != 0 ]; then 
	fasta_list=`echo *.fas`
else
	fasta_list=`echo *.fasta`
fi



#########################################################################


echo "Performing WGA. Running ProgressiveMauve. Please wait..."
progressiveMauve --output=$PREFIX.xmfa $fasta_list
echo ""
tput setaf 2; echo "Done!"; tput sgr0


#########################################################################


echo "Converting xmfa alignment into multifasta. Please wait..."
perl convertAlignment.pl -i $PREFIX.xmfa -o $PREFIX.conv.fas -f fasta -g xmfa -c
echo ""
tput setaf 2; echo "Done!"; tput sgr0


#########################################################################


echo "Converting MSA to VCF format. Please wait..."
java -jar msa2vcf.jar $PREFIX.conv.fas > $PREFIX.vcf
echo ""
tput setaf 2; echo "Done!"; tput sgr0


#########################################################################


echo "Computing all motif positions in genomes. Please wait..."
echo $PREFIX.conv.fas | awk '{if($0~">") print $0; else print toupper($0)}' |\
seqkit locate --ignore-case -p $motif_up $PREFIX.conv.fas |\
sed 's/\// /g' | sed 's/.fa/ /g' |\
awk -v OFS='\t' -v var1=$motif_up -v var2=$rev_com 'NR>1{if($6=="+") print "Motifs",$7,$8,$1"-"var1; else print "Motifs",$7,$8,$1"-"var2}' |\
sort > "All_"$motif_up.txt
echo ""
tput setaf 2; echo "Done!"; tput sgr0


#########################################################################


echo "Computing orthologous conserved / variable motif positions. Please wait..."

awk '{print $2,$3}' "All_"$motif_up.txt | sort | uniq -c |\
awk -v OFS='\t' -v var=$2 '{if($1==var) print $2,$3}' > temp1


vcftools --vcf $PREFIX.vcf --remove-indels --recode --recode-INFO-all --stdout |\
awk -v OFS='\t' 'NR>6{print "Motifs",$2,$2,$4"-"$5}' | sort |\
bedtools intersect -wa -wb -a - -b "All_"$motif_up.txt |\
awk -v OFS='\t' '{print $6,$7}' | sort | uniq | sort > temp2


vcftools --vcf $PREFIX.vcf --keep-only-indels --recode --recode-INFO-all --stdout |\
awk -v OFS='\t' 'NR>6{print "Motifs",$2,$2,$4"-"$5}' | sort |\
bedtools intersect -wa -wb -a - -b "All_"$motif_up.txt |\
awk -v OFS='\t' '{print $6,$7}' | sort | uniq | sort > temp3


tmp1=`comm -12 temp2 temp3`
tmp2=`comm -12 temp1 temp3`
fgrep -x "$tmp1" -v temp3 > temp4
fgrep -x "$tmp2" -v temp4 | awk -v OFS='\t' '{print "Indel",$1,$2}' | sort -k2,2n > $PREFIX.Indels.txt
awk -v OFS='\t' '{print "SNP",$1,$2}' temp2 | sort -k2,2n > $PREFIX.SNPs.txt
awk -v OFS='\t' '{print "Orthologous_Conserved",$1,$2}' temp1 | sort -k2,2n > $PREFIX.Conserved.txt


tmp3=`cat $PREFIX.Conserved.txt $PREFIX.SNPs.txt $PREFIX.Indels.txt | awk '{print $2,$3}' | sort`
awk '{print $2,$3}' "All_"$motif_up.txt | sort | uniq | sort > temp5
fgrep -x "$tmp3" -v temp5 | awk -v OFS='\t' '{print "Non-Orthologous",$1,$2}' | sort -k2,2n > $PREFIX.NonOrthologous.txt

rm *.sslist *.log temp* *.backbone *.bbcols

echo ""
tput setaf 2; echo "Done!"; tput sgr0


#########################################################################

