#!/bin/bash

#########################################################################
# File Name: GO_GetCounts.sh
# Author(s): Pedro H. Oliveira
# Institution: Mount Sinai School of Medicine, NY, USA
# Mail: pcphco@gmail.com
# Created Time: Thu 23 May 2019 12:51:54 AM (EDT)
#########################################################################


if [[ "$1" == "-h" || "$1" == "-help" ]];
then
	echo "#############################################################################################################################"
	echo "Performs RNA-seq paired-end read cleaning and mapping for differential expression analysis."
	echo ""
	echo "Requires: Java, bwa, AdapterRemoval, Trimmomatic, SortMeRNA, Samtools, featureCounts"
	echo ""
	echo "Usage: GO_GetCounts.sh <*.fastq1 file> <*.fastq2 file> <*.fasta reference file> <*.SAF annotation file> <PATH to SILVA rRNA_databases folder> <PATH to SILVA index folder> <PATH to Trimmomatic.jar> <PATH to adapters *.fa> <file prefix>"
	echo "#############################################################################################################################"
   exit 0
fi


filein1="$1" # Paired-end forward *.fastq1 file
filein2="$2" # Paired-end reverse *.fastq2 file
filein3="$3" # Reference *.fasta file
filein4="$4" # Annotation file in Simplified Annotation Format (SAF)
DBPATH="$5" #PATH to SILVA fasta database files (rRNA_databases folder).
IDXPATH="$6" #PATH where indexed SILVA database files will be generated.
TRIMPATH="$7" #PATH to Trimmomatic *.jar file.
ADAPATH="$8" #PATH to adapters *.fa file (obtained using AdapterRemoval).
prefix="$9" #file prefix.


######################################################


echo "Running Trimmomatic to eliminate adapters and low quality reads..."
java -jar "$TRIMPATH" PE -phred33 -threads 8 "$filein1" "$filein2" "$prefix"_R1_001_P.fastq "$prefix"_R1_001_U.fastq "$prefix"_R2_001_P.fastq "$prefix"_R2_001_U.fastq ILLUMINACLIP:"$ADAPATH":2:30:10:8:True SLIDINGWINDOW:4:15 LEADING:20 TRAILING:20 MINLEN:50
echo ""
tput setaf 2; echo "Done!"; tput sgr0


######################################################


echo "Running SortMeRNA to eliminate rRNA reads..."
indexdb_rna --ref "$DBPATH"/silva-bac-16s-id90.fasta,"$IDXPATH"/silva-bac-16s-db:\
"$DBPATH"/silva-bac-23s-id98.fasta,"$IDXPATH"/silva-bac-23s-db:\
"$DBPATH"/rfam-5s-database-id98.fasta,"$IDXPATH"/rfam-5s-db -v

./merge-paired-reads.sh "$prefix"_R1_001_P.fastq "$prefix"_R2_001_P.fastq "$prefix"_outfile.fastq #Script merge-paired-reads.sh is included in the SortMeRNA package.

sortmerna --ref "$DBPATH"/silva-bac-16s-id90.fasta,"$IDXPATH"/silva-bac-16s-db:\
"$DBPATH"/silva-bac-23s-id98.fasta,"$IDXPATH"/silva-bac-23s-db:\
"$DBPATH"/rfam-5s-database-id98.fasta,"$IDXPATH"/rfam-5s-db:\
 --reads "$prefix"_outfile.fastq --num_alignments 1 -a 8 --paired_in --fastx --aligned "$prefix"_outfile_rRNA\
 --other "$prefix"_outfile_non_rRNA --log -v
 
./unmerge-paired-reads.sh "$prefix"_outfile_non_rRNA.fastq "$prefix"_R1_001_P_rRNA.fastq "$prefix"_R2_001_P_rRNA.fastq #Script unmerge-paired-reads.sh is included in the SortMeRNA package. 

echo ""
tput setaf 2; echo "Done!"; tput sgr0


######################################################


echo "Read mapping..."
bwa index "$filein3"
bwa mem -t 8 "$filein3" "$prefix"_R1_001_P_rRNA.fastq "$prefix"_R2_001_P_rRNA.fastq > "$prefix".sam
samtools view -bS "$prefix".sam > "$prefix".bam
samtools sort -T temp -o "$prefix".sorted.bam "$prefix".bam
samtools index "$prefix".sorted.bam
featureCounts -p -t exon -g gene_id -T 8 -F SAF -a "$filein4" -o "$prefix".count.txt "$prefix".sorted.bam #Does not consider multi-mapping and multi-overlapping reads.

rm *._U.fastq
rm *_outfile.fastq
rm *.sam
rm $prefix.bam
rm $prefix_outfile_non_rRNA.fastq
rm $prefix_outfile_rRNA.fastq


echo ""
tput setaf 2; echo "Done!"; tput sgr0


######################################################


