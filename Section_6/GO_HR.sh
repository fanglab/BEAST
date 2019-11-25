#!/bin/csh -f


#########################################################################
# File Name: GO_HR.sh
# Author(s): Pedro H. Oliveira
# Institution: Mount Sinai School of Medicine, NY, USA
# Mail: pcphco@gmail.com
# Created Time: Thu 24 Sep 2019 17:49:53 PM (EDT)
#########################################################################


if ("$1" == "-h" || "$1" == "-help") then
	echo "#############################################################################################################################"
	echo "Runs ClonalFrameML (Didelot et al. 2015) given an ordered core genome alignment and corresponding phylogenetic tree in Newick format"
	echo ""
	echo "Requires ClonalFrameML, R"
	echo ""
	echo "Usage: GO_HR.sh <ordered_core_alignment> <Newick_tree> <Species_prefix>"
	echo "#############################################################################################################################"
    exit 0
endif


set rootdir = `dirname -- "$0"`
set wrk_dir = `cd $rootdir && pwd`
set alignment = "$1"
set tree = "$2"
set prefix = "$3"
set chr_length = `awk '/^>/ {if (seqlen) print seqlen;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' $alignment | uniq`


#########################################################################


echo "Computing Mean Branch Length and Transition-Transversion ratio of ordered alignment"
echo ""
mkdir align_folder
set align_folder = "$wrk_dir/align_folder"
cp "$1" "$align_folder"
Rscript TTR-MBL.R $align_folder $tree
rm -r $align_folder
echo ""
tput setaf 2; echo "Done!"; tput sgr0


#########################################################################


echo "Running CFML. It may take a while..."
echo "Starting Standard model:"
echo ""
set M = `sed 's/,/ /g' MBL.csv | awk 'NR>1{print $2}'`
set kappa = `sed 's/,/ /g' TTR.csv | awk 'NR>1{print $8}'`
ClonalFrameML $tree $alignment $prefix.smout -kappa $kappa -prior_mean "0.1 0.001 0.1 $M" -emsim 100 > $prefix.smout.log.txt
echo ""
echo "Standard Model finished."
echo ""


echo "Starting Per-branch model:"
echo ""
set outfile = `echo $prefix.smout.em.txt`
set Rtheta = `awk '{if($1=="R/theta") print $2}' $outfile`
set Delta = `awk '{if($1=="1/delta") print $2}' $outfile`
set nu = `awk '{if($1=="nu") print $2}' $outfile`
ClonalFrameML $tree $alignment $prefix.pbmout -kappa $kappa -embranch true -embranch_dispersion 0.1 -initial_values "$Rtheta $Delta $nu" > $prefix.pbmout.log.txt
echo ""
echo "Per-branch model finished."
echo ""
tput setaf 2; echo "Done!"; tput sgr0


#########################################################################


echo "Plotting CFML results:"
echo ""
awk -v var="$chr_length" 'BEGIN {for(i=1;i<=var;i++) printf "%d ",i; print}' |\
tr ' ' '\n' | awk 'NF' > Core_Sites.txt
Rscript cfml_results_2.R $prefix.pbmout Core_Sites.txt
echo ""
tput setaf 2; echo "Done!"; tput sgr0


#########################################################################

