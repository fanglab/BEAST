#!/bin/csh -f

#########################################################################
# File Name: GO_Integrons.sh
# Author(s): Pedro H. Oliveira
# Institution: Mount Sinai School of Medicine, NY, USA
# Mail: pcphco@gmail.com
# Created Time: Thu 12 Sep 2019 16:19:10 PM (EDT)
#########################################################################


#Computes Integrons from a FASTA file using Integron Finder.
#Requires: Integron Finder (Cury et al, 2016. Nucleic Acids Res. 44(10): 4539–4550)



if ("$1" == "-h" || "$1" == "-help") then
	echo ""
	echo "Usage: ./GO_Integrons.sh <*.fasta>"
	echo ""
    exit 0
endif

set fasta_file = "$1"
set rootdir = `dirname -- "$0"``
set wrk_dir = `cd $rootdir && pwd`


#########################################################################

echo "Running IntegronFinder..."

integron_finder $fasta_file --local_max --func_annot

tput setaf 2; echo "Done!"; tput sgr0
