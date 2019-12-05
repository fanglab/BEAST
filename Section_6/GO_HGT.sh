#!/bin/bash

#########################################################################
# File Name: GO_HGT.sh
# Author(s): Pedro H. Oliveira
# Institution: Mount Sinai School of Medicine, NY, USA
# Mail: pcphco@gmail.com
# Created Time: Thu 26 Sep 2019 10:54:04 AM (EDT)
#########################################################################


if [[ "$1" == "-h" || "$1" == "-help" ]];
then
	echo "#############################################################################################################################"
	echo "Runs Count (Csuros, 2010) to perform ancestral reconstruction and infer family and lineage specific characteristic along the evolutionary tree."
	echo ""
	echo "Requires: Count, Java"
	echo ""
	echo "Usage: GO_HGT.sh <Pan_genome_matrix.csv> <Newick_tree> <Species_prefix> <Posterior_gain_probability>"
	echo "#############################################################################################################################"
   exit 0
fi


matrix=$1 #Pan_genome_matrix
tree=$2 #Tree in Newick format
prefix=$3 #Species prefix
post_prob=$4 #Posterior gain probability
prefixuc=$(echo "$prefix" | tr '[a-z]' '[A-Z]') #Uppercase species prefix
prefixlc=$(echo "$prefix" | tr '[A-Z]' '[a-z]') #Uppercase species prefix


#########################################################################


echo "Running COUNT under asymmetric Wagner Parsimony. Please wait..."
awk '{for(i=1;++i<=NF-1;) printf $i"\t"; print $(NF)}' $matrix |\
awk '{print "Gene"NR-1, "\t", $0}' | sed 's/Gene0/Family/g' > temp1
java -Xmx2048M -cp Count.jar ca.umontreal.iro.evolution.genecontent.AsymmetricWagner -gain 1 $tree temp1 |\
sed -e '1,/_gain/d' | awk '{if($0 !~ "noname" && $0 !~ "total") print $0}' |\
awk '{print $3, $4}' | sort > $prefix.wagner.out
echo ""
tput setaf 2; echo "Done!"; tput sgr0


#########################################################################


echo "Ancestral reconstruction using posterior probabilities and assuming uniform duplication. Please wait..."
java -Xmx2048M -cp Count.jar ca.umontreal.iro.evolution.genecontent.ML -uniform_duplication true -opt_rounds 100 $tree temp1 > $prefixuc.ML1.r
java -Xmx2048M -cp Count.jar ca.umontreal.iro.evolution.genecontent.Posteriors $tree temp1 $prefixuc.ML1.r > $prefixuc.Posterior.out
rm temp*
echo ""
tput setaf 2; echo "Done!"; tput sgr0


#########################################################################


echo "Parsing Count output. Computing gains at tree tips. Please wait..."


tmp1=$(awk '
   FNR==7{
        for(i=1;i<=NF;i++) 
            if ($i ~ str1 || $i ~ "Family") {
              h=(h)?h FS $i:$i
              f=(f)?f FS i:i
            }
        print h
        nf=split(f,fA,FS);next
   }
   {
       for(i=1;i<=nf;i++) 
           printf("%s%c",$fA[i], (i==nf)?ORS:FS)
   }' str1=$prefixlc FS="\t" $prefixuc.Posterior.out
)



tmp2=$(
echo "$tmp1" |\
awk '
   FNR==1{
        for(i=1;i<=NF;i++) 
            if ($i ~ ":loss" || $i ~ ":gain" || $i ~ "Family") {
              h=(h)?h FS $i:$i
              f=(f)?f FS i:i
            }
        print h
        nf=split(f,fA,FS);next
   }
   {
       for(i=1;i<=nf;i++) 
           printf("%s%c",$fA[i], (i==nf)?ORS:FS)
   }' FS="\t"
)


tmp3=$(echo "$tmp2" | awk '{if($1 !~ "ABSENT") print $0}')



#---Choosing only lines with at least one entry >=post_prob.---
#---Counting the number of fields >post_prob for each line.---
tmp4=$(echo "$tmp3" | tail -n +2 | awk -v var=$post_prob '{ for(i=2; i<=NF; i++) if($i>=var && $i!~"NaN") print $0}' |\
uniq | awk -v var=$post_prob '{ c=1; for(i=2; i<=NF; i++) if($i>=var) c++; print c-1 "\t" $0}')



#---Deleting columns of losses (even columns)---
#---Deleting lines having no gains---
tmp5=$(echo "$tmp4" | awk '{for (i=3; i<NF; i++) printf $i " "; print $NF}' |\
awk '{for (i=1;i<=NF;i+=2) printf("%.5f ", $i); print ""}' |\
awk -v var=$post_prob '{for(i=1; i<=NF; i++)  $i=($i>=var) ? 1 : $i}1' |\
awk -v var=$post_prob '{for(i=1; i<=NF; i++)  $i=($i<var) ? 0 : $i}1')
tmp6=$(echo "$tmp4" | awk '{print $2}')
tmp7=$(paste <(echo "$tmp6") <(echo "$tmp5") | awk '{if($0 !~ "ABSENT") print $0}' |\
awk -v var=$post_prob '{ for(i=2; i<=NF; i++) if($i>=var && $i!~"NaN") print $0}' | uniq)


#---Printing header (only gains) n times---
tmp8=$(echo "$tmp3" |
awk '
   FNR==1{
        for(i=1;i<=NF;i++) 
            if ($i ~ ":gain") {
              h=(h)?h FS $i:$i
              f=(f)?f FS i:i
            }
        print h
        nf=split(f,fA,FS);next
   }
   {
       for(i=1;i<=nf;i++) 
           printf("%s%c",$fA[i], (i==nf)?ORS:FS)
   }' FS="\t" | head -1 | sed 's/:gain//g' | awk '{ for(i=1; i<=NF; i++) print $i}')


#---Computing the sum of gene families acquired at each tree tip---
tmp9=$(echo "$tmp7" | awk '{$1="";print}')
tmp10=$(echo "$tmp9" | awk '{for (i=1;i<=NF;i++) a[i]+=$i}END{for (i=1;i<=NF;i++) printf a[i] OFS; printf "\n"}' | tr ' ' '\n')
paste <(echo "$tmp8") <(echo "$tmp10") | sort > $prefixuc.Gains.out
tput setaf 2; echo "Done!"; tput sgr0


#########################################################################







