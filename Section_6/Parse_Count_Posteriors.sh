#!/bin/bash

for z in *.Posterior.out
do 
	
a=$(echo $z)
b=${file/.Posterior.out/}
c=${b:0:4}							 #e.g.: cdif
d="$(echo $c | tr '[a-z]' '[A-Z]')"  #e.g.: CDIF
fam=$(echo Family)
 
tmp1=$(
awk '
   FNR==7{
        for(i=1;i<=NF;i++) 
            if ($i ~ str1 || $i ~ str2) {
              h=(h)?h FS $i:$i
              f=(f)?f FS i:i
            }
        print h
        nf=split(f,fA,FS);next
   }
   {
       for(i=1;i<=nf;i++) 
           printf("%s%c",$fA[i], (i==nf)?ORS:FS)
   }' str1=$c str2=$fam FS="\t" $a
)


tmp2=$(echo "$tmp1" | 
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
   
tmp3=$(echo "$tmp2" | awk '{if($1 != "ABSENT") print $0}')

#---Choosing only lines with at least one entry >=0.2.---
#---Counting the number of fields >0.2 for each line.---
tmp4=$(echo "$tmp3" | tail -n +2 | awk '{ for(i=2; i<=NF; i++) if($i>=0.2 && $i!~"NaN") print $0}' | uniq | awk '{ c=1; for(i=2; i<=NF; i++) if($i>=0.2) c++; print c-1 "\t" $0}')

#---Deleting columns of losses (even columns)---
#---Deleting lines having no gains---
tmp5=$(echo "$tmp4" | awk '{for (i=3; i<NF; i++) printf $i " "; print $NF}' | awk '{for (i=1;i<=NF;i+=2) printf("%.5f ", $i); print ""}' | awk '{for(i=1; i<=NF; i++)  $i=($i>=0.2) ? 1 : $i}1' | awk '{for(i=1; i<=NF; i++)  $i=($i<0.2) ? 0 : $i}1')
tmp6=$(echo "$tmp4" | awk '{print $2}')
tmp7=$(paste <(echo "$tmp6") <(echo "$tmp5") | awk '{if($0 !~ "ABSENT") print $0}' | awk '{ for(i=2; i<=NF; i++) if($i>=0.2 && $i!~"NaN") print $0}' | uniq)


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


#---Computing the sum of each column of genes acquired---
tmp9=$(echo "$tmp7" | awk '{$1="";print}')
tmp10=$(echo "$tmp9" | awk '{for (i=1;i<=NF;i++) a[i]+=$i}END{for (i=1;i<=NF;i++) printf a[i] OFS; printf "\n"}' | tr ' ' '\n')
paste <(echo "$tmp8") <(echo "$tmp10") | sort | tail +2 > $d.Posteriors.out

done

tmp11=$(echo "$tmp8" | while read x; do echo -n "$x " ; done)
tmp12=$(echo "$tmp11" | awk '{for(i=0;i<2580;i++)print}')
tmp13=$(paste <(echo "$tmp7") <(echo "$tmp12"))
tmp14=$(echo "$tmp13" | awk '{print $1,$2":"$47,$3":"$48,$4":"$49,$5":"$50,$6":"$51,$7":"$52,$8":"$53,$9":"$54,$10":"$55,$11":"$56,$12":"$57,$13":"$58,$14":"$59,$15":"$60,$16":"$61,$17":"$62,$18":"$63,$19":"$64,$20":"$65,$21":"$66,$22":"$67,$23":"$68,$24":"$69,$25":"$70,$26":"$71,$27":"$72,$28":"$73,$29":"$74,$30":"$75,$31":"$76,$32":"$77,$33":"$78,$34":"$79,$35":"$80,$36":"$81,$37":"$82,$38":"$83,$39":"$84,$40":"$85,$41":"$86,$42":"$87,$43":"$88,$44":"$89,$45":"$90,$46":"$91}' | awk '{for(i=1;i<=NF;++i)if($i~/1:/)print $1,$i}' | sed 's/1://g' | sort | awk '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x,a[x]}' | sort | awk '{$1=""; print $0}')
tmp15=$(echo "$tmp14" | awk '{print NF}')
tmp16=$(paste <(echo "$tmp15") <(echo "$tmp14") | awk '{if($1>1) print $0}' | awk '{$1=""; print $0}')


echo "$tmp16" | 
while read f; do
	echo $f
		for a in $f; do
    		for b in $f; do
        		printf "%s - %s\n" "$a" "$b"
    		done
		done
done | awk '{if($2=="-" && $3>$1) print $1 "\t" $3}' | sort | uniq -c | awk '{print $2,$3,$1}' | sort > Genes_exchanged_bw_branches.txt


