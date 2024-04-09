#!/bin/bash

curl -O https://s3.amazonaws.com/bedtools-tutorials/web/make-matrix.py
pigs=(AW BKS BMA LW PTR) # pig individuals (breeds)

for pig in ${pigs[*]};
do
  cd .../${pig}/all_callsets
  vcfs=(`ls ./ | grep vcf`) # all short- and long-read-based callsets
  let len_vcfs_minus2=${#vcfs[*]}-2
  let len_vcfs_minus1=${#vcfs[*]}-1
  for i in $(seq 0 1 $len_vcfs_minus2);
  do
  	vcf1=${vcfs[i]}
  	bedtools jaccard -a $vcf1 -b $vcf1 | awk 'NR>1' | cut -f 3 > ${vcf1}.${vcf1}.jaccard
  	let i_plus1=$i+1
  	for j in $(seq $i_plus1 1 $len_vcfs_minus1);
  	do
  		vcf2=${vcfs[j]}
  		bedtools jaccard -a $vcf1 -b $vcf2 | awk 'NR>1' | cut -f 3 > ${vcf1}.${vcf2}.jaccard
  		cp ${vcf1}.${vcf2}.jaccard ${vcf2}.${vcf1}.jaccard
  	done
  done
  vcf1=${vcfs[$len_vcfs_minus1]}
  bedtools jaccard -a $vcf1 -b $vcf1 | awk 'NR>1' | cut -f 3 > ${vcf1}.${vcf1}.jaccard
  
  find . \
  	| grep jaccard \
  	| xargs grep "" \
  	| sed -e s"/\.\///" \
  	| perl -p -e "s/.vcf./.vcf\t/" \
  	| perl -p -e "s/.jaccard:/\t/" \
  	> pairwise.dnase.txt
  	
  cat pairwise.dnase.txt | sed -e 's/_sv.vcf//g' -e 's/.vcf//g' > pairwise.dnase.shortnames.txt
  
  awk 'NF==3' pairwise.dnase.shortnames.txt | python2 ~/make-matrix.py > dnase.shortnames.distance.matrix
  rm .../${pig}/all_callsets/*jaccard
done
