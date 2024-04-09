#!/bin/bash

ref_path=... # reference genome path
pig=... # the pig individual under consideration
ppn=... # the number of threads
chrs=`cut -f 1 $ref_path/autosomes_X.txt` # there are two columns in autosomes_X.txt, 1st column being chromosomes, 2nd column being lengths
lens=`cut -f 2 $ref_path/autosomes_X.txt`

# HiFi or ONT
for platform in HiFi ONT;
do
  cd .../${pig}/${platform}
  if [ -d depth_coverage ]; then rm -rf depth_coverage; fi
  mkdir depth_coverage
  parallel "samtools depth -r {} alignment_minimap2_${platform}.sorted.bam > depth_coverage/{}_depth.txt " ::: ${chrs[*]}
  :> depth_coverage/chrs_coverage_depth.txt
  cat $ref_path/autosomes_X.txt | while read i; 
  do 
  	line=($i)
  	chr=${line[0]}
  	len=${line[1]}
  	awk -v chr=$chr -v len=$len 'BEGIN{sum=0;Nr_row=0} {sum+=$3;Nr_row++} END {print chr,sum/Nr_row,Nr_row/len}' depth_coverage/${chr}_depth.txt >> depth_coverage/chrs_coverage_depth.txt
  done
done

# WGS
cd .../${pig}/WGS
if [ -d depth_coverage ]; then rm -rf depth_coverage; fi
mkdir depth_coverage
parallel "samtools depth -r {} ${pig}.sorted.markdup.bam > depth_coverage/{}_depth.txt " ::: ${chrs[*]}
:> depth_coverage/chrs_coverage_depth.txt
cat $ref_path/autosomes_X.txt | while read i;
do 
	line=($i)
	chr=${line[0]}
	len=${line[1]}
	awk -v chr=$chr -v len=$len 'BEGIN{sum=0;Nr_row=0} {sum+=$3;Nr_row++} END {print chr,sum/Nr_row,Nr_row/len}' depth_coverage/${chr}_depth.txt >> depth_coverage/chrs_coverage_depth.txt
done

# Hi-C
cd .../${pig}/Hi-C
bam_path=HiC-Pro_work_dir/output/bowtie_results/bwt2/${pig}/
if [ -d depth_coverage ]; then rm -rf depth_coverage; fi
mkdir depth_coverage
samtools sort -@ $ppn -o $bam_path/alignment.sorted.bam -O bam $bam_path/${pig}_GCF_000003025.6_Sscrofa11.1_genomic.bwt2pairs.bam
samtools index -@ $ppn $bam_path/alignment.sorted.bam
parallel "samtools depth -r {} $bam_path/alignment.sorted.bam > depth_coverage/{}_depth.txt " ::: ${chrs[*]}
:> depth_coverage/chrs_coverage_depth.txt
cat $ref_path/autosomes_X.txt | while read i;
do 
	line=($i)
	chr=${line[0]}
	len=${line[1]}
	awk -v chr=$chr -v len=$len 'BEGIN{sum=0;Nr_row=0} {sum+=$3;Nr_row++} END {print chr,sum/Nr_row,Nr_row/len}' depth_coverage/${chr}_depth.txt >> depth_coverage/chrs_coverage_depth.txt
done
