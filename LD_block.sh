#!/bin/bash

pig=... # the pig individual under consideration
SNPdata=... # the plink .bim .bam .fam files of SNP data

cd .../$pig
plink --bfile $SNPdata \
 --chr 1,4,18,23 \
 --blocks no-pheno-req \
 --blocks-max-kb 1000 \
 --blocks-strong-lowci 0.5001 \
 --blocks-strong-highci 1 \
 --out LD
