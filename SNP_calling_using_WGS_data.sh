#!/bin/bash

ref_path=...  # reference genome path
ppn=... # the number of threads
pig=... # the pig individual under consideration

###### bwa-mem2 alignment for WGS adapter sequences removed reads
bwa-mem2 mem -t $ppn -R '@RG\tID:'${pig}'\tPL:illumina\tSM:'${pig}'' ${ref_path}/bwa-mem2/GCF_000003025.6_Sscrofa11.1_genomic ${pig}_1_adprm.fq.gz ${pig}_2_adprm.fq.gz | samtools view -bS - > ${pig}.bam

###### samtools sort and index mapped files
samtools sort -@ $ppn -O bam -o ${pig}.sorted.bam ${pig}.bam

###### mark duplicated reads
picard MarkDuplicates -I ${pig}.sorted.bam -O ${pig}.sorted.markdup.bam -M ${pig}.markdup_metrics.txt
samtools index -@ $ppn ${pig}.sorted.markdup.bam

###### produce gvcf file to call SNPs across pig breeds
gatk --java-options "-Xmx5g" HaplotypeCaller -I ${pig}.sorted.markdup.bam -O ${pig}.g.vcf.gz -R ${ref_path}/GCF_000003025.6_Sscrofa11.1_genomic.fna -ERC GVCF

###### SNP calling across pig breeds based on gvcf files
chrs=(chr1 chr2 ... chr18) # autosomes
for chr in ${chrs[@]};
do
	gatk --java-options "-Xmx30g" GenomicsDBImport \
	 -V AW.g.vcf.gz \
	 -V BKS.g.vcf.gz \
	 -V BMA.g.vcf.gz \
	 -V LW.g.vcf.gz \
	 -V PTR.g.vcf.gz \
	 --genomicsdb-workspace-path my_database \
	 --intervals $chr
	gatk GenotypeGVCFs \
	 -R ${ref_path}/GCF_000003025.6_Sscrofa11.1_genomic.fna \
	 -V gendb://my_database \
	 -G StandardAnnotation \
	 -O fivebreeds_${chr}.vcf
	rm -rf my_database
done

###### make plink files
for chr in ${chrs[@]};
do
  plink --vcf fivebreeds_${chr}.vcf --allow-extra-chr --biallelic-only strict --chr $chr --recode --out fivebreeds_${chr}
done

chrnum=0
for chr in ${chrs[@]};
do
  let chrnum++
  awk -v chrnum=$chrnum '{print chrnum,chrnum":"$4,0,$4}' fivebreeds_${chr}.map > fivebreeds_${chr}_recoded.map
done

echo fivebreeds_chr2.ped fivebreeds_chr2_recoded.map > allchrslist.txt
for i in $(seq 3 1 18);
do
	echo fivebreeds_chr${i}.ped fivebreeds_chr${i}_recoded.map >> allchrslist.txt
done
plink --ped fivebreeds_chr1.ped --map fivebreeds_chr1_recoded.map --merge-list allchrslist.txt --make-bed --out fivebreeds_allchrs
plink --bfile fivebreeds_allchrs --geno 0 --make-bed --out fivebreeds_allchrs_geno0
plink --bfile fivebreeds_allchrs_geno0 --indep-pairwise 50 5 0.1
plink --bfile fivebreeds_allchrs_geno0 --extract plink.prune.in --make-bed --out fivebreeds_allchrs_geno0_LDpruned
plink --bfile fivebreeds_allchrs_geno0_LDpruned --recodeA --out fivebreeds_allchrs_geno0_LDpruned
