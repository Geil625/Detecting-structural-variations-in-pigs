#!/bin/bash

ref_path=... # reference genome path
pigs=(AW BKS BMA LW PTR) # pig individuals (breeds)
platforms=(HiFi ONT) # long-read sequencing platforms
callers=(DeBreak Sniffles SVIM cuteSV dysgu Picky) # callers
aligners=(minimap2 NGMLR lra) # aligners

###### make benchmark set for each pig individual (breed) ######
for pig in ${pigs[*]};
do
	for platform in ${platforms[*]};
	do
		cd .../$pig/$platform
		# merge alignment-based callsets
		vcfs=$(ls | grep vcf.gz) # all alignment-based callsets
		:> all_callsets.txt
		for vcf in ${vcfs[*]};
		do
			combination=${vcf%.vcf.gz*}
			cp ${combination}.vcf.gz temp.vcf.gz
			bgzip -d temp.vcf.gz
			mv temp.vcf ${combination}.vcf	
			echo `pwd`/${combination}.vcf >> all_callsets.txt
		done
		SURVIVOR merge all_callsets.txt 1000 4 1 1 0 50 SURVIVOR_merged_alignment_callset.vcf
		bcftools annotate --set-id '%CHROM\_%ID\_%POS' SURVIVOR_merged_alignment_callset.vcf > SURVIVOR_merged_alignment_callset_uniqID.vcf
		rm SURVIVOR_merged_alignment_callset.vcf
		ls | grep vcf | grep -v gz | grep -v SURVIVOR | xargs rm
		rm all_callsets.txt
	done

	# merge assembly-based callsets
	cd .../$pig/HiFi
	:> all_callsets.txt
	cp MUMandCo.vcf.gz temp.vcf.gz
	bgzip -d temp.vcf.gz
	mv temp.vcf MUMandCo.vcf
	echo `pwd`/MUMandCo.vcf >> all_callsets.txt
	cp svimasm.vcf.gz temp.vcf.gz
	bgzip -d temp.vcf.gz
	mv temp.vcf svimasm.vcf
	echo `pwd`/svimasm.vcf >> all_callsets.txt
	SURVIVOR merge all_callsets.txt 1000 1 1 1 0 50 SURVIVOR_merged_assembly_callset.vcf
	bcftools annotate --set-id '%CHROM\_%ID\_%POS' SURVIVOR_merged_assembly_callset.vcf > SURVIVOR_merged_assembly_callset_uniqID.vcf
	rm SURVIVOR_merged_assembly_callset.vcf
	ls | grep vcf | grep -v gz | grep -v SURVIVOR | xargs rm
	rm all_callsets.txt

	# merge HiFi and ONT alignment-based sets
	cd .../$pig
	:> all_callsets.txt
	for platform in ${platforms[*]};
	do
		echo `pwd`/${platform}/SURVIVOR_merged_alignment_callset_uniqID.vcf >> `pwd`/$pig/all_callsets.txt
	done
	SURVIVOR merge all_callsets.txt 1000 1 1 1 0 50 SURVIVOR_merged_across_platforms_alignment_set.vcf
	rm all_callsets.txt

	# merge assembly-based and alignment-based sets
	:> all_callsets.txt
	echo `pwd`/HiFi/SURVIVOR_merged_assembly_callset_uniqID.vcf >> all_callsets.txt
	echo `pwd`/SURVIVOR_merged_across_platforms_alignment_set.vcf >> all_callsets.txt
	SURVIVOR merge all_callsets.txt 1000 2 1 1 0 50 SURVIVOR_merged_assembly_alignment_set.vcf
	rm all_callsets.txt
done

# svanalyzer to cluster SVs and remove duplicate SVs
for pig in ${pigs[*]};
do
	cd .../$pig
	svanalyzer merge --ref $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna --variants SURVIVOR_merged_assembly_alignment_set.vcf --prefix svanalyzer_filtered_set
done

# sort and index
conda activate bcftools
for pig in ${pigs[*]};
do
	cd .../$pig
	bcftools sort svanalyzer_filtered_set.clustered.vcf -o benchmark_set.vcf.gz -O z
	tabix -p vcf benchmark_set.vcf.gz
done

###### make hotspot and abundant sets ######
cd ...
:> all_callsets.txt
for pig in ${pigs[*]};
do
	echo `pwd`/$pig/svanalyzer_filtered_set.clustered.vcf >> all_callsets.txt
done
SURVIVOR merge all_callsets.txt 1000 1 1 1 0 50 SURVIVOR_union_acrossPigs.vcf
SURVIVOR merge all_callsets.txt 1000 5 1 1 0 50 SURVIVOR_overlap_acrossPigs.vcf
svanalyzer merge --ref $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna --variants SURVIVOR_union_acrossPigs.vcf --prefix SURVIVOR_union_acrossPigs
svanalyzer merge --ref $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna --variants SURVIVOR_overlap_acrossPigs.vcf --prefix SURVIVOR_overlap_acrossPigs
rm all_callsets.txt

# sort and index
bcftools sort SURVIVOR_union_acrossPigs.clustered.vcf -o hotspot.vcf.gz -O z
bcftools sort SURVIVOR_overlap_acrossPigs.clustered.vcf -o abundant.vcf.gz -O z
tabix -p vcf hotspot.vcf.gz
tabix -p vcf abundant.vcf.gz
