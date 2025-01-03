#!/bin/bash

ref_path=...  # reference genome path
ppn=... # the number of threads
aligner=... # the aligner under consideration
platform=... # the platform under consideration

############ SV calling based on HiFi or ONT long reads ############

###### integrate HiFi bam data and transform to fastq.gz
bam2fastq -o integrated hifi_reads1.bam ... hifi_readsN.bam

###### integrate ONT fastq.gz data
cat fastq_pass/* > fastq_pass/combined.fastq.gz # only use reads with quality > 7

###### index reference genome for minimap2 and lra
minimap2 -d ${ref_path}/GCF_000003025.6_Sscrofa11.1_genomic.mmi ${ref_path}/minimap2/GCF_000003025.6_Sscrofa11.1_genomic.fna
lra index -ONT GCF_000003025.6_Sscrofa11.1_genomic.fna # [-CCS] output is identical

###### minimap2 alignment respectively for HiFi and ONT reads
minimap2 -t $ppn -ax map-hifi $ref_path/minimap2/GCF_000003025.6_Sscrofa11.1_genomic.mmi integrated.fastq.gz | samtools view -bS - > alignment_minimap2_HiFi.bam
minimap2 -t $ppn -ax map-ont $ref_path/minimap2/GCF_000003025.6_Sscrofa11.1_genomic.mmi combined.fastq.gz | samtools view -bS - > alignment_minimap2_ONT.bam

###### NGMLR alignment respectively for HiFi and ONT reads
ngmlr -t $ppn -r $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna -q integrated.fastq.gz | samtools view -bS - > alignment_ngmlr_HiFi.bam
ngmlr -x ont -t $ppn -r $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna -q combined.fastq.gz | samtools view -bS - > alignment_ngmlr_ONT.bam

###### lra alignment respectively for HiFi and ONT reads
zcat integrated.fastq.gz | lra align -CCS $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna /dev/stdin -t $ppn -p s | samtools view -bS - > alignment_lra_HiFi.bam
zcat combined.fastq.gz | lra align -ONT $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna /dev/stdin -t $ppn -p s | samtools view -bS - > alignment_lra_ONT.bam

###### samtools sort and index reads mapped files
samtools sort -o alignment_${aligner}_${platform}.sorted.bam -@ $ppn alignment_${aligner}_${platform}.bam
samtools sort -n -o alignment_${aligner}_${platform}.sorted_byName.bam -@ $ppn alignment_${aligner}_${platform}.bam
samtools index alignment_${aligner}_${platform}.sorted.bam

###### SV calling by cuteSV, DeBreak, dysgu, Picky, Sniffles, and SVIM respectively using HiFi and ONT reads
# cuteSV
cuteSV alignment_${aligner}_HiFi.sorted.bam $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna cuteSV_${aligner}_HiFi.vcf ./ -s 1 -l 50 --genotype \
 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --threads $ppn
cuteSV alignment_${aligner}_ONT.sorted.bam $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna cuteSV_${aligner}_ONT.vcf ./ -s 5 -l 50 --genotype \
 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --threads $ppn

# debreak
debreak --bam alignment_${aligner}_HiFi.sorted.bam -o debreak_${aligner}_HiFi/ -t $ppn --min_size 50 --min_quality 20 -m 1 --rescue_large_ins --rescue_dup --poa -r $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna
debreak --bam alignment_${aligner}_ONT.sorted.bam -o debreak_${aligner}_ONT/ -t $ppn --min_size 50 --min_quality 20 -m 5 --rescue_large_ins --rescue_dup --poa -r $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna

# dysgu
dysgu call -p $ppn --mode pacbio --min-support 1 --min-size 50 --mq 20 --overwrite --clean $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna ${aligner}_HiFi_temp_dir alignment_${aligner}_HiFi.sorted.bam > dysgu_${aligner}_HiFi.vcf
dysgu call -p $ppn --mode nanopore --min-support 5 --min-size 50 --mq 20 --overwrite --clean $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna ${aligner}_ONT_temp_dir alignment_${aligner}_ONT.sorted.bam > dysgu_${aligner}_ONT.vcf

# Picky
if [ ! -d "picky_${aligner}_${platform}" ]; then
    mkdir picky_${aligner}_${platform}
fi
samtools view -h alignment_${aligner}_${platform}.sorted_byName.bam | picky.pl sam2align | picky.pl callSV --oprefix picky_${aligner}_${platform}
picky.pl xls2vcf --xls picky_${aligner}_${platform}.profile.DEL.xls \
 --xls picky_${aligner}_${platform}.profile.INS.xls \
 --xls picky_${aligner}_${platform}.profile.INDEL.xls \
 --xls picky_${aligner}_${platform}.profile.INV.xls \
 --xls picky_${aligner}_${platform}.profile.TTLC.xls \
 --xls picky_${aligner}_${platform}.profile.TDSR.xls \
 --xls picky_${aligner}_${platform}.profile.TDC.xls > picky_${aligner}_${platform}.allsv.vcf
mv picky_${aligner}_${platform}*.xls picky_${aligner}_${platform}*.vcf picky_${aligner}_${platform}/

# sniffles
sniffles --minsvlen 50 --mapq 20 --minsupport 1 --input alignment_${aligner}_HiFi.sorted.bam --vcf sniffles_${aligner}_HiFi.vcf --threads $ppn --reference $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna
sniffles --minsvlen 50 --mapq 20 --minsupport 5 --input alignment_${aligner}_ONT.sorted.bam --vcf sniffles_${aligner}_ONT.vcf --threads $ppn --reference $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna

# svim
svim alignment --min_mapq 20 --min_sv_size 50 svim_${aligner}_HiFi alignment_${aligner}_HiFi.sorted.bam $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna
svim alignment --min_mapq 20 --min_sv_size 50 svim_${aligner}_ONT alignment_${aligner}_ONT.sorted.bam $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna

###### SVIM-asm using assembled pig genome based on total HiFi data
minimap2 -a -x asm5 --cs -r2k -t $ppn $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna assembled_genome.fa | samtools view -bS - > assembly_alignment_minimap2.bam
samtools sort -@ $ppn -o assembly_alignment_minimap2.sorted.bam assembly_alignment_minimap2.bam
samtools index assembly_alignment_minimap2.sorted.bam
svim-asm haploid svim-asm assembly_alignment_minimap2.sorted.bam $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna

###### MUM&Co using assembled pig genome based on total HiFi data
bash mumandco_v3.8.sh -r $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna -q assembled_genome.fa -g 2501912388 -t $ppn -o MUMandCo


############ SV calling based on Hi-C short reads ############

###### HiC-Pro processes Hi-C data and produces contact matrices
pig_path=... # pig raw Hi-C data path
idx_path=... # bowtie2 index path
sft_path=... # HiC-Pro software path
pig=... # the pig individual under consideration
if [ -d 'HiC-Pro_work_dir' ]; then rm -rf HiC-Pro_work_dir; fi
mkdir HiC-Pro_work_dir
cd HiC-Pro_work_dir
mkdir rawdata
mkdir rawdata/${pig}

cp ${pig_path}/${pig}_1_adprm.fq.gz rawdata/${pig}/${pig}_R1.fastq.gz # using reads with adapter sequences removed
cp ${pig_path}/${pig}_2_adprm.fq.gz rawdata/${pig}/${pig}_R2.fastq.gz # using reads with adapter sequences removed

cp $sft_path/config-hicpro.txt ./

sed -i '/^N_CPU/s#[0-9]\+$#'${ppn}'#' config-hicpro.txt
sed -i '/^MIN_MAPQ/s/[0-9]\+$/20/' config-hicpro.txt
sed -i 's/JOB_NAME =/JOB_NAME = '${pig}'_HiC-Pro/g' config-hicpro.txt
sed -i 's/JOB_MEM =/JOB_MEM = 30gb/g' config-hicpro.txt
sed -i 's/JOB_WALLTIME =/JOB_WALLTIME = 800:00:00/g' config-hicpro.txt
sed -i 's#BOWTIE2_IDX_PATH =#BOWTIE2_IDX_PATH = '${idx_path}'#g' config-hicpro.txt
sed -i 's/REFERENCE_GENOME = hg19/REFERENCE_GENOME = GCF_000003025.6_Sscrofa11.1_genomic/g' config-hicpro.txt
sed -i 's#GENOME_SIZE = chrom_hg19.sizes#GENOME_SIZE = '${ref_path}'/autosomes_X.txt#g' config-hicpro.txt
sed -i 's#GENOME_FRAGMENT = HindIII_resfrag_hg19.bed#GENOME_FRAGMENT = '${idx_path}'/GCF_000003025.6_Sscrofa11.1_genomic.bed#g' config-hicpro.txt
sed -i 's/LIGATION_SITE = AAGCTAGCTT/LIGATION_SITE = GATCGATC/g' config-hicpro.txt
sed -i 's/MIN_FRAG_SIZE =/MIN_FRAG_SIZE = 100/g' config-hicpro.txt
sed -i 's/MAX_FRAG_SIZE =/MAX_FRAG_SIZE = 100000/g' config-hicpro.txt
sed -i 's/MIN_INSERT_SIZE =/MIN_INSERT_SIZE = 100/g' config-hicpro.txt
sed -i 's/MAX_INSERT_SIZE =/MAX_INSERT_SIZE = 1000/g' config-hicpro.txt
sed -i 's/BIN_SIZE = 20000 40000 150000 500000 1000000/BIN_SIZE = 5000 10000 50000 100000 500000 1000000/g' config-hicpro.txt

$sft_path/bin/HiC-Pro -i rawdata -o output -c config-hicpro.txt

###### convert matrix format to cool format
raw_path=.../$pig/HiC/HiC-Pro_work_dir/output/hic_results/matrix/${pig}/raw
for reso in 5000 10000 50000;
do
	hicConvertFormat -m $raw_path/$reso/${pig}_${reso}.matrix \
	 --bedFileHicpro $raw_path/$reso/${pig}_${reso}_abs.bed \
	 --inputFormat hicpro --outputFormat cool -o $raw_path/$reso/matrix.cool
	cooler balance -p $ppn $raw_path/$reso/matrix.cool
done

###### SV calling by EagleC
predictSV --hic-5k $raw_path/5000/matrix.cool \
          --hic-10k $raw_path/10000/matrix.cool \
          --hic-50k $raw_path/50000/matrix.cool \
          -O eaglec_chrnum -g other --balance-type ICE --output-format full


############ SV calling based on WGS short reads ############

###### index reference genome for bwa-mem2
bwa-mem2 index -p ${ref_path}/GCF_000003025.6_Sscrofa11.1_genomic ${ref_path}/GCF_000003025.6_Sscrofa11.1_genomic.fna

###### bwa-mem2 alignment for WGS reads
bwa-mem2 mem -t $ppn -R '@RG\tID:'${pig}'\tPL:illumina\tSM:'${pig}'' ${ref_path}/bwa-mem2/GCF_000003025.6_Sscrofa11.1_genomic ${pig}_1_adprm.fq.gz ${pig}_2_adprm.fq.gz | samtools view -bS - > ${pig}.bam

###### samtools sort and index mapped files
samtools view -@ $ppn -q 20 -O bam -o ${pig}_mapq20.bam ${pig}.bam
samtools sort -@ $ppn -O bam -o ${pig}_mapq20.sorted.bam ${pig}_mapq20.bam

###### mark duplicated reads
picard MarkDuplicates -I ${pig}_mapq20.sorted.bam -O ${pig}_mapq20.sorted.markdup.bam -M ${pig}_mapq20.markdup_metrics.txt
samtools index -@ $ppn ${pig}_mapq20.sorted.markdup.bam
samtools index -@ $ppn ${pig}_mapq20.sorted.bam

###### SV calling respectively by DELLY, LUMPY, Manta, and Wham
# DELLY
delly call -g $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna ${pig}_mapq20.sorted.markdup.bam > delly_mapq20_markdup.vcf

# LUMPY
extractSplitReads_BwaMem_path=.../lumpy-sv-master/scripts
samtools view -@ $ppn -b -F 1294 ${pig}_mapq20.sorted.markdup.bam | samtools sort -@ $ppn - > ${pig}_mapq20.discordants.sorted.markdup.bam
samtools view -@ $ppn -h ${pig}_mapq20.sorted.markdup.bam | $extractSplitReads_BwaMem_path/extractSplitReads_BwaMem -i stdin | samtools view -@ $ppn -Sb - | samtools sort -@ $ppn - > ${pig}_mapq20.splitters.sorted.markdup.bam
lumpyexpress -B ${pig}_mapq20.sorted.markdup.bam -S ${pig}_mapq20.splitters.sorted.markdup.bam -D ${pig}_mapq20.discordants.sorted.markdup.bam -o lumpy_mapq20_markdup.vcf

# Manta
manta_exe_path=.../manta/bin
$manta_exe_path/configManta.py --referenceFasta $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna --bam ${pig}_mapq20.sorted.markdup.bam --runDir manta_dir_mapq20_markdup
./manta_dir_mapq20_markdup/runWorkflow.py -j $ppn

# Wham
whamg -x $ppn -a $ref_path/GCF_000003025.6_Sscrofa11.1_genomic.fna -f ${pig}_mapq20.sorted.markdup.bam > wham_mapq20_markdup.vcf  2> wham_mapq20_markdup.err
