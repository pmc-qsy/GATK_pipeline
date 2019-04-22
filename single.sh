#!/usr/bin/bash

# path of certain software and tools
trimmomatic=/home/jiangbx/biosoft/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar
bwa=/home/jiangbx/biosoft/bwa/bwa-0.7.17/bwa
samtools=/home/jiangbx/biosoft/samtools/samtools-1.9/
gatk=/home/jiangbx/biosoft/GATK/gatk-4.1.0.0/gatk

# reference
reference=/your_path_to_reference/hg38/Homo_sampiens_assembly38.fasta
GATK_bundle=/your_path_to_GATK_bundle/hg38

## this step is not neccessary, all depends to whether you have index or not
#$gatk IndexFeatureFile --feature-file $GATK_bundle/hapmap_3.3.hg38.vcf
#$gatk IndexFeatureFile --feature-file $GATK_bundle/1000_omni2.5.hg38.vcf
#$gatk IndexFeatureFile --feature-file $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf
#$gatk IndexFeatureFile --feature-file $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf
#$gatk IndexFeatureFile --feature-file $GATK_bundle/dbsnp_146.hg38.vcf

## shell parameters
fq1=$1
fq2=$2
RGID=$3 ## Read Group, which means Lane ID normally
library=$4 ## ID of seq library
sample=$5 ## sample ID
outdir=$6 ## output path

## set output path according to your sample
outdir=${outdir}/${sample}

## get fastq suffix name from fastq1 file;
## example: *.1.fq.gz
fq_file_name='basement $fq1'
fq_file_name=${fq1_file_name%%.1.fq.gz}

# output path
if [ ! -d $outdir/cleanfq ]
then mkdir -p $outdir/cleanfq
fi

if [ ! -d $outdir/bwa ]
then mkdir -p $outdir/bwa
fi

if [ ! -d $outdir/gatk]
then mkdir -p $outdir/gatk
fi

## quality control by Trimmomatic
## a key parameter in ILLUMINACLIP is keepBothReads=True
time java -jar ${trimmomatic} PE \
	$fq1 $fq2 \
	$outdir/cleanfq/${fq_file_name}.paired.1.fq.gz ${fq_file_name}.unpaired.1.fq.gz \
	$outdir/cleanfq/${fq_file_name}.paired.2.fq.gz ${fq_file_name}.unpaired.2.fq.gz \
	ILLUMINACLIP:/home/jiangbx/biosoft/trimmomatic/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:True \
	SLIDINGWINDOW:5:15 LEADING:5 TRAILING:5 MINLEN:50 && echo "** fq QC done **"

## alignment by bwa mem
time $bwa mem -t 8 -M -Y -R "@RG\tID:$RGID\tPL:ILLUMINA\tPU:$PU\tLB:$library\tSM:$sample" $reference/Homo_sampiens_assembly38.fasta \
	$outdir/cleanfq/${fq_file_name}.paired.1.fq.gz $outdir/cleanfq/${fq_file_name}.paired.2.fq.gz | $samtools view -Sb - > $outdir/bwa/${sample}.bam && \
	echo "** BWA MEM done **" && \
time $samtools sort -@ 4 -m 4G -O bam -o $outdir/bwa/${sample}.sorted.bam  $ outdir/bwa/${sample}.bam && echo "** sorted raw bam file done **"

## this step is not so neccessary 
# time $samtools index $outdir/bwa/${sample}.sorted.bam && echo "** ${sample}.sorted.bam index done **"

## mark duplications
$gatk MarkDuplicates \
	-I $outdir/bwa/${sample}.sorted.bam \
	-M $outdir/bwa/${sample}.markdup_metrics.txt \
	-O $outdir/bwa/${sample}.sorted.markdup.bam && echo "** ${sample}.sorted.bam MarkDuplicates done **"

## make index for ${sample}.sorted.markdup.bam
time $samtools index $outdir/bwa/${sample}.sorted.markdup.bam && echo "** ${sample}.sorted.markdup.bam index done **"

## do BQSR
## p.s. Does your vcf file have an index? GATK4 does not support on the fly indexing of VCFs anymore.
time $gatk BaseRecalibrator \
	-R $reference/Homo_sampiens_assembly38.fasta \
	-I $outdir/bwa/${sample}.sorted.markdup.bam \
	--known-sites $GATK_bundle/1000G_phase1.indels.hg38.vcf \
	--known-sites $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf \
	--known-sites $GATK_bundle/dbsnp_146.hg38.vcf \
	-O $outdir/bwa/${sample}.sorted.markdup.recal_data.table && echo "** ${sample}.sorted.markdup.recal_data.table done **"

time $gatk ApplyBQSR \
	--bqsr-recal-file $outdir/bwa/${sample}.sorted.markdup.recal_data.table \
	-R $reference/Homo_sampiens_assembly38.fasta \
	-I $outdir/bwa/${sample}.sorted.markdup.bam \
	-O $outdir/bwa/${sample}.sorted.markdup.BQSR.bam && echo "** ApplyBQSR done **"

## make index for ${sample}.sorted.markdup.BQSR.bam
time $samtools index $outdir/bwa/${sample}.sorted.markdup.BQSR.bam && echo "** ${sample}.sorted.markdup.BQSR.bam index done **"

## for single sample, we have four methods to call snp with the same results, so just pick one of them.
## Methods one, use HaplotypeCaller to output VCF of your sample, but it can be a bit slow when the input is a bit large.
time $gatk HaplotypeCaller \
	-R $reference/Homo_sampiens_assembly38.fasta \
	-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
	-O $outdir/gatk/${sample}.HC.vcf.gz && echo "** ${sample}.HC.vcf.gz done **"

: '
## Methods two, output every single VCF file for each chromosome, then merge the VCF files. It can speed up a little, but not so neccessary.
chrom={ chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM }
for i in ${chrom[@]}; do
	time $gatk HaplotypeCaller \
	-R $reference/Homo_sampiens_assembly38.fasta \
	-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
	-L $i \
	-O $outdir/gatk/${sample}.HC.${i}.vcf.gz && echo "** ${sample}.HC.${i}.vcf.gz done **"
done && wait
merge_vcfs=""
for i in ${chrom[@]}; do
	merge_vcfs=${merge_vcfs}" -I $outdir/gatk/${sample}.HC.${i}.vcf.gz \\"\n
done && time $gatk MergeVcfs ${merge_vcfs} -O $outdir/gatk/${sample}.HC.vcf.gz && echo "** MergeVcfs done **"

## Methods three, output all samples' gVCF file, then do GenotypeGVCFs. It's not neccessary in single sample, but is a golden standard in multi sample.
time $gatk HaplotypeCaller \
	--emit-ref-confidence GVCF \
	-R $reference/Homo_sampiens_assembly38.fasta \
	-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
	-O $outdir/gatk/${sample}.HC.g.vcf.gz && echo "** GVCF ${sample}.HC.g.vcf.gz done **"
time $gatk GenotypeGVCFs \
	-R $reference/Homo_sampiens_assembly38.fasta \
	-V $outdir/gatk/${sample}.HC.g.vcf.gz \
	-O $outdir/gatk/${sample}.HC.vcf.gz && echo "** ${sample}.HC.vcf.gz done **"

## Methods four, output gvcf of every chromosome, then do GenotypeGVCFs for each of them. This can speed up.
chrom={ chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM }
for i in ${chrom[@]}; do
	time $gatk HaplotypeCaller \
		--emit-ref-confidence GVCF \
		-R $reference/Homo_sampiens_assembly38.fasta \
		-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
		-L $i \
		-O $outdir/gatk/${sample}.HC.${i}.g.vcf.gz && \
	time $gatk GenotypeGVCFs \
		-R $reference/Homo_sampiens_assembly38.fasta \
		-V $outdir/gatk/${sample}.HC.${i}.g.vcf.gz \
		-O $outdir/gatk/${sample}.Hc.${i}.vcf.gz && echo "** ${sample}.Hc.${i}.vcf.gz done **"
'

## VQSR. Seperate two different modes of SNP and Indel due to their different standard for quality.
## SNP mode
time $gatk VariantRecalibrator \
	-R $reference/Homo_sampiens_assembly38.fasta \
	-V $outdir/gatk/${sample}.HC.vcf.gz \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATK_bundle/hapmap_3.3.hg38.vcf \
	-resource:omini,known=false,training=true,truth=false,prior=12.0 $GATK_bundle/1000G_omni2.5.hg38.vcf \
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $GATK_bundle/dbsnp_146.hg38.vcf \
	-an DP -an QD -an FS -an SOR ReadPosRankSum -an MQRankSum \
	-mode SNP \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
	-rscriptFile $outdir/gatk/${sample}.HC.snps.plots.R \
	--tranches-file $outdir/gatk/${sample}.HC.snps.tranches \
	-O $outdir/gatk/${sample}.HC.snps.recal && \
time $gatk ApplyVQSR \
	-R $reference/Homo_sampiens_assembly38.fasta \
	-V $outdir/gatk/${sample}.HC.vcf.gz \
	--ts_filter_level 99.0 \
	--tranches-file $outdir/gatk/${sample}.HC.snps.tranches \
	-recalFile $outdir/gatk/${sample}.HC.snps.recal \
	-mode SNP \
	-O $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz && echo "** SNPs VQSR done **"

## Indel mode
time $gatk VariantRecalibrator \
	-R $reference/Homo_sampiens_assembly38.fasta \
	-input $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz \
	-resource:mills,known=true,training=true,truth=true,prior=12.0 $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf \
	-an DP -an QD -an FS -an SOR ReadPosRankSum -an MQRankSum \
	-mode INDEL \
	--max-gaussians 6 \
	-rscriptFile $outdir/gatk/${sample}.HC.snps.indels.plots.R \
	--tranches-file $outdir/gatk/${sample}.HC.snps.indels.transches \
	-O $outdir/gatk/${sample}.HC.snps.indels.recal && \
time $gatk ApplyVQSR \
	-R $reference/Homo_sampiens_assembly38.fasta \
	-input $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz \
	--ts_filter_level 99.0 \
	--tranches-file $outdir/gatk/${sample}.HC.snps.indels.transches \
	-recalFile $outdir/gatk/${sample}.HC.snps.indels.recal \
	-model INDEL \
	-O $outdir/gatk/${sample}.HC.VQSR.vcf.gz && echo "** SNPs and Indels VQSR done **"

	## remove certain files of middle steps
	# rm -f $outdir/bwa/${sample}.bam $outdir/bwa/${sample}.sorted.bam $outdir/bwa/${sample}.sorted.markdup.bam*
	# rm -f $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz && rm -f $outdir/gatk/${sample}.HC.*.recal

