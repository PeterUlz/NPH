# Check whether we can make use of plasma-DNA signal coming from protection of nucleosomes

#external programs:

# .) samtools
# .) BWA
# .) R
# .) wgsim
# .) ASCP (Aspera)
# .) fastq-dump (from SRA)
# .) Picard
# .) NPred from NuPOP 

#Steps
# 1) Align plasma-seq samples for several pools
# 2) Trim down plasma-seq alignments to 60bp
# 3) Calculate coverage around various transcription start site (TSS) combinations
# 4) Calculate coverage around various transcription end site (TES) combinations
# 5) Create single-gene wig files for analysis in UCSC Genome Browser

single_genes: output/Single_Genes_Wig/Giant_Plasma_Merge_ERBB2.wig output/Single_Genes_Wig/Merged_Female_Controls_ERBB2.wig \
              output/Single_Genes_Wig/Merged_Male_Controls_ERBB2.wig output/Single_Genes_Wig/ERBB2_positive_ERBB2.wig \
              output/Single_Genes_Wig/Breast_Cancer_ERBB2.wig output/Single_Genes_Wig/Giant_Plasma_Merge_MYC.wig \
              output/Single_Genes_Wig/Merged_Female_Controls_MYC.wig output/Single_Genes_Wig/Merged_Male_Controls_MYC.wig \
              output/Single_Genes_Wig/MYC_positive_MYC.wig output/Single_Genes_Wig/Breast_Cancer_MYC.wig \
              output/Single_Genes_Wig/Prostate_Cancer_MYC.wig output/Single_Genes_Wig/Giant_Plasma_Merge_CCND1.wig \
              output/Single_Genes_Wig/Merged_Female_Controls_CCND1.wig output/Single_Genes_Wig/Merged_Male_Controls_CCND1.wig \
              output/Single_Genes_Wig/CCND1_positive_CCND1.wig output/Single_Genes_Wig/Breast_Cancer_CCND1.wig \
              output/Single_Genes_Wig/Giant_Plasma_Merge_AR.wig output/Single_Genes_Wig/Merged_Female_Controls_AR.wig \
              output/Single_Genes_Wig/Merged_Male_Controls_AR.wig output/Single_Genes_Wig/Breast_Cancer_AR.wig \
              output/Single_Genes_Wig/Prostate_Cancer_AR.wig

single_genes_tes: output/Single_Genes_Wig_TES/Giant_Plasma_Merge_ERBB2.wig output/Single_Genes_Wig_TES/Merged_Female_Controls_ERBB2.wig \
              output/Single_Genes_Wig_TES/Merged_Male_Controls_ERBB2.wig output/Single_Genes_Wig_TES/ERBB2_positive_ERBB2.wig \
              output/Single_Genes_Wig_TES/Breast_Cancer_ERBB2.wig output/Single_Genes_Wig_TES/Giant_Plasma_Merge_MYC.wig \
              output/Single_Genes_Wig_TES/Merged_Female_Controls_MYC.wig output/Single_Genes_Wig_TES/Merged_Male_Controls_MYC.wig \
              output/Single_Genes_Wig_TES/MYC_positive_MYC.wig output/Single_Genes_Wig_TES/Breast_Cancer_MYC.wig \
              output/Single_Genes_Wig_TES/Prostate_Cancer_MYC.wig output/Single_Genes_Wig_TES/Giant_Plasma_Merge_CCND1.wig \
              output/Single_Genes_Wig_TES/Merged_Female_Controls_CCND1.wig output/Single_Genes_Wig_TES/Merged_Male_Controls_CCND1.wig \
              output/Single_Genes_Wig_TES/CCND1_positive_CCND1.wig output/Single_Genes_Wig_TES/Breast_Cancer_CCND1.wig \
              output/Single_Genes_Wig_TES/Giant_Plasma_Merge_AR.wig output/Single_Genes_Wig_TES/Merged_Female_Controls_AR.wig \
              output/Single_Genes_Wig_TES/Merged_Male_Controls_AR.wig output/Single_Genes_Wig_TES/Breast_Cancer_AR.wig \
              output/Single_Genes_Wig_TES/Prostate_Cancer_AR.wig

tissue_specific: output/TSS_coverage/Breast_specific/MergedMale_Controls_Breast_Specific_tss.txt output/TSS_coverage/Breast_specific/MergedFemale_Controls_Breast_Specific_tss.txt \
                 output/TSS_coverage/Breast_specific/Breast_Cancer_Breast_Specific_tss.txt output/TSS_coverage/Breast_specific/Prostate_Cancer_Breast_Specific_tss.txt \
                 output/TSS_coverage/Breast_specific/MergedMale_Controls_Prostate_Specific_tss.txt output/TSS_coverage/Breast_specific/MergedFemale_Controls_Prostate_Specific_tss.txt \
                 output/TSS_coverage/Breast_specific/Breast_Cancer_Prostate_Specific_tss.txt output/TSS_coverage/Breast_specific/Prostate_Cancer_Prostate_Specific_tss.txt

####################################################################################################################################
# 1
# Align Plasma Seq samples
#  -) align to hg19 (bwa mem)
#  -) sample-level duplicate removal using samtools
#  -) merge individual samples using samtools

#Pool of male controls
output/alignments/Merged_Male_Controls_rmdup.bam:
	for i in ./data/MaleControls/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`;bwa mem -t 10 ~/RefSeq/hg19_070510/hg19 $$i > intermediate/MaleControls/$$sample.sam; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/MaleControls/$$sample.sam OUTPUT=intermediate/MaleControls/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/MaleControls/$$sample.sam; \
	samtools rmdup -s intermediate/MaleControls/$$sample.bam intermediate/MaleControls/$$sample.rmdup.bam;rm intermediate/MaleControls/$$sample.bam;rm intermediate/MaleControls/$$sample.bai;samtools index intermediate/MaleControls/$$sample.rmdup.bam; done
	samtools merge -f output/alignments/Merged_Male_Controls_rmdup.bam intermediate/MaleControls/*.rmdup.bam; samtools index output/alignments/Merged_Male_Controls_rmdup.bam

#Pool of female controls
output/alignments/Merged_Female_Controls_rmdup.bam:
	for i in ./data/FemaleControls/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`;bwa mem -t 10 ~/RefSeq/hg19_070510/hg19 $$i > intermediate/FemaleControls/$$sample.sam; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/FemaleControls/$$sample.sam OUTPUT=intermediate/FemaleControls/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/FemaleControls/$$sample.sam; \
	samtools rmdup -s intermediate/FemaleControls/$$sample.bam intermediate/FemaleControls/$$sample.rmdup.bam;rm intermediate/FemaleControls/$$sample.bam;rm intermediate/FemaleControls/$$sample.bai;samtools index intermediate/FemaleControls/$$sample.rmdup.bam; done
	samtools merge -f output/alignments/Merged_Female_Controls_rmdup.bam intermediate/FemaleControls/*.rmdup.bam; samtools index output/alignments/Merged_Female_Controls_rmdup.bam

#Large pool of ~320 samples
output/alignments/Merged_Giant_Plasma_rmdup.bam:
	for i in ./data/Plasma_Giant_Merge/*R1_001.fastq.gz;do base=`basename $$i`;echo "working on $$i" ; \
	sample=`echo $$base | sed s/R1_001\.fastq\.gz//`;bwa mem -t 4 ~/RefSeq/hg19_070510/hg19 $$i> intermediate/Giant_Merge/$$sample.sam; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/Giant_Merge/$$sample.sam OUTPUT=intermediate/Giant_Merge/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/Giant_Merge/$$sample.sam; \
	samtools rmdup -s intermediate/Giant_Merge/$$sample.bam intermediate/Giant_Merge/$$sample.rmdup.bam; \
	rm intermediate/Giant_Merge/$$sample.bam;rm intermediate/Giant_Merge/$$sample.bai;samtools index intermediate/Giant_Merge/$$sample.rmdup.bam; done
	samtools merge -f output/alignments/Merged_Plasma_rmdup.bam intermediate/Giant_Merge/*.rmdup.bam; samtools index output/alignments/Merged_Plasma_rmdup.bam
	samtools merge -f output/alignments/Merged_Giant_Plasma_rmdup.bam output/alignments/Merged_Plasma_rmdup.bam output/alignments/Merged_Paired_rmdup.bam; samtools index output/alignments/Merged_Giant_Plasma_rmdup.bam

#Breast cancer samples
output/alignments/BreastCancer_rmdup.bam:
	for i in ./data/Breast/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`;bwa mem -t 10 ~/RefSeq/hg19_070510/hg19 $$i > intermediate/Breast/$$sample.sam; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/Breast/$$sample.sam OUTPUT=intermediate/Breast/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/Breast/$$sample.sam; \
	samtools rmdup -s intermediate/Breast/$$sample.bam intermediate/Breast/$$sample.rmdup.bam;rm intermediate/Breast/$$sample.bam;rm intermediate/Breast/$$sample.bai;samtools index intermediate/Breast/$$sample.rmdup.bam; done
	samtools merge -f output/alignments/BreastCancer_rmdup.bam intermediate/Breast/*.rmdup.bam; samtools index output/alignments/BreastCancer_rmdup.bam

#Prostate cancer samples
output/alignments/ProstateCancer_rmdup.bam:
	for i in ./data/Prostate/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`;bwa mem -t 10 ~/RefSeq/hg19_070510/hg19 $$i > intermediate/Prostate/$$sample.sam; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/Prostate/$$sample.sam OUTPUT=intermediate/Prostate/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/Prostate/$$sample.sam; \
	samtools rmdup -s intermediate/Prostate/$$sample.bam intermediate/Prostate/$$sample.rmdup.bam;rm intermediate/Prostate/$$sample.bam;rm intermediate/Prostate/$$sample.bai;samtools index intermediate/Prostate/$$sample.rmdup.bam; done
	samtools merge -f output/alignments/ProstateCancer_rmdup.bam intermediate/Prostate/*.rmdup.bam; samtools index output/alignments/ProstateCancer_rmdup.bam

#CNV samples (low-coverage WGS high-molecular-weight DNA)
output/alignments/CNV_samples_rmdup.bam:
	for i in ./data/CNV_samples/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`;bwa mem -t 6 ~/RefSeq/hg19_070510/hg19 $$i > intermediate/CNV_samples/$$sample.sam; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/CNV_samples/$$sample.sam OUTPUT=intermediate/CNV_samples/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/CNV_samples/$$sample.sam; \
	samtools rmdup -s intermediate/CNV_samples/$$sample.bam intermediate/CNV_samples/$$sample.rmdup.bam;rm intermediate/CNV_samples/$$sample.bam;rm intermediate/CNV_samples/$$sample.bai;samtools index intermediate/CNV_samples/$$sample.rmdup.bam; done
	samtools merge -f output/alignments/CNV_samples_rmdup.bam intermediate/CNV_samples/*.rmdup.bam; samtools index output/alignments/CNV_samples_rmdup.bam

#HiSeq NA12878 reads (not plasma-associated; consitutional DNA from HiSeq)
output/alignments/NA12878_HiSeq_rmdup.bam: data/HiSeq_150bp_NA12878/sorted_S1_L001_R1_001.fastq.gz data/HiSeq_150bp_NA12878/sorted_S1_L001_R1_002.fastq.gz 
	cat data/HiSeq_150bp_NA12878/sorted_S1_L001_R1_001.fastq.gz data/HiSeq_150bp_NA12878/sorted_S1_L001_R1_002.fastq.gz > data/HiSeq_150bp_NA12878/NA12878_HiSeq.fastq.gz
	bwa mem -t 10 ~/RefSeq/hg19_070510/hg19 data/HiSeq_150bp_NA12878/NA12878_HiSeq.fastq.gz > intermediate/NA12878_HiSeq/NA12878_HiSeq.sam
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/NA12878_HiSeq/NA12878_HiSeq.sam OUTPUT=intermediate/NA12878_HiSeq/NA12878_HiSeq.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
	rm intermediate/NA12878_HiSeq/NA12878_HiSeq.sam
	samtools rmdup -s intermediate/NA12878_HiSeq/NA12878_HiSeq.bam output/alignments/NA12878_HiSeq_rmdup.bam
	rm intermediate/NA12878_HiSeq/$$sample.bam
	rm intermediate/NA12878_HiSeq/$$sample.bai
	samtools index output/alignments/NA12878_HiSeq_rmdup.bam

#Monophasic breast cancer samples
output/alignments/Monophasic_Breast_rmdup.bam: 
	for i in ./data/Monophasic_Breast/*R1_001.fastq.gz;do base=`basename $$i`;echo "working on $$i" ; \
	sample=`echo $$base | sed s/R1_001\.fastq\.gz//`;bwa mem -t 4 ~/RefSeq/hg19_070510/hg19 $$i > intermediate/Monophasic_Breast/$$sample.sam; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/Monophasic_Breast/$$sample.sam OUTPUT=intermediate/Monophasic_Breast/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/Monophasic_Breast/$$sample.sam; \
	samtools rmdup -s intermediate/Monophasic_Breast/$$sample.bam intermediate/Monophasic_Breast/$$sample.rmdup.bam;rm intermediate/Monophasic_Breast/$$sample.bam;rm intermediate/Monophasic_Breast/$$sample.bai;samtools index intermediate/Monophasic_Breast/$$sample.rmdup.bam; done
	samtools merge -f output/alignments/Monophasic_Breast_rmdup.bam intermediate/Monophasic_Breast/*.rmdup.bam; samtools index output/alignments/Monophasic_Breast_rmdup.bam

#Biphasic breast cancer samples
output/alignments/Biphasic_Breast_rmdup.bam: 
	for i in ./data/Biphasic_Breast/*R1_001.fastq.gz;do base=`basename $$i`;echo "working on $$i" ; \
	sample=`echo $$base | sed s/R1_001\.fastq\.gz//`;bwa mem -t 4 ~/RefSeq/hg19_070510/hg19 $$i > intermediate/Biphasic_Breast/$$sample.sam; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/Biphasic_Breast/$$sample.sam OUTPUT=intermediate/Biphasic_Breast/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/Biphasic_Breast/$$sample.sam; \
	samtools rmdup -s intermediate/Biphasic_Breast/$$sample.bam intermediate/Biphasic_Breast/$$sample.rmdup.bam; \
	rm intermediate/Biphasic_Breast/$$sample.bam;rm intermediate/Biphasic_Breast/$$sample.bai;samtools index intermediate/Biphasic_Breast/$$sample.rmdup.bam; done
	samtools merge -f output/alignments/Biphasic_Breast_rmdup.bam intermediate/Biphasic_Breast/*.rmdup.bam; samtools index output/alignments/Biphasic_Breast_rmdup.bam

#8q overrepresented breast cancer samples
output/alignments/Breast_8q_gain_rmdup.bam: 
	for i in ./data/Breast_8q/*R1_001.fastq.gz;do base=`basename $$i`;echo "working on $$i" ; \
	sample=`echo $$base | sed s/R1_001\.fastq\.gz//`;bwa mem -t 4 ~/RefSeq/hg19_070510/hg19 $$i > intermediate/Breast_8q/$$sample.sam; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/Breast_8q/$$sample.sam OUTPUT=intermediate/Breast_8q/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/Breast_8q/$$sample.sam; \
	samtools rmdup -s intermediate/Breast_8q/$$sample.bam intermediate/Breast_8q/$$sample.rmdup.bam; \
	rm intermediate/Breast_8q/$$sample.bam;rm intermediate/Breast_8q/$$sample.bai;samtools index intermediate/Breast_8q/$$sample.rmdup.bam; done
	samtools merge -f output/alignments/Breast_8q_gain_rmdup.bam intermediate/Breast_8q/*.rmdup.bam; samtools index output/alignments/Breast_8q_gain_rmdup.bam

#ERBB2 amplified samples
output/FocalAmps/ERBB2_rmdup.bam: 
	for i in ./data/FocalAmps/ERBB2/*R1_001.fastq.gz;do base=`basename $$i`;echo "working on $$i" ; \
	sample=`echo $$base | sed s/R1_001\.fastq\.gz//`;bwa mem -t 4 ~/RefSeq/hg19_070510/hg19 $$i > intermediate/FocalAmps/ERBB2/$$sample.sam; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/FocalAmps/ERBB2/$$sample.sam OUTPUT=intermediate/FocalAmps/ERBB2/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/FocalAmps/ERBB2/$$sample.sam; \
	samtools rmdup -s intermediate/FocalAmps/ERBB2/$$sample.bam intermediate/FocalAmps/ERBB2/$$sample.rmdup.bam;rm intermediate/FocalAmps/ERBB2/$$sample.bam;rm intermediate/FocalAmps/ERBB2/$$sample.bai;samtools index intermediate/FocalAmps/ERBB2/$$sample.rmdup.bam; done
	samtools merge -f output/FocalAmps/ERBB2_rmdup.bam intermediate/FocalAmps/ERBB2/*.rmdup.bam; samtools index output/FocalAmps/ERBB2_rmdup.bam

#CCND1 amplified samples
output/FocalAmps/CCND1_rmdup.bam: 
	for i in ./data/FocalAmps/CCND1/*R1_001.fastq.gz;do base=`basename $$i`;echo "working on $$i" ; \
	sample=`echo $$base | sed s/R1_001\.fastq\.gz//`;bwa mem -t 4 ~/RefSeq/hg19_070510/hg19 $$i > intermediate/FocalAmps/CCND1/$$sample.sam; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/FocalAmps/CCND1/$$sample.sam OUTPUT=intermediate/FocalAmps/CCND1/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/FocalAmps/CCND1/$$sample.sam; \
	samtools rmdup -s intermediate/FocalAmps/CCND1/$$sample.bam intermediate/FocalAmps/CCND1/$$sample.rmdup.bam;rm intermediate/FocalAmps/CCND1/$$sample.bam;rm intermediate/FocalAmps/CCND1/$$sample.bai;samtools index intermediate/FocalAmps/CCND1/$$sample.rmdup.bam; done
	samtools merge -f output/FocalAmps/CCND1_rmdup.bam intermediate/FocalAmps/CCND1/*.rmdup.bam; samtools index output/FocalAmps/CCND1_rmdup.bam

#MYC amplified samples
output/FocalAmps/MYC_rmdup.bam: 
	for i in ./data/FocalAmps/MYC/*R1_001.fastq.gz;do base=`basename $$i`;echo "working on $$i" ; \
	sample=`echo $$base | sed s/R1_001\.fastq\.gz//`;bwa mem -t 4 ~/RefSeq/hg19_070510/hg19 $$i > intermediate/FocalAmps/MYC/$$sample.sam; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/FocalAmps/MYC/$$sample.sam OUTPUT=intermediate/FocalAmps/MYC/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/FocalAmps/MYC/$$sample.sam; \
	samtools rmdup -s intermediate/FocalAmps/MYC/$$sample.bam intermediate/FocalAmps/MYC/$$sample.rmdup.bam;rm intermediate/FocalAmps/MYC/$$sample.bam;rm intermediate/FocalAmps/MYC/$$sample.bai;samtools index intermediate/FocalAmps/MYC/$$sample.rmdup.bam; done
	samtools merge -f output/FocalAmps/MYC_rmdup.bam intermediate/FocalAmps/MYC/*.rmdup.bam; samtools index output/FocalAmps/MYC_rmdup.bam


####################################################################################################################################
# 2
# Trim down primary alignments to 60bp to get cleaner nucleosome signal
#
#
output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam:
	samtools view -h output/alignments/Merged_Giant_Plasma_rmdup.bam | ./scripts/trim_alignments_to_60bp.py | samtools view -S -b -o output/trimmed_reads/Merged_Giant_Plasma_trim60bp.bam -
	samtools sort output/trimmed_reads/Merged_Giant_Plasma_trim60bp.bam output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted
	samtools index output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam
	rm output/trimmed_reads/Merged_Giant_Plasma_trim60bp.bam
output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam:
	samtools view -h output/alignments/Merged_Male_Controls_rmdup.bam | ./scripts/trim_alignments_to_60bp.py | samtools view -S -b -o output/trimmed_reads/Merged_Male_Controls_trim60bp.bam -
	samtools sort output/trimmed_reads/Merged_Male_Controls_trim60bp.bam output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted
	samtools index output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam
	rm output/trimmed_reads/Merged_Male_Controls_trim60bp.bam
output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam:
	samtools view -h output/alignments/Merged_Female_Controls_rmdup.bam | ./scripts/trim_alignments_to_60bp.py | samtools view -S -b -o output/trimmed_reads/Merged_Female_Controls_trim60bp.bam -
	samtools sort output/trimmed_reads/Merged_Female_Controls_trim60bp.bam output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted
	samtools index output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam
	rm output/trimmed_reads/Merged_Female_Controls_trim60bp.bam
output/trimmed_reads/BreastCancer_trim60bp.sorted.bam:
	samtools view -h output/alignments/BreastCancer_rmdup.bam | ./scripts/trim_alignments_to_60bp.py | samtools view -S -b -o output/trimmed_reads/BreastCancer_trim60bp.bam -
	samtools sort output/trimmed_reads/BreastCancer_trim60bp.bam output/trimmed_reads/BreastCancer_trim60bp.sorted
	samtools index output/trimmed_reads/BreastCancer_trim60bp.sorted.bam
	rm output/trimmed_reads/BreastCancer_trim60bp.bam
output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam:
	samtools view -h output/alignments/ProstateCancer_rmdup.bam | ./scripts/trim_alignments_to_60bp.py | samtools view -S -b -o output/trimmed_reads/ProstateCancer_trim60bp.bam -
	samtools sort output/trimmed_reads/ProstateCancer_trim60bp.bam output/trimmed_reads/ProstateCancer_trim60bp.sorted
	samtools index output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam
	rm output/trimmed_reads/ProstateCancer_trim60bp.bam
output/trimmed_reads/CNV_samples_trim60bp.sorted.bam:
	samtools view -h output/alignments/CNV_samples_rmdup.bam | ./scripts/trim_alignments_to_60bp.py | samtools view -S -b -o output/trimmed_reads/CNV_samples_trim60bp.bam -
	samtools sort output/trimmed_reads/CNV_samples_trim60bp.bam output/trimmed_reads/CNV_samples_trim60bp.sorted
	samtools index output/trimmed_reads/CNV_samples_trim60bp.sorted.bam
	rm output/trimmed_reads/CNV_samples_trim60bp.bam
output/trimmed_reads/Monophasic_Breast_trim60bp.sorted.bam:
	samtools view -h output/alignments/Monophasic_Breast_rmdup.bam | ./scripts/trim_alignments_to_60bp.py | samtools view -S -b -o output/trimmed_reads/Monophasic_Breast_trim60bp.bam -
	samtools sort output/trimmed_reads/Monophasic_Breast_trim60bp.bam output/trimmed_reads/Monophasic_Breast_trim60bp.sorted
	samtools index output/trimmed_reads/Monophasic_Breast_trim60bp.sorted.bam
	rm output/trimmed_reads/Monophasic_Breast_trim60bp.bam
output/trimmed_reads/Biphasic_Breast_trim60bp.sorted.bam: 
	samtools view -h output/alignments/Biphasic_Breast_rmdup.bam | ./scripts/trim_alignments_to_60bp.py | samtools view -S -b -o output/trimmed_reads/Biphasic_Breast_trim60bp.bam -
	samtools sort output/trimmed_reads/Biphasic_Breast_trim60bp.bam output/trimmed_reads/Biphasic_Breast_trim60bp.sorted
	samtools index output/trimmed_reads/Biphasic_Breast_trim60bp.sorted.bam
	rm output/trimmed_reads/Biphasic_Breast_trim60bp.bam
output/trimmed_reads/Breast_8q_gain_trim60bp.sorted.bam: 
	samtools view -h output/alignments/Breast_8q_gain_rmdup.bam | ./scripts/trim_alignments_to_60bp.py | samtools view -S -b -o output/trimmed_reads/Breast_8q_gain_trim60bp.bam -
	samtools sort output/trimmed_reads/Breast_8q_gain_trim60bp.bam output/trimmed_reads/Breast_8q_gain_trim60bp.sorted
	samtools index output/trimmed_reads/Breast_8q_gain_trim60bp.sorted.bam
	rm output/trimmed_reads/Breast_8q_gain_trim60bp.bam
output/trimmed_reads/ERBB2_trim60bp.sorted.bam: 
	samtools view -h output/FocalAmps/ERBB2_rmdup.bam | ./scripts/trim_alignments_to_60bp.py | samtools view -S -b -o output/trimmed_reads/ERBB2_trim60bp.bam -
	samtools sort output/trimmed_reads/ERBB2_trim60bp.bam output/trimmed_reads/ERBB2_trim60bp.sorted
	samtools index output/trimmed_reads/ERBB2_trim60bp.sorted.bam
	rm output/trimmed_reads/ERBB2_trim60bp.bam
output/trimmed_reads/CCND1_trim60bp.sorted.bam: 
	samtools view -h output/FocalAmps/CCND1_rmdup.bam | ./scripts/trim_alignments_to_60bp.py | samtools view -S -b -o output/trimmed_reads/CCND1_trim60bp.bam -
	samtools sort output/trimmed_reads/CCND1_trim60bp.bam output/trimmed_reads/CCND1_trim60bp.sorted
	samtools index output/trimmed_reads/CCND1_trim60bp.sorted.bam
	rm output/trimmed_reads/CCND1_trim60bp.bam
output/trimmed_reads/MYC_trim60bp.sorted.bam:
	samtools view -h output/FocalAmps/MYC_rmdup.bam | ./scripts/trim_alignments_to_60bp.py | samtools view -S -b -o output/trimmed_reads/MYC_trim60bp.bam -
	samtools sort output/trimmed_reads/MYC_trim60bp.bam output/trimmed_reads/MYC_trim60bp.sorted
	samtools index output/trimmed_reads/MYC_trim60bp.sorted.bam
	rm output/trimmed_reads/MYC_trim60bp.bam
####################################################################################################################################
# 3
# Analyze mean coverage around Transcription start sites (TSS)
#   -) for 10000 genes (first one in refseq files)
#   -) Blood FPKM>5
#   -) Blood Top20
#   -) Blood Top100
#   -) Brain FPKM>5 else FPKM<1
#   
#3.1 10000 genes
output/TSS_coverage/10000genes/Merged_Giant_Plasma_10000genes_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam -t 10 > output/TSS_coverage/10000genes/Merged_Giant_Plasma_10000genes_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/10000genes/Merged_Giant_Plasma_10000genes_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/10000genes/Merged_Giant_Plasma_10000genes_tss.txt
output/TSS_coverage/10000genes/Merged_Male_Controls_10000genes_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -t 10 > output/TSS_coverage/10000genes/Merged_Male_Controls_10000genes_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/10000genes/Merged_Male_Controls_10000genes_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/10000genes/Merged_Male_Controls_10000genes_tss.txt
output/TSS_coverage/10000genes/Merged_Female_Controls_10000genes_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -t 10 > output/TSS_coverage/10000genes/Merged_Female_Controls_10000genes_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/10000genes/Merged_Female_Controls_10000genes_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/10000genes/Merged_Female_Controls_10000genes_tss.txt
output/TSS_coverage/10000genes/BreastCancer_10000genes_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -t 10 > output/TSS_coverage/10000genes/BreastCancer_10000genes_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/10000genes/BreastCancer_10000genes_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/10000genes/BreastCancer_10000genes_tss.txt
output/TSS_coverage/10000genes/ProstateCancer_10000genes_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam -t 10 > output/TSS_coverage/10000genes/ProstateCancer_10000genes_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/10000genes/ProstateCancer_10000genes_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/10000genes/ProstateCancer_10000genes_tss.txt
output/TSS_coverage/10000genes/CNV_samples_10000genes_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/CNV_samples_trim60bp.sorted.bam -t 10 > output/TSS_coverage/10000genes/CNV_samples_10000genes_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/10000genes/CNV_samples_10000genes_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/10000genes/CNV_samples_10000genes_tss.txt
output/TSS_coverage/10000genes/Monophasic_Breast_10000genes_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/Monophasic_Breast_trim60bp.sorted.bam -t 10 > output/TSS_coverage/10000genes/Monophasic_Breast_10000genes_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/10000genes/Monophasic_Breast_10000genes_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/10000genes/Monophasic_Breast_10000genes_tss.txt
output/TSS_coverage/10000genes/Biphasic_Breast_10000genes_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/Biphasic_Breast_trim60bp.sorted.bam -t 10 > output/TSS_coverage/10000genes/Biphasic_Breast_10000genes_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/10000genes/Biphasic_Breast_10000genes_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/10000genes/Biphasic_Breast_10000genes_tss.txt
output/TSS_coverage/10000genes/Breast_8q_gain_10000genes_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/Breast_8q_gain_trim60bp.sorted.bam -t 10 > output/TSS_coverage/10000genes/Breast_8q_gain_10000genes_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/10000genes/Breast_8q_gain_10000genes_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/10000genes/Breast_8q_gain_10000genes_tss.txt
output/TSS_coverage/10000genes/ERBB2_10000genes_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/ERBB2_trim60bp.sorted.bam -t 10 > output/TSS_coverage/10000genes/ERBB2_10000genes_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/10000genes/ERBB2_10000genes_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/10000genes/ERBB2_10000genes_tss.txt
output/TSS_coverage/10000genes/CCND1_10000genes_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/CCND1_trim60bp.sorted.bam -t 10 > output/TSS_coverage/10000genes/CCND1_10000genes_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/10000genes/CCND1_10000genes_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/10000genes/CCND1_10000genes_tss.txt
output/TSS_coverage/10000genes/MYC_10000genes_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/MYC_trim60bp.sorted.bam -t 10 > output/TSS_coverage/10000genes/MYC_10000genes_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/10000genes/MYC_10000genes_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/10000genes/MYC_10000genes_tss.txt

#3.2 Blood_FPKM5
output/TSS_coverage/Blood_FPKM5/Merged_Giant_Plasma_Blood_FPKM5_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_FPKM5/Merged_Giant_Plasma_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_FPKM5/Merged_Giant_Plasma_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_FPKM5/Merged_Giant_Plasma_Blood_FPKM5_tss.txt
output/TSS_coverage/Blood_FPKM5/Merged_Male_Controls_Blood_FPKM5_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_FPKM5/Merged_Male_Controls_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_FPKM5/Merged_Male_Controls_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_FPKM5/Merged_Male_Controls_Blood_FPKM5_tss.txt
output/TSS_coverage/Blood_FPKM5/Merged_Female_Controls_Blood_FPKM5_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_FPKM5/Merged_Female_Controls_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_FPKM5/Merged_Female_Controls_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_FPKM5/Merged_Female_Controls_Blood_FPKM5_tss.txt
output/TSS_coverage/Blood_FPKM5/BreastCancer_Blood_FPKM5_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_FPKM5/BreastCancer_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_FPKM5/BreastCancer_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_FPKM5/BreastCancer_Blood_FPKM5_tss.txt
output/TSS_coverage/Blood_FPKM5/ProstateCancer_Blood_FPKM5_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_FPKM5/ProstateCancer_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_FPKM5/ProstateCancer_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_FPKM5/ProstateCancer_Blood_FPKM5_tss.txt
output/TSS_coverage/Blood_FPKM5/CNV_samples_Blood_FPKM5_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/CNV_samples_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_FPKM5/CNV_samples_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_FPKM5/CNV_samples_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_FPKM5/CNV_samples_Blood_FPKM5_tss.txt
output/TSS_coverage/Blood_FPKM5/Monophasic_Breast_Blood_FPKM5_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Monophasic_Breast_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_FPKM5/Monophasic_Breast_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_FPKM5/Monophasic_Breast_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_FPKM5/Monophasic_Breast_Blood_FPKM5_tss.txt
output/TSS_coverage/Blood_FPKM5/Biphasic_Breast_Blood_FPKM5_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Biphasic_Breast_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_FPKM5/Biphasic_Breast_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_FPKM5/Biphasic_Breast_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_FPKM5/Biphasic_Breast_Blood_FPKM5_tss.txt
output/TSS_coverage/Blood_FPKM5/Breast_8q_gain_Blood_FPKM5_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Breast_8q_gain_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_FPKM5/Breast_8q_gain_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_FPKM5/Breast_8q_gain_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_FPKM5/Breast_8q_gain_Blood_FPKM5_tss.txt
output/TSS_coverage/Blood_FPKM5/ERBB2_Blood_FPKM5_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/ERBB2_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_FPKM5/ERBB2_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_FPKM5/ERBB2_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_FPKM5/ERBB2_Blood_FPKM5_tss.txt
output/TSS_coverage/Blood_FPKM5/CCND1_Blood_FPKM5_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/CCND1_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_FPKM5/CCND1_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_FPKM5/CCND1_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_FPKM5/CCND1_Blood_FPKM5_tss.txt
output/TSS_coverage/Blood_FPKM5/MYC_Blood_FPKM5_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/MYC_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_FPKM5/MYC_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_FPKM5/MYC_Blood_FPKM5_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_FPKM5/MYC_Blood_FPKM5_tss.txt

#3.3 Blood_Top20
output/TSS_coverage/Blood_Top20/Merged_Giant_Plasma_Blood_Top20_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top20/Merged_Giant_Plasma_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top20/Merged_Giant_Plasma_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top20/Merged_Giant_Plasma_Blood_Top20_tss.txt
output/TSS_coverage/Blood_Top20/Merged_Male_Controls_Blood_Top20_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top20/Merged_Male_Controls_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top20/Merged_Male_Controls_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top20/Merged_Male_Controls_Blood_Top20_tss.txt
output/TSS_coverage/Blood_Top20/Merged_Female_Controls_Blood_Top20_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top20/Merged_Female_Controls_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top20/Merged_Female_Controls_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top20/Merged_Female_Controls_Blood_Top20_tss.txt
output/TSS_coverage/Blood_Top20/BreastCancer_Blood_Top20_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top20/BreastCancer_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top20/BreastCancer_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top20/BreastCancer_Blood_Top20_tss.txt
output/TSS_coverage/Blood_Top20/ProstateCancer_Blood_Top20_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top20/ProstateCancer_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top20/ProstateCancer_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top20/ProstateCancer_Blood_Top20_tss.txt
output/TSS_coverage/Blood_Top20/CNV_samples_Blood_Top20_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/CNV_samples_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top20/CNV_samples_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top20/CNV_samples_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top20/CNV_samples_Blood_Top20_tss.txt
output/TSS_coverage/Blood_Top20/Monophasic_Breast_Blood_Top20_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Monophasic_Breast_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top20/Monophasic_Breast_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top20/Monophasic_Breast_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top20/Monophasic_Breast_Blood_Top20_tss.txt
output/TSS_coverage/Blood_Top20/Biphasic_Breast_Blood_Top20_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Biphasic_Breast_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top20/Biphasic_Breast_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top20/Biphasic_Breast_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top20/Biphasic_Breast_Blood_Top20_tss.txt
output/TSS_coverage/Blood_Top20/Breast_8q_gain_Blood_Top20_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Breast_8q_gain_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top20/Breast_8q_gain_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top20/Breast_8q_gain_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top20/Breast_8q_gain_Blood_Top20_tss.txt
output/TSS_coverage/Blood_Top20/ERBB2_Blood_Top20_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/ERBB2_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top20/ERBB2_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top20/ERBB2_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top20/ERBB2_Blood_Top20_tss.txt
output/TSS_coverage/Blood_Top20/CCND1_Blood_Top20_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/CCND1_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top20/CCND1_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top20/CCND1_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top20/CCND1_Blood_Top20_tss.txt
output/TSS_coverage/Blood_Top20/MYC_Blood_Top20_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/MYC_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top20/MYC_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top20/MYC_Blood_Top20_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top20/MYC_Blood_Top20_tss.txt

#3.4 Blood_Top100
output/TSS_coverage/Blood_Top100/Merged_Giant_Plasma_Blood_Top100_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top100/Merged_Giant_Plasma_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top100/Merged_Giant_Plasma_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top100/Merged_Giant_Plasma_Blood_Top100_tss.txt
output/TSS_coverage/Blood_Top100/Merged_Male_Controls_Blood_Top100_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top100/Merged_Male_Controls_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top100/Merged_Male_Controls_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top100/Merged_Male_Controls_Blood_Top100_tss.txt
output/TSS_coverage/Blood_Top100/Merged_Female_Controls_Blood_Top100_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top100/Merged_Female_Controls_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top100/Merged_Female_Controls_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top100/Merged_Female_Controls_Blood_Top100_tss.txt
output/TSS_coverage/Blood_Top100/BreastCancer_Blood_Top100_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top100/BreastCancer_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top100/BreastCancer_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top100/BreastCancer_Blood_Top100_tss.txt
output/TSS_coverage/Blood_Top100/ProstateCancer_Blood_Top100_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top100/ProstateCancer_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top100/ProstateCancer_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top100/ProstateCancer_Blood_Top100_tss.txt
output/TSS_coverage/Blood_Top100/CNV_samples_Blood_Top100_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/CNV_samples_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top100/CNV_samples_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top100/CNV_samples_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top100/CNV_samples_Blood_Top100_tss.txt
output/TSS_coverage/Blood_Top100/Monophasic_Breast_Blood_Top100_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Monophasic_Breast_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top100/Monophasic_Breast_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top100/Monophasic_Breast_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top100/Monophasic_Breast_Blood_Top100_tss.txt
output/TSS_coverage/Blood_Top100/Biphasic_Breast_Blood_Top100_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Biphasic_Breast_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top100/Biphasic_Breast_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top100/Biphasic_Breast_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top100/Biphasic_Breast_Blood_Top100_tss.txt
output/TSS_coverage/Blood_Top100/Breast_8q_gain_Blood_Top100_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Breast_8q_gain_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top100/Breast_8q_gain_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top100/Breast_8q_gain_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top100/Breast_8q_gain_Blood_Top100_tss.txt
output/TSS_coverage/Blood_Top100/ERBB2_Blood_Top100_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/ERBB2_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top100/ERBB2_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top100/ERBB2_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top100/ERBB2_Blood_Top100_tss.txt
output/TSS_coverage/Blood_Top100/CCND1_Blood_Top100_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/CCND1_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top100/CCND1_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top100/CCND1_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top100/CCND1_Blood_Top100_tss.txt
output/TSS_coverage/Blood_Top100/MYC_Blood_Top100_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/MYC_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Blood_Top100/MYC_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Blood_Top100/MYC_Blood_Top100_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Blood_Top100/MYC_Blood_Top100_tss.txt

####################################################################################################################################
# 4
# Analyze mean coverage around Transcription end sites (TES)
#   -) for 10000 genes (first one in refseq files)
#   -) Blood 
#   -) Blood Top20
#   -) Blood Top100
#   -) Brain FPKM>5 else FPKM<1
#3.1 10000 genes
output/TES_coverage/10000genes/Merged_Giant_Plasma_10000genes_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam -t 10 > output/TES_coverage/10000genes/Merged_Giant_Plasma_10000genes_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/10000genes/Merged_Giant_Plasma_10000genes_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/10000genes/Merged_Giant_Plasma_10000genes_TES.txt
output/TES_coverage/10000genes/Merged_Male_Controls_10000genes_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -t 10 > output/TES_coverage/10000genes/Merged_Male_Controls_10000genes_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/10000genes/Merged_Male_Controls_10000genes_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/10000genes/Merged_Male_Controls_10000genes_TES.txt
output/TES_coverage/10000genes/Merged_Female_Controls_10000genes_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -t 10 > output/TES_coverage/10000genes/Merged_Female_Controls_10000genes_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/10000genes/Merged_Female_Controls_10000genes_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/10000genes/Merged_Female_Controls_10000genes_TES.txt
output/TES_coverage/10000genes/BreastCancer_10000genes_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -t 10 > output/TES_coverage/10000genes/BreastCancer_10000genes_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/10000genes/BreastCancer_10000genes_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/10000genes/BreastCancer_10000genes_TES.txt
output/TES_coverage/10000genes/CNV_samples_10000genes_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/CNV_samples_trim60bp.sorted.bam -t 10 > output/TES_coverage/10000genes/CNV_samples_10000genes_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/10000genes/CNV_samples_10000genes_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/10000genes/CNV_samples_10000genes_TES.txt
output/TES_coverage/10000genes/ProstateCancer_10000genes_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam -t 10 > output/TES_coverage/10000genes/ProstateCancer_10000genes_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/10000genes/ProstateCancer_10000genes_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/10000genes/ProstateCancer_10000genes_TES.txt
output/TES_coverage/10000genes/Monophasic_Breast_10000genes_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/Monophasic_Breast_trim60bp.sorted.bam -t 10 > output/TES_coverage/10000genes/Monophasic_Breast_10000genes_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/10000genes/Monophasic_Breast_10000genes_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/10000genes/Monophasic_Breast_10000genes_TES.txt
output/TES_coverage/10000genes/Biphasic_Breast_10000genes_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/Biphasic_Breast_trim60bp.sorted.bam -t 10 > output/TES_coverage/10000genes/Biphasic_Breast_10000genes_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/10000genes/Biphasic_Breast_10000genes_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/10000genes/Biphasic_Breast_10000genes_TES.txt
output/TES_coverage/10000genes/Breast_8q_gain_10000genes_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/Breast_8q_gain_trim60bp.sorted.bam -t 10 > output/TES_coverage/10000genes/Breast_8q_gain_10000genes_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/10000genes/Breast_8q_gain_10000genes_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/10000genes/Breast_8q_gain_10000genes_TES.txt
output/TES_coverage/10000genes/ERBB2_10000genes_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/ERBB2_trim60bp.sorted.bam -t 10 > output/TES_coverage/10000genes/ERBB2_10000genes_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/10000genes/ERBB2_10000genes_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/10000genes/ERBB2_10000genes_TES.txt
output/TES_coverage/10000genes/CCND1_10000genes_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/CCND1_trim60bp.sorted.bam -t 10 > output/TES_coverage/10000genes/CCND1_10000genes_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/10000genes/CCND1_10000genes_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/10000genes/CCND1_10000genes_TES.txt
output/TES_coverage/10000genes/MYC_10000genes_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 10000 -b output/trimmed_reads/MYC_trim60bp.sorted.bam -t 10 > output/TES_coverage/10000genes/MYC_10000genes_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/10000genes/MYC_10000genes_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/10000genes/MYC_10000genes_TES.txt

#3.2 Blood_FPKM5
output/TES_coverage/Blood_FPKM5/Merged_Giant_Plasma_Blood_FPKM5_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_FPKM5/Merged_Giant_Plasma_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_FPKM5/Merged_Giant_Plasma_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_FPKM5/Merged_Giant_Plasma_Blood_FPKM5_TES.txt
output/TES_coverage/Blood_FPKM5/Merged_Male_Controls_Blood_FPKM5_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_FPKM5/Merged_Male_Controls_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_FPKM5/Merged_Male_Controls_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_FPKM5/Merged_Male_Controls_Blood_FPKM5_TES.txt
output/TES_coverage/Blood_FPKM5/Merged_Female_Controls_Blood_FPKM5_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_FPKM5/Merged_Female_Controls_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_FPKM5/Merged_Female_Controls_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_FPKM5/Merged_Female_Controls_Blood_FPKM5_TES.txt
output/TES_coverage/Blood_FPKM5/BreastCancer_Blood_FPKM5_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_FPKM5/BreastCancer_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_FPKM5/BreastCancer_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_FPKM5/BreastCancer_Blood_FPKM5_TES.txt
output/TES_coverage/Blood_FPKM5/ProstateCancer_Blood_FPKM5_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_FPKM5/ProstateCancer_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_FPKM5/ProstateCancer_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_FPKM5/ProstateCancer_Blood_FPKM5_TES.txt
output/TES_coverage/Blood_FPKM5/CNV_samples_Blood_FPKM5_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/CNV_samples_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_FPKM5/CNV_samples_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_FPKM5/CNV_samples_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_FPKM5/CNV_samples_Blood_FPKM5_TES.txt
output/TES_coverage/Blood_FPKM5/Monophasic_Breast_Blood_FPKM5_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Monophasic_Breast_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_FPKM5/Monophasic_Breast_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_FPKM5/Monophasic_Breast_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_FPKM5/Monophasic_Breast_Blood_FPKM5_TES.txt
output/TES_coverage/Blood_FPKM5/Biphasic_Breast_Blood_FPKM5_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Biphasic_Breast_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_FPKM5/Biphasic_Breast_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_FPKM5/Biphasic_Breast_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_FPKM5/Biphasic_Breast_Blood_FPKM5_TES.txt
output/TES_coverage/Blood_FPKM5/Breast_8q_gain_Blood_FPKM5_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Breast_8q_gain_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_FPKM5/Breast_8q_gain_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_FPKM5/Breast_8q_gain_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_FPKM5/Breast_8q_gain_Blood_FPKM5_TES.txt
output/TES_coverage/Blood_FPKM5/ERBB2_Blood_FPKM5_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/ERBB2_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_FPKM5/ERBB2_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_FPKM5/ERBB2_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_FPKM5/ERBB2_Blood_FPKM5_TES.txt
output/TES_coverage/Blood_FPKM5/CCND1_Blood_FPKM5_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/CCND1_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_FPKM5/CCND1_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_FPKM5/CCND1_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_FPKM5/CCND1_Blood_FPKM5_TES.txt
output/TES_coverage/Blood_FPKM5/MYC_Blood_FPKM5_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/MYC_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_FPKM5/MYC_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_FPKM5/MYC_Blood_FPKM5_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_FPKM5/MYC_Blood_FPKM5_TES.txt

#3.3 Blood_Top20
output/TES_coverage/Blood_Top20/Merged_Giant_Plasma_Blood_Top20_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top20/Merged_Giant_Plasma_Blood_Top20_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top20/Merged_Giant_Plasma_Blood_Top20_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top20/Merged_Giant_Plasma_Blood_Top20_TES.txt
output/TES_coverage/Blood_Top20/Merged_Male_Controls_Blood_Top20_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top20/Merged_Male_Controls_Blood_Top20_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top20/Merged_Male_Controls_Blood_Top20_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top20/Merged_Male_Controls_Blood_Top20_TES.txt
output/TES_coverage/Blood_Top20/Merged_Female_Controls_Blood_Top20_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top20/Merged_Female_Controls_Blood_Top20_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top20/Merged_Female_Controls_Blood_Top20_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top20/Merged_Female_Controls_Blood_Top20_TES.txt
output/TES_coverage/Blood_Top20/BreastCancer_Blood_Top20_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top20/BreastCancer_Blood_Top20_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top20/BreastCancer_Blood_Top20_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top20/BreastCancer_Blood_Top20_TES.txt
output/TES_coverage/Blood_Top20/ProstateCancer_Blood_Top20_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top20/ProstateCancer_Blood_Top20_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top20/ProstateCancer_Blood_Top20_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top20/ProstateCancer_Blood_Top20_TES.txt
output/TES_coverage/Blood_Top20/CNV_samples_Blood_Top20_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/CNV_samples_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top20/CNV_samples_Blood_Top20_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top20/CNV_samples_Blood_Top20_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top20/CNV_samples_Blood_Top20_TES.txt
output/TES_coverage/Blood_Top20/Monophasic_Breast_Blood_Top20_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Monophasic_Breast_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top20/Monophasic_Breast_Blood_Top20_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top20/Monophasic_Breast_Blood_Top20_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top20/Monophasic_Breast_Blood_Top20_TES.txt
output/TES_coverage/Blood_Top20/Biphasic_Breast_Blood_Top20_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Biphasic_Breast_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top20/Biphasic_Breast_Blood_Top20_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top20/Biphasic_Breast_Blood_Top20_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top20/Biphasic_Breast_Blood_Top20_TES.txt
output/TES_coverage/Blood_Top20/Breast_8q_gain_Blood_Top20_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Breast_8q_gain_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top20/Breast_8q_gain_Blood_Top20_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top20/Breast_8q_gain_Blood_Top20_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top20/Breast_8q_gain_Blood_Top20_TES.txt
output/TES_coverage/Blood_Top20/ERBB2_Blood_Top20_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/ERBB2_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top20/ERBB2_Blood_Top20_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top20/ERBB2_Blood_Top20_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top20/ERBB2_Blood_Top20_TES.txt
output/TES_coverage/Blood_Top20/CCND1_Blood_Top20_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/CCND1_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top20/CCND1_Blood_Top20_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top20/CCND1_Blood_Top20_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top20/CCND1_Blood_Top20_TES.txt
output/TES_coverage/Blood_Top20/MYC_Blood_Top20_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top20Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/MYC_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top20/MYC_Blood_Top20_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top20/MYC_Blood_Top20_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top20/MYC_Blood_Top20_TES.txt

#3.4 Blood_Top100
output/TES_coverage/Blood_Top100/Merged_Giant_Plasma_Blood_Top100_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top100/Merged_Giant_Plasma_Blood_Top100_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top100/Merged_Giant_Plasma_Blood_Top100_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top100/Merged_Giant_Plasma_Blood_Top100_TES.txt
output/TES_coverage/Blood_Top100/Merged_Male_Controls_Blood_Top100_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top100/Merged_Male_Controls_Blood_Top100_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top100/Merged_Male_Controls_Blood_Top100_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top100/Merged_Male_Controls_Blood_Top100_TES.txt
output/TES_coverage/Blood_Top100/Merged_Female_Controls_Blood_Top100_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top100/Merged_Female_Controls_Blood_Top100_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top100/Merged_Female_Controls_Blood_Top100_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top100/Merged_Female_Controls_Blood_Top100_TES.txt
output/TES_coverage/Blood_Top100/BreastCancer_Blood_Top100_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top100/BreastCancer_Blood_Top100_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top100/BreastCancer_Blood_Top100_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top100/BreastCancer_Blood_Top100_TES.txt
output/TES_coverage/Blood_Top100/ProstateCancer_Blood_Top100_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top100/ProstateCancer_Blood_Top100_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top100/ProstateCancer_Blood_Top100_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top100/ProstateCancer_Blood_Top100_TES.txt
output/TES_coverage/Blood_Top100/CNV_samples_Blood_Top100_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/CNV_samples_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top100/CNV_samples_Blood_Top100_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top100/CNV_samples_Blood_Top100_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top100/CNV_samples_Blood_Top100_TES.txt
output/TES_coverage/Blood_Top100/Monophasic_Breast_Blood_Top100_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Monophasic_Breast_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top100/Monophasic_Breast_Blood_Top100_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top100/Monophasic_Breast_Blood_Top100_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top100/Monophasic_Breast_Blood_Top100_TES.txt
output/TES_coverage/Blood_Top100/Biphasic_Breast_Blood_Top100_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Biphasic_Breast_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top100/Biphasic_Breast_Blood_Top100_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top100/Biphasic_Breast_Blood_Top100_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top100/Biphasic_Breast_Blood_Top100_TES.txt
output/TES_coverage/Blood_Top100/Breast_8q_gain_Blood_Top100_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Breast_8q_gain_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top100/Breast_8q_gain_Blood_Top100_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top100/Breast_8q_gain_Blood_Top100_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top100/Breast_8q_gain_Blood_Top100_TES.txt
output/TES_coverage/Blood_Top100/ERBB2_Blood_Top100_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/ERBB2_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top100/ERBB2_Blood_Top100_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top100/ERBB2_Blood_Top100_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top100/ERBB2_Blood_Top100_TES.txt
output/TES_coverage/Blood_Top100/CCND1_Blood_Top100_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/CCND1_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top100/CCND1_Blood_Top100_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top100/CCND1_Blood_Top100_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top100/CCND1_Blood_Top100_TES.txt
output/TES_coverage/Blood_Top100/MYC_Blood_Top100_TES.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/GTEx/Top100Blood_noChrM.csv -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/MYC_trim60bp.sorted.bam -t 10 > output/TES_coverage/Blood_Top100/MYC_Blood_Top100_TES.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/Blood_Top100/MYC_Blood_Top100_TES.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/Blood_Top100/MYC_Blood_Top100_TES.txt


####################################################################################################################################
# 5
# Create WIG files for single genes from various BAM files for TSS
#   

#ERBB2
output/Single_Genes_Wig/Giant_Plasma_Merge_ERBB2.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Giant_Plasma_Merge_ERBB2.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam -g ERBB2
output/Single_Genes_Wig/Merged_Female_Controls_ERBB2.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Female_Controls_ERBB2.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -g ERBB2
output/Single_Genes_Wig/Merged_Male_Controls_ERBB2.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Male_Controls_ERBB2.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -g ERBB2
output/Single_Genes_Wig/Merged_Male_Controls_ERBB2_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Male_Controls_ERBB2_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -g ERBB2
output/Single_Genes_Wig/Merged_Male_Controls_ERBB2_large_smoothed_400bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Male_Controls_ERBB2_large.wig -o output/Single_Genes_Wig/Merged_Male_Controls_ERBB2_large_smoothed_400bp.wig -w 400
output/Single_Genes_Wig/Merged_Male_Controls_ERBB2_large_smoothed_1000bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Male_Controls_ERBB2_large.wig -o output/Single_Genes_Wig/Merged_Male_Controls_ERBB2_large_smoothed_1000bp.wig -w 1000

output/Single_Genes_Wig/ERBB2_positive_ERBB2.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/ERBB2_positive_ERBB2.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/ERBB2_trim60bp.sorted.bam -g ERBB2
output/Single_Genes_Wig/ERBB2_positive_ERBB2_smoothed_10window.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/ERBB2_positive_ERBB2.wig -o output/Single_Genes_Wig/ERBB2_positive_ERBB2_smoothed_10window.wig
output/Single_Genes_Wig/ERBB2_positive_ERBB2_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/ERBB2_positive_ERBB2_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/ERBB2_trim60bp.sorted.bam -g ERBB2
output/Single_Genes_Wig/ERBB2_positive_ERBB2_large_smoothed_400bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/ERBB2_positive_ERBB2_large.wig -o output/Single_Genes_Wig/ERBB2_positive_ERBB2_large_smoothed_400bp.wig -w 400
output/Single_Genes_Wig/ERBB2_positive_ERBB2_large_smoothed_1000bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/ERBB2_positive_ERBB2_large.wig -o output/Single_Genes_Wig/ERBB2_positive_ERBB2_large_smoothed_1000bp.wig -w 1000

output/Single_Genes_Wig/Breast_Cancer_ERBB2.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Breast_Cancer_ERBB2.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -g ERBB2

#MYC
output/Single_Genes_Wig/Giant_Plasma_Merge_MYC.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Giant_Plasma_Merge_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam -g MYC
output/Single_Genes_Wig/Merged_Female_Controls_MYC.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Female_Controls_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -g MYC
output/Single_Genes_Wig/Merged_Male_Controls_MYC.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Male_Controls_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -g MYC
output/Single_Genes_Wig/MYC_positive_MYC.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/MYC_positive_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/MYC_trim60bp.sorted.bam -g MYC
output/Single_Genes_Wig/Breast_Cancer_MYC.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Breast_Cancer_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -g MYC
output/Single_Genes_Wig/Prostate_Cancer_MYC.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Prostate_Cancer_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam -g MYC

#CCND1
output/Single_Genes_Wig/Giant_Plasma_Merge_CCND1.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Giant_Plasma_Merge_CCND1.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam -g CCND1
output/Single_Genes_Wig/Merged_Female_Controls_CCND1.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Female_Controls_CCND1.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -g CCND1
output/Single_Genes_Wig/Merged_Male_Controls_CCND1.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Male_Controls_CCND1.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -g CCND1
output/Single_Genes_Wig/Merged_Male_Controls_CCND1_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Male_Controls_CCND1_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -g CCND1
output/Single_Genes_Wig/Merged_Male_Controls_CCND1_large_smoothed_400bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Male_Controls_CCND1_large.wig -o output/Single_Genes_Wig/Merged_Male_Controls_CCND1_large_smoothed_400bp.wig -w 400
output/Single_Genes_Wig/Merged_Male_Controls_CCND1_large_smoothed_1000bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Male_Controls_CCND1_large.wig -o output/Single_Genes_Wig/Merged_Male_Controls_CCND1_large_smoothed_1000bp.wig -w 1000

output/Single_Genes_Wig/CCND1_positive_CCND1.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/CCND1_positive_CCND1.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/CCND1_trim60bp.sorted.bam -g CCND1
output/Single_Genes_Wig/CCND1_positive_CCND1_smoothed_10window.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/CCND1_positive_CCND1.wig -o output/Single_Genes_Wig/CCND1_positive_CCND1_smoothed_10window.wig  
output/Single_Genes_Wig/CCND1_positive_CCND1_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/CCND1_positive_CCND1_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/CCND1_trim60bp.sorted.bam -g CCND1
output/Single_Genes_Wig/CCND1_positive_CCND1_large_smoothed_400bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/CCND1_positive_CCND1_large.wig -o output/Single_Genes_Wig/CCND1_positive_CCND1_large_smoothed_400bp.wig -w 400
output/Single_Genes_Wig/CCND1_positive_CCND1_large_smoothed_1000bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/CCND1_positive_CCND1_large.wig -o output/Single_Genes_Wig/CCND1_positive_CCND1_large_smoothed_1000bp.wig -w 1000

output/Single_Genes_Wig/Breast_Cancer_CCND1.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Breast_Cancer_CCND1.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -g CCND1

#AR
output/Single_Genes_Wig/Giant_Plasma_Merge_AR.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Giant_Plasma_Merge_AR.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam -g AR
output/Single_Genes_Wig/Merged_Female_Controls_AR.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Female_Controls_AR.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -g AR
output/Single_Genes_Wig/Merged_Female_Controls_AR_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Female_Controls_AR_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -g AR
output/Single_Genes_Wig/Merged_Female_Controls_AR_large_smoothed_400bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Female_Controls_AR_large.wig -o output/Single_Genes_Wig/Merged_Female_Controls_AR_large_smoothed_400bp.wig -w 400
output/Single_Genes_Wig/Merged_Female_Controls_AR_large_smoothed_1000bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Female_Controls_AR_large.wig -o output/Single_Genes_Wig/Merged_Female_Controls_AR_large_smoothed_1000bp.wig -w 1000
output/Single_Genes_Wig/Merged_Male_Controls_AR.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Male_Controls_AR.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -g AR
output/Single_Genes_Wig/Breast_Cancer_AR.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Breast_Cancer_AR.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -g AR
output/Single_Genes_Wig/Prostate_Cancer_AR.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Prostate_Cancer_AR.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam -g AR
output/Single_Genes_Wig/Prostate_Cancer_AR_smoothed_10window.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Prostate_Cancer_AR.wig -o output/Single_Genes_Wig/Prostate_Cancer_AR_smoothed_10window.wig
output/Single_Genes_Wig/Prostate_Cancer_AR_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Prostate_Cancer_AR_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam -g AR
output/Single_Genes_Wig/Prostate_Cancer_AR_large_smoothed_400bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Prostate_Cancer_AR_large.wig -o output/Single_Genes_Wig/Prostate_Cancer_AR_large_smoothed_400bp.wig -w 400
output/Single_Genes_Wig/Prostate_Cancer_AR_large_smoothed_1000bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Prostate_Cancer_AR_large.wig -o output/Single_Genes_Wig/Prostate_Cancer_AR_large_smoothed_1000bp.wig -w 1000


#HBB
output/Single_Genes_Wig/Merged_Female_Controls_HBB_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Female_Controls_HBB_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -g HBB
output/Single_Genes_Wig/Merged_Female_Controls_HBB_large_smoothed_400bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Female_Controls_HBB_large.wig -o output/Single_Genes_Wig/Merged_Female_Controls_HBB_large_smoothed_400bp.wig -w 400
output/Single_Genes_Wig/Merged_Female_Controls_HBB_large_smoothed_1000bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Female_Controls_HBB_large.wig -o output/Single_Genes_Wig/Merged_Female_Controls_HBB_large_smoothed_1000bp.wig -w 1000

#HBA2
output/Single_Genes_Wig/Merged_Female_Controls_HBA2_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Female_Controls_HBA2_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -g HBA2
output/Single_Genes_Wig/Merged_Female_Controls_HBA2_large_smoothed_400bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Female_Controls_HBA2_large.wig -o output/Single_Genes_Wig/Merged_Female_Controls_HBA2_large_smoothed_400bp.wig -w 400
output/Single_Genes_Wig/Merged_Female_Controls_HBA2_large_smoothed_1000bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Female_Controls_HBA2_large.wig -o output/Single_Genes_Wig/Merged_Female_Controls_HBA2_large_smoothed_1000bp.wig -w 1000

#HBA1
output/Single_Genes_Wig/Merged_Female_Controls_HBA1_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Female_Controls_HBA1_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -g HBA1
output/Single_Genes_Wig/Merged_Female_Controls_HBA1_large_smoothed_400bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Female_Controls_HBA1_large.wig -o output/Single_Genes_Wig/Merged_Female_Controls_HBA1_large_smoothed_400bp.wig -w 400
output/Single_Genes_Wig/Merged_Female_Controls_HBA1_large_smoothed_1000bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Female_Controls_HBA1_large.wig -o output/Single_Genes_Wig/Merged_Female_Controls_HBA1_large_smoothed_1000bp.wig -w 1000

#S100A9
output/Single_Genes_Wig/Merged_Female_Controls_S100A9_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Female_Controls_S100A9_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -g S100A9
output/Single_Genes_Wig/Merged_Female_Controls_S100A9_large_smoothed_400bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Female_Controls_S100A9_large.wig -o output/Single_Genes_Wig/Merged_Female_Controls_S100A9_large_smoothed_400bp.wig -w 400
output/Single_Genes_Wig/Merged_Female_Controls_S100A9_large_smoothed_1000bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Female_Controls_S100A9_large.wig -o output/Single_Genes_Wig/Merged_Female_Controls_S100A9_large_smoothed_1000bp.wig -w 1000

#FTL
output/Single_Genes_Wig/Merged_Female_Controls_FTL_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Female_Controls_FTL_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -g FTL
output/Single_Genes_Wig/Merged_Female_Controls_FTL_large_smoothed_400bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Female_Controls_FTL_large.wig -o output/Single_Genes_Wig/Merged_Female_Controls_FTL_large_smoothed_400bp.wig -w 400
output/Single_Genes_Wig/Merged_Female_Controls_FTL_large_smoothed_1000bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Female_Controls_FTL_large.wig -o output/Single_Genes_Wig/Merged_Female_Controls_FTL_large_smoothed_1000bp.wig -w 1000


####################################################################################################################################
# 6
# Create WIG files for single genes from various BAM files for TES
# also create additional file smoothed with 20bp-window

#ERBB2
output/Single_Genes_Wig_TES/Giant_Plasma_Merge_ERBB2.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Giant_Plasma_Merge_ERBB2.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam -g ERBB2
output/Single_Genes_Wig_TES/Merged_Female_Controls_ERBB2.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Merged_Female_Controls_ERBB2.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -g ERBB2
output/Single_Genes_Wig_TES/Merged_Male_Controls_ERBB2.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Merged_Male_Controls_ERBB2.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -g ERBB2
output/Single_Genes_Wig_TES/ERBB2_positive_ERBB2.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/ERBB2_positive_ERBB2.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/ERBB2_trim60bp.sorted.bam -g ERBB2
output/Single_Genes_Wig_TES/Breast_Cancer_ERBB2.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Breast_Cancer_ERBB2.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -g ERBB2

#MYC
output/Single_Genes_Wig_TES/Giant_Plasma_Merge_MYC.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Giant_Plasma_Merge_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam -g MYC
output/Single_Genes_Wig_TES/Merged_Female_Controls_MYC.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Merged_Female_Controls_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -g MYC
output/Single_Genes_Wig_TES/Merged_Male_Controls_MYC.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Merged_Male_Controls_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -g MYC
output/Single_Genes_Wig_TES/MYC_positive_MYC.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/MYC_positive_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/MYC_trim60bp.sorted.bam -g MYC
output/Single_Genes_Wig_TES/Breast_Cancer_MYC.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Breast_Cancer_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -g MYC
output/Single_Genes_Wig_TES/Prostate_Cancer_MYC.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Prostate_Cancer_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam -g MYC

#CCND1
output/Single_Genes_Wig_TES/Giant_Plasma_Merge_CCND1.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Giant_Plasma_Merge_CCND1.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam -g CCND1
output/Single_Genes_Wig_TES/Merged_Female_Controls_CCND1.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Merged_Female_Controls_CCND1.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -g CCND1
output/Single_Genes_Wig_TES/Merged_Male_Controls_CCND1.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Merged_Male_Controls_CCND1.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -g CCND1
output/Single_Genes_Wig_TES/CCND1_positive_CCND1.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/CCND1_positive_CCND1.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/CCND1_trim60bp.sorted.bam -g CCND1
output/Single_Genes_Wig_TES/Breast_Cancer_CCND1.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Breast_Cancer_CCND1.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -g CCND1

#AR
output/Single_Genes_Wig_TES/Giant_Plasma_Merge_AR.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Giant_Plasma_Merge_AR.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Giant_Plasma_trim60bp.sorted.bam -g AR
output/Single_Genes_Wig_TES/Merged_Female_Controls_AR.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Merged_Female_Controls_AR.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -g AR
output/Single_Genes_Wig_TES/Merged_Male_Controls_AR.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Merged_Male_Controls_AR.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -g AR
output/Single_Genes_Wig_TES/Breast_Cancer_AR.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Breast_Cancer_AR.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -g AR
output/Single_Genes_Wig_TES/Prostate_Cancer_AR.wig:
	./scripts/single_gene_wig_from_bam_TES.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig_TES/Prostate_Cancer_AR.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam -g AR


####################################################################################################################################
#
# Check whether TSS profile of breast and prostate specific genes is higher in breast cancer samples than controls
# 
#Breast specific genes
output/TSS_coverage/Breast_specific/MergedMale_Controls_Breast_Specific_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Breast_FPKMover5_Blood_unexpressed.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Breast_specific/MergedMale_Controls_Breast_Specific_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Breast_specific/MergedMale_Controls_Breast_Specific_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Breast_specific/MergedMale_Controls_Breast_Specific_tss.txt
output/TSS_coverage/Breast_specific/MergedFemale_Controls_Breast_Specific_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Breast_FPKMover5_Blood_unexpressed.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Breast_specific/MergedFemale_Controls_Breast_Specific_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Breast_specific/MergedFemale_Controls_Breast_Specific_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Breast_specific/MergedFemale_Controls_Breast_Specific_tss.txt
output/TSS_coverage/Breast_specific/Breast_Cancer_Breast_Specific_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Breast_FPKMover5_Blood_unexpressed.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Breast_specific/Breast_Cancer_Breast_Specific_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Breast_specific/Breast_Cancer_Breast_Specific_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Breast_specific/Breast_Cancer_Breast_Specific_tss.txt
output/TSS_coverage/Breast_specific/Prostate_Cancer_Breast_Specific_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Breast_FPKMover5_Blood_unexpressed.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Breast_specific/Prostate_Cancer_Breast_Specific_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Breast_specific/Prostate_Cancer_Breast_Specific_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Breast_specific/Prostate_Cancer_Breast_Specific_tss.txt

#Prostate specific genes
output/TSS_coverage/Prostate_specific/MergedMale_Controls_Prostate_Specific_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Prostate_FPKMover5_Blood_unexpressed.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Prostate_specific/MergedMale_Controls_Prostate_Specific_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Prostate_specific/MergedMale_Controls_Prostate_Specific_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Prostate_specific/MergedMale_Controls_Prostate_Specific_tss.txt
output/TSS_coverage/Prostate_specific/MergedFemale_Controls_Prostate_Specific_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Prostate_FPKMover5_Blood_unexpressed.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Prostate_specific/MergedFemale_Controls_Prostate_Specific_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Prostate_specific/MergedFemale_Controls_Prostate_Specific_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Prostate_specific/MergedFemale_Controls_Prostate_Specific_tss.txt
output/TSS_coverage/Prostate_specific/Breast_Cancer_Prostate_Specific_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Prostate_FPKMover5_Blood_unexpressed.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/BreastCancer_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Prostate_specific/Breast_Cancer_Prostate_Specific_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Prostate_specific/Breast_Cancer_Prostate_Specific_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Prostate_specific/Breast_Cancer_Prostate_Specific_tss.txt
output/TSS_coverage/Prostate_specific/Prostate_Cancer_Prostate_Specific_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/GTEx/Prostate_FPKMover5_Blood_unexpressed.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam -t 10 > output/TSS_coverage/Prostate_specific/Prostate_Cancer_Prostate_Specific_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Prostate_specific/Prostate_Cancer_Prostate_Specific_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Prostate_specific/Prostate_Cancer_Prostate_Specific_tss.txt

####################################################################################################################################
# Check coverage around TSS between expressed and unexpressed genes
output/PredictActiveGenes/MergedMale_Expression.txt:
	./scripts/analyze_active_genes_by_TSS_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam \
    -t 4 -egl ./ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -uegl ./ref/GTEx/Genes_not_expressed_in_blood_but_other_tissue.txt  \
    > output/PredictActiveGenes/MergedMale_Expression.txt
output/PredictActiveGenes/MergedFemale_Expression.txt:
	./scripts/analyze_active_genes_by_TSS_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Female_Controls_trim60bp.sorted.bam \
    -t 4 -egl ./ref/Expression_atlas_GTEx_whole_blood_FPKM5_names_only.csv -uegl ./ref/GTEx/Genes_not_expressed_in_blood_but_other_tissue.txt  \
    > output/PredictActiveGenes/MergedFemale_Expression.txt
####################################################################################################################################
#
# Check paired end plasma-samples for
#   -) coverage at ordered arrays
#   -) insert sizes
#   -) amount of 147bp fragments

#merge samples and align
output/alignments/Merged_Paired_rmdup.bam:
	for i in ./data/Paired_End_Plasma/*R1_001.fastq.gz;do base=`basename $$i`;echo "working on $$i" ; \
	fastqgz2=`echo $$i | sed s/R1/R2/`;sample=`echo $$base | sed s/R1_001\.fastq\.gz//`;bwa mem -t 4 ~/RefSeq/hg19_070510/hg19 $$i $$fastqgz2> intermediate/Paired_End/$$sample.sam; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/Paired_End/$$sample.sam OUTPUT=intermediate/Paired_End/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/Paired_End/$$sample.sam; \
	samtools rmdup -s intermediate/Paired_End/$$sample.bam intermediate/Paired_End/$$sample.rmdup.bam;rm intermediate/Paired_End/$$sample.bam;rm intermediate/Paired_End/$$sample.bai;samtools index intermediate/Paired_End/$$sample.rmdup.bam; done
	samtools merge -f output/alignments/Merged_Paired_rmdup.bam intermediate/Paired_End/*.rmdup.bam; samtools index output/alignments/Merged_Paired_rmdup.bam

#make insert size graph
output/alignments/Merged_Paired.insertsize.pdf:
	~/Downloads/jre1.7.0/bin/java -jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar CollectInsertSizeMetrics I=output/alignments/Merged_Paired.bam H=output/alignments/Merged_Paired.insertsize.pdf O=output/alignments/Merged_Paired.stats


####################################################################################################################################
#
# Try to infer genes having 5' nucleosome depleted region and calculate TSS profile of these genes
ref/refSeq_TSS_with_5prime_NDR.txt: ref/wgEncodeSydhNsomeGm12878Sig.bigWig scripts/identify_genes_with_5primer_NDR.py
	./scripts/identify_genes_with_5prime_NDR.py

output/TSS_coverage/5prime_NDR/Merged_Male_Controls_tss_coverage_5prime_NDR_genes.txt.png: output/alignments/Merged_Male_Controls_rmdup.bam ./scripts/analyze_TSS_read_depth_of_gene_list.py scripts/create_TSS_plot.R
	scripts/analyze_TSS_read_depth_of_gene_list.py -rg ref/refSeq_extended_names_strand.bed -gl ref/refSeq_TSS_with_5prime_NDR.txt \
    -b output/alignments/Merged_Male_Controls_rmdup.bam > output/TSS_coverage/5prime_NDR/GMerged_Male_Controls_tss_coverage_5prime_NDR_genes.txt 2> output/TSS_coverage/5prime_NDR/Merged_Male_Controls_tss_coverage_5prime_NDR_genes.txt.log
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/5prime_NDR/Merged_Male_Controls_tss_coverage_5prime_NDR_genes.txt

output/TSS_coverage/5prime_NDR/GiantPlasmaMerge_tss_coverage_5prime_NDR_genes.txt.png: output/alignments/Merged_Giant_Plasma_rmdup.bam ./scripts/analyze_TSS_read_depth_of_gene_list.py scripts/create_TSS_plot.R
	scripts/analyze_TSS_read_depth_of_gene_list.py -rg ref/refSeq_extended_names_strand.bed -gl ref/refSeq_TSS_with_5prime_NDR.txt \
    -b output/alignments/Merged_Giant_Plasma_rmdup.bam > output/TSS_coverage/5prime_NDR/GiantPlasmaMerge_tss_coverage_5prime_NDR_genes.txt 2> output/TSS_coverage/5prime_NDR/GiantPlasmaMerge_tss_coverage_5prime_NDR_genes.txt.log
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/5prime_NDR/GiantPlasmaMerge_tss_coverage_5prime_NDR_genes.txt

####################################################################################################################################
#
# Analyze Nucleosome tracks from ENCODE MNAse Seq Experiments 
#
# Download ENCODE Nucleosome tracks
ref/wgEncodeSydhNsome/wgEncodeSydhNsomeGm12878Sig.bigWig:
	wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhNsome/wgEncodeSydhNsomeGm12878Sig.bigWig
	mv wgEncodeSydhNsomeGm12878Sig.bigWig ./ref/wgEncodeSydhNsomeGm12878Sig.bigWig

ref/wgEncodeSydhNsome/wgEncodeSydhNsomeK562Sig.bigWig:
	wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhNsome/wgEncodeSydhNsomeK562Sig.bigWig
	mv wgEncodeSydhNsomeK562Sig.bigWig ./ref/wgEncodeSydhNsomeK562Sig.bigWig

output/TSS_coverage/10000random/ENCODE_GM12878_tss_coverage.txt.png: ref/wgEncodeSydhNsomeGm12878Sig.bigWig scripts/analyze_TSS_read_depth_fromBigWig.py
	scripts/analyze_TSS_read_depth_fromBigWig.py -rg ref/refseq_genes.bed -bw ref/wgEncodeSydhNsomeGm12878Sig.bigWig > output/TSS_coverage/10000random/ENCODE_GM12878_tss_coverage.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/10000random/ENCODE_GM12878_tss_coverage.txt

output/TSS_coverage/10000random/ENCODE_K562_tss_coverage.txt.png: ref/wgEncodeSydhNsomeK562Sig.bigWig scripts/analyze_TSS_read_depth_fromBigWig.py
	scripts/analyze_TSS_read_depth_fromBigWig.py -rg ref/refseq_genes.bed -bw ref/wgEncodeSydhNsomeK562Sig.bigWig > output/TSS_coverage/10000random/ENCODE_K562_tss_coverage.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/10000random/ENCODE_K562_tss_coverage.txt

output/TSS_coverage/10000bp/ENCODE_GM12878_tss_coverage_10000bp.txt.png: ref/wgEncodeSydhNsomeGm12878Sig.bigWig scripts/analyze_TSS_read_depth_fromBigWig.py
	scripts/analyze_TSS_read_depth_fromBigWig.py -s 10000 -e 10000 -rg ref/refseq_genes.bed -bw ref/wgEncodeSydhNsomeGm12878Sig.bigWig > output/TSS_coverage/10000bp/ENCODE_GM12878_tss_coverage_10000bp.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/10000bp/ENCODE_GM12878_tss_coverage_10000bp.txt

output/TSS_coverage/10000bp/ENCODE_K562_tss_coverage_10000bp.txt.png: ref/wgEncodeSydhNsomeK562Sig.bigWig scripts/analyze_TSS_read_depth_fromBigWig.py
	scripts/analyze_TSS_read_depth_fromBigWig.py -s 10000 -e 10000  -rg ref/refseq_genes.bed -bw ref/wgEncodeSydhNsomeK562Sig.bigWig > output/TSS_coverage/10000bp/ENCODE_K562_tss_coverage_10000bp.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/10000bp/ENCODE_K562_tss_coverage_10000bp.txt

output/TES_coverage/10000random/ENCODE_GM12878_tes_coverage.txt.png: ref/wgEncodeSydhNsomeGm12878Sig.bigWig scripts/analyze_TES_read_depth_fromBigWig.py
	scripts/analyze_TES_read_depth_fromBigWig.py -rg ref/refseq_genes.bed -bw ref/wgEncodeSydhNsomeGm12878Sig.bigWig > output/TES_coverage/10000random/ENCODE_GM12878_tes_coverage.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TES_coverage/10000random/ENCODE_GM12878_tes_coverage.txt

output/TES_coverage/10000random/ENCODE_K562_tes_coverage.txt.png: ref/wgEncodeSydhNsomeK562Sig.bigWig scripts/analyze_TES_read_depth_fromBigWig.py
	scripts/analyze_TES_read_depth_fromBigWig.py -rg ref/refseq_genes.bed -bw ref/wgEncodeSydhNsomeK562Sig.bigWig > output/TES_coverage/10000random/ENCODE_K562_tes_coverage.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TES_coverage/10000random/ENCODE_K562_tes_coverage.txt

####################################################################################################################################
#
# Analyze 8q overrepresented genes in Breast Cancers having 8q overrepresentation 
#
#
output/TSS_coverage/8q_gain/Breast_8q_gain_8q_genes.txt:
	./scripts/analyze_TSS_read_depth_of_gene_list.py -rg ref/refSeq_extended_names_strand.bed -b output/alignments/Breast_8q_gain_rmdup.bam -gl ref/refGenes_8q.txt > output/TSS_coverage/8q_gain/Breast_8q_gain_8q_genes.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/8q_gain/Breast_8q_gain_8q_genes.txt

output/TSS_coverage/8q_gain/Merged_Male_8q_genes.txt:
	./scripts/analyze_TSS_read_depth_of_gene_list.py -rg ref/refSeq_extended_names_strand.bed -b output/alignments/Merged_Male_Controls_rmdup.bam -gl ref/refGenes_8q.txt > output/TSS_coverage/8q_gain/Merged_Male_8q_genes.txt

####################################################################################################################################
#
# Analyze Top 20 expressed genes in wholeblood
#
output/Top20/Top100_tss.txt:
	./scripts/analyze_TSS_read_depth_of_gene_list.py -rg ref/refSeq_extended_names_strand.bed -b output/alignments/Merged_Male_Controls_rmdup.bam -gl ref/GTEx/Top100Blood_noChrM.csv -p > output/Top20/Top100_tss.txt

output/Top20/Top20_tss.txt:
	./scripts/analyze_TSS_read_depth_of_gene_list.py -rg ref/refSeq_extended_names_strand.bed -b output/alignments/Merged_Male_Controls_rmdup.bam -gl ref/GTEx/Top20Blood_noChrM.csv -p > output/Top20/Top20_tss.txt

output/Top20/Bottom5000_tss.txt:
	./scripts/analyze_TSS_read_depth_of_gene_list.py -rg ref/refSeq_extended_names_strand.bed -b output/alignments/Merged_Male_Controls_rmdup.bam -gl ref/GTEx/NotExpressedInBloodButOtherTissueFPKMover5.txt > output/Top20/Bottom5000_tss.txt

output/Top20/Top20_tes.txt:
	./scripts/analyze_TES_read_depth_of_gene_list.py -rg ./ref/refSeq_extended_names_strand.bed -b output/alignments/Merged_Male_Controls_rmdup.bam -gl ./ref/GTEx/Top20Blood_noChrM.csv > output/Top20/Top20_tes.txt

