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
#  -) (trim reads to 60bp)
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

output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam:
	for i in ./data/MaleControls/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`; \
    zcat $$i | ./scripts/fastx_trimmer -Q33 -f 53 -l 113 -z -o ./intermediate/trimmed/MergedMale/$${base}.trimmed.fastq.gz; \
    bwa mem -t 10 ~/RefSeq/hg19_070510/hg19 ./intermediate/trimmed/MergedMale/$${base}.trimmed.fastq.gz > intermediate/trimmed/MergedMale/$$sample.sam; \
    rm ./intermediate/trimmed/MergedMale/$${base}.trimmed.fastq.gz; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/trimmed/MergedMale/$$sample.sam OUTPUT=intermediate/trimmed/MergedMale/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/trimmed/MergedMale/$$sample.sam; \
	samtools rmdup -s intermediate/trimmed/MergedMale/$$sample.bam intermediate/trimmed/MergedMale/$$sample.rmdup.bam;rm intermediate/trimmed/MergedMale/$$sample.bam;rm intermediate/trimmed/MergedMale/$$sample.bai;samtools index intermediate/trimmed/MergedMale/$$sample.rmdup.bam; done
	samtools merge -f output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam intermediate/trimmed/MergedMale/*.rmdup.bam; samtools index output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam

#Pool of female controls
output/alignments/Merged_Female_Controls_rmdup.bam:
	for i in ./data/FemaleControls/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`;bwa mem -t 10 ~/RefSeq/hg19_070510/hg19 $$i > intermediate/FemaleControls/$$sample.sam; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/FemaleControls/$$sample.sam OUTPUT=intermediate/FemaleControls/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/FemaleControls/$$sample.sam; \
	samtools rmdup -s intermediate/FemaleControls/$$sample.bam intermediate/FemaleControls/$$sample.rmdup.bam;rm intermediate/FemaleControls/$$sample.bam;rm intermediate/FemaleControls/$$sample.bai;samtools index intermediate/FemaleControls/$$sample.rmdup.bam; done
	samtools merge -f output/alignments/Merged_Female_Controls_rmdup.bam intermediate/FemaleControls/*.rmdup.bam; samtools index output/alignments/Merged_Female_Controls_rmdup.bam

output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam:
	for i in ./data/FemaleControls/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`; \
    zcat $$i | ./scripts/fastx_trimmer -Q33 -f 53 -l 113 -m 113 -z -o ./intermediate/trimmed/MergedFemale/$${base}.trimmed.fastq.gz; \
    bwa mem -t 10 ~/RefSeq/hg19_070510/hg19 ./intermediate/trimmed/MergedFemale/$${base}.trimmed.fastq.gz > intermediate/trimmed/MergedFemale/$$sample.sam; \
    rm ./intermediate/trimmed/MergedFemale/$${base}.trimmed.fastq.gz; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/trimmed/MergedFemale/$$sample.sam OUTPUT=intermediate/trimmed/MergedFemale/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/trimmed/MergedFemale/$$sample.sam; \
	samtools rmdup -s intermediate/trimmed/MergedFemale/$$sample.bam intermediate/trimmed/MergedFemale/$$sample.rmdup.bam;rm intermediate/trimmed/MergedFemale/$$sample.bam;rm intermediate/trimmed/MergedFemale/$$sample.bai;samtools index intermediate/trimmed/MergedFemale/$$sample.rmdup.bam; done
	samtools merge -f output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam intermediate/trimmed/MergedFemale/*.rmdup.bam; samtools index output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam

output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam:
	samtools merge -f output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam
	samtools index output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam

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

output/trimmed_reads/BreastCancer_rmdup_trimmed.bam:
	for i in ./data/Breast/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`; \
    zcat $$i | ./scripts/fastx_trimmer -Q33 -f 53 -l 113 -m 113 -z -o ./intermediate/trimmed/Breast/$${base}.trimmed.fastq.gz; \
    bwa mem -t 10 ~/RefSeq/hg19_070510/hg19 ./intermediate/trimmed/Breast/$${base}.trimmed.fastq.gz > intermediate/trimmed/Breast/$$sample.sam; \
    rm ./intermediate/trimmed/Breast/$${base}.trimmed.fastq.gz; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/trimmed/Breast/$$sample.sam OUTPUT=intermediate/trimmed/Breast/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/trimmed/Breast/$$sample.sam; \
	samtools rmdup -s intermediate/trimmed/Breast/$$sample.bam intermediate/trimmed/Breast/$$sample.rmdup.bam;rm intermediate/trimmed/Breast/$$sample.bam;rm intermediate/trimmed/Breast/$$sample.bai;samtools index intermediate/trimmed/Breast/$$sample.rmdup.bam; done
	samtools merge -f output/trimmed_reads/BreastCancer_rmdup_trimmed.bam intermediate/trimmed/Breast/*.rmdup.bam; samtools index output/trimmed_reads/BreastCancer_rmdup_trimmed.bam

#Prostate cancer samples
output/alignments/ProstateCancer_rmdup.bam:
	for i in ./data/Prostate/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`;bwa mem -t 10 ~/RefSeq/hg19_070510/hg19 $$i > intermediate/Prostate/$$sample.sam; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/Prostate/$$sample.sam OUTPUT=intermediate/Prostate/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/Prostate/$$sample.sam; \
	samtools rmdup -s intermediate/Prostate/$$sample.bam intermediate/Prostate/$$sample.rmdup.bam;rm intermediate/Prostate/$$sample.bam;rm intermediate/Prostate/$$sample.bai;samtools index intermediate/Prostate/$$sample.rmdup.bam; done
	samtools merge -f output/alignments/ProstateCancer_rmdup.bam intermediate/Prostate/*.rmdup.bam; samtools index output/alignments/ProstateCancer_rmdup.bam

output/trimmed_reads/ProstateCancer_rmdup_trimmed.bam:
	for i in ./data/Prostate/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`; \
    zcat $$i | ./scripts/fastx_trimmer -Q33 -f 53 -l 113 -m 113 -z -o ./intermediate/trimmed/Prostate/$${base}.trimmed.fastq.gz; \
    bwa mem -t 10 ~/RefSeq/hg19_070510/hg19 ./intermediate/trimmed/Breast/$${base}.trimmed.fastq.gz > intermediate/trimmed/Prostate/$$sample.sam; \
    rm ./intermediate/trimmed/Prostate/$${base}.trimmed.fastq.gz; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/trimmed/Prostate/$$sample.sam OUTPUT=intermediate/trimmed/Prostate/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/trimmed/Prostate/$$sample.sam; \
	samtools rmdup -s intermediate/trimmed/Prostate/$$sample.bam intermediate/trimmed/Prostate/$$sample.rmdup.bam;rm intermediate/trimmed/Prostate/$$sample.bam;rm intermediate/trimmed/Prostate/$$sample.bai;samtools index intermediate/trimmed/Prostate/$$sample.rmdup.bam; done
	samtools merge -f output/trimmed_reads/ProstateCancer_rmdup_trimmed.bam intermediate/trimmed/Prostate/*.rmdup.bam; samtools index output/trimmed_reads/ProstateCancer_rmdup_trimmed.bam

#CNV samples (low-coverage WGS high-molecular-weight DNA)
output/alignments/CNV_samples_rmdup.bam:
	for i in ./data/CNV_samples/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`;bwa mem -t 6 ~/RefSeq/hg19_070510/hg19 $$i > intermediate/CNV_samples/$$sample.sam; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/CNV_samples/$$sample.sam OUTPUT=intermediate/CNV_samples/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/CNV_samples/$$sample.sam; \
	samtools rmdup -s intermediate/CNV_samples/$$sample.bam intermediate/CNV_samples/$$sample.rmdup.bam;rm intermediate/CNV_samples/$$sample.bam;rm intermediate/CNV_samples/$$sample.bai;samtools index intermediate/CNV_samples/$$sample.rmdup.bam; done
	samtools merge -f output/alignments/CNV_samples_rmdup.bam intermediate/CNV_samples/*.rmdup.bam; samtools index output/alignments/CNV_samples_rmdup.bam

output/trimmed_reads/CNV_samples_rmdup_trimmed.bam:
	for i in ./data/CNV_samples/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`; \
    zcat $$i | ./scripts/fastx_trimmer -Q33 -f 53 -l 113 -m 113 -z -o ./intermediate/trimmed/CNVsamples/$${base}.trimmed.fastq.gz; \
    bwa mem -t 10 ~/RefSeq/hg19_070510/hg19 ./intermediate/trimmed/CNVsamples/$${base}.trimmed.fastq.gz > intermediate/trimmed/CNVsamples/$$sample.sam; \
    rm ./intermediate/trimmed/CNVsamples/$${base}.trimmed.fastq.gz; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/trimmed/CNVsamples/$$sample.sam OUTPUT=intermediate/trimmed/CNVsamples/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/trimmed/CNVsamples/$$sample.sam; \
	samtools rmdup -s intermediate/trimmed/CNVsamples/$$sample.bam intermediate/trimmed/CNVsamples/$$sample.rmdup.bam;rm intermediate/trimmed/CNVsamples/$$sample.bam;rm intermediate/trimmed/CNVsamples/$$sample.bai;samtools index intermediate/trimmed/CNVsamples/$$sample.rmdup.bam; done
	samtools merge -f output/trimmed_reads/CNV_samples_rmdup_trimmed.bam intermediate/trimmed/CNVsamples/*.rmdup.bam; samtools index output/trimmed_reads/CNV_samples_rmdup_trimmed.bam

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

output/trimmed_reads/Biphasic_Breast_rmdup_trimmed.bam:
	for i in ./data/Biphasic_Breast/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`; \
    zcat $$i | ./scripts/fastx_trimmer -Q33 -f 53 -l 113 -m 113 -z -o ./intermediate/trimmed/BiphasicBreast/$${base}.trimmed.fastq.gz; \
    bwa mem -t 10 ~/RefSeq/hg19_070510/hg19 ./intermediate/trimmed/BiphasicBreast/$${base}.trimmed.fastq.gz > intermediate/trimmed/BiphasicBreast/$$sample.sam; \
    rm ./intermediate/trimmed/BiphasicBreast/$${base}.trimmed.fastq.gz; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/trimmed/BiphasicBreast/$$sample.sam OUTPUT=intermediate/trimmed/BiphasicBreast/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/trimmed/BiphasicBreast/$$sample.sam; \
	samtools rmdup -s intermediate/trimmed/BiphasicBreast/$$sample.bam intermediate/trimmed/BiphasicBreast/$$sample.rmdup.bam;rm intermediate/trimmed/BiphasicBreast/$$sample.bam;rm intermediate/trimmed/BiphasicBreast/$$sample.bai;samtools index intermediate/trimmed/BiphasicBreast/$$sample.rmdup.bam; done
	samtools merge -f output/trimmed_reads/Biphasic_Breast_rmdup_trimmed.bam intermediate/trimmed/BiphasicBreast/*.rmdup.bam; samtools index output/trimmed_reads/Biphasic_Breast_rmdup_trimmed.bam

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

output/trimmed_reads/ERBB2_rmdup_trimmed.bam:
	for i in ./data/FocalAmps/ERBB2/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`; \
    zcat $$i | ./scripts/fastx_trimmer -Q33 -f 53 -l 113 -m 113 -z -o ./intermediate/trimmed/ERBB2/$${base}.trimmed.fastq.gz; \
    bwa mem -t 10 ~/RefSeq/hg19_070510/hg19 ./intermediate/trimmed/ERBB2/$${base}.trimmed.fastq.gz > intermediate/trimmed/ERBB2/$$sample.sam; \
    rm ./intermediate/trimmed/ERBB2/$${base}.trimmed.fastq.gz; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/trimmed/ERBB2/$$sample.sam OUTPUT=intermediate/trimmed/ERBB2/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/trimmed/ERBB2/$$sample.sam; \
	samtools rmdup -s intermediate/trimmed/ERBB2/$$sample.bam intermediate/trimmed/ERBB2/$$sample.rmdup.bam;rm intermediate/trimmed/ERBB2/$$sample.bam;rm intermediate/trimmed/ERBB2/$$sample.bai;samtools index intermediate/trimmed/ERBB2/$$sample.rmdup.bam; done
	samtools merge -f output/trimmed_reads/ERBB2_rmdup_trimmed.bam intermediate/trimmed/ERBB2/*.rmdup.bam; samtools index output/trimmed_reads/ERBB2_rmdup_trimmed.bam

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

output/trimmed_reads/CCND1_rmdup_trimmed.bam:
	for i in ./data/FocalAmps/CCND1/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`; \
    zcat $$i | ./scripts/fastx_trimmer -Q33 -f 53 -l 113 -m 113 -z -o ./intermediate/trimmed/CCND1/$${base}.trimmed.fastq.gz; \
    bwa mem -t 10 ~/RefSeq/hg19_070510/hg19 ./intermediate/trimmed/CCND1/$${base}.trimmed.fastq.gz > intermediate/trimmed/CCND1/$$sample.sam; \
    rm ./intermediate/trimmed/CCND1/$${base}.trimmed.fastq.gz; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/trimmed/CCND1/$$sample.sam OUTPUT=intermediate/trimmed/CCND1/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/trimmed/CCND1/$$sample.sam; \
	samtools rmdup -s intermediate/trimmed/CCND1/$$sample.bam intermediate/trimmed/CCND1/$$sample.rmdup.bam;rm intermediate/trimmed/CCND1/$$sample.bam;rm intermediate/trimmed/CCND1/$$sample.bai;samtools index intermediate/trimmed/CCND1/$$sample.rmdup.bam; done
	samtools merge -f output/trimmed_reads/CCND1_rmdup_trimmed.bam intermediate/trimmed/CCND1/*.rmdup.bam; samtools index output/trimmed_reads/CCND1_rmdup_trimmed.bam

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

output/trimmed_reads/MYC_rmdup_trimmed.bam:
	for i in ./data/FocalAmps/MYC/*.fastq.gz;do base=`basename $$i`;echo "working on $$i" ;sample=`echo $$base | sed s/\.fastq\.gz//`; \
    zcat $$i | ./scripts/fastx_trimmer -Q33 -f 53 -l 113 -m 113 -z -o ./intermediate/trimmed/MYC/$${base}.trimmed.fastq.gz; \
    bwa mem -t 8 ~/RefSeq/hg19_070510/hg19 ./intermediate/trimmed/MYC/$${base}.trimmed.fastq.gz > intermediate/trimmed/MYC/$$sample.sam; \
    rm ./intermediate/trimmed/MYC/$${base}.trimmed.fastq.gz; \
	~/Downloads/jre1.7.0/bin/java -Xmx4g -Djava.io.tmpdir=/tmp \
	-jar ~/Software/SNP_calling/picard-tools-1.128/picard.jar SortSam \
	SO=coordinate INPUT=intermediate/trimmed/MYC/$$sample.sam OUTPUT=intermediate/trimmed/MYC/$$sample.bam \
	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true; rm intermediate/trimmed/MYC/$$sample.sam; \
	samtools rmdup -s intermediate/trimmed/MYC/$$sample.bam intermediate/trimmed/MYC/$$sample.rmdup.bam;rm intermediate/trimmed/MYC/$$sample.bam;rm intermediate/trimmed/MYC/$$sample.bai;samtools index intermediate/trimmed/MYC/$$sample.rmdup.bam; done
	samtools merge -f output/trimmed_reads/MYC_rmdup_trimmed.bam intermediate/trimmed/MYC/*.rmdup.bam; samtools index output/trimmed_reads/MYC_rmdup_trimmed.bam

####################################################################################################################################
# 3
# Analyze mean coverage around Transcription start sites (TSS)
#   -) All genes
#   -) Apoptosis Genes
#   -) Plasma RNASeq Top1000 (only Genes with NM accession in RefSeq)
#   -) Plasma RNASeq Bottom1000 (only Genes with NM accession in RefSeq)
#   -) Plasma RNASeq Top20 (only Genes with NM accession in RefSeq)
#   -) Plasma RNASeq Top50 (only Genes with NM accession in RefSeq)
#   -) Plasma RNASeq Top100 (only Genes with NM accession in RefSeq)
#   -) Plasma RNASeq Top500 (only Genes with NM accession in RefSeq)
#   

#3.1 All genes
output/TSS_coverage/All/Merged_Male_Controls_all_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/All/Merged_Male_Controls_all_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/All/Merged_Male_Controls_all_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/All/Merged_Male_Controls_all_tss.txt
output/TSS_coverage/All/Merged_Female_Controls_all_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/All/Merged_Female_Controls_all_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/All/Merged_Female_Controls_all_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/All/Merged_Female_Controls_all_tss.txt
output/TSS_coverage/All/Merged_Controls_all_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/All/Merged_Controls_all_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/All/Merged_Controls_all_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/All/Merged_Controls_all_tss.txt

#3.2 Apoptosis genes
output/TSS_coverage/Apoptosis_genes/Merged_Male_Controls_Apoptosis_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Apoptosis_genes.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/Apoptosis_genes/Merged_Male_Controls_Apoptosis_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Apoptosis_genes/Merged_Male_Controls_Apoptosis_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Apoptosis_genes/Merged_Male_Controls_Apoptosis_tss.txt
output/TSS_coverage/Apoptosis_genes/Merged_Female_Controls_Apoptosis_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Apoptosis_genes.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/Apoptosis_genes/Merged_Female_Controls_Apoptosis_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Apoptosis_genes/Merged_Female_Controls_Apoptosis_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Apoptosis_genes/Merged_Female_Controls_Apoptosis_tss.txt
output/TSS_coverage/Apoptosis_genes/Merged_Controls_Apoptosis_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Apoptosis_genes.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/Apoptosis_genes/Merged_Controls_Apoptosis_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/Apoptosis_genes/Merged_Controls_Apoptosis_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/Apoptosis_genes/Merged_Controls_Apoptosis_tss.txt

#3.3 Plasma RNASeq Top1000
output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Top1000_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top1000_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Top1000_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Top1000_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Top1000_tss.txt
output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Top1000_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top1000_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Top1000_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Top1000_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Top1000_tss.txt
output/TSS_coverage/PlasmaRNASeq/BreastCancer_Top1000_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top1000_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/BreastCancer_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/BreastCancer_Top1000_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/BreastCancer_Top1000_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/BreastCancer_Top1000_tss.txt
output/TSS_coverage/PlasmaRNASeq/ProstateCancer_Top1000_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top1000_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/ProstateCancer_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/ProstateCancer_Top1000_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/ProstateCancer_Top1000_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/ProstateCancer_Top1000_tss.txt
output/TSS_coverage/PlasmaRNASeq/CNV_samples_Top1000_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top1000_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/CNV_samples_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/CNV_samples_Top1000_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/CNV_samples_Top1000_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/CNV_samples_Top1000_tss.txt
output/TSS_coverage/PlasmaRNASeq/Biphasic_Breast_Top1000_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top1000_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/BreastCancer_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/Biphasic_Breast_Top1000_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Biphasic_Breast_Top1000_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Biphasic_Breast_Top1000_tss.txt
output/TSS_coverage/PlasmaRNASeq/ERBB2_Top1000_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top1000.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/ERBB2_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/ERBB2_Top1000_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/ERBB2_Top1000_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/ERBB2_Top1000_tss.txt
output/TSS_coverage/PlasmaRNASeq/CCND1_Top1000_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top1000_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/CCND1_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/CCND1_Top1000_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/CCND1_Top1000_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/CCND1_Top1000_tss.txt
output/TSS_coverage/PlasmaRNASeq/MYC_Top1000_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top1000_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/MYC_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/MYC_Top1000_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/MYC_Top1000_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/MYC_Top1000_tss.txt
output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top1000_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top1000_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top1000_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top1000_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top1000_tss.txt

#3.4 Bottom1000
output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Bottom1000_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Bottom1000_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Bottom1000_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Bottom1000_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Bottom1000_tss.txt
output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Bottom1000_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Bottom1000_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Bottom1000_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Bottom1000_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Bottom1000_tss.txt
output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Bottom1000_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Bottom1000_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Bottom1000_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Bottom1000_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Bottom1000_tss.txt

#Top20
output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top20_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top20_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top20_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top20_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top20_tss.txt
output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Top20_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top20_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Top20_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Top20_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Top20_tss.txt
output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Top20_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top20_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Top20_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Top20_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Top20_tss.txt

#Top20 but without any Hemoglobin gene
output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top20_tss_noHemoglobin.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top20_NMonly_noHemoglobin.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top20_tss_noHemoglobin.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top20_tss_noHemoglobin.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top20_tss_noHemoglobin.txt

#Top50
output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top50_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top50_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top50_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top50_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top50_tss.txt
output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Top50_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top50_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Top50_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Top50_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Top50_tss.txt
output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Top50_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top50_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Top50_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Top50_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Top50_tss.txt
#Top100
output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top100_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top100_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top100_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top100_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top100_tss.txt
output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Top100_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top100_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Top100_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Top100_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Top100_tss.txt
output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Top100_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top100_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Top100_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Top100_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Top100_tss.txt
#Top500
output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top500_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top500_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top500_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top500_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/MergedControls_Plasma_Top500_tss.txt
output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Top500_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top500_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Top500_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Top500_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Male_Controls_Plasma_Top500_tss.txt
output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Top500_tss.txt: ./scripts/analyze_TSS_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TSS_coverage.py -gl ref/Plasma-RNASeq/Top500_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam -t 10 > output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Top500_tss.txt
	cat scripts/create_TSS_plot_extended.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Top500_tss.txt
	cat scripts/create_TSS_plot.R | R --slave --args output/TSS_coverage/PlasmaRNASeq/Merged_Female_Controls_Plasma_Top500_tss.txt

####################################################################################################################################
# 4
# Analyze mean coverage around Transcription end sites (TES)
#   -) Plasma RNASeq Top1000
#

#4.1 Plasma RNASeq Top1000 and Bottom1000 only entries with NM number
output/TES_coverage/PlasmaRNASeq/MergedControls_Plasma_Bottom1000_NMonly_tes.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/Plasma-RNASeq/Bottom1000_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam -t 10 > output/TES_coverage/PlasmaRNASeq/MergedControls_Plasma_Bottom1000_NMonly_tes.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/PlasmaRNASeq/MergedControls_Plasma_Bottom1000_NMonly_tes.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/PlasmaRNASeq/MergedControls_Plasma_Bottom1000_NMonly_tes.txt
output/TES_coverage/PlasmaRNASeq/MergedControls_Plasma_Top1000_NMonly_tes.txt: ./scripts/analyze_TES_coverage.py ref/refSeq_extended_names_strand.bed
	./scripts/analyze_TES_coverage.py -gl ref/Plasma-RNASeq/Top1000_NMonly.txt -rg ref/refSeq_extended_names_strand.bed -m 0 -b output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam -t 10 > output/TES_coverage/PlasmaRNASeq/MergedControls_Plasma_Top1000_NMonly_tes.txt
	cat scripts/create_TES_plot_extended.R | R --slave --args output/TES_coverage/PlasmaRNASeq/MergedControls_Plasma_Top1000_NMonly_tes.txt
	cat scripts/create_TES_plot.R | R --slave --args output/TES_coverage/PlasmaRNASeq/MergedControls_Plasma_Top1000_NMonly_tes.txt

####################################################################################################################################
# 5
# Create WIG files for single genes from various BAM files for TSS
#   

#ERBB2
output/Single_Genes_Wig/Merged_Female_Controls_ERBB2.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Female_Controls_ERBB2.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam -g ERBB2
output/Single_Genes_Wig/Merged_Male_Controls_ERBB2.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Male_Controls_ERBB2.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam -g ERBB2
output/Single_Genes_Wig/Merged_Male_Controls_ERBB2_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Male_Controls_ERBB2_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam -g ERBB2
output/Single_Genes_Wig/Merged_Male_Controls_ERBB2_large_smoothed_100bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Male_Controls_ERBB2_large.wig -o output/Single_Genes_Wig/Merged_Male_Controls_ERBB2_large_smoothed_100bp.wig -w 100
output/Single_Genes_Wig/Merged_Controls_ERBB2_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Controls_ERBB2_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam -g ERBB2
output/Single_Genes_Wig/Merged_Controls_ERBB2_large_smoothed_100bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Controls_ERBB2_large.wig -o output/Single_Genes_Wig/Merged_Controls_ERBB2_large_smoothed_100bp.wig -w 100

output/Single_Genes_Wig/ERBB2_positive_ERBB2.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/ERBB2_positive_ERBB2.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/ERBB2_rmdup_trimmed.bam -g ERBB2
output/Single_Genes_Wig/ERBB2_positive_ERBB2_smoothed_10window.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/ERBB2_positive_ERBB2.wig -o output/Single_Genes_Wig/ERBB2_positive_ERBB2_smoothed_10window.wig
output/Single_Genes_Wig/ERBB2_positive_ERBB2_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/ERBB2_positive_ERBB2_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/ERBB2_rmdup_trimmed.bam -g ERBB2
output/Single_Genes_Wig/ERBB2_positive_ERBB2_large_smoothed_100bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/ERBB2_positive_ERBB2_large.wig -o output/Single_Genes_Wig/ERBB2_positive_ERBB2_large_smoothed_100bp.wig -w 100

output/Single_Genes_Wig/Breast_Cancer_ERBB2.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Breast_Cancer_ERBB2.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/BreastCancer_rmdup_trimmed.bam -g ERBB2

#MYC
output/Single_Genes_Wig/Merged_Female_Controls_MYC.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Female_Controls_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam -g MYC
output/Single_Genes_Wig/Merged_Male_Controls_MYC.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Male_Controls_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam -g MYC
output/Single_Genes_Wig/MYC_positive_MYC.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/MYC_positive_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/MYC_trim60bp.sorted.bam -g MYC
output/Single_Genes_Wig/MYC_positive_MYC_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/MYC_positive_MYC.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/MYC_trim60bp.sorted.bam -g MYC
output/Single_Genes_Wig/MYC_positive_MYC_large_smoothed_100bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/MYC_positive_MYC_large.wig -o output/Single_Genes_Wig/MYC_positive_MYC_large_smoothed_100bp.wig -w 100
output/Single_Genes_Wig/Merged_Controls_MYC_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Controls_MYC_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam -g MYC
output/Single_Genes_Wig/Merged_Controls_MYC_large_smoothed_100bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Controls_MYC_large.wig -o output/Single_Genes_Wig/Merged_Controls_MYC_large_smoothed_100bp.wig -w 100

output/Single_Genes_Wig/Breast_Cancer_MYC.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Breast_Cancer_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/BreastCancer_rmdup_trimmed.bam -g MYC
output/Single_Genes_Wig/Prostate_Cancer_MYC.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Prostate_Cancer_MYC.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/ProstateCancer_trim60bp.sorted.bam -g MYC

#CCND1
output/Single_Genes_Wig/Merged_Female_Controls_CCND1.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Female_Controls_CCND1.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam -g CCND1
output/Single_Genes_Wig/Merged_Male_Controls_CCND1.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Male_Controls_CCND1.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam -g CCND1
output/Single_Genes_Wig/Merged_Male_Controls_CCND1_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Male_Controls_CCND1_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/Merged_Male_Controls_trim60bp.sorted.bam -g CCND1
output/Single_Genes_Wig/Merged_Male_Controls_CCND1_large_smoothed_100bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Male_Controls_CCND1_large.wig -o output/Single_Genes_Wig/Merged_Male_Controls_CCND1_large_smoothed_100bp.wig -w 100
output/Single_Genes_Wig/Merged_Controls_CCND1_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Controls_CCND1_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam -g CCND1
output/Single_Genes_Wig/Merged_Controls_CCND1_large_smoothed_100bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Controls_CCND1_large.wig -o output/Single_Genes_Wig/Merged_Controls_CCND1_large_smoothed_100bp.wig -w 100

output/Single_Genes_Wig/CCND1_positive_CCND1.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/CCND1_positive_CCND1.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/CCND1_rmdup_trimmed.bam -g CCND1
output/Single_Genes_Wig/CCND1_positive_CCND1_smoothed_10window.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/CCND1_positive_CCND1.wig -o output/Single_Genes_Wig/CCND1_positive_CCND1_smoothed_10window.wig  
output/Single_Genes_Wig/CCND1_positive_CCND1_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/CCND1_positive_CCND1_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/CCND1_rmdup_trimmed.bam -g CCND1
output/Single_Genes_Wig/CCND1_positive_CCND1_large_smoothed_100bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/CCND1_positive_CCND1_large.wig -o output/Single_Genes_Wig/CCND1_positive_CCND1_large_smoothed_100bp.wig -w 100

output/Single_Genes_Wig/Breast_Cancer_CCND1.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Breast_Cancer_CCND1.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/BreastCancer_rmdup_trimmed.bam -g CCND1

#AR
output/Single_Genes_Wig/Merged_Female_Controls_AR.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Female_Controls_AR.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam -g AR
output/Single_Genes_Wig/Merged_Female_Controls_AR_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Female_Controls_AR_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam -g AR
output/Single_Genes_Wig/Merged_Female_Controls_AR_large_smoothed_100bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Female_Controls_AR_large.wig -o output/Single_Genes_Wig/Merged_Female_Controls_AR_large_smoothed_100bp.wig -w 100
output/Single_Genes_Wig/Merged_Male_Controls_AR.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Male_Controls_AR.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam -g AR
output/Single_Genes_Wig/Merged_Controls_AR_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Merged_Controls_AR_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam -g AR
output/Single_Genes_Wig/Merged_Controls_AR_large_smoothed_100bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Merged_Controls_AR_large.wig -o output/Single_Genes_Wig/Merged_Controls_AR_large_smoothed_100bp.wig -w 100

output/Single_Genes_Wig/Breast_Cancer_AR.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Breast_Cancer_AR.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/BreastCancer_rmdup_trimmed.bam -g AR
output/Single_Genes_Wig/Prostate_Cancer_AR.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Prostate_Cancer_AR.wig \
    -s 10000 -e 10000 -b output/trimmed_reads/ProstateCancer_rmdup_trimmed.bam -g AR
output/Single_Genes_Wig/Prostate_Cancer_AR_smoothed_10window.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Prostate_Cancer_AR.wig -o output/Single_Genes_Wig/Prostate_Cancer_AR_smoothed_10window.wig
output/Single_Genes_Wig/Prostate_Cancer_AR_large.wig:
	./scripts/single_gene_wig_from_bam.py -rg ./ref/refSeq_extended_names_strand.bed -w output/Single_Genes_Wig/Prostate_Cancer_AR_large.wig \
    -s 100000 -e 100000 -b output/trimmed_reads/ProstateCancer_rmdup_trimmed.bam -g AR
output/Single_Genes_Wig/Prostate_Cancer_AR_large_smoothed_100bp.wig:
	./scripts/smooth_wig_file.py -wig output/Single_Genes_Wig/Prostate_Cancer_AR_large.wig -o output/Single_Genes_Wig/Prostate_Cancer_AR_large_smoothed_100bp.wig -w 100

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
#  Try to predict whether single genes are expressed by analysis of Top1000 and Bottom1000 genes at four parameters:
#    1) Mean Coverage +/- 1000bp around Transcription start site
#    2) Slope of 400bp-smoothed coverage data in
#       2.1) 1000bp before TSS and TSS (5' Slope, should be negative for expressed genes)
#       2.2) TSS and 1000bp into gene (3' Slope, should be positive for expressed genes)
#    3) Coverage in smaller region around TSS (-150;+50); corrected by mean coverage of +/-1000bp around TSS

#Merged Controls
output/PredictActiveGenes/Merged_Controls/MergedControls_TSSCoverage_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top1000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom1000_NMonly.txt \
    > output/PredictActiveGenes/Merged_Controls/MergedControls_TSSCoverage_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes/Merged_Controls/MergedControls_TSSCoverage_small_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage_small.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top1000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom1000_NMonly.txt -ws 1000 -we 1000 -ss 150 -se 50 \
    > output/PredictActiveGenes/Merged_Controls/MergedControls_TSSCoverage_small_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes/Merged_Controls/MergedControls_LowPass_Slope_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_active_genes_by_low_pass_slope_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top1000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom1000_NMonly.txt -tmp ./intermediate/slope/ \
    > output/PredictActiveGenes/Merged_Controls/MergedControls_LowPass_Slope_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes/Merged_Controls/MergedControls_prediction_data.csv:
	./scripts/combine_prediction_parameters.py -o output/PredictActiveGenes/Merged_Controls/MergedControls_prediction_data.csv \
    -bcov output/PredictActiveGenes/Merged_Controls/MergedControls_TSSCoverage_PlasmaRNASeq_NMonly.txt \
    -scov output/PredictActiveGenes/Merged_Controls/MergedControls_TSSCoverage_small_PlasmaRNASeq_NMonly.txt \
    -slope output/PredictActiveGenes/Merged_Controls/MergedControls_LowPass_Slope_PlasmaRNASeq_NMonly.txt

#Male Controls 
output/PredictActiveGenes/Male_Controls/Male_Controls_TSSCoverage_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top1000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom1000_NMonly.txt \
    > output/PredictActiveGenes/Male_Controls/Male_Controls_TSSCoverage_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes/Male_Controls/Male_Controls_TSSCoverage_small_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage_small.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top1000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom1000_NMonly.txt -ws 1000 -we 1000 -ss 150 -se 50 \
    > output/PredictActiveGenes/Male_Controls/Male_Controls_TSSCoverage_small_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes/Male_Controls/Male_Controls_LowPass_Slope_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_active_genes_by_low_pass_slope_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top1000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom1000_NMonly.txt -tmp ./intermediate/slope/ \
    > output/PredictActiveGenes/Male_Controls/Male_Controls_LowPass_Slope_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes/Male_Controls/Male_Controls_prediction_data.csv:
	./scripts/combine_prediction_parameters.py -o output/PredictActiveGenes/Male_Controls/Male_Controls_prediction_data.csv \
    -bcov output/PredictActiveGenes/Male_Controls/Male_Controls_TSSCoverage_PlasmaRNASeq_NMonly.txt \
    -scov output/PredictActiveGenes/Male_Controls/Male_Controls_TSSCoverage_small_PlasmaRNASeq_NMonly.txt \
    -slope output/PredictActiveGenes/Male_Controls/Male_Controls_LowPass_Slope_PlasmaRNASeq_NMonly.txt

#Female Controls 
output/PredictActiveGenes/Female_Controls/Female_Controls_TSSCoverage_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top1000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom1000_NMonly.txt \
    > output/PredictActiveGenes/Female_Controls/Female_Controls_TSSCoverage_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes/Female_Controls/Female_Controls_TSSCoverage_small_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage_small.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top1000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom1000_NMonly.txt -ws 1000 -we 1000 -ss 150 -se 50 \
    > output/PredictActiveGenes/Female_Controls/Female_Controls_TSSCoverage_small_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes/Female_Controls/Female_Controls_LowPass_Slope_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_active_genes_by_low_pass_slope_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top1000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom1000_NMonly.txt -tmp ./intermediate/slope/ \
    > output/PredictActiveGenes/Female_Controls/Female_Controls_LowPass_Slope_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes/Female_Controls/Female_Controls_prediction_data.csv:
	./scripts/combine_prediction_parameters.py -o output/PredictActiveGenes/Female_Controls/Female_Controls_prediction_data.csv \
    -bcov output/PredictActiveGenes/Female_Controls/Female_Controls_TSSCoverage_PlasmaRNASeq_NMonly.txt \
    -scov output/PredictActiveGenes/Female_Controls/Female_Controls_TSSCoverage_small_PlasmaRNASeq_NMonly.txt \
    -slope output/PredictActiveGenes/Female_Controls/Female_Controls_LowPass_Slope_PlasmaRNASeq_NMonly.txt

#BreastCancer
output/PredictActiveGenes/Breast/Breast_Cancer_TSSCoverage_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/BreastCancer_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top1000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom1000_NMonly.txt \
    > output/PredictActiveGenes/Breast/Breast_Cancer_TSSCoverage_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes/Breast/Breast_Cancer_TSSCoverage_small_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage_small.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/BreastCancer_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top1000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom1000_NMonly.txt -ws 1000 -we 1000 -ss 150 -se 50 \
    > output/PredictActiveGenes/Breast/Breast_Cancer_TSSCoverage_small_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes/Breast/Breast_Cancer_LowPass_Slope_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_active_genes_by_low_pass_slope_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/BreastCancer_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top1000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom1000_NMonly.txt -tmp ./intermediate/slope/ \
    > output/PredictActiveGenes/Breast/Breast_Cancer_LowPass_Slope_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes/Breast/Breast_Cancer_prediction_data.csv:
	./scripts/combine_prediction_parameters.py -o output/PredictActiveGenes/Breast/Breast_Cancer_prediction_data.csv \
    -bcov output/PredictActiveGenes/Breast/Breast_Cancer_TSSCoverage_PlasmaRNASeq_NMonly.txt \
    -scov output/PredictActiveGenes/Breast/Breast_Cancer_TSSCoverage_small_PlasmaRNASeq_NMonly.txt \
    -slope output/PredictActiveGenes/Breast/Breast_Cancer_LowPass_Slope_PlasmaRNASeq_NMonly.txt

#CNV samples 
output/PredictActiveGenes/CNVsamples/CNV_samples_TSSCoverage_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/CNV_samples_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top1000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom1000_NMonly.txt \
    > output/PredictActiveGenes/CNVsamples/CNV_samples_TSSCoverage_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes/CNVsamples/CNV_samples_TSSCoverage_small_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage_small.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/CNV_samples_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top1000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom1000_NMonly.txt -ws 1000 -we 1000 -ss 150 -se 50 \
    > output/PredictActiveGenes/CNVsamples/CNV_samples_TSSCoverage_small_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes/CNVsamples/CNV_samples_LowPass_Slope_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_active_genes_by_low_pass_slope_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/CNV_samples_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top1000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom1000_NMonly.txt -tmp ./intermediate/slope/ \
    > output/PredictActiveGenes/CNVsamples/CNV_samples_LowPass_Slope_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes/CNVsamples/CNV_samples_prediction_data.csv:
	./scripts/combine_prediction_parameters.py -o output/PredictActiveGenes/CNVsamples/CNV_samples_prediction_data.csv \
    -bcov output/PredictActiveGenes/CNVsamples/CNV_samples_TSSCoverage_PlasmaRNASeq_NMonly.txt \
    -scov output/PredictActiveGenes/CNVsamples/CNV_samples_TSSCoverage_small_PlasmaRNASeq_NMonly.txt \
    -slope output/PredictActiveGenes/CNVsamples/CNV_samples_LowPass_Slope_PlasmaRNASeq_NMonly.txt

#"Dilution" series adding more and more genes and see wether there is an upward trend in coverage
output/PredictActiveGenes/Dilution/MergedControls_TSSCoverage_PlasmaRNASeq_Top20_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam \
    -t 1 -egl ./ref/Plasma-RNASeq/Top20_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom50_NM.txt \
    > output/PredictActiveGenes/Dilution/MergedControls_TSSCoverage_PlasmaRNASeq_Top20_NMonly.txt
output/PredictActiveGenes/Dilution/MergedControls_TSSCoverage_PlasmaRNASeq_Top50_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam \
    -t 1 -egl ./ref/Plasma-RNASeq/Top50_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom50_NM.txt \
    > output/PredictActiveGenes/Dilution/MergedControls_TSSCoverage_PlasmaRNASeq_Top50_NMonly.txt
output/PredictActiveGenes/Dilution/MergedControls_TSSCoverage_PlasmaRNASeq_Top100_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam \
    -t 2 -egl ./ref/Plasma-RNASeq/Top100_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom50_NM.txt \
    > output/PredictActiveGenes/Dilution/MergedControls_TSSCoverage_PlasmaRNASeq_Top100_NMonly.txt
output/PredictActiveGenes/Dilution/MergedControls_TSSCoverage_PlasmaRNASeq_Top500_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam \
    -t 3 -egl ./ref/Plasma-RNASeq/Top500_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom50_NM.txt \
    > output/PredictActiveGenes/Dilution/MergedControls_TSSCoverage_PlasmaRNASeq_Top500_NMonly.txt
output/PredictActiveGenes/Dilution/MergedControls_TSSCoverage_PlasmaRNASeq_Top1000_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top1000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom1000_NMonly.txt \
    > output/PredictActiveGenes/Dilution/MergedControls_TSSCoverage_PlasmaRNASeq_Top1000_NMonly.txt
output/PredictActiveGenes/Dilution/MergedControls_TSSCoverage_PlasmaRNASeq_Top5000_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top5000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom50_NM.txt \
    > output/PredictActiveGenes/Dilution/MergedControls_TSSCoverage_PlasmaRNASeq_Top5000_NMonly.txt
output/PredictActiveGenes/Dilution/MergedControls_TSSCoverage_PlasmaRNASeq_Top10000_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top10000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom50_NM.txt \
    > output/PredictActiveGenes/Dilution/MergedControls_TSSCoverage_PlasmaRNASeq_Top10000_NMonly.txt
output/PredictActiveGenes/Dilution/MergedControls_TSSCoverage_PlasmaRNASeq_Top15000_NMonly.txt:
	./scripts/analyze_active_genes_by_TSS_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam \
    -t 10 -egl ./ref/Plasma-RNASeq/Top15000_NMonly.txt -uegl ./ref/Plasma-RNASeq/Bottom50_NM.txt \
    > output/PredictActiveGenes/Dilution/MergedControls_TSSCoverage_PlasmaRNASeq_Top15000_NMonly.txt

####################################################################################################################################
#  Try to predict whether single genes are expressed by analysis of
#    1) Mean Coverage +/- 1000bp around Transcription start site
#    2) Slope of 400bp-smoothed coverage data in
#       2.1) 1000bp before TSS and TSS (5' Slope, should be negative for expressed genes)
#       2.2) TSS and 1000bp into gene (3' Slope, should be positive for expressed genes)
#    3) Coverage in smaller region around TSS (-150;+50); corrected by mean coverage of +/-1000bp around TSS

#Merged Controls
output/PredictActiveGenes_All/Merged_Controls/MergedControls_TSSCoverage_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_all_genes_by_TSS_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam \
    -t 10 -gl ./ref/Plasma-RNASeq/AllGenes_NMonly.txt \
    > output/PredictActiveGenes_All/Merged_Controls/MergedControls_TSSCoverage_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes_All/Merged_Controls/MergedControls_TSSCoverage_small_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_all_genes_by_TSS_coverage_small.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam \
    -t 10 -gl ./ref/Plasma-RNASeq/AllGenes_NMonly.txt -ws 1000 -we 1000 -ss 150 -se 50 \
    > output/PredictActiveGenes_All/Merged_Controls/MergedControls_TSSCoverage_small_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes_All/Merged_Controls/MergedControls_LowPass_Slope_PlasmaRNASeq_NMonly.txt:
	./scripts/analyze_all_genes_by_low_pass_slope_coverage.py -rg ./ref/refSeq_extended_names_strand.bed -b ./output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam \
    -t 10 -gl ./ref/Plasma-RNASeq/AllGenes_NMonly.txt -tmp ./intermediate/slope/ \
    > output/PredictActiveGenes_All/Merged_Controls/MergedControls_LowPass_Slope_PlasmaRNASeq_NMonly.txt
output/PredictActiveGenes_All/Merged_Controls/MergedControls_prediction_data.csv:
	./scripts/combine_prediction_parameters.py -o output/PredictActiveGenes/Merged_Controls/MergedControls_prediction_data.csv \
    -bcov output/PredictActiveGenes_All/Merged_Controls/MergedControls_TSSCoverage_PlasmaRNASeq_NMonly.txt \
    -scov output/PredictActiveGenes_All/Merged_Controls/MergedControls_TSSCoverage_small_PlasmaRNASeq_NMonly.txt \
    -slope output/PredictActiveGenes_All/Merged_Controls/MergedControls_LowPass_Slope_PlasmaRNASeq_NMonly.txt

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
# Create wig file for ordered array at Gaffney site: Hg19:chr12:34,484,733-34,560,733 
output/Gaffney_nucleosome_array/MergedMaleControls_Gaffney_site.counts:
	./scripts/single_region_wig_from_bam.py -chr chr12 -s 34484000 -e 34561000 -b output/trimmed_reads/Merged_Male_Controls_rmdup_trimmed.bam -w output/Gaffney_nucleosome_array/MergedMaleControls_Gaffney_site.wig
	./scripts/wiggleToCounts.py -w output/Gaffney_nucleosome_array/MergedMaleControls_Gaffney_site.wig -o output/Gaffney_nucleosome_array/MergedMaleControls_Gaffney_site.counts -s 34484733 -e 34560733
output/Gaffney_nucleosome_array/MergedControls_Gaffney_site.counts:
	./scripts/single_region_wig_from_bam.py -chr chr12 -s 34484000 -e 34561000 -b output/trimmed_reads/Merged_Controls_rmdup_trimmed.bam -w output/Gaffney_nucleosome_array/MergedControls_Gaffney_site.wig
	./scripts/wiggleToCounts.py -w output/Gaffney_nucleosome_array/MergedControls_Gaffney_site.wig -o output/Gaffney_nucleosome_array/MergedControls_Gaffney_site.counts -s 34484733 -e 34560733
output/Gaffney_nucleosome_array/MergedFemaleControls_Gaffney_site.counts:
	./scripts/single_region_wig_from_bam.py -chr chr12 -s 34484000 -e 34561000 -b output/trimmed_reads/Merged_Female_Controls_rmdup_trimmed.bam -w output/Gaffney_nucleosome_array/MergedFemaleControls_Gaffney_site.wig
	./scripts/wiggleToCounts.py -w output/Gaffney_nucleosome_array/MergedFemaleControls_Gaffney_site.wig -o output/Gaffney_nucleosome_array/MergedFemaleControls_Gaffney_site.counts -s 34484733 -e 34560733
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

