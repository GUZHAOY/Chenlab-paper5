#!/usr/bin/env bash

# This script demonstrates TRIBE analysis steps.
# Please modify paths, indexes, and file names according to your environment.

bowtie_indexes="/home/data/genome/hg38.bybowtie2/hg38.bybowtie2"

for i in CGS622
do
    mkdir -p ${i}/trim
    trim_galore -q 25 --phred33 --fastqc --stringency 3 --gzip --length 50 --trim-n --max_n 2 \
    --output_dir ${i}/trim --cores 32 --paired \
    ${i}/raw_data/${i}_1.fq.gz \
    ${i}/raw_data/${i}_2.fq.gz

    mkdir -p ${i}/aligned
    bowtie2 --threads 32 -q --phred33 -D 15 -R 2 -N 1 -L 22 -i S,1,1.15 -x $bowtie_indexes \
    -1 ${i}/trim/${i}_1_val_1.fq.gz -2 ${i}/trim/${i}_2_val_2.fq.gz \
    | samtools view -@ 32 -S -h -q 10 -b -o ${i}/aligned/${i}.bam

    samtools sort -@ 32 -n ${i}/aligned/${i}.bam -o ${i}/aligned/${i}.sort.bam
    samtools fixmate -@ 32 -r -m ${i}/aligned/${i}.sort.bam ${i}/aligned/${i}.sort.fixmate.bam
    samtools sort -@ 32 ${i}/aligned/${i}.sort.fixmate.bam -o ${i}/aligned/${i}.sort.fixmate.sort.bam
    samtools markdup -@ 32 -r --output-fmt BAM ${i}/aligned/${i}.sort.fixmate.sort.bam ${i}/aligned/${i}.sort.fixmate.sort.rmdup.bam

    mkdir -p ${i}/samfiles
    samtools view -@ 32 ${i}/aligned/${i}.sort.fixmate.sort.rmdup.bam > ${i}/samfiles/${i}.sam
done

conda activate TRIBEenv
export PATH="/home/bioinfo/zhaoyu/CODE/try:$PATH"

# load_table.sh calls (example, please adjust according to your needs)
load_table.sh /home/bioinfo/zhaoyu/CGS621/samfiles/CGS621.sam ADAR1 rnalibs 1
load_table.sh /home/bioinfo/zhaoyu/CGS622/samfiles/CGS622.sam ADAR1 rnalibs 3
load_table.sh CGS223/samfiles/CGS223.sam T1 rnalibs 1
load_table.sh CGS224/samfiles/CGS224.sam T1 rnalibs 2
load_table.sh CGS391/samfiles/CGS391.sam T1 rnalibs 4
load_table.sh CGS132/samfiles/CGS132.sam T1 rnalibs 3
load_table.sh CGS424/samfiles/CGS424.sam T1 rnalibs 5
load_table.sh CGS425/samfiles/CGS425.sam T1 rnalibs 6
load_table.sh CGS426/samfiles/CGS426.sam T1 rnalibs 7
load_table.sh CGS381.sam T2 rnalibs 4
load_table.sh ../samfiles/CGS094.sam ADAR713 rnalibs 4
load_table.sh /home/bioinfo/zhaoyu/tribe/CGS199/samfiles/CGS199.sam B1 rnalibs 2
load_table.sh /home/bioinfo/zhaoyu/tribe/CGS200/samfiles/CGS200.sam B1 rnalibs 3
load_table.sh /home/bioinfo/zhaoyu/tribe/CGS207/samfiles/CGS207.sam B1 rnalibs 4
load_table.sh /home/bioinfo/zhaoyu/tribe/CGS208/samfiles/CGS208.sam B1 rnalibs 5
load_table.sh /home/bioinfo/zhaoyu/tribe/CGS209/samfiles/CGS209.sam B1 rnalibs 6
load_table.sh /home/bioinfo/zhaoyu/tribe/CGS210/samfiles/CGS210.sam B1 rnalibs 7
load_table.sh /home/bioinfo/zhaoyu/tribe/CGS253/samfiles/CGS253.sam B1 rnalibs 8
load_table.sh /home/bioinfo/zhaoyu/tribe/CGS254/samfiles/CGS254.sam B1 rnalibs 9
load_table.sh /home/bioinfo/zhaoyu/tribe/CGS189/samfiles/CGS189.sam B1 rnalibs 10
load_table.sh /home/bioinfo/zhaoyu/tribe/CGS190/samfiles/CGS190.sam B1 rnalibs 11
load_table.sh /home/bioinfo/zhaoyu/tribe/CGS216/samfiles/CGS216.sam B1 rnalibs 12
load_table.sh /home/bioinfo/zhaoyu/tribe/CGS217/samfiles/CGS217.sam B1 rnalibs 13
load_table.sh /home/bioinfo/zhaoyu/CGS418-428/analysis/CGS420/samfiles/CGS420.sam C1 rnalibs 3
load_table.sh /home/bioinfo/zhaoyu/CGS418-428/analysis/CGS421/samfiles/CGS421.sam C1 rnalibs 4
load_table.sh /home/storage_2/CGS382-389/analysis/fastQC_1/zhaoyu/data/CGS380-381/samfiles/CGS380.sam C1 rnalibs 5
load_table.sh /home/storage_2/CGS382-389/analysis/fastQC_1/zhaoyu/data/CGS380-381/samfiles/CGS381.sam B1 rnalibs 6
load_table.sh /home/storage_2/CGS382-389/analysis/fastQC_1/zhaoyu/data/CGS380-381/samfiles/CGS381.sam B1 rnalibs 8
load_table.sh /home/bioinfo/zhaoyu/CGS431-432/analysis/CGS431/samfiles/CGS431.sam C3 rnalibs 1
load_table.sh /home/bioinfo/zhaoyu/CGS431-432/analysis/CGS432/samfiles/CGS432.sam C3 rnalibs 2
load_table.sh /home/storage_2/CGS382-389/analysis/fastQC_1/zhaoyu/tribe/CGS199/samfiles/CGS199.sam C3 rnalibs 3
load_table.sh /home/storage_2/CGS382-389/analysis/fastQC_1/zhaoyu/tribe/CGS200/samfiles/CGS200.sam C3 rnalibs 4

###################################
# STEP4: Find RNA editing sites
###################################
conda activate TRIBEenv
rnaedit_wtRNA_RNA.sh

# Example bedtools operations (modify file names accordingly)
bedtools intersect -wa -wb -f 0.9 -r -a rnalibs_3_1_A2G_1%.bedgraph -b rnalibs_3_2_A2G_1%.bedgraph > 1.bedgraph
bedtools intersect -wa -v -a 1.bedgraph -b /home/bioinfo/zhaoyu/TRIBE/ADAR/rnalibs_3_1_A2G_1%.bedgraph > 2.bedgraph
bedtools intersect -wa -v -a 2.bedgraph -b /home/bioinfo/zhaoyu/TRIBE/ADAR/rnalibs_3_3_A2G_1%.bedgraph > D1-ADAR-hPGCLC-rep1.bedgraph
summarize_results.pl D1-ADAR-hPGCLC-rep1.bedgraph > D1-ADAR-hPGCLC_rep1.xls

# Add any additional TRIBE-related scripts or analysis steps here
