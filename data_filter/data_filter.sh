#!/usr/bin/env bash

# Example data filtering code
# This code runs cutadapt, trim_galore, and STAR/featureCounts for initial QC and filtering.
# Adjust paths and parameters as needed.

for i in CGS431 CGS432
do
    mkdir ${i}
    mkdir ${i}/01_cutted
    cutadapt -j 20 -u 20 -U 20 -o ${i}/01_cutted/${i}_1_cut.fq.gz -p ${i}/01_cutted/${i}_2_cut.fq.gz ../raw_data/${i}_1.fq.gz ../raw_data/${i}_2.fq.gz
    mkdir ${i}/02_trimmed
    mkdir ${i}/03_mapped
    mkdir ${i}/04_counts

    trim_galore -q 25 --phred33 -e 0.1 --stringency 4 --gzip --trim-n -o ${i}/02_trimmed/ --paired ${i}/01_cutted/${i}_1_cut.fq.gz ${i}/01_cutted/${i}_2_cut.fq.gz

    STAR --runThreadN 15 --runMode alignReads --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM SortedByCoordinate --genomeDir /home/data/genome/STAR_RNA_UCSCfa_GENCODEgtf_hg38_index --readFilesIn ${i}/02_trimmed/${i}_1_cut_val_1.fq.gz ${i}/02_trimmed/${i}_2_cut_val_2.fq.gz --outFileNamePrefix ${i}/03_mapped/

    samtools sort -@ 20 ${i}/03_mapped/Aligned.sortedByCoord.out.bam -o ${i}/03_mapped/${i}.sort.bam

    featureCounts -T 10 \
    -a /home/data/genome/GENCODE/gencode_v41_gtf/gencode.v41.annotation.gtf \
    -p -C -t exon -g gene_name \
    -o ${i}/04_counts/${i}.txt ${i}/03_mapped/${i}.sort.bam
done
