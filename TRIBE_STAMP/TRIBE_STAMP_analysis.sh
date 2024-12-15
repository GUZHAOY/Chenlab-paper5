#!/usr/bin/env bash

# This script demonstrates TRIBE-STAMP analysis using Bullseye.
# Please adjust paths, annotation files, and sample names accordingly.

conda activate bullseye

for i in CGS611 CGS612 CGS613 CGS614 CGS615 CGS616 CGS617 CGS618 CGS619 CGS620
do
    cd ${i}/03_mapped
    perl /home/storage_3/CGS418-428/analysis/Bullseye/Bullseye/Bul/pl/parseBAM.pl \
    --input Aligned.sortedByCoord.out.bam \
    --output output.matrix \
    --cpu 60 \
    --minCoverage 10 \
    --removeDuplicates
    cd ../..
done

# Finding editing sites with Bullseye
perl /home/bioinfo_2/zhaoyu/Bul/pl/Find_edit_site.pl \
--annotationFile /home/bioinfo_2/zhaoyu/Bul/index/annotation.refFlat \
--EditedMatrix /home/storage_3/CGS530-536/analysis/CGS534/03_mapped/output.matrix.gz \
--controlMatrix /home/bioinfo_2/zhaoyu/CGS611-620/analysis/03_mapped/CGS611/03_mapped/output.matrix.gz \
--minEdit 5 \
--maxEdit 90 \
--editFoldThreshold 3 \
--MinEditSites 3 \
--cpu 60 \
--outfile /home/bioinfo_2/zhaoyu/CGS611-620/analysis/03_mapped/CE-D1-ADAR-hESC/CGS534-WT1.bed \
--fallback /home/bioinfo_2/zhaoyu/Bul/index/ \
--verbose \
--editType A2G

# Repeat the Find_edit_site.pl command for other samples as needed
perl /home/bioinfo_2/zhaoyu/Bul/pl/Find_edit_site.pl \
--annotationFile /home/bioinfo_2/zhaoyu/Bul/index/annotation.refFlat \
--EditedMatrix /home/storage_3/CGS530-536/analysis/CGS535/03_mapped/output.matrix.gz \
--controlMatrix /home/bioinfo_2/zhaoyu/CGS611-620/analysis/03_mapped/CGS611/03_mapped/output.matrix.gz \
--minEdit 5 \
--maxEdit 90 \
--editFoldThreshold 3 \
--MinEditSites 3 \
--cpu 60 \
--outfile /home/bioinfo_2/zhaoyu/CGS611-620/analysis/03_mapped/CE-D1-ADAR-hESC/CGS353-WT1.bed \
--fallback /home/bioinfo_2/zhaoyu/Bul/index/ \
--verbose \
--editType A2G

# ... (other Find_edit_site.pl calls)

perl /home/bioinfo_2/zhaoyu/Bul/pl/summarize_sites.pl --MinRep 2 --mut 3 CGS253-WT1.bed CGS254-WT1.bed > D1-rAPOBEC1-1.bed
perl /home/bioinfo_2/zhaoyu/Bul/pl/summarize_sites.pl --MinRep 2 --mut 3 CGS253-WT2.bed CGS254-WT2.bed > D1-rAPOBEC1-2.bed

# Run R script for summarizing results (example)
Rscript -e "
library(tidyverse)
library(openxlsx)
setwd('D:/暑研/editing/新/bed')
sites <- read_delim('D1-rAPOBEC1-2.bed', delim = '\t',
  col_names = c('chr','start','end','name','score','strand','control','control_cov','score2','cov','type','n_rep')) %>%
  separate(name, into=c('Gene_name','biotype'), sep='([|])') %>%
  mutate(site_id=paste(chr,start, sep = '_'), score=round(score*100,digits=2)) %>%
  group_by(Gene_name) %>%
  summarise(Num_edit_sites=n(),
            Avg_edting=round(mean(score),digits=2),
            Features=paste(biotype,collapse=';'),
            Edit_percent_read_str=paste(score,collapse=';'),
            Identifier_str=paste(site_id,collapse=';'))
write.xlsx(sites, 'D1-rAPOBEC1-2.xlsx', colNames=TRUE)
"

# Co-editing analysis
perl /home/bioinfo_2/zhaoyu/Bul/pl/co_editing.pl \
--bed1 /home/bioinfo_2/zhaoyu/CGS611-620/analysis/03_mapped/CE-D1-rAPOBEC1-hESC/CE-D1-rAPOBEC1-1.bed \
--bed2 /home/bioinfo_2/zhaoyu/CGS611-620/analysis/03_mapped/CE-N3-ADAR-hESC/CE-N3-ADAR-1.bed \
--bam /home/bioinfo_2/zhaoyu/Bul/D1-N3-hESC/CGS380/Aligned.sortedByCoord.out.bam \
--outfile D1rAPO+N3ADAR-1.tsv \
--refFlat /home/bioinfo_2/zhaoyu/Bul/index/annotation.refFlat \
--minCov 20 \
--distance 150 \
--removeDup \
--removeMultiMapped

perl /home/bioinfo_2/zhaoyu/Bul/pl/co_editing.pl \
--bed1 /home/bioinfo_2/zhaoyu/CGS611-620/analysis/03_mapped/CE-D1-rAPOBEC1-hESC/CE-D1-rAPOBEC1-1.bed \
--bed2 /home/bioinfo_2/zhaoyu/CGS611-620/analysis/03_mapped/CE-N3-ADAR-hESC/CE-N3-ADAR-1.bed \
--bam /home/bioinfo_2/zhaoyu/Bul/D1-N3-hESC/CGS381/Aligned.sortedByCoord.out.bam \
--outfile D1rAPO+N3ADAR-2.tsv \
--refFlat /home/bioinfo_2/zhaoyu/Bul/index/annotation.refFlat \
--minCov 20 \
--distance 150 \
--removeDup \
--removeMultiMapped
