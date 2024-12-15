# Set the working directory (please modify to your actual path)
setwd("/Users/yaoshuo/Desktop/CGS617-620")

# Load required packages
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(tidyr)

########################################
# Basic RNA-seq data import and PCA analysis
########################################

# Example file paths (modify to actual file names/paths)
WThesc_1 <- "CGS223counts.txt"
WThesc_2 <- "CGS224counts.txt"
endo_1 <- "CGS243.txt"
endo_2 <- "CGS244.txt"
meso_1 <- "CGS247.txt"
meso_2 <- "CGS248.txt"
ecto_1 <- "CGS266.txt"
ecto_2 <- "CGS267.txt"
rac1_1 <- "CGS617_counts.txt"
rac1_2 <- "CGS618_counts.txt"
diaph1_1 <- "CGS619_counts.txt"
diaph1_2 <- "CGS620_counts.txt"

# Read data and extract count columns
WThesc_1 <- read.table(WThesc_1, header = TRUE, row.names = 1)[, -4:-1]
WThesc_2 <- read.table(WThesc_2, header = TRUE, row.names = 1)[, -4:-1]
endo_1 <- read.table(endo_1, header = TRUE, row.names = 1)[, -4:-1]
endo_2 <- read.table(endo_2, header = TRUE, row.names = 1)[, -4:-1]
meso_1 <- read.table(meso_1, header = TRUE, row.names = 1)[, -4:-1]
meso_2 <- read.table(meso_2, header = TRUE, row.names = 1)[, -4:-1]
ecto_1 <- read.table(ecto_1, header = TRUE, row.names = 1)[, -4:-1]
ecto_2 <- read.table(ecto_2, header = TRUE, row.names = 1)[, -4:-1]
rac1_1 <- read.table(rac1_1, header = TRUE, row.names = 1)[, -4:-1]
rac1_2 <- read.table(rac1_2, header = TRUE, row.names = 1)[, -4:-1]
diaph1_1 <- read.table(diaph1_1, header = TRUE, row.names = 1)[, -4:-1]
diaph1_2 <- read.table(diaph1_2, header = TRUE, row.names = 1)[, -4:-1]

# Combine count data
count <- cbind(WThesc_1, WThesc_2, endo_1, endo_2, meso_1, meso_2, ecto_1, ecto_2,
               rac1_1, rac1_2, diaph1_1, diaph1_2)

# Remove the Length column if present
column_index <- which(names(count) == "Length")
if(length(column_index) > 0) {
  count <- count[, -column_index]
}

# Rename columns
colnames(count) <- c("control hESC_1", "control hESC_2", "Endoderm_1", "Endoderm_2",
                     "Mesoderm_1", "Mesoderm_2", "Ectoderm_1", "Ectoderm_2",
                     "RAC1-KD_1", "RAC1-KD_2", "DIAPH1-KD_1", "DIAPH1-KD_2")

# Filter out genes with no expression
count <- count[apply(count, 1, function(x) sum(x > 1) > 0), ] %>% as.matrix()

# Create metadata for samples
info <- data.frame(row.names = colnames(count),
                   celltype = c(rep("control hESC", 2), rep("Endoderm", 2),
                                rep("Mesoderm", 2), rep("Ectoderm", 2),
                                rep("RAC1-KD", 2), rep("DIAPH1-KD", 2)),
                   sample = colnames(count))
info$celltype <- as.factor(info$celltype)

# DESeq2 analysis and PCA
ddsOrig <- DESeqDataSetFromMatrix(countData = count, colData = info, design = ~ celltype)
dds <- DESeq(ddsOrig)
vstd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vstd, intgroup = c('celltype'), returnData = TRUE)
pca <- prcomp(t(assay(vstd)))
summary(pca)

# Custom colors for cell types
custom_colors <- c("control hESC" = "skyblue", "Endoderm" = "brown",
                   "Mesoderm" = "pink", "Ectoderm" = "cyan",
                   "RAC1-KD" = "red", "DIAPH1-KD" = "green")

# Plot PCA
ggplot(pcaData, aes(PC1, PC2, color = celltype)) +
  theme_classic(base_size = 15) +
  geom_point(size = 7, alpha = 0.5) +
  geom_text_repel(aes(label = celltype), show.legend = FALSE, max.overlaps = 20) +
  scale_color_manual(values = custom_colors, breaks = names(custom_colors)) +
  labs(x = 'PCA1:49%', y = 'PCA2:41%')


########################################
# Differential expression and GO enrichment analysis
########################################

# Set working directory (modify to your path)
setwd("/Users/yaoshuo/Desktop/CGS621-626/DEG+GO正常")

# Read count data for DEG analysis
CGS623 <- read.table("CGS623_counts.txt", header = TRUE, sep = "\t", row.names = 1)
CGS624 <- read.table("CGS624_counts.txt", header = TRUE, sep = "\t", row.names = 1)
CGS625 <- read.table("CGS625_counts.txt", header = TRUE, sep = "\t", row.names = 1)
CGS626 <- read.table("CGS626_counts.txt", header = TRUE, sep = "\t", row.names = 1)

CGS623_counts <- as.numeric(CGS623[, ncol(CGS623)])
CGS624_counts <- as.numeric(CGS624[, ncol(CGS624)])
CGS625_counts <- as.numeric(CGS625[, ncol(CGS625)])
CGS626_counts <- as.numeric(CGS626[, ncol(CGS626)])

count_data <- cbind(CGS623_counts, CGS624_counts, CGS625_counts, CGS626_counts)
row.names(count_data) <- row.names(CGS623)
count_data_clean <- na.omit(count_data)

colnames(count_data_clean) <- c("control_rep1", "control_rep2", "treat_rep1", "treat_rep2")
condition <- factor(c("control", "control", "treat", "treat"))
col_data <- data.frame(condition = condition)

dds <- DESeqDataSetFromMatrix(countData = count_data_clean, colData = col_data, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
write.csv(as.data.frame(res), file = "BFP_ACTA_DEGs.csv", row.names = TRUE)
resOrdered <- res[order(res$padj),]
resSig <- subset(resOrdered, padj < 0.05)

# Export up- and down-regulated genes
write.csv(subset(resSig, log2FoldChange > 2), file = "BFP_ACTA_upDEGs.csv", row.names = TRUE)
write.csv(subset(resSig, log2FoldChange < -2), file = "BFP_ACTA_downDEGs.csv", row.names = TRUE)

# Generate logTPM data (gene_lengths required)
gene_lengths <- CGS623$Length
rpk <- count_data_clean[, 1:4] / gene_lengths
scaling_factors <- colSums(rpk) / 1e6
tpm <- rpk / scaling_factors
logTPM <- log(tpm + 1)
write.csv(logTPM, file = "logTPM_data.csv")

# Mark up/down genes for scatter plot
upDEGs <- read.csv("BFP_ACTA_upDEGs.csv", header = TRUE, row.names = 1)
downDEGs <- read.csv("BFP_ACTA_downDEGs.csv", header = TRUE, row.names = 1)
test <- read.csv("logTPM_data.csv", header = TRUE, row.names = 1)

Up <- rownames(upDEGs)
Down <- rownames(downDEGs)
test$ect.DEG0.05 <- 'no diff'
test[Up, "ect.DEG0.05"] <- "Up(0.05)"
test[Down, "ect.DEG0.05"] <- "Down(0.05)"

test$wt_mean <- (test$control_rep1 + test$control_rep2)/2
test$d1_mean <- (test$treat_rep1 + test$treat_rep2)/2

pdf("BFP_ACTA_scatter_plot123.pdf", width = 7, height = 7)
col0 <- rep(rgb(244/255, 243/255, 183/255, 0.22), nrow(test))
col0[test$ect.DEG0.05 == "Down(0.05)"] <- rgb(90/255, 124/255, 175/255, 0.98)
col0[test$ect.DEG0.05 == "Up(0.05)"] <- rgb(201/255, 100/255, 85/255, 0.98)
plot(test$wt_mean, test$d1_mean, xlab="", ylab="", col=NA, bg=col0, pch=21, main="BFP_ACTA (|FC|>2;padj<0.05)")
title(xlab="Control mean exp(log(TPM)+1)", ylab="Treatment mean exp(log(TPM)+1)", line=2.3, cex.lab=1.5)
legend("topleft", c("upDEG", "downDEG"), pt.bg=c("#FF000090","#0000FF90"), pch=c(21,21), cex=1.0)
dev.off()

# GO enrichment analysis
resSig <- read.csv("BFP_ACTA_DEGs.csv", header = TRUE, row.names = 1)

# GO for up-regulated genes
up <- subset(resSig, log2FoldChange > 2)
up$id <- mapIds(org.Hs.eg.db, keys = row.names(up), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
id2 <- unique(na.omit(up$id))

pdf("BFP_ACTA_Upregulated_Genes_GOterm_Enrichment.pdf")
go.up.bp <- enrichGO(gene = id2, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                     ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(go.up.bp, title = "Up-regulated Genes (BP) GO Enrichment")

go.up.mf <- enrichGO(gene = id2, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                     ont = "MF", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(go.up.mf, title = "Up-regulated Genes (MF) GO Enrichment")

go.up.cc <- enrichGO(gene = id2, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                     ont = "CC", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(go.up.cc, title = "Up-regulated Genes (CC) GO Enrichment")
dev.off()

# GO for down-regulated genes
down <- subset(resSig, log2FoldChange < -2)
down$id <- mapIds(org.Hs.eg.db, keys = row.names(down), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
id_down <- unique(na.omit(down$id))

pdf("BFP_ACTA_Downregulated_Genes_GOterm_Enrichment.pdf")
go.down.bp <- enrichGO(gene = id_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                       ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(go.down.bp, title = "Down-regulated Genes (BP) GO Enrichment")

go.down.mf <- enrichGO(gene = id_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                       ont = "MF", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(go.down.mf, title = "Down-regulated Genes (MF) GO Enrichment")

go.down.cc <- enrichGO(gene = id_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                       ont = "CC", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(go.down.cc, title = "Down-regulated Genes (CC) GO Enrichment")
dev.off()


########################################
# Boxplot example
########################################
setwd("D:/暑研/editing/新/bed/hPGCLC/boxplot/1206/2")
data <- read.csv("filtered_boxplot-2.csv")
names(data)[1] <- "Geneid"
data_long <- data %>%
  gather(key = "sample", value = "expression", -Geneid) %>%
  separate(sample, into = c("group", "replicate"), sep = "_", remove = FALSE) %>%
  group_by(Geneid, group) %>% summarise(mean_expression = mean(expression, na.rm = TRUE)) %>% 
  spread(key = group, value = mean_expression)

compare_groups <- c("CE.hESC")
p_values <- data_long %>% summarise(across(all_of(compare_groups), ~ t.test(.x, data_long$N3.ADAR.hESC)$p.value))
p_values <- c(CE_hESC = p_values$CE.hESC)

new_column_names <- sapply(names(p_values), function(name) paste0(name, " (P=", sprintf("%.2f", p_values[[name]]), ")"))
names(data_long)[names(data_long) %in% compare_groups] <- new_column_names
data_long <- data_long %>% gather(key = "group", value = "expression", -Geneid)
data_long$group <- factor(data_long$group, levels = c("N3.ADAR.hESC", new_column_names))

pdf(file = "FigureS5H-short-1.pdf", width = 4, height = 7)
ggplot(data_long, aes(x = group, y = expression, color = group)) +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  theme_minimal() + labs(x = "", y = "log2(tpm+1)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), 
        legend.position = "right") +
  geom_text(data = data.frame(group = new_column_names, expression = rep(12.5, length(new_column_names))),
            aes(x = group, y = expression, label = paste("P =", sprintf("%.2f", p_values))),
            color = "black", vjust = -0.5 )
dev.off()


########################################
# Heatmap example
########################################
setwd("D:/暑研/editing/新/D1_genelist_heatmap/Figure1J/1")

files <- list.files(path = ".", pattern = "\\.txt$", full.names = TRUE)
expr <- lapply(files, function(x) {
  tmp <- read.table(file = x, header = TRUE, skip = 1)
  tmp <- tmp[, c(1, 7)]
  colnames(tmp) <- c("Geneid", basename(x))
  return(tmp)
})

df <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), expr)

poluripotency_genes <- c("NANOG", "POU5F1", "SOX2")
germline_genes <- c("SOX17", "TFAP2C", "PRDM1", "NANOS3", "DAZL", "DDX4")
genes <- c(poluripotency_genes, germline_genes)

df <- column_to_rownames(df, var = "Geneid")
filtered_data <- df[rownames(df) %in% genes, ]
log_transformed_data <- log2(filtered_data + 1)
ordered_data <- log_transformed_data[genes, , drop = FALSE]

row_split <- factor(c(rep("Poluripotency", length(poluripotency_genes)),
                      rep("Germline", length(germline_genes))),
                    levels = c("Poluripotency", "Germline"))

colnames(ordered_data) <- c("Control_hESC1", "Control_hESC2", "D1-KO_hESC1", "D1-KO_hESC2",
                            "Control_hPGCLC1", "Control_hPGCLC2", "D1-KO_hPGCLC1", "D1-KO_hPGCLC2")

pdf(file="Figure1L-heatmap-6.pdf", width=8, height=8)
Heatmap(as.matrix(ordered_data),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        row_split = row_split,
        column_names_rot = 60,
        row_gap = unit(3, "mm"),
        name = "log2(TPM+1)",
        rect_gp = gpar(col = "white", lwd = 6),
        col = colorRamp2(c(0, 8, 16), c("#3A70A2", "#FFFFCC", "#d93506")),
        row_title_gp = gpar(fontsize = 15),
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 12))
dev.off()
