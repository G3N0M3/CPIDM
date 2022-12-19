#######
# Analysis for Immune cells (CD45+)
#######


# Create Directory --------------------------------------------------------
if (!dir.exists(paste0("./output/", ver, "/Islet"))) {
  dir.create(paste0("./output/", ver, "/Islet"))
  print(paste("Output GO directory for version", ver, "created"))
} else {
  print("DIRECTORY ALREADY EXISTS!!!")
  print("Further analysis will OVERWRITE existing files")
}


# Subset ------------------------------------------------------------------
# islet cells
table_Is <- sampleTable[sampleTable$cellType == "Islet", ]


# DESeqDataSet ------------------------------------------------------------------
dds_Is <- DESeqDataSetFromHTSeqCount(sampleTable = table_Is,
                                     directory = "./data/count/",
                                     design = ~ condition)


# Pre-processing ----------------------------------------------------------
# removing rows that have only 0 or 1 read
dds_filt_Is <- dds_Is[rowSums(counts(dds_Is)) > 1, ]
# factoring
dds_filt_Is$condition <- relevel(dds_filt_Is$condition, ref = "Islet_CTR")


# DESeq2 ------------------------------------------------------------------
DE_Is <- DESeq(dds_filt_Is)


# Results -----------------------------------------------------------------
res_Is <- results(DE_Is)


# regularized log transformation ------------------------------------------
rld_Is <- rlog(DE_Is, blind=TRUE)


# Heatmap -----------------------------------------------------------------
# select genes
select <- order(rowMeans(counts(DE_Is, normalized = TRUE)),
                decreasing = TRUE)[1:30]
# plot heatmap
png(filename = paste0("./output/", ver, "/Islet/heatmap_Is.png"),
    width = 5000, height = 4000, units = "px",
    bg = "white", res = 600)
pheatmap(assay(rld_Is)[select, ],
         cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE,
         annotation_col = as.data.frame(colData(DE_Is))[, c("age", "treatment")])
dev.off()


# Gene symbol -------------------------------------------------------------
# Refine
DE_annot_Is <- DE_Is
row.names(DE_annot_Is) <- sapply(row.names(DE_annot_Is),
                                 FUN = function(x) {strsplit(x, "\\.")[[1]][1]},
                                 USE.NAMES = FALSE)
sampNames_Is <- c("IsCTR_1", "IsCTR_2", "IsCTR_3", "IsCTR_4",
                  "IsCTLA4_1", "IsCTLA4_2", "IsCTLA4_3", "IsCTLA4_4",
                  "IsPDL1_1", "IsPDL1_2", "IsPDL1_3")
colnames(DE_annot_Is) <- sampNames_Is
# Annotation
rowData(DE_annot_Is)$symbols <- mapIds(org.Mm.eg.db,  # AnnotationDb object
                                       keys = row.names(DE_annot_Is),  # keys to select from db
                                       column = "SYMBOL",  # the column to search on (match)
                                       keytype = "ENSEMBL",  # annotation type of keys
                                       multiVals = "first")  # action upon multiple mapping
row.names(DE_annot_Is)[!is.na(rowData(DE_annot_Is)$symbols)] <-
  rowData(DE_annot_Is)$symbols[!is.na(rowData(DE_annot_Is)$symbols)]


# Heatmap with Symbol -----------------------------------------------------
# select genes
select <- order(rowMeans(counts(DE_annot_Is, normalized = TRUE)),
                decreasing = TRUE)[1:30]
# rlogTransform
rld_annot_Is <- rlog(DE_annot_Is, blind=TRUE)
png(filename = paste0("./output/", ver, "/Islet/heatmap_annot_Is.png"),
    width = 5000, height = 4000, units = "px",
    bg = "white", res = 600)
pheatmap(assay(rld_annot_Is)[select, ],
         cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(DE_annot_Is))[, c("age", "treatment")])
dev.off()


# sample-to-sample --------------------------------------------------------
# calculate sample distance
sampleDists <- dist(t(assay(rld_Is)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld_Is$condition
colnames(sampleDistMatrix) <- NULL

## save plot
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(25)
# clustered
png(filename = paste0("./output/", ver, "/Islet/sample_annot_Is.png"),
    width = 6000, height = 5000, units = "px",
    bg = "white", res = 600)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, fontsize = 15)
dev.off()
# not-clustered
png(filename = paste0("./output/", ver, "/Islet/sample_annot_unclust_Is.png"),
    width = 6000, height = 5000, units = "px",
    bg = "white", res = 600)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize = 15)
dev.off()


# PCA Plot ----------------------------------------------------------------
ntop_ratio <- 0.05
ntop_Is <- floor(nrow(rld_Is) * ntop_ratio)
pca_Is <- plotPCA(rld_Is,
                  intgroup = c("treatment"),
                  ntop = ntop_Is)

pca_data_Is <- plotPCA(rld_Is, intgroup = c("treatment"), ntop = ntop_Is)

pca_Is <-
  pca_Is +
  theme_light() +
  labs(color = "Treatment",
       title = paste("ntop =", ntop_Is)) +
  scale_color_discrete(breaks = c("PDL1", "CTLA4", "CTR")) +
  theme(text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))

ggsave("pca_Is.png",
       #plot = last_plot(),
       plot = pca_Is,
       device = "png",
       path = paste0("./output/", ver, "/Islet"),
       width = 3500, height = 2000, units = "px", dpi = 600,
       scale = 1.25
)


# DEG List ----------------------------------------------------------------
res_annot_Is <- results(DE_annot_Is)
csv_Is <- as.data.frame(res_annot_Is[abs(res_annot_Is$log2FoldChange) > 1 &
                                       !is.na(res_annot_Is$padj) &
                                       res_annot_Is$padj < 0.1, ])
write.csv(x = csv_Is[order(csv_Is$pvalue, decreasing = FALSE), ],
          file = paste0("./output/", ver, "/Islet/DEG_Is.csv"))
