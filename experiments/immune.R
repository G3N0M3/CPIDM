#######
# Analysis for Immune cells (CD45+)
#######


# Create Directory --------------------------------------------------------
if (!dir.exists(paste0("./output/", ver, "/Immune"))) {
  dir.create(paste0("./output/", ver, "/Immune"))
  print(paste("Output GO directory for version", ver, "created"))
} else {
  print("DIRECTORY ALREADY EXISTS!!!")
  print("Further analysis will OVERWRITE existing files")
}


# Subset ------------------------------------------------------------------
# immune cells
table_Im <- sampleTable[sampleTable$cellType == "Immune", ]


# DESeqDataSet ------------------------------------------------------------------
dds_Im <- DESeqDataSetFromHTSeqCount(sampleTable = table_Im,
                                     directory = "./data/count/",
                                     design = ~ condition)


# Pre-processing ----------------------------------------------------------
# removing rows that have only 0 or 1 read
dds_filt_Im <- dds_Im[rowSums(counts(dds_Im)) > 1, ]
# factoring
dds_filt_Im$condition <- relevel(dds_filt_Im$condition, ref = "Immune_CTR")


# DESeq2 ------------------------------------------------------------------
DE_Im <- DESeq(dds_filt_Im)


# Results -----------------------------------------------------------------
res_Im <- results(DE_Im)


# regularized log transformation ------------------------------------------
rld_Im <- rlog(DE_Im, blind=TRUE)


# Heatmap -----------------------------------------------------------------
# select genes
select <- order(rowMeans(counts(DE_Im, normalized = TRUE)),
                decreasing = TRUE)[1:30]
# plot heatmap
png(filename = paste0("./output/", ver, "/Immune/heatmap_Im.png"),
    width = 5000, height = 4000, units = "px",
    bg = "white", res = 600)
pheatmap(assay(rld_Im)[select, ],
         cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE,
         annotation_col = as.data.frame(colData(DE_Im))[, c("age", "treatment")],
         fontsize = 12)
dev.off()


# Gene symbol -------------------------------------------------------------
# Refine
DE_annot_Im <- DE_Im
row.names(DE_annot_Im) <- sapply(row.names(DE_annot_Im),
                                 FUN = function(x) {strsplit(x, "\\.")[[1]][1]},
                                 USE.NAMES = FALSE)
sampNames_Im <- c("ImCTR_1", "ImCTR_2", "ImCTR_3", "ImCTR_4",
                  "ImCTLA4_1", "ImCTLA4_2", "ImCTLA4_3", "ImCTLA4_4",
                  "ImPDL1_1", "ImPDL1_2", "ImPDL1_3")
colnames(DE_annot_Im) <- sampNames_Im
# Annotation
rowData(DE_annot_Im)$symbols <- mapIds(org.Mm.eg.db,  # AnnotationDb object
                                       keys = row.names(DE_annot_Im),  # keys to select from db
                                       column = "SYMBOL",  # the column to search on (match)
                                       keytype = "ENSEMBL",  # annotation type of keys
                                       multiVals = "first")  # action upon multiple mapping
row.names(DE_annot_Im)[!is.na(rowData(DE_annot_Im)$symbols)] <-
  rowData(DE_annot_Im)$symbols[!is.na(rowData(DE_annot_Im)$symbols)]


# Heatmap with Symbol -----------------------------------------------------
# select genes
select <- order(rowMeans(counts(DE_annot_Im, normalized = TRUE)),
                decreasing = TRUE)[1:30]
# rlogTransform
rld_annot_Im <- rlog(DE_annot_Im, blind=TRUE)
png(filename = paste0("./output/", ver, "/Immune/heatmap_annot_Im.png"),
    width = 5000, height = 4000, units = "px",
    bg = "white", res = 600)
pheatmap(assay(rld_annot_Im)[select, ],
         cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 12,
         annotation_col = as.data.frame(colData(DE_annot_Im))[, c("age", "treatment")])
dev.off()


# sample-to-sample --------------------------------------------------------
# calculate sample distance
sampleDists <- dist(t(assay(rld_Im)))
sampleDistMatrix <- as.matrix(sampleDists)
rld_Im$samples <- sampNames_Im
rownames(sampleDistMatrix) <- rld_Im$samples
colnames(sampleDistMatrix) <- NULL

## save plot
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(25)
# clustered
png(filename = paste0("./output/", ver, "/Immune/sample_annot_Im.png"),
    width = 6000, height = 5000, units = "px",
    bg = "white", res = 600)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         cex = 1.25)
dev.off()
# not-clustered
png(filename = paste0("./output/", ver, "/Immune/sample_annot_unclust_Im.png"),
    width = 6000, height = 5000, units = "px",
    bg = "white", res = 600)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         cluster_rows = FALSE, cluster_cols = FALSE,
         cex = 1.25)
dev.off()


# PCA Plot ----------------------------------------------------------------
ntop_ratio_Im <- 0.03
ntop_Im <- floor(nrow(rld_Im) * ntop_ratio)
pca_Im <- plotPCA(rld_Im,
               intgroup = c("treatment"),
               ntop = ntop_Im)

pca_data <- plotPCA(rld_Im, intgroup = c("cellType", "treatment"),
                    ntop = ntop, returnData = TRUE)

pca_Im <- pca_Im +
  theme_light() +
  labs(color = "Treatment",
       title = paste("ntop =", ntop_Im)) +
  scale_color_discrete(breaks = c("PDL1", "CTLA4", "CTR")) +
  theme(text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))

ggsave("pca_Im.png",
       #plot = last_plot(),
       plot = pca_Im,
       device = "png",
       path = paste0("./output/", ver, "/Immune"),
       width = 3200, height = 2400, units = "px", dpi = 600,
       scale = 1.25
)


# DEG List ----------------------------------------------------------------
res_annot_Im <- results(DE_annot_Im)
csv_Im <- as.data.frame(res_annot_Im[abs(res_annot_Im$log2FoldChange) > 1 &
                                       !is.na(res_annot_Im$padj) &
                                       res_annot_Im$padj < 0.05, ])
write.csv(x = csv_Im[order(csv_Im$pvalue, decreasing = FALSE), ],
          file = paste0("./output/", ver, "/Immune/DEG_Im.csv"))
