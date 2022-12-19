  #######
# Analysis for Total
#######


# Create Directory --------------------------------------------------------
if (!dir.exists(paste0("./output/", ver, "/multi"))) {
  dir.create(paste0("./output/", ver, "/multi"))
  print(paste("Output GO directory for version", ver, "created"))
} else {
  print("DIRECTORY ALREADY EXISTS!!!")
  print("Further analysis will OVERWRITE existing files")
}

  
# DESeqDataSet ------------------------------------------------------------------
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = "./data/count/",
                                  design = ~ condition)


# Pre-processing ----------------------------------------------------------
# removing rows that have only 0 or 1 read
dds_filt <- dds[rowSums(counts(dds)) > 1, ]
# factoring
dds_filt$condition <- relevel(dds_filt$condition, ref = "Immune_CTR")


# DESeq2 ------------------------------------------------------------------
DE <- DESeq(dds_filt)


# Results -----------------------------------------------------------------
res <- results(DE)


# normTransform -----------------------------------------------------------
ntd <- normTransform(DE)
#meanSdPlot(assay(ntd))


# Heatmap -----------------------------------------------------------------
# select genes
select <- order(rowMeans(counts(DE, normalized = TRUE)),
                decreasing = TRUE)[1:30]
# regularized log transformation
rld <- rlog(DE, blind=TRUE)
png(filename = paste0("./output/", ver, "/multi/heatmap_multi.png"),
    width = 5000, height = 4000, units = "px",
    bg = "white", res = 600)
pheatmap(
  assay(rld)[select, ],
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_col = as.data.frame(colData(DE))[, c("age", "cellType", "treatment")])
dev.off()


# Gene symbol -------------------------------------------------------------
# Refine
DE_annot <- DE
row.names(DE_annot) <- sapply(row.names(DE_annot),
                        FUN = function(x) {strsplit(x, "\\.")[[1]][1]},
                        USE.NAMES = FALSE)
sampNames <- c("ImCTR_1", "ImCTR_2", "ImCTR_3", "ImCTR_4",
               "ImCTLA4_1", "ImCTLA4_2", "ImCTLA4_3", "ImCTLA4_4",
               "ImPDL1_1", "ImPDL1_2", "ImPDL1_3",
               "IsCTR_1", "IsCTR_2", "IsCTR_3", "IsCTR_4",
               "IsCTLA4_1", "IsCTLA4_2", "IsCTLA4_3", "IsCTLA4_4",
               "IsPDL1_1", "IsPDL1_2", "IsPDL1_3")
colnames(DE_annot) <- sampNames

# Annotation
rowData(DE_annot)$symbols <- mapIds(org.Mm.eg.db,  # AnnotationDb object
                                    keys = row.names(DE_annot),  # keys to select from db
                                    column = "SYMBOL",  # the column to search on (match)
                                    keytype = "ENSEMBL",  # annotation type of keys
                                    multiVals = "first")  # action upon multiple mapping
row.names(DE_annot)[!is.na(rowData(DE_annot)$symbols)] <-
  rowData(DE_annot)$symbols[!is.na(rowData(DE_annot)$symbols)]


# Heatmap with Symbol -----------------------------------------------------
# regularized log transformation
rld_annot <- rlog(DE_annot, blind=TRUE)
png(filename = paste0("./output/", ver, "/multi/heatmap_annot_multi.png"),
    width = 5000, height = 4000, units = "px",
    bg = "white", res = 600)
pheatmap(assay(rld_annot)[select, ],
         cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = FALSE,
         annotation_col = as.data.frame(colData(DE_annot))[, c("age", "cellType", "treatment")]
         )
dev.off()


# sample-to-sample --------------------------------------------------------
# calculate sample distance
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rld$samples <- sampNames
rownames(sampleDistMatrix) <- rld$samples
colnames(sampleDistMatrix) <- NULL

## save plot
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(25)
# clustered
png(filename = paste0("./output/", ver, "/multi/sample_annot_multi.png"),
    width = 6000, height = 5000, units = "px",
    bg = "white", res = 600)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         cex = 1.25)
dev.off()
# not clustered
png(filename = paste0("./output/", ver, "/multi/sample_annot_unclust_multi.png"),
    width = 6000, height = 5000, units = "px",
    bg = "white", res = 600)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         cluster_cols = FALSE, cluster_rows = FALSE,
         cex = 1.25)
dev.off()


# PCA Plot ----------------------------------------------------------------
ntop_ratio <- 0.03
ntop <- floor(nrow(rld) * ntop_ratio)
pca <- plotPCA(rld,
               intgroup = c("cellType", "treatment"),
               ntop = ntop)

pca_data <- plotPCA(rld, intgroup = c("cellType", "treatment"),
                    ntop = ntop, returnData = TRUE)

pca <- pca +
  theme_light() +
  labs(color = "Treatment",
       title = paste("ntop =", ntop)) +
  scale_color_discrete(breaks = c("Immune:PDL1", "Immune:CTLA4", "Immune:CTR",
                                  "Islet:PDL1", "Islet:CTLA4", "Islet:CTR")) +
  theme(text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))

ggsave(paste0("pca_", ntop, ".png"),
       #plot = last_plot(),
       plot = pca,
       device = "png",
       path = paste0("./output/", ver, "/multi"),
       width = 5000, height = 2000, units = "px", dpi = 600,
       scale = 1.25
)
