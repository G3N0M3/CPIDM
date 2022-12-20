#######
# Plotting with DEGs for Immune cells
#######


# Input Data --------------------------------------------------------------
# DEG list from venn.R
rld_DEG_Im <- rld_Im[row.names(rld_annot_Im) %in% DEG_ImPC, ]


# PCA ---------------------------------------------------------------------
degPCA_Im <- plotPCA(rld_DEG_Im,
                  intgroup = c("treatment"))

data_degPCA_Im <- plotPCA(rld_DEG_Im, intgroup = c("treatment"),
                          returnData = TRUE)

degPCA_Im <- degPCA_Im +
  theme_light() +
  labs(color = "Treatment") +
  scale_color_discrete(breaks = c("PDL1", "CTLA4", "CTR"))

ggsave("degPCA_Im.png",
       #plot = last_plot(),
       plot = degPCA_Im,
       device = "png",
       path = paste0("./output/", ver, "/Immune"),
       width = 2500, height = 1000, units = "px", dpi = 600,
       scale = 1.25
)


# Heatmap -----------------------------------------------------------------
rld_annot_DEG_Im <- rld_annot_Im[row.names(rld_annot_Im) %in% DEG_ImPC, ]

png(filename = paste0("./output/", ver, "/Immune/DEG_heatmap_Im.png"),
    width = 5000, height = 7000, units = "px",
    bg = "white", res = 600)
pheatmap(assay(rld_annot_DEG_Im),
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, fontsize = 20,
         scale = "row",
         treeheight_row = 0, treeheight_col = 0,
         col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
dev.off()


# Sample ------------------------------------------------------------------
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(25)
# calculate sample distance
sampleDists <- dist(t(assay(rld_annot_DEG_Im)))
sampleDistMatrix <- as.matrix(sampleDists)
rld_annot_DEG_Im$samples <- c("ImCTR_1", "ImCTR_2", "ImCTR_3", "ImCTR_4",
                        "ImCTLA4_1", "ImCTLA4_2", "ImCTLA4_3", "ImCTLA4_4",
                        "ImPDL1_1", "ImPDL1_2", "ImPDL1_3")
rownames(sampleDistMatrix) <- rld_annot_DEG_Im$samples
colnames(sampleDistMatrix) <- NULL

## save plots
# clustered
png(filename = paste0("./output/", ver, "/Immune/DEG_sample_Im.png"),
    width = 6000, height = 5000, units = "px",
    bg = "white", res = 600)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, fontsize = 15)
dev.off()
# not-clustered
png(filename = paste0("./output/", ver, "/Immune/DEG_sample_unclust_Im.png"),
    width = 6000, height = 5000, units = "px",
    bg = "white", res = 600)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, fontsize = 15,
         cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()


# Paper DEGs --------------------------------------------------------------
paperDEGs_Im <- c("Cd8a", "Cd160", "Serna4a", "Fasl", "Oasl1", "Klrk1",
                  "Bst1", "Cxcl9", "Ly6a", "Il18bp", "Isg15", "Irf7",
                  "Cd40", "Tnfrs1a", "Il27", "Gzmb", "Ccl7", "Xaf1",
                  "Ccl5", "Gzma", "Pdcd1", "Nrp", "Btnl2", "Pf4", "Arid2",
                  "Tnfrsf26", "Sarm1", "Homer2", "Rbpjl", "H2-Ob")


# Heatmap with paperDEGs --------------------------------------------------
rld_annot_paperDEG_Im <- rld_annot_Im[row.names(rld_annot_Im) %in% paperDEGs_Im, ]

png(filename = paste0("./output/", ver, "/Immune/paperDEGs_heatmap_Im.png"),
    width = 5000, height = 7000, units = "px",
    bg = "white", res = 600)
pheatmap(assay(rld_annot_paperDEG_Im),
         cluster_rows = TRUE, cluster_cols = FALSE,
         show_rownames = TRUE, fontsize = 20,
         scale = "row",
         treeheight_row = 0,
         col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(25))
dev.off()
