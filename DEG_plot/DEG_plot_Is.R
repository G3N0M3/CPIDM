#######
# Plotting with DEGs for Islet cells
#######


# Input Data --------------------------------------------------------------
# DEG list from venn.R
rld_DEG_Is <- rld_Is[row.names(rld_annot_Is) %in% DEG_IsPC, ]

degPCA_Is <- plotPCA(rld_DEG_Is,
                     intgroup = c("treatment"))

degPCA_Is <- degPCA_Is +
  theme_light() +
  labs(color = "Treatment") +
  scale_color_discrete(breaks = c("PDL1", "CTLA4", "CTR"))

ggsave("degPCA_Is.png",
       #plot = last_plot(),
       plot = degPCA_Is,
       device = "png",
       path = paste0("./output/", ver, "/Islet"),
       width = 2500, height = 1000, units = "px", dpi = 600,
       scale = 1.25
)


# Heatmap -----------------------------------------------------------------
rld_annot_DEG_Is <- rld_annot_Is[row.names(rld_annot_Is) %in% DEG_IsPC, ]

png(filename = paste0("./output/", ver, "/Islet/DEG_heatmap_Is.png"),
    width = 5000, height = 7000, units = "px",
    bg = "white", res = 600)
pheatmap(assay(rld_annot_DEG_Is),
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = FALSE, fontsize = 20)
dev.off()


# Sample ------------------------------------------------------------------
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(25)
# calculate sample distance
sampleDists <- dist(t(assay(rld_annot_DEG_Is)))
sampleDistMatrix <- as.matrix(sampleDists)
rld_annot_DEG_Is$samples <- c("IsCTR_1", "IsCTR_2", "IsCTR_3", "IsCTR_4",
                              "IsCTLA4_1", "IsCTLA4_2", "IsCTLA4_3", "IsCTLA4_4",
                              "IsPDL1_1", "IsPDL1_2", "IsPDL1_3")
rownames(sampleDistMatrix) <- rld_annot_DEG_Is$samples
colnames(sampleDistMatrix) <- NULL

## save plots
# clustered
png(filename = paste0("./output/", ver, "/Islet/DEG_sample_Is.png"),
    width = 6000, height = 5000, units = "px",
    bg = "white", res = 600)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, fontsize = 15)
dev.off()
# not-clustered
png(filename = paste0("./output/", ver, "/Islet/DEG_sample_unclust_Is.png"),
    width = 6000, height = 5000, units = "px",
    bg = "white", res = 600)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, fontsize = 15,
         cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()


# Paper DEGs --------------------------------------------------------------
paperDEG_Is <- c("Slc2a2", "Fam83d", "Lin7c", "Gm37013", "Mras", "Cenpj", "Mlkl",
                 "Sprr1a", "Maff", "Isg20", "Ifit1", "Usp18", "Gbp4",
                 "Irf7", "Ly6e", "Cxcl10", "Ifit3", "Gbp2", "Oasl2", "Oas1a",
                 "Ifi44", "Rsad2", "Znfx1", "Fcgbp", "Ddit3", "Ifit2",
                 "Batf2", "Irf1", "Cd274", "Gbp7", "Ankrd13a")


# Heatmap with paperDEGs --------------------------------------------------
rld_annot_paperDEG_Is <- rld_annot_Is[row.names(rld_annot_Is) %in% paperDEG_Is, ]

png(filename = paste0("./output/", ver, "/Islet/paperDEGs_heatmap_Is.png"),
    width = 5000, height = 7000, units = "px",
    bg = "white", res = 600)
pheatmap(assay(rld_annot_DEG_Is),
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, fontsize = 20,
         scale = "row",
         col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(25))
dev.off()


# Sample ------------------------------------------------------------------
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(25)
# calculate sample distance
sampleDists <- dist(t(assay(rld_annot_DEG_Is)))
sampleDistMatrix <- as.matrix(sampleDists)
rld_annot_DEG_Is$samples <- c("IsCTR_1", "IsCTR_2", "IsCTR_3", "IsCTR_4",
                              "IsCTLA4_1", "IsCTLA4_2", "IsCTLA4_3", "IsCTLA4_4",
                              "IsPDL1_1", "IsPDL1_2", "IsPDL1_3")
rownames(sampleDistMatrix) <- rld_annot_DEG_Is$samples
colnames(sampleDistMatrix) <- NULL

## save plots
# clustered
png(filename = paste0("./output/", ver, "/Islet/paperDEGs_sample_Is.png"),
    width = 6000, height = 5000, units = "px",
    bg = "white", res = 600)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, fontsize = 15)
dev.off()
# not-clustered
png(filename = paste0("./output/", ver, "/Islet/paperDEGs_sample_unclust_Is.png"),
    width = 6000, height = 5000, units = "px",
    bg = "white", res = 600)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, fontsize = 15,
         cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()


# PCA ---------------------------------------------------------------------
rld_DEG_Is <- rld_Is[row.names(rld_annot_Is) %in% paperDEG_Is, ]

degPCA_Is <- plotPCA(rld_DEG_Is,
                     intgroup = c("treatment"))

degPCA_Is <- degPCA_Is +
  theme_light() +
  labs(color = "Treatment") +
  scale_color_discrete(breaks = c("PDL1", "CTLA4", "CTR"))

ggsave("paper_degPCA_Is.png",
       #plot = last_plot(),
       plot = degPCA_Is,
       device = "png",
       path = paste0("./output/", ver, "/Islet"),
       width = 3000, height = 1000, units = "px", dpi = 600,
       scale = 1.25
)

