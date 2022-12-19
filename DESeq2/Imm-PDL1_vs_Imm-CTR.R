##### Immune_PDL1 vs Immune_CTR #####



#### Subset sampleTable ####
table_ImPctr <- sampleTable[sampleTable$condition %in%
                              c("Immune_PDL1", "Immune_CTR"),]



#### DESeqDataSet ####
dds_ImPctr <- DESeqDataSetFromHTSeqCount(sampleTable = table_ImPctr,
                                         directory = "./data/count/",
                                         design = ~ condition)



#### Pre-filtering and Factoring ####
## Pre-filtering
# removing rows that have only 0 or 1 read
dds_filt_ImPctr <- dds_ImPctr[rowSums(counts(dds_ImPctr)) > 1, ]
## Factoring
dds_filt_ImPctr$condition <- relevel(dds_filt_ImPctr$condition,
                                     ref = "Immune_CTR")



#### DESeq2 ####
DE_ImPctr <- DESeq(dds_filt_ImPctr)



#### Results ####
res_ImPctr <- results(DE_ImPctr)



#### MA-plot ####
res_rev_ImPctr <- res_ImPctr[order(res_ImPctr$padj, decreasing = TRUE), ]
png(filename = paste0("./output/", ver, "/MAplot/MA_ImPctr.png"),
    width = 10000, height = 10000, units = "px",
    bg = "white", res = 1500)
plotMA(res_rev_ImPctr,
       xlim = c(1e-01, 1e+07), ylim = c(-10, 10),
       alpha = 0.05, colSig = "red")
title("anti-PD-L1 treated vs Control (CD45+)")
abline(h=1, col = "grey40", lwd = 2, lty = 2)
abline(h=-1, col = "grey40", lwd = 2, lty = 2)
text(x = 1e6 * 3 , y = 2, "log2FC = 1")
text(x = 1e6 * 3 , y = -2, "log2FC = -1")
dev.off()



#### lfcShrink ####
shrink_ImPctr <- lfcShrink(DE_ImPctr, coef = 2, res = res_ImPctr)
# applying reverse ordered by padj 
shrink_rev_ImPctr <- shrink_ImPctr[order(shrink_ImPctr$padj, decreasing = TRUE), ]
png(filename = paste0("./output/", ver, "/MAplot/MA_shrink_ImPctr.png"),
    width = 10000, height = 10000, units = "px",
    bg = "white", res = 1500)
plotMA(shrink_rev_ImPctr,
       xlim = c(1e-01, 1e+07), ylim = c(-10, 10), alpha = 0.05,
       colSig = "red")
title("anti-PD-L1 treated vs Control (CD45+)")
dev.off()