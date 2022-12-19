##### Immune_CTLA4 vs Immune_CTR #####



#### Subset sampleTable ####
table_ImCctr <- sampleTable[sampleTable$condition %in%
                              c("Immune_CTLA4", "Immune_CTR"),]



#### DESeqDataSet ####
dds_ImCctr <- DESeqDataSetFromHTSeqCount(sampleTable = table_ImCctr,
                                         directory = "./data/count/",
                                         design = ~ condition)



#### Pre-filtering and Factoring ####
## Pre-filtering
# removing rows that have only 0 or 1 read
dds_filt_ImCctr <- dds_ImCctr[rowSums(counts(dds_ImCctr)) > 1, ]
## Factoring
dds_filt_ImCctr$condition <- relevel(dds_filt_ImCctr$condition,
                                     ref = "Immune_CTR")



#### DESeq2 ####
DE_ImCctr <- DESeq(dds_filt_ImCctr)



#### Results ####
res_ImCctr <- results(DE_ImCctr)



#### MA-plot ####
res_rev_ImCctr <- res_ImCctr[order(res_ImCctr$padj, decreasing = TRUE), ]
png(filename = paste0("./output/", ver, "/MAplot/MA_ImCctr.png"),
    width = 10000, height = 10000, units = "px",
    bg = "white", res = 1500)
plotMA(res_rev_ImCctr,
       xlim = c(1e-01, 1e+07), ylim = c(-10, 10),
       alpha = 0.05, colSig = "red")
title("anti-PD-L1 treated vs Control (CD45+)")
abline(h=1, col = "grey40", lwd = 2, lty = 2)
abline(h=-1, col = "grey40", lwd = 2, lty = 2)
text(x = 1e6 * 3 , y = 2, "log2FC = 1")
text(x = 1e6 * 3 , y = -2, "log2FC = -1")
dev.off()



#### lfcShrink ####
shrink_ImCctr <- lfcShrink(DE_ImCctr, coef = 2, res = res_ImCctr)
shrink_rev_ImCctr <- shrink_ImCctr[order(shrink_ImCctr$padj, decreasing = TRUE), ]
png(filename = paste0("./output/", ver, "/MAplot/MA_shrink_ImCctr.png"),
    width = 10000, height = 10000, units = "px",
    bg = "white", res = 1500)
plotMA(shrink_rev_ImCctr,
       xlim = c(1e-01, 1e+07), ylim = c(-10, 10), alpha = 0.05,
       colSig = "red")
title("anti-CTLA-4 treated vs Control (CD45+)")
dev.off()
