##### Islet_PDL1 vs Islet_CTLA4 #####



#### Subset sampleTable ####
table_IsPC <- sampleTable[sampleTable$condition %in%
                            c("Islet_PDL1", "Islet_CTLA4"),]



#### DESeqDataSet ####
dds_IsPC <- DESeqDataSetFromHTSeqCount(sampleTable = table_IsPC,
                                       directory = "./data/count/",
                                       design = ~ condition)



#### Pre-filtering and Factoring ####
## Pre-filtering
# removing rows that have only 0 or 1 read
dds_filt_IsPC <- dds_IsPC[rowSums(counts(dds_IsPC)) > 1, ]
## Factoring
dds_filt_IsPC$condition <- relevel(dds_filt_IsPC$condition,
                                   ref = "Islet_CTLA4")



#### DESeq2 ####
DE_IsPC <- DESeq(dds_filt_IsPC)



#### Results ####
res_IsPC <- results(DE_IsPC)



#### MA-plot ####
res_rev_IsPC <- res_IsPC[order(res_IsPC$padj, decreasing = TRUE), ]
png(filename = paste0("./output/", ver, "/MAplot/MA_IsPC.png"),
    width = 4000, height = 4000, units = "px",
    bg = "white", res = 600)
plotMA(res_rev_IsPC,
       xlim = c(1e-01, 1e+07), ylim = c(-10, 10),
       alpha = 0.1, colSig = "red")
title("anti-PD-L1 treated vs anti-CTLA-4 treated (CD45-)")
abline(h=1, col = "grey40", lwd = 2, lty = 2)
abline(h=-1, col = "grey40", lwd = 2, lty = 2)
text(x = 1e6 * 3 , y = 2, "log2FC = 1")
text(x = 1e6 * 3 , y = -2, "log2FC = -1")
dev.off()



#### lfcShrink ####
shrink_IsPC <- lfcShrink(DE_IsPC, coef = 2, res = res_IsPC)
shrink_rev_IsPC <- shrink_IsPC[order(shrink_IsPC$padj, decreasing = TRUE), ]
png(filename = paste0("./output/", ver, "/MAplot/MA_shrink_IsPC.png"),
    width = 4000, height = 4000, units = "px",
    bg = "white", res = 600)
plotMA(shrink_rev_IsPC,
       xlim = c(1e-01, 1e+07), ylim = c(-10, 10),
       alpha = 0.1, colSig = "red")
title("anti-PD-L1 treated vs anti-CTLA-4 treated (CD45-)")
dev.off()

