##### Immune_PDL1 vs Immune_CTLA4 #####



#### Subset sampleTable ####
table_ImPC <- sampleTable[sampleTable$condition %in%
                              c("Immune_PDL1", "Immune_CTLA4"), ]



#### DESeqDataSet ####
dds_ImPC <- DESeqDataSetFromHTSeqCount(sampleTable = table_ImPC,
                                         directory = "./data/count/",
                                         design = ~ condition)



#### Pre-filtering and Factoring ####
## Pre-filtering
# removing rows that have only 0 or 1 read
dds_filt_ImPC <- dds_ImPC[rowSums(counts(dds_ImPC)) > 1, ]
## Factoring
dds_filt_ImPC$condition <- relevel(dds_filt_ImPC$condition, ref = "Immune_CTLA4")



#### DESeq2 ####
DE_ImPC <- DESeq(dds_filt_ImPC)



#### Results ####
res_ImPC <- results(DE_ImPC)



#### MA-plot ####
res_rev_ImPC <- res_ImPC[order(res_ImPC$padj, decreasing = TRUE), ]
png(filename = paste0("./output/", ver, "/MAplot/MA_ImPC.png"),
    width = 4000, height = 4000, units = "px",
    bg = "white", res = 600)
plotMA(res_rev_ImPC,
       xlim = c(1e-01, 1e+07), ylim = c(-10, 10),
       alpha = 0.1, colSig = "red")
title("anti-PD-L1 treated vs anti-CTLA-4 treated (CD45+)")
abline(h=1, col = "grey40", lwd = 2, lty = 2)
abline(h=-1, col = "grey40", lwd = 2, lty = 2)
text(x = 1e6 * 3 , y = 2, "log2FC = 1")
text(x = 1e6 * 3 , y = -2, "log2FC = -1")
dev.off()



#### lfcShrink ####
shrink_ImPC <- lfcShrink(DE_ImPC, coef = 2, res = res_ImPC)
shrink_rev_ImPC <- shrink_ImPC[order(shrink_ImPC$padj, decreasing = TRUE), ]
png(filename = paste0("./output/", ver, "/MAplot/MA_shrink_ImPC.png"),
    width = 4000, height = 4000, units = "px",
    bg = "white", res = 600)
plotMA(shrink_rev_ImPC,
       xlim = c(1e-01, 1e+07), ylim = c(-10, 10), alpha = 0.1,
       colSig = "red")
title("anti-PD-L1 treated vs anti-CTLA-4 treated (CD45+)")
dev.off()

