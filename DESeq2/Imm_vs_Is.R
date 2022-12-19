#######
# Immune vs Islet, based on CD45
#######


# sampleTable -------------------------------------------------------------
table_ImIs <- sampleTable
condition <- sapply(table_ImIs$condition,
                    FUN = function(x) {strsplit(x, "_")[[1]][1]},
                    USE.NAMES = FALSE)
table_ImIs$condition <- as.factor(condition)


# DESeqDataSet ------------------------------------------------------------
dds_ImIs <- DESeqDataSetFromHTSeqCount(sampleTable = table_ImIs,
                                       directory = "./data/count/",
                                       design = ~ condition)


# Pre-processing ----------------------------------------------------------
# removing rows that have only 0 or 1 read
dds_filt_ImIs <- dds_ImIs[rowSums(counts(dds_ImIs)) > 1, ]

# Factoring
dds_filt_ImIs$condition <- relevel(dds_filt_ImIs$condition, ref = "Islet")


# DESeq2 ------------------------------------------------------------------
DE_ImIs <- DESeq(dds_filt_ImIs)



# Results -----------------------------------------------------------------
res_ImIs <- results(DE_ImIs)


# MA-plot -----------------------------------------------------------------
png(filename = paste0("./output/", ver, "/MAplot/MA_ImIs.png"),
    width = 10000, height = 10000, units = "px",
    bg = "white", res = 1500)
plotMA(res_ImIs[order(res_ImIs$padj, decreasing = TRUE), ],
       xlim = c(1e-01, 1e+07), ylim = c(-10, 10),
       alpha = 0.05, colSig = "red")
title("Immune (CD45+) vs Islet (CD45-)")
abline(h=1, col = "grey40", lwd = 2, lty = 2)
abline(h=-1, col = "grey40", lwd = 2, lty = 2)
text(x = 1e6 * 3 , y = 2, "log2FC = 1")
text(x = 1e6 * 3 , y = -2, "log2FC = -1")
dev.off()


# lfcShrink ---------------------------------------------------------------
shrink_ImIs <- lfcShrink(DE_ImIs, coef = 2, res = res_ImIs)
png(filename = paste0("./output/", ver, "/MAplot/MA_shrink_ImIs.png"),
    width = 10000, height = 10000, units = "px",
    bg = "white", res = 1500)
plotMA(shrink_ImIs[order(shrink_ImIs$padj, decreasing = TRUE), ],
       xlim = c(1e-01, 1e+07), ylim = c(-10, 10), alpha = 0.05,
       colSig = "red")
title("Immune (CD45+) vs Islet (CD45-)")
dev.off()
