#######
# Creates boxplots of DEG in ImPC
#######


# Input Data --------------------------------------------------------------
# expression values
count_Im <- as.data.frame(counts(DE_Im))
# Annotation
row.names(count_Im) <- sapply(row.names(count_Im),
                              FUN = function(x) {strsplit(x, "\\.")[[1]][1]},
                              USE.NAMES = FALSE)
count_Im$symbol <- mapIds(org.Mm.eg.db,  # AnnotationDb object
                          keys = row.names(count_Im),  # keys to select from db
                          column = "SYMBOL",  # the column to search on (match)
                          keytype = "ENSEMBL",  # annotation type of keys
                          multiVals = "first")  # action upon multiple mapping
count_Im$symbol[is.na(count_Im$symbol)] <- row.names(count_Im)[is.na(count_Im$symbol)]
# for t.test
comparisons <- list(c("anti-CTLA-4", "anti-PD-L1"),
                    c("anti-PD-L1", "Control"),
                    c("anti-CTLA-4", "Control"))

#gene
#compare_means(count ~ treatment, data = df_Im)


# Plotting ----------------------------------------------------------------
# genes of interest
gene_list <- c("Cd8a", "Cd274", "Cxcl9", "Cxcl10", "Fasl", "Gzma", "Gzmb",
               "Ifng", "Irf1", "Pdcd1")
for (gene in gene_list) {
  
  # data.frame of input
  df_Im <-
    data.frame(run = colnames(count_Im)[1:length(colnames(count_Im))-1],
               treatment = rep(c("Control", "anti-CTLA-4", "anti-PD-L1"),
                               times = c(4, 4, 3)),
               count = as.numeric(count_Im[count_Im$symbol == gene, ][1:ncol(count_Im)-1])
               )
  # plotting
  boxplot_Im <-
    ggplot(df_Im,
           aes(x = as.factor(treatment),
               y = count)) +
    geom_boxplot(fill = "#ECECEC") +
    geom_point(aes(color = treatment), size = 5) +
    scale_color_manual(values = c("#4971C3", "#DC3912", "#FF9900")) +
    theme_light() +
    labs(x = "Treatment", y = "Count",
         title = gene,
         color = "Treatment") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 12)) +
    # t.test
    stat_compare_means(comparisons = comparisons,
                       p.adjust.methods = "fdr",
                       aes(label = format.pval(..p.adj.., digits = 3)))
  gene
  compare_means(count ~ treatment, data = df_Im)
  # Saving
  ggsave(paste0("Im_boxplot_", gene, ".png"),
         #plot = last_plot(),
         plot = boxplot_Im,
         device = "png",
         path = paste0("./output/", ver, "/expression"),
         width = 2500, height = 2000, units = "px", dpi = 600,
         scale = 1.25)
}
