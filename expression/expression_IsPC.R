#######
# Creates boxplots of DEG in ImPC
#######



# Import ------------------------------------------------------------------
library(ggpubr)

# Input Data --------------------------------------------------------------
# expression values
count_Is <- as.data.frame(counts(DE_Is))
# Annotation
row.names(count_Is) <- sapply(row.names(count_Is),
                              FUN = function(x) {strsplit(x, "\\.")[[1]][1]},
                              USE.NAMES = FALSE)
count_Is$symbol <- mapIds(org.Mm.eg.db,  # AnnotationDb object
                          keys = row.names(count_Is),  # keys to select from db
                          column = "SYMBOL",  # the column to search on (match)
                          keytype = "ENSEMBL",  # annotation type of keys
                          multiVals = "first")  # action upon multiple mapping
count_Is$symbol[is.na(count_Is$symbol)] <- row.names(count_Is)[is.na(count_Is$symbol)]
# for t.test
comparisons <- list(c("anti-CTLA-4", "anti-PD-L1"),
                    c("anti-PD-L1", "Control"),
                    c("anti-CTLA-4", "Control"))


# Plotting ----------------------------------------------------------------
# genes of interest
gene_list <- c("Cd274", "Cxcl10", "Irf1")
for (gene in gene_list) {
  # data.frame of input
  df_Is <-
    data.frame(run = colnames(count_Is)[1:length(colnames(count_Is))-1],
               treatment = rep(c("Control", "anti-CTLA-4", "anti-PD-L1"),
                               times = c(4, 4, 3)),
               count = as.numeric(count_Is[count_Is$symbol == gene, ][1:ncol(count_Is)-1])
               )
  # plotting
  boxplot_Is <-
    ggplot(df_Is,
           aes(x = as.factor(treatment),
               y = count)) +
    geom_boxplot(fill = "#ECECEC") +
    geom_point(aes(color = treatment), size = 5, position = "jitter") +
    scale_color_manual(values = c("#4971C3", "#DC3912", "#FF9900")) +
    theme_light() +
    labs(x = "Treatment", y = "Count",
         title = gene,
         color = "Treatment") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 12)) +
    # t.test
    stat_compare_means(comparisons = comparisons)
  
  
  # Saving ------------------------------------------------------------------
  ggsave(paste0("Is_boxplot_", gene, ".png"),
         #plot = last_plot(),
         plot = boxplot_Is,
         device = "png",
         path = paste0("./output/", ver, "/expression"),
         width = 2500, height = 2000, units = "px", dpi = 600,
         scale = 1.25)
}

