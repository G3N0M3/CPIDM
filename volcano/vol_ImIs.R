####
# script for plotting volcano plot (ImIs)
####



# Modify Input ------------------------------------------------------------
## from DESeqResults object
# which is generated in each experiment script
vol <- as.data.frame(res[, c("log2FoldChange", "padj")])
# remove NA values in padj column
vol <- vol[!is.na(vol$padj), ]
# columns for aes mapping purposes
vol$alpha <- ifelse(vol$padj < 0.1,
                         ifelse(abs(vol$log2FoldChange) > 1,
                                1, .5),
                         .2)
vol$color <- ifelse(vol$padj < 0.1,
                         ifelse(abs(vol$log2FoldChange) > 1,
                                ifelse(vol$log2FoldChange > 1,
                                       "Sig (UP)",
                                       "Sig (DOWN)"),
                                "Sig (NOT)"),
                         "not-Sig")
vol$color <- factor(vol$color,
                         levels = c("Sig (UP)", "Sig (DOWN)",
                                    "Sig (NOT)", "not-Sig"))


# Plotting ----------------------------------------------------------------
volPlot <-
  ggplot(vol) +
  # basic
  geom_point(aes(x = log2FoldChange,
                 y = -log10(padj),
                 alpha = alpha,
                 color = color),
             stroke = NA) +
  scale_color_manual(values = c("red", "blue", "green4", "black")) +
  # background
  theme_classic() +
  # additional lines
  geom_hline(yintercept = -log10(0.1), linewidth = .75,
             color = "red3", linetype = "dotted") +
  geom_vline(xintercept = 1, linewidth = .75, color = "#bababa", linetype = "dotted") +
  geom_vline(xintercept = -1, linewidth = .75, color = "#bababa", linetype = "dotted") +
  geom_vline(xintercept = 0, linewidth = .75, color = "#bababa") +
  # legends
  guides(alpha = "none") +  # remove legend for alpha
  # labels
  labs(title = "Immune (CD45+) vs Islet (Cd45-)",  # title
       color = "padj < 0.1"  # legend for color
       ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  # coordination limits
  coord_cartesian(xlim = c(-12, 12), ylim = c(0, NA)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(.01, 0)) +
  # annotations
  annotate("text", label = "LFC = 1",
           x = 1, y = 52, hjust = -0.25,
           size = 3) +  # LFC=1
  annotate("text", label = "LFC = -1",
           x = -1, y = 52, hjust = 1.25,
           size = 3) +  # LFC=-1
  annotate("text", label = "padj = 0.1",
           x = 10, y = -log10(0.1), vjust = -1,
           size = 3)


# Save --------------------------------------------------------------------
ggsave("volcano_ImIs.png",
       plot = volPlot,
       device = "png",
       path = paste0("./output/", ver, "/volcano"),
       width = 4000, height = 3500, units = "px", dpi = 600
)
