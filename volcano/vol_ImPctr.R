####
# script for plotting volcano plot (ImPctr)
####



# Modify Input ------------------------------------------------------------
# original data from DESeqResults object
# which is generated in each experiment script
vol_ImPctr <- as.data.frame(res_ImPctr[, c("log2FoldChange", "padj")])
# remove NA values in padj column
vol_ImPctr <- vol_ImPctr[!is.na(vol_ImPctr$padj), ]
# columns for aes mapping purposes
vol_ImPctr$alpha <- ifelse(vol_ImPctr$padj < 0.05,
                         ifelse(abs(vol_ImPctr$log2FoldChange) > 1,
                                1, .5),
                         .2)
vol_ImPctr$color <- ifelse(vol_ImPctr$padj < 0.05,
                         ifelse(abs(vol_ImPctr$log2FoldChange) > 1,
                                ifelse(vol_ImPctr$log2FoldChange > 1,
                                       "Sig (UP)",
                                       "Sig (DOWN)"),
                                "Sig (NOT)"),
                         "not-Sig")
vol_ImPctr$color <- factor(vol_ImPctr$color,
                           levels = c("Sig (UP)", "Sig (DOWN)",
                                      "Sig (NOT)", "not-Sig"))


# Plotting ----------------------------------------------------------------
volPlot_ImPctr <-
  ggplot(vol_ImPctr) +
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
  geom_hline(yintercept = -log10(0.05), linewidth = .75,
             color = "red2", linetype = "dotted") +
  geom_vline(xintercept = 1, linewidth = .75, color = "#bababa", linetype = "dotted") +
  geom_vline(xintercept = -1, linewidth = .75, color = "#bababa", linetype = "dotted") +
  geom_vline(xintercept = 0, linewidth = .75, color = "#babaab") +
  # legends
  guides(color = guide_legend(title = ),  # change legend name
         alpha = "none"  # remove legend for alpha
  ) +
  # plot title
  labs(title = "anti-PD-L1 treated vs Control (CD45+)",
       color = "padj < 0.05"  # legend for color
       ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  # coordination limits
  coord_cartesian(xlim = c(-12, 12), ylim = c(0, 10)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  # annotations
  annotate("text", label = "LFC = 1",
           x = 1, y = 9, hjust = -0.25,
           size = 3) +  # LFC=1
  annotate("text", label = "LFC = -1",
           x = -1, y = 9, hjust = 1.25,
           size = 3) +  # LFC=-1
  annotate("text", label = "padj = 0.05",
           x = 10, y = -log10(0.05), vjust = -1,
           size = 3)


# Save --------------------------------------------------------------------
ggsave("volcano_ImPctr.png",
       plot = volPlot_ImPctr,
       device = "png",
       path = paste0("./output/", ver, "/volcano"),
       width = 10000, height = 10000, units = "px", dpi = 1750
)
