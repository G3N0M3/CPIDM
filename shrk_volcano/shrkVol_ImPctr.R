####
# script for plotting volcano plot (ImPctr)
####



# Modify Input ------------------------------------------------------------
## from "lfcShrink"ed DESeqResults object
# generated in each experiment script
vol_shr_ImPctr <- as.data.frame(shrink_ImPctr[, c("log2FoldChange", "padj")])
# remove NA values in padj column
vol_shr_ImPctr <- vol_shr_ImPctr[!is.na(vol_shr_ImPctr$padj), ]
# columns for aes mapping purposes
vol_shr_ImPctr$alpha <- ifelse(vol_shr_ImPctr$padj < 0.05, 1, .2)
vol_shr_ImPctr$color <- ifelse(vol_shr_ImPctr$padj < 0.05, "Sig", "not-Sig")


# Plotting ----------------------------------------------------------------
volPlot_shr_ImPctr <-
  ggplot(vol_shr_ImPctr) +
  # basic
  geom_point(aes(x = log2FoldChange,
                 y = -log10(padj),
                 alpha = alpha,
                 color = color),
             stroke = NA) +
  scale_color_manual(values = c("black", "red")) +
  # background
  theme_classic() +
  # additional lines
  geom_hline(yintercept = -log10(0.05), linewidth = .75,
             color = "red3", linetype = "dotted") +
  geom_vline(xintercept = 1, linewidth = .75, color = "#bababa", linetype = "dotted") +
  geom_vline(xintercept = -1, linewidth = .75, color = "#bababa", linetype = "dotted") +
  geom_vline(xintercept = 0, linewidth = .75, color = "#bababa") +
  # legends
  guides(alpha = "none") +  # remove legend for alpha
  # labels
  labs(title = "anti-PD-L1 treated vs Control (CD45+)",  # title
       x = "log2FoldChange (shrinked)",  # x label
       color = "padj < 0.05"  # legend for color
       ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  # coordination limits
  coord_cartesian(xlim = c(-12, 12), ylim = c(0, 10)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))


# Save --------------------------------------------------------------------
ggsave("volcano_shr_ImPctr.png",
       plot = volPlot_shr_ImPctr,
       device = "png",
       path = paste0("./output/", ver, "/shrkVolcano"),
       width = 10000, height = 10000, units = "px", dpi = 1750
)
