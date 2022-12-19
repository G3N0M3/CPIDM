####
# script for plotting shrinked volcano plot (ImPC)
####



# Modify Input ------------------------------------------------------------
## from "lfcShrink"ed DESeqResults object
# generated in each experiment script
vol_shr_ImPC <- as.data.frame(shrink_ImPC[, c("log2FoldChange", "padj")])
# remove NA values in padj column
vol_shr_ImPC <- vol_shr_ImPC[!is.na(vol_shr_ImPC$padj), ]
# columns for aes mapping purposes
vol_shr_ImPC$alpha <- ifelse(vol_shr_ImPC$padj < 0.05, 1, .2)
vol_shr_ImPC$color <- ifelse(vol_shr_ImPC$padj < 0.05, "Sig", "not-Sig")


# Plotting ----------------------------------------------------------------
volPlot_shr_ImPC <-
  ggplot(vol_shr_ImPC) +
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
  labs(title = "anti-PD-L1 treated vs anti-CTLA-4 (CD45+)",  # title
       x = "log2FoldChange (shrinked)",  # x label
       color = "padj < 0.05"  # legend for color
       ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  # coordination limits
  coord_cartesian(xlim = c(-12, 12), ylim = c(0, 10)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))


# Save --------------------------------------------------------------------
ggsave("volcano_shr_ImPC.png",
       plot = volPlot_shr_ImPC,
       device = "png",
       path = paste0("./output/", ver, "/shrkVolcano"),
       width = 4000, height = 3000, units = "px", dpi = 600
)
