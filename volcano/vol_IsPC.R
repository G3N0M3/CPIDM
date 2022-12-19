####
# script for plotting volcano plot (IsPC)
####



# Modify Input ------------------------------------------------------------
# original data from DESeqResults object
# which is generated in each experiment script
vol_IsPC <- as.data.frame(res_IsPC[, c("log2FoldChange", "padj")])
# remove NA values in padj column
vol_IsPC <- vol_IsPC[!is.na(vol_IsPC$padj), ]
# columns for aes mapping purposes
vol_IsPC$alpha <- ifelse(vol_IsPC$padj < 0.1,
                         ifelse(abs(vol_IsPC$log2FoldChange) > 1,
                                1, .5),
                         .2)
vol_IsPC$color <- ifelse(vol_IsPC$padj < 0.1,
                         ifelse(abs(vol_IsPC$log2FoldChange) > 1,
                                ifelse(vol_IsPC$log2FoldChange > 1,
                                       "Sig (UP)",
                                       "Sig (DOWN)"),
                                "Sig (NOT)"),
                         "not-Sig")
vol_IsPC$color <- factor(vol_IsPC$color,
                         levels = c("Sig (UP)", "Sig (DOWN)",
                                    "Sig (NOT)", "not-Sig"))


# Plotting ----------------------------------------------------------------
volPlot_IsPC <-
  ggplot(vol_IsPC) +
  # basic
  geom_point(aes(x = log2FoldChange,
                 y = -log10(padj),
                 alpha = alpha,
                 color = color),
             stroke = NA) +
  scale_color_manual(values = c("red", "blue", "black", "black")) +
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
  # plot title
  labs(title = "anti-PD-L1 treated vs anti-CTLA-4 treated (CD45-)",  # title
       color = "padj < 0.1"  # legend for color
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
  annotate("text", label = "padj = 0.1",
           x = 10, y = -log10(0.1), vjust = -1,
           size = 3)


# Save --------------------------------------------------------------------
ggsave("volcano_IsPC.png",
       plot = volPlot_IsPC,
       device = "png",
       path = paste0("./output/", ver, "/volcano"),
       width = 4000, height = 3000, units = "px", dpi = 600
)
