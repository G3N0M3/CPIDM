####
# script for plotting volcano plot (ImPC)
####



# Modify Input ------------------------------------------------------------
## from DESeqResults object
# which is generated in each experiment script
vol_ImPC <- as.data.frame(res_ImPC[, c("log2FoldChange", "padj")])
# remove NA values in padj column
vol_ImPC <- vol_ImPC[!is.na(vol_ImPC$padj), ]
# columns for aes mapping purposes
vol_ImPC$alpha <- ifelse(vol_ImPC$padj < 0.1,
                         ifelse(abs(vol_ImPC$log2FoldChange) > 1,
                                1, .5),
                         .2)
vol_ImPC$color <- ifelse(vol_ImPC$padj < 0.1,
                         ifelse(abs(vol_ImPC$log2FoldChange) > 1,
                                ifelse(vol_ImPC$log2FoldChange > 1,
                                       "Sig (UP)",
                                       "Sig (DOWN)"),
                                "Sig (NOT)"),
                         "not-Sig")
vol_ImPC$color <- factor(vol_ImPC$color,
                         levels = c("Sig (UP)", "Sig (DOWN)",
                                    "Sig (NOT)", "not-Sig"))


# Plotting ----------------------------------------------------------------
volPlot_ImPC <-
  ggplot(vol_ImPC) +
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
  labs(title = "anti-PD-L1 treated vs anti-CTLA-4 (CD45+)",  # title
       color = "padj < 0.1"  # legend for color
       ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
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
ggsave("volcano_ImPC.png",
       plot = volPlot_ImPC,
       device = "png",
       path = paste0("./output/", ver, "/volcano"),
       width = 4000, height = 3500, units = "px", dpi = 600
)
