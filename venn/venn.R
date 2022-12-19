#######
# Script that creates Venn diagram for Immune cell experiments
#######


# DEG list ----------------------------------------------------------------
# ImPC
cut_ImPC <- ord_annot_ImPC[!is.na(ord_annot_ImPC$padj) &
                             as.numeric(ord_annot_ImPC$padj) < 0.1 &
                             abs(as.numeric(ord_annot_ImPC$log2FoldChange)) > 1, ]
DEG_ImPC <- row.names(cut_ImPC)
write.csv(x = cut_ImPC[order(cut_ImPC$padj, decreasing = FALSE), ],
          file = paste0("./output/", ver, "/venn/DEG_ImPC.csv"))
# ImPctr
cut_ImPctr <- ord_annot_ImPctr[!is.na(ord_annot_ImPctr$padj) &
                                 as.numeric(ord_annot_ImPctr$padj) < 0.1 &
                                 abs(as.numeric(ord_annot_ImPctr$log2FoldChange)) > 1, ]
DEG_ImPctr <- row.names(cut_ImPctr)
write.csv(x = cut_ImPctr[order(cut_ImPC$padj, decreasing = FALSE), ],
          file = paste0("./output/", ver, "/venn/DEG_ImPctr.csv"))
# ImCctr
cut_ImCctr <- ord_annot_ImCctr[!is.na(ord_annot_ImCctr$padj) &
                                 as.numeric(ord_annot_ImCctr$padj) < 0.1 &
                                 abs(as.numeric(ord_annot_ImCctr$log2FoldChange)) > 1, ]
DEG_ImCctr <- row.names(cut_ImCctr)
write.csv(x = cut_ImCctr[order(cut_ImPC$padj, decreasing = FALSE), ],
          file = paste0("./output/", ver, "/venn/DEG_ImCctr.csv"))
# IsPC
cut_IsPC <- ord_annot_IsPC[!is.na(ord_annot_IsPC$padj) &
                             as.numeric(ord_annot_IsPC$padj) < 0.1 &
                             abs(as.numeric(ord_annot_IsPC$log2FoldChange)) > 1, ]
DEG_IsPC <- row.names(cut_IsPC)
write.csv(x = cut_IsPC[order(cut_IsPC$padj, decreasing = FALSE), ],
          file = paste0("./output/", ver, "/venn/DEG_IsPC.csv"))


# Plotting ----------------------------------------------------------------
venn_pal <- c("#d8e2f2", "#d6d4d4", "#fae4d5")
venn.diagram(
  
  # basic
  x = list(DEG_ImPC, DEG_ImPctr, DEG_ImCctr),
  category.names = c("anti-PDL1 vs anti-CTLA4",
                     "anti-CTLA4 vs CTR",
                     "anti-PDL1 vs CTR"),
  
  # output features
  output = TRUE,
  filename = paste0("./output/", ver, "/venn/venn.png"),
  imagetype = "png",
  height = 5000, width = 5000,
  resolution = 600,
  compression = "lzw",
  margin = 0.03,
  
  # circles
  lwd = 1,
  lty = 1, # 0: blank, 1: solid
  fill = venn_pal,
  
  # numbers
  cex = 1.5,
  fontface = "plain",
  fontfamily = "sans",
  print.mode = "raw",  # "raw", "percent", c("raw", "percent")
  #sigdigs = 2,
  
  # category names
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  
  # category layout
  rotation = 1,
  # diagram rotation
  rotation.degree = -90,  # positive rotates counter-clockwise
  
  ## category location (1: top-left, 2: top-right, 3: bottom)
  # category name location rotation in degrees, 0 -> 12h
  cat.pos = c(0, 180, -70),
  # category name location distance, 0 -> edge, minus to outer, plus to inner
  cat.dist = c(-0.07, -0.07, -0.09)
)
