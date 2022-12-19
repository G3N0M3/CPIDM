#######
# Pctr_PC GO
#######


# DEG ---------------------------------------------------------------------
# ImPC
go_deg_ImPC <- res_ImPC[!is.na(res_ImPC$padj) &
                          res_ImPC$padj < 0.1 &
                          abs(res_ImPC$log2FoldChange) > 1, ]
# ImPctr
go_deg_ImPctr <- res_ImPctr[!is.na(res_ImPctr$padj) &
                              res_ImPctr$padj < 0.1 &
                              abs(res_ImPctr$log2FoldChange) > 1, ]
# ImCctr
go_deg_ImCctr <- res_ImCctr[!is.na(res_ImCctr$padj) &
                              res_ImCctr$padj < 0.1 &
                              abs(res_ImCctr$log2FoldChange) > 1, ]
# DEG list
deg_list_ImPC <- row.names(go_deg_ImPC)
deg_list_ImPctr <- row.names(go_deg_ImPctr)
deg_list_ImCctr <- row.names(go_deg_ImCctr)


# Targeted Location -------------------------------------------------------
tar_loc <- intersect(deg_list_ImPctr, deg_list_ImPC)
tar_loc <- setdiff(tar_loc, deg_list_ImCctr)
length(tar_loc)  # has to be 116


ref_tar_loc <- sapply(tar_loc,
                      FUN = function(x) {strsplit(x, "\\.")[[1]][1]},
                      USE.NAMES = FALSE)


# GO ----------------------------------------------------------------------
bgGenes <- sapply(row.names(res_ImPC),
                  FUN = function(x) {strsplit(x, "\\.")[[1]][1]},
                  USE.NAMES = FALSE)
#UP_ImPC <- ref_tar_loc
allUP <- factor(as.integer(bgGenes %in% ref_tar_loc))
names(allUP) <- bgGenes


onts <- c("MF", "BP", "CC")


# UP ----------------------------------------------------------------------
tab_UP <- as.list(onts)
names(tab_UP) <- onts
for (i in 1:length(onts)) {
  tgd <- new("topGOdata",
             ontology = onts[i],
             allGenes = allUP,
             nodeSize = 5,
             annot = annFUN.org,
             mapping = "org.Mm.eg.db",
             ID = "ensembl")
  
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher")
  
  tab_UP[[i]] <- GenTable(tgd,
                               Fisher.elim = resultTopGO.elim,
                               orderBy = "Fisher.classic",
                               topNodes = 200)
}
upGOres <- rbind.fill(tab_UP)


# Save --------------------------------------------------------------------
write.csv(upGOres[order(as.numeric(upGOres$Fisher.elim)), ],
          file = paste0("./output/", ver, "/GO/UP.csv"),
          row.names = FALSE)


# Input data --------------------------------------------------------------
go_up <- upGOres[order(as.numeric(upGOres$Fisher.elim)), ]

go_up$pval <- -log10(as.numeric(go_up$Fisher.elim))
go_up <- go_up[order(go_up$pval, decreasing = TRUE), ]
go_up <- go_up[1:10, ]


# Plot --------------------------------------------------------------------
goBar <- ggplot(go_up) +
  geom_bar(aes(x = pval, y = GO.ID),
           stat = 'identity') +
  scale_y_discrete(limits = rev(go_up$GO.ID)) +
  theme_light() +
  coord_cartesian(xlim = c(0, 11)) +
  scale_x_continuous(expand = c(0, 0)) +
  geom_text(aes(x = .15, y = GO.ID, label = Term),
            color = "white", hjust = 0, fontface = "bold") +
  labs(x = "-Log10(p value)", y = "") +
  theme(axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12))

goBar

ggsave("goBarplot.png",
       plot = goBar,
       device = "png",
       path = paste0("./output/", ver, "/GO_barplot"),
       width = 5000, height = 3000, units = "px", dpi = 600
)


