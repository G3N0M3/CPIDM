#######
# Script for GO (ImPC)
#######


# GO ----------------------------------------------------------------------
bgGenes_ImPC <- row.names(annot_ImPC)
UP_ImPC <- row.names(annot_ImPC[annot_ImPC$padj < 0.1 &
                                  !is.na(annot_ImPC$padj) &
                                  annot_ImPC$log2FoldChange > 1, ]
)
allUP_ImPC <- factor(as.integer(bgGenes_ImPC %in% UP_ImPC))
names(allUP_ImPC) <- bgGenes_ImPC


onts <- c("MF", "BP", "CC")


# UP ----------------------------------------------------------------------
tab_UP_ImPC <- as.list(onts)
names(tab_UP_ImPC) <- onts
for (i in 1:length(onts)) {
  tgd <- new("topGOdata",
             ontology = onts[i],
             allGenes = allUP_ImPC,
             nodeSize = 5,
             annot = annFUN.org,
             mapping = "org.Mm.eg.db",
             ID = "ensembl")
  
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher")
  
  tab_UP_ImPC[[i]] <- GenTable(tgd,
                               Fisher.elim = resultTopGO.elim,
                               orderBy = "Fisher.classic",
                               topNodes = 200)
}
upGOres_ImPC <- rbind.fill(tab_UP_ImPC)


# Save --------------------------------------------------------------------
write.csv(upGOres_ImPC[order(upGOres_ImPC$Fisher.elim), ],
          file = paste0("./output/", ver, "/GO/UP_ImPC.csv"),
          row.names = FALSE)
