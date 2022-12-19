#######
# Script for GO (ImPctr)
#######


# GO ----------------------------------------------------------------------
bgGenes_ImPctr <- row.names(annot_ImPctr)
UP_ImPctr <- row.names(annot_ImPctr[annot_ImPctr$padj < 0.1 &
                                      !is.na(annot_ImPctr$padj) &
                                      annot_ImPctr$log2FoldChange > 1, ]
                       )
DOWN_ImPctr <- row.names(annot_ImPctr[annot_ImPctr$padj < 0.1 &
                                        !is.na(annot_ImPctr$padj) &
                                        annot_ImPctr$log2FoldChange < 1, ]
                         )
allUP_ImPctr <- factor(as.integer(bgGenes_ImPctr %in% UP_ImPctr))
allDOWN_ImPctr <- factor(as.integer(bgGenes_ImPctr %in% DOWN_ImPctr))
names(allUP_ImPctr) <- bgGenes_ImPctr
names(allDOWN_ImPctr) <- bgGenes_ImPctr


onts <- c("MF", "BP", "CC")


# UP ----------------------------------------------------------------------
tab_UP_ImPctr <- as.list(onts)
names(tab_UP_ImPctr) <- onts
for (i in 1:length(onts)) {
  tgd <- new("topGOdata",
             ontology = onts[i],
             allGenes = allUP_ImPctr,
             nodeSize = 5,
             annot = annFUN.org,
             mapping = "org.Mm.eg.db",
             ID = "ensembl")
  
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher")
  
  tab_UP_ImPctr[[i]] <- GenTable(tgd,
                                 Fisher.elim = resultTopGO.elim,
                                 orderBy = "Fisher.classic",
                                 topNodes = 200)
}
upGOres_ImPctr <- rbind.fill(tab_UP_ImPctr)


# DOWN --------------------------------------------------------------------
tab_DOWN_ImPctr <- as.list(onts)
names(tab_DOWN_ImPctr) <- onts
for (i in 1:length(onts)) {
  tgd <- new("topGOdata",
             ontology = onts[i],
             allGenes = allDOWN_ImPctr,
             nodeSize = 5,
             annot = annFUN.org,
             mapping = "org.Mm.eg.db",
             ID = "ensembl")
  
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher")
  
  tab_DOWN_ImPctr[[i]] <- GenTable(tgd,
                                   Fisher.elim = resultTopGO.elim,
                                   orderBy = "Fisher.classic",
                                   topNodes = 200)
}
downGOres_ImPctr <- rbind.fill(tab_DOWN_ImPctr)


# Save --------------------------------------------------------------------
write.csv(upGOres_ImPctr[order(upGOres_ImPctr$Fisher.elim), ],
          file = paste0("./output/", ver, "/GO/UP_ImPctr.csv"),
          row.names = FALSE)
write.csv(downGOres_ImPctr[order(downGOres_ImPctr$Fisher.elim), ],
          file = paste0("./output/", ver, "/GO/DOWN_ImPctr.csv"),
          row.names = FALSE)
