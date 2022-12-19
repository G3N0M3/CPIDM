#######
# Script for GO (ImCctr)
#######


# GO ----------------------------------------------------------------------
bgGenes_ImCctr <- row.names(annot_ImCctr)
UP_ImCctr <- row.names(annot_ImCctr[annot_ImCctr$padj < 0.1 &
                                      !is.na(annot_ImCctr$padj) &
                                      annot_ImCctr$log2FoldChange > 1, ]
                       )
DOWN_ImCctr <- row.names(annot_ImCctr[annot_ImCctr$padj < 0.1 &
                                        !is.na(annot_ImCctr$padj) &
                                        annot_ImCctr$log2FoldChange < 1, ]
                         )
allUP_ImCctr <- factor(as.integer(bgGenes_ImCctr %in% UP_ImCctr))
allDOWN_ImCctr <- factor(as.integer(bgGenes_ImCctr %in% DOWN_ImCctr))
names(allUP_ImCctr) <- bgGenes_ImCctr
names(allDOWN_ImCctr) <- bgGenes_ImCctr


onts <- c("MF", "BP", "CC")


# UP ----------------------------------------------------------------------
tab_UP_ImCctr <- as.list(onts)
names(tab_UP_ImCctr) <- onts
for (i in 1:length(onts)) {
  tgd <- new("topGOdata",
             ontology = onts[i],
             allGenes = allUP_ImCctr,
             nodeSize = 5,
             annot = annFUN.org,
             mapping = "org.Mm.eg.db",
             ID = "ensembl")
  
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher")
  
  tab_UP_ImCctr[[i]] <- GenTable(tgd,
                                 Fisher.elim = resultTopGO.elim,
                                 orderBy = "Fisher.classic",
                                 topNodes = 200)
}
upGOres_ImCctr <- rbind.fill(tab_UP_ImCctr)


# DOWN --------------------------------------------------------------------
tab_DOWN_ImCctr <- as.list(onts)
names(tab_DOWN_ImCctr) <- onts
for (i in 1:length(onts)) {
  tgd <- new("topGOdata",
             ontology = onts[i],
             allGenes = allDOWN_ImCctr,
             nodeSize = 5,
             annot = annFUN.org,
             mapping = "org.Mm.eg.db",
             ID = "ensembl")
  
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher")
  
  tab_DOWN_ImCctr[[i]] <- GenTable(tgd,
                                   Fisher.elim = resultTopGO.elim,
                                   orderBy = "Fisher.classic",
                                   topNodes = 200)
}
downGOres_ImCctr <- rbind.fill(tab_DOWN_ImCctr)


# Save --------------------------------------------------------------------
write.csv(upGOres_ImCctr[order(upGOres_ImCctr$Fisher.elim), ],
          file = paste0("./output/", ver, "/GO/UP_ImCctr.csv"),
          row.names = FALSE)
write.csv(downGOres_ImCctr[order(downGOres_ImCctr$Fisher.elim), ],
          file = paste0("./output/", ver, "/GO/DOWN_ImCctr.csv"),
          row.names = FALSE)
