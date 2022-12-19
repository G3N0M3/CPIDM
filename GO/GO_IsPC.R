#######
# Script for GO (IsPC)
# There are no significantly down-regulated genes in IsPC
#######


# GO ----------------------------------------------------------------------
bgGenes_IsPC <- row.names(annot_IsPC)
UP_IsPC <- row.names(annot_IsPC[annot_IsPC$padj < 0.1 &
                                  !is.na(annot_IsPC$padj) &
                                  annot_IsPC$log2FoldChange > 1, ]
)
allUP_IsPC <- factor(as.integer(bgGenes_IsPC %in% UP_IsPC))
names(allUP_IsPC) <- bgGenes_IsPC


onts <- c("MF", "BP", "CC")


# UP ----------------------------------------------------------------------
tab_UP_IsPC <- as.list(onts)
names(tab_UP_IsPC) <- onts
for (i in 1:length(onts)) {
  tgd <- new("topGOdata",
             ontology = onts[i],
             allGenes = allUP_IsPC,
             nodeSize = 5,
             annot = annFUN.org,
             mapping = "org.Mm.eg.db",
             ID = "ensembl")
  
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher")
  
  tab_UP_IsPC[[i]] <- GenTable(tgd,
                               Fisher.elim = resultTopGO.elim,
                               orderBy = "Fisher.classic",
                               topNodes = 200)
}
upGOres_IsPC <- rbind.fill(tab_UP_IsPC)


# DOWN --------------------------------------------------------------------
# There are no significantly down-regulated genes


# Save --------------------------------------------------------------------
write.csv(upGOres_IsPC[order(upGOres_IsPC$Fisher.elim), ],
          file = paste0("./output/", ver, "/GO/UP_IsPC.csv"),
          row.names = FALSE)
