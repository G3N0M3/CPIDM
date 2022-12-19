#######
# Annotating Gene list with gene symbol based on Ensembl ID (ImPctr)
#######

# Refine ------------------------------------------------------------------
# removing numbers behind periods
# ref: http://asia.ensembl.org/Help/Faq?id=488
if ("annot_ImPctr" %in% ls()) {rm(annot_ImPctr)}
annot_ImPctr <- res_ImPctr
row.names(annot_ImPctr) <- sapply(row.names(annot_ImPctr),
                                FUN = function(x) {strsplit(x, "\\.")[[1]][1]},
                                USE.NAMES = FALSE)

if (length(unique(res_ImPctr)) == length(unique(annot_ImPctr))) {
  print("No genes overlap after refining row names")
} else {
  print("WARNING: SOME GENES HAVE OVERLAPPING NAMES!!")
}


# Annotation --------------------------------------------------------------
annot_ImPctr$symbol <- mapIds(org.Mm.eg.db,  # AnnotationDb object
                              keys = row.names(annot_ImPctr),  # keys to select from db
                              column = "SYMBOL",  # the column to search on (match)
                              keytype = "ENSEMBL",  # annotation type of keys
                              multiVals = "first"  # action upon multiple mapping
                              )


# Save List af Ordering ---------------------------------------------------
ord_annot_ImPctr <- annot_ImPctr[order(annot_ImPctr$padj), ]
row.names(ord_annot_ImPctr)[!is.na(ord_annot_ImPctr$symbol)] <-
  ord_annot_ImPctr$symbol[!is.na(ord_annot_ImPctr$symbol)]
ord_annot_ImPctr <- ord_annot_ImPctr[order(ord_annot_ImPctr$padj, decreasing = FALSE), ]
# full
write.csv(as.data.frame(ord_annot_ImPctr[!is.na(annot_ImPctr$padj), ]),
          file = paste0("./output/", ver, "/annotation/annot_ImPctr.csv"))
# cut-off (padj < 0.05 and abs(FC) > 2)
write.csv(as.data.frame(ord_annot_ImPctr[!is.na(ord_annot_ImPctr$padj) &
                                           as.numeric(ord_annot_ImPctr$padj) < 0.1 &
                                           abs(as.numeric(ord_annot_ImPctr$log2FoldChange)) > 1, ]),
          file = paste0("./output/", ver, "/annotation/annot_cut_ImPctr.csv"))
