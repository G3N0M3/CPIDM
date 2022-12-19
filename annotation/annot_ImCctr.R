#######
# Annotating Gene list with gene symbol based on Ensembl ID (ImCctr)
#######

# Refine ------------------------------------------------------------------
# removing numbers behind periods
# ref: http://asia.ensembl.org/Help/Faq?id=488
if ("annot_ImCctr" %in% ls()) {rm(annot_ImCctr)}
annot_ImCctr <- res_ImCctr
row.names(annot_ImCctr) <- sapply(row.names(annot_ImCctr),
                                FUN = function(x) {strsplit(x, "\\.")[[1]][1]},
                                USE.NAMES = FALSE)

if (length(unique(res_ImCctr)) == length(unique(annot_ImCctr))) {
  print("No genes overlap after refining row names")
} else {
  print("WARNING: SOME GENES HAVE OVERLAPPING NAMES!!")
}


# Annotation --------------------------------------------------------------
annot_ImCctr$symbol <- mapIds(org.Mm.eg.db,  # AnnotationDb object
                              keys = row.names(annot_ImCctr),  # keys to select from db
                              column = "SYMBOL",  # the column to search on (match)
                              keytype = "ENSEMBL",  # annotation type of keys
                              multiVals = "first"  # action upon multiple mapping
                              )


# Save List af Ordering ---------------------------------------------------
ord_annot_ImCctr <- annot_ImCctr[order(annot_ImCctr$padj), ]
row.names(ord_annot_ImCctr)[!is.na(ord_annot_ImCctr$symbol)] <-
  ord_annot_ImCctr$symbol[!is.na(ord_annot_ImCctr$symbol)]
ord_annot_ImCctr <- ord_annot_ImCctr[order(ord_annot_ImCctr$padj, decreasing = FALSE), ]
# full
write.csv(as.data.frame(ord_annot_ImCctr[!is.na(annot_ImCctr$padj), ]),
          file = paste0("./output/", ver, "/annotation/annot_ImCctr.csv"))
# cut-off (padj < 0.05 and abs(FC) > 2)
write.csv(as.data.frame(ord_annot_ImPC[!is.na(ord_annot_ImCctr$padj) &
                                         as.numeric(ord_annot_ImCctr$padj) < 0.1 &
                                         abs(as.numeric(ord_annot_ImCctr$log2FoldChange)) > 1, ]),
          file = paste0("./output/", ver, "/annotation/annot_cut_ImCctr.csv"))
