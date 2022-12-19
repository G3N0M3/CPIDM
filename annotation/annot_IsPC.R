#######
# Annotating Gene list with gene symbol based on Ensembl ID (IsPC)
#######

# Refine ------------------------------------------------------------------
# removing numbers behind periods
# ref: http://asia.ensembl.org/Help/Faq?id=488
if ("annot_IsPC" %in% ls()) {rm(annot_IsPC)}
annot_IsPC <- res_IsPC
row.names(annot_IsPC) <- sapply(row.names(annot_IsPC),
                                FUN = function(x) {strsplit(x, "\\.")[[1]][1]},
                                USE.NAMES = FALSE)

if (length(unique(res_IsPC)) == length(unique(annot_IsPC))) {
  print("No genes overlap after refining row names")
} else {
  print("WARNING: SOME GENES HAVE OVERLAPPING NAMES!!")
}


# Annotation --------------------------------------------------------------
annot_IsPC$symbol <- mapIds(org.Mm.eg.db,  # AnnotationDb object
                            keys = row.names(annot_IsPC),  # keys to select from db
                            column = "SYMBOL",  # the column to search on (match)
                            keytype = "ENSEMBL",  # annotation type of keys
                            multiVals = "first"  # action upon multiple mapping
                            )


# Save List af Ordering ---------------------------------------------------
ord_annot_IsPC <- annot_IsPC[order(annot_IsPC$padj), ]
row.names(ord_annot_IsPC)[!is.na(ord_annot_IsPC$symbol)] <-
  ord_annot_IsPC$symbol[!is.na(ord_annot_IsPC$symbol)]
ord_annot_IsPC <- ord_annot_IsPC[order(ord_annot_IsPC$padj, decreasing = FALSE), ]
# full
write.csv(as.data.frame(ord_annot_IsPC[!is.na(annot_IsPC$padj), ]),
          file = paste0("./output/", ver, "/annotation/annot_IsPC.csv"))
# cut-off (padj < 0.05 and abs(FC) > 2)
write.csv(as.data.frame(ord_annot_IsPC[!is.na(ord_annot_IsPC$padj) &
                                         as.numeric(ord_annot_IsPC$padj) < 0.1 &
                                         abs(as.numeric(ord_annot_IsPC$log2FoldChange)) > 1, ]),
          file = paste0("./output/", ver, "/annotation/annot_cut_IsPC.csv"))
