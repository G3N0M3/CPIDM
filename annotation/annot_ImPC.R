#######
# Annotating Gene list with gene symbol based on Ensembl ID (ImPC)
#######

# Refine ------------------------------------------------------------------
# removing numbers behind periods
# ref: http://asia.ensembl.org/Help/Faq?id=488
if ("annot_ImPC" %in% ls()) {rm(annot_ImPC)}
annot_ImPC <- res_ImPC
row.names(annot_ImPC) <- sapply(row.names(annot_ImPC),
                                FUN = function(x) {strsplit(x, "\\.")[[1]][1]},
                                USE.NAMES = FALSE)

if (length(unique(res_ImPC)) == length(unique(annot_ImPC))) {
  print("No genes overlap after refining row names")
} else {
  print("WARNING: SOME GENES HAVE OVERLAPPING NAMES!!")
}


# Annotation --------------------------------------------------------------
annot_ImPC$symbol <- mapIds(org.Mm.eg.db,  # AnnotationDb object
                            keys = row.names(annot_ImPC),  # keys to select from db
                            column = "SYMBOL",  # the column to search on (match)
                            keytype = "ENSEMBL",  # annotation type of keys
                            multiVals = "first"  # action upon multiple mapping
                            )


# Save List af Ordering ---------------------------------------------------
ord_annot_ImPC <- annot_ImPC[order(annot_ImPC$padj), ]
row.names(ord_annot_ImPC)[!is.na(ord_annot_ImPC$symbol)] <-
  ord_annot_ImPC$symbol[!is.na(ord_annot_ImPC$symbol)]
ord_annot_ImPC <- ord_annot_ImPC[order(ord_annot_ImPC$padj, decreasing = FALSE), ]
# full
write.csv(as.data.frame(ord_annot_ImPC[!is.na(ord_annot_ImPC$padj), ]),
          file = paste0("./output/", ver, "/annotation/annot_ImPC.csv"))
# cut-off (padj < 0.05 and abs(FC) > 2)
write.csv(as.data.frame(ord_annot_ImPC[!is.na(ord_annot_ImPC$padj) &
                                         as.numeric(ord_annot_ImPC$padj) < 0.1 &
                                         abs(as.numeric(ord_annot_ImPC$log2FoldChange)) > 1, ]),
          file = paste0("./output/", ver, "/annotation/annot_cut_ImPC.csv"))
