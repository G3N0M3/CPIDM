#######
# script that creates 3D (or 2D) PCA plots for Immune set
#######


# Input data --------------------------------------------------------------
choice_ntop = 0  # 1: use ntop, 0: use DEG

if (choice_ntop) {
  PC_data_Im <- assay(rld_Im)[order(rowVars(assay(rld_Im)), decreasing = TRUE), ][1:ntop_Im, ]
} else {
  # rld_annot_Im is used since DEG_ImPC is annotated
  select <- row.names(rld_annot_Im) %in% DEG_ImPC
  PC_data_Im <- assay(rld_annot_Im)[select, ]
}


# Calculate PCs -----------------------------------------------------------
prin_comp <- prcomp(t(PC_data_Im), rank. = 3)

components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components <- cbind(components, rld_Im$samples)
rownames(components) <- rld_Im$samples
head(components)

summary(prin_comp)
tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]["Cumulative Proportion", 3]

tit <- paste0("Total Explained Variance = ",
             round(tot_explained_variance_ratio * 100,
                   digits = 2), "%")
tit


# Plot --------------------------------------------------------------------
fig <- plot_ly(components,
               x = ~PC1,
               y = ~PC2,
               z = ~PC3,
               color = ~rld_Im$condition,
               colors = c('#00BC48','#fa7970','#6A96FA')) %>%
  add_markers(size = 12)
fig <- fig %>%
  layout(
    title = tit,
    scene = list(bgcolor = "#ffffff")
  )
fig
