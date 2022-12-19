#######
# Creates GO barplot for ImPC
#######


# Input data --------------------------------------------------------------
go_up_ImPC <- read.csv(file = paste0("./output/", ver, "/GO/UP_ImPC.csv"))
go_up_ImPC$pval <- as.numeric(-log10(go_up_ImPC$Fisher.elim))
go_up_ImPC <- go_up_ImPC[order(go_up_ImPC$pval,
                               decreasing = TRUE), ]
go_up_ImPC <- go_up_ImPC[1:10, ]


# Plot --------------------------------------------------------------------
goBar_ImPC <- ggplot(go_up_ImPC) +
  geom_bar(aes(x = pval, y = GO.ID),
           stat = 'identity') +
  scale_y_discrete(limits = rev(go_up_ImPC$GO.ID)) +
  theme_light() +
  coord_cartesian(xlim = c(0, 11)) +
  scale_x_continuous(expand = c(0, 0)) +
  geom_text(aes(x = .15, y = GO.ID, label = Term),
            color = "white", hjust = 0, fontface = "bold") +
  labs(x = "-Log10(p value)", y = "") +
  theme(axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12))
goBar_ImPC
ggsave("goBarplot_ImPC.png",
       plot = goBar_ImPC,
       device = "png",
       path = paste0("./output/", ver, "/GO_barplot"),
       width = 5000, height = 3000, units = "px", dpi = 600
)
