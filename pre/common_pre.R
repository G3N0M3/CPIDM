# Sample Annotation -------------------------------------------------------
# Read file
colData <- read.csv("./data/annotation/annotation.csv")
head(colData)


# Create sampleTable ------------------------------------------------------
sampleFiles <- grep("htseqcount", list.files("./data/count/"), value = TRUE)
head(sampleFiles)
sampleNames <- sapply(sampleFiles,
                      FUN = function(x){strsplit(x, "\\.")[[1]][1]},
                      USE.NAMES = FALSE)
sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          age = colData$Age,
                          cellType = ifelse(colData$CD45 == "+", "Immune", "Islet"),
                          treatment = colData$Ab_mod,
                          condition = colData$Annot)