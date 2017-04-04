#### plotting geneaology

########### Library phylotate
library(phylotate)
library(ape)
pardefault <- par(no.readonly = T)
par(mfrow=c(2,2))

tree <- read_annotated(filename = paste0(analysis.dir,"out.trees.txt"), format = "newick")

plot(tree, type = "phylogram", show.tip.label = FALSE)
plot(tree, show.tip.label = FALSE)
