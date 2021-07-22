---
title: "pHMMer_genetrees"
output: pdf_document
---

#Using this tutorial (https://rforbiochemists.blogspot.com/2017/01/playing-with-phylogenetic-tree.html)

library(ape)
library(ggtree)

#load tree
mytree_crtBa <- '/bigdata/stajichlab/shared/projects/Cyanobacteria/Cyanobacteria_pigments/genetree_pHMMer/crtBa/crtBa.newick'
tree_crtBa <- read.newick(mytree_crtBa)
tree_crtBa

#look at tree data
#plot(tree_crtBa, "u", 
# use.edge.length = FALSE,
#edge.color = "black")
#class(tree_crtBa)
#str(tree_crtBa)
#tree_crtBa$tip.label

#load in habitat data
habitat <- read.table('/bigdata/stajichlab/shared/projects/Cyanobacteria/Cyanobacteria_pigments/genetree_pHMMer/JasonTest/species_habitat.tsv'
, sep = '\t', col.names = c("Species","Habitat"))

#set up tip colors 
colors <- habitat$Habitat
colors[] <- c("cornflowerblue", "forest green", "darkorange1", "saddlebrown")
#Warning message:
#number of items to replace is not a multiple of replacement length

colors[habitat$Habitat == c("aqutic")] <- "cornflowerblue"
colors[habitat$Habitat == c("hydroterrestrial")] <- "forest green"
colors[habitat$Habitat == c("soil")] <- "darkorange1"
colors[habitat$Habitat == c("subaerial")] <- "saddlebrown"
colors[habitat$Habitat == c("unknown")] <- "grey"
tipcol <- rep('black', length(tree_crtBa$tip.label))

#loop to replace tip colors of tip labels
for(i in 1:length(habitat$Species)){
  tipcol[grep(habitat$Species[i],tree_crtBa$tip.label)] <- colors[i]
}

#gene tree with colored taxa
plot(tree_crtBa,
     use.edge.length = FALSE,
     tip.color = tipcol,
     cex = 0.25)