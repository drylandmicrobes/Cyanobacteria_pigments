#library(ggplot2)
library(ggtree)
library(phytools)
#library(tidyverse)

treeext = "_FT_LG_GAMMA_long.tre"
infolder = "genetree_pHMMer"
outfolder = ""
# args would be
# Rscripts generate_gene_trees.R INFOLDER OUTFOLDER

args <- commandArgs(trailingOnly = TRUE)
if( length(args) > 0 ) {
    infolder = args[1]
    if( length(args) > 1 ) {
        outfolder = args[2]
    }
}

if(outfolder == "" ) {
    outfolder = infolder
}


printTree <- function(treefile,odir) {
    gene_name = sub(pattern = "(.*)_long.tre$", replacement = "\\1",
                    basename(treefile))
    outfile <- file.path(odir,sprintf("%s.pdf",gene_name))
    tree <- read.tree(treefile)
#    rootedtree <- root(tree, outgroup)
    rootedtree <- midpoint.root(tree)
#    rootedtree <- tree
    pt <- ggtree(rootedtree,layout="circular") + 
        geom_tiplab(font.size=1)
    ggsave(outfile,pt,width=30,height=30) 
}

treefiles <- list.files(infolder,pattern=treeext,full.names=TRUE)

lapply(treefiles,printTree,outfolder)
