#library(ggplot2)
library(ggtree)
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
    pt <- ggtree(tree) + geom_tiplab(font.size=2)
    ggsave(outfile,pt,width=15,height=15) 
}


treefiles <- list.files(infolder,pattern=treeext,full.names=TRUE)

lapply(treefiles,printTree,outfolder)
