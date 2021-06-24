#library(ggplot2)
library(ggtree)
#library(tidyverse)

treeext = "_long.tre"
infolder = "genetree_BLASTP"
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
    print(outfile)
    tree <- read.tree(treefile)
}


treefiles <- list.files(infolder,pattern="_FT_LG_GAMMA_long.tre$",full.names=TRUE)

printTree(treefiles[1],outfolder)
