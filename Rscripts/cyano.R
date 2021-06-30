
#load the necessary packages 
library(ggtree)
library(phytools)
library(Biostrings)
library(readr)
library(seqinr)




#read crtI gene tree from Newick file  
crtI_gene_tree<-read.tree("crtI_duplication_neuro.genetree-labels.newick")


#build tree with MSA visualization 
p <- ggtree(crtI_gene_tree) + geom_tiplab(size=2)
msaplot(p,"crtI_msa.fasta", offset=3, width=2)

#build tree with a small portion of MSA data 
p <- ggtree(crtI_gene_tree) + 
    geom_tiplab(offset=4, align=TRUE) + xlim(NA, 12)
msaplot(p,"crtI_msa.fasta", window=c(120, 200))  

#could be really good 
p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(330,360))

#around 350 there is some good solid differences in the aa sequences but the rest
#still looks pretty conserved with lots of blank spaces 
p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(350,400))

#also good, but too many blank spaces so it looks like the alignment may not be the best... 
#p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
#msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(250,300))

#really good just alot of gaps 
p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
p
msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(200,250))

#also really good, I think this will be the one to use! 
p <- ggtree(crtI_gene_tree) + geom_tiplab(size=2.75)
msaplot(p,"crtI_msa.fasta", offset=1, width=2, window=c(210,240))


#some good conserved areas 
p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(150,200))

#highly conserved, start of alignments (dont go lower)
#p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
#msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(100,150))

#one good row showing one aa defining the two diff clades but is it only 1 aa causing a major diff in the phylogenetics of the gene? 
p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(400,450))


#set working directory to folder with newick tree and msa fasta files

#load the necessary packages 
library(ggtree)
library(phytools)
library(Biostrings)
library(readr)
library(seqinr)

#read crtW gene tree from Newick file  
crtW_gene_tree<-read.newick("crtW.genetree-labels.newick")
habitat_info <- read_csv("Habitat_info - Sheet1.csv")

v <- ggtree(crtW_gene_tree) %<+% habitat_info + geom_tippoint((aes(color="Habitat")))

#build tree with MSA visualization 
v <- ggtree(crtW_gene_tree) %<+% habitat_info + geom_tippoint((aes(color=c("forest green","orange","orange","light blue","dark blue","orange","dark blue","dark blue", "orange", "orange","dark blue","light blue","dark blue","dark blue","orange","orange","orange")))) + geom_tiplab(size=3)
msaplot(v,"crtW_msa.fasta", offset=2, width=2)

q <- ggtree(crtW_gene_tree) + geom_tiplab(size=3)
msaplot(q,"crtW_msa.fasta", offset=1, width=2, window=c(100,200))

q <- ggtree(crtW_gene_tree) + geom_tiplab(size=3)
msaplot(q,"crtW_msa.fasta", offset=3, width=2, window=c(100,150))


#Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

crtW_paper_tree<-read.newick("crtW_paper.genetree-labels.newick")

w <- ggtree(crtW_paper_tree) + geom_tiplab(size=2)
msaplot(w,"crtw_paper_msa.fasta", offset=1, width=2)
