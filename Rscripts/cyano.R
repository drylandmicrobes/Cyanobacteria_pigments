
#load the necessary packages 
library(ggtree)
library(phytools)
library(Biostrings)
library(readr)
library(seqinr)




#read crtI gene tree from Newick file  
crtI_gene_tree<-read.tree("crtI_duplication_neuro.genetree-labels.newick")


#build tree with MSA visualization 
p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
p
msaplot(p,"crtI_msa.fasta", offset=3, width=2)

#build tree with a small portion of MSA data 
p <- ggtree(crtI_gene_tree) + 
    geom_tiplab(offset=4, align=TRUE) + xlim(NA, 12)
p
msaplot(p,"crtI_msa.fasta", window=c(120, 200))  

#could be really good 
p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
p
msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(330,360))

#around 350 there is some good solid differences in the aa sequences but the rest
#still looks pretty conserved with lots of blank spaces 
p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
p
msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(350,400))

#also good, but too many blank spaces so it looks like the alignment may not be the best... 
#p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
#msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(250,300))

#really good just alot of gaps 
p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
p
msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(200,250))

#also really good, I think this will be the one to use! 
p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(210,240))


#some good conserved areas 
p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
p
msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(150,200))

#highly conserved, start of alignments (dont go lower)
#p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
#msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(100,150))

#one good row showing one aa defining the two diff clades but is it only 1 aa causing a major diff in the phylogenetics of the gene? 
p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
p
msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(400,450))

#not good for figure
#p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
#msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(450,500))

#not the best
#p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
#msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(550,605))

#no good 
#p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
#msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(550,600))

#no good 
#p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
#msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(600,650))

#no good 
#p <- ggtree(crtI_gene_tree) + geom_tiplab(size=3)
#msaplot(p,"crtI_msa.fasta", offset=3, width=2, window=c(650,700))

#set working directory to folder with newick tree and msa fasta files

#load the necessary packages 
library(ggtree)
library(phytools)
library(Biostrings)
library(readr)
library(seqinr)

#read crtW gene tree from Newick file  
crtW_gene_tree<-read.newick("crtW.genetree-labels.newick")


#build tree with MSA visualization 
q <- ggtree(crtW_gene_tree) + geom_tiplab(size=3)
q
msaplot(q,"crtW_msa.fasta", offset=3, width=2)


#Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

