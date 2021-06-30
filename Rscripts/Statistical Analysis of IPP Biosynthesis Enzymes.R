---
title: "Statistical Analysis of IPP Biosynthesis Enzymes"
author: "Dionne Martin"
date: "6/28/2021"
output: pdf_document
---
  
#################### PREFACE ##########################
#Code for PCA, NMDS, and PERMANOVA/ simper test 
#of raw BLASTp hits of IPP biosynthesis enzymes: 
#dxs, hmgs, ispC, ispD, ispE, ispG, phaA. 
#These enzymes are from two distinct pathways (mevolonate
# and non-mevoloante pathway) to create the precursor 
#for all carotenoids, IPP. 
#######################################################

```{pca}
#load packages and blastp raw hits data 
library(devtools)
library(ggbiplot)
library(readr)
IPP <- read_csv("/bigdata/stajichlab/shared/projects/Cyanobacteria/Terresterial_Cyano/Cyano Carotenoids BLASTp Results - IPP dataset.csv")
#organize data to just raw hits and strain names 
row.names(IPP) = IPP$Strain
IPP_matrix <- data.matrix(IPP[c(3:9)])
range(IPP_matrix) #0 6 

#analyze pca statistics to find proportion of variatiance explained by each gene 
pca_IPP= prcomp(IPP_matrix, center= TRUE, scale.= TRUE)
summary(pca_IPP)
str(pca_IPP)
#dxs-0.32 hmgs-0.19 ispC-0.14 ispD-0.12 
#ispE-0.09 ispG-0.08 phaA-0.02 
 

#make groups based on habitat
habitat= IPP$Habitat

#try to make group colours, this should be added to plot code 
scale_colour_manual(name= "Habitat", values= c("brown", "dark orange", "blue", "pink", "red", "dark blue", "forest green", "dark yellow", "red", "black"))

#make the PCA plot with strains as points 
#png("PCA Plot of Carotenoid BLAST Hit .png", width=1000, height=1000)
ggbiplot(pca_IPP,ellipse=FALSE,obs.scale = 3,var.scale = 1, groups=hits.habitat) +
  scale_colour_manual(name="Habitat", values= c("saddlebrown", "darkorange1", "cornflowerblue", "deeppink", "red", "dark blue", "forest green", "gold", "darkmagenta", "black"))+
  ggtitle("PCA of Raw BLASTp Hits of IPP Enzymes")+
  theme_minimal()+
  theme(legend.position = "bottom")

#make the PCA plot with strain names (looks too messy with all names...)
#png("PCA Plot of Carotenoid BLAST Hit .png", width=1000, height=1000)
ggbiplot(pca_IPP,labels= rownames(IPP),ellipse=FALSE,obs.scale = 3,var.scale = 1, groups=hits.habitat) +
  scale_colour_manual(name="Habitat", values= c("saddlebrown", "darkorange1", "cornflowerblue", "deeppink", "red", "dark blue", "forest green", "gold", "darkmagenta", "black"))+
  ggtitle("PCA of Raw BLASTp Hits of IPP Enzymes")+
  theme_minimal()+
  theme(legend.position = "bottom")

```
```{nmds}
#load packages 
library(vegan)
library(ggplot2)
library(readr)

#prep data
IPP <- read_csv("/bigdata/stajichlab/shared/projects/Cyanobacteria/Terresterial_Cyano/Cyano Carotenoids BLASTp Results - IPP dataset.csv")
IPP_matrix <- data.matrix(IPP[c(3:9)])
nmds_IPP_habitat= as.data.frame(IPP[,1:2])
#Do you need to transform your data? range should be 0-10 if not, use root function
range(IPP_matrix) # 0-6 

#make a distance matrix
IPP_matrix_ds=vegdist(IPP_matrix, method = 'euclidian')

#Run metaMDS on distance matrix
IPP_nms=metaMDS(IPP_matrix_ds, distance="euclidian",k=2,trymax=500)
IPP_nms
#stress=0.038

#assign colors to habitat
cols <- nmds_IPP_habitat$Strain
cols[] <- c("saddlebrown", "darkorange1", "cornflowerblue", "deeppink", "red", "dark blue", "forest green", "gold", "darkmagenta", "black")
#Warning message:
#number of items to replace is not a multiple of replacement length

cols[nmds_IPP_habitat$Habitat == c("freshwater")] <- "cornflowerblue"
cols[nmds_IPP_habitat$Habitat == c("saltwater")] <- "dark blue"
cols[nmds_IPP_habitat$Habitat == c("hot spring")] <- "deeppink"
cols[nmds_IPP_habitat$Habitat == c("soil")] <- "forest green"
cols[nmds_IPP_habitat$Habitat == c("cave")] <- "saddlebrown"
cols[nmds_IPP_habitat$Habitat == c("desert")] <- "darkorange1"
cols[nmds_IPP_habitat$Habitat == c("tropical")] <- "black"
cols[nmds_IPP_habitat$Habitat == c("temperate")] <- "darkmagenta"
cols[nmds_IPP_habitat$Habitat == c("symbiont")] <- "gold"
cols[nmds_IPP_habitat$Habitat == c("rock")] <- "red"
pairs(character(nrow(nmds_IPP_habitat$Strain, col=cols)
                
#habitatvec = factor(small$Habitat, levels = c("freshwater", "saltwater", "hot spring", "soil", "cave", "desert", "tropical", "temperate", "symbiont", "rock"))
#colors <- c("cornflowerblue", "dark blue", "deeppink", "forest green", "saddlebrown", "darkorange1", "black", "darkmagenta", "gold", "red")
habitatinfo= c(rep(IPP$Habitat))
                
#NMDS Plot With Species Names 
ordiplot(IPP_nms, type="n")
orditorp(IPP_nms,display="sites",labels=IPP$Strain, col=cols ,air=0.00025)
par(xpd=TRUE)
#legend(0.6,0.25,legend =unique(nmds_small_habitat$Strain), fill=cols, border= "black")
#change coordinates of the legend for each specific graph
#legend(-0.4,-0.2, "stress = 0.13") 
#change coordinates and stress value for each specific graph
title(main="NMDS of Raw BLASTp")
scores <-
  scores(IPP_nms, display = "sites", "species")
cent <-
  aggregate(scores ~ habitatinfo, data = nmds_IPP_habitat, FUN = "mean")
names(cent) [-1] <- colnames(scores)
points(cent [,-1],
      pch = c( 8 , 8 , 8, 8),
      col = unique(colors),
      bg = c("black"),
      lwd = 3.0,
      cex = 2.0
)
                
#NMDS Plot with Points
plot(IPP_nms)
points(IPP_nms,display= "sites",pch=21, bg=cols) 
par(xpd=TRUE)
#legend(0.6,0.25,legend =habitatvec, fill=colvec, border= "black")
#legend(-0.4,-0.2, "stress = 0.13") 
title(main="NMDS of Raw BLASTp")
scores <-
  scores(IPP_nms, display = "sites", "species")
cent <-
  aggregate(scores ~ habitatinfo, data = nmds_IPP_habitat, FUN = "mean")
names(cent) [-1] <- colnames(scores)
points(cent [,-1],
      pch = c( 8 , 8 , 8, 8),
      col = unique(colors),
      bg = c("black"),
      lwd = 3.0,
      cex = 2.0
)
```
                
```{permanova}
                
#data prepartion 
rownames(IPP_matrix) = IPP$Strain
str(IPP_matrix)
nmds_IPP_habitat= as.data.frame(IPP[,1:2])
str(nmds_IPP_habitat)

#make a distance matrix
IPP_matrix_ds=vegdist(IPP_matrix, method = 'euclidian')
                
#Run metaMDS on distance matrix
IPP_nms=metaMDS(IPP_matrix_ds, distance="euclidian",k=2,trymax=500)
IPP_nms
#stress=0.038
                
#Run PERMANOVA
pmv_IPP = adonis(IPP_matrix ~ habitat, data = nmds_IPP_habitat, permutations = 999, method = 'euclidean')
pmv_IPP
#PERMANOVA Results: Based off of the carotenoid genes being compared
#cyanobacterial strains differ statistically significant 
#between their habitats (p-value:0.027).
#Habitats explian 13.66% of the variance in the copy numbers of critical 
#carotenoid biosynthesis enzymes between cyanobacterial species. 
                
#Visualise the permanova
densityplot(permustats(pmv_IPP))
                
#Distance Based Dispersion test 
bd_IPP = betadisper(IPP_matrix_ds, IPP$Habitat)
boxplot(bd_IPP)
anova(bd_IPP) 
#F-Test: pvalue 0.4245
permutest(bd_IPP)
#permutation test: pvalue 0.392
#the significance of a permutation test is represented by its P-value. 
#The P-value is the probability of obtaining a result at least 
#as extreme as the test statistic given that the null hypothesis is true. 
                
#both the F-Test and permutation test show that 
#the avergae distance to the median for the 
#habitat data is not statistically significant. 
                
#simper test: determine which genes are responsible for the difference between 2 compared habiatat groups 
simp_IPP= simper(IPP_matrix, group= IPP$Habitat)
summary(simp_IPP)

#Results of simper test
Contrast: saltwater_symbiont 

average      sd  ratio     ava    avb cumsum
dxs  0.076152 0.06269 1.2148 2.50000 3.6667 0.4780
phaA 0.043808 0.03119 1.4045 0.08333 0.6667 0.7529
hmgs 0.026323 0.03246 0.8108 0.20833 0.3333 0.9181
ispE 0.010192 0.02300 0.4431 1.16667 1.0000 0.9821
ispD 0.002852 0.01380 0.2067 0.95833 1.0000 1.0000
ispC 0.000000 0.00000    NaN 1.00000 1.0000 1.0000
ispG 0.000000 0.00000    NaN 1.00000 1.0000 1.0000

Contrast: saltwater_freshwater 

average      sd  ratio     ava     avb cumsum
dxs  0.059540 0.04972 1.1974 2.50000 3.03571 0.4774
phaA 0.031087 0.03905 0.7960 0.08333 0.42857 0.7267
hmgs 0.016171 0.02814 0.5746 0.20833 0.07143 0.8564
ispE 0.012572 0.02604 0.4828 1.16667 1.03571 0.9572
ispD 0.005335 0.01891 0.2822 0.95833 1.03571 1.0000
ispC 0.000000 0.00000    NaN 1.00000 1.00000 1.0000
ispG 0.000000 0.00000    NaN 1.00000 1.00000 1.0000

Contrast: saltwater_cave 

average      sd  ratio     ava avb cumsum
dxs  0.078977 0.05560 1.4205 2.50000 3.3 0.5155
ispE 0.018173 0.03283 0.5536 1.16667 0.9 0.6342
hmgs 0.017741 0.02964 0.5985 0.20833 0.1 0.7500
phaA 0.016780 0.03230 0.5196 0.08333 0.2 0.8595
ispD 0.009669 0.02479 0.3900 0.95833 0.9 0.9226
ispC 0.005926 0.01784 0.3322 1.00000 1.1 0.9613
ispG 0.005926 0.01784 0.3322 1.00000 1.1 1.0000

Contrast: saltwater_desert 

average      sd  ratio     ava     avb cumsum
dxs  0.050717 0.04334 1.1701 2.50000 3.03704 0.3454
phaA 0.046487 0.05304 0.8764 0.08333 0.66667 0.6621
hmgs 0.016648 0.02886 0.5768 0.20833 0.07407 0.7755
ispE 0.015109 0.02891 0.5227 1.16667 1.00000 0.8784
ispC 0.007623 0.02184 0.3490 1.00000 0.96296 0.9303
ispD 0.007563 0.02170 0.3485 0.95833 1.00000 0.9818
ispG 0.002671 0.01366 0.1956 1.00000 0.96296 1.0000

Contrast: saltwater_tropical 

average      sd  ratio     ava    avb cumsum
dxs  0.057726 0.04841 1.1925 2.50000 3.0000 0.5760
phaA 0.014084 0.03191 0.4413 0.08333 0.1429 0.7165
hmgs 0.013997 0.02744 0.5102 0.20833 0.0000 0.8561
ispE 0.011237 0.02526 0.4448 1.16667 1.0000 0.9682
ispD 0.003183 0.01535 0.2074 0.95833 1.0000 1.0000
ispC 0.000000 0.00000    NaN 1.00000 1.0000 1.0000
ispG 0.000000 0.00000    NaN 1.00000 1.0000 1.0000

Contrast: saltwater_soil 

average      sd  ratio     ava avb cumsum
dxs  0.040574 0.03684 1.1014 2.50000 2.8 0.4065
phaA 0.031001 0.03732 0.8307 0.08333 0.4 0.7171
hmgs 0.013909 0.02725 0.5104 0.20833 0.0 0.8565
ispE 0.011166 0.02510 0.4449 1.16667 1.0 0.9683
ispD 0.003159 0.01522 0.2076 0.95833 1.0 1.0000
ispC 0.000000 0.00000    NaN 1.00000 1.0 1.0000
ispG 0.000000 0.00000    NaN 1.00000 1.0 1.0000

Contrast: saltwater_rock 

average      sd  ratio     ava avb cumsum
dxs  0.042705 0.03704 1.1530 2.50000   3 0.5581
hmgs 0.014087 0.02807 0.5019 0.20833   0 0.7422
ispE 0.011310 0.02585 0.4376 1.16667   1 0.8900
phaA 0.005208 0.02552 0.2041 0.08333   0 0.9581
ispD 0.003205 0.01570 0.2041 0.95833   1 1.0000
ispC 0.000000 0.00000    NaN 1.00000   1 1.0000
ispG 0.000000 0.00000    NaN 1.00000   1 1.0000

Contrast: saltwater_hot spring 

average      sd  ratio     ava avb cumsum
dxs  0.100941 0.10191 0.9905 2.50000 2.6 0.4909
phaA 0.028856 0.05408 0.5336 0.08333 0.4 0.6313
ispC 0.022633 0.04570 0.4952 1.00000 0.8 0.7414
ispG 0.022633 0.04570 0.4952 1.00000 0.8 0.8514
hmgs 0.015015 0.03025 0.4964 0.20833 0.0 0.9245
ispE 0.012059 0.02784 0.4332 1.16667 1.0 0.9831
ispD 0.003475 0.01729 0.2010 0.95833 1.0 1.0000

Contrast: saltwater_temperate 

average      sd  ratio     ava avb cumsum
dxs  0.123968 0.09280 1.3359 2.50000 4.5 0.6700
phaA 0.035822 0.03525 1.0163 0.08333 0.5 0.8636
hmgs 0.012452 0.02460 0.5063 0.20833 0.0 0.9309
ispE 0.009993 0.02264 0.4414 1.16667 1.0 0.9849
ispD 0.002790 0.01355 0.2059 0.95833 1.0 1.0000
ispC 0.000000 0.00000    NaN 1.00000 1.0 1.0000
ispG 0.000000 0.00000    NaN 1.00000 1.0 1.0000

Contrast: symbiont_freshwater 

average      sd  ratio    ava     avb cumsum
dxs  0.067010 0.05403 1.2402 3.6667 3.03571 0.5190
phaA 0.035249 0.03215 1.0963 0.6667 0.42857 0.7921
hmgs 0.022671 0.03080 0.7361 0.3333 0.07143 0.9676
ispE 0.002149 0.01125 0.1910 1.0000 1.03571 0.9843
ispD 0.002027 0.01061 0.1910 1.0000 1.03571 1.0000
ispC 0.000000 0.00000    NaN 1.0000 1.00000 1.0000
ispG 0.000000 0.00000    NaN 1.0000 1.00000 1.0000

Contrast: symbiont_cave 

average      sd  ratio    ava avb cumsum
dxs  0.063906 0.06438 0.9926 3.6667 3.3 0.4320
phaA 0.037090 0.03118 1.1897 0.6667 0.2 0.6827
hmgs 0.023397 0.03151 0.7426 0.3333 0.1 0.8408
ispD 0.006405 0.01958 0.3271 1.0000 0.9 0.8841
ispE 0.006405 0.01958 0.3271 1.0000 0.9 0.9274
ispC 0.005370 0.01641 0.3273 1.0000 1.1 0.9637
ispG 0.005370 0.01641 0.3273 1.0000 1.1 1.0000

Contrast: symbiont_desert 

average      sd  ratio    ava     avb cumsum
dxs  0.051741 0.05432 0.9526 3.6667 3.03704 0.3919
phaA 0.039642 0.04007 0.9893 0.6667 0.66667 0.6922
hmgs 0.022730 0.03074 0.7395 0.3333 0.07407 0.8644
ispC 0.006795 0.01953 0.3478 1.0000 0.96296 0.9159
ispE 0.004475 0.01598 0.2801 1.0000 1.00000 0.9498
ispD 0.004260 0.01528 0.2787 1.0000 1.00000 0.9820
ispG 0.002372 0.01219 0.1946 1.0000 0.96296 1.0000

Contrast: symbiont_tropical 

average      sd  ratio    ava    avb cumsum
dxs  0.06515 0.05630 1.1572 3.6667 3.0000 0.5164
phaA 0.03893 0.03154 1.2341 0.6667 0.1429 0.8250
hmgs 0.02208 0.03207 0.6885 0.3333 0.0000 1.0000
ispC 0.00000 0.00000    NaN 1.0000 1.0000 1.0000
ispD 0.00000 0.00000    NaN 1.0000 1.0000 1.0000
ispE 0.00000 0.00000    NaN 1.0000 1.0000 1.0000
ispG 0.00000 0.00000    NaN 1.0000 1.0000 1.0000

Contrast: symbiont_soil 

average      sd  ratio    ava avb cumsum
dxs  0.05159 0.06209 0.8309 3.6667 2.8 0.4811
phaA 0.03371 0.03275 1.0291 0.6667 0.4 0.7954
hmgs 0.02194 0.03214 0.6828 0.3333 0.0 1.0000
ispC 0.00000 0.00000    NaN 1.0000 1.0 1.0000
ispD 0.00000 0.00000    NaN 1.0000 1.0 1.0000
ispE 0.00000 0.00000    NaN 1.0000 1.0 1.0000
ispG 0.00000 0.00000    NaN 1.0000 1.0 1.0000

Contrast: symbiont_rock 

average      sd  ratio    ava avb cumsum
phaA 0.04183 0.03644 1.1480 0.6667   0 0.4051
dxs  0.03922 0.06792 0.5774 3.6667   3 0.7848
hmgs 0.02222 0.03849 0.5774 0.3333   0 1.0000
ispC 0.00000 0.00000    NaN 1.0000   1 1.0000
ispD 0.00000 0.00000    NaN 1.0000   1 1.0000
ispE 0.00000 0.00000    NaN 1.0000   1 1.0000
ispG 0.00000 0.00000    NaN 1.0000   1 1.0000

Contrast: symbiont_hot spring 

average      sd  ratio    ava avb cumsum
dxs  0.10252 0.13246 0.7740 3.6667 2.6 0.4745
phaA 0.05210 0.03694 1.4103 0.6667 0.4 0.7157
hmgs 0.02364 0.03573 0.6618 0.3333 0.0 0.8251
ispC 0.01889 0.03927 0.4810 1.0000 0.8 0.9126
ispG 0.01889 0.03927 0.4810 1.0000 0.8 1.0000
ispD 0.00000 0.00000    NaN 1.0000 1.0 1.0000
ispE 0.00000 0.00000    NaN 1.0000 1.0 1.0000

Contrast: symbiont_temperate 

average      sd  ratio    ava avb cumsum
dxs  0.08241 0.07701 1.0700 3.6667 4.5 0.6335
phaA 0.02801 0.03094 0.9054 0.6667 0.5 0.8488
hmgs 0.01968 0.03056 0.6438 0.3333 0.0 1.0000
ispC 0.00000 0.00000    NaN 1.0000 1.0 1.0000
ispD 0.00000 0.00000    NaN 1.0000 1.0 1.0000
ispE 0.00000 0.00000    NaN 1.0000 1.0 1.0000
ispG 0.00000 0.00000    NaN 1.0000 1.0 1.0000

Contrast: freshwater_cave 

average      sd  ratio     ava avb cumsum
dxs  0.070732 0.05854 1.2082 3.03571 3.3 0.5046
phaA 0.030036 0.03516 0.8544 0.42857 0.2 0.7189
hmgs 0.009735 0.02271 0.4287 0.07143 0.1 0.7883
ispE 0.009199 0.02403 0.3828 1.03571 0.9 0.8540
ispD 0.009058 0.02362 0.3834 1.03571 0.9 0.9186
ispC 0.005705 0.01719 0.3319 1.00000 1.1 0.9593
ispG 0.005705 0.01719 0.3319 1.00000 1.1 1.0000

Contrast: freshwater_desert 

average      sd  ratio     ava     avb cumsum
dxs  0.052166 0.04348 1.1996 3.03571 3.03704 0.4074
phaA 0.044163 0.04705 0.9386 0.42857 0.66667 0.7523
hmgs 0.008382 0.02139 0.3919 0.07143 0.07407 0.8178
ispC 0.007293 0.02091 0.3488 1.00000 0.96296 0.8747
ispE 0.006922 0.02038 0.3397 1.03571 1.00000 0.9288
ispD 0.006565 0.01939 0.3386 1.03571 1.00000 0.9801
ispG 0.002552 0.01306 0.1954 1.00000 0.96296 1.0000

Contrast: freshwater_tropical 

average      sd  ratio     ava    avb cumsum
dxs  0.059840 0.04981 1.2012 3.03571 3.0000 0.6070
phaA 0.029983 0.03601 0.8326 0.42857 0.1429 0.9111
hmgs 0.004177 0.01512 0.2763 0.07143 0.0000 0.9535
ispE 0.002366 0.01234 0.1916 1.03571 1.0000 0.9775
ispD 0.002218 0.01157 0.1917 1.03571 1.0000 1.0000
ispC 0.000000 0.00000    NaN 1.00000 1.0000 1.0000
ispG 0.000000 0.00000    NaN 1.00000 1.0000 1.0000

Contrast: freshwater_soil 

average      sd  ratio     ava avb cumsum
dxs  0.049810 0.04086 1.2189 3.03571 2.8 0.5371
phaA 0.034217 0.03577 0.9566 0.42857 0.4 0.9061
hmgs 0.004155 0.01504 0.2763 0.07143 0.0 0.9509
ispE 0.002351 0.01227 0.1917 1.03571 1.0 0.9762
ispD 0.002206 0.01151 0.1917 1.03571 1.0 1.0000
ispC 0.000000 0.00000    NaN 1.00000 1.0 1.0000
ispG 0.000000 0.00000    NaN 1.00000 1.0 1.0000

Contrast: freshwater_rock 

average      sd  ratio     ava avb cumsum
dxs  0.046118 0.03646 1.2650 3.03571   3 0.5571
phaA 0.027855 0.03665 0.7600 0.42857   0 0.8935
hmgs 0.004202 0.01543 0.2724 0.07143   0 0.9443
ispE 0.002381 0.01260 0.1890 1.03571   1 0.9730
ispD 0.002232 0.01181 0.1890 1.03571   1 1.0000
ispC 0.000000 0.00000    NaN 1.00000   1 1.0000
ispG 0.000000 0.00000    NaN 1.00000   1 1.0000

Contrast: freshwater_hot spring 

average      sd  ratio     ava avb cumsum
dxs  0.102827 0.11482 0.8955 3.03571 2.6 0.5184
phaA 0.043950 0.04821 0.9117 0.42857 0.4 0.7400
ispC 0.021129 0.04278 0.4939 1.00000 0.8 0.8465
ispG 0.021129 0.04278 0.4939 1.00000 0.8 0.9531
hmgs 0.004417 0.01625 0.2718 0.07143 0.0 0.9753
ispE 0.002533 0.01351 0.1876 1.03571 1.0 0.9881
ispD 0.002359 0.01254 0.1882 1.03571 1.0 1.0000

Contrast: freshwater_temperate 

average      sd  ratio     ava avb cumsum
dxs  0.107461 0.08057 1.3338 3.03571 4.5 0.7256
phaA 0.032776 0.03264 1.0042 0.42857 0.5 0.9469
hmgs 0.003770 0.01374 0.2745 0.07143 0.0 0.9723
ispE 0.002108 0.01107 0.1904 1.03571 1.0 0.9866
ispD 0.001990 0.01045 0.1904 1.03571 1.0 1.0000
ispC 0.000000 0.00000    NaN 1.00000 1.0 1.0000
ispG 0.000000 0.00000    NaN 1.00000 1.0 1.0000

Contrast: cave_desert 

average      sd  ratio ava     avb cumsum
dxs  0.057372 0.05727 1.0018 3.3 3.03704 0.3741
phaA 0.042737 0.04895 0.8731 0.2 0.66667 0.6528
ispC 0.012632 0.02641 0.4782 1.1 0.96296 0.7352
ispE 0.011114 0.02560 0.4341 0.9 1.00000 0.8077
ispD 0.010869 0.02502 0.4343 0.9 1.00000 0.8786
hmgs 0.010393 0.02404 0.4324 0.1 0.07407 0.9463
ispG 0.008229 0.02146 0.3834 1.1 0.96296 1.0000

Contrast: cave_tropical 

average      sd  ratio ava    avb cumsum
dxs  0.068828 0.06026 1.1421 3.3 3.0000 0.5751
phaA 0.018332 0.02935 0.6247 0.2 0.1429 0.7283
ispD 0.007096 0.02148 0.3303 0.9 1.0000 0.7876
ispE 0.007096 0.02148 0.3303 0.9 1.0000 0.8469
hmgs 0.006624 0.02005 0.3304 0.1 0.0000 0.9023
ispC 0.005847 0.01769 0.3305 1.1 1.0000 0.9511
ispG 0.005847 0.01769 0.3305 1.1 1.0000 1.0000

Contrast: cave_soil 

average      sd  ratio ava avb cumsum
dxs  0.058452 0.05713 1.0231 3.3 2.8 0.4864
phaA 0.029415 0.03379 0.8705 0.2 0.4 0.7311
ispD 0.007048 0.02137 0.3298 0.9 1.0 0.7898
ispE 0.007048 0.02137 0.3298 0.9 1.0 0.8484
hmgs 0.006583 0.01996 0.3299 0.1 0.0 0.9032
ispC 0.005817 0.01763 0.3299 1.1 1.0 0.9516
ispG 0.005817 0.01763 0.3299 1.1 1.0 1.0000

Contrast: cave_rock 

average      sd  ratio ava avb cumsum
dxs  0.048599 0.05931 0.8194 3.3   3 0.5178
phaA 0.012549 0.02652 0.4732 0.2   0 0.6514
ispD 0.007143 0.02259 0.3162 0.9   1 0.7275
ispE 0.007143 0.02259 0.3162 0.9   1 0.8036
hmgs 0.006667 0.02108 0.3162 0.1   0 0.8747
ispC 0.005882 0.01860 0.3162 1.1   1 0.9373
ispG 0.005882 0.01860 0.3162 1.1   1 1.0000

Contrast: cave_hot spring 

average      sd  ratio ava avb cumsum
dxs  0.107816 0.13267 0.8127 3.3 2.6 0.4939
phaA 0.033092 0.04838 0.6840 0.2 0.4 0.6455
ispC 0.027485 0.04778 0.5752 1.1 0.8 0.7714
ispG 0.027485 0.04778 0.5752 1.1 0.8 0.8973
ispD 0.007663 0.02389 0.3208 0.9 1.0 0.9324
ispE 0.007663 0.02389 0.3208 0.9 1.0 0.9675
hmgs 0.007093 0.02201 0.3223 0.1 0.0 1.0000

Contrast: cave_temperate 

average      sd  ratio ava avb cumsum
dxs  0.101761 0.08619 1.1806 3.3 4.5 0.6264
phaA 0.031693 0.03295 0.9618 0.2 0.5 0.8215
ispD 0.006275 0.01935 0.3242 0.9 1.0 0.8601
ispE 0.006275 0.01935 0.3242 0.9 1.0 0.8987
hmgs 0.005903 0.01820 0.3243 0.1 0.0 0.9350
ispC 0.005278 0.01627 0.3244 1.1 1.0 0.9675
ispG 0.005278 0.01627 0.3244 1.1 1.0 1.0000

Contrast: desert_tropical 

average      sd  ratio     ava    avb cumsum
dxs  0.048825 0.04354 1.1214 3.03704 3.0000 0.4156
phaA 0.043828 0.05045 0.8688 0.66667 0.1429 0.7887
ispC 0.007503 0.02153 0.3485 0.96296 1.0000 0.8526
hmgs 0.005081 0.01806 0.2814 0.07407 0.0000 0.8959
ispE 0.004929 0.01754 0.2809 1.00000 1.0000 0.9378
ispD 0.004674 0.01673 0.2793 1.00000 1.0000 0.9776
ispG 0.002628 0.01346 0.1952 0.96296 1.0000 1.0000

Contrast: desert_soil 

average      sd  ratio     ava avb cumsum
phaA 0.043524 0.04701 0.9258 0.66667 0.4 0.4349
dxs  0.031899 0.03825 0.8339 3.03704 2.8 0.7536
ispC 0.007454 0.02138 0.3487 0.96296 1.0 0.8281
hmgs 0.005049 0.01793 0.2815 0.07407 0.0 0.8786
ispE 0.004898 0.01743 0.2810 1.00000 1.0 0.9275
ispD 0.004646 0.01663 0.2794 1.00000 1.0 0.9739
ispG 0.002610 0.01336 0.1953 0.96296 1.0 1.0000

Contrast: desert_rock 

average      sd  ratio     ava avb cumsum
phaA 0.044143 0.05294 0.8339 0.66667   0 0.4796
dxs  0.022923 0.03323 0.6898 3.03704   3 0.7286
ispC 0.007552 0.02198 0.3435 0.96296   1 0.8107
hmgs 0.005115 0.01844 0.2774 0.07407   0 0.8663
ispE 0.004960 0.01791 0.2769 1.00000   1 0.9202
ispD 0.004703 0.01709 0.2752 1.00000   1 0.9713
ispG 0.002646 0.01375 0.1925 0.96296   1 1.0000

Contrast: desert_hot spring 

average      sd  ratio     ava avb cumsum
dxs  0.088349 0.11982 0.7374 3.03704 2.6 0.4250
phaA 0.056344 0.05864 0.9608 0.66667 0.4 0.6960
ispC 0.025423 0.04322 0.5883 0.96296 0.8 0.8183
ispG 0.022033 0.04206 0.5239 0.96296 0.8 0.9243
hmgs 0.005465 0.01991 0.2745 0.07407 0.0 0.9506
ispE 0.005285 0.01927 0.2743 1.00000 1.0 0.9760
ispD 0.004991 0.01830 0.2727 1.00000 1.0 1.0000

Contrast: desert_temperate 

average      sd  ratio     ava avb cumsum
dxs  0.095230 0.08253 1.1539 3.03704 4.5 0.6094
phaA 0.038979 0.04042 0.9644 0.66667 0.5 0.8588
ispC 0.006661 0.01921 0.3468 0.96296 1.0 0.9014
hmgs 0.004510 0.01614 0.2795 0.07407 0.0 0.9303
ispE 0.004388 0.01572 0.2791 1.00000 1.0 0.9584
ispD 0.004180 0.01505 0.2778 1.00000 1.0 0.9851
ispG 0.002324 0.01199 0.1939 0.96296 1.0 1.0000

Contrast: tropical_soil 

average      sd  ratio    ava avb cumsum
dxs  0.04587 0.04085 1.1230 3.0000 2.8 0.6089
phaA 0.02946 0.03463 0.8507 0.1429 0.4 1.0000
hmgs 0.00000 0.00000    NaN 0.0000 0.0 1.0000
ispC 0.00000 0.00000    NaN 1.0000 1.0 1.0000
ispD 0.00000 0.00000    NaN 1.0000 1.0 1.0000
ispE 0.00000 0.00000    NaN 1.0000 1.0 1.0000
ispG 0.00000 0.00000    NaN 1.0000 1.0 1.0000

Contrast: tropical_rock 

average     sd ratio    ava avb cumsum
dxs  0.041026 0.0386 1.063 3.0000   3 0.8116
phaA 0.009524 0.0252 0.378 0.1429   0 1.0000
hmgs 0.000000 0.0000   NaN 0.0000   0 1.0000
ispC 0.000000 0.0000   NaN 1.0000   1 1.0000
ispD 0.000000 0.0000   NaN 1.0000   1 1.0000
ispE 0.000000 0.0000   NaN 1.0000   1 1.0000
ispG 0.000000 0.0000   NaN 1.0000   1 1.0000

Contrast: tropical_hot spring 

average      sd  ratio    ava avb cumsum
dxs  0.10210 0.12166 0.8392 3.0000 2.6 0.5741
phaA 0.03162 0.05022 0.6296 0.1429 0.4 0.7519
ispC 0.02206 0.04502 0.4901 1.0000 0.8 0.8759
ispG 0.02206 0.04502 0.4901 1.0000 0.8 1.0000
hmgs 0.00000 0.00000    NaN 0.0000 0.0 1.0000
ispD 0.00000 0.00000    NaN 1.0000 1.0 1.0000
ispE 0.00000 0.00000    NaN 1.0000 1.0 1.0000

Contrast: tropical_temperate 

average      sd  ratio    ava avb cumsum
dxs  0.10783 0.08547 1.2616 3.0000 4.5 0.7677
phaA 0.03262 0.03407 0.9576 0.1429 0.5 1.0000
hmgs 0.00000 0.00000    NaN 0.0000 0.0 1.0000
ispC 0.00000 0.00000    NaN 1.0000 1.0 1.0000
ispD 0.00000 0.00000    NaN 1.0000 1.0 1.0000
ispE 0.00000 0.00000    NaN 1.0000 1.0 1.0000
ispG 0.00000 0.00000    NaN 1.0000 1.0 1.0000

Contrast: soil_rock 

average      sd  ratio ava avb cumsum
phaA 0.02762 0.03786 0.7296 0.4   0 0.6591
dxs  0.01429 0.03194 0.4472 2.8   3 1.0000
hmgs 0.00000 0.00000    NaN 0.0   0 1.0000
ispC 0.00000 0.00000    NaN 1.0   1 1.0000
ispD 0.00000 0.00000    NaN 1.0   1 1.0000
ispE 0.00000 0.00000    NaN 1.0   1 1.0000
ispG 0.00000 0.00000    NaN 1.0   1 1.0000

Contrast: soil_hot spring 

average      sd  ratio ava avb cumsum
dxs  0.08494 0.11928 0.7121 2.8 2.6 0.4910
phaA 0.04451 0.04765 0.9341 0.4 0.4 0.7482
ispC 0.02178 0.04450 0.4894 1.0 0.8 0.8741
ispG 0.02178 0.04450 0.4894 1.0 0.8 1.0000
hmgs 0.00000 0.00000    NaN 0.0 0.0 1.0000
ispD 0.00000 0.00000    NaN 1.0 1.0 1.0000
ispE 0.00000 0.00000    NaN 1.0 1.0 1.0000

Contrast: soil_temperate 

average      sd  ratio ava avb cumsum
dxs  0.09980 0.09510 1.0494 2.8 4.5 0.7605
phaA 0.03144 0.03333 0.9433 0.4 0.5 1.0000
hmgs 0.00000 0.00000    NaN 0.0 0.0 1.0000
ispC 0.00000 0.00000    NaN 1.0 1.0 1.0000
ispD 0.00000 0.00000    NaN 1.0 1.0 1.0000
ispE 0.00000 0.00000    NaN 1.0 1.0 1.0000
ispG 0.00000 0.00000    NaN 1.0 1.0 1.0000

Contrast: rock_hot spring 

average      sd  ratio ava avb cumsum
dxs  0.08000 0.14453 0.5535   3 2.6 0.5353
phaA 0.02500 0.05590 0.4472   0 0.4 0.7026
ispC 0.02222 0.04969 0.4472   1 0.8 0.8513
ispG 0.02222 0.04969 0.4472   1 0.8 1.0000
hmgs 0.00000 0.00000    NaN   0 0.0 1.0000
ispD 0.00000 0.00000    NaN   1 1.0 1.0000
ispE 0.00000 0.00000    NaN   1 1.0 1.0000

Contrast: rock_temperate 

average      sd  ratio ava avb cumsum
dxs  0.08824 0.12478 0.7071   3 4.5 0.7258
phaA 0.03333 0.04714 0.7071   0 0.5 1.0000
hmgs 0.00000 0.00000    NaN   0 0.0 1.0000
ispC 0.00000 0.00000    NaN   1 1.0 1.0000
ispD 0.00000 0.00000    NaN   1 1.0 1.0000
ispE 0.00000 0.00000    NaN   1 1.0 1.0000
ispG 0.00000 0.00000    NaN   1 1.0 1.0000

Contrast: hot spring_temperate 

average      sd  ratio ava avb cumsum
dxs  0.14844 0.15745 0.9428 2.6 4.5 0.6423
phaA 0.04599 0.04240 1.0846 0.4 0.5 0.8413
ispC 0.01833 0.03885 0.4719 0.8 1.0 0.9207
ispG 0.01833 0.03885 0.4719 0.8 1.0 1.0000
hmgs 0.00000 0.00000    NaN 0.0 0.0 1.0000
ispD 0.00000 0.00000    NaN 1.0 1.0 1.0000
ispE 0.00000 0.00000    NaN 1.0 1.0 1.0000
Permutation: free
Number of permutations: 0
