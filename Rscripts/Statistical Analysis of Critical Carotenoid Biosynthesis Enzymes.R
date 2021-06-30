---
title: "Statistical Analysis of Critical Carotenoid Biosynthesis Enzymes"
author: "Dionne Martin"
date: "6/25/2021"
output: pdf_document
---
  
#################### PREFACE ##########################
#Code for PCA, NMDS, and PERMANOVA/ simper test 
#of raw BLASTp hits of critical carotenoid 
#biosynthesis enzymes: crtE, crtIa, crtIb, crtBa, 
#crtBb, cruA, cruP, and crtL. 
#these enzymes are necesary for the production of 
# beta carotene for all downstream carotenoids 
#######################################################

```{pca}
#load packages and blastp raw hits data 
library(devtools)
library(ggbiplot)
library(readr)
small <- read_csv("/bigdata/stajichlab/shared/projects/Cyanobacteria/Terresterial_Cyano/Cyano Carotenoids BLASTp Results - small dataset (1).csv")
#organize data to just raw hits and strain names 
row.names(small) = small$Strain
small_matrix <- data.matrix(small[c(3:10)])
pca= prcomp(small_matrix, center= TRUE, scale.= TRUE)


#analyze pca statistics to find proportion of variatiance explained by each gene 
summary(pca)
str(pca)
#crtE-0.4522 crtBa-0.1573 crtBb-0.1263 crtIa-0.1196
#crtIb-0.096 cruA-0.036 cruP-0.009 crtL-0.002

#make groups based on habitat
habitat= small$Habitat

#try to make group colours, this should be added to plot code 
scale_colour_manual(name= "Habitat", values= c("brown", "dark orange", "blue", "pink", "red", "dark blue", "forest green", "dark yellow", "red", "black"))

#make the PCA plot with strains as points 
png("PCA Plot of Carotenoid BLAST Hit .png", width=1000, height=1000)
ggbiplot(pca,ellipse=FALSE,obs.scale = 3,var.scale = 1, groups=hits.habitat) +
  scale_colour_manual(name="Habitat", values= c("saddlebrown", "darkorange1", "cornflowerblue", "deeppink", "red", "dark blue", "forest green", "gold", "darkmagenta", "black"))+
  ggtitle("PCA of Raw BLASTp Hits of Critical Enzymes")+
  theme_minimal()+
  theme(legend.position = "bottom")

#make the PCA plot with strain names
png("PCA Plot of Carotenoid BLAST Hit .png", width=1000, height=1000)
ggbiplot(pca,labels= rownames(small),ellipse=FALSE,obs.scale = 3,var.scale = 1, groups=hits.habitat) +
  scale_colour_manual(name="Habitat", values= c("saddlebrown", "darkorange1", "cornflowerblue", "deeppink", "red", "dark blue", "forest green", "gold", "darkmagenta", "black"))+
  ggtitle("PCA of Raw BLASTp Hits of Critical Enzymes")+
  theme_minimal()+
  theme(legend.position = "bottom")

```
```{nmds}
#load packages 
library(vegan)
library(ggplot2)
library(readr)

#prep data
small <- read_csv("/bigdata/stajichlab/shared/projects/Cyanobacteria/Terresterial_Cyano/Cyano Carotenoids BLASTp Results - small dataset (1).csv")
small_matrix <- subset(small, select = c(crtE, crtBa, crtBb, crtIa, crtIb, cruA, cruP, crtL))
nmds_small_habitat= as.data.frame(small[,1:2])
#Do you need to transform your data? range should be 0-10 if not, use root function
range(small_matrix) # 0-7 no transformation needed

#make a distance matrix
small_matrix_ds=vegdist(small_matrix, method = 'euclidian')

#Run metaMDS on distance matrix
small_nms=metaMDS(small_matrix_ds, distance="euclidian",k=2,trymax=500)
small_nms
#stress=0.091

#assign colors to habitat
cols <- nmds_small_habitat$Strain
cols[] <- c("saddlebrown", "darkorange1", "cornflowerblue", "deeppink", "red", "dark blue", "forest green", "gold", "darkmagenta", "black")
#Warning message:
#number of items to replace is not a multiple of replacement length

cols[nmds_small_habitat$Habitat == c("freshwater")] <- "cornflowerblue"
cols[nmds_small_habitat$Habitat == c("saltwater")] <- "dark blue"
cols[nmds_small_habitat$Habitat == c("hot spring")] <- "deeppink"
cols[nmds_small_habitat$Habitat == c("soil")] <- "forest green"
cols[nmds_small_habitat$Habitat == c("cave")] <- "saddlebrown"
cols[nmds_small_habitat$Habitat == c("desert")] <- "darkorange1"
cols[nmds_small_habitat$Habitat == c("tropical")] <- "black"
cols[nmds_small_habitat$Habitat == c("temperate")] <- "darkmagenta"
cols[nmds_small_habitat$Habitat == c("symbiont")] <- "gold"
cols[nmds_small_habitat$Habitat == c("rock")] <- "red"
pairs(character(nrow(nmds_small_habitat$Strain, col=cols)

#habitatvec = factor(small$Habitat, levels = c("freshwater", "saltwater", "hot spring", "soil", "cave", "desert", "tropical", "temperate", "symbiont", "rock"))
#colors <- c("cornflowerblue", "dark blue", "deeppink", "forest green", "saddlebrown", "darkorange1", "black", "darkmagenta", "gold", "red")
habitatinfo= c(rep(small$Habitat))

#Correct Plot With Species Names 
ordiplot(small_nms, type="n")
orditorp(small_nms,display="sites",labels=small$Strain, col=cols ,air=0.00025)
par(xpd=TRUE)
#legend(0.6,0.25,legend =unique(nmds_small_habitat$Strain), fill=cols, border= "black")
#change coordinates of the legend for each specific graph
#legend(-0.4,-0.2, "stress = 0.13") 
#change coordinates and stress value for each specific graph
title(main="NMDS of Raw BLASTp")
scores <-
  scores(small_nms, display = "sites", "species")
cent <-
  aggregate(scores ~ habitatinfo, data = nmds_small_habitat, FUN = "mean")
names(cent) [-1] <- colnames(scores)
points(cent [,-1],
       pch = c( 8 , 8 , 8, 8),
       col = unique(colors),
       bg = c("black"),
       lwd = 3.0,
       cex = 2.0
)

#NMDS Plot 
plot(small_nms)
points(small_nms,display= "sites",pch=21, bg=cols) 
par(xpd=TRUE)
#legend(0.6,0.25,legend =habitatvec, fill=colvec, border= "black")
#legend(-0.4,-0.2, "stress = 0.13") 
title(main="NMDS of Raw BLASTp")
scores <-
  scores(small_nms, display = "sites", "species")
cent <-
  aggregate(scores ~ habitatinfo, data = nmds_small_habitat, FUN = "mean")
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
rownames(small_matrix) = small$Strain
#habitat= small$Habitat
str(small_matrix)
nmds_small_habitat= as.data.frame(small[,1:2])
str(nmds_small_habitat)
#habitatinfo= c(rep(small$Habitat))

#make a distance matrix
small_matrix_ds=vegdist(small_matrix, method = 'euclidian')

#Run metaMDS on distance matrix
small_nms=metaMDS(small_matrix_ds, distance="euclidian",k=2,trymax=500)

#Run PERMANOVA
pmv_small = adonis(small_matrix ~ habitat, data = nmds_small_habitat, permutations = 999, method = 'euclidean')
pmv_small
#PERMANOVA Results: Based off of the carotenoid genes being compared
#cyanobacterial strains differ statistically significant 
#between their habitats (p-value:0.001).
#Habitats explian 31.91% of the variance in the copy numbers of critical 
#carotenoid biosynthesis enzymes between cyanobacterial species. 


#Visualise the permanova
densityplot(permustats(pmv_small))

#Distance Based Dispersion test (supposed to be non-euclidean distance)

bd_small = betadisper(small_matrix_ds, small$Habitat)
boxplot(bd_small)
anova(bd_small) 
#F-Test: pvalue 0.5072
permutest(bd_small)
#permutation test: pvalue 0.49

#both the F-Test and permutation test show that 
#the avergae distance to the median for the 
#habitat data is not statistically significant. 

#simper test: determine which genes are responsible for the difference between 2 compared habiatat groups 
simp_small= simper(small_matrix, group= small$Habitat)
summary(simp_small)
#Results of simper test
Contrast: saltwater_symbiont 

average      sd  ratio   ava   avb cumsum
cruA  0.05134 0.03026 1.6965 0.625 1.667 0.1890
cruP  0.05134 0.03026 1.6965 0.625 1.667 0.3779
crtL  0.04251 0.03318 1.2813 1.083 0.000 0.5344
crtBb 0.04095 0.03093 1.3240 2.333 2.333 0.6851
crtIa 0.03512 0.03234 1.0860 2.542 3.000 0.8143
crtIb 0.02329 0.02362 0.9859 2.583 2.667 0.9000
crtE  0.01624 0.02132 0.7618 2.083 1.667 0.9598
crtBa 0.01092 0.01568 0.6965 1.000 1.333 1.0000

Contrast: saltwater_freshwater 

average      sd  ratio   ava    avb cumsum
cruP  0.05099 0.03211 1.5879 0.625 1.8571 0.2012
cruA  0.05024 0.03124 1.6079 0.625 1.8214 0.3995
crtL  0.03757 0.02928 1.2830 1.083 0.1071 0.5478
crtIb 0.03465 0.02914 1.1892 2.583 3.1786 0.6846
crtIa 0.02980 0.02661 1.1200 2.542 2.9643 0.8022
crtBb 0.02931 0.02607 1.1243 2.333 2.7143 0.9179
crtBa 0.01310 0.01657 0.7905 1.000 1.3929 0.9696
crtE  0.00771 0.01430 0.5391 2.083 2.1786 1.0000

Contrast: saltwater_cave 

average      sd  ratio   ava avb cumsum
cruA  0.05424 0.03318 1.6344 0.625 1.9 0.1779
cruP  0.05424 0.03318 1.6344 0.625 1.9 0.3558
crtIb 0.05029 0.03869 1.3000 2.583 3.1 0.5208
crtIa 0.04379 0.03541 1.2368 2.542 3.1 0.6645
crtL  0.04128 0.03221 1.2816 1.083 0.0 0.7999
crtBb 0.03488 0.02659 1.3120 2.333 1.5 0.9143
crtBa 0.01442 0.01818 0.7935 1.000 1.4 0.9616
crtE  0.01169 0.01686 0.6934 2.083 2.3 1.0000

Contrast: saltwater_desert 

average      sd  ratio   ava    avb cumsum
cruA  0.04899 0.03248 1.5080 0.625 1.5926 0.1699
cruP  0.04899 0.03248 1.5080 0.625 1.5926 0.3399
crtBb 0.04644 0.02917 1.5923 2.333 1.0741 0.5010
crtIb 0.04245 0.03499 1.2130 2.583 3.2222 0.6483
crtIa 0.04166 0.03492 1.1930 2.542 3.1852 0.7928
crtL  0.04120 0.03320 1.2408 1.083 0.1111 0.9358
crtBa 0.01260 0.01842 0.6838 1.000 1.2593 0.9795
crtE  0.00592 0.01296 0.4568 2.083 2.1111 1.0000

Contrast: saltwater_tropical 

average      sd  ratio   ava   avb cumsum
cruA  0.057027 0.03521 1.6194 0.625 2.000 0.1965
cruP  0.057027 0.03521 1.6194 0.625 2.000 0.3930
crtBb 0.046377 0.02194 2.1138 2.333 1.000 0.5529
crtL  0.041774 0.03236 1.2907 1.083 0.000 0.6968
crtIb 0.040161 0.03342 1.2017 2.583 3.286 0.8352
crtIa 0.036252 0.03278 1.1059 2.542 3.143 0.9601
crtE  0.007041 0.01417 0.4968 2.083 2.143 0.9844
crtBa 0.004530 0.01119 0.4046 1.000 1.143 1.0000

Contrast: saltwater_soil 

average      sd  ratio   ava avb cumsum
cruP  0.056070 0.03384 1.6568 0.625 2.2 0.1876
cruA  0.055895 0.03342 1.6727 0.625 2.2 0.3746
crtIb 0.054464 0.03409 1.5976 2.583 4.0 0.5568
crtIa 0.049066 0.03254 1.5077 2.542 3.8 0.7210
crtBb 0.039036 0.02671 1.4615 2.333 2.8 0.8515
crtL  0.036946 0.02876 1.2848 1.083 0.0 0.9752
crtE  0.007427 0.01300 0.5716 2.083 2.2 1.0000
crtBa 0.000000 0.00000    NaN 1.000 1.0 1.0000

Contrast: saltwater_rock 

average       sd  ratio   ava avb cumsum
crtIa 0.134593 0.033283 4.0439 2.542   7 0.3345
crtIb 0.103491 0.030331 3.4121 2.583   6 0.5917
cruA  0.045459 0.028034 1.6215 0.625   2 0.7047
cruP  0.045459 0.028034 1.6215 0.625   2 0.8177
crtBb 0.037783 0.018829 2.0067 2.333   1 0.9116
crtL  0.033413 0.026014 1.2844 1.083   0 0.9947
crtE  0.002142 0.007269 0.2947 2.083   2 1.0000
crtBa 0.000000 0.000000    NaN 1.000   1 1.0000

Contrast: saltwater_hot spring 

average      sd  ratio   ava avb cumsum
crtBb 0.04525 0.04082 1.1085 2.333 3.4 0.1649
crtIa 0.04463 0.02862 1.5594 2.542 3.6 0.3276
cruA  0.04417 0.03348 1.3195 0.625 1.6 0.4886
cruP  0.04417 0.03348 1.3195 0.625 1.6 0.6496
crtL  0.04003 0.03247 1.2325 1.083 0.0 0.7955
crtIb 0.03787 0.02763 1.3706 2.583 3.4 0.9335
crtE  0.01229 0.02173 0.5655 2.083 1.8 0.9783
crtBa 0.00596 0.01204 0.4951 1.000 1.2 1.0000

Contrast: saltwater_temperate 

average      sd ratio   ava avb cumsum
cruA  0.05512 0.03136 1.758 0.625 1.5 0.1992
cruP  0.05512 0.03136 1.758 0.625 1.5 0.3983
crtBb 0.05302 0.02403 2.206 2.333 1.0 0.5899
crtL  0.03839 0.03201 1.199 1.083 0.5 0.7286
crtIb 0.02639 0.02388 1.105 2.583 2.0 0.8240
crtIa 0.02449 0.02410 1.016 2.542 2.0 0.9125
crtE  0.02423 0.02352 1.030 2.083 1.5 1.0000
crtBa 0.00000 0.00000   NaN 1.000 1.0 1.0000

Contrast: symbiont_freshwater 

average      sd  ratio   ava    avb cumsum
crtBb 0.040777 0.02688 1.5168 2.333 2.7143 0.2468
crtIa 0.030030 0.02592 1.1585 3.000 2.9643 0.4285
crtIb 0.028859 0.02359 1.2234 2.667 3.1786 0.6031
crtE  0.017720 0.02149 0.8246 1.667 2.1786 0.7103
crtBa 0.014866 0.01621 0.9169 1.333 1.3929 0.8003
cruA  0.014628 0.01959 0.7465 1.667 1.8214 0.8888
cruP  0.014211 0.01952 0.7280 1.667 1.8571 0.9748
crtL  0.004162 0.01215 0.3425 0.000 0.1071 1.0000

Contrast: symbiont_cave 

average      sd  ratio   ava avb cumsum
crtIb 0.04492 0.03503 1.2824 2.667 3.1 0.2323
crtIa 0.04400 0.03248 1.3550 3.000 3.1 0.4599
crtBb 0.03994 0.03389 1.1786 2.333 1.5 0.6665
crtE  0.02247 0.02426 0.9262 1.667 2.3 0.7827
crtBa 0.01574 0.01744 0.9028 1.333 1.4 0.8642
cruA  0.01313 0.01772 0.7410 1.667 1.9 0.9321
cruP  0.01313 0.01772 0.7410 1.667 1.9 1.0000
crtL  0.00000 0.00000    NaN 0.000 0.0 1.0000

Contrast: symbiont_desert 

average      sd  ratio   ava    avb cumsum
crtBb 0.045238 0.04065 1.1129 2.333 1.0741 0.2251
crtIa 0.039015 0.03343 1.1671 3.000 3.1852 0.4193
crtIb 0.036195 0.02932 1.2346 2.667 3.2222 0.5994
cruA  0.021585 0.02754 0.7839 1.667 1.5926 0.7068
cruP  0.021585 0.02754 0.7839 1.667 1.5926 0.8143
crtE  0.016628 0.02143 0.7760 1.667 2.1111 0.8970
crtBa 0.016459 0.01925 0.8548 1.333 1.2593 0.9789
crtL  0.004236 0.01224 0.3460 0.000 0.1111 1.0000

Contrast: symbiont_tropical 

average      sd  ratio   ava   avb cumsum
crtBb 0.04327 0.03897 1.1105 2.333 1.000 0.2614
crtIa 0.03479 0.02914 1.1939 3.000 3.143 0.4715
crtIb 0.03385 0.02760 1.2265 2.667 3.286 0.6759
crtE  0.01747 0.02226 0.7849 1.667 2.143 0.7815
cruA  0.01210 0.01762 0.6869 1.667 2.000 0.8546
cruP  0.01210 0.01762 0.6869 1.667 2.000 0.9277
crtBa 0.01197 0.01566 0.7644 1.333 1.143 1.0000
crtL  0.00000 0.00000    NaN 0.000 0.000 1.0000

Contrast: symbiont_soil 

average      sd  ratio   ava avb cumsum
crtIb 0.049082 0.02663 1.8429 2.667 4.0 0.2599
crtBb 0.042995 0.03351 1.2830 2.333 2.8 0.4875
crtIa 0.037573 0.03132 1.1998 3.000 3.8 0.6865
crtE  0.016833 0.02013 0.8362 1.667 2.2 0.7756
cruP  0.016637 0.01991 0.8354 1.667 2.2 0.8637
cruA  0.016469 0.01961 0.8400 1.667 2.2 0.9509
crtBa 0.009268 0.01364 0.6795 1.333 1.0 1.0000
crtL  0.000000 0.00000    NaN 0.000 0.0 1.0000

Contrast: symbiont_rock 

average      sd  ratio   ava avb cumsum
crtIa 0.114105 0.03318 3.4389 3.000   7 0.4024
crtIb 0.095457 0.02301 4.1484 2.667   6 0.7390
crtBb 0.035742 0.03875 0.9224 2.333   1 0.8651
crtE  0.010101 0.01750 0.5774 1.667   2 0.9007
cruA  0.009804 0.01698 0.5774 1.667   2 0.9353
cruP  0.009804 0.01698 0.5774 1.667   2 0.9699
crtBa 0.008547 0.01480 0.5774 1.333   1 1.0000
crtL  0.000000 0.00000    NaN 0.000   0 1.0000

Contrast: symbiont_hot spring 

average      sd  ratio   ava avb cumsum
crtBb 0.05460 0.03883 1.4062 2.333 3.4 0.2818
crtIa 0.03360 0.03131 1.0732 3.000 3.6 0.4552
crtIb 0.03222 0.02142 1.5044 2.667 3.4 0.6214
cruA  0.02340 0.03161 0.7404 1.667 1.6 0.7422
cruP  0.02340 0.03161 0.7404 1.667 1.6 0.8629
crtE  0.01445 0.01870 0.7723 1.667 1.8 0.9375
crtBa 0.01211 0.01556 0.7782 1.333 1.2 1.0000
crtL  0.00000 0.00000    NaN 0.000 0.0 1.0000

Contrast: symbiont_temperate 

average      sd  ratio   ava avb cumsum
crtBb 0.04898 0.04646 1.0542 2.333 1.0 0.2390
crtIa 0.03927 0.03745 1.0486 3.000 2.0 0.4307
crtIb 0.02538 0.01992 1.2741 2.667 2.0 0.5546
crtE  0.01994 0.02205 0.9043 1.667 1.5 0.6519
cruA  0.01994 0.02205 0.9043 1.667 1.5 0.7493
cruP  0.01994 0.02205 0.9043 1.667 1.5 0.8466
crtL  0.01994 0.02205 0.9043 0.000 0.5 0.9439
crtBa 0.01149 0.01781 0.6455 1.333 1.0 1.0000

Contrast: freshwater_cave 

average      sd  ratio    ava avb cumsum
crtIb 0.048607 0.03353 1.4495 3.1786 3.1 0.2726
crtBb 0.041448 0.02612 1.5868 2.7143 1.5 0.5051
crtIa 0.039780 0.02976 1.3369 2.9643 3.1 0.7283
crtBa 0.015352 0.01635 0.9390 1.3929 1.4 0.8144
crtE  0.011811 0.01562 0.7562 2.1786 2.3 0.8806
cruA  0.009134 0.01873 0.4878 1.8214 1.9 0.9319
cruP  0.008102 0.01818 0.4455 1.8571 1.9 0.9773
crtL  0.004043 0.01178 0.3433 0.1071 0.0 1.0000

Contrast: freshwater_desert 

average      sd  ratio    ava    avb cumsum
crtBb 0.054301 0.02787 1.9487 2.7143 1.0741 0.2745
crtIb 0.038013 0.03118 1.2190 3.1786 3.2222 0.4666
crtIa 0.035720 0.02949 1.2113 2.9643 3.1852 0.6472
cruA  0.019147 0.02843 0.6735 1.8214 1.5926 0.7439
cruP  0.018733 0.02893 0.6476 1.8571 1.5926 0.8386
crtBa 0.016649 0.01888 0.8820 1.3929 1.2593 0.9228
crtE  0.008087 0.01428 0.5665 2.1786 2.1111 0.9636
crtL  0.007192 0.01490 0.4826 0.1071 0.1111 1.0000

Contrast: freshwater_tropical 

average      sd  ratio    ava   avb cumsum
crtBb 0.054315 0.02162 2.5119 2.7143 1.000 0.3411
crtIb 0.035137 0.02809 1.2510 3.1786 3.286 0.5618
crtIa 0.031367 0.02660 1.1790 2.9643 3.143 0.7588
crtBa 0.013270 0.01566 0.8471 1.3929 1.143 0.8422
crtE  0.008710 0.01449 0.6011 2.1786 2.143 0.8969
cruA  0.006837 0.01814 0.3769 1.8214 2.000 0.9398
cruP  0.005491 0.01716 0.3200 1.8571 2.000 0.9743
crtL  0.004091 0.01189 0.3441 0.1071 0.000 1.0000

Contrast: freshwater_soil 

average      sd  ratio    ava avb cumsum
crtIb 0.040489 0.03017 1.3421 3.1786 4.0 0.2582
crtIa 0.037819 0.02717 1.3922 2.9643 3.8 0.4994
crtBb 0.033250 0.02380 1.3973 2.7143 2.8 0.7115
cruA  0.011459 0.01948 0.5881 1.8214 2.2 0.7846
crtBa 0.011077 0.01398 0.7925 1.3929 1.0 0.8553
cruP  0.010417 0.01899 0.5484 1.8571 2.2 0.9217
crtE  0.008655 0.01318 0.6567 2.1786 2.2 0.9769
crtL  0.003623 0.01056 0.3430 0.1071 0.0 1.0000

Contrast: freshwater_rock 

average       sd  ratio    ava avb cumsum
crtIa 0.110690 0.030491 3.6303 2.9643   7 0.4227
crtIb 0.078035 0.031179 2.5028 3.1786   6 0.7207
crtBb 0.045052 0.018474 2.4387 2.7143   1 0.8928
crtBa 0.010197 0.012965 0.7865 1.3929   1 0.9317
cruA  0.005479 0.014683 0.3732 1.8214   2 0.9527
crtE  0.004717 0.010364 0.4551 2.1786   2 0.9707
cruP  0.004397 0.013880 0.3168 1.8571   2 0.9875
crtL  0.003281 0.009645 0.3401 0.1071   0 1.0000

Contrast: freshwater_hot spring 

average      sd  ratio    ava avb cumsum
crtBb 0.03993 0.03050 1.3091 2.7143 3.4 0.2289
crtIa 0.03412 0.02490 1.3704 2.9643 3.6 0.4245
crtIb 0.03077 0.02647 1.1623 3.1786 3.4 0.6009
cruA  0.01966 0.03206 0.6133 1.8214 1.6 0.7136
cruP  0.01908 0.03276 0.5825 1.8571 1.6 0.8230
crtE  0.01385 0.02153 0.6436 2.1786 1.8 0.9024
crtBa 0.01311 0.01536 0.8535 1.3929 1.2 0.9775
crtL  0.00392 0.01162 0.3375 0.1071 0.0 1.0000

Contrast: freshwater_temperate 

average      sd  ratio    ava avb cumsum
crtBb 0.06129 0.02360 2.5965 2.7143 1.0 0.2637
crtIb 0.04046 0.03072 1.3173 3.1786 2.0 0.4378
crtIa 0.03284 0.02671 1.2295 2.9643 2.0 0.5791
crtE  0.02502 0.02358 1.0612 2.1786 1.5 0.6868
cruA  0.02020 0.02111 0.9570 1.8214 1.5 0.7737
cruP  0.02020 0.02111 0.9570 1.8571 1.5 0.8606
crtL  0.01858 0.01896 0.9799 0.1071 0.5 0.9406
crtBa 0.01381 0.01745 0.7910 1.3929 1.0 1.0000

Contrast: cave_desert 

average      sd  ratio ava    avb cumsum
crtIb 0.054774 0.03957 1.3842 3.1 3.2222 0.2845
crtIa 0.047669 0.03654 1.3045 3.1 3.1852 0.5321
crtBb 0.019855 0.02515 0.7893 1.5 1.0741 0.6352
cruA  0.018397 0.02941 0.6256 1.9 1.5926 0.7308
cruP  0.018397 0.02941 0.6256 1.9 1.5926 0.8263
crtBa 0.017658 0.02030 0.8700 1.4 1.2593 0.9180
crtE  0.011668 0.01644 0.7095 2.3 2.1111 0.9786
crtL  0.004115 0.01186 0.3469 0.0 0.1111 1.0000

Contrast: cave_tropical 

average      sd  ratio ava   avb cumsum
crtIb 0.052589 0.03671 1.4325 3.1 3.286 0.3583
crtIa 0.043620 0.03377 1.2916 3.1 3.143 0.6555
crtBb 0.016839 0.02267 0.7426 1.5 1.000 0.7702
crtBa 0.014294 0.01691 0.8451 1.4 1.143 0.8676
crtE  0.011900 0.01625 0.7322 2.3 2.143 0.9486
cruA  0.003769 0.01143 0.3297 1.9 2.000 0.9743
cruP  0.003769 0.01143 0.3297 1.9 2.000 1.0000
crtL  0.000000 0.00000    NaN 0.0 0.000 1.0000

Contrast: cave_soil 

average      sd  ratio ava avb cumsum
crtIb 0.055910 0.03797 1.4724 3.1 4.0 0.2948
crtIa 0.047233 0.03232 1.4612 3.1 3.8 0.5438
crtBb 0.045332 0.03079 1.4724 1.5 2.8 0.7828
crtBa 0.012023 0.01514 0.7940 1.4 1.0 0.8461
crtE  0.011205 0.01462 0.7664 2.3 2.2 0.9052
cruP  0.009070 0.01548 0.5857 1.9 2.2 0.9530
cruA  0.008909 0.01520 0.5863 1.9 2.2 1.0000
crtL  0.000000 0.00000    NaN 0.0 0.0 1.0000

Contrast: cave_rock 

average       sd  ratio ava avb cumsum
crtIa 0.111544 0.048109 2.3186 3.1   7 0.4749
crtIb 0.084304 0.053257 1.5830 3.1   6 0.8338
crtBb 0.013816 0.019411 0.7117 1.5   1 0.8926
crtBa 0.011000 0.014324 0.7679 1.4   1 0.9395
crtE  0.008158 0.013191 0.6185 2.3   2 0.9742
cruA  0.003030 0.009583 0.3162 1.9   2 0.9871
cruP  0.003030 0.009583 0.3162 1.9   2 1.0000
crtL  0.000000 0.000000    NaN 0.0   0 1.0000

Contrast: cave_hot spring 

average      sd  ratio ava avb cumsum
crtBb 0.05933 0.03908 1.5184 1.5 3.4 0.2612
crtIb 0.05024 0.03233 1.5541 3.1 3.4 0.4823
crtIa 0.04695 0.02887 1.6262 3.1 3.6 0.6890
cruA  0.01917 0.03403 0.5632 1.9 1.6 0.7733
cruP  0.01917 0.03403 0.5632 1.9 1.6 0.8577
crtE  0.01829 0.02450 0.7467 2.3 1.8 0.9382
crtBa 0.01404 0.01663 0.8441 1.4 1.2 1.0000
crtL  0.00000 0.00000    NaN 0.0 0.0 1.0000

Contrast: cave_temperate 

average      sd  ratio ava avb cumsum
crtIb 0.05317 0.04137 1.2852 3.1 2.0 0.2397
crtIa 0.04560 0.03916 1.1645 3.1 2.0 0.4453
crtE  0.03059 0.02665 1.1482 2.3 1.5 0.5832
cruA  0.01936 0.02012 0.9619 1.9 1.5 0.6704
cruP  0.01936 0.02012 0.9619 1.9 1.5 0.7577
crtL  0.01936 0.02012 0.9619 0.0 0.5 0.8449
crtBb 0.01915 0.02613 0.7327 1.5 1.0 0.9312
crtBa 0.01525 0.01948 0.7831 1.4 1.0 1.0000

Contrast: desert_tropical 

average      sd  ratio    ava   avb cumsum
crtIb 0.041895 0.03425 1.2231 3.2222 3.286 0.2841
crtIa 0.039961 0.03305 1.2090 3.1852 3.143 0.5551
cruA  0.017452 0.03095 0.5639 1.5926 2.000 0.6735
cruP  0.017452 0.03095 0.5639 1.5926 2.000 0.7918
crtBa 0.013570 0.01783 0.7611 1.2593 1.143 0.8838
crtE  0.007543 0.01430 0.5274 2.1111 2.143 0.9350
crtBb 0.005422 0.01340 0.4047 1.0741 1.000 0.9718
crtL  0.004164 0.01197 0.3478 0.1111 0.000 1.0000

Contrast: desert_soil 

average      sd  ratio    ava avb cumsum
crtBb 0.054524 0.03601 1.5140 1.0741 2.8 0.2650
crtIb 0.045488 0.03666 1.2408 3.2222 4.0 0.4861
crtIa 0.041729 0.03393 1.2297 3.1852 3.8 0.6889
cruP  0.021152 0.02978 0.7103 1.5926 2.2 0.7917
cruA  0.020980 0.02950 0.7111 1.5926 2.2 0.8937
crtBa 0.010403 0.01507 0.6905 1.2593 1.0 0.9442
crtE  0.007785 0.01306 0.5961 2.1111 2.2 0.9821
crtL  0.003689 0.01062 0.3473 0.1111 0.0 1.0000

Contrast: desert_rock 

average       sd  ratio    ava avb cumsum
crtIa 0.112468 0.043884 2.5629 3.1852   7 0.4639
crtIb 0.082614 0.041258 2.0024 3.2222   6 0.8046
cruA  0.013619 0.024100 0.5651 1.5926   2 0.8608
cruP  0.013619 0.024100 0.5651 1.5926   2 0.9170
crtBa 0.009479 0.013796 0.6871 1.2593   1 0.9561
crtBb 0.004357 0.010797 0.4035 1.0741   1 0.9740
crtL  0.003342 0.009695 0.3447 0.1111   0 0.9878
crtE  0.002956 0.008529 0.3466 2.1111   2 1.0000

Contrast: desert_hot spring 

average      sd  ratio    ava avb cumsum
crtBb 0.073595 0.04081 1.8033 1.0741 3.4 0.3168
crtIa 0.039707 0.03283 1.2096 3.1852 3.6 0.4877
crtIb 0.037457 0.03220 1.1633 3.2222 3.4 0.6489
cruA  0.025576 0.03413 0.7493 1.5926 1.6 0.7590
cruP  0.025576 0.03413 0.7493 1.5926 1.6 0.8691
crtBa 0.013698 0.01767 0.7750 1.2593 1.2 0.9280
crtE  0.012728 0.02174 0.5853 2.1111 1.8 0.9828
crtL  0.003992 0.01172 0.3407 0.1111 0.0 1.0000

Contrast: desert_temperate 

average      sd  ratio    ava avb cumsum
crtIb 0.048374 0.03616 1.3377 3.2222 2.0 0.2276
crtIa 0.046822 0.03548 1.3197 3.1852 2.0 0.4478
cruA  0.026534 0.02929 0.9060 1.5926 1.5 0.5727
cruP  0.026534 0.02929 0.9060 1.5926 1.5 0.6975
crtE  0.024356 0.02364 1.0304 2.1111 1.5 0.8121
crtL  0.020325 0.02105 0.9655 0.1111 0.5 0.9077
crtBa 0.013353 0.01949 0.6850 1.2593 1.0 0.9705
crtBb 0.006266 0.01564 0.4006 1.0741 1.0 1.0000

Contrast: tropical_soil 

average      sd  ratio   ava avb cumsum
crtBb 0.053572 0.03364 1.5926 1.000 2.8 0.3369
crtIb 0.041919 0.03227 1.2990 3.286 4.0 0.6005
crtIa 0.039895 0.02810 1.4198 3.143 3.8 0.8513
crtE  0.008371 0.01351 0.6198 2.143 2.2 0.9040
cruP  0.005785 0.01177 0.4916 2.000 2.2 0.9404
cruA  0.005622 0.01143 0.4916 2.000 2.2 0.9757
crtBa 0.003864 0.00964 0.4008 1.143 1.0 1.0000
crtL  0.000000 0.00000    NaN 0.000 0.0 1.0000

Contrast: tropical_rock 

average       sd ratio   ava avb cumsum
crtIa 0.110088 0.036295 3.033 3.143   7 0.5627
crtIb 0.078012 0.035844 2.176 3.286   6 0.9615
crtE  0.003968 0.010499 0.378 2.143   2 0.9817
crtBa 0.003571 0.009449 0.378 1.143   1 1.0000
crtBb 0.000000 0.000000   NaN 1.000   1 1.0000
cruA  0.000000 0.000000   NaN 2.000   2 1.0000
cruP  0.000000 0.000000   NaN 2.000   2 1.0000
crtL  0.000000 0.000000   NaN 0.000   0 1.0000

Contrast: tropical_hot spring 

average      sd  ratio   ava avb cumsum
crtBb 0.073345 0.03541 2.0712 1.000 3.4 0.3622
crtIa 0.037597 0.02635 1.4268 3.143 3.6 0.5479
crtIb 0.034157 0.02922 1.1690 3.286 3.4 0.7166
cruA  0.017778 0.03628 0.4900 2.000 1.6 0.8044
cruP  0.017778 0.03628 0.4900 2.000 1.6 0.8923
crtE  0.013529 0.02221 0.6092 2.143 1.8 0.9591
crtBa 0.008287 0.01339 0.6187 1.143 1.2 1.0000
crtL  0.000000 0.00000    NaN 0.000 0.0 1.0000

Contrast: tropical_temperate 

average      sd  ratio   ava avb cumsum
crtIb 0.047066 0.03716 1.2664 3.286 2.0 0.2655
crtIa 0.041571 0.03497 1.1888 3.143 2.0 0.5001
crtE  0.025084 0.02470 1.0156 2.143 1.5 0.6416
cruA  0.019589 0.02047 0.9568 2.000 1.5 0.7521
cruP  0.019589 0.02047 0.9568 2.000 1.5 0.8626
crtL  0.019589 0.02047 0.9568 0.000 0.5 0.9731
crtBa 0.004762 0.01210 0.3934 1.143 1.0 1.0000
crtBb 0.000000 0.00000    NaN 1.000 1.0 1.0000

Contrast: soil_rock 

average      sd  ratio ava avb cumsum
crtIa 0.083588 0.03675 2.2745 3.8   7 0.4265
crtIb 0.053059 0.03803 1.3953 4.0   6 0.6973
crtBb 0.044924 0.03115 1.4423 2.8   1 0.9265
cruP  0.004878 0.01091 0.4472 2.2   2 0.9514
crtE  0.004762 0.01065 0.4472 2.2   2 0.9757
cruA  0.004762 0.01065 0.4472 2.2   2 1.0000
crtBa 0.000000 0.00000    NaN 1.0   1 1.0000
crtL  0.000000 0.00000    NaN 0.0   0 1.0000

Contrast: soil_hot spring 

average      sd  ratio ava avb cumsum
crtBb 0.04443 0.03025 1.4685 2.8 3.4 0.2634
crtIb 0.03638 0.03327 1.0933 4.0 3.4 0.4791
crtIa 0.02768 0.03449 0.8025 3.8 3.6 0.6432
cruP  0.02103 0.03436 0.6121 2.2 1.6 0.7679
cruA  0.02088 0.03411 0.6121 2.2 1.6 0.8916
crtE  0.01315 0.01977 0.6650 2.2 1.8 0.9696
crtBa 0.00513 0.01051 0.4881 1.0 1.2 1.0000
crtL  0.00000 0.00000    NaN 0.0 0.0 1.0000

Contrast: soil_temperate 

average      sd ratio ava avb cumsum
crtIb 0.06588 0.03864 1.705 4.0 2.0 0.2409
crtBb 0.05999 0.03857 1.555 2.8 1.0 0.4603
crtIa 0.05943 0.03514 1.691 3.8 2.0 0.6776
cruP  0.02375 0.02255 1.053 2.2 1.5 0.7645
crtE  0.02355 0.02215 1.063 2.2 1.5 0.8506
cruA  0.02355 0.02215 1.063 2.2 1.5 0.9367
crtL  0.01730 0.01844 0.938 0.0 0.5 1.0000
crtBa 0.00000 0.00000   NaN 1.0 1.0 1.0000

Contrast: rock_hot spring 

average      sd  ratio ava avb cumsum
crtIa 0.095003 0.04343 2.1874   7 3.6 0.3536
crtIb 0.073339 0.03868 1.8963   6 3.4 0.6265
crtBb 0.061118 0.03339 1.8302   1 3.4 0.8539
cruA  0.013793 0.03084 0.4472   2 1.6 0.9053
cruP  0.013793 0.03084 0.4472   2 1.6 0.9566
crtE  0.006897 0.01542 0.4472   2 1.8 0.9823
crtBa 0.004762 0.01065 0.4472   1 1.2 1.0000
crtL  0.000000 0.00000    NaN   0 0.0 1.0000

Contrast: rock_temperate 

average     sd  ratio ava avb cumsum
crtIa 0.15625 0.0000    Inf   7 2.0 0.4545
crtIb 0.12500 0.0000    Inf   6 2.0 0.8182
crtE  0.01562 0.0221 0.7071   2 1.5 0.8636
cruA  0.01562 0.0221 0.7071   2 1.5 0.9091
cruP  0.01562 0.0221 0.7071   2 1.5 0.9545
crtL  0.01562 0.0221 0.7071   0 0.5 1.0000
crtBa 0.00000 0.0000    NaN   1 1.0 1.0000
crtBb 0.00000 0.0000    NaN   1 1.0 1.0000

Contrast: hot spring_temperate 

average      sd  ratio ava avb cumsum
crtBb 0.08253 0.03976 2.0758 3.4 1.0 0.2891
crtIa 0.05401 0.02877 1.8771 3.6 2.0 0.4782
crtIb 0.04660 0.02723 1.7114 3.4 2.0 0.6415
cruA  0.02929 0.03308 0.8855 1.6 1.5 0.7441
cruP  0.02929 0.03308 0.8855 1.6 1.5 0.8467
crtE  0.01877 0.02062 0.9101 1.8 1.5 0.9124
crtL  0.01877 0.02062 0.9101 0.0 0.5 0.9781
crtBa 0.00625 0.01318 0.4743 1.2 1.0 1.0000
Permutation: free
Number of permutations: 0
```