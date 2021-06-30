---
title: "Statistical Analysis of Accessory Carotenoid Biosynthesis Enzymes"
author: "Dionne Martin"
date: "6/28/2021"
output: pdf_document
---
  
#################### PREFACE ##########################
#Code for PCA, NMDS, and PERMANOVA/ simper test 
#of raw BLASTp hits of accessory carotenoid 
#biosynthesis enzymes: crtRa, crtRb, crtG, crtX, 
#crtW, crtOa, crtN,crtP, crtQa, crtQb,
#cruF, and crtM. Both crtOb and crtZ have lots of 0's in
#their columns and can not be used when performing PCA. 
#These enzymes are create accessory carotenoid pigments 
#that could be important for photoprotection. 
#######################################################

```{pca}
#load packages and blastp raw hits data 
library(devtools)
library(ggbiplot)
library(readr)
accessory <- read_csv("/bigdata/stajichlab/shared/projects/Cyanobacteria/Terresterial_Cyano/Cyano Carotenoids BLASTp Results - accesory dataset.csv")
#organize data to just raw hits and strain names 
#row.names(accessory) = accessory$Strain
#accessory_matrix <- data.matrix(accessory[c(3:16)])
accessory_matrix = subset(accessory, select=c(crtRa, crtRb, crtG, crtX, crtW, crtOa, crtN,crtP, crtQa, crtQb, cruF, crtM),row.names=1)
pca_accessory= prcomp(accessory_matrix, center= TRUE, scale.= TRUE)

#analyze pca statistics to find proportion of variatiance explained by each gene 
summary(pca_accessory)
str(pca_accessory)
#crtRa-0.488 crtRb-0.123 crtG-0.104 crtX-0.095
#crtW-0.078 crtOa-0.054 crtN-0.05 crtP-0.04 
#crtQa-0.03 crtQb-0.014 cruF- 0.009 crtM-0.004 

#make groups based on habitat
habitat= accessory$Habitat

#try to make group colours, this should be added to plot code 
scale_colour_manual(name= "Habitat", values= c("brown", "dark orange", "blue", "pink", "red", "dark blue", "forest green", "dark yellow", "red", "black"))

#make the PCA plot with strains as points 
#png("PCA Plot of Carotenoid BLAST Hit .png", width=1000, height=1000)
ggbiplot(pca_accessory,ellipse=FALSE,obs.scale = 3,var.scale = 1, groups=hits.habitat) +
  scale_colour_manual(name="Habitat", values= c("saddlebrown", "darkorange1", "cornflowerblue", "deeppink", "red", "dark blue", "forest green", "gold", "darkmagenta", "black"))+
  ggtitle("PCA of Raw BLASTp Hits of Accessory Enzymes")+
  theme_minimal()+
  theme(legend.position = "bottom")

#make the PCA plot with strain names (looks too messy with all names...)
#png("PCA Plot of Carotenoid BLAST Hit .png", width=1000, height=1000)
ggbiplot(pca_accessory,labels= rownames(small),ellipse=FALSE,obs.scale = 3,var.scale = 1, groups=hits.habitat) +
  scale_colour_manual(name="Habitat", values= c("saddlebrown", "darkorange1", "cornflowerblue", "deeppink", "red", "dark blue", "forest green", "gold", "darkmagenta", "black"))+
  ggtitle("PCA of Raw BLASTp Hits of Accessory Enzymes")+
  theme_minimal()+
  theme(legend.position = "bottom")

```
```{nmds}
#load packages 
library(vegan)
library(ggplot2)
library(readr)

#prep data
accessory <- read_csv("/bigdata/stajichlab/shared/projects/Cyanobacteria/Terresterial_Cyano/Cyano Carotenoids BLASTp Results - accesory dataset.csv")
accessory_matrix = subset(accessory, select=c(crtRa, crtRb, crtG, crtX, crtW, crtOa, crtN,crtP, crtQa, crtQb, cruF, crtM))
nmds_accessory_habitat= as.data.frame(small[,1:2])
#Do you need to transform your data? range should be 0-10 if not, use root function
range((accessory_matrix)^0.5) # 0-11 should be 0-10 but lets see if the data looks fine 

#make a distance matrix
accessory_matrix_ds=vegdist(accessory_matrix, method = 'euclidian')

#Run metaMDS on distance matrix
accessory_nms=metaMDS(accessory_matrix_ds, distance="euclidian",k=2,trymax=500)
accessory_nms
#stress=0.101

#assign colors to habitat
cols <- nmds_accessory_habitat$Strain
cols[] <- c("saddlebrown", "darkorange1", "cornflowerblue", "deeppink", "red", "dark blue", "forest green", "gold", "darkmagenta", "black")
#Warning message:
#number of items to replace is not a multiple of replacement length

cols[nmds_accessory_habitat$Habitat == c("freshwater")] <- "cornflowerblue"
cols[nmds_accessory_habitat$Habitat == c("saltwater")] <- "dark blue"
cols[nmds_accessory_habitat$Habitat == c("hot spring")] <- "deeppink"
cols[nmds_accessory_habitat$Habitat == c("soil")] <- "forest green"
cols[nmds_accessory_habitat$Habitat == c("cave")] <- "saddlebrown"
cols[nmds_accessory_habitat$Habitat == c("desert")] <- "darkorange1"
cols[nmds_accessory_habitat$Habitat == c("tropical")] <- "black"
cols[nmds_accessory_habitat$Habitat == c("temperate")] <- "darkmagenta"
cols[nmds_accessory_habitat$Habitat == c("symbiont")] <- "gold"
cols[nmds_accessory_habitat$Habitat == c("rock")] <- "red"
pairs(character(nrow(nmds_accessory_habitat$Strain, col=cols)
                
#habitatvec = factor(small$Habitat, levels = c("freshwater", "saltwater", "hot spring", "soil", "cave", "desert", "tropical", "temperate", "symbiont", "rock"))
#colors <- c("cornflowerblue", "dark blue", "deeppink", "forest green", "saddlebrown", "darkorange1", "black", "darkmagenta", "gold", "red")
habitatinfo= c(rep(accessory$Habitat))
                
#NMDS Plot With Species Names 
ordiplot(accessory_nms, type="n")
orditorp(accessory_nms,display="sites",labels=accessory$Strain, col=cols ,air=0.00025)
par(xpd=TRUE)
#legend(0.6,0.25,legend =unique(nmds_small_habitat$Strain), fill=cols, border= "black")
#change coordinates of the legend for each specific graph
#legend(-0.4,-0.2, "stress = 0.13") 
#change coordinates and stress value for each specific graph
title(main="NMDS of Raw BLASTp")
scores <-
  scores(accessory_nms, display = "sites", "species")
cent <-
  aggregate(scores ~ habitatinfo, data = nmds_accessory_habitat, FUN = "mean")
names(cent) [-1] <- colnames(scores)
      points(cent [,-1],
      pch = c( 8 , 8 , 8, 8),
      col = unique(colors),
      bg = c("black"),
      lwd = 3.0,
      cex = 2.0
)
                
#NMDS Plot with Points
ordiplot(accessory_nms, choices= c(1,2))
points(accessory_nms,display= "sites",pch=21, bg=cols) 
par(xpd=TRUE)
#legend(0.6,0.25,legend =habitat, fill=cols, border= "black")
#legend(-0.4,-0.2, "stress = 0.101") 
title(main="NMDS of Raw BLASTp")
scores <-
  scores(accessory_nms, display = "sites", "species")
cent <-
  aggregate(scores ~ habitatinfo, data = nmds_accessory_habitat, FUN = "mean")
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
rownames(accessory_matrix) = accessory$Strain
#habitat= small$Habitat
str(accessory_matrix)
nmds_accessory_habitat= as.data.frame(accessory[,1:2])
str(nmds_accessory_habitat)
#habitatinfo= c(rep(small$Habitat))

#make a distance matrix
accessory_matrix_ds=vegdist(accessory_matrix^0.5, method = 'euclidian')

#Run metaMDS on distance matrix
accessory_nms=metaMDS(accessory_matrix_ds, distance="euclidian",k=2,trymax=500)
#stress=0.161

#Run PERMANOVA
pmv_accessory = adonis(accessory_matrix ~ habitat, data = nmds_accessory_habitat, permutations = 999, method = 'euclidean')
pmv_accessory
#PERMANOVA Results: Based off of the carotenoid genes being compared
#cyanobacterial strains differ statistically significant 
#between their habitats (p-value:0.001).
#Habitats explian 21.76% of the variance in the copy numbers of critical 
#carotenoid biosynthesis enzymes between cyanobacterial species. 


#Visualise the permanova
densityplot(permustats(pmv_accessory))

#Distance Based Dispersion test 
bd_accessory = betadisper(accessory_matrix_ds, accessory$Habitat)
boxplot(bd_accessory)
anova(bd_accessory) 
#F-Test: pvalue 0.1009
permutest(bd_accessory)
#permutation test: pvalue 0.111

#both the F-Test and permutation test show that 
#the avergae distance to the median for the 
#habitat data is not statistically significant. 

#simper test: determine which genes are responsible for the difference between 2 compared habiatat groups 
simp_accessory= simper(accessory_matrix, group= accessory$Habitat)
summary(simp_accessory)
#Results of simper test

Contrast: saltwater_symbiont 

average       sd  ratio    ava    avb cumsum
crtN  0.055210 0.039622 1.3934 3.6667 5.3333 0.1880
crtOa 0.042306 0.030137 1.4038 2.2917 4.0000 0.3321
crtP  0.040316 0.030942 1.3029 2.9583 4.0000 0.4694
crtQb 0.032894 0.024414 1.3473 2.3333 3.6667 0.5814
crtX  0.032585 0.031968 1.0193 0.2917 1.6667 0.6924
crtQa 0.026979 0.019676 1.3712 2.2917 3.3333 0.7843
cruF  0.019485 0.011764 1.6564 0.2083 1.0000 0.8507
crtW  0.018929 0.015799 1.1981 0.2500 1.0000 0.9151
crtG  0.013969 0.010933 1.2777 0.1667 0.6667 0.9627
crtRb 0.007101 0.009628 0.7376 1.0833 1.3333 0.9869
crtM  0.002115 0.007595 0.2784 1.0833 1.0000 0.9941
crtRa 0.001734 0.008568 0.2024 1.0833 1.0000 1.0000

Contrast: saltwater_freshwater 

average      sd  ratio    ava    avb cumsum
crtN  0.043784 0.03799 1.1524 3.6667 4.3214 0.1666
crtP  0.043478 0.03353 1.2967 2.9583 3.7500 0.3319
crtOa 0.035681 0.02808 1.2706 2.2917 3.0000 0.4677
crtX  0.025094 0.01938 1.2945 0.2917 1.0714 0.5631
crtG  0.022438 0.01809 1.2404 0.1667 0.8571 0.6485
crtQa 0.021374 0.01932 1.1063 2.2917 2.8214 0.7298
crtQb 0.021248 0.01924 1.1042 2.3333 2.8214 0.8106
cruF  0.017591 0.01244 1.4142 0.2083 0.8214 0.8775
crtW  0.011796 0.01445 0.8161 0.2500 0.3929 0.9224
crtRb 0.008613 0.01478 0.5828 1.0833 1.3214 0.9552
crtM  0.006059 0.01147 0.5284 1.0833 1.1786 0.9782
crtRa 0.005723 0.01179 0.4854 1.0833 1.1786 1.0000

Contrast: saltwater_cave 

average      sd  ratio    ava avb cumsum
crtX  0.049939 0.01815 2.7516 0.2917 2.3 0.1717
crtN  0.039800 0.03861 1.0307 3.6667 4.8 0.3086
crtP  0.036361 0.02775 1.3104 2.9583 3.8 0.4336
crtOa 0.030022 0.02325 1.2915 2.2917 3.2 0.5369
crtQa 0.026458 0.02155 1.2278 2.2917 2.9 0.6279
crtQb 0.024917 0.01746 1.4275 2.3333 2.6 0.7136
crtW  0.023015 0.02094 1.0992 0.2500 1.1 0.7927
cruF  0.018344 0.01156 1.5868 0.2083 0.9 0.8558
crtRb 0.012827 0.01313 0.9767 1.0833 1.1 0.8999
crtG  0.010995 0.01670 0.6583 0.1667 0.4 0.9377
crtRa 0.009479 0.01460 0.6491 1.0833 0.9 0.9703
crtM  0.008637 0.01250 0.6908 1.0833 1.3 1.0000

Contrast: saltwater_desert 

average      sd  ratio    ava    avb cumsum
crtN  0.055400 0.04030 1.3747 3.6667 5.4074 0.1770
crtX  0.054684 0.03965 1.3793 0.2917 2.7037 0.3518
crtP  0.044070 0.03188 1.3822 2.9583 4.2593 0.4926
crtOa 0.042484 0.03099 1.3710 2.2917 3.9259 0.6284
crtQa 0.024900 0.02210 1.1265 2.2917 3.1111 0.7080
crtQb 0.022610 0.02023 1.1179 2.3333 3.0370 0.7802
crtW  0.018134 0.01357 1.3366 0.2500 1.0000 0.8382
cruF  0.017391 0.01117 1.5568 0.2083 0.9259 0.8938
crtG  0.009482 0.01431 0.6624 0.1667 0.3333 0.9241
crtRa 0.008709 0.01635 0.5326 1.0833 1.2593 0.9519
crtRb 0.008073 0.01107 0.7289 1.0833 1.2593 0.9777
crtM  0.006978 0.01112 0.6273 1.0833 1.1852 1.0000

Contrast: saltwater_tropical 

average       sd  ratio    ava   avb cumsum
crtX  0.049920 0.036170 1.3801 0.2917 2.714 0.1464
crtN  0.048437 0.033332 1.4532 3.6667 5.571 0.2885
crtP  0.041310 0.029018 1.4236 2.9583 4.571 0.4097
crtQb 0.039908 0.032756 1.2183 2.3333 4.286 0.5267
crtOa 0.035462 0.024796 1.4302 2.2917 3.857 0.6307
crtQa 0.034632 0.030336 1.1416 2.2917 4.000 0.7323
crtG  0.021980 0.016635 1.3213 0.1667 1.143 0.7968
crtW  0.019008 0.013913 1.3662 0.2500 1.143 0.8525
crtRa 0.018521 0.017331 1.0686 1.0833 2.000 0.9069
cruF  0.017081 0.009513 1.7956 0.2083 1.000 0.9570
crtRb 0.007684 0.011220 0.6848 1.0833 1.429 0.9795
crtM  0.006986 0.010408 0.6712 1.0833 1.286 1.0000

Contrast: saltwater_soil 

average       sd  ratio    ava avb cumsum
crtOa 0.067590 0.037188 1.8175 2.2917 5.6 0.2023
crtN  0.065197 0.043766 1.4897 3.6667 6.8 0.3975
crtP  0.052459 0.037698 1.3916 2.9583 5.4 0.5545
crtX  0.035742 0.017397 2.0545 0.2917 2.0 0.6615
crtQa 0.030337 0.018658 1.6260 2.2917 3.6 0.7523
crtQb 0.026516 0.018250 1.4529 2.3333 3.4 0.8316
cruF  0.017135 0.009544 1.7953 0.2083 1.0 0.8829
crtRa 0.015189 0.013028 1.1659 1.0833 1.8 0.9284
crtW  0.011804 0.012605 0.9365 0.2500 0.6 0.9637
crtG  0.006207 0.011286 0.5500 0.1667 0.2 0.9823
crtRb 0.004081 0.007165 0.5696 1.0833 1.2 0.9945
crtM  0.001837 0.006359 0.2889 1.0833 1.0 1.0000

Contrast: saltwater_rock 

average       sd  ratio    ava avb cumsum
crtQa 0.043147 0.017614 2.4496 2.2917   4 0.1565
crtW  0.040665 0.010900 3.7306 0.2500   2 0.3041
crtP  0.027807 0.018582 1.4965 2.9583   3 0.4050
crtN  0.025794 0.031514 0.8185 3.6667   4 0.4986
crtOa 0.025253 0.019113 1.3212 2.2917   3 0.5902
crtRa 0.023308 0.003058 7.6216 1.0833   2 0.6748
crtX  0.022030 0.007347 2.9986 0.2917   1 0.7547
crtQb 0.021825 0.012569 1.7364 2.3333   3 0.8339
crtRb 0.021689 0.007247 2.9927 1.0833   2 0.9126
cruF  0.019312 0.010305 1.8740 0.2083   1 0.9826
crtG  0.002709 0.010022 0.2703 0.1667   0 0.9925
crtM  0.002080 0.007220 0.2880 1.0833   1 1.0000

Contrast: saltwater_hot spring 

average       sd  ratio    ava avb cumsum
crtN  0.086033 0.055500 1.5501 3.6667 4.2 0.2221
crtP  0.077555 0.042073 1.8433 2.9583 4.0 0.4222
crtOa 0.068025 0.030534 2.2278 2.2917 3.6 0.5978
crtQb 0.030610 0.027753 1.1029 2.3333 3.4 0.6768
crtX  0.028445 0.027995 1.0161 0.2917 1.4 0.7502
crtQa 0.026581 0.024422 1.0884 2.2917 3.2 0.8188
cruF  0.023202 0.016362 1.4180 0.2083 1.0 0.8787
crtW  0.016083 0.017940 0.8965 0.2500 0.6 0.9202
crtM  0.011570 0.022808 0.5073 1.0833 0.8 0.9501
crtRb 0.010780 0.019708 0.5470 1.0833 0.8 0.9779
crtG  0.006587 0.012589 0.5232 0.1667 0.2 0.9949
crtRa 0.001975 0.009999 0.1975 1.0833 1.0 1.0000

Contrast: saltwater_temperate 

average       sd  ratio    ava avb cumsum
crtP  0.047095 0.036293 1.2976 2.9583 4.5 0.1697
crtN  0.046281 0.042537 1.0880 3.6667 5.0 0.3366
crtW  0.030982 0.017132 1.8084 0.2500 1.5 0.4482
crtOa 0.026928 0.020405 1.3197 2.2917 3.0 0.5453
crtRa 0.024781 0.003417 7.2523 1.0833 2.0 0.6346
cruF  0.020580 0.010895 1.8889 0.2083 1.0 0.7088
crtQb 0.019475 0.017595 1.1069 2.3333 1.5 0.7790
crtQa 0.018559 0.017450 1.0636 2.2917 1.5 0.8459
crtX  0.014700 0.014350 1.0244 0.2917 0.5 0.8989
crtG  0.013581 0.013751 0.9877 0.1667 0.5 0.9478
crtRb 0.012252 0.012596 0.9727 1.0833 1.5 0.9920
crtM  0.002226 0.007669 0.2902 1.0833 1.0 1.0000

Contrast: symbiont_freshwater 

average       sd  ratio    ava    avb cumsum
crtN  0.043130 0.032402 1.3311 5.3333 4.3214 0.1784
crtOa 0.036474 0.025426 1.4345 4.0000 3.0000 0.3293
crtP  0.034394 0.025629 1.3420 4.0000 3.7500 0.4716
crtX  0.029570 0.023647 1.2505 1.6667 1.0714 0.5939
crtQb 0.026590 0.018376 1.4470 3.6667 2.8214 0.7039
crtQa 0.020842 0.015949 1.3067 3.3333 2.8214 0.7901
crtW  0.017065 0.014257 1.1970 1.0000 0.3929 0.8607
crtG  0.013318 0.014481 0.9196 0.6667 0.8571 0.9158
crtRb 0.009160 0.012140 0.7545 1.3333 1.3214 0.9537
cruF  0.004051 0.009081 0.4461 1.0000 0.8214 0.9705
crtM  0.003709 0.008356 0.4439 1.0000 1.1786 0.9858
crtRa 0.003423 0.007596 0.4507 1.0000 1.1786 1.0000

Contrast: symbiont_cave 

average       sd  ratio    ava avb cumsum
crtN  0.038280 0.024540 1.5599 5.3333 4.8 0.1568
crtX  0.035844 0.017841 2.0091 1.6667 2.3 0.3036
crtOa 0.030178 0.012362 2.4412 4.0000 3.2 0.4272
crtQb 0.030096 0.019954 1.5083 3.6667 2.6 0.5505
crtP  0.028459 0.018259 1.5586 4.0000 3.8 0.6671
crtQa 0.023815 0.018613 1.2795 3.3333 2.9 0.7646
crtW  0.018304 0.017305 1.0577 1.0000 1.1 0.8396
crtRb 0.012294 0.012205 1.0073 1.3333 1.1 0.8899
crtG  0.012250 0.011192 1.0945 0.6667 0.4 0.9401
crtM  0.006371 0.010268 0.6204 1.0000 1.3 0.9662
crtRa 0.006334 0.010215 0.6200 1.0000 0.9 0.9922
cruF  0.001915 0.005964 0.3211 1.0000 0.9 1.0000

Contrast: symbiont_desert 

average       sd  ratio    ava    avb cumsum
crtN  0.043237 0.032551 1.3283 5.3333 5.4074 0.1817
crtX  0.040216 0.033452 1.2022 1.6667 2.7037 0.3507
crtP  0.031844 0.026081 1.2210 4.0000 4.2593 0.4845
crtOa 0.029966 0.025755 1.1635 4.0000 3.9259 0.6104
crtQb 0.025066 0.017315 1.4476 3.6667 3.0370 0.7157
crtQa 0.020191 0.017008 1.1871 3.3333 3.1111 0.8005
crtW  0.014745 0.012645 1.1661 1.0000 1.0000 0.8625
crtG  0.010567 0.009832 1.0747 0.6667 0.3333 0.9069
crtRb 0.008322 0.009642 0.8631 1.3333 1.2593 0.9419
crtRa 0.006147 0.012865 0.4778 1.0000 1.2593 0.9677
crtM  0.004799 0.008406 0.5709 1.0000 1.1852 0.9879
cruF  0.002891 0.007287 0.3968 1.0000 0.9259 1.0000

Contrast: symbiont_tropical 

average       sd  ratio    ava   avb cumsum
crtX  0.037175 0.030774 1.2080 1.6667 2.714 0.1544
crtN  0.035387 0.027271 1.2976 5.3333 5.571 0.3013
crtQb 0.031077 0.026634 1.1668 3.6667 4.286 0.4304
crtP  0.027715 0.023646 1.1720 4.0000 4.571 0.5455
crtOa 0.027019 0.018819 1.4357 4.0000 3.857 0.6577
crtQa 0.025058 0.026029 0.9627 3.3333 4.000 0.7617
crtRa 0.015420 0.015418 1.0001 1.0000 2.000 0.8257
crtG  0.014801 0.012981 1.1402 0.6667 1.143 0.8872
crtW  0.013504 0.012723 1.0614 1.0000 1.143 0.9433
crtRb 0.008463 0.009696 0.8729 1.3333 1.429 0.9784
crtM  0.005194 0.008650 0.6005 1.0000 1.286 1.0000
cruF  0.000000 0.000000    NaN 1.0000 1.000 1.0000

Contrast: symbiont_soil 

average       sd  ratio    ava avb cumsum
crtN  0.045012 0.033595 1.3398 5.3333 6.8 0.1967
crtOa 0.038066 0.032975 1.1544 4.0000 5.6 0.3631
crtP  0.037385 0.028958 1.2910 4.0000 5.4 0.5265
crtX  0.028469 0.019785 1.4389 1.6667 2.0 0.6509
crtQb 0.021143 0.017905 1.1809 3.6667 3.4 0.7434
crtQa 0.015551 0.019990 0.7779 3.3333 3.6 0.8113
crtW  0.015078 0.013158 1.1459 1.0000 0.6 0.8772
crtRa 0.012513 0.011408 1.0968 1.0000 1.8 0.9319
crtG  0.009504 0.008374 1.1349 0.6667 0.2 0.9735
crtRb 0.006074 0.007773 0.7815 1.3333 1.2 1.0000
cruF  0.000000 0.000000    NaN 1.0000 1.0 1.0000
crtM  0.000000 0.000000    NaN 1.0000 1.0 1.0000

Contrast: symbiont_rock 

average       sd  ratio    ava avb cumsum
crtN  0.03582 0.027014 1.3258 5.3333   4 0.1577
crtOa 0.03017 0.005013 6.0184 4.0000   3 0.2905
crtP  0.03017 0.005013 6.0184 4.0000   3 0.4233
crtQb 0.02461 0.008062 3.0530 3.6667   3 0.5317
crtX  0.02452 0.024590 0.9973 1.6667   1 0.6397
crtW  0.02172 0.024823 0.8752 1.0000   2 0.7353
crtRa 0.01915 0.004540 4.2180 1.0000   2 0.8196
crtQa 0.01626 0.028163 0.5774 3.3333   4 0.8912
crtRb 0.01369 0.012465 1.0979 1.3333   2 0.9515
crtG  0.01102 0.009545 1.1546 0.6667   0 1.0000
cruF  0.00000 0.000000    NaN 1.0000   1 1.0000
crtM  0.00000 0.000000    NaN 1.0000   1 1.0000

Contrast: symbiont_hot spring 

average      sd  ratio    ava avb cumsum
crtN  0.077018 0.05946 1.2953 5.3333 4.2 0.2331
crtP  0.061744 0.04249 1.4531 4.0000 4.0 0.4200
crtOa 0.055579 0.04420 1.2576 4.0000 3.6 0.5883
crtX  0.033374 0.03178 1.0503 1.6667 1.4 0.6893
crtQb 0.029778 0.02694 1.1055 3.6667 3.4 0.7794
crtQa 0.024285 0.02138 1.1358 3.3333 3.2 0.8530
crtW  0.018314 0.01648 1.1115 1.0000 0.6 0.9084
crtRb 0.012400 0.01700 0.7295 1.3333 0.8 0.9459
crtG  0.011393 0.01027 1.1096 0.6667 0.2 0.9804
crtM  0.006467 0.01429 0.4525 1.0000 0.8 1.0000
crtRa 0.000000 0.00000    NaN 1.0000 1.0 1.0000
cruF  0.000000 0.00000    NaN 1.0000 1.0 1.0000

Contrast: symbiont_temperate 

average       sd  ratio    ava avb cumsum
crtN  0.04035 0.026609 1.5163 5.3333 5.0 0.1569
crtQb 0.03866 0.023946 1.6145 3.6667 1.5 0.3072
crtQa 0.03296 0.019073 1.7283 3.3333 1.5 0.4353
crtOa 0.03165 0.004416 7.1682 4.0000 3.0 0.5583
crtP  0.02739 0.030955 0.8848 4.0000 4.5 0.6648
crtX  0.02718 0.027674 0.9821 1.6667 0.5 0.7705
crtRa 0.02016 0.004519 4.4597 1.0000 2.0 0.8488
crtW  0.01879 0.019510 0.9632 1.0000 1.5 0.9219
crtG  0.01007 0.011436 0.8806 0.6667 0.5 0.9610
crtRb 0.01002 0.011397 0.8792 1.3333 1.5 1.0000
cruF  0.00000 0.000000    NaN 1.0000 1.0 1.0000
crtM  0.00000 0.000000    NaN 1.0000 1.0 1.0000

Contrast: freshwater_cave 

average       sd  ratio    ava avb cumsum
crtN  0.029010 0.028253 1.0268 4.3214 4.8 0.1294
crtX  0.028329 0.019348 1.4642 1.0714 2.3 0.2557
crtOa 0.024652 0.022188 1.1110 3.0000 3.2 0.3656
crtP  0.024528 0.024026 1.0209 3.7500 3.8 0.4750
crtQa 0.022015 0.018017 1.2219 2.8214 2.9 0.5732
crtQb 0.021347 0.017486 1.2208 2.8214 2.6 0.6684
crtW  0.019993 0.017842 1.1206 0.3929 1.1 0.7575
crtG  0.017125 0.014480 1.1827 0.8571 0.4 0.8339
crtRb 0.014091 0.015047 0.9364 1.3214 1.1 0.8968
crtRa 0.009627 0.012585 0.7650 1.1786 0.9 0.9397
crtM  0.008146 0.010819 0.7529 1.1786 1.3 0.9760
cruF  0.005379 0.009656 0.5571 0.8214 0.9 1.0000

Contrast: freshwater_desert 

average       sd  ratio    ava    avb cumsum
crtN  0.042344 0.032555 1.3007 4.3214 5.4074 0.1715
crtX  0.037669 0.033620 1.1204 1.0714 2.7037 0.3240
crtOa 0.034428 0.027397 1.2566 3.0000 3.9259 0.4634
crtP  0.033073 0.026910 1.2290 3.7500 4.2593 0.5974
crtQa 0.018949 0.016662 1.1373 2.8214 3.1111 0.6741
crtQb 0.017871 0.015889 1.1247 2.8214 3.0370 0.7465
crtW  0.015965 0.011868 1.3453 0.3929 1.0000 0.8111
crtG  0.015310 0.013351 1.1468 0.8571 0.3333 0.8731
crtRb 0.009721 0.012293 0.7908 1.3214 1.2593 0.9125
crtRa 0.008431 0.013565 0.6215 1.1786 1.2593 0.9466
crtM  0.007179 0.010243 0.7009 1.1786 1.1852 0.9757
cruF  0.005998 0.009886 0.6067 0.8214 0.9259 1.0000

Contrast: freshwater_tropical 

average       sd  ratio    ava   avb cumsum
crtN  0.034667 0.026897 1.2889 4.3214 5.571 0.1368
crtX  0.034214 0.030895 1.1074 1.0714 2.714 0.2718
crtQb 0.032307 0.026365 1.2254 2.8214 4.286 0.3992
crtOa 0.028608 0.022430 1.2755 3.0000 3.857 0.5121
crtP  0.028245 0.024016 1.1761 3.7500 4.571 0.6235
crtQa 0.027056 0.025673 1.0539 2.8214 4.000 0.7303
crtW  0.016386 0.011981 1.3677 0.3929 1.143 0.7949
crtRa 0.015997 0.014887 1.0746 1.1786 2.000 0.8580
crtG  0.015928 0.012955 1.2295 0.8571 1.143 0.9209
crtRb 0.009610 0.012048 0.7976 1.3214 1.429 0.9588
crtM  0.006854 0.009312 0.7360 1.1786 1.286 0.9858
cruF  0.003593 0.007879 0.4561 0.8214 1.000 1.0000

Contrast: freshwater_soil 

average       sd  ratio    ava avb cumsum
crtOa 0.050382 0.035059 1.4371 3.0000 5.6 0.2048
crtN  0.047262 0.037837 1.2491 4.3214 6.8 0.3970
crtP  0.034648 0.033765 1.0261 3.7500 5.4 0.5378
crtX  0.021736 0.017261 1.2593 1.0714 2.0 0.6262
crtQa 0.020515 0.015237 1.3464 2.8214 3.6 0.7096
crtQb 0.018386 0.014780 1.2440 2.8214 3.4 0.7844
crtG  0.014478 0.012717 1.1385 0.8571 0.2 0.8432
crtRa 0.012827 0.011206 1.1446 1.1786 1.8 0.8954
crtW  0.011707 0.011963 0.9786 0.3929 0.6 0.9430
crtRb 0.007092 0.010925 0.6492 1.3214 1.2 0.9718
cruF  0.003604 0.007907 0.4558 0.8214 1.0 0.9865
crtM  0.003331 0.007346 0.4534 1.1786 1.0 1.0000

Contrast: freshwater_rock 

average       sd  ratio    ava avb cumsum
crtW  0.033234 0.013135 2.5301 0.3929   2 0.1584
crtQa 0.025261 0.017783 1.4205 2.8214   4 0.2788
crtP  0.025101 0.018731 1.3401 3.7500   3 0.3984
crtN  0.021128 0.022785 0.9273 4.3214   4 0.4991
crtOa 0.020665 0.020851 0.9911 3.0000   3 0.5976
crtG  0.017678 0.014393 1.2282 0.8571   0 0.6819
crtRa 0.017019 0.008381 2.0306 1.1786   2 0.7630
crtRb 0.016454 0.009188 1.7909 1.3214   2 0.8414
crtX  0.012783 0.011802 1.0831 1.0714   1 0.9024
crtQb 0.012742 0.010659 1.1954 2.8214   3 0.9631
cruF  0.004035 0.008855 0.4557 0.8214   1 0.9823
crtM  0.003706 0.008197 0.4521 1.1786   1 1.0000

Contrast: freshwater_hot spring 

average       sd  ratio    ava avb cumsum
crtN  0.076091 0.050249 1.5143 4.3214 4.2 0.2161
crtP  0.068636 0.043397 1.5816 3.7500 4.0 0.4110
crtOa 0.059851 0.035812 1.6712 3.0000 3.6 0.5810
crtX  0.028652 0.022689 1.2628 1.0714 1.4 0.6624
crtQb 0.027387 0.020965 1.3063 2.8214 3.4 0.7402
crtQa 0.023680 0.019520 1.2131 2.8214 3.2 0.8074
crtG  0.018785 0.018304 1.0263 0.8571 0.2 0.8608
crtW  0.014968 0.015565 0.9616 0.3929 0.6 0.9033
crtRb 0.014201 0.021724 0.6537 1.3214 0.8 0.9436
crtM  0.011291 0.018737 0.6026 1.1786 0.8 0.9757
cruF  0.004717 0.011033 0.4276 0.8214 1.0 0.9891
crtRa 0.003847 0.008765 0.4390 1.1786 1.0 1.0000

Contrast: freshwater_temperate 

average       sd  ratio    ava avb cumsum
crtN  0.032178 0.031065 1.0358 4.3214 5.0 0.1393
crtQa 0.027109 0.017950 1.5103 2.8214 1.5 0.2566
crtQb 0.027109 0.017950 1.5103 2.8214 1.5 0.3740
crtP  0.026935 0.029549 0.9116 3.7500 4.5 0.4905
crtW  0.025630 0.015451 1.6587 0.3929 1.5 0.6015
crtOa 0.021825 0.021902 0.9965 3.0000 3.0 0.6959
crtRa 0.017965 0.008805 2.0402 1.1786 2.0 0.7737
crtX  0.017563 0.015210 1.1547 1.0714 0.5 0.8497
crtG  0.014596 0.013744 1.0620 0.8571 0.5 0.9129
crtRb 0.011931 0.012454 0.9580 1.3214 1.5 0.9645
cruF  0.004280 0.009312 0.4596 0.8214 1.0 0.9831
crtM  0.003914 0.008591 0.4556 1.1786 1.0 1.0000

Contrast: cave_desert 

average       sd  ratio ava    avb cumsum
crtN  0.035402 0.025441 1.3915 4.8 5.4074 0.1605
crtX  0.029115 0.021696 1.3420 2.3 2.7037 0.2926
crtP  0.026764 0.019018 1.4073 3.8 4.2593 0.4139
crtOa 0.025647 0.018776 1.3659 3.2 3.9259 0.5302
crtQa 0.022445 0.018439 1.2173 2.9 3.1111 0.6320
crtQb 0.021377 0.018025 1.1859 2.6 3.0370 0.7289
crtW  0.012821 0.014413 0.8895 1.1 1.0000 0.7871
crtRb 0.012444 0.012481 0.9970 1.1 1.2593 0.8435
crtRa 0.011315 0.015421 0.7337 0.9 1.2593 0.8948
crtG  0.010388 0.012706 0.8175 0.4 0.3333 0.9419
crtM  0.008457 0.010728 0.7884 1.3 1.1852 0.9802
cruF  0.004358 0.008352 0.5218 0.9 0.9259 1.0000

Contrast: cave_tropical 

average       sd  ratio ava   avb cumsum
crtQb 0.035243 0.027394 1.2865 2.6 4.286 0.1534
crtQa 0.029419 0.025675 1.1458 2.9 4.000 0.2815
crtN  0.027890 0.016728 1.6673 4.8 5.571 0.4028
crtX  0.027185 0.018327 1.4834 2.3 2.714 0.5212
crtP  0.021587 0.015597 1.3840 3.8 4.571 0.6151
crtRa 0.019571 0.017004 1.1510 0.9 2.000 0.7003
crtOa 0.019056 0.014524 1.3120 3.2 3.857 0.7833
crtG  0.017406 0.013902 1.2520 0.4 1.143 0.8590
crtRb 0.012561 0.011750 1.0690 1.1 1.429 0.9137
crtW  0.010508 0.012725 0.8258 1.1 1.143 0.9594
crtM  0.007584 0.009253 0.8196 1.3 1.286 0.9924
cruF  0.001736 0.005304 0.3273 0.9 1.000 1.0000

Contrast: cave_soil 

average       sd  ratio ava avb cumsum
crtOa 0.040037 0.027929 1.4335 3.2 5.6 0.1755
crtN  0.037938 0.029155 1.3012 4.8 6.8 0.3419
crtP  0.028907 0.026519 1.0901 3.8 5.4 0.4686
crtQa 0.022676 0.017553 1.2918 2.9 3.6 0.5680
crtQb 0.021743 0.017971 1.2099 2.6 3.4 0.6633
crtW  0.017648 0.015362 1.1489 1.1 0.6 0.7407
crtX  0.016598 0.012323 1.3469 2.3 2.0 0.8135
crtRa 0.016558 0.013822 1.1980 0.9 1.8 0.8861
crtRb 0.010100 0.009952 1.0149 1.1 1.2 0.9304
crtG  0.008423 0.011454 0.7354 0.4 0.2 0.9673
crtM  0.005720 0.008992 0.6361 1.3 1.0 0.9924
cruF  0.001741 0.005331 0.3266 0.9 1.0 1.0000

Contrast: cave_rock 

average       sd  ratio ava avb cumsum
crtQa 0.026705 0.021099 1.2657 2.9   4 0.1456
crtX  0.025190 0.008394 3.0009 2.3   1 0.2829
crtRa 0.022186 0.013200 1.6808 0.9   2 0.4039
crtW  0.021375 0.010352 2.0648 1.1   2 0.5204
crtRb 0.018645 0.016285 1.1449 1.1   2 0.6221
crtQb 0.016430 0.017262 0.9518 2.6   3 0.7116
crtP  0.014946 0.014415 1.0368 3.8   3 0.7931
crtN  0.014389 0.020300 0.7088 4.8   4 0.8716
crtOa 0.007643 0.009893 0.7725 3.2   3 0.9132
crtG  0.007620 0.013378 0.5696 0.4   0 0.9548
crtM  0.006370 0.010291 0.6189 1.3   1 0.9895
cruF  0.001923 0.006081 0.3162 0.9   1 1.0000

Contrast: cave_hot spring 

average       sd  ratio ava avb cumsum
crtN  0.074679 0.052306 1.4277 4.8 4.2 0.2060
crtP  0.064858 0.036894 1.7579 3.8 4.0 0.3849
crtOa 0.057596 0.028392 2.0286 3.2 3.6 0.5438
crtX  0.037993 0.027779 1.3677 2.3 1.4 0.6486
crtQb 0.031289 0.021258 1.4719 2.6 3.4 0.7349
crtQa 0.027230 0.021966 1.2396 2.9 3.2 0.8100
crtW  0.021242 0.020068 1.0585 1.1 0.6 0.8686
crtRb 0.014446 0.016059 0.8995 1.1 0.8 0.9084
crtM  0.013809 0.019994 0.6906 1.3 0.8 0.9465
crtG  0.009984 0.014673 0.6804 0.4 0.2 0.9741
crtRa 0.007249 0.012104 0.5989 0.9 1.0 0.9941
cruF  0.002148 0.006802 0.3157 0.9 1.0 1.0000

Contrast: cave_temperate 

average       sd  ratio ava avb cumsum
crtX  0.036701 0.013239 2.7722 2.3 0.5 0.1690
crtQa 0.029458 0.021073 1.3979 2.9 1.5 0.3046
crtQb 0.025780 0.017341 1.4866 2.6 1.5 0.4233
crtRa 0.023377 0.013645 1.7132 0.9 2.0 0.5310
crtN  0.022427 0.020895 1.0734 4.8 5.0 0.6342
crtP  0.018992 0.015579 1.2190 3.8 4.5 0.7217
crtW  0.016389 0.014152 1.1580 1.1 1.5 0.7972
crtRb 0.014980 0.014772 1.0141 1.1 1.5 0.8661
crtG  0.012297 0.012312 0.9988 0.4 0.5 0.9227
crtOa 0.008028 0.010118 0.7934 3.2 3.0 0.9597
crtM  0.006730 0.010588 0.6356 1.3 1.0 0.9907
cruF  0.002020 0.006219 0.3249 0.9 1.0 1.0000

Contrast: desert_tropical 

average       sd  ratio    ava   avb cumsum
crtX  0.033497 0.027769 1.2063 2.7037 2.714 0.1441
crtN  0.033177 0.025319 1.3103 5.4074 5.571 0.2868
crtQb 0.030107 0.024173 1.2455 3.0370 4.286 0.4163
crtP  0.026161 0.020634 1.2679 4.2593 4.571 0.5289
crtQa 0.025730 0.023629 1.0889 3.1111 4.000 0.6395
crtOa 0.025379 0.018888 1.3436 3.9259 3.857 0.7487
crtG  0.016375 0.012670 1.2925 0.3333 1.143 0.8191
crtRa 0.016168 0.014837 1.0897 1.2593 2.000 0.8887
crtRb 0.008995 0.009900 0.9086 1.2593 1.429 0.9274
crtM  0.007256 0.009306 0.7796 1.1852 1.286 0.9586
crtW  0.007018 0.009917 0.7077 1.0000 1.143 0.9888
cruF  0.002606 0.006439 0.4048 0.9259 1.000 1.0000

Contrast: desert_soil 

average       sd  ratio    ava avb cumsum
crtN  0.040313 0.032180 1.2527 5.4074 6.8 0.1807
crtOa 0.036790 0.028467 1.2924 3.9259 5.6 0.3455
crtP  0.033217 0.026161 1.2697 4.2593 5.4 0.4944
crtX  0.027830 0.024946 1.1156 2.7037 2.0 0.6191
crtQa 0.018585 0.014757 1.2594 3.1111 3.6 0.7024
crtQb 0.017723 0.013899 1.2751 3.0370 3.4 0.7818
crtW  0.014467 0.010406 1.3903 1.0000 0.6 0.8466
crtRa 0.013450 0.012169 1.1053 1.2593 1.8 0.9069
crtG  0.007118 0.009627 0.7393 0.3333 0.2 0.9388
crtRb 0.006667 0.008358 0.7977 1.2593 1.2 0.9687
crtM  0.004374 0.007544 0.5799 1.1852 1.0 0.9883
cruF  0.002614 0.006460 0.4046 0.9259 1.0 1.0000

Contrast: desert_rock 

average       sd  ratio    ava avb cumsum
crtN  0.035615 0.022956 1.5515 5.4074   4 0.1645
crtX  0.031032 0.032078 0.9674 2.7037   1 0.3078
crtP  0.028628 0.018155 1.5769 4.2593   3 0.4400
crtOa 0.024180 0.018241 1.3256 3.9259   3 0.5517
crtQa 0.020717 0.015614 1.3268 3.1111   4 0.6473
crtW  0.018997 0.011138 1.7056 1.0000   2 0.7351
crtRa 0.016652 0.009631 1.7290 1.2593   2 0.8120
crtRb 0.013882 0.009969 1.3925 1.2593   2 0.8761
crtQb 0.012764 0.012023 1.0616 3.0370   3 0.9350
crtG  0.006357 0.010614 0.5989 0.3333   0 0.9644
crtM  0.004821 0.008347 0.5776 1.1852   1 0.9866
cruF  0.002895 0.007194 0.4024 0.9259   1 1.0000

Contrast: desert_hot spring 

average      sd  ratio    ava avb cumsum
crtN  0.075903 0.05801 1.3084 5.4074 4.2 0.2178
crtP  0.062386 0.04363 1.4299 4.2593 4.0 0.3969
crtOa 0.055854 0.04104 1.3610 3.9259 3.6 0.5571
crtX  0.044991 0.04044 1.1127 2.7037 1.4 0.6863
crtQb 0.026288 0.02031 1.2944 3.0370 3.4 0.7617
crtQa 0.024141 0.02079 1.1613 3.1111 3.2 0.8310
crtW  0.017000 0.01337 1.2710 1.0000 0.6 0.8798
crtRb 0.012451 0.01708 0.7291 1.2593 0.8 0.9155
crtM  0.010806 0.01557 0.6942 1.1852 0.8 0.9465
crtG  0.008475 0.01232 0.6877 0.3333 0.2 0.9708
crtRa 0.006884 0.01479 0.4656 1.2593 1.0 0.9906
cruF  0.003280 0.00857 0.3828 0.9259 1.0 1.0000

Contrast: desert_temperate 

average       sd  ratio    ava avb cumsum
crtX  0.041446 0.034649 1.1962 2.7037 0.5 0.1678
crtN  0.035449 0.026722 1.3266 5.4074 5.0 0.3113
crtQa 0.030491 0.020719 1.4716 3.1111 1.5 0.4347
crtQb 0.028519 0.019093 1.4937 3.0370 1.5 0.5502
crtP  0.026311 0.021483 1.2248 4.2593 4.5 0.6567
crtOa 0.025321 0.018878 1.3413 3.9259 3.0 0.7592
crtRa 0.017480 0.010029 1.7429 1.2593 2.0 0.8300
crtW  0.012935 0.012759 1.0139 1.0000 1.5 0.8824
crtG  0.010476 0.010682 0.9808 0.3333 0.5 0.9248
crtRb 0.010472 0.010648 0.9835 1.2593 1.5 0.9672
crtM  0.005059 0.008680 0.5828 1.1852 1.0 0.9877
cruF  0.003049 0.007517 0.4056 0.9259 1.0 1.0000

Contrast: tropical_soil 

average       sd  ratio   ava avb cumsum
crtOa 0.032935 0.025959 1.2688 3.857 5.6 0.1482
crtN  0.030839 0.027731 1.1121 5.571 6.8 0.2870
crtP  0.027275 0.021103 1.2925 4.571 5.4 0.4098
crtQb 0.026802 0.021723 1.2338 4.286 3.4 0.5304
crtX  0.025821 0.022656 1.1397 2.714 2.0 0.6466
crtQa 0.022088 0.022827 0.9676 4.000 3.6 0.7460
crtG  0.015810 0.012284 1.2871 1.143 0.2 0.8171
crtW  0.014603 0.010466 1.3953 1.143 0.6 0.8828
crtRa 0.014491 0.011719 1.2366 2.000 1.8 0.9481
crtRb 0.006797 0.008997 0.7555 1.429 1.2 0.9786
crtM  0.004745 0.007729 0.6140 1.286 1.0 1.0000
cruF  0.000000 0.000000    NaN 1.000 1.0 1.0000

Contrast: tropical_rock 

average       sd  ratio   ava avb cumsum
crtN  0.029902 0.014737 2.0291 5.571   4 0.1377
crtQb 0.028245 0.023120 1.2217 4.286   3 0.2677
crtX  0.027398 0.032063 0.8545 2.714   1 0.3938
crtP  0.024865 0.019383 1.2829 4.571   3 0.5083
crtQa 0.022445 0.024213 0.9270 4.000   4 0.6116
crtG  0.018491 0.014819 1.2478 1.143   0 0.6967
crtOa 0.017762 0.014329 1.2396 3.857   3 0.7785
crtRb 0.014812 0.007046 2.1023 1.429   2 0.8467
crtW  0.014448 0.006768 2.1346 1.143   2 0.9132
crtRa 0.013639 0.009671 1.4103 2.000   2 0.9760
crtM  0.005222 0.008952 0.5834 1.286   1 1.0000
cruF  0.000000 0.000000    NaN 1.000   1 1.0000

Contrast: tropical_hot spring 

average      sd  ratio   ava avb cumsum
crtN  0.06966 0.05147 1.3534 5.571 4.2 0.1933
crtP  0.05762 0.04154 1.3874 4.571 4.0 0.3532
crtOa 0.05031 0.03366 1.4944 3.857 3.6 0.4929
crtX  0.04125 0.03662 1.1265 2.714 1.4 0.6074
crtQb 0.03506 0.03099 1.1313 4.286 3.4 0.7046
crtQa 0.03062 0.02792 1.0967 4.000 3.2 0.7896
crtG  0.01902 0.01626 1.1696 1.143 0.2 0.8424
crtRa 0.01688 0.01705 0.9905 2.000 1.0 0.8893
crtW  0.01684 0.01364 1.2352 1.143 0.6 0.9360
crtRb 0.01187 0.01573 0.7547 1.429 0.8 0.9690
crtM  0.01119 0.01607 0.6962 1.286 0.8 1.0000
cruF  0.00000 0.00000    NaN 1.000 1.0 1.0000

Contrast: tropical_temperate 

average       sd  ratio   ava avb cumsum
crtQb 0.045092 0.031222 1.4442 4.286 1.5 0.1851
crtQa 0.039893 0.028995 1.3759 4.000 1.5 0.3489
crtX  0.037594 0.033179 1.1331 2.714 0.5 0.5032
crtN  0.027222 0.016960 1.6051 5.571 5.0 0.6150
crtP  0.019381 0.013779 1.4066 4.571 4.5 0.6945
crtOa 0.018505 0.014279 1.2960 3.857 3.0 0.7705
crtG  0.016244 0.012809 1.2681 1.143 0.5 0.8372
crtRa 0.014224 0.009667 1.4714 2.000 2.0 0.8955
crtRb 0.010834 0.010104 1.0722 1.429 1.5 0.9400
crtW  0.009136 0.009676 0.9442 1.143 1.5 0.9775
crtM  0.005474 0.009020 0.6069 1.286 1.0 1.0000
cruF  0.000000 0.000000    NaN 1.000 1.0 1.0000

Contrast: soil_rock 

average       sd  ratio ava avb cumsum
crtN  0.043865 0.031095 1.4107 6.8   4 0.2029
crtOa 0.042488 0.026996 1.5738 5.6   3 0.3993
crtP  0.037591 0.027804 1.3520 5.4   3 0.5732
crtW  0.026082 0.017395 1.4994 0.6   2 0.6938
crtX  0.015316 0.016734 0.9153 2.0   1 0.7646
crtRb 0.014708 0.008362 1.7590 1.2   2 0.8327
crtQb 0.014077 0.008296 1.6969 3.4   3 0.8978
crtRa 0.010203 0.009669 1.0552 1.8   2 0.9449
crtQa 0.008333 0.018634 0.4472 3.6   4 0.9835
crtG  0.003571 0.007986 0.4472 0.2   0 1.0000
cruF  0.000000 0.000000    NaN 1.0   1 1.0000
crtM  0.000000 0.000000    NaN 1.0   1 1.0000

Contrast: soil_hot spring 

average       sd  ratio ava avb cumsum
crtN  0.082307 0.061665 1.3348 6.8 4.2 0.2439
crtP  0.065627 0.047982 1.3677 5.4 4.0 0.4384
crtOa 0.063432 0.055916 1.1344 5.6 3.6 0.6264
crtX  0.029667 0.021762 1.3633 2.0 1.4 0.7143
crtQb 0.025711 0.018286 1.4060 3.4 3.4 0.7905
crtQa 0.023750 0.020287 1.1707 3.6 3.2 0.8608
crtW  0.013750 0.013540 1.0155 0.6 0.6 0.9016
crtRa 0.013680 0.012597 1.0860 1.8 1.0 0.9421
crtRb 0.008290 0.012688 0.6533 1.2 0.8 0.9667
crtG  0.005819 0.009088 0.6403 0.2 0.2 0.9839
crtM  0.005423 0.011358 0.4775 1.0 0.8 1.0000
cruF  0.000000 0.000000    NaN 1.0 1.0 1.0000

Contrast: soil_temperate 

average       sd  ratio ava avb cumsum
crtOa 0.044336 0.026453 1.6760 5.6 3.0 0.1784
crtQa 0.036674 0.017167 2.1363 3.6 1.5 0.3259
crtN  0.035770 0.029686 1.2050 6.8 5.0 0.4698
crtQb 0.033196 0.017234 1.9262 3.4 1.5 0.6033
crtX  0.024966 0.018044 1.3837 2.0 0.5 0.7038
crtP  0.023743 0.023824 0.9966 5.4 4.5 0.7993
crtW  0.021120 0.015240 1.3858 0.6 1.5 0.8842
crtRa 0.010672 0.009567 1.1154 1.8 2.0 0.9272
crtG  0.009061 0.009740 0.9302 0.2 0.5 0.9636
crtRb 0.009045 0.009713 0.9313 1.2 1.5 1.0000
cruF  0.000000 0.000000    NaN 1.0 1.0 1.0000
crtM  0.000000 0.000000    NaN 1.0 1.0 1.0000

Contrast: rock_hot spring 

average       sd  ratio ava avb cumsum
crtN  0.066603 0.044747 1.4884   4 4.2 0.1881
crtP  0.062107 0.021419 2.8997   3 4.0 0.3636
crtOa 0.055884 0.022676 2.4645   3 3.6 0.5214
crtW  0.031338 0.022330 1.4034   2 0.6 0.6100
crtQa 0.029480 0.026701 1.1041   4 3.2 0.6932
crtRb 0.027725 0.019899 1.3933   2 0.8 0.7716
crtQb 0.025088 0.014389 1.7436   3 3.4 0.8424
crtX  0.024546 0.017550 1.3986   1 1.4 0.9118
crtRa 0.021475 0.006922 3.1026   2 1.0 0.9724
crtM  0.006250 0.013975 0.4472   1 0.8 0.9901
crtG  0.003509 0.007846 0.4472   0 0.2 1.0000
cruF  0.000000 0.000000    NaN   1 1.0 1.0000

Contrast: rock_temperate 

average      sd  ratio ava avb cumsum
crtQa 0.05041 0.01356 3.7161   4 1.5 0.2947
crtP  0.03020 0.01385 2.1802   3 4.5 0.4714
crtQb 0.03020 0.01385 2.1802   3 1.5 0.6480
crtN  0.02000 0.02828 0.7071   4 5.0 0.7649
crtRb 0.01020 0.01443 0.7071   2 1.5 0.8246
crtG  0.01000 0.01414 0.7071   0 0.5 0.8831
crtX  0.01000 0.01414 0.7071   1 0.5 0.9415
crtW  0.01000 0.01414 0.7071   2 1.5 1.0000
crtRa 0.00000 0.00000    NaN   2 2.0 1.0000
crtOa 0.00000 0.00000    NaN   3 3.0 1.0000
cruF  0.00000 0.00000    NaN   1 1.0 1.0000
crtM  0.00000 0.00000    NaN   1 1.0 1.0000

Contrast: hot spring_temperate 

average       sd  ratio ava avb cumsum
crtN  0.080122 0.061657 1.2995 4.2 5.0 0.2077
crtP  0.068835 0.054625 1.2601 4.0 4.5 0.3862
crtOa 0.059381 0.023830 2.4919 3.6 3.0 0.5401
crtQb 0.035837 0.025906 1.3833 3.4 1.5 0.6330
crtQa 0.032167 0.022378 1.4374 3.2 1.5 0.7164
crtW  0.025391 0.020232 1.2550 0.6 1.5 0.7823
crtX  0.025046 0.022323 1.1220 1.4 0.5 0.8472
crtRa 0.022812 0.007371 3.0948 1.0 2.0 0.9063
crtRb 0.018045 0.021175 0.8522 0.8 1.5 0.9531
crtG  0.011297 0.012933 0.8735 0.2 0.5 0.9824
crtM  0.006782 0.014299 0.4743 0.8 1.0 1.0000
cruF  0.000000 0.000000    NaN 1.0 1.0 1.0000
Permutation: free
Number of permutations: 0