---
title: "Statistical Analysis of Percent Identities"
author: "Dionne Martin"
date: "6/28/2021"
output: pdf_document
---
  
#################### PREFACE ##########################
#Code for PCA, NMDS, and PERMANOVA/ simper test 
#of the highest percent identity BLASTp hit of criticl 
#carotenoid biosynthesis enzymes: crtE, crtBa, crtBb,
#crtIa, crtIb, cruA, cruP, and crtL.
#Can percent idenitities give us an insight on 
#disimilarities in the pathway between habitats?  
#######################################################

```{pca}
#load packages and blastp raw hits data 
library(devtools)
library(ggbiplot)
library(readr)
percent <- read_csv("/bigdata/stajichlab/shared/projects/Cyanobacteria/Terresterial_Cyano/Cyano Carotenoids BLASTp Results - percent data.csv")
#organize data to just raw hits and strain names 
row.names(percent) = percent$Strain
percent_matrix <- data.matrix(percent[c(3:10)])
range(percent_matrix) #0 100 

#analyze pca statistics to find proportion of variatiance explained by each gene 
pca_percent= prcomp(percent_matrix, center= TRUE, scale.= TRUE)
summary(pca_percent)
str(pca_percent)
#crtE-0.45 crtBa-0.26 crtBb-0.09 crtIa-0.07
#crtIb-0.04 cruA-0.03 cruP-0.01 crtL-0.006


#make groups based on habitat
habitat= percent$Habitat


#try to make group colours, this should be added to plot code 
scale_colour_manual(name= "Habitat", values= c("brown", "dark orange", "blue", "pink", "red", "dark blue", "forest green", "dark yellow", "red", "black"))

#make the PCA plot with strains as points 
#png("PCA Plot of Carotenoid BLAST Hit .png", width=1000, height=1000)
ggbiplot(pca_percent,ellipse=TRUE,obs.scale = 3,var.scale = 1, groups=habitat) +
  scale_colour_manual(name="Habitat", values= c("saddlebrown", "darkorange1", "cornflowerblue", "deeppink", "red", "dark blue", "forest green", "gold", "darkmagenta", "black"))+
  ggtitle("PCA of BLASTp Percent Identities")+
  theme_minimal()+
  theme(legend.position = "bottom")

#make the PCA plot with strain names (looks too messy with all names...)
#png("PCA Plot of Carotenoid BLAST Hit .png", width=1000, height=1000)
ggbiplot(pca_percent,labels= rownames(percent),ellipse=FALSE,obs.scale = 3,var.scale = 1, groups=habitat) +
  scale_colour_manual(name="Habitat", values= c("saddlebrown", "darkorange1", "cornflowerblue", "deeppink", "red", "dark blue", "forest green", "gold", "darkmagenta", "black"))+
  ggtitle("PCA of BLASTp Percent Identities")+
  theme_minimal()+
  theme(legend.position = "bottom")

```
```{nmds}
#load packages 
library(vegan)
library(ggplot2)
library(readr)

#prep data
percent <- read_csv("/bigdata/stajichlab/shared/projects/Cyanobacteria/Terresterial_Cyano/Cyano Carotenoids BLASTp Results - percent data.csv")
percent_matrix <- data.matrix(percent[c(3:10)])
nmds_percent_habitat= as.data.frame(percent[,1:2])
#Do you need to transform your data? range should be 0-10 if not, use root function
range((percent_matrix)^0.5) # 0-10

#make a distance matrix
percent_matrix_ds=vegdist(percent_matrix^0.5, method = 'euclidian')

#Run metaMDS on distance matrix
percent_nms=metaMDS(percent_matrix_ds, distance="euclidian",k=2,trymax=500)
percent_nms
#stress=0.071

#assign colors to habitat
cols <- nmds_percent_habitat$Strain
cols[] <- c("saddlebrown", "darkorange1", "cornflowerblue", "deeppink", "red", "dark blue", "forest green", "gold", "darkmagenta", "black")
#Warning message:
#number of items to replace is not a multiple of replacement length

cols[nmds_percent_habitat$Habitat == c("freshwater")] <- "cornflowerblue"
cols[nmds_percent_habitat$Habitat == c("saltwater")] <- "dark blue"
cols[nmds_percent_habitat$Habitat == c("hot spring")] <- "deeppink"
cols[nmds_percent_habitat$Habitat == c("soil")] <- "forest green"
cols[nmds_percent_habitat$Habitat == c("cave")] <- "saddlebrown"
cols[nmds_percent_habitat$Habitat == c("desert")] <- "darkorange1"
cols[nmds_percent_habitat$Habitat == c("tropical")] <- "black"
cols[nmds_percent_habitat$Habitat == c("temperate")] <- "darkmagenta"
cols[nmds_percent_habitat$Habitat == c("symbiont")] <- "gold"
cols[nmds_percent_habitat$Habitat == c("rock")] <- "red"
pairs(character(nrow(nmds_percent_habitat$Strain, col=cols)
                
#habitatvec = factor(small$Habitat, levels = c("freshwater", "saltwater", "hot spring", "soil", "cave", "desert", "tropical", "temperate", "symbiont", "rock"))
#colors <- c("cornflowerblue", "dark blue", "deeppink", "forest green", "saddlebrown", "darkorange1", "black", "darkmagenta", "gold", "red")
habitatinfo= c(rep(percent$Habitat))
                
#NMDS Plot With Species Names 
ordiplot(percent_nms, type="n")
orditorp(percent_nms,display="sites",labels=percent$Strain, col=cols ,air=0.00025)
par(xpd=TRUE)
#legend(0.6,0.25,legend =unique(nmds_small_habitat$Strain), fill=cols, border= "black")
#change coordinates of the legend for each specific graph
#legend(-0.4,-0.2, "stress = 0.13") 
#change coordinates and stress value for each specific graph
title(main="NMDS of Raw BLASTp")
scores <-
  scores(percent_nms, display = "sites", "species")
cent <-
  aggregate(scores ~ habitatinfo, data = nmds_percent_habitat, FUN = "mean")
names(cent) [-1] <- colnames(scores)
points(cent [,-1],
pch = c( 8 , 8 , 8, 8),
col = unique(colors),
bg = c("black"),
lwd = 3.0,
cex = 2.0
)
                
#NMDS Plot with Points
plot(percent_nms)
points(percent_nms,display= "sites",pch=21, bg=cols) 
par(xpd=TRUE)
#legend(0.6,0.25,legend =habitatvec, fill=colvec, border= "black")
#legend(-0.4,-0.2, "stress = 0.13") 
title(main="NMDS of Raw BLASTp")
scores <-
  scores(percent_nms, display = "sites", "species")
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
rownames(percent_matrix) = percent$Strain
str(percent_matrix)
nmds_percent_habitat= as.data.frame(percent[,1:2])
str(nmds_percent_habitat)
                
#make a distance matrix
percent_matrix_ds=vegdist(percent_matrix^0.5, method = 'euclidian')
                
#Run metaMDS on distance matrix
percent_nms=metaMDS(percent_matrix_ds, distance="euclidian",k=2,trymax=500)
percent_nms
#stress=0.0713
                
#Run PERMANOVA
pmv_percent = adonis(percent_matrix^0.5 ~ habitat, data = nmds_percent_habitat, permutations = 999, method = 'euclidean')
pmv_percent
#PERMANOVA Results: Based off of the carotenoid genes being compared
#cyanobacterial strains differ statistically significant 
#between their habitats (p-value:0.001).
#Habitats explian 34.65% of the variance in the copy numbers of critical 
#carotenoid biosynthesis enzymes between cyanobacterial species. 
                
#Visualise the permanova
densityplot(permustats(pmv_percent))
                
#Distance Based Dispersion test 
bd_percent = betadisper(percent_matrix_ds, percent$Habitat)
boxplot(bd_percent)
anova(bd_percent) 
#F-Test: pvalue 0.1353
permutest(bd_percent)
#permutation test: pvalue 0.392
#the significance of a permutation test is represented by its P-value. 
#The P-value is the probability of obtaining a result at least 
#as extreme as the test statistic given that the null hypothesis is true. 
                
#both the F-Test and permutation test show that 
#the average distance to the median for the 
#habitat data is not statistically significant. 
                
#simper test: determine which genes are responsible for the difference between 2 compared habiatat groups 
simp_percent= simper(percent_matrix^0.5, group= percent$Habitat)
summary(simp_percent)
                
#Results of simper test
Contrast: saltwater_symbiont 

average       sd ratio   ava   avb cumsum
cruP  0.055115 0.030760 1.792 2.289 7.883 0.3026
crtL  0.047584 0.030816 1.544 5.147 0.000 0.5638
cruA  0.044204 0.034086 1.297 3.608 7.740 0.8065
crtBa 0.012597 0.004927 2.557 7.603 8.128 0.8756
crtE  0.012465 0.006994 1.782 7.944 8.252 0.9441
crtBb 0.004084 0.005657 0.722 8.290 8.454 0.9665
crtIa 0.003666 0.002264 1.619 8.473 8.818 0.9866
crtIb 0.002440 0.001395 1.749 8.784 8.547 1.0000

Contrast: saltwater_freshwater 

average       sd  ratio   ava    avb cumsum
cruP  0.052112 0.032343 1.6112 2.289 7.6211 0.3160
crtL  0.046050 0.030700 1.5000 5.147 0.9875 0.5953
cruA  0.042564 0.036561 1.1642 3.608 7.0578 0.8534
crtBa 0.008428 0.004699 1.7934 7.603 8.4828 0.9045
crtBb 0.005247 0.007257 0.7230 8.290 8.5940 0.9363
crtE  0.004981 0.003518 1.4157 7.944 8.4192 0.9665
crtIa 0.002885 0.001620 1.7808 8.473 8.7000 0.9840
crtIb 0.002632 0.002662 0.9888 8.784 8.7609 1.0000

Contrast: saltwater_cave 

average       sd  ratio   ava   avb cumsum
cruP  0.054816 0.031729 1.7276 2.289 7.919 0.2913
crtL  0.047822 0.030935 1.5459 5.147 0.000 0.5455
cruA  0.045052 0.036831 1.2232 3.608 8.207 0.7849
crtBa 0.009962 0.005692 1.7500 7.603 8.624 0.8379
crtIb 0.008993 0.013270 0.6777 8.784 7.918 0.8857
crtIa 0.008639 0.010469 0.8251 8.473 8.204 0.9316
crtE  0.006880 0.004680 1.4699 7.944 8.616 0.9681
crtBb 0.005994 0.006129 0.9779 8.290 7.909 1.0000

Contrast: saltwater_desert 

average       sd  ratio   ava    avb cumsum
cruP  0.052263 0.034722 1.5052 2.289 7.3333 0.2928
crtL  0.045471 0.032339 1.4061 5.147 0.8498 0.5476
cruA  0.043701 0.034830 1.2547 3.608 7.2010 0.7924
crtBa 0.009365 0.005295 1.7686 7.603 8.4536 0.8449
crtBb 0.008210 0.019395 0.4233 8.290 7.7479 0.8909
crtE  0.006904 0.004973 1.3882 7.944 8.4929 0.9295
crtIb 0.006289 0.010495 0.5993 8.784 8.2232 0.9648
crtIa 0.006288 0.009187 0.6844 8.473 8.3781 1.0000

Contrast: saltwater_tropical 

average       sd  ratio   ava   avb cumsum
cruP  0.057427 0.033332 1.7229 2.289 8.480 0.3210
crtL  0.047103 0.030473 1.5457 5.147 0.000 0.5844
cruA  0.045425 0.037398 1.2146 3.608 7.488 0.8383
crtBa 0.011324 0.005160 2.1946 7.603 8.841 0.9016
crtE  0.007667 0.003732 2.0545 7.944 8.784 0.9445
crtBb 0.004372 0.004640 0.9422 8.290 8.029 0.9689
crtIa 0.003266 0.001772 1.8438 8.473 8.777 0.9872
crtIb 0.002289 0.001640 1.3959 8.784 8.592 1.0000

Contrast: saltwater_soil 

average       sd  ratio   ava   avb cumsum
cruP  0.053841 0.032596 1.6518 2.289 8.103 0.3204
crtL  0.046712 0.030170 1.5483 5.147 0.000 0.5985
cruA  0.045331 0.037317 1.2148 3.608 8.459 0.8682
crtBa 0.008259 0.004421 1.8682 7.603 8.492 0.9174
crtE  0.004539 0.002320 1.9568 7.944 8.416 0.9444
crtIa 0.004095 0.002204 1.8584 8.473 8.903 0.9688
crtBb 0.003689 0.005905 0.6246 8.290 8.628 0.9907
crtIb 0.001556 0.001234 1.2602 8.784 8.737 1.0000

Contrast: saltwater_rock 

average        sd  ratio   ava   avb cumsum
cruP  0.053340 0.0332472 1.6043 2.289 7.976 0.3154
crtL  0.047163 0.0309877 1.5220 5.147 0.000 0.5943
cruA  0.042818 0.0375971 1.1389 3.608 8.081 0.8475
crtBa 0.012606 0.0047631 2.6465 7.603 8.981 0.9221
crtE  0.004339 0.0023072 1.8808 7.944 8.391 0.9477
crtIb 0.003447 0.0015644 2.2033 8.784 8.407 0.9681
crtBb 0.003164 0.0051503 0.6144 8.290 8.244 0.9868
crtIa 0.002226 0.0007649 2.9100 8.473 8.604 1.0000

Contrast: saltwater_hot spring 

average       sd  ratio   ava   avb cumsum
cruP  0.049678 0.036552 1.3591 2.289 6.534 0.2845
crtL  0.049195 0.031903 1.5420 5.147 0.000 0.5662
cruA  0.043654 0.038655 1.1293 3.608 6.677 0.8161
crtBb 0.011014 0.012829 0.8585 8.290 8.037 0.8792
crtBa 0.006024 0.004886 1.2329 7.603 7.909 0.9137
crtE  0.005330 0.003731 1.4286 7.944 8.402 0.9442
crtIb 0.005153 0.005190 0.9930 8.784 8.269 0.9737
crtIa 0.004585 0.002761 1.6610 8.473 8.601 1.0000

Contrast: saltwater_temperate 

average       sd  ratio   ava   avb cumsum
cruP  0.054084 0.033396 1.6195 2.289 7.960 0.2967
cruA  0.041035 0.029709 1.3812 3.608 6.904 0.5219
crtL  0.035835 0.032575 1.1001 5.147 3.914 0.7185
crtIb 0.019281 0.016984 1.1352 8.784 6.719 0.8243
crtIa 0.017768 0.016296 1.0903 8.473 6.663 0.9218
crtBa 0.007125 0.004103 1.7365 7.603 8.331 0.9609
crtE  0.003853 0.002436 1.5815 7.944 8.237 0.9820
crtBb 0.003275 0.005227 0.6266 8.290 8.233 1.0000

Contrast: symbiont_freshwater 

average       sd  ratio   ava    avb cumsum
cruA  0.017551 0.024300 0.7222 7.740 7.0578 0.2650
cruP  0.012403 0.014082 0.8808 7.883 7.6211 0.4523
crtE  0.010732 0.008235 1.3033 8.252 8.4192 0.6144
crtBa 0.008486 0.008705 0.9748 8.128 8.4828 0.7425
crtL  0.008199 0.023814 0.3443 0.000 0.9875 0.8664
crtBb 0.004123 0.004919 0.8381 8.454 8.5940 0.9286
crtIb 0.002630 0.002666 0.9866 8.547 8.7609 0.9683
crtIa 0.002097 0.001655 1.2665 8.818 8.7000 1.0000

Contrast: symbiont_cave 

average       sd  ratio   ava   avb cumsum
crtE  0.011372 0.008908 1.2767 8.252 8.616 0.1853
cruP  0.011093 0.009873 1.1236 7.883 7.919 0.3661
crtBa 0.009241 0.009574 0.9652 8.128 8.624 0.5167
cruA  0.008833 0.008828 1.0006 7.740 8.207 0.6607
crtIa 0.007725 0.011690 0.6609 8.818 8.204 0.7866
crtIb 0.007473 0.011958 0.6249 8.547 7.918 0.9084
crtBb 0.005622 0.004923 1.1419 8.454 7.909 1.0000
crtL  0.000000 0.000000    NaN 0.000 0.000 1.0000

Contrast: symbiont_desert 

average       sd  ratio   ava    avb cumsum
cruP  0.017363 0.025254 0.6876 7.883 7.3333 0.2182
cruA  0.016113 0.021763 0.7404 7.740 7.2010 0.4206
crtE  0.011442 0.009397 1.2176 8.252 8.4929 0.5644
crtBa 0.009031 0.009695 0.9315 8.128 8.4536 0.6779
crtBb 0.007837 0.018095 0.4331 8.454 7.7479 0.7763
crtL  0.007457 0.021257 0.3508 0.000 0.8498 0.8700
crtIa 0.005538 0.010034 0.5519 8.818 8.3781 0.9396
crtIb 0.004805 0.009327 0.5152 8.547 8.2232 1.0000

Contrast: symbiont_tropical 

average        sd  ratio   ava   avb cumsum
cruA  0.018315 0.0249013 0.7355 7.740 7.488 0.3430
crtE  0.010671 0.0093441 1.1420 8.252 8.784 0.5428
cruP  0.008714 0.0103803 0.8394 7.883 8.480 0.7060
crtBa 0.008387 0.0105204 0.7972 8.128 8.841 0.8630
crtBb 0.004010 0.0024142 1.6610 8.454 8.029 0.9381
crtIa 0.001975 0.0014567 1.3555 8.818 8.777 0.9751
crtIb 0.001329 0.0009588 1.3857 8.547 8.592 1.0000
crtL  0.000000 0.0000000    NaN 0.000 0.000 1.0000

Contrast: symbiont_soil 

average       sd  ratio   ava   avb cumsum
crtE  0.010260 0.008162 1.2571 8.252 8.416 0.2453
cruP  0.009028 0.007924 1.1393 7.883 8.103 0.4612
cruA  0.008102 0.009914 0.8172 7.740 8.459 0.6549
crtBa 0.008084 0.008834 0.9151 8.128 8.492 0.8482
crtBb 0.002551 0.002297 1.1103 8.454 8.628 0.9092
crtIa 0.001924 0.001622 1.1859 8.818 8.903 0.9552
crtIb 0.001872 0.001032 1.8146 8.547 8.737 1.0000
crtL  0.000000 0.000000    NaN 0.000 0.000 1.0000

Contrast: symbiont_rock 

average        sd  ratio   ava   avb cumsum
crtE  0.010348 0.0091886 1.1262 8.252 8.391 0.2450
cruP  0.009633 0.0078913 1.2207 7.883 7.976 0.4732
cruA  0.008046 0.0092156 0.8730 7.740 8.081 0.6637
crtBa 0.007580 0.0132189 0.5734 8.128 8.981 0.8432
crtBb 0.003080 0.0013147 2.3428 8.454 8.244 0.9161
crtIa 0.002327 0.0014392 1.6170 8.818 8.604 0.9712
crtIb 0.001214 0.0009905 1.2261 8.547 8.407 1.0000
crtL  0.000000 0.0000000    NaN 0.000 0.000 1.0000

Contrast: symbiont_hot spring 

average       sd  ratio   ava   avb cumsum
cruP  0.022239 0.027562 0.8069 7.883 6.534 0.2674
cruA  0.021201 0.028369 0.7473 7.740 6.677 0.5223
crtE  0.011360 0.008625 1.3171 8.252 8.402 0.6589
crtBa 0.011347 0.006683 1.6979 8.128 7.909 0.7954
crtBb 0.009608 0.012399 0.7749 8.454 8.037 0.9109
crtIa 0.004006 0.003646 1.0988 8.818 8.601 0.9591
crtIb 0.003403 0.004475 0.7604 8.547 8.269 1.0000
crtL  0.000000 0.000000    NaN 0.000 0.000 1.0000

Contrast: symbiont_temperate 

average       sd  ratio   ava   avb cumsum
crtL  0.034534 0.036950 0.9346 0.000 3.914 0.2971
crtIa 0.018960 0.017390 1.0903 8.818 6.663 0.4602
crtIb 0.016399 0.016692 0.9824 8.547 6.719 0.6012
cruA  0.013318 0.011451 1.1631 7.740 6.904 0.7158
crtE  0.011087 0.008007 1.3847 8.252 8.237 0.8111
cruP  0.009852 0.007299 1.3498 7.883 7.960 0.8959
crtBa 0.008924 0.008046 1.1090 8.128 8.331 0.9727
crtBb 0.003179 0.001550 2.0503 8.454 8.233 1.0000

Contrast: freshwater_cave 

average       sd  ratio    ava   avb cumsum
cruA  0.014771 0.025599 0.5770 7.0578 8.207 0.2236
cruP  0.009977 0.013913 0.7171 7.6211 7.919 0.3747
crtIb 0.008725 0.012318 0.7083 8.7609 7.918 0.5068
crtL  0.008236 0.023867 0.3451 0.9875 0.000 0.6315
crtBb 0.007459 0.005970 1.2494 8.5940 7.909 0.7444
crtIa 0.007312 0.011161 0.6551 8.7000 8.204 0.8551
crtE  0.004971 0.003447 1.4419 8.4192 8.616 0.9304
crtBa 0.004597 0.003238 1.4196 8.4828 8.624 1.0000

Contrast: freshwater_desert 

average       sd  ratio    ava    avb cumsum
cruA  0.021405 0.028827 0.7425 7.0578 7.2010 0.2602
cruP  0.016596 0.026242 0.6324 7.6211 7.3333 0.4619
crtL  0.014178 0.028757 0.4930 0.9875 0.8498 0.6342
crtBb 0.009629 0.018235 0.5281 8.5940 7.7479 0.7512
crtIb 0.006181 0.009788 0.6315 8.7609 8.2232 0.8263
crtIa 0.004981 0.009676 0.5148 8.7000 8.3781 0.8869
crtE  0.004903 0.004916 0.9974 8.4192 8.4929 0.9465
crtBa 0.004403 0.004408 0.9989 8.4828 8.4536 1.0000

Contrast: freshwater_tropical 

average       sd  ratio    ava   avb cumsum
cruA  0.022496 0.031617 0.7115 7.0578 7.488 0.3853
cruP  0.008848 0.014022 0.6310 7.6211 8.480 0.5369
crtL  0.008124 0.023558 0.3449 0.9875 0.000 0.6760
crtBb 0.006114 0.004102 1.4905 8.5940 8.029 0.7807
crtE  0.004438 0.002824 1.5715 8.4192 8.784 0.8568
crtBa 0.004134 0.002660 1.5539 8.4828 8.841 0.9276
crtIb 0.002563 0.002737 0.9367 8.7609 8.592 0.9715
crtIa 0.001666 0.001344 1.2398 8.7000 8.777 1.0000

Contrast: freshwater_soil 

average       sd  ratio    ava   avb cumsum
cruA  0.013907 0.026230 0.5302 7.0578 8.459 0.3380
crtL  0.008064 0.023392 0.3447 0.9875 0.000 0.5340
cruP  0.006975 0.013165 0.5298 7.6211 8.103 0.7035
crtBb 0.003069 0.004826 0.6359 8.5940 8.628 0.7781
crtE  0.002679 0.001877 1.4275 8.4192 8.416 0.8432
crtBa 0.002333 0.001799 1.2972 8.4828 8.492 0.8999
crtIb 0.002210 0.002548 0.8672 8.7609 8.737 0.9536
crtIa 0.001910 0.001855 1.0296 8.7000 8.903 1.0000

Contrast: freshwater_rock 

average        sd  ratio    ava   avb cumsum
cruA  0.013383 0.0256456 0.5219 7.0578 8.081 0.3004
crtL  0.008135 0.0239425 0.3398 0.9875 0.000 0.4830
cruP  0.006707 0.0133175 0.5036 7.6211 7.976 0.6335
crtBb 0.004582 0.0043667 1.0493 8.5940 8.244 0.7363
crtBa 0.004419 0.0025088 1.7614 8.4828 8.981 0.8355
crtIb 0.003278 0.0028617 1.1455 8.7609 8.407 0.9091
crtE  0.002570 0.0019794 1.2983 8.4192 8.391 0.9668
crtIa 0.001480 0.0008161 1.8141 8.7000 8.604 1.0000

Contrast: freshwater_hot spring 

average       sd  ratio    ava   avb cumsum
cruA  0.024322 0.033831 0.7189 7.0578 6.677 0.2922
cruP  0.020459 0.028461 0.7188 7.6211 6.534 0.5379
crtBb 0.009921 0.013190 0.7521 8.5940 8.037 0.6571
crtL  0.008448 0.024526 0.3445 0.9875 0.000 0.7586
crtBa 0.007116 0.004105 1.7336 8.4828 7.909 0.8440
crtIb 0.005026 0.005236 0.9601 8.7609 8.269 0.9044
crtE  0.004133 0.003125 1.3227 8.4192 8.402 0.9541
crtIa 0.003825 0.003014 1.2690 8.7000 8.601 1.0000

Contrast: freshwater_temperate 

average       sd  ratio    ava   avb cumsum
crtL  0.035469 0.034704 1.0220 0.9875 3.914 0.3254
cruA  0.020040 0.021170 0.9466 7.0578 6.904 0.5093
crtIb 0.018143 0.015814 1.1473 8.7609 6.719 0.6758
crtIa 0.017976 0.016121 1.1151 8.7000 6.663 0.8407
cruP  0.006767 0.013403 0.5049 7.6211 7.960 0.9028
crtBb 0.004759 0.004426 1.0751 8.5940 8.233 0.9464
crtE  0.003508 0.002545 1.3783 8.4192 8.237 0.9786
crtBa 0.002331 0.001828 1.2750 8.4828 8.331 1.0000

Contrast: cave_desert 

average       sd  ratio   ava    avb cumsum
cruP  0.016310 0.024974 0.6531 7.919 7.3333 0.2142
cruA  0.014305 0.022297 0.6416 8.207 7.2010 0.4020
crtIb 0.009675 0.012862 0.7522 7.918 8.2232 0.5290
crtIa 0.009619 0.012948 0.7429 8.204 8.3781 0.6553
crtL  0.007492 0.021310 0.3516 0.000 0.8498 0.7537
crtBb 0.007234 0.017188 0.4209 7.909 7.7479 0.8487
crtE  0.005863 0.005473 1.0712 8.616 8.4929 0.9257
crtBa 0.005659 0.005187 1.0911 8.624 8.4536 1.0000

Contrast: cave_tropical 

average       sd  ratio   ava   avb cumsum
cruA  0.015583 0.025550 0.6099 8.207 7.488 0.3058
cruP  0.007721 0.007320 1.0547 7.919 8.480 0.4574
crtIb 0.007645 0.011817 0.6469 7.918 8.592 0.6074
crtIa 0.007348 0.011434 0.6426 8.204 8.777 0.7516
crtE  0.004751 0.003441 1.3807 8.616 8.784 0.8448
crtBa 0.004642 0.003841 1.2085 8.624 8.841 0.9359
crtBb 0.003264 0.003823 0.8536 7.909 8.029 1.0000
crtL  0.000000 0.000000    NaN 0.000 0.000 1.0000

Contrast: cave_soil 

average       sd  ratio   ava   avb cumsum
crtIb 0.007996 0.012265 0.6519 7.918 8.737 0.1961
crtIa 0.007506 0.011853 0.6333 8.204 8.903 0.3802
cruP  0.006785 0.005676 1.1953 7.919 8.103 0.5466
crtBb 0.006380 0.004947 1.2898 7.909 8.628 0.7031
crtBa 0.004195 0.002839 1.4775 8.624 8.492 0.8060
crtE  0.004180 0.002890 1.4466 8.616 8.416 0.9085
cruA  0.003731 0.003327 1.1216 8.207 8.459 1.0000
crtL  0.000000 0.000000    NaN 0.000 0.000 1.0000

Contrast: cave_rock 

average       sd  ratio   ava   avb cumsum
crtIb 0.007271 0.011721 0.6203 7.918 8.407 0.1943
crtIa 0.007265 0.011246 0.6460 8.204 8.604 0.3885
cruP  0.006814 0.005713 1.1928 7.919 7.976 0.5707
crtE  0.004279 0.002976 1.4380 8.616 8.391 0.6850
crtBa 0.004260 0.004323 0.9854 8.624 8.981 0.7989
crtBb 0.003884 0.004364 0.8900 7.909 8.244 0.9027
cruA  0.003640 0.001967 1.8506 8.207 8.081 1.0000
crtL  0.000000 0.000000    NaN 0.000 0.000 1.0000

Contrast: cave_hot spring 

average       sd  ratio   ava   avb cumsum
cruP  0.020778 0.027530 0.7547 7.919 6.534 0.2493
cruA  0.018470 0.029960 0.6165 8.207 6.677 0.4709
crtBb 0.012442 0.009509 1.3084 7.909 8.037 0.6202
crtIa 0.009013 0.011065 0.8145 8.204 8.601 0.7283
crtIb 0.008919 0.011252 0.7927 7.918 8.269 0.8353
crtBa 0.008282 0.005530 1.4975 8.624 7.909 0.9347
crtE  0.005444 0.003931 1.3847 8.616 8.402 1.0000
crtL  0.000000 0.000000    NaN 0.000 0.000 1.0000

Contrast: cave_temperate 

average       sd  ratio   ava   avb cumsum
crtL  0.034697 0.035665 0.9728 0.000 3.914 0.3364
crtIa 0.018535 0.015631 1.1858 8.204 6.663 0.5161
crtIb 0.016698 0.015242 1.0956 7.918 6.719 0.6780
cruA  0.012888 0.010047 1.2829 8.207 6.904 0.8029
cruP  0.006913 0.005647 1.2243 7.919 7.960 0.8699
crtE  0.005103 0.003761 1.3568 8.616 8.237 0.9194
crtBa 0.004426 0.003173 1.3951 8.624 8.331 0.9623
crtBb 0.003889 0.004369 0.8901 7.909 8.233 1.0000

Contrast: desert_tropical 

average       sd  ratio    ava   avb cumsum
cruA  0.022412 0.029555 0.7583 7.2010 7.488 0.3237
cruP  0.013973 0.026656 0.5242 7.3333 8.480 0.5255
crtL  0.007385 0.021016 0.3514 0.8498 0.000 0.6322
crtBb 0.005384 0.017227 0.3125 7.7479 8.029 0.7099
crtIb 0.005067 0.009219 0.5496 8.2232 8.592 0.7831
crtIa 0.005058 0.009860 0.5130 8.3781 8.777 0.8562
crtBa 0.005021 0.005223 0.9615 8.4536 8.841 0.9287
crtE  0.004937 0.005358 0.9215 8.4929 8.784 1.0000

Contrast: desert_soil 

average       sd  ratio    ava   avb cumsum
cruA  0.013517 0.023030 0.5869 7.2010 8.459 0.2231
cruP  0.012628 0.025424 0.4967 7.3333 8.103 0.4316
crtBb 0.008507 0.017885 0.4757 7.7479 8.628 0.5720
crtL  0.007327 0.020855 0.3513 0.8498 0.000 0.6929
crtIb 0.005472 0.009539 0.5737 8.2232 8.737 0.7832
crtIa 0.005354 0.010176 0.5262 8.3781 8.903 0.8716
crtE  0.003945 0.004627 0.8526 8.4929 8.416 0.9367
crtBa 0.003833 0.004265 0.8986 8.4536 8.492 1.0000

Contrast: desert_rock 

average       sd  ratio    ava   avb cumsum
cruA  0.013382 0.021980 0.6088 7.2010 8.081 0.2332
cruP  0.012963 0.025564 0.5071 7.3333 7.976 0.4590
crtL  0.007395 0.021369 0.3460 0.8498 0.000 0.5879
crtBb 0.005446 0.018045 0.3018 7.7479 8.244 0.6828
crtBa 0.004931 0.005687 0.8670 8.4536 8.981 0.7687
crtIa 0.004769 0.009500 0.5019 8.3781 8.604 0.8518
crtIb 0.004534 0.008979 0.5049 8.2232 8.407 0.9308
crtE  0.003974 0.004715 0.8430 8.4929 8.391 1.0000

Contrast: desert_hot spring 

average       sd  ratio    ava   avb cumsum
cruP  0.025035 0.033553 0.7461 7.3333 6.534 0.2523
cruA  0.024642 0.031734 0.7765 7.2010 6.677 0.5007
crtBb 0.014810 0.019112 0.7749 7.7479 8.037 0.6500
crtBa 0.008033 0.005124 1.5677 8.4536 7.909 0.7310
crtL  0.007697 0.021938 0.3508 0.8498 0.000 0.8085
crtIa 0.006978 0.009635 0.7243 8.3781 8.601 0.8789
crtIb 0.006594 0.009147 0.7209 8.2232 8.269 0.9453
crtE  0.005422 0.005274 1.0280 8.4929 8.402 1.0000

Contrast: desert_temperate 

average       sd  ratio    ava   avb cumsum
crtL  0.034945 0.035321 0.9894 0.8498 3.914 0.3019
crtIa 0.018281 0.015730 1.1621 8.3781 6.663 0.4599
cruA  0.017637 0.019415 0.9084 7.2010 6.904 0.6123
crtIb 0.016827 0.015164 1.1097 8.2232 6.719 0.7576
cruP  0.013259 0.025710 0.5157 7.3333 7.960 0.8722
crtBb 0.005658 0.018181 0.3112 7.7479 8.233 0.9211
crtE  0.005035 0.004779 1.0536 8.4929 8.237 0.9646
crtBa 0.004100 0.004015 1.0211 8.4536 8.331 1.0000

Contrast: tropical_soil 

average       sd  ratio   ava   avb cumsum
cruA  0.013933 0.026473 0.5263 7.488 8.459 0.4130
crtBb 0.005045 0.001635 3.0861 8.029 8.628 0.5625
cruP  0.004138 0.002038 2.0299 8.480 8.103 0.6851
crtBa 0.003706 0.002057 1.8018 8.841 8.492 0.7949
crtE  0.003433 0.002330 1.4736 8.784 8.416 0.8967
crtIb 0.001775 0.001295 1.3701 8.592 8.737 0.9493
crtIa 0.001711 0.001455 1.1764 8.777 8.903 1.0000
crtL  0.000000 0.000000    NaN 0.000 0.000 1.0000

Contrast: tropical_rock 

average        sd  ratio   ava   avb cumsum
cruA  0.015389 0.0263808 0.5834 7.488 8.081 0.4884
cruP  0.004611 0.0030054 1.5342 8.480 7.976 0.6348
crtE  0.003604 0.0024595 1.4655 8.784 8.391 0.7491
crtBa 0.002338 0.0024009 0.9738 8.841 8.981 0.8233
crtBb 0.002040 0.0006205 3.2881 8.029 8.244 0.8881
crtIb 0.001834 0.0009064 2.0240 8.592 8.407 0.9463
crtIa 0.001692 0.0010970 1.5419 8.777 8.604 1.0000
crtL  0.000000 0.0000000    NaN 0.000 0.000 1.0000

Contrast: tropical_hot spring 

average       sd  ratio   ava   avb cumsum
cruA  0.025616 0.034881 0.7344 7.488 6.677 0.3279
cruP  0.019516 0.029642 0.6584 8.480 6.534 0.5777
crtBb 0.011534 0.009050 1.2746 8.029 8.037 0.7254
crtBa 0.009083 0.005438 1.6704 8.841 7.909 0.8417
crtE  0.004878 0.003531 1.3814 8.784 8.402 0.9041
crtIa 0.003781 0.003292 1.1487 8.777 8.601 0.9525
crtIb 0.003708 0.004427 0.8376 8.592 8.269 1.0000
crtL  0.000000 0.000000    NaN 0.000 0.000 1.0000

Contrast: tropical_temperate 

average        sd  ratio   ava   avb cumsum
crtL  0.034201 0.0355444 0.9622 0.000 3.914 0.3158
cruA  0.022751 0.0211064 1.0779 7.488 6.904 0.5258
crtIa 0.018401 0.0166793 1.1032 8.777 6.663 0.6957
crtIb 0.016609 0.0159932 1.0385 8.592 6.719 0.8491
crtE  0.004997 0.0033085 1.5103 8.784 8.237 0.8952
cruP  0.004819 0.0029548 1.6310 8.480 7.960 0.9397
crtBa 0.004534 0.0025552 1.7742 8.841 8.331 0.9816
crtBb 0.001997 0.0009879 2.0214 8.029 8.233 1.0000

Contrast: soil_rock 

average        sd ratio   ava   avb cumsum
crtBa 0.0041277 0.0015130 2.728 8.492 8.981 0.2342
crtBb 0.0032441 0.0013454 2.411 8.628 8.244 0.4183
cruA  0.0031883 0.0019269 1.655 8.459 8.081 0.5993
crtIb 0.0027858 0.0012648 2.203 8.737 8.407 0.7574
crtIa 0.0025213 0.0014874 1.695 8.903 8.604 0.9004
cruP  0.0010691 0.0004520 2.365 8.103 7.976 0.9611
crtE  0.0006855 0.0002865 2.392 8.416 8.391 1.0000
crtL  0.0000000 0.0000000   NaN 0.000 0.000 1.0000

Contrast: soil_hot spring 

average       sd  ratio   ava   avb cumsum
cruP  0.017708 0.028572 0.6198 8.103 6.534 0.2847
cruA  0.017466 0.031023 0.5630 8.459 6.677 0.5656
crtBb 0.008512 0.013122 0.6487 8.628 8.037 0.7024
crtBa 0.006950 0.003770 1.8436 8.492 7.909 0.8142
crtIb 0.004391 0.004773 0.9200 8.737 8.269 0.8848
crtIa 0.003905 0.003838 1.0174 8.903 8.601 0.9476
crtE  0.003261 0.002316 1.4076 8.416 8.402 1.0000
crtL  0.000000 0.000000    NaN 0.000 0.000 1.0000

Contrast: soil_temperate 

average        sd  ratio   ava   avb cumsum
crtL  0.033932 0.0357682 0.9487 0.000 3.914 0.3660
crtIa 0.019373 0.0168360 1.1507 8.903 6.663 0.5749
crtIb 0.017515 0.0164139 1.0671 8.737 6.719 0.7639
cruA  0.013447 0.0113958 1.1800 8.459 6.904 0.9089
crtBb 0.003395 0.0015692 2.1635 8.628 8.233 0.9455
crtE  0.002359 0.0018059 1.3062 8.416 8.237 0.9709
crtBa 0.001463 0.0013439 1.0886 8.492 8.331 0.9867
cruP  0.001231 0.0006338 1.9423 8.103 7.960 1.0000

Contrast: rock_hot spring 

average       sd  ratio   ava   avb cumsum
cruP  0.017632 0.031120 0.5666 7.976 6.534 0.2715
cruA  0.016903 0.032761 0.5160 8.081 6.677 0.5319
crtBb 0.010422 0.011558 0.9018 8.244 8.037 0.6924
crtBa 0.009724 0.006538 1.4872 8.981 7.909 0.8421
crtIa 0.003916 0.002471 1.5849 8.604 8.601 0.9024
crtE  0.003198 0.002549 1.2549 8.391 8.402 0.9517
crtIb 0.003136 0.004271 0.7342 8.407 8.269 1.0000
crtL  0.000000 0.000000    NaN 0.000 0.000 1.0000

Contrast: rock_temperate 

average        sd   ratio   ava   avb cumsum
crtL  0.0342450 0.0484298  0.7071 0.000 3.914 0.3942
crtIa 0.0169692 0.0226762  0.7483 8.604 6.663 0.5895
crtIb 0.0156844 0.0209023  0.7504 8.407 6.719 0.7700
cruA  0.0107323 0.0145739  0.7364 8.081 6.904 0.8935
crtBa 0.0056206 0.0001301 43.2022 8.981 8.331 0.9582
crtE  0.0023778 0.0019182  1.2396 8.391 8.237 0.9856
crtBb 0.0008205 0.0001465  5.6016 8.244 8.233 0.9950
cruP  0.0004308 0.0002086  2.0657 7.976 7.960 1.0000

Contrast: hot spring_temperate 

average       sd  ratio   ava   avb cumsum
crtL  0.035644 0.037649 0.9467 0.000 3.914 0.2658
cruA  0.023468 0.024927 0.9415 6.677 6.904 0.4407
crtIa 0.019178 0.016307 1.1761 8.601 6.663 0.5837
cruP  0.017997 0.029684 0.6063 6.534 7.960 0.7179
crtIb 0.016861 0.015067 1.1191 8.269 6.719 0.8436
crtBb 0.010651 0.011016 0.9668 8.037 8.233 0.9230
crtBa 0.006222 0.003328 1.8694 7.909 8.331 0.9694
crtE  0.004101 0.002750 1.4912 8.402 8.237 1.0000
Permutation: free
Number of permutations: 0

                