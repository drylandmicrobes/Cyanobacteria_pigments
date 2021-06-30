---
title: "Statistical Analysis of All pHHMer Results"
author: "Dionne Martin"
date: "6/28/2021"
output: pdf_document
---
  
#################### PREFACE ##########################
#Code for PCA, NMDS, and PERMANOVA/ simper test 
#of pHMMer results. ALL ENZYMES TEST RUN 
#carotenoid biosynthesis enzymes:
#Can percent idenitities give us an insight on 
#disimilarities in the pathway between habitats?  
#######################################################

```{pca}
#load packages and blastp raw hits data 
library(devtools)
library(ggbiplot)
library(readr)
pHMMER <- read_csv("/bigdata/stajichlab/shared/projects/Cyanobacteria/Terresterial_Cyano/Cyano Carotenoids pHMMER Results - phmmer (1).csv")
#organize data to just raw hits and strain names 
#row.names(pHMMER) = pHMMER$Strain
#pHMMER_matrix <- data.matrix(pHMMER[c(3:31)])
pHMMER_matrix <- subset(pHMMER, select=c(crtE, crtBa, crtBb, crtIa, crtIb, cruA, cruP, crtL, crtRa, crtRb, crtG, crtX, crtW, crtOa, crtN,crtP, crtQa, crtQb, cruF, crtM),row.names=1)
range(pHMMER_matrix^0.5) #0 23

#analyze pca statistics to find proportion of variatiance explained by each gene 
pca_pHMMER= prcomp(pHMMER_matrix^0.5, center= TRUE, scale.= TRUE)
summary(pca_pHMMER)
str(pca_pHMMER)
#crtE-0.45 crtBa-0.26 crtBb-0.09 crtIa-0.07
#crtIb-0.04 cruA-0.03 cruP-0.01 crtL-0.006


#make groups based on habitat
habitat= pHMMER$Habitat


#try to make group colours, this should be added to plot code 
scale_colour_manual(name= "Habitat", values= c("brown", "dark orange", "blue", "pink", "red", "dark blue", "forest green", "dark yellow", "red", "black"))

#make the PCA plot with strains as points NEED TO DO PCA OF SEPERATE SETS OF GENES 
#png("PCA Plot of Carotenoid BLAST Hit .png", width=1000, height=1000)
ggbiplot(pca_pHMMER,ellipse=TRUE,obs.scale = 3,var.scale = 1, groups=habitat) +
  scale_colour_manual(name="Habitat", values= c("saddlebrown", "darkorange1", "cornflowerblue", "deeppink", "red", "dark blue", "forest green", "gold", "darkmagenta", "black"))+
  ggtitle("PCA of pHMMer Results")+
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
pHMMER <- read_csv("/bigdata/stajichlab/shared/projects/Cyanobacteria/Terresterial_Cyano/Cyano Carotenoids pHMMER Results - phmmer (1).csv")
pHMMER_matrix <- subset(pHMMER, select=c(crtE, crtBa, crtBb, crtIa, crtIb, cruA, cruP, crtL, crtRa, crtRb, crtG, crtX, crtW, crtOa, crtN,crtP, crtQa, crtQb, cruF, crtM),row.names=1)
nmds_pHMMER_habitat= as.data.frame(pHMMER[,1:2])
#Do you need to transform your data? range should be 0-10 if not, use root function
range(pHMMER_matrix^0.5) # 0-23

#make a distance matrix
pHMMER_matrix_ds=vegdist(pHMMER_matrix^0.5, method = 'euclidian')

#Run metaMDS on distance matrix
pHMMER_nms=metaMDS(pHMMER_matrix_ds, distance="euclidian",k=2,trymax=500)
pHMMER_nms
#stress=0.129

#assign colors to habitat
cols <- nmds_pHMMER_habitat$Strain
cols[] <- c("saddlebrown", "darkorange1", "cornflowerblue", "deeppink", "red", "dark blue", "forest green", "gold", "darkmagenta", "black")
#Warning message:
#number of items to replace is not a multiple of replacement length

cols[nmds_pHMMER_habitat$Habitat == c("freshwater")] <- "cornflowerblue"
cols[nmds_pHMMER_habitat$Habitat == c("saltwater")] <- "dark blue"
cols[nmds_pHMMER_habitat$Habitat == c("hot spring")] <- "deeppink"
cols[nmds_pHMMER_habitat$Habitat == c("soil")] <- "forest green"
cols[nmds_pHMMER_habitat$Habitat == c("cave")] <- "saddlebrown"
cols[nmds_pHMMER_habitat$Habitat == c("desert")] <- "darkorange1"
cols[nmds_pHMMER_habitat$Habitat == c("tropical")] <- "black"
cols[nmds_pHMMER_habitat$Habitat == c("temperate")] <- "darkmagenta"
cols[nmds_pHMMER_habitat$Habitat == c("symbiont")] <- "gold"
cols[nmds_pHMMER_habitat$Habitat == c("rock")] <- "red"
pairs(character(nrow(nmds_pHMMER_habitat$Strain, col=cols)
                
#habitatvec = factor(small$Habitat, levels = c("freshwater", "saltwater", "hot spring", "soil", "cave", "desert", "tropical", "temperate", "symbiont", "rock"))
#colors <- c("cornflowerblue", "dark blue", "deeppink", "forest green", "saddlebrown", "darkorange1", "black", "darkmagenta", "gold", "red")
habitatinfo= c(rep(pHMMER$Habitat))
                
#NMDS Plot With Species Names 
ordiplot(pHMMER_nms, type="n")
orditorp(pHMMER_nms,display="sites",labels=pHMMER$Strain, col=cols ,air=0.00025)
par(xpd=TRUE)
#legend(0.6,0.25,legend =unique(nmds_small_habitat$Strain), fill=cols, border= "black")
#change coordinates of the legend for each specific graph
#legend(-0.4,-0.2, "stress = 0.13") 
#change coordinates and stress value for each specific graph
title(main="NMDS of pHMMer Results")
scores <-
  scores(pHMMER_nms, display = "sites", "species")
cent <-
  aggregate(scores ~ habitatinfo, data = nmds_pHMMER_habitat, FUN = "mean")
names(cent) [-1] <- colnames(scores)
points(cent [,-1],
pch = c( 8 , 8 , 8, 8),
col = unique(colors),
bg = c("black"),
lwd = 3.0,
cex = 2.0
)
                
                #NMDS Plot with Points
                plot(pHMMER_nms)
                points(pHMMER_nms,display= "sites",pch=21, bg=cols) 
                par(xpd=TRUE)
                #legend(0.6,0.25,legend =habitatvec, fill=colvec, border= "black")
                #legend(-0.4,-0.2, "stress = 0.13") 
                title(main="NMDS of Raw BLASTp")
                scores <-
                  scores(pHMMER_nms, display = "sites", "species")
                cent <-
                  aggregate(scores ~ habitatinfo, data = nmds_pHMMER_habitat, FUN = "mean")
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
                rownames(pHMMER_matrix) = pHMMER$Strain
                str(pHMMER_matrix)
                nmds_pHMMER_habitat= as.data.frame(pHMMER[,1:2])
                str(nmds_pHMMER_habitat)
                
                #make a distance matrix
                pHMMER_matrix_ds=vegdist(pHMMER_matrix, method = 'euclidian')
                
                #Run metaMDS on distance matrix
                pHMMER_nms=metaMDS(pHMMER_matrix_ds, distance="euclidian",k=2,trymax=500)
                pHMMER_nms
                #stress=0.102
                
                #Run PERMANOVA
                pmv_pHMMER = adonis(pHMMER_matrix ~ habitat, data = nmds_pHMMER_habitat, permutations = 999, method = 'euclidean')
                pmv_pHMMER
                #PERMANOVA Results: Based off of the carotenoid genes being compared
                #cyanobacterial strains differ statistically significant 
                #between their habitats (p-value:0.001).
                #Habitats explian 20.65% of the variance in the copy numbers of critical 
                #carotenoid biosynthesis enzymes between cyanobacterial species. 
                
                #Visualise the permanova
                densityplot(permustats(pmv_pHMMER))
                
                #Distance Based Dispersion test 
                bd_pHMMER = betadisper(pHMMER_matrix_ds, pHMMER$Habitat)
                boxplot(bd_pHMMER)
                anova(bd_pHMMER) 
                #F-Test: pvalue 0.002305
                permutest(bd_pHMMER)
                #permutation test: pvalue 0.007
                #the significance of a permutation test is represented by its P-value. 
                #The P-value is the probability of obtaining a result at least 
                #as extreme as the test statistic given that the null hypothesis is true. 
                
                #both the F-Test and permutation test show that 
                #the avergae distance to the median for the 
                #habitat data is not statistically significant. 
                
                #simper test: determine which genes are responsible for the difference between 2 compared habiatat groups 
                simp_pHMMER= simper(pHMMER_matrix, group= pHMMER$Habitat)
                summary(simp_pHMMER)
                
                #Results of simper test
                crtIb 0.019606 0.018979 1.0330 2.9583  4.7037 0.6051
                crtQb 0.019163 0.015623 1.2266 3.0000  4.6667 0.6662
                cruF  0.019163 0.015623 1.2266 3.0000  4.6667 0.7273
                crtBb 0.013183 0.009606 1.3724 3.0000  1.6296 0.7694
                cruA  0.011671 0.008102 1.4406 0.6250  1.6667 0.8066
                cruP  0.011125 0.007658 1.4527 0.6250  1.5926 0.8420
                crtRb 0.010584 0.011236 0.9419 1.2500  2.2222 0.8758
                crtM  0.009586 0.008623 1.1117 1.2917  2.0741 0.9064
                crtL  0.009205 0.007530 1.2225 1.0833  0.1111 0.9357
                crtW  0.008147 0.006691 1.2175 0.2083  1.1111 0.9617
                crtBa 0.003757 0.005805 0.6472 1.0833  1.3333 0.9737
                crtRa 0.003475 0.008123 0.4277 1.0833  1.2963 0.9847
                crtG  0.002980 0.004700 0.6340 0.1250  0.2963 0.9942
                crtE  0.001809 0.003774 0.4794 2.0833  2.1852 1.0000
                
                Contrast: saltwater_tropical 
                
                average       sd  ratio    ava    avb cumsum
                crtIa 0.035511 0.020838 1.7041 3.7500  8.143 0.1051
                crtOa 0.031380 0.022511 1.3939 8.0000 10.429 0.1979
                crtN  0.028064 0.022693 1.2367 7.0833 10.000 0.2810
                crtQa 0.027599 0.020371 1.3548 3.0000  6.571 0.3627
                crtP  0.025463 0.021629 1.1773 7.2500  9.286 0.4380
                crtQb 0.024514 0.019980 1.2269 3.0000  6.143 0.5106
                cruF  0.024514 0.019980 1.2269 3.0000  6.143 0.5831
                crtX  0.023368 0.013381 1.7464 0.3333  3.429 0.6523
                crtIb 0.019057 0.015481 1.2310 2.9583  5.143 0.7086
                cruA  0.015641 0.011372 1.3754 0.6250  2.429 0.7549
                cruP  0.013808 0.008571 1.6109 0.6250  2.286 0.7958
                crtRb 0.012512 0.011812 1.0592 1.2500  2.857 0.8328
                crtBb 0.010669 0.008063 1.3232 3.0000  2.000 0.8644
                crtW  0.009781 0.007399 1.3218 0.2083  1.429 0.8933
                crtG  0.008788 0.007842 1.1207 0.1250  1.286 0.9193
                crtL  0.008741 0.006993 1.2499 1.0833  0.000 0.9452
                crtRa 0.006645 0.008248 0.8057 1.0833  2.000 0.9649
                crtM  0.006638 0.007657 0.8669 1.2917  1.857 0.9845
                crtBa 0.003654 0.005649 0.6469 1.0833  1.429 0.9953
                crtE  0.001577 0.003231 0.4881 2.0833  2.143 1.0000
                
                Contrast: saltwater_soil 
                
                average       sd  ratio    ava  avb cumsum
                crtOa 0.0509971 0.030107 1.6939 8.0000 14.6 0.1470
                crtP  0.0359999 0.038151 0.9436 7.2500 12.0 0.2507
                crtN  0.0355600 0.023230 1.5308 7.0833 11.6 0.3532
                crtIa 0.0355355 0.020881 1.7018 3.7500  8.0 0.4556
                crtQa 0.0289835 0.016270 1.7814 3.0000  7.0 0.5391
                crtIb 0.0240755 0.014439 1.6674 2.9583  6.0 0.6084
                crtQb 0.0237570 0.012222 1.9437 3.0000  6.2 0.6769
                cruF  0.0237570 0.012222 1.9437 3.0000  6.2 0.7454
                crtBb 0.0196028 0.011313 1.7328 3.0000  4.8 0.8018
                crtX  0.0159825 0.009772 1.6355 0.3333  2.4 0.8479
                cruA  0.0125976 0.007979 1.5789 0.6250  2.2 0.8842
                cruP  0.0123768 0.007454 1.6604 0.6250  2.2 0.9199
                crtL  0.0082336 0.006543 1.2583 1.0833  0.0 0.9436
                crtW  0.0067817 0.009407 0.7209 0.2083  1.0 0.9631
                crtG  0.0030201 0.003597 0.8396 0.1250  0.4 0.9718
                crtRa 0.0029310 0.005116 0.5729 1.0833  1.4 0.9803
                crtRb 0.0026097 0.004847 0.5384 1.2500  1.2 0.9878
                crtM  0.0020140 0.005544 0.3633 1.2917  1.0 0.9936
                crtE  0.0017437 0.003079 0.5664 2.0833  2.2 0.9986
                crtBa 0.0004757 0.002307 0.2062 1.0833  1.0 1.0000
                
                Contrast: saltwater_rock 
                
                average       sd  ratio    ava avb cumsum
                crtOa 0.0831233 0.067196 1.2370 8.0000 4.0 0.1342
                crtP  0.0730738 0.068864 1.0611 7.2500 3.5 0.2521
                crtN  0.0674548 0.053509 1.2606 7.0833 5.0 0.3610
                crtIa 0.0652559 0.021404 3.0488 3.7500 6.0 0.4663
                crtIb 0.0467702 0.015038 3.1102 2.9583 4.0 0.5417
                crtQa 0.0427491 0.017917 2.3860 3.0000 3.5 0.6107
                crtQb 0.0427357 0.016382 2.6086 3.0000 3.5 0.6797
                cruF  0.0427357 0.016382 2.6086 3.0000 3.5 0.7487
                crtBb 0.0341236 0.022240 1.5343 3.0000 0.5 0.8038
                crtM  0.0205283 0.021579 0.9513 1.2917 1.5 0.8369
                crtRa 0.0183746 0.008088 2.2718 1.0833 1.5 0.8665
                crtL  0.0161355 0.016203 0.9958 1.0833 0.0 0.8926
                crtRb 0.0160117 0.012085 1.3249 1.2500 1.0 0.9184
                cruA  0.0133329 0.012244 1.0889 0.6250 1.5 0.9399
                crtE  0.0111678 0.012172 0.9175 2.0833 1.5 0.9580
                cruP  0.0098352 0.010614 0.9266 0.6250 1.0 0.9738
                crtW  0.0087751 0.008172 1.0738 0.2083 1.0 0.9880
                crtX  0.0056571 0.007522 0.7521 0.3333 0.5 0.9971
                crtG  0.0010671 0.003869 0.2758 0.1250 0.0 0.9989
                crtBa 0.0007118 0.003625 0.1964 1.0833 1.0 1.0000
                
                Contrast: saltwater_hot spring 
                
                average       sd  ratio    ava avb cumsum
                crtOa 0.068508 0.038145 1.7960 8.0000 9.8 0.1705
                crtN  0.053791 0.035392 1.5199 7.0833 7.6 0.3044
                crtP  0.051242 0.039001 1.3139 7.2500 7.0 0.4320
                crtQa 0.033927 0.022263 1.5239 3.0000 7.0 0.5164
                crtBb 0.027105 0.015855 1.7096 3.0000 6.0 0.5839
                crtQb 0.027105 0.015855 1.7096 3.0000 6.0 0.6514
                cruF  0.027105 0.015855 1.7096 3.0000 6.0 0.7188
                crtIb 0.020219 0.012368 1.6348 2.9583 5.0 0.7692
                crtIa 0.019085 0.013655 1.3976 3.7500 4.6 0.8167
                crtX  0.011279 0.010370 1.0877 0.3333 1.6 0.8448
                cruA  0.011222 0.009379 1.1965 0.6250 1.6 0.8727
                cruP  0.011222 0.009379 1.1965 0.6250 1.6 0.9006
                crtL  0.010927 0.010244 1.0667 1.0833 0.0 0.9278
                crtM  0.010062 0.011640 0.8644 1.2917 1.6 0.9529
                crtRb 0.005274 0.009393 0.5615 1.2500 0.8 0.9660
                crtW  0.004090 0.005928 0.6899 0.2083 0.4 0.9762
                crtE  0.003847 0.007148 0.5382 2.0833 1.8 0.9858
                crtG  0.002106 0.003822 0.5510 0.1250 0.2 0.9910
                crtRa 0.001891 0.003508 0.5391 1.0833 1.2 0.9957
                crtBa 0.001725 0.003515 0.4907 1.0833 1.2 1.0000
                
                Contrast: saltwater_temperate 
                
                average       sd  ratio    ava avb cumsum
                crtOa 0.0362490 0.028342 1.2790 8.0000 9.5 0.1328
                crtIa 0.0303865 0.020765 1.4634 3.7500 6.5 0.2441
                crtN  0.0280663 0.023410 1.1989 7.0833 9.5 0.3469
                crtP  0.0178560 0.023530 0.7589 7.2500 8.5 0.4123
                cruA  0.0164163 0.011806 1.3905 0.6250 2.0 0.4725
                crtIb 0.0162838 0.009876 1.6489 2.9583 4.5 0.5321
                crtRb 0.0162285 0.015897 1.0208 1.2500 3.0 0.5915
                crtW  0.0120135 0.007145 1.6814 0.2083 1.5 0.6356
                crtBb 0.0119329 0.009403 1.2691 3.0000 1.5 0.6793
                cruP  0.0118127 0.007696 1.5349 0.6250 1.5 0.7225
                crtQa 0.0111809 0.009612 1.1633 3.0000 3.0 0.7635
                crtQb 0.0108174 0.008724 1.2400 3.0000 3.0 0.8031
                cruF  0.0108174 0.008724 1.2400 3.0000 3.0 0.8427
                crtBa 0.0097621 0.010151 0.9617 1.0833 2.0 0.8785
                crtM  0.0085372 0.008697 0.9817 1.2917 2.0 0.9098
                crtL  0.0082550 0.007364 1.1210 1.0833 0.5 0.9400
                crtX  0.0058682 0.006067 0.9673 0.3333 0.5 0.9615
                crtE  0.0054697 0.005327 1.0268 2.0833 1.5 0.9816
                crtG  0.0043205 0.004322 0.9998 0.1250 0.5 0.9974
                crtRa 0.0007152 0.002420 0.2956 1.0833 1.0 1.0000
                
                Contrast: symbiont_freshwater 
                
                average       sd  ratio  ava    avb cumsum
                crtOa 0.058452 0.047815 1.2225 7.50 8.4643 0.1293
                crtN  0.055741 0.041835 1.3324 7.50 8.0714 0.2526
                crtP  0.045336 0.039154 1.1579 6.25 7.2143 0.3529
                crtIa 0.034963 0.031186 1.1211 3.75 5.4643 0.4302
                crtQa 0.031943 0.028103 1.1367 3.75 4.8571 0.5008
                crtIb 0.028467 0.027510 1.0348 2.75 4.7143 0.5638
                crtQb 0.027504 0.027580 0.9972 3.00 4.6429 0.6246
                cruF  0.027504 0.027580 0.9972 3.00 4.6429 0.6855
                crtBb 0.026728 0.025655 1.0418 2.50 4.3571 0.7446
                crtRb 0.024135 0.022722 1.0622 2.75 2.7857 0.7980
                crtM  0.016554 0.016422 1.0080 1.25 2.5000 0.8346
                crtX  0.014864 0.010441 1.4237 1.75 1.2143 0.8675
                crtE  0.011412 0.013779 0.8282 1.25 2.1429 0.8927
                crtG  0.009746 0.008202 1.1883 1.00 0.8571 0.9143
                cruA  0.009389 0.012406 0.7568 1.25 1.8571 0.9350
                cruP  0.009115 0.012311 0.7404 1.25 1.8214 0.9552
                crtW  0.007451 0.009568 0.7787 0.75 0.5714 0.9717
                crtBa 0.006020 0.010328 0.5829 0.75 1.2143 0.9850
                crtRa 0.005761 0.009364 0.6152 0.75 1.2143 0.9977
                crtL  0.001019 0.003180 0.3206 0.00 0.1071 1.0000
                
                Contrast: symbiont_cave 
                
                average       sd  ratio  ava avb cumsum
                crtN  0.058564 0.049094 1.1929 7.50 9.4 0.1325
                crtOa 0.057364 0.045333 1.2654 7.50 8.9 0.2623
                crtP  0.052548 0.051325 1.0238 6.25 9.2 0.3812
                crtIa 0.033718 0.028828 1.1696 3.75 5.3 0.4575
                crtQa 0.031447 0.025985 1.2102 3.75 4.9 0.5286
                crtQb 0.024984 0.022040 1.1336 3.00 4.1 0.5851
                cruF  0.024984 0.022040 1.1336 3.00 4.1 0.6417
                crtIb 0.023187 0.021748 1.0662 2.75 3.9 0.6941
                crtX  0.020237 0.016520 1.2250 1.75 2.9 0.7399
                crtRb 0.017682 0.013398 1.3197 2.75 1.4 0.7799
                crtBb 0.016770 0.016285 1.0298 2.50 2.8 0.8178
                crtM  0.015697 0.014666 1.0704 1.25 2.5 0.8534
                crtE  0.012601 0.014319 0.8801 1.25 2.3 0.8819
                cruA  0.011052 0.014403 0.7674 1.25 2.2 0.9069
                crtW  0.009755 0.011530 0.8460 0.75 1.1 0.9289
                cruP  0.008889 0.012148 0.7317 1.25 1.9 0.9491
                crtBa 0.008523 0.011834 0.7202 0.75 1.5 0.9683
                crtG  0.007961 0.008446 0.9425 1.00 0.4 0.9863
                crtRa 0.006035 0.008089 0.7461 0.75 1.0 1.0000
                crtL  0.000000 0.000000    NaN 0.00 0.0 1.0000
                
                Contrast: symbiont_desert 
                
                average       sd  ratio  ava     avb cumsum
                crtN  0.062621 0.057251 1.0938 7.50 10.6667 0.1361
                crtOa 0.061909 0.052681 1.1752 7.50  9.9259 0.2706
                crtP  0.046265 0.041471 1.1156 6.25  7.7778 0.3711
                crtIa 0.036031 0.032215 1.1185 3.75  5.8148 0.4494
                crtQa 0.034075 0.030879 1.1035 3.75  5.2593 0.5234
                crtIb 0.028882 0.028353 1.0186 2.75  4.7037 0.5862
                crtQb 0.027863 0.026229 1.0623 3.00  4.6667 0.6467
                cruF  0.027863 0.026229 1.0623 3.00  4.6667 0.7073
                crtX  0.025201 0.026538 0.9496 1.75  3.1481 0.7620
                crtRb 0.020474 0.016529 1.2387 2.75  2.2222 0.8065
                crtBb 0.015515 0.011199 1.3854 2.50  1.6296 0.8402
                crtM  0.013163 0.013484 0.9762 1.25  2.0741 0.8688
                crtE  0.011436 0.013995 0.8172 1.25  2.1852 0.8937
                cruA  0.009291 0.011029 0.8424 1.25  1.6667 0.9139
                crtW  0.009000 0.009495 0.9478 0.75  1.1111 0.9334
                cruP  0.008667 0.010653 0.8136 1.25  1.5926 0.9523
                crtG  0.007212 0.007596 0.9494 1.00  0.2963 0.9679
                crtBa 0.006889 0.009276 0.7427 0.75  1.3333 0.9829
                crtRa 0.006831 0.012128 0.5632 0.75  1.2963 0.9978
                crtL  0.001033 0.003212 0.3217 0.00  0.1111 1.0000
                
                Contrast: symbiont_tropical 
                
                average       sd  ratio  ava    avb cumsum
                crtOa 0.054600 0.045766 1.1930 7.50 10.429 0.1197
                crtN  0.053448 0.044739 1.1947 7.50 10.000 0.2369
                crtP  0.046391 0.046778 0.9917 6.25  9.286 0.3386
                crtIa 0.044766 0.039363 1.1373 3.75  8.143 0.4368
                crtQa 0.035272 0.029103 1.2120 3.75  6.571 0.5142
                crtQb 0.030806 0.029654 1.0388 3.00  6.143 0.5817
                cruF  0.030806 0.029654 1.0388 3.00  6.143 0.6493
                crtIb 0.026454 0.023556 1.1230 2.75  5.143 0.7073
                crtX  0.021456 0.018271 1.1743 1.75  3.429 0.7543
                crtRb 0.019905 0.015815 1.2586 2.75  2.857 0.7980
                crtBb 0.013777 0.010332 1.3334 2.50  2.000 0.8282
                cruA  0.012057 0.015516 0.7770 1.25  2.429 0.8546
                crtM  0.010077 0.010736 0.9386 1.25  1.857 0.8767
                crtRa 0.010017 0.011543 0.8678 0.75  2.000 0.8987
                cruP  0.009934 0.012176 0.8159 1.25  2.286 0.9205
                crtE  0.009873 0.012303 0.8025 1.25  2.143 0.9421
                crtW  0.009849 0.010058 0.9792 0.75  1.429 0.9637
                crtG  0.009744 0.009376 1.0392 1.00  1.286 0.9851
                crtBa 0.006810 0.009626 0.7074 0.75  1.429 1.0000
                crtL  0.000000 0.000000    NaN 0.00  0.000 1.0000
                
                Contrast: symbiont_soil 
                
                average       sd  ratio  ava  avb cumsum
                crtOa 0.069884 0.063480 1.1009 7.50 14.6 0.1519
                crtP  0.058888 0.053662 1.0974 6.25 12.0 0.2798
                crtN  0.054220 0.051761 1.0475 7.50 11.6 0.3977
                crtIa 0.042120 0.035382 1.1904 3.75  8.0 0.4892
                crtQa 0.034643 0.028283 1.2249 3.75  7.0 0.5645
                crtIb 0.030116 0.024701 1.2192 2.75  6.0 0.6300
                crtQb 0.029518 0.025393 1.1625 3.00  6.2 0.6941
                cruF  0.029518 0.025393 1.1625 3.00  6.2 0.7583
                crtBb 0.026144 0.022612 1.1562 2.50  4.8 0.8151
                crtX  0.015630 0.012294 1.2713 1.75  2.4 0.8490
                crtRb 0.014351 0.012135 1.1825 2.75  1.2 0.8802
                crtE  0.009359 0.011026 0.8489 1.25  2.2 0.9006
                cruA  0.008715 0.010968 0.7946 1.25  2.2 0.9195
                cruP  0.008444 0.010446 0.8084 1.25  2.2 0.9379
                crtW  0.008258 0.009829 0.8402 0.75  1.0 0.9558
                crtG  0.006272 0.006344 0.9887 1.00  0.4 0.9694
                crtM  0.005655 0.006015 0.9401 1.25  1.0 0.9817
                crtRa 0.005532 0.007898 0.7004 0.75  1.4 0.9937
                crtBa 0.002881 0.005320 0.5416 0.75  1.0 1.0000
                crtL  0.000000 0.000000    NaN 0.00  0.0 1.0000
                
                Contrast: symbiont_rock 
                
                average       sd  ratio  ava avb cumsum
                crtN  0.108987 0.100831 1.0809 7.50 5.0 0.1542
                crtOa 0.072570 0.058114 1.2488 7.50 4.0 0.2568
                crtP  0.061050 0.050535 1.2081 6.25 3.5 0.3432
                crtIa 0.060823 0.043845 1.3872 3.75 6.0 0.4293
                crtE  0.048386 0.115415 0.4192 1.25 1.5 0.4977
                crtBa 0.043137 0.117329 0.3677 0.75 1.0 0.5587
                crtIb 0.042350 0.029032 1.4587 2.75 4.0 0.6187
                crtQa 0.042282 0.033507 1.2619 3.75 3.5 0.6785
                crtQb 0.039002 0.028472 1.3699 3.00 3.5 0.7336
                cruF  0.039002 0.028472 1.3699 3.00 3.5 0.7888
                crtRb 0.025226 0.023311 1.0821 2.75 1.0 0.8245
                crtBb 0.024584 0.020886 1.1770 2.50 0.5 0.8593
                cruA  0.018628 0.018407 1.0120 1.25 1.5 0.8856
                crtM  0.016834 0.013704 1.2284 1.25 1.5 0.9095
                crtX  0.015676 0.015453 1.0144 1.75 0.5 0.9316
                crtRa 0.015595 0.011083 1.4071 0.75 1.5 0.9537
                cruP  0.014661 0.019013 0.7711 1.25 1.0 0.9744
                crtW  0.009835 0.009835 1.0000 0.75 1.0 0.9883
                crtG  0.008236 0.011488 0.7169 1.00 0.0 1.0000
                crtL  0.000000 0.000000    NaN 0.00 0.0 1.0000
                
                Contrast: symbiont_hot spring 
                
                average       sd  ratio  ava avb cumsum
                crtOa 0.070926 0.054560 1.3000 7.50 9.8 0.1353
                crtN  0.058652 0.046846 1.2520 7.50 7.6 0.2471
                crtP  0.053351 0.040522 1.3166 6.25 7.0 0.3489
                crtQa 0.047910 0.039722 1.2061 3.75 7.0 0.4402
                crtBb 0.040880 0.036128 1.1315 2.50 6.0 0.5182
                crtQb 0.040252 0.036432 1.1049 3.00 6.0 0.5950
                cruF  0.040252 0.036432 1.1049 3.00 6.0 0.6717
                crtIa 0.034393 0.034584 0.9945 3.75 4.6 0.7373
                crtIb 0.032880 0.035646 0.9224 2.75 5.0 0.8000
                crtRb 0.017987 0.017466 1.0298 2.75 0.8 0.8343
                crtX  0.015015 0.013559 1.1074 1.75 1.6 0.8630
                crtE  0.011640 0.017787 0.6544 1.25 1.8 0.8852
                crtM  0.011596 0.010895 1.0644 1.25 1.6 0.9073
                cruA  0.010241 0.014457 0.7084 1.25 1.6 0.9268
                cruP  0.010241 0.014457 0.7084 1.25 1.6 0.9464
                crtRa 0.007512 0.015676 0.4792 0.75 1.2 0.9607
                crtBa 0.007385 0.015598 0.4735 0.75 1.2 0.9748
                crtG  0.007020 0.008388 0.8369 1.00 0.2 0.9882
                crtW  0.006211 0.006960 0.8923 0.75 0.4 1.0000
                crtL  0.000000 0.000000    NaN 0.00 0.0 1.0000
                
                Contrast: symbiont_temperate 
                
                average       sd  ratio  ava avb cumsum
                crtOa 0.063785 0.056010 1.1388 7.50 9.5 0.1404
                crtN  0.062886 0.056079 1.1214 7.50 9.5 0.2789
                crtP  0.051511 0.053657 0.9600 6.25 8.5 0.3923
                crtIa 0.043132 0.039177 1.1009 3.75 6.5 0.4873
                crtIb 0.025432 0.029309 0.8677 2.75 4.5 0.5433
                crtRb 0.025102 0.021329 1.1769 2.75 3.0 0.5986
                crtQa 0.025095 0.017795 1.4102 3.75 3.0 0.6538
                crtQb 0.020263 0.017517 1.1568 3.00 3.0 0.6984
                cruF  0.020263 0.017517 1.1568 3.00 3.0 0.7430
                crtBb 0.015422 0.011520 1.3386 2.50 1.5 0.7770
                crtBa 0.014960 0.018474 0.8098 0.75 2.0 0.8099
                cruA  0.014556 0.017177 0.8474 1.25 2.0 0.8420
                crtX  0.013448 0.009228 1.4574 1.75 0.5 0.8716
                crtM  0.012800 0.013458 0.9512 1.25 2.0 0.8998
                crtW  0.011930 0.012565 0.9495 0.75 1.5 0.9261
                cruP  0.009056 0.012300 0.7363 1.25 1.5 0.9460
                crtE  0.008500 0.009558 0.8892 1.25 1.5 0.9647
                crtG  0.007761 0.007533 1.0303 1.00 0.5 0.9818
                crtL  0.004307 0.005105 0.8437 0.00 0.5 0.9913
                crtRa 0.003960 0.007470 0.5300 0.75 1.0 1.0000
                
                Contrast: freshwater_cave 
                
                average       sd  ratio    ava avb  cumsum
                crtOa 0.0215391 0.016363 1.3163 8.4643 8.9 0.09814
                crtN  0.0198386 0.015628 1.2694 8.0714 9.4 0.18852
                crtP  0.0181472 0.015167 1.1965 7.2143 9.2 0.27120
                crtIa 0.0169787 0.013742 1.2355 5.4643 5.3 0.34856
                crtQa 0.0166798 0.013601 1.2263 4.8571 4.9 0.42456
                crtQb 0.0160249 0.011020 1.4541 4.6429 4.1 0.49757
                cruF  0.0160249 0.011020 1.4541 4.6429 4.1 0.57058
                crtIb 0.0157367 0.011341 1.3876 4.7143 3.9 0.64228
                crtBb 0.0152360 0.009840 1.5484 4.3571 2.8 0.71170
                crtRb 0.0132181 0.014212 0.9301 2.7857 1.4 0.77192
                crtX  0.0131183 0.007850 1.6712 1.2143 2.9 0.83169
                crtW  0.0085000 0.008102 1.0492 0.5714 1.1 0.87042
                crtG  0.0060221 0.005036 1.1957 0.8571 0.4 0.89785
                crtM  0.0055790 0.006592 0.8463 2.5000 2.5 0.92327
                crtBa 0.0045106 0.006180 0.7299 1.2143 1.5 0.94382
                crtRa 0.0040398 0.004969 0.8130 1.2143 1.0 0.96223
                cruA  0.0030483 0.004788 0.6366 1.8571 2.2 0.97612
                crtE  0.0025736 0.003524 0.7304 2.1429 2.3 0.98784
                cruP  0.0018735 0.003761 0.4982 1.8214 1.9 0.99638
                crtL  0.0007945 0.002309 0.3441 0.1071 0.0 1.00000
                
                Contrast: freshwater_desert 
                
                average       sd  ratio    ava     avb cumsum
                crtOa 0.030067 0.023142 1.2992 8.4643  9.9259 0.1171
                crtN  0.028385 0.021126 1.3436 8.0714 10.6667 0.2277
                crtBb 0.020462 0.012979 1.5765 4.3571  1.6296 0.3074
                crtQa 0.019079 0.014877 1.2825 4.8571  5.2593 0.3817
                crtIb 0.018708 0.014479 1.2921 4.7143  4.7037 0.4545
                crtIa 0.017298 0.014231 1.2155 5.4643  5.8148 0.5219
                crtQb 0.017190 0.012948 1.3277 4.6429  4.6667 0.5889
                cruF  0.017190 0.012948 1.3277 4.6429  4.6667 0.6558
                crtP  0.016741 0.014615 1.1454 7.2143  7.7778 0.7210
                crtX  0.016635 0.017027 0.9769 1.2143  3.1481 0.7858
                crtRb 0.014297 0.013022 1.0980 2.7857  2.2222 0.8415
                crtW  0.008017 0.006842 1.1717 0.5714  1.1111 0.8727
                crtM  0.007461 0.007906 0.9437 2.5000  2.0741 0.9018
                crtG  0.005783 0.005011 1.1540 0.8571  0.2963 0.9243
                cruA  0.004564 0.005864 0.7783 1.8571  1.6667 0.9421
                crtBa 0.003918 0.005596 0.7000 1.2143  1.3333 0.9573
                cruP  0.003863 0.005730 0.6742 1.8214  1.5926 0.9724
                crtRa 0.003702 0.007407 0.4998 1.2143  1.2963 0.9868
                crtE  0.001958 0.003539 0.5534 2.1429  2.1852 0.9944
                crtL  0.001426 0.002946 0.4842 0.1071  0.1111 1.0000
                
                Contrast: freshwater_tropical 
                
                average       sd  ratio    ava    avb  cumsum
                crtOa 0.0253884 0.017506 1.4502 8.4643 10.429 0.09968
                crtIa 0.0211082 0.015333 1.3767 5.4643  8.143 0.18255
                crtQa 0.0211081 0.016826 1.2545 4.8571  6.571 0.26542
                crtN  0.0207664 0.015971 1.3002 8.0714 10.000 0.34695
                crtP  0.0205407 0.014351 1.4313 7.2143  9.286 0.42759
                crtQb 0.0190307 0.016582 1.1477 4.6429  6.143 0.50231
                cruF  0.0190307 0.016582 1.1477 4.6429  6.143 0.57702
                crtBb 0.0172942 0.011536 1.4991 4.3571  2.000 0.64492
                crtIb 0.0172187 0.011952 1.4407 4.7143  5.143 0.71252
                crtX  0.0157396 0.011105 1.4174 1.2143  3.429 0.77432
                crtRb 0.0141507 0.012141 1.1655 2.7857  2.857 0.82987
                crtW  0.0090892 0.006859 1.3252 0.5714  1.429 0.86556
                crtM  0.0074635 0.006683 1.1169 2.5000  1.857 0.89486
                crtG  0.0068550 0.005757 1.1907 0.8571  1.286 0.92177
                crtRa 0.0062384 0.007334 0.8506 1.2143  2.000 0.94627
                cruA  0.0046012 0.006487 0.7092 1.8571  2.429 0.96433
                crtBa 0.0036854 0.005141 0.7169 1.2143  1.429 0.97880
                cruP  0.0029853 0.004317 0.6915 1.8214  2.286 0.99052
                crtE  0.0016892 0.003023 0.5587 2.1429  2.143 0.99715
                crtL  0.0007255 0.002125 0.3413 0.1071  0.000 1.00000
                
                Contrast: freshwater_soil 
                
                average       sd  ratio    ava  avb cumsum
                crtOa 0.0408807 0.023529 1.7375 8.4643 14.6 0.1586
                crtP  0.0308832 0.031823 0.9705 7.2143 12.0 0.2783
                crtIa 0.0255606 0.014528 1.7594 5.4643  8.0 0.3775
                crtN  0.0249698 0.017494 1.4274 8.0714 11.6 0.4743
                crtQa 0.0189509 0.014255 1.3294 4.8571  7.0 0.5478
                crtIb 0.0172761 0.012362 1.3975 4.7143  6.0 0.6148
                crtQb 0.0149505 0.010670 1.4012 4.6429  6.2 0.6728
                cruF  0.0149505 0.010670 1.4012 4.6429  6.2 0.7308
                crtBb 0.0148340 0.010162 1.4598 4.3571  4.8 0.7883
                crtRb 0.0109463 0.013161 0.8317 2.7857  1.2 0.8308
                crtX  0.0101044 0.008052 1.2549 1.2143  2.4 0.8700
                crtM  0.0095947 0.005327 1.8012 2.5000  1.0 0.9072
                crtW  0.0073204 0.009303 0.7869 0.5714  1.0 0.9356
                crtG  0.0046566 0.004386 1.0617 0.8571  0.4 0.9536
                crtRa 0.0031961 0.005052 0.6327 1.2143  1.4 0.9660
                cruA  0.0025720 0.004032 0.6379 1.8571  2.2 0.9760
                cruP  0.0022683 0.003739 0.6067 1.8214  2.2 0.9848
                crtE  0.0018059 0.002897 0.6234 2.1429  2.2 0.9918
                crtBa 0.0014265 0.003784 0.3770 1.2143  1.0 0.9973
                crtL  0.0006905 0.002020 0.3418 0.1071  0.0 1.0000
                
                Contrast: freshwater_rock 
                
                average       sd  ratio    ava avb cumsum
                crtOa 0.067601 0.058426 1.1570 8.4643 4.0 0.1171
                crtIa 0.059644 0.024882 2.3970 5.4643 6.0 0.2204
                crtN  0.055842 0.047602 1.1731 8.0714 5.0 0.3171
                crtP  0.055373 0.048536 1.1409 7.2143 3.5 0.4130
                crtIb 0.043446 0.027723 1.5671 4.7143 4.0 0.4882
                crtQa 0.041861 0.030354 1.3791 4.8571 3.5 0.5607
                crtQb 0.040568 0.028343 1.4313 4.6429 3.5 0.6309
                cruF  0.040568 0.028343 1.4313 4.6429 3.5 0.7012
                crtBb 0.040390 0.025571 1.5795 4.3571 0.5 0.7711
                crtRb 0.024701 0.025940 0.9522 2.7857 1.0 0.8139
                crtM  0.019896 0.020113 0.9892 2.5000 1.5 0.8484
                cruA  0.017190 0.011898 1.4448 1.8571 1.5 0.8781
                crtRa 0.014844 0.006862 2.1632 1.2143 1.5 0.9038
                cruP  0.013817 0.014406 0.9591 1.8214 1.0 0.9278
                crtX  0.010762 0.012204 0.8819 1.2143 0.5 0.9464
                crtG  0.009343 0.008831 1.0580 0.8571 0.0 0.9626
                crtW  0.009336 0.010407 0.8971 0.5714 1.0 0.9787
                crtE  0.008678 0.008854 0.9801 2.1429 1.5 0.9938
                crtBa 0.002467 0.007325 0.3368 1.2143 1.0 0.9980
                crtL  0.001128 0.003533 0.3192 0.1071 0.0 1.0000
                
                Contrast: freshwater_hot spring 
                
                average       sd  ratio    ava avb cumsum
                crtOa 0.0583896 0.034016 1.7166 8.4643 9.8 0.1691
                crtN  0.0453119 0.033339 1.3591 8.0714 7.6 0.3004
                crtP  0.0421817 0.028375 1.4866 7.2143 7.0 0.4225
                crtQa 0.0279561 0.018614 1.5019 4.8571 7.0 0.5035
                crtBb 0.0219891 0.014431 1.5237 4.3571 6.0 0.5672
                crtQb 0.0218503 0.015148 1.4424 4.6429 6.0 0.6305
                cruF  0.0218503 0.015148 1.4424 4.6429 6.0 0.6938
                crtIb 0.0173135 0.014554 1.1896 4.7143 5.0 0.7439
                crtIa 0.0166816 0.017389 0.9593 5.4643 4.6 0.7922
                crtRb 0.0160821 0.019201 0.8375 2.7857 0.8 0.8388
                crtM  0.0119808 0.012629 0.9487 2.5000 1.6 0.8735
                crtX  0.0110244 0.008875 1.2422 1.2143 1.6 0.9054
                crtG  0.0064561 0.006310 1.0231 0.8571 0.2 0.9241
                cruA  0.0057148 0.009609 0.5947 1.8571 1.6 0.9407
                crtW  0.0055493 0.008690 0.6386 0.5714 0.4 0.9568
                cruP  0.0054771 0.009521 0.5753 1.8214 1.6 0.9726
                crtE  0.0035370 0.005961 0.5933 2.1429 1.8 0.9829
                crtBa 0.0025832 0.005098 0.5067 1.2143 1.2 0.9904
                crtRa 0.0024752 0.004649 0.5324 1.2143 1.2 0.9975
                crtL  0.0008532 0.002611 0.3268 0.1071 0.0 1.0000
                
                Contrast: freshwater_temperate 
                
                average       sd  ratio    ava avb cumsum
                crtOa 0.031065 0.020527 1.5134 8.4643 9.5 0.1273
                crtBb 0.020654 0.012936 1.5966 4.3571 1.5 0.2119
                crtIa 0.020346 0.014478 1.4053 5.4643 6.5 0.2952
                crtN  0.018173 0.014531 1.2506 8.0714 9.5 0.3697
                crtQa 0.017381 0.014119 1.2310 4.8571 3.0 0.4409
                crtRb 0.017348 0.013852 1.2524 2.7857 3.0 0.5119
                crtQb 0.016086 0.012616 1.2751 4.6429 3.0 0.5778
                cruF  0.016086 0.012616 1.2751 4.6429 3.0 0.6437
                crtP  0.013122 0.012980 1.0109 7.2143 8.5 0.6975
                crtIb 0.012588 0.008395 1.4994 4.7143 4.5 0.7491
                crtW  0.010704 0.006573 1.6286 0.5714 1.5 0.7929
                crtBa 0.008113 0.008001 1.0140 1.2143 2.0 0.8262
                cruA  0.007971 0.004007 1.9896 1.8571 2.0 0.8588
                crtM  0.007946 0.007626 1.0419 2.5000 2.0 0.8914
                crtX  0.007081 0.006497 1.0900 1.2143 0.5 0.9204
                crtG  0.005288 0.005138 1.0292 0.8571 0.5 0.9420
                crtE  0.005170 0.004987 1.0367 2.1429 1.5 0.9632
                cruP  0.003855 0.003993 0.9656 1.8214 1.5 0.9790
                crtL  0.003559 0.003639 0.9780 0.1071 0.5 0.9936
                crtRa 0.001563 0.004240 0.3687 1.2143 1.0 1.0000
                
                Contrast: cave_desert 
                
                average       sd  ratio ava     avb cumsum
                crtN  0.0241303 0.017156 1.4065 9.4 10.6667 0.1079
                crtOa 0.0226029 0.019243 1.1746 8.9  9.9259 0.2090
                crtIa 0.0188671 0.014416 1.3087 5.3  5.8148 0.2933
                crtQa 0.0181460 0.013749 1.3198 4.9  5.2593 0.3745
                crtP  0.0180696 0.016402 1.1016 9.2  7.7778 0.4553
                crtIb 0.0170814 0.014952 1.1424 3.9  4.7037 0.5317
                crtQb 0.0158799 0.012221 1.2993 4.1  4.6667 0.6027
                cruF  0.0158799 0.012221 1.2993 4.1  4.6667 0.6737
                crtX  0.0144819 0.012290 1.1784 2.9  3.1481 0.7384
                crtBb 0.0098881 0.007588 1.3031 2.8  1.6296 0.7826
                crtRb 0.0097822 0.008606 1.1367 1.4  2.2222 0.8264
                crtM  0.0073781 0.007837 0.9415 2.5  2.0741 0.8594
                crtW  0.0053306 0.006220 0.8570 1.1  1.1111 0.8832
                crtBa 0.0051829 0.006215 0.8339 1.5  1.3333 0.9064
                cruA  0.0050736 0.006592 0.7697 2.2  1.6667 0.9291
                crtRa 0.0048767 0.007245 0.6731 1.0  1.2963 0.9509
                crtG  0.0039119 0.004993 0.7835 0.4  0.2963 0.9684
                cruP  0.0035075 0.005526 0.6347 1.9  1.5926 0.9840
                crtE  0.0027617 0.003789 0.7289 2.3  2.1852 0.9964
                crtL  0.0008059 0.002308 0.3491 0.0  0.1111 1.0000
                
                Contrast: cave_tropical 
                
                average       sd  ratio ava    avb cumsum
                crtIa 0.023351 0.015884 1.4701 5.3  8.143 0.1052
                crtQa 0.020429 0.015926 1.2828 4.9  6.571 0.1973
                crtQb 0.019827 0.015463 1.2822 4.1  6.143 0.2867
                cruF  0.019827 0.015463 1.2822 4.1  6.143 0.3761
                crtOa 0.017817 0.011663 1.5276 8.9 10.429 0.4564
                crtN  0.017137 0.012262 1.3975 9.4 10.000 0.5336
                crtP  0.016378 0.012675 1.2922 9.2  9.286 0.6074
                crtIb 0.016222 0.012511 1.2967 3.9  5.143 0.6805
                crtRb 0.011349 0.009605 1.1816 1.4  2.857 0.7317
                crtX  0.011234 0.007690 1.4607 2.9  3.429 0.7823
                crtBb 0.007521 0.006570 1.1447 2.8  2.000 0.8162
                crtG  0.007414 0.006420 1.1549 0.4  1.286 0.8496
                crtRa 0.007337 0.007301 1.0049 1.0  2.000 0.8827
                crtM  0.007318 0.006666 1.0978 2.5  1.857 0.9157
                crtW  0.005507 0.006044 0.9111 1.1  1.429 0.9405
                crtBa 0.004685 0.005458 0.8583 1.5  1.429 0.9616
                cruA  0.003723 0.005088 0.7318 2.2  2.429 0.9784
                cruP  0.002405 0.003500 0.6873 1.9  2.286 0.9892
                crtE  0.002387 0.003303 0.7228 2.3  2.143 1.0000
                crtL  0.000000 0.000000    NaN 0.0  0.000 1.0000
                
                Contrast: cave_soil 
                
                average       sd  ratio ava  avb cumsum
                crtOa 0.034961 0.017572 1.9896 8.9 14.6 0.1456
                crtIa 0.027000 0.015008 1.7991 5.3  8.0 0.2581
                crtP  0.026309 0.027282 0.9643 9.2 12.0 0.3677
                crtIb 0.018670 0.012084 1.5450 3.9  6.0 0.4455
                crtQa 0.018237 0.013336 1.3676 4.9  7.0 0.5215
                crtN  0.017465 0.013722 1.2728 9.4 11.6 0.5942
                crtBb 0.016791 0.008143 2.0621 2.8  4.8 0.6642
                crtQb 0.016273 0.009745 1.6699 4.1  6.2 0.7320
                cruF  0.016273 0.009745 1.6699 4.1  6.2 0.7998
                crtM  0.009285 0.005266 1.7631 2.5  1.0 0.8385
                crtX  0.009149 0.004985 1.8351 2.9  2.4 0.8766
                crtW  0.008494 0.007319 1.1605 1.1  1.0 0.9120
                crtRa 0.004358 0.004852 0.8982 1.0  1.4 0.9301
                crtRb 0.003824 0.003854 0.9921 1.4  1.2 0.9460
                crtG  0.003531 0.004071 0.8675 0.4  0.4 0.9608
                crtBa 0.003238 0.005251 0.6167 1.5  1.0 0.9743
                crtE  0.002386 0.003134 0.7612 2.3  2.2 0.9842
                cruA  0.002070 0.003092 0.6693 2.2  2.2 0.9928
                cruP  0.001725 0.002909 0.5930 1.9  2.2 1.0000
                crtL  0.000000 0.000000    NaN 0.0  0.0 1.0000
                
                Contrast: cave_rock 
                
                average       sd  ratio ava avb cumsum
                crtP  0.070500 0.060197 1.1711 9.2 3.5 0.1232
                crtOa 0.065830 0.060134 1.0947 8.9 4.0 0.2381
                crtN  0.062923 0.055897 1.1257 9.4 5.0 0.3481
                crtIa 0.057339 0.023648 2.4247 5.3 6.0 0.4482
                crtQa 0.041379 0.027851 1.4857 4.9 3.5 0.5205
                crtIb 0.039697 0.019363 2.0501 3.9 4.0 0.5898
                crtQb 0.037644 0.022074 1.7053 4.1 3.5 0.6556
                cruF  0.037644 0.022074 1.7053 4.1 3.5 0.7214
                crtX  0.025924 0.015482 1.6745 2.9 0.5 0.7666
                crtBb 0.025244 0.016557 1.5247 2.8 0.5 0.8107
                crtM  0.018876 0.018057 1.0454 2.5 1.5 0.8437
                cruA  0.018328 0.015381 1.1916 2.2 1.5 0.8757
                cruP  0.013789 0.014332 0.9622 1.9 1.0 0.8998
                crtRa 0.013671 0.007643 1.7887 1.0 1.5 0.9237
                crtRb 0.011986 0.009738 1.2308 1.4 1.0 0.9446
                crtW  0.011858 0.012617 0.9399 1.1 1.0 0.9653
                crtE  0.010090 0.010037 1.0053 2.3 1.5 0.9830
                crtBa 0.005344 0.009388 0.5693 1.5 1.0 0.9923
                crtG  0.004403 0.008251 0.5337 0.4 0.0 1.0000
                crtL  0.000000 0.000000    NaN 0.0 0.0 1.0000
                
                Contrast: cave_hot spring 
                
                average       sd  ratio ava avb cumsum
                crtOa 0.057701 0.028831 2.0013 8.9 9.8 0.1655
                crtP  0.046398 0.038265 1.2126 9.2 7.0 0.2986
                crtN  0.046044 0.040418 1.1392 9.4 7.6 0.4307
                crtQa 0.027184 0.017158 1.5843 4.9 7.0 0.5087
                crtBb 0.023172 0.011991 1.9324 2.8 6.0 0.5752
                crtQb 0.021613 0.013084 1.6519 4.1 6.0 0.6372
                cruF  0.021613 0.013084 1.6519 4.1 6.0 0.6992
                crtIa 0.017253 0.016626 1.0377 5.3 4.6 0.7487
                crtIb 0.016886 0.011260 1.4997 3.9 5.0 0.7972
                crtX  0.015639 0.013135 1.1907 2.9 1.6 0.8420
                crtM  0.011554 0.011782 0.9807 2.5 1.6 0.8752
                crtW  0.008886 0.008623 1.0306 1.1 0.4 0.9007
                cruA  0.006504 0.011140 0.5838 2.2 1.6 0.9193
                crtRb 0.006229 0.007069 0.8811 1.4 0.8 0.9372
                cruP  0.005129 0.009464 0.5419 1.9 1.6 0.9519
                crtE  0.004701 0.006718 0.6997 2.3 1.8 0.9654
                crtBa 0.004436 0.006559 0.6763 1.5 1.2 0.9781
                crtRa 0.003858 0.004459 0.8652 1.0 1.2 0.9892
                crtG  0.003768 0.005733 0.6573 0.4 0.2 1.0000
                crtL  0.000000 0.000000    NaN 0.0 0.0 1.0000
                
                Contrast: cave_temperate 
                
                average       sd  ratio ava avb cumsum
                crtOa 0.025888 0.011621 2.2277 8.9 9.5 0.1182
                crtIa 0.022643 0.014313 1.5820 5.3 6.5 0.2216
                crtX  0.017380 0.005332 3.2592 2.9 0.5 0.3010
                crtQa 0.016583 0.013068 1.2689 4.9 3.0 0.3767
                crtRb 0.014837 0.012141 1.2220 1.4 3.0 0.4444
                crtQb 0.012940 0.010333 1.2523 4.1 3.0 0.5035
                cruF  0.012940 0.010333 1.2523 4.1 3.0 0.5626
                crtIb 0.012536 0.006611 1.8962 3.9 4.5 0.6199
                crtN  0.012085 0.009241 1.3078 9.4 9.5 0.6750
                crtP  0.010847 0.008706 1.2459 9.2 8.5 0.7246
                crtBb 0.010172 0.005207 1.9535 2.8 1.5 0.7710
                crtBa 0.007839 0.007630 1.0274 1.5 2.0 0.8068
                crtM  0.007799 0.007686 1.0148 2.5 2.0 0.8424
                cruA  0.007399 0.003473 2.1302 2.2 2.0 0.8762
                crtW  0.006345 0.005694 1.1143 1.1 1.5 0.9052
                crtE  0.006291 0.005630 1.1173 2.3 1.5 0.9339
                crtG  0.004454 0.004727 0.9422 0.4 0.5 0.9543
                cruP  0.003519 0.003647 0.9649 1.9 1.5 0.9703
                crtL  0.003457 0.003573 0.9675 0.0 0.5 0.9861
                crtRa 0.003042 0.003846 0.7908 1.0 1.0 1.0000
                
                Contrast: desert_tropical 
                
                average       sd  ratio     ava    avb  cumsum
                crtOa 0.0238725 0.019604 1.2178  9.9259 10.429 0.09698
                crtN  0.0232158 0.016944 1.3702 10.6667 10.000 0.19128
                crtQa 0.0218365 0.016097 1.3566  5.2593  6.571 0.27999
                crtIa 0.0205316 0.016532 1.2420  5.8148  8.143 0.36340
                crtQb 0.0201092 0.015919 1.2632  4.6667  6.143 0.44508
                cruF  0.0201092 0.015919 1.2632  4.6667  6.143 0.52677
                crtP  0.0198931 0.016060 1.2387  7.7778  9.286 0.60758
                crtIb 0.0185207 0.014919 1.2414  4.7037  5.143 0.68282
                crtX  0.0162294 0.013237 1.2261  3.1481  3.429 0.74875
                crtRb 0.0117605 0.009703 1.2120  2.2222  2.857 0.79652
                crtG  0.0072170 0.006507 1.1091  0.2963  1.286 0.82584
                crtM  0.0070138 0.006988 1.0037  2.0741  1.857 0.85433
                crtRa 0.0068817 0.008224 0.8368  1.2963  2.000 0.88229
                crtBb 0.0066037 0.005904 1.1184  1.6296  2.000 0.90911
                cruA  0.0063355 0.007794 0.8128  1.6667  2.429 0.93485
                cruP  0.0046681 0.005959 0.7834  1.5926  2.286 0.95381
                crtW  0.0044065 0.005507 0.8002  1.1111  1.429 0.97171
                crtBa 0.0043466 0.005238 0.8298  1.3333  1.429 0.98937
                crtE  0.0018803 0.003382 0.5560  2.1852  2.143 0.99701
                crtL  0.0007366 0.002125 0.3466  0.1111  0.000 1.00000
                
                Contrast: desert_soil 
                
                average       sd  ratio     ava  avb cumsum
                crtOa 0.0340684 0.024778 1.3749  9.9259 14.6 0.1315
                crtP  0.0303446 0.030941 0.9807  7.7778 12.0 0.2486
                crtIa 0.0245038 0.016564 1.4794  5.8148  8.0 0.3432
                crtBb 0.0214829 0.011903 1.8049  1.6296  4.8 0.4261
                crtN  0.0212307 0.016099 1.3187 10.6667 11.6 0.5081
                crtIb 0.0197858 0.013568 1.4583  4.7037  6.0 0.5844
                crtQa 0.0190533 0.013672 1.3936  5.2593  7.0 0.6580
                crtQb 0.0159685 0.010586 1.5084  4.6667  6.2 0.7196
                cruF  0.0159685 0.010586 1.5084  4.6667  6.2 0.7812
                crtX  0.0132954 0.013147 1.0113  3.1481  2.4 0.8325
                crtW  0.0081539 0.006314 1.2915  1.1111  1.0 0.8640
                crtRb 0.0078055 0.008077 0.9664  2.2222  1.2 0.8941
                crtM  0.0070356 0.006068 1.1594  2.0741  1.0 0.9213
                cruA  0.0043131 0.005522 0.7810  1.6667  2.2 0.9380
                crtRa 0.0040219 0.006664 0.6035  1.2963  1.4 0.9535
                cruP  0.0038696 0.005305 0.7295  1.5926  2.2 0.9684
                crtG  0.0029650 0.003361 0.8823  0.2963  0.4 0.9799
                crtBa 0.0025309 0.004161 0.6083  1.3333  1.0 0.9896
                crtE  0.0019871 0.003219 0.6173  2.1852  2.2 0.9973
                crtL  0.0007015 0.002020 0.3473  0.1111  0.0 1.0000
                
                Contrast: desert_rock 
                
                average       sd  ratio     ava avb cumsum
                crtOa 0.075219 0.062399 1.2054  9.9259 4.0 0.1295
                crtN  0.073404 0.058786 1.2487 10.6667 5.0 0.2558
                crtIa 0.058807 0.028247 2.0819  5.8148 6.0 0.3570
                crtP  0.058334 0.049362 1.1818  7.7778 3.5 0.4574
                crtQa 0.043657 0.034395 1.2693  5.2593 3.5 0.5325
                crtIb 0.043256 0.028569 1.5141  4.7037 4.0 0.6070
                crtQb 0.039935 0.027410 1.4570  4.6667 3.5 0.6757
                cruF  0.039935 0.027410 1.4570  4.6667 3.5 0.7444
                crtX  0.027969 0.030993 0.9024  3.1481 0.5 0.7926
                crtRb 0.019026 0.017122 1.1112  2.2222 1.0 0.8253
                crtM  0.016888 0.015589 1.0833  2.0741 1.5 0.8544
                cruA  0.015518 0.010962 1.4157  1.6667 1.5 0.8811
                crtRa 0.015387 0.010694 1.4389  1.2963 1.5 0.9076
                crtBb 0.013415 0.013673 0.9811  1.6296 0.5 0.9307
                cruP  0.012065 0.012152 0.9928  1.5926 1.0 0.9514
                crtW  0.011173 0.010040 1.1129  1.1111 1.0 0.9707
                crtE  0.008767 0.009121 0.9612  2.1852 1.5 0.9857
                crtBa 0.004154 0.007445 0.5580  1.3333 1.0 0.9929
                crtG  0.002989 0.005762 0.5188  0.2963 0.0 0.9980
                crtL  0.001143 0.003576 0.3196  0.1111 0.0 1.0000
                
                Contrast: desert_hot spring 
                
                average       sd  ratio     ava avb cumsum
                crtOa 0.058317 0.039802 1.4652  9.9259 9.8 0.1579
                crtN  0.050907 0.045827 1.1108 10.6667 7.6 0.2958
                crtP  0.042623 0.030728 1.3871  7.7778 7.0 0.4112
                crtBb 0.029113 0.014087 2.0666  1.6296 6.0 0.4900
                crtQa 0.027922 0.019423 1.4376  5.2593 7.0 0.5656
                crtQb 0.021979 0.015407 1.4265  4.6667 6.0 0.6251
                cruF  0.021979 0.015407 1.4265  4.6667 6.0 0.6846
                crtIb 0.019939 0.016329 1.2211  4.7037 5.0 0.7386
                crtX  0.019660 0.021124 0.9307  3.1481 1.6 0.7919
                crtIa 0.019090 0.018347 1.0405  5.8148 4.6 0.8436
                crtRb 0.011621 0.012929 0.8988  2.2222 0.8 0.8750
                crtM  0.010359 0.010725 0.9659  2.0741 1.6 0.9031
                crtW  0.008528 0.006576 1.2967  1.1111 0.4 0.9262
                cruA  0.006585 0.008875 0.7419  1.6667 1.6 0.9440
                cruP  0.006051 0.008632 0.7010  1.5926 1.6 0.9604
                crtE  0.003702 0.006207 0.5965  2.1852 1.8 0.9704
                crtBa 0.003663 0.005304 0.6906  1.3333 1.2 0.9803
                crtRa 0.003529 0.007649 0.4613  1.2963 1.2 0.9899
                crtG  0.002872 0.004285 0.6702  0.2963 0.2 0.9977
                crtL  0.000865 0.002619 0.3303  0.1111 0.0 1.0000
                
                Contrast: desert_temperate 
                
                average       sd  ratio     ava avb cumsum
                crtOa 0.031059 0.023034 1.3484  9.9259 9.5 0.1256
                crtN  0.022815 0.014632 1.5592 10.6667 9.5 0.2178
                crtIa 0.021699 0.014713 1.4748  5.8148 6.5 0.3055
                crtQa 0.019454 0.015787 1.2323  5.2593 3.0 0.3842
                crtX  0.018879 0.018248 1.0346  3.1481 0.5 0.4605
                crtIb 0.016573 0.011035 1.5018  4.7037 4.5 0.5275
                crtQb 0.016223 0.013031 1.2449  4.6667 3.0 0.5931
                cruF  0.016223 0.013031 1.2449  4.6667 3.0 0.6587
                crtRb 0.014638 0.012058 1.2140  2.2222 3.0 0.7179
                crtP  0.013913 0.014887 0.9346  7.7778 8.5 0.7741
                cruA  0.008610 0.006608 1.3030  1.6667 2.0 0.8090
                crtBa 0.008287 0.008170 1.0143  1.3333 2.0 0.8425
                crtM  0.007870 0.007976 0.9867  2.0741 2.0 0.8743
                crtBb 0.006083 0.005341 1.1389  1.6296 1.5 0.8989
                crtE  0.005312 0.005264 1.0090  2.1852 1.5 0.9203
                crtW  0.005085 0.004996 1.0180  1.1111 1.5 0.9409
                cruP  0.004601 0.005076 0.9063  1.5926 1.5 0.9595
                crtG  0.003830 0.003922 0.9766  0.2963 0.5 0.9750
                crtL  0.003496 0.003598 0.9715  0.1111 0.5 0.9891
                crtRa 0.002690 0.007225 0.3723  1.2963 1.0 1.0000
                
                Contrast: tropical_soil 
                
                average       sd  ratio    ava  avb cumsum
                crtOa 0.027564 0.019816 1.3910 10.429 14.6 0.1174
                crtP  0.027333 0.025687 1.0641  9.286 12.0 0.2339
                crtQa 0.019401 0.013638 1.4225  6.571  7.0 0.3165
                crtBb 0.018681 0.010505 1.7784  2.000  4.8 0.3961
                crtIa 0.018183 0.016327 1.1136  8.143  8.0 0.4735
                crtN  0.017668 0.012709 1.3903 10.000 11.6 0.5488
                crtIb 0.017587 0.012317 1.4279  5.143  6.0 0.6237
                crtQb 0.016233 0.011987 1.3542  6.143  6.2 0.6929
                cruF  0.016233 0.011987 1.3542  6.143  6.2 0.7620
                crtX  0.011490 0.008830 1.3013  3.429  2.4 0.8110
                crtRb 0.009647 0.009261 1.0416  2.857  1.2 0.8521
                crtW  0.008672 0.006019 1.4408  1.429  1.0 0.8890
                crtG  0.006114 0.005855 1.0442  1.286  0.4 0.9150
                crtRa 0.005822 0.006605 0.8814  2.000  1.4 0.9398
                crtM  0.004617 0.005577 0.8279  1.857  1.0 0.9595
                cruA  0.003229 0.004522 0.7141  2.429  2.2 0.9733
                crtBa 0.002543 0.004248 0.5986  1.429  1.0 0.9841
                cruP  0.002014 0.002728 0.7382  2.286  2.2 0.9927
                crtE  0.001717 0.002797 0.6140  2.143  2.2 1.0000
                crtL  0.000000 0.000000    NaN  0.000  0.0 1.0000
                
                Contrast: tropical_rock 
                
                average       sd  ratio    ava avb cumsum
                crtOa 0.067847 0.054669 1.2411 10.429 4.0 0.1192
                crtP  0.063441 0.052623 1.2056  9.286 3.5 0.2306
                crtIa 0.059847 0.042029 1.4240  8.143 6.0 0.3357
                crtN  0.056820 0.052112 1.0903 10.000 5.0 0.4355
                crtQa 0.043619 0.032317 1.3497  6.571 3.5 0.5121
                crtQb 0.041953 0.030111 1.3933  6.143 3.5 0.5858
                cruF  0.041953 0.030111 1.3933  6.143 3.5 0.6595
                crtIb 0.038657 0.023857 1.6204  5.143 4.0 0.7274
                crtX  0.025496 0.019957 1.2775  3.429 0.5 0.7722
                crtRb 0.020041 0.016449 1.2184  2.857 1.0 0.8074
                cruA  0.017930 0.017618 1.0177  2.429 1.5 0.8389
                crtRa 0.015517 0.010792 1.4378  2.000 1.5 0.8662
                cruP  0.014480 0.014631 0.9897  2.286 1.0 0.8916
                crtBb 0.014330 0.012192 1.1753  2.000 0.5 0.9168
                crtM  0.014192 0.011569 1.2267  1.857 1.5 0.9417
                crtW  0.011549 0.010830 1.0664  1.429 1.0 0.9620
                crtG  0.010089 0.010197 0.9895  1.286 0.0 0.9797
                crtE  0.007563 0.008504 0.8894  2.143 1.5 0.9930
                crtBa 0.003999 0.007157 0.5588  1.429 1.0 1.0000
                crtL  0.000000 0.000000    NaN  0.000 0.0 1.0000
                
                Contrast: tropical_hot spring 
                
                average       sd  ratio    ava avb cumsum
                crtOa 0.052598 0.034166 1.5395 10.429 9.8 0.1432
                crtN  0.043329 0.037168 1.1658 10.000 7.6 0.2612
                crtP  0.042417 0.035316 1.2011  9.286 7.0 0.3766
                crtQa 0.027769 0.019560 1.4196  6.571 7.0 0.4522
                crtIa 0.027597 0.023015 1.1991  8.143 4.6 0.5274
                crtBb 0.025100 0.013270 1.8916  2.000 6.0 0.5957
                crtQb 0.023180 0.017312 1.3390  6.143 6.0 0.6588
                cruF  0.023180 0.017312 1.3390  6.143 6.0 0.7219
                crtIb 0.018013 0.013698 1.3150  5.143 5.0 0.7710
                crtX  0.017910 0.014598 1.2268  3.429 1.6 0.8197
                crtRb 0.013592 0.013173 1.0318  2.857 0.8 0.8567
                crtW  0.009690 0.007269 1.3331  1.429 0.4 0.8831
                crtM  0.008490 0.008648 0.9818  1.857 1.6 0.9062
                cruA  0.007739 0.012099 0.6397  2.429 1.6 0.9273
                crtG  0.007718 0.007695 1.0029  1.286 0.2 0.9483
                crtRa 0.006210 0.007812 0.7949  2.000 1.2 0.9652
                cruP  0.006050 0.009708 0.6232  2.286 1.6 0.9817
                crtBa 0.003499 0.005148 0.6798  1.429 1.2 0.9912
                crtE  0.003231 0.005636 0.5733  2.143 1.8 1.0000
                crtL  0.000000 0.000000    NaN  0.000 0.0 1.0000
                
                Contrast: tropical_temperate 
                
                average       sd  ratio    ava avb cumsum
                crtOa 0.025489 0.017454 1.4604 10.429 9.5 0.1041
                crtQa 0.023670 0.018817 1.2579  6.571 3.0 0.2008
                crtQb 0.021055 0.018562 1.1343  6.143 3.0 0.2868
                cruF  0.021055 0.018562 1.1343  6.143 3.0 0.3728
                crtIa 0.020610 0.014939 1.3796  8.143 6.5 0.4570
                crtX  0.018392 0.012037 1.5279  3.429 0.5 0.5321
                crtRb 0.014798 0.010978 1.3480  2.857 3.0 0.5925
                crtP  0.014642 0.011690 1.2525  9.286 8.5 0.6524
                crtIb 0.014350 0.009926 1.4457  5.143 4.5 0.7110
                crtN  0.012906 0.012113 1.0654 10.000 9.5 0.7637
                cruA  0.007754 0.004933 1.5717  2.429 2.0 0.7953
                crtBa 0.007168 0.006867 1.0439  1.429 2.0 0.8246
                crtG  0.006893 0.006775 1.0173  1.286 0.5 0.8528
                crtM  0.006723 0.007103 0.9464  1.857 2.0 0.8802
                crtBb 0.006003 0.004377 1.3713  2.000 1.5 0.9048
                crtRa 0.005898 0.007949 0.7420  2.000 1.0 0.9288
                cruP  0.004960 0.004362 1.1371  2.286 1.5 0.9491
                crtE  0.004782 0.004983 0.9597  2.143 1.5 0.9686
                crtW  0.004505 0.004358 1.0336  1.429 1.5 0.9870
                crtL  0.003173 0.003355 0.9456  0.000 0.5 1.0000
                
                Contrast: soil_rock 
                
                average       sd  ratio  ava avb cumsum
                crtOa 0.096645 0.065606 1.4731 14.6 4.0 0.1680
                crtP  0.074587 0.062933 1.1852 12.0 3.5 0.2977
                crtN  0.064479 0.055944 1.1526 11.6 5.0 0.4098
                crtIa 0.053349 0.041243 1.2935  8.0 6.0 0.5026
                crtQa 0.041968 0.033544 1.2511  7.0 3.5 0.5755
                crtIb 0.038001 0.029414 1.2920  6.0 4.0 0.6416
                crtBb 0.037406 0.024856 1.5049  4.8 0.5 0.7066
                crtQb 0.036579 0.031251 1.1705  6.2 3.5 0.7702
                cruF  0.036579 0.031251 1.1705  6.2 3.5 0.8338
                crtX  0.016526 0.014421 1.1459  2.4 0.5 0.8625
                cruA  0.014501 0.011694 1.2400  2.2 1.5 0.8877
                cruP  0.012465 0.012645 0.9858  2.2 1.0 0.9094
                crtRa 0.012066 0.006720 1.7956  1.4 1.5 0.9304
                crtM  0.011300 0.002195 5.1477  1.0 1.5 0.9500
                crtW  0.009693 0.010667 0.9086  1.0 1.0 0.9669
                crtRb 0.008812 0.005368 1.6418  1.2 1.0 0.9822
                crtE  0.007221 0.007589 0.9515  2.2 1.5 0.9948
                crtG  0.003009 0.004183 0.7194  0.4 0.0 1.0000
                crtBa 0.000000 0.000000    NaN  1.0 1.0 1.0000
                crtL  0.000000 0.000000    NaN  0.0 0.0 1.0000
                
                Contrast: soil_hot spring 
                
                average       sd  ratio  ava avb cumsum
                crtOa 0.0551789 0.054821 1.0065 14.6 9.8 0.1696
                crtP  0.0516907 0.044494 1.1617 12.0 7.0 0.3284
                crtN  0.0435085 0.043896 0.9912 11.6 7.6 0.4621
                crtIa 0.0299964 0.020382 1.4717  8.0 4.6 0.5543
                crtQa 0.0244066 0.017366 1.4054  7.0 7.0 0.6293
                crtBb 0.0190346 0.012548 1.5169  4.8 6.0 0.6878
                crtIb 0.0183375 0.015166 1.2091  6.0 5.0 0.7441
                crtQb 0.0172706 0.013893 1.2432  6.2 6.0 0.7972
                cruF  0.0172706 0.013893 1.2432  6.2 6.0 0.8502
                crtX  0.0123405 0.010528 1.1722  2.4 1.6 0.8882
                crtW  0.0067962 0.009186 0.7399  1.0 0.4 0.9090
                crtM  0.0061805 0.005307 1.1647  1.0 1.6 0.9280
                cruA  0.0052329 0.008908 0.5874  2.2 1.6 0.9441
                cruP  0.0050267 0.008547 0.5881  2.2 1.6 0.9596
                crtE  0.0032820 0.005197 0.6315  2.2 1.8 0.9697
                crtRb 0.0031244 0.004853 0.6439  1.2 0.8 0.9793
                crtRa 0.0030226 0.004774 0.6332  1.4 1.2 0.9885
                crtG  0.0027367 0.003318 0.8248  0.4 0.2 0.9969
                crtBa 0.0009926 0.002039 0.4869  1.0 1.2 1.0000
                crtL  0.0000000 0.000000    NaN  0.0 0.0 1.0000
                
                Contrast: soil_temperate 
                
                average       sd  ratio  ava avb cumsum
                crtOa 0.037283 0.024308 1.5338 14.6 9.5 0.1357
                crtP  0.024862 0.031276 0.7949 12.0 8.5 0.2261
                crtQa 0.024752 0.015336 1.6140  7.0 3.0 0.3162
                crtIa 0.023692 0.018583 1.2750  8.0 6.5 0.4024
                crtBb 0.022488 0.011209 2.0064  4.8 1.5 0.4842
                crtQb 0.020291 0.010898 1.8619  6.2 3.0 0.5580
                cruF  0.020291 0.010898 1.8619  6.2 3.0 0.6319
                crtIb 0.016072 0.008819 1.8225  6.0 4.5 0.6903
                crtN  0.014127 0.011893 1.1879 11.6 9.5 0.7417
                crtRb 0.012192 0.012139 1.0044  1.2 3.0 0.7861
                crtX  0.011768 0.009104 1.2927  2.4 0.5 0.8289
                crtW  0.010011 0.005685 1.7610  1.0 1.5 0.8653
                crtBa 0.007010 0.007539 0.9298  1.0 2.0 0.8909
                cruA  0.006441 0.003045 2.1153  2.2 2.0 0.9143
                crtM  0.006059 0.006482 0.9347  1.0 2.0 0.9363
                crtE  0.004810 0.004747 1.0134  2.2 1.5 0.9538
                cruP  0.004145 0.003826 1.0832  2.2 1.5 0.9689
                crtG  0.003191 0.003404 0.9374  0.4 0.5 0.9805
                crtL  0.003030 0.003241 0.9347  0.0 0.5 0.9916
                crtRa 0.002322 0.004907 0.4731  1.4 1.0 1.0000
                
                Contrast: rock_hot spring 
                
                average       sd  ratio ava avb cumsum
                crtOa 0.075364 0.056737 1.3283 4.0 9.8 0.1201
                crtIa 0.063678 0.028465 2.2371 6.0 4.6 0.2216
                crtQa 0.060243 0.042345 1.4227 3.5 7.0 0.3177
                crtN  0.059096 0.040637 1.4542 5.0 7.6 0.4119
                crtBb 0.058543 0.033548 1.7450 0.5 6.0 0.5052
                crtP  0.053616 0.043101 1.2440 3.5 7.0 0.5907
                crtQb 0.052788 0.040471 1.3043 3.5 6.0 0.6748
                cruF  0.052788 0.040471 1.3043 3.5 6.0 0.7590
                crtIb 0.050698 0.035883 1.4129 4.0 5.0 0.8398
                crtRa 0.018107 0.014328 1.2637 1.5 1.2 0.8686
                cruA  0.015951 0.014263 1.1183 1.5 1.6 0.8941
                crtM  0.014862 0.012676 1.1725 1.5 1.6 0.9178
                cruP  0.012525 0.015471 0.8096 1.0 1.6 0.9377
                crtX  0.012347 0.011912 1.0366 0.5 1.6 0.9574
                crtRb 0.009688 0.007027 1.3787 1.0 0.8 0.9729
                crtW  0.007725 0.008685 0.8895 1.0 0.4 0.9852
                crtE  0.006262 0.007735 0.8096 1.5 1.8 0.9952
                crtG  0.001662 0.003705 0.4485 0.0 0.2 0.9978
                crtBa 0.001377 0.003026 0.4551 1.0 1.2 1.0000
                crtL  0.000000 0.000000    NaN 0.0 0.0 1.0000
                
                Contrast: rock_temperate 
                
                average       sd  ratio ava avb cumsum
                crtOa 0.078816 0.069481 1.1344 4.0 9.5 0.1361
                crtP  0.068473 0.068627 0.9978 3.5 8.5 0.2543
                crtN  0.065122 0.072462 0.8987 5.0 9.5 0.3667
                crtIa 0.065072 0.039267 1.6572 6.0 6.5 0.4791
                crtIb 0.045219 0.025341 1.7844 4.0 4.5 0.5572
                crtQa 0.035081 0.013079 2.6823 3.5 3.0 0.6177
                crtQb 0.035081 0.013079 2.6823 3.5 3.0 0.6783
                cruF  0.035081 0.013079 2.6823 3.5 3.0 0.7389
                crtRb 0.026666 0.024955 1.0685 1.0 3.0 0.7849
                cruA  0.019428 0.022912 0.8480 1.5 2.0 0.8185
                crtM  0.017477 0.015654 1.1164 1.5 2.0 0.8486
                crtRa 0.014253 0.002362 6.0334 1.5 1.0 0.8732
                crtBb 0.013735 0.015150 0.9066 0.5 1.5 0.8970
                cruP  0.013489 0.015298 0.8817 1.0 1.5 0.9202
                crtW  0.013489 0.015298 0.8817 1.0 1.5 0.9435
                crtBa 0.012369 0.016595 0.7453 1.0 2.0 0.9649
                crtX  0.005939 0.008272 0.7179 0.5 0.5 0.9751
                crtE  0.004963 0.006150 0.8070 1.5 1.5 0.9837
                crtL  0.004717 0.006050 0.7798 0.0 0.5 0.9919
                crtG  0.004717 0.006050 0.7798 0.0 0.5 1.0000
                
                Contrast: hot spring_temperate 
                
                average       sd  ratio ava avb cumsum
                crtOa 0.059844 0.042020 1.4242 9.8 9.5 0.1561
                crtN  0.047318 0.045855 1.0319 7.6 9.5 0.2796
                crtP  0.046507 0.037204 1.2500 7.0 8.5 0.4009
                crtBb 0.029420 0.014832 1.9836 6.0 1.5 0.4776
                crtQa 0.028818 0.020776 1.3871 7.0 3.0 0.5528
                crtIa 0.023847 0.022185 1.0749 4.6 6.5 0.6150
                crtQb 0.022804 0.014648 1.5569 6.0 3.0 0.6745
                cruF  0.022804 0.014648 1.5569 6.0 3.0 0.7340
                crtRb 0.017230 0.018295 0.9418 0.8 3.0 0.7789
                crtIb 0.013641 0.011683 1.1676 5.0 4.5 0.8145
                crtW  0.011300 0.008495 1.3302 0.4 1.5 0.8440
                cruA  0.011037 0.011576 0.9535 1.6 2.0 0.8728
                crtX  0.010337 0.007504 1.3775 1.6 0.5 0.8998
                crtM  0.010216 0.010655 0.9588 1.6 2.0 0.9264
                crtBa 0.008866 0.010208 0.8686 1.2 2.0 0.9495
                cruP  0.006569 0.008828 0.7441 1.6 1.5 0.9667
                crtE  0.004118 0.004622 0.8910 1.8 1.5 0.9774
                crtG  0.003764 0.004259 0.8836 0.2 0.5 0.9872
                crtL  0.003670 0.004190 0.8759 0.0 0.5 0.9968
                crtRa 0.001225 0.002590 0.4730 1.2 1.0 1.0000
                Permutation: free
                Number of permutations: 0