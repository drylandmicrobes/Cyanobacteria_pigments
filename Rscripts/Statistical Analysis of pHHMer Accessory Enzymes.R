---
title: "Statistical Analysis of pHHMer Accessory Enzymes"
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
phmmer_accessory <- read_csv("Cyano Carotenoids pHMMER Results - phmmer_accessory.csv")
#organize data to just raw hits and strain names 
#row.names(phmmer_accessory) = phmmer_accessory$Strain
#phmmer_accessory_matrix <- data.matrix(phmmer_accessory[c(3:16)])
phmmer_accessory_matrix = subset(phmmer_accessory, select=c(crtRa, crtRb, crtG, crtX, crtW, crtOa, crtN,crtP, crtQa, crtQb, cruF, crtM),row.names=1)
range(phmmer_accessory_matrix) #0 23

#analyze pca statistics to find proportion of variatiance explained by each gene 
pca_phmmer_accessory= prcomp(phmmer_accessory_matrix, center= TRUE, scale.= TRUE)
summary(pca_phmmer_accessory)
str(pca_phmmer_accessory)
#


#make groups based on habitat
habitat= phmmer_accessory$Habitat


#try to make group colours, this should be added to plot code 
scale_colour_manual(name= "Habitat", values= c("brown", "dark orange", "blue", "pink", "red", "dark blue", "forest green", "dark yellow", "red", "black"))

#make the PCA plot with strains as points NEED TO DO PCA OF SEPERATE SETS OF GENES 
#png("PCA Plot of Carotenoid BLAST Hit .png", width=1000, height=1000)
ggbiplot(pca_phmmer_accessory,ellipse=TRUE,obs.scale = 3,var.scale = 1, groups=habitat) +
  scale_colour_manual(name="Habitat", values= c("saddlebrown", "darkorange1", "cornflowerblue", "deeppink", "red", "dark blue", "forest green", "gold", "darkmagenta", "black"))+
  ggtitle("PCA of pHMMer Accessory Enzymes")+
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
nmds_phmmer_accessory_habitat= as.data.frame(phmmer_accessory[,1:2])
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
                rownames(phmmer_accessory_matrix) = phmmer_accessory$Strain
                str(phmmer_accessory_matrix)
                nmds_phmmer_accessory_habitat= as.data.frame(phmmer_accessory[,1:2])
                str(nmds_phmmer_accessory_habitat)
                
                #make a distance matrix
                phmmer_accessory_matrix_ds=vegdist(phmmer_accessory_matrix, method = 'euclidian')
                
                #Run metaMDS on distance matrix
                phmmer_accessory_nms=metaMDS(phmmer_accessory_matrix_ds, distance="euclidian",k=2,trymax=500)
                phmmer_accessory_nms
                #stress=0.08
                
                #Run PERMANOVA
                pmv_phmmer_accessory = adonis(phmmer_accessory_matrix ~ habitat, data = nmds_phmmer_small_habitat, permutations = 999, method = 'euclidean')
                pmv_phmmer_accessory
                #PERMANOVA Results: Based off of the carotenoid genes being compared
                #cyanobacterial strains differ statistically significant 
                #between their habitats (p-value:0.001).
                #Habitats explian 19.13% of the variance in the copy numbers of critical 
                #carotenoid biosynthesis enzymes between cyanobacterial species. 
                
                #Visualise the permanova
                densityplot(permustats(pmv_phmmer_small))
                
                #Distance Based Dispersion test 
                bd_phmmer_accessory = betadisper(phmmer_accessory_matrix_ds, phmmer_small$Habitat)
                boxplot(bd_phmmer_accessory)
                anova(bd_phmmer_accessory) 
                #F-Test: pvalue 0.003
                permutest(bd_phmmer_accessory)
                #permutation test: pvalue 0.008
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
                