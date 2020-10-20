
################################### pRDA ####################################
library(ggplot2)
library(tidyr)
library(dplyr)
library(vegan) #only library for pRDA
library(ade4)
library(RVAideMemoire)

assoc <- read.table("association.txt", header = TRUE)
coord <- read.table("coord.txt", header = TRUE)
flo <- read.table("flo.txt", header = TRUE)
mil <- read.table("mil.txt", header = TRUE)

miltri <- mil[, -10]
miltri <- as.data.frame(scale(miltri, center = TRUE, scale = TRUE))


venn <- varpart(flo, miltri, coord) #variance partitioning of the different explanatory variables
venn
plot(venn, bg=2:3, Xnames = c("mil","coord"))

flo <- decostand(flo, "hellinger") #RDA works better

rda_VG_vegan <- rda(flo, miltri)
rda_VG_vegan
anova.cca(rda_VG_vegan) #meme si c'est contre intuitif d'utilser cca c'est normal : on voit que c'est significatif --> on peut continuer
#test colinearity between parameters
test_rda <- vif.cca(rda_VG_vegan) 
test_rda
#we have some values greater than ten -> select some variables
rda_VG_vegan <- rda(flo~., miltri)
ordistep(rda_VG_vegan, perm.max=500)
#we keep : Mg, K+, Capa_Reten, Altitude
selected <- miltri[,c("Mg..", "K.", "Capa_Reten", "Altitude")]
rda_VG_vegan <- rda(flo,selected)
test_rda <- vif.cca(rda_VG_vegan) 
test_rda
#no colinearity . C'est good

#now, let's do it with ade4
pca_VG_ade4 <- dudi.pca(flo, scale = F, scannf = FALSE, nf = 2)
rda_VG_ade4 <- pcaiv(pca_VG_ade4, selected, scannf = FALSE, nf = 2)

#now we have our 2 RDA -> let's analyse !
rda_VG_vegan
rda_VG_ade4
MVA.synt(rda_VG_vegan)



#calculate with ADE4 the percentage of species explained by env:
sum(rda_VG_ade4$eig)/sum(pca_VG_ade4$eig) # See later on 

# 15,72  % of the flo community is explained by the environnement represented by Mg..", "K.", "Capa_Reten", "Altitude (p<0,01)    (21,9 % avant de faire le tri des variables corélées)
# The first axis represents 69,61 % of the constrained inertia (8,87/81,05 =10,9% du total) and is explained by and 


randtest(rda_VG_ade4) # permutation test to see if constrained inertia is significant (i.e, impact of environment on communities is significant)
plot(randtest(rda_VG_ade4)) # graphic representation of the randtest. Result : SIGNIFICANT

# same as vegan : decomposition of the inertia for each axis
inertia.dudi(rda_VG_ade4, col = TRUE, row = FALSE)

plot(rda_VG_ade4, boxes = FALSE)
