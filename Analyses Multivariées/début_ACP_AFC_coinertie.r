library(ade4)
library(vegan)
library(ggplot2)
library(corrplot)
library(tidyr)
library(dplyr)
library(RVAideMemoire)
setwd("C:/Users/serva/Google Drive/1-Partage Ordis/M2/ASA/Projet")

assoc <- read.table("association.txt", header = TRUE)
coord <- read.table("coord.txt", header = TRUE)
flo <- read.table("flo.txt",header = TRUE)
mil <- read.table("mil.txt",header = TRUE)

############## ACP sur les variabels environnementales #########
# Matrice des corrÃ©lations
mcor <- cor(mil)
corrplot(mcor, method = "number", type = "upper", order = "hclust", tl.col = "black")

data_ACP <- mil[,-10]
View(data_ACP)

# RÃ©alisation de l'ACP
acp_env <- dudi.pca(data_ACP, scannf = FALSE, nf = 2)
acp_env$eig # Valeurs propres brutes 
acp_env$eig*100/sum(acp_env$eig) # Valeurs propres en pourcentage
s.label(acp_env$li, boxes = FALSE) # Affichage des stations dans le nouveau plan formÃ© par l'acp
s.corcircle(acp_env$co) ## Affichage des variables environnementales dans le nouveau plan de l'acp

iner_acp_env <- inertia.dudi(acp_env, col.inertia = TRUE, row.inertia = TRUE)
iner_acp_env$col.abs # RÃ©cupÃ©ration des contributions absolues pour les colonnes (ie les variables environnementales)

seuil = 1/length(iner_acp_env$col.abs[,1])*100
seuil

row.names(iner_acp_env$col.abs)[which(iner_acp_env$col.abs[,1] > seuil)] # RÃ©cupÃ©ration des variables qui contribuent Ã  l'axe 1
row.names(iner_acp_env$col.abs)[which(iner_acp_env$col.abs[,2] > seuil)] # RÃ©cupÃ©ration des variables qui contribuent Ã  l'axe 2


# RÃ©alisation de la CAH sur l'ACP
dist_env <- dist.quant(acp_env$li,1)
dendo <- hclust(dist_env,"ward.D2")
plot(dendo$height, type="s")
abline(h = 6.7)

plot(dendo)
abline(h = 6.7)

decoup <- cutree(dendo, 6)
decoup

s.class(acp_env$li,as.factor(decoup), clabel = 0.5) # Remarque : pour le cluster 6, seulement prÃ©sence de deux stations donc pas bien visible sous l'Ã©tiquette. Changer le clabel Ã  0 pour pouvoir voir les points.

# Affichage des Ã©tiquettes qui correspondent 
indice_axe1 <- which(iner_acp_env$col.abs[,1] > seuil)
s.label(acp_env$co[indice_axe1,], boxes = FALSE, add.plot = TRUE)

indice_axe2 <- which(iner_acp_env$col.abs[,2] > seuil)
s.label(acp_env$co[indice_axe2,], boxes=TRUE, add.plot = TRUE, clabel = 0.8)

#  decoupe des stations en fonction des associations #
asso <- assoc %>%
  mutate(
    group = dense_rank(association) # crée une colonne si la station "C1" alors donne le rang "1"
  ) %>%
  select(group)
decoup_asso <- (t(asso)) # passe les lignes deviennent colonnes et colonne deviennent lignes
s.label(acp_env$li) # Station with CAP
s.class(acp_env$li, as.factor(decoup_asso), col = 1:7, label = c("","","","","","","") , add.plot = TRUE) # add group of associations
title(main = "CAP with association group") # add title
legend(x= 4.5,y = 4,legend = c("C1","C2","C3","C4","C5","C6","C7"),col = 1:7, pch = 1) # add legend




#################### AFC sur le tableau des occurences d'espÃ¨ces#################### 

afc_flo <- dudi.coa(flo, scannf = FALSE, nf = 2)
summary(afc_flo)

s.label(afc_flo$li)
s.label(afc_flo$co)

iner_afc_flo <- inertia.dudi(afc_flo, col.inertia = TRUE, row.inertia = TRUE)
a = format(iner_afc_flo$col.abs,scientific = FALSE)
seuil_AFC<-1/min(c(nrow(flo),ncol(afc_flo)))*100

which(iner_afc_flo$col.abs[,1] > seuil_AFC)

dist_flo <- dist.quant(afc_flo$li,1)

dendo <- hclust(dist_flo,"ward.D2")

plot(dendo$height, type="s")
abline(h = 1.94)

plot(dendo)
abline(h = 1.94)

decoup <- cutree(dendo, 7)
decoup

s.class(afc_flo$li,as.factor(decoup))

# Affichage des Ã©tiquettes qui correspondent 
indice_axe1 <- which(iner_afc_flo$col.abs[,1] > seuil)
s.label(afc_flo$co[indice_axe1,], boxes = FALSE)

indice_axe2 <- which(iner_afc_flo$col.abs[,2] > seuil)
s.label(afc_flo$co[indice_axe2,], boxes=TRUE, add.plot = TRUE)

#  decoupe des stations en fonction des associations #
asso <- assoc %>%
  mutate(
    group = dense_rank(association) # crée une colonne si la station "C1" alors donne le rang "1"
  ) %>%
  select(group)
decoup_asso <- (t(asso)) # passe les lignes deviennent colonnes et colonne deviennent lignes

s.label(afc_flo$li) # station with AFC
s.class(afc_flo$li,as.factor(decoup_asso), col = 1:7, label = c("","","","","","",""), add.plot=T) # add group of association
title(main = "CA with association group") # add title
legend(x= 1,y = 1,legend = c("C1","C2","C3","C4","C5","C6","C7"),col = 1:7, pch = 1) # add legend

#############coinertia####################
afc_flo #already done before
acp_coinertia_env <- dudi.pca(data_ACP, row.w = afc_flo$lw, scannf = FALSE, nf = 2) #weight the PCA row with the CA row weights
coi <- coinertia(acp_coinertia_env, afc_flo, scannf = FALSE, nf=2) #selection of 2 axis
coi

randtest(coi,nrepet=1000)
plot(randtest(coi), main = "Monte-Carlo test") #the coinertia analysis is significant


iner=inertia.dudi(coi,col.inertia=T,row.inertia=T)
iner

#In particular : 
abscoiE=iner$col.abs #get the absolute contribution of each environmental variable on each axis
abscoiE
selectE=rbind(abscoiE[abscoiE[,1]>(100/10),], abscoiE[abscoiE[,2]>(100/10),]) #variable which contributes the most to the axis (threshold of 10%)
selectE

abscoiS=iner$row.abs #get the absolute contribution of each species on each axis
abscoiS
selectS=rbind(abscoiS[abscoiS[,1]>5,], abscoiS[abscoiS[,2]>5,]) #species which contributes the most to the axis (threshold of 5%)
selectS

plot(coi)

#représentations plus précises (PAS FORCEMENT BIEN FAITES !!!!!!!!!)
#plot the correlation circle of the first table:
MVA.plot(coi,"corr",space=1) 
s.corcircle(coi$c1,box = FALSE, grid = TRUE,  possub = "bottomleft",cgrid = 0, fullcircle = TRUE, add.plot = FALSE)
#plot the correlation circle of the second table:
MVA.plot(coi,"corr",space=2) 
s.arrow(coi$l1,boxes=FALSE,add.plot=F)
#correlation between pairs of axes (of the 2 correlation circles): 2 points clouds lead to 2 correlation circles: 2 first axes more or less correlated and so on
MVA.synt(coi)
#plot of the stations
MVA.plot(coi,space=1) #to see stations from the first table
MVA.plot(coi,space=2) #to see stations from the second table
s.match(coi$mX, coi$mY) #to see changes associated with changing table, permet de faire des classes de stations ! 
fac=cbind(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5)
fac=as.factor(sapply(fac, factor))
s.match.class(coi$mX, coi$mY,fac,col1 = rep("BLUE",nlevels(fac)),label =levels(fac))
