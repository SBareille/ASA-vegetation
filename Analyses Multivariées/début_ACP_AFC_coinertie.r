library(ade4)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(tidyr)
library(dplyr)
library(RVAideMemoire)
setwd("C:/Users/serva/Google Drive/1-Partage Ordis/M2/ASA/Projet")

assoc <- read.table("association.txt", header = TRUE)
coord <- read.table("coord.txt", header = TRUE)
flo <- read.table("flo.txt",header = TRUE)
env <- read.table("mil.txt",header = TRUE)

# standardization of data before statistical analysis
flo <- decostand(flo, "hellinger") # Lots of zero
env= as.data.frame(scale(env, center=TRUE, scale=TRUE)) # variables with different units

#################### CA on the Species occurrence table (with stations, and then associations) #################### 

# We have already used the hellinger transformation on the data

CA_flo <- dudi.coa(flo, scannf = FALSE, nf = 2) # we keep two axis
summary(CA_flo)
CA_flo$eig # Gross eigen values of axis
CA_pourcent <- CA_flo$eig*100/sum(CA_flo$eig) # Pourcentages of eigen values of axis
barplot(CA_pourcent, main = "Histogram of Pourcentages of eigen values", ylim=c(0,50), col = heat.colors(10))
s.label(CA_flo$li) # Stations display
s.label(CA_flo$co) # Species display

iner_CA_flo <- inertia.dudi(CA_flo, col.inertia = TRUE, row.inertia = TRUE)
iner_CA_flo$col.abs # Absolute contribution of species on axis

threshold_CA<-1/min(c(nrow(flo),ncol(CA_flo)))*100 # ordinarily used for CA
threshold_CA

row.names(iner_CA_flo$col.abs)[which(iner_CA_flo$col.abs[,1] > threshold_CA)] # Species which contribute the most to axis 1
row.names(iner_CA_flo$col.abs)[which(iner_CA_flo$col.abs[,2] > threshold_CA)] # Species which contribute the most to axis 2

# dist_flo <- dist.quant(CA_flo$li,1)
# 
# dendo <- hclust(dist_flo,"ward.D2")
# 
# plot(dendo$height, type="s")
# abline(h = 1.94)
# 
# plot(dendo)
# abline(h = 1.94)
# 
# decoup <- cutree(dendo, 7)
# decoup
# 
# s.class(CA_flo$li,as.factor(decoup))
# 
# # Affichage des Ã©tiquettes qui correspondent 
# indice_axe1 <- which(iner_CA_flo$col.abs[,1] > threshold)
# s.label(CA_flo$co[indice_axe1,], boxes = FALSE)
# 
# indice_axe2 <- which(iner_CA_flo$col.abs[,2] > threshold)
# s.label(CA_flo$co[indice_axe2,], boxes=TRUE, add.plot = TRUE)

#  we want to see vegetal associations on the graph to study the species repartition on them
asso <- assoc %>%
  mutate(
    group = dense_rank(association) # creates a column : if station is "C1", then ots rank will be "1"
  ) %>%
  select(group)
decoup_asso <- (t(asso)) # ligns become colomns and colomns become ligns

#s.label(CA_flo$li) # station with CA
s.class(CA_flo$li,as.factor(decoup_asso), col=brewer.pal(n=7, name = 'Set2'), label = c("","","","","","",""), add.plot=T) # add group of associations.
title(main = "CA with association groups") # add title
legend(x= 1,y = 1,legend = c("C1","C2","C3","C4","C5","C6","C7"),col=brewer.pal(n=7, name = 'Set2'), pch = 1) # add legend

############## PCA on environmental variables (with stations, and then associations) #########

# We already have centered and reduce the data

# Correlation matrix
mcor <- cor(env)
corrplot(mcor, method = "number", type = "upper", order = "hclust", tl.col = "black")
# Na+/L and Conductivity are highly correlated (0.96). Na+/L is also correlated with NA+/100g (0.71), as we may have guessed 
# We decide to get rid of Na+/L
data_PCA <- env[,-10]
View(data_PCA)

# PCA
PCA_env <- dudi.pca(data_PCA, scannf = FALSE, nf = 2) # we keep 2 axis
PCA_env$eig # Gross eigen values of axis
PCA_pourcent <- PCA_env$eig*100/sum(PCA_env$eig) # Pourcentages of eigen values of axis
barplot(PCA_pourcent, main = "Histogram of Pourcentages of eigen values", ylim=c(0,50), col = heat.colors(10))
s.label(PCA_env$li, boxes = FALSE) # Stations display
s.corcircle(PCA_env$co) # Environmental variables display

iner_PCA_env <- inertia.dudi(PCA_env, col.inertia = TRUE, row.inertia = TRUE)
iner_PCA_env$col.abs # Absolute contribution of environmental variables on axis

threshold_PCA = 1/length(iner_PCA_env$col.abs[,1])*100
threshold_PCA

row.names(iner_PCA_env$col.abs)[which(iner_PCA_env$col.abs[,1] > threshold_PCA)] # Environmental variables which contributes the most to axis 1
row.names(iner_PCA_env$col.abs)[which(iner_PCA_env$col.abs[,2] > threshold_PCA)] # Environmental variables which contributes the most to axis 2


# RÃ©alisation de la CAH sur l'PCA     (j'ai tout mis en commentaire temps qu'on n'a pas décidé ce qu'on en faisait)
# dist_env <- dist.quant(PCA_env$li,1)
# dendo <- hclust(dist_env,"ward.D2")
# plot(dendo$height, type="s")
# abline(h = 6.7)
# 
# plot(dendo)
# abline(h = 6.7)
# 
# decoup <- cutree(dendo, 6)
# decoup
# 
# s.class(PCA_env$li,as.factor(decoup), clabel = 0.5) # Remarque : pour le cluster 6, seulement prÃ©sence de deux stations donc pas bien visible sous l'Ã©tiquette. Changer le clabel Ã  0 pour pouvoir voir les points.
# 
# # Affichage des Ã©tiquettes qui correspondent 
# indice_axe1 <- which(iner_PCA_env$col.abs[,1] > threshold_PCA)
# s.label(PCA_env$co[indice_axe1,], boxes = FALSE, add.plot = TRUE)
# 
# indice_axe2 <- which(iner_PCA_env$col.abs[,2] > threshold_PCA)
# s.label(PCA_env$co[indice_axe2,], boxes=TRUE, add.plot = TRUE, clabel = 0.8)

#  we want to see see vegetal associations on the graph to study the effect of environmental variables on them

s.label(PCA_env$li, ) # station with PCA
s.class(PCA_env$li, as.factor(decoup_asso), col=brewer.pal(n=7, name = 'Set2'), label = c("","","","","","","") , add.plot = TRUE) # add group of associations. decoup_asso has been created for CA
title(main = "PCA with association group") # add title
legend(x= 4.5,y = -0.5,legend = c("C1","C2","C3","C4","C5","C6","C7"),col = brewer.pal(n=7, name = 'Set2'), pch = 1) # add legend



#############  coinertia  ####################
CA_flo #already done before
PCA_coinertia_env <- dudi.pca(data_PCA, row.w = CA_flo$lw, scannf = FALSE, nf = 2) #weight the PCA row with the CA row weights
coi <- coinertia(PCA_coinertia_env, CA_flo, scannf = FALSE, nf=2) #selection of 2 axis
coi # we read the RV coefficient :  0.289

randtest(coi,nrepet=1000)
plot(randtest(coi), main = "Monte-Carlo test") # the coinertia analysis is significant


iner=inertia.dudi(coi,col.inertia=T,row.inertia=T)
iner
# The first axis represents 64.43 % of the coinertia, and the second axis 12.37 %

#In particular : 
abscoiE = iner$col.abs #get the absolute contribution of each environmental variable on each axis
abscoiE
selectE = rbind(abscoiE[abscoiE[,1]>(100/10),], abscoiE[abscoiE[,2]>(100/10),]) #variable which contributes the most to the axis (threshold of 10%)
selectE = selectE[-7,] # we get rid of the seventh line because K+ contributes to axis 1 and 2 and appears two times in the selected varaibles
selectE 
# Six variables contribute to the first axis with a threshold of 10% : K2O, Na+/100g, K+,  Conduc  on the right part, and Altitude on the left part
# Four variables contribute to this second axis with a threshold of 10% :K+ on the top part, and Limon, Capa_Reten on the bottom part


abscoiS=iner$row.abs #get the absolute contribution of each species on each axis
abscoiS
selectS=rbind(abscoiS[abscoiS[,1]>5,], abscoiS[abscoiS[,2]>5,]) #species which contributes the most to the axis (threshold of 5%)
selectS
# Six species contribute to the first axis with a threshold of 5% : E1 E2 E11  on the right part, and E21 E43 E47 on the left part
# Six species contribute to the second axis with a threshold of 5% : E1 E2 E21  on the top part, and E11 E9 E19 on the bottom part


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

#  coinertia with association #
asso <- assoc %>%
  mutate(
    group = dense_rank(association) # crée une colonne si la station "C1" alors donne le rang "1"
  ) %>%
  select(group)
decoup_asso <- (t(asso)) # passe les lignes deviennent colonnes et colonne deviennent lignes
s.match.class(coi$mX, coi$mY,as.factor(decoup_asso), col1 = c(1:7), col2 = c(1:7), label = c("","","","","","",""))
title(main = "Coinertia with association groups")
