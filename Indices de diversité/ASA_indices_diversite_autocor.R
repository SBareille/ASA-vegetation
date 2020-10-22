coord <- read.table("coord.txt", header = TRUE)
flo <- read.table("flo.txt",header = TRUE)
mil <- read.table("mil.txt",header = TRUE)





library(ade4)
library(spdep)
library(tripack)
library(nlme)
library(TSA)
library(gwrr)
library(tidyverse)
require(devtools)
install_version("TSA", version = "1.2.1", repos = "http://cran.us.r-project.org")
install_version("TSA", version = "1.2.1", repos = "https://pbil.univ-lyon1.fr/CRAN/")
library(TSA) 



data=cbind(flo,mil,coord)
plot(data$x,data$y) #representation des coordonnees des stations
xy=data[,68:69] #extraction des coordonnees (x et y)
knear4<-knearneigh(as.matrix(xy),4) # estimation des connexions voisines de type tour par le critere des 4 plus proches voisins
knear4
knn2nb(knear4)




space_mat.bin<- nb2listw(knn2nb(knear4), style="B",zero.policy=TRUE) # matrice des poids spatiale avec la formule des poids binaire
#space_mat.stan<- nb2listw(knn2nb(knear4), style="W", zero.policy=TRUE)

# Calcul de l'indice de diversite de Simpson
Simpson=function(l){
  N=sum(l)
  Si = 0
  for(i in l){
    if(i>0){
    Si=Si+(i*(i-1))
  }}
  return(1-Si/(N*(N-1)))
}

SimFlo=apply(flo, 1, Simpson) # calcul de l'indice de Simpson pour les especes du jeu de donnes
data=cbind(data,SimFlo)



# test de Moran sur l'indice de Simpson
moran.test(data[,70] , space_mat.bin)


#autocorrelogramme de l'indice de Simpson
cor4_flo<-sp.correlogram(knn2nb(knear4), data[,70], order=10, method="I",zero.policy=TRUE, style="B")
print(cor4_flo,p.adjust.method="bonferonni")
plot(cor4_flo,main="Moran autocorrelogram for Simpson index (4 neighbors criterion)")

#c'est bon 



#calculation of Shannon index
Shannon=function(l){
  H=rep(0,nrow(l))
  for (i in 1:nrow(l)){
    for(j in 1:ncol(l)){
      n=l[i,j]
      if(n>0){
        p=n/sum(l[,j])
        H[i]=H[i]-(p)*log2(p)
      }}} 
  return(H)
}

sh_flo=Shannon(flo)

Pielou=function(f){#calcul de l'indice de Pielou
  H=Shannon(f)
  I=rep(0,length(H))
  for(i in 1:length(H)){
    I[i]=H[i]/max(H)
  }
  return(I)
}

P=Pielou(flo)

data=cbind(data,P)

# test de Moran sur l'indice de Shannon
moran.test(data[,70] , space_mat.bin)

cor4_piel<-sp.correlogram(knn2nb(knear4), data[,71], order=10, method="I",zero.policy=TRUE, style="B")
print(cor4_piel,p.adjust.method="bonferonni")
plot(cor4_piel,main="Moran autocorrelogram for Pielou idex (4 neighbors criterion)")





################################################
##BROUILLON, TENTATIVES DE TRUCS
################################################


moddddd=lm(association~Altitude+Argile+Sable+K2O+Conduc+SimFlo+Mg..+Capa_Reten+Na..l, data = div35_2)



moddd=lm(SimFlo~association, data = div35_2)

##test de Moran sur les associations##
#PAS POSSIBLE
assoc <- read.table("association.txt", header = TRUE)
data=cbind(data,assoc)
moran.test(data[,71] , space_mat.bin)

cor4_asso<-sp.correlogram(knn2nb(knear4), div35_2[,71], order=10, method="I",zero.policy=TRUE, style="B")
print(cor4_asso,p.adjust.method="bonferonni")
plot(cor4_asso,main="Moran autocorrelogram for Simpson index on associations (4 neighbors criterion)")