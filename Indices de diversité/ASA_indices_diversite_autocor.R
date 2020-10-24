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

coord <- read.table("coord.txt", header = TRUE)
flo <- read.table("flo.txt",header = TRUE)
mil <- read.table("mil.txt",header = TRUE)


data=cbind(flo,mil,coord)
plot(data$x,data$y) #representation of station coordinates
xy=data[,68:69] #extraction of the coordinates
knear4<-knearneigh(as.matrix(xy),4) # rook type of connection approximate by the 4th nearest neighbor criteria
knear4
knn2nb(knear4)

# spatial weighted matrix with binary weight formulation
space_mat.bin<- nb2listw(knn2nb(knear4), style="B",zero.policy=TRUE) 

# Calculation if Simpson Index
Simpson=function(l){
  N=sum(l)
  Si = 0
  for(i in l){
    if(i>0){
    Si=Si+(i*(i-1))
  }}
  return(1-Si/(N*(N-1)))
}

# Calculation of Simpson Index for plant species
SimFlo=apply(flo, 1, Simpson)
data=cbind(data,SimFlo)



# Moran test on Simpson Index
moran.test(data[,70] , space_mat.bin)


#autocorrelogramm of Simson Index
cor4_flo<-sp.correlogram(knn2nb(knear4), data[,70], order=10, method="I",zero.policy=TRUE, style="B")
print(cor4_flo,p.adjust.method="bonferonni")
plot(cor4_flo,main="Moran autocorrelogram for Simpson Index (4 neighbors criterion)")


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

#calculation of Pielou Index
Pielou=function(f){
  H=Shannon(f)
  I=rep(0,length(H))
  for(i in 1:length(H)){
    I[i]=H[i]/max(H)
  }
  return(I)
}

P=Pielou(flo)

data=cbind(data,P)

# Moran test on Pielou Index
moran.test(data[,71] , space_mat.bin)

cor4_piel<-sp.correlogram(knn2nb(knear4), data[,71], order=10, method="I",zero.policy=TRUE, style="B")
print(cor4_piel,p.adjust.method="bonferonni")
plot(cor4_piel,main="Moran autocorrelogram for Pielou Index (4 neighbors criterion)")






# specifique fichness calculation
richness=function(station){
  counter = 0
  for (i in station){
    if (i == 0){
      counter = counter + 1
    }
  }
  specific_richness = length(station)-counter
  return(specific_richness)
}

# Calculation of Specific richness for plant species
RichnessFlo=apply(flo, 1, richness)
data=cbind(data,RichnessFlo)


# Moran test on Specific richness
moran.test(data[,72] , space_mat.bin)

cor4_piel<-sp.correlogram(knn2nb(knear4), data[,72], order=10, method="I",zero.policy=TRUE, style="B")
print(cor4_piel,p.adjust.method="bonferonni")
plot(cor4_piel,main="Moran autocorrelogram for the specific richness (4 neighbors criterion)")

