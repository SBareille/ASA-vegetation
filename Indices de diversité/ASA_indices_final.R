########################################################################
######## Diversity and richness Index and spatial autocorrelation ######
########################################################################

library(ade4)
library(spdep)
library(tripack)
library(nlme)
library(TSA)
library(gwrr)
library(tidyverse)
require(devtools)
#install_version("TSA", version = "1.2.1", repos = "http://cran.us.r-project.org")
#install_version("TSA", version = "1.2.1", repos = "https://pbil.univ-lyon1.fr/CRAN/")
library(TSA) 

# Data extraction
coord <- read.table("coord.txt", header = TRUE)
flo <- read.table("flo.txt",header = TRUE)
mil <- read.table("mil.txt",header = TRUE)


data=cbind(flo,mil,coord)
plot(data$x,data$y) #representation of station coordinates
xy=data$x #extraction of the coordinates
xy=cbind(xy,data$y)
# The spatial configurations of the stations is irregular, so we will compute neighbor relationships with an irregular configuration
# Considering we are studying plant species, we will use Gabriel configuration
flo.gab<-gabrielneigh(as.matrix(xy))# We compute a neighbor relationship using a Gabriel graph method

# Spatial representation
s.label(xy,clabel=0, cpoint=1,neig=nb2neig(graph2nb(flo.gab)))
title("Stations with Gabriel graph method")

# Look into Gabriel configuration
graph2nb(flo.gab)
# Neighbour list object:
# Number of regions: 97 
# Number of nonzero links: 225 
# Percentage nonzero weights: 2.391327 
# Average number of links: 2.319588 
# 9 regions with no links:
#  5 15 31 38 52 62 73 88 97
# Non-symmetric neighbours list

flogab.bin <- nb2listw(graph2nb(flo.gab), style="B", zero.policy=TRUE) # Gabriel method with binary formula
flogab.std <- nb2listw(graph2nb(flo.gab), style="W", zero.policy=TRUE) # Gabriel method with standardized formula


# Calculation of Simpson Index
# Index which estimate 2 individuals taken randomly are from the same species
# We take 1-Simpson Index to have an intuitive measure of the diversity (when it is equal to 1, there is a maximum of 
# diversity (equal repartition of individuals between the species)
# when it's equal to 0, the diversity is minimum (1 or a little number of species have all the effectives, when the other species
# have only 1 or some individuals)
# Index directly representative of the heterogeneity
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


# Test of the global spatial aurocorrelation on Simpson Index with Moran test
moran.test(data$SimFlo, flogab.bin, zero.policy=TRUE)

# Moran I statistic standard deviate = 6.5241, p-value = 3.42e-11
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
#         0.397659906      -0.011494253       0.003933059 

# p<0.05 => there is significant global spatial autocorrelation on Simpson Index

moran.test(data$SimFlo, flogab.std, zero.policy=TRUE)

# Moran I statistic standard deviate = 5.156, p-value = 1.262e-07
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
#         0.368299074      -0.011494253       0.005425926 

# Moran I statistic with binary method is greater than standardized method (0.397>0.368) => we will keep binary method


# Due to the large extend of the data set (scale of Mafragh plain) we might expect non stationnary issue
# we will investigate this via the local Moran index
locm.simpson <- localmoran(data$SimFlo, flogab.bin, p.adjust.method="bonferroni")

locm.simpson[lm.simpson[,5]<0.05,]
#             Ii        E.Ii   Var.Ii      Z.Ii    Pr(z > 0)
# 18    5.317127 -0.05208333 4.453161  2.544346 3.284474e-02
# 23    6.788113 -0.04166667 3.598454  3.600383 7.943715e-04
# 24    7.305899 -0.04166667 3.598454  3.873339 2.683858e-04
# 43    8.585015 -0.04166667 3.598454  4.547637 1.356290e-05
# 44   13.465873 -0.04166667 3.598454  7.120628 2.685914e-12
# 45   18.979867 -0.02083333 1.835152 14.025991 1.621364e-44
# 46    4.904896 -0.02083333 1.835152  3.636089 4.152131e-04
# 53    4.454540 -0.02083333 1.835152  3.303644 1.431554e-03

# 8 stations have a significant positive value of local index (Pr(z>0)<0.05):
# 18, 23, 24, 43, 44, 45, 46, 53
# These 8 stations have positive spatial autocorrelation

#calculation of Shannon index
# Index usually calculate to estimate the diversity
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

# calculation of Pielou Index
# The Index has values from 0 to 1, when Shannon Index has values from 0 to log(number of species)
# So 1 is a mean of a maximum diversity
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
moran.test(data$P , flogab.bin, zero.policy=TRUE)
# Moran I statistic standard deviate = 5.8232, p-value = 2.886e-09
# alternative hypothesis: greater
# sample estimates:
# Moran I statistic       Expectation          Variance 
#       0.362850124      -0.011494253       0.004132523

# Moran I statistic>0 and p<0.05 => there is significant positive global spatial autocorrelation on Pielou Index

# Moran local test on Pielou Index
locm.pielou <- localmoran(data$P, flogab.bin, p.adjust.method="bonferroni")

locm.pielou[locm.pielou[,5]<0.05,]
#             Ii        E.Ii    Var.Ii      Z.Ii    Pr(z > 0)
# 24    4.619587 -0.04166667 3.7525129  2.406255 4.029251e-02
# 43    4.852428 -0.04166667 3.7525129  2.526453 2.880516e-02
# 44    5.521656 -0.04166667 3.7525129  2.871926 1.019948e-02
# 45    3.281437 -0.02083333 1.9156436  2.385916 2.555494e-02
# 53    3.609980 -0.02083333 1.9156436  2.623291 1.306271e-02
# 65   16.384299 -0.02083333 1.9156436 11.852839 3.120332e-32
# 66    7.843992 -0.04166667 3.7525129  4.070773 1.171436e-04
# 89    2.595506 -0.01041667 0.9676686  2.649099 8.070673e-03
# 90    3.766791 -0.02083333 1.9156436  2.736589 9.311969e-03

# 9 stations have a significant positive value of local index (Pr(z>0)<0.05):
# 24, 43, 44, 45, 53, 65, 66, 89, 90


# Specific richness calculation
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
moran.test(data$RichnessFlo , flogab.bin, zero.policy=TRUE)
# Moran I statistic standard deviate = 5.346, p-value = 4.497e-08
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
#         0.334862922      -0.011494253       0.004197566

# Moran I statistic>0 and p<0.05 => there is significant positive global spatial autocorrelation on Specific Richness


# Moran local test on Specific Richness
locm.rich <- localmoran(data$RichnessFlo, flogab.bin, p.adjust.method="bonferroni")

locm.rich[locm.rich[,5]<0.05,]
#             Ii        E.Ii    Var.Ii     Z.Ii    Pr(z > 0)
# 43    4.781059 -0.04166667 3.8027496 2.473112 3.348557e-02
# 44    5.800318 -0.04166667 3.8027496 2.995791 6.843341e-03
# 45    4.127342 -0.02083333 1.9418909 2.976766 4.369590e-03
# 53    3.343874 -0.02083333 1.9418909 2.414543 2.363245e-02
# 65   13.234222 -0.02083333 1.9418909 9.511940 2.807080e-21
# 66    4.764983 -0.04166667 3.8027496 2.464868 3.426598e-02
# 89    2.327864 -0.01041667 0.9810745 2.360726 1.823918e-02

# 10 stations have a significant positive value of local index (Pr(z>0)<0.05):
# 6, 24, 43, 44, 45, 46, 53, 54, 55, 65


# Calculation of Hill Index
# While Shannon Index is sensible to rare species, and Simpson Index to abundant species, Hill Index is a
# mathematical relation between these index.
# So It "neutralizes" the sensibility of each index : It's a more synthetic index of the diversity
Hill=function(f){
  sh=Shannon(f)
  si=Simpson(f)
  hi=rep(0,length(sh_flo))
  for(i in 1:length(sh_flo)){
    hi[i]=1/((1+si[i])*exp(sh[i]))
  }
  return(1-hi)
}

H=Hill(flo)

data=cbind(data,H)

# Moran test on Hill Index
moran.test(data$H , flogab.bin, zero.policy=TRUE)
# Moran I statistic standard deviate = 7.4001, p-value = 6.803e-14
# alternative hypothesis: greater
# sample estimates:
#  Moran I statistic       Expectation          Variance 
#        0.468616625      -0.011494253       0.004209256

# Moran I statistic>0 and p<0.05 => there is significant positive global spatial autocorrelation on Specific Richness

# Moran local test on Hill Index
locm.hill <- localmoran(data$H, flogab.bin, p.adjust.method="bonferroni")

locm.hill[locm.hill[,5]<0.05,]
#             Ii        E.Ii   Var.Ii     Z.Ii    Pr(z > 0)
# 6     4.025601 -0.02083333 1.946608 2.900236 5.593235e-03
# 24    9.414943 -0.04166667 3.811778 4.843638 3.187081e-06
# 43    8.957230 -0.04166667 3.811778 4.609199 1.010559e-05
# 44   11.241319 -0.04166667 3.811778 5.779100 1.877531e-08
# 45    7.950316 -0.02083333 1.946608 5.713230 1.662777e-08
# 46    4.587768 -0.02083333 1.946608 3.303162 1.434016e-03
# 53    9.746351 -0.02083333 1.946608 7.000517 3.825296e-12
# 54    5.395024 -0.02083333 1.946608 3.881753 1.555590e-04
# 55    6.083711 -0.03125000 2.889373 3.597426 6.427637e-04
# 65    3.859100 -0.02083333 1.946608 2.780898 8.131321e-03

# 10 stations have a significant positive value of local index (Pr(z>0)<0.05):
# 6, 24, 43, 44, 45, 46, 53, 54, 55, 65

# Spatial representation of Specific Richness
stations_rich=cbind(coord, RichnessFlo)
rich_map <- ggplot(stations_rich, aes(x=x, y=y, label=rownames(stations_rich))) + 
  geom_point(shape=21, aes(size=abs(RichnessFlo), color = I("black"), fill = as.factor(sign(RichnessFlo)))) +
  scale_fill_manual(values = c("white", "black")) +
  ggtitle("Specific Richness")

rich_map

# Spatial representation of Hill Index
stations_Hill=cbind(coord, H)
Hill_map <- ggplot(stations_Hill, aes(x=x, y=y, label=rownames(stations_Hill))) + 
  geom_point(shape=21, aes(size=abs(H), color = I("black"), fill = as.factor(sign(H)))) +
  scale_fill_manual(values = c("white", "black")) +
  ggtitle("Hill Index")

Hill_map


# Graphic representation of the stations with a significant positive value of local index for Hill Index and Specific Richness
stations <- c('6','24','43','44','45','46','53','54','55','65')
stations_name <- as.data.frame(stations)
stations_assoc <- cbind(stations_name, coord[c(6,24,43,44,45,46,53,54,55,65),])
stations_assoc
#map
stations_assoc_map <- ggplot(stations_assoc, aes(x=x, y=y, label=rownames(stations_assoc))) + 
  geom_text() + 
  scale_color_brewer(palette = "Set2") +
  geom_label(size=8) +
  labs(title = "Spatial representation of significative local index stations") +
  theme(plot.title = element_text(size = 30), 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))
stations_assoc_map

