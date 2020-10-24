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

#Autocorrelogram of Pielou Index
cor4_piel<-sp.correlogram(knn2nb(knear4), data[,71], order=10, method="I",zero.policy=TRUE, style="B")
print(cor4_piel,p.adjust.method="bonferonni")
plot(cor4_piel,main="Moran autocorrelogram for Pielou Index (4 neighbors criterion)")






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
moran.test(data[,72] , space_mat.bin)

cor4_piel<-sp.correlogram(knn2nb(knear4), data[,72], order=10, method="I",zero.policy=TRUE, style="B")
print(cor4_piel,p.adjust.method="bonferonni")
plot(cor4_piel,main="Moran autocorrelogram for the specific richness (4 neighbors criterion)")



#######################################################################################################
#Moran Mais en Local (j'ai l'impression que c'est mieux car notre domaine d'étude est vraiment large)

# avec pielou
localI.pp1=localmoran(data[,71], space_mat.bin, p.adjust.method="bonferroni")
localI.pp1[localI.pp1[,5]<0.05,]
s.value(xy,localI.pp1[,1], ylim=c(min(data[,69])-15,max(data[,69])+15), xlim=c(min(data[,68])-15,max(data[,68])+15), sub="")


# avec richesse
localI.pp1=localmoran(data[,72], space_mat.bin, p.adjust.method="bonferroni")
localI.pp1[localI.pp1[,5]<0.05,]
s.value(xy,localI.pp1[,1], ylim=c(min(data[,69])-15,max(data[,69])+15), xlim=c(min(data[,68])-15,max(data[,68])+15), sub="")

# Globalement y'a des zones plus homogene en terme de richesse et diversité mais en soit je vois pas trop ce que ça nous apporte de savoir ça...




# location data 
xy=data[,68:69]

# Our points are not regularly distributed through space.We computed a neighbor relationship using several methods

# example with Gabriel method:
data.gab<-gabrielneigh(as.matrix(xy))
s.label(xy,clabel=0, cpoint=1,neig=nb2neig(graph2nb(data.gab)))

# truc écrit dans le script originel mais interessant : 
# the part nb2neig(graph2nb(data.gab)) is necessary to translate the object from one package to another 
# s.lable is a function from ade4 and yo need to translate your object from tripack->spdep->ade4



# synthetic view of the neighbor relationship:
graph2nb(data.gab)


# Methode de Delaunay:
# to compute a neighbor relationship using the delaunay triangulation method you should use this function using dirrectly the coordinates:
data.tri<-tri2nb (xy)

# graphic representation:
s.label(xy,clabel=0,cpoint=1,neig=nb2neig(data.tri))

# synthetic view:
data.tri


# encore un autre truc du même genre:

# to compute a neighbor relationship using a distance criteria you should use this function 
# using dirrectly the coordinates as a matrix and indicating the min and max of distance.
# the max of distance is defined depending on your knowledge of the system or for different hypothesis you want to test:
data.dnear<-dnearneigh(as.matrix(xy),0,30)
# graphic representation:
s.label(xy,clabel=0, cpoint=1,neig=nb2neig(data.dnear))
# synthetic view:
data.dnear

# je capte pas trop pk on fait 3 trucs comme ça..


# when comparing gabriel, triangulation and distance results: what would you say? 
# va niker ta mere?

# what kind of biological process does each criteria can represent?
# et bah pk tu nous filles pas la reponse ? tas cru j'ai le temps d'ecrire en cours ? 

# C'est le poids spatial avec diverses methodes de calcules (renommées)
ponddatatri.bin <- nb2listw(data.tri, style="B",zero.policy=TRUE)
ponddatagab.bin <- nb2listw(graph2nb(data.gab), style="B", zero.policy=TRUE)
ponddatadnear.bin <- nb2listw(data.dnear, style="B",zero.policy=TRUE)

# On teste avec Moran pour voir si c'est Pvalue < 0.05
PVALUES = c()
for (i in 57:67 ){
  a = moran.test(data[,i], ponddatagab.bin, zero.policy=TRUE)$p.value
  b = moran.test(data[,i], ponddatatri.bin, zero.policy=TRUE)$p.value
  c = moran.test(data[,i], ponddatadnear.bin, zero.policy=TRUE)$p.value
  PVALUES = c(PVALUES,names(data)[i],a,b,c)
}
for (i in 70:72 ){
  a = moran.test(data[,i], ponddatagab.bin, zero.policy=TRUE)$p.value
  b=moran.test(data[,i], ponddatatri.bin, zero.policy=TRUE)$p.value
  c= moran.test(data[,i], ponddatadnear.bin, zero.policy=TRUE)$p.value
  PVALUES = c(PVALUES,names(data)[i],a,b,c)
}


PVALUES
# j'ai l'implression qu'il y a une forte autocorrelation spaciale pour tout le jeux de données. A la fois
 # pour les fleurs et aussi pour les parametre phisico chimiques...


# alors je tente un autre truc nouveaux là
# Je crois que c'est pour chercher en un seul coup si y'a correlation dans le jeu de données

# Dans le ''''cours'''' j'ai trouvé ça : 
# Is there any similary between environmental variables of ponds explained by their geographical proximity?
#   Mantel test: function mantel.rtest
# m'enfin c'est pas clair , c'est pas détaillé, c'est ce foutre de notre geule sérieux

# en gros on regarde les distance entre les stations, puis on compute les ressemblances phisico chimique et on regarde si les 
# stations les plus proches se ressembles etc..

geo = dist(xy)     # matrices des distances
chim = dist(scale(data[,57:67]))   # matrice des ressemblance phisico chimiques
# faut peut etre sortir les variables correlées qui doivent rajouter de l'auto cor ?

r1 = mantel.rtest(geo,chim,nrepet=1000)  # le fameux teste ...

r1  # p-value: 0.02297702
plot(r1, main = "Mantel's test")

# what can you conclude?
# je sais pas. 






##############Spatial modelling: practice######################

# we have seen just before taht there is autocorrelation at different scale in the data of planta nd bird richnesses in Ille & vilaine
# we are now interested in the relationship that could occur between plant richness and bird richness taking into account the spatial effect
# we will use the queen criteria with standardized weight (why? because this is my choice.... you can try with other criteria): 

# Ici je pense qu'on a pas comme la prof 2 jeux de données à comparer. Elle elle a genre richesse oiseaux et richesse plante
# Elle regarde si ya même auto cor ou un truc dans le genre. Nous on a Richesse spécifique, Mais faudrait identifié  1  seul param pour tester avec

div.knear8<-knearneigh(as.matrix(xy),8)
pond8.stand<- nb2listw(knn2nb(div.knear8), style="W",zero.policy=TRUE)
moran.test(data$RichnessFlo, pond8.stand)


# Start with a linear model that explain bird richness with plant richness and check if residuals are spatially autocorrelated
# Je choisis K car les Modeles lineraires de tout à l'heure le faisait ressortir
lmO=lm(data$RichnessFlo~data$K.)   # j'ai plein de Warning je les comprends pas
summary(lmO)
moran.test(lmO$residuals,pond8.stand)  #p-value = 0.0002419 

# what is the relationship estimated between plant and bird?  

# are the residuals autocorrelated?  Je crois... 



# to define which spatial models you have to do you have to perform the Lagrange multiplier diagnostic
# using the linear model, and the spatial criteria
lm.LMtests(lmO,pond8.stand, test="all")
# which tests are significant?

# so which model you must do?


#### Je me suis arrété là



# the function to performe the model is the one below, you use the same formula as for linear model 
# but you add an argument concerning the spatial connection of your points
lmsarD=sacsarlm(rich_ois~rich_pp,data=div35,pond8.stand)

# before analysing the results : check if the residuals are always spatially autocorrelated
ressar=lmsarD$residuals
moran.test(ressar, pond8.stand, alternative="two.sided")
# are they?

# to analyze the result of your model:
summary(lmsarD)

# can you repace the estimated coefficients in the equation related to this model: y = ? W y + ?X + ? W e + ? ?
# can you repace the estimated coefficients in the equation related to this model: y = 0.63 W y + ?X + -0.65643 W e + ? ?
#??? c'est rrho , lambda       peut etre que beta c'est estimate
# you can plot the results, and comapre with teh line obtained with the non spatial model:
x11()
plot(div35$rich_ois~div35$rich_pp, main="Relationship between plant and bird richnesses", xlab="plant richness", ylab="bird richness")
abline(lm(lmsarD$fitted.values~div35$rich_pp),col="red")
abline(lm(lmO$fitted.values~div35$rich_pp))

# what to conclude?


# however we have seen with local Moran index that the spatial autocorrelation was non stationary
# now we will investigated if the relationship between plant and bird is the same everywhere in Ille & Vilaine

# for that we used geographically weighted regression in order to estiamte a regression coefficient per point and map its value
# for this function we do not use connection critria, 
# R will use a kernel function to estimate the spatial weight of the connection between points (like in geastatistical tool of analysis)
lmgwrD=gwr.est(rich_ois~rich_pp,data=div35,locs=xy)
summary(lmgwrD)
lmgwrD$beta



# Phi: kernel bandwidth
# RMSPE:  root square of error of the prediction based on bandwidth estimated
# Beta:  matrix with the coefficients for each point (intercept line [1,] and regression coefficient line [2,])
# Yhat:  model estimate
# RMSE:  root square of error estimation 
# Rsquare:  model fitting

# what is the range of variation of the regression between bird richness and plant richness?

#you can map the value of the coefficient:

plot(xy)
s.value(xy,lmgwrD$beta[2,],add.plot=T)

