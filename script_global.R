library(ade4)
library(vegan)
library(ggplot2)
library(corrplot)
library(tidyr)
library(dplyr)
library(RVAideMemoire)
library(tripack)
library(spdep)
library(adespatial)
setwd("C:/Users/serva/Google Drive/1-Partage Ordis/M2/ASA/Projet")

# Data extraction
assoc <- read.table("association.txt", header = TRUE)
coord <- read.table("coord.txt", header = TRUE)
flo <- read.table("flo.txt",header = TRUE)
env <- read.table("mil.txt",header = TRUE) # to have english names for variables we changed the .txt : Clay	Silt	Sand	K2O	Mg++	Na+/100g	K+	Conduc	Reten_Capa	Na+/l	Altitude

#######################################################################
#################### PART 1 : data inspection #########################
#######################################################################

#################### flo INSPECTION #########################

names(flo)
str(flo)
effectifs=apply(X = flo, 2, sum)
# get the most represented species
effectifs_ordered <- cbind(colnames(flo),as.data.frame(effectifs))
effectifs_ordered <- effectifs_ordered[order(effectifs_ordered[,2],decreasing=T), ]
effectifs_ordered

par(mfrow=c(1,1))
hist(effectifs
     ,breaks=20
     , col='red'
     , xlab='Species effectif'
     , ylab='Number of Species')

############# Environmental variables inspection ########

names(env)
# We want to see distribution of environmental variables 
# First we transform data to be able to use ggplot2
env_data <- env[-10] %>% 
  pivot_longer(cols = 1:10,  names_to = "Variable", values_to = "Y" ) %>%
  arrange(Variable)

# Distribution of environmental variables 
ggplot(env_data, aes(x = Variable, y = Y, colour = Variable))+
  geom_point(position = position_jitterdodge(dodge.width = 0.7),size = 2)+
  geom_boxplot(alpha = 0.5)
  
########### Spatial representations to understand data #########

##### Spatial representation of associations
#create a data frame with stations_name, coord and assoc 
stations_name <- as.data.frame(row.names(coord))
stations_assoc <- cbind(stations_name, coord, assoc)
stations_assoc
#map
stations_assoc_map <- ggplot(stations_assoc, aes(x=x, y=y, label=rownames(stations_assoc), colour=factor(association))) + 
  geom_text() + 
  scale_color_manual( values=c(2,6,3,4,5,1,8)) +
  geom_label(size=8)+
  labs(title = "Spatial representation of associations", colour="Vegetal associations") +
  theme(plot.title = element_text(size = 30), 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))
stations_assoc_map

##### Spatial representation of each environemental variables
#create a data frame with stations_name, coord and the deviations to the mean of each observations of each variables (for more lisibility)
mean_env <- apply(X = env, c(2), mean)
env_deviated <- data.frame(matrix(ncol=ncol(env),nrow=nrow(env), dimnames=list(NULL, colnames(env))))
for (j in 1:ncol(env)) {
  for (i in 1:nrow(env)) {
    env_deviated[i,j] <- env[i,j] - mean_env[j]}
}
stations_env_deviated <- cbind(stations_name, coord, env_deviated)

#map the environmental variables
clay_map <- ggplot(stations_env_deviated, aes(x=x, y=y, label=rownames(stations_env_deviated))) + 
  geom_point(shape=21, aes(size=abs(Clay), color = I("black"), fill = as.factor(sign(Clay)))) +
  scale_fill_manual(values = c("white", "black")) +
  labs(title = "Clay (%granulometry)",
       subtitle = "Deviation from the mean (mean = 60.38 %)")
clay_map

silt_map <- ggplot(stations_env_deviated, aes(x=x, y=y, label=rownames(stations_env_deviated))) + 
  geom_point(shape=21, aes(size=abs(Silt), color = I("black"), fill = as.factor(sign(Silt)))) +
  scale_fill_manual(values = c("white", "black")) +
  labs(title = "Silt (%granulometry)",
       subtitle = "Deviation from the mean (mean = 30,02 %)")
silt_map

sand_map <- ggplot(stations_env_deviated, aes(x=x, y=y, label=rownames(stations_env_deviated))) + 
  geom_point(shape=21, aes(size=abs(Sand), color = I("black"), fill = as.factor(sign(Sand)))) +
  scale_fill_manual(values = c("white", "black")) +
  labs(title = "Sand (%granulometry)",
       subtitle = "Deviation from the mean (mean = 8,49 %)")
sand_map

K2O_map <- ggplot(stations_env_deviated, aes(x=x, y=y, label=rownames(stations_env_deviated))) + 
  geom_point(shape=21, aes(size=abs(K2O), color = I("black"), fill = as.factor(sign(K2O)))) +
  scale_fill_manual(values = c("white", "black")) +
  labs(title = "K20 (%granulometry)",
       subtitle = "Deviation from the mean (mean = 1,32)")
K2O_map

Mg_map <- ggplot(stations_env_deviated, aes(x=x, y=y, label=rownames(stations_env_deviated))) + 
  geom_point(shape=21, aes(size=abs(Mg..), color = I("black"), fill = as.factor(sign(Mg..)))) +
  scale_fill_manual(values = c("white", "black")) +
  labs(title = "Mg++",
       subtitle = "Deviation from the mean (mean = 14,60 mg/L)")
Mg_map

Na100_map <- ggplot(stations_env_deviated, aes(x=x, y=y, label=rownames(stations_env_deviated))) + 
  geom_point(shape=21, aes(size=abs(Na..100g), color = I("black"), fill = as.factor(sign(Na..100g)))) +
  scale_fill_manual(values = c("white", "black")) +
  labs(title = "Na+/100g",
       subtitle = "Deviation from the mean (mean = 7,56 mg/L)")
Na100_map

K_map <- ggplot(stations_env_deviated, aes(x=x, y=y, label=rownames(stations_env_deviated))) + 
  geom_point(shape=21, aes(size=abs(K.), color = I("black"), fill = as.factor(sign(K.)))) +
  scale_fill_manual(values = c("white", "black")) +
  labs(title = "K+",
       subtitle = "Deviation from the mean (mean = 1,50 mg/L)")
K_map

Conduc_map <- ggplot(stations_env_deviated, aes(x=x, y=y, label=rownames(stations_env_deviated))) + 
  geom_point(shape=21, aes(size=abs(Conduc), color = I("black"), fill = as.factor(sign(Conduc)))) +
  scale_fill_manual(values = c("white", "black")) +
  labs(title = "Conductivity",
       subtitle = "Deviation from the mean (mean = 11,40)")
Conduc_map

Reten_map <- ggplot(stations_env_deviated, aes(x=x, y=y, label=rownames(stations_env_deviated))) + 
  geom_point(shape=21, aes(size=abs(Reten_Capa), color = I("black"), fill = as.factor(sign(Reten_Capa)))) +
  scale_fill_manual(values = c("white", "black")) +
  labs(title = "Retention Capacity",
       subtitle = "Deviation from the mean (mean = 38,55)")
Reten_map

NaL_map <- ggplot(stations_env_deviated, aes(x=x, y=y, label=rownames(stations_env_deviated))) + 
  geom_point(shape=21, aes(size=abs(Na..l), color = I("black"), fill = as.factor(sign(Na..l)))) +
  scale_fill_manual(values = c("white", "black")) +
  labs(title = "Na+/l",
       subtitle = "Deviation from the mean (mean = 63,94)")
NaL_map

Altitude_map <- ggplot(stations_env_deviated, aes(x=x, y=y, label=rownames(stations_env_deviated))) + 
  geom_point(shape=21, aes(size=abs(Altitude), color = I("black"), fill = as.factor(sign(Altitude)))) +
  scale_fill_manual(values = c("white", "black")) +
  labs(title = "Altitude",
       subtitle = "Deviation from the mean (mean = 3,20)")
Altitude_map


#########################################################################
################# PART 2 : statistical analysis #########################
#########################################################################


# standardization of data before statistical analysis
env= as.data.frame(scale(env, center=TRUE, scale=TRUE)) # variables with different units
# first question : choosing the right transformation for our species values. We will try two approaches and compare them to see the one to keep
# First approach : hellinger followed by an CA        VS      2nd approach : chi square followed by a PCA

######### First approach : hellinger followed by an CA on the Species occurrence table (with stations, and then associations) ############# 

flo_helli <- decostand(flo, "hellinger") # Lots of zero

CA_flo <- dudi.coa(flo_helli, scannf = FALSE, nf = 2) # we keep two axis
summary(CA_flo) # the first axis represents 10,51% of the total inertia, and the second axis 7,48%
CA_flo$eig # Gross eigen values of axis
CA_pourcent <- CA_flo$eig*100/sum(CA_flo$eig) # Pourcentages of eigen values of axis
barplot(CA_pourcent, main = "Histogram of Pourcentages of eigen values", ylim=c(0,50), col = heat.colors(10))
s.label(CA_flo$li) # Stations display : there is a Guttman effect
        # shows that there is a strong spatial species gradient, but it also mean that our representation is not optimal
s.label(CA_flo$co) # Species display

# Now we want to see which species contributes to the axis
iner_CA_flo <- inertia.dudi(CA_flo, col.inertia = TRUE, row.inertia = TRUE)
iner_CA_flo$col.abs # Absolute contribution of species on axis
threshold_CA<-1/min(c(nrow(flo_helli),ncol(flo_helli)))*100 # ordinarily used for CA
threshold_CA
row.names(iner_CA_flo$col.abs)[which(iner_CA_flo$col.abs[,1] > threshold_CA)] # Species which contribute the most to axis 1
row.names(iner_CA_flo$col.abs)[which(iner_CA_flo$col.abs[,2] > threshold_CA)] # Species which contribute the most to axis 2
# Almost all species have a  contribution

#  we want to see vegetal associations on the graph to study the species repartition on them
asso <- assoc %>%
  mutate(
    group = dense_rank(association) # creates a column : if station is "C1", then ots rank will be "1"
  ) %>%
  select(group)
decoup_asso <- (t(asso)) # ligns become colomns and colomns become ligns

#s.label(CA_flo$li) # station with CA
s.label(CA_flo$co[c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,17,18,19,23,26,37,38,40,41,43,44,45,46,47,48),]) # species that contributes significantly to axis1 and 2
s.class(CA_flo$li,as.factor(decoup_asso), col=c(2,6,3,4,5,1,8), label = c("C1","C2","C3","C4","C5","C6","C7"), add.plot=T) # add group of associations.
title(main = "CA on species occurence with association groups (after Hellinger transformation)") # add title
# We can see which species characterize which associations ! 
# For example : specie E2 caracterize association n°1. 
# It is in adequation with the classification made by the phytosociologists : E2=Sciprus littoralis,   Class 1= facies Sciprus littoralis
# We see that C1 and C2 are close to each other (both are associations of Scirpetum maritimi, the other are Suadetum fructoasae)

######### Second approach : chi square followed by a PCA on the Species occurrence table (with stations, and then associations) ###########

flo_chi <- decostand(flo, "chi.square") # known for a better representation of rare species

PCA_flo <- dudi.pca(flo_chi, scannf = FALSE, nf = 2) # we keep two axis
summary(PCA_flo) # The first axis represents 8,60% of the total inertia, and the second axis 7,43%. This is less than with hellinger+CA
PCA_flo$eig # Gross eigen values of axis
PCA_pourcent <- PCA_flo$eig*100/sum(PCA_flo$eig) # Pourcentages of eigen values of axis
barplot(PCA_pourcent, main = "Histogram of Pourcentages of eigen values", ylim=c(0,50), col = heat.colors(10))
s.label(PCA_flo$li) # Stations display
s.corcircle(PCA_flo$co) # Species display

# Now we want to see which species contributes to the axis
iner_PCA_flo <- inertia.dudi(PCA_flo, col.inertia = TRUE, row.inertia = TRUE)
iner_PCA_flo$col.abs # Absolute contribution of species on axis
threshold_PCA<-1/length(iner_PCA_flo$col.abs[,1])*100 # ordinarily used for PCA
threshold_PCA
row.names(iner_PCA_flo$col.abs)[which(iner_PCA_flo$col.abs[,1] > threshold_PCA)] # Species which contribute the most to axis 1
row.names(iner_PCA_flo$col.abs)[which(iner_PCA_flo$col.abs[,2] > threshold_PCA)] # Species which contribute the most to axis 2
# Similarly to CA, we have a lot of species which contributes to both axis 
s.corcircle(PCA_flo$co[c(1,3,4,5,6,8,10,11,12,13,15,17,18,20,23,24,25,26,27,33,36,37,38,41,42,44,45,47,48),]) # species that contributes significantly to axis1 and 2

#  we want to see vegetal associations on the graph to study the species repartition on them
s.label(PCA_flo$li) # station with PCA
s.class(PCA_flo$li,as.factor(decoup_asso), col=c(2,6,3,4,5,1,8), label = c("C1","C2","C3","C4","C5","C6","C7"), add.plot=T) # add group of associations.
title(main = "PCA on species occurence table with association groups (after chi-square transformation)") # add title
# Similarly to we can see which species characterize which associations. 
# But because it is cocircle and not label on the same graph it is a bit more difficult to read
# The good thing is that the guttman effect has disappeared
# We can see that C3 seams to be close to C1 and C2, but also to C4.

# Because chi-square transformation enables to remove the Guttman effect, we choose with transformation !

# Now we want to see stations display from an environmental point of view 

################ PCA on environmental variables (with stations, and then associations) ############

# We already have centered and reduce the data

# Correlation matrix
mcor <- cor(env)
corrplot(mcor, method = "number", type = "upper", order = "hclust", tl.col = "black")
# Na+/L and Conductivity are highly correlated (0.96). Na+/L is also correlated with NA+/100g (0.71), as we may have guessed 
# We decide to get rid of Na+/L
env_sorted <- env[,-10]
View(env_sorted)

# PCA
PCA_env <- dudi.pca(env_sorted, scannf = FALSE, nf = 2) # we keep 2 axis
summary(PCA_env) # The first axis represents 39.97% of inertia, and the second axis 15.53%
PCA_env$eig # Gross eigen values of axis
PCA_pourcent <- PCA_env$eig*100/sum(PCA_env$eig) # Pourcentages of eigen values of axis : 
barplot(PCA_pourcent, main = "Histogram of Pourcentages of eigen values", ylim=c(0,50), col = heat.colors(10))
s.label(PCA_env$li, boxes = FALSE) # Stations display
s.corcircle(PCA_env$co) # Environmental variables display

iner_PCA_env <- inertia.dudi(PCA_env, col.inertia = TRUE, row.inertia = TRUE)
iner_PCA_env$col.abs # Absolute contribution of environmental variables on axis
threshold_PCA = 1/length(iner_PCA_env$col.abs[,1])*100
threshold_PCA
row.names(iner_PCA_env$col.abs)[which(iner_PCA_env$col.abs[,1] > threshold_PCA)] # Environmental variables which contributes the most to axis 1
row.names(iner_PCA_env$col.abs)[which(iner_PCA_env$col.abs[,2] > threshold_PCA)] # Environmental variables which contributes the most to axis 2

#  we want to see see vegetal associations on the graph to study the effect of environmental variables on them

s.class(PCA_env$li, as.factor(decoup_asso), col=c(2,6,3,4,5,1,8), label = c("C1","C2","C3","C4","C5","C6","C7")) 
title(main = "PCA on environmental variables with association groups") # add title
# we are gonna group different class according to their possition on the first axis
s.class(PCA_env$li, as.factor(decoup_asso), col=c("blue","blue","blue","green","green","orange","orange"), label = c("C1","C2","C3","C4","C5","C6","C7")) # color = soil type
title(main = "PCA on environmental variables with association groups") # add title
# We distinguish 2 ensembles : in blue C1, C2, C3 in humid zone, with a soil of low granulometry, VS in yellow C6, C7 with a sandy soil
# Class C4 and C5 make the transition

########################################################################
############ Multivariate  spatial autocorrelation #####################
########################################################################

##### Spatial analysis ----
xy <- as.matrix(coord) # transformation of the data frame of coordinates in a matrix type for compatibility with package functions

## Spatial connectivity criteria
nbtri <- tri2nb(xy) # Delaunay triangulation connectivity criterion
nbtri # Output : """ Neighbour list object:
# Number of regions: 97 
# Number of nonzero links: 554 
# Percentage nonzero weights: 5.88798 
# Average number of links: 5.71134 """

nbgab <- graph2nb(gabrielneigh(xy), sym = TRUE) # Gabriel graph connectivity criterion
nbgab # Output : """  Neighbour list object:
# Number of regions: 97 
# Number of nonzero links: 450 
# Percentage nonzero weights: 4.782655 
# Average number of links: 4.639175 """

# Graphical representation functions from the adegraphics package. Takes as first argument the coordinates of the stations
# and as second argument the connectivity criterion. Other arguments are graphical parameters.
delaunaygraph <- adegraphics::s.label(xy, nb = nbtri, main = "Delaunay", pnb.edge.col = "blue", plot = FALSE) # Delaunay connectivity criterion
gabgraph <- adegraphics::s.label(xy, nb = nbgab, main = "Gabriel", pnb.edge.col = "blue", plot = FALSE) # Gabriel connectivity criterion
adegraphics::ADEgS(c(delaunaygraph, gabgraph), layout = c(2,1)) # Graphical representation of both connectivity criteria in one figure with ADEgS function

# The gabriel graph connectivity criterion is the most suited to represent the biological information we want to show.
# More specifically, since stations are well distributed and the species studied are plants, we just want to get the close neighbors that are around each station
# The delaunay graph has a large bias since it connects some stations that are far away and are not really connected except if we would consider water or air dispersion elements

# Obtaining the weight matrix for the gabriel graph connectivity criterion
listwgab <- nb2listw(nbgab, style = "W") # Style W for row standardised method
listwgab # printing the information of the weight matrix

# Moran random test on the line coordinates of the pca. 
# What this does is that it checks if stations that are connected spatially also have close coordinates in the pca respectively for each axis, 
# ie that they have close community composition (since the coordinates in the ordination plan depends on the species composition)
moran.randtest(PCA_flo$li, listw = listwgab) # Output :
# """ class: krandtest lightkrandtest 
# Monte-Carlo tests
# Call: moran.randtest(x = pca_flo$li, listw = listwgab)
# Number of tests:   2 
# Adjustment method for multiple comparisons:   none 
# Permutation number:   999 
# Test       Obs   Std.Obs   Alter Pvalue
# 1 Axis1 0.7335086 10.836005 greater  0.001
# 2 Axis2 0.3160862  4.939626 greater  0.001 """

# Both axes have a significant moran.randtest value (p-value < 0.05) and the moran index is positive in both cases (0.73 for axis 1 ; 0.32 for axis 2). 
# This means that the community composition is spatially structured and more specifically that there is a positive autocorrelation
# This can be visualised thanks to the graphical representation of the values of the coordinates for each axis of the ordination on a spatial map of the stations 
# The first argument of s.value from the adegraphics package is the spatial coordinates of the stations 
# and the second parameter is the coordinates of the stations on the first or the second axis of the pca
g1 <- adegraphics::s.value(coord, PCA_flo$li[,1], symbol = "circle", col = c("black", "white"), ppoint.cex = 0.9, plot = FALSE, main = "axis 1") # First axis
g2 <- adegraphics::s.value(coord, PCA_flo$li[,2], symbol = "circle", col = c("black", "white"), ppoint.cex = 1, plot = FALSE, main = "axis 2") # Second axis
adegraphics::ADEgS(c(g1,g2),layout = c(2,1)) # General figure of the 2 plots

########################################################################
######## Diversity and richness Index and spatial autocorrelation ######
########################################################################
env <- read.table("mil.txt",header = TRUE) # we reload env data to don't have centered and reduced result
data=cbind(flo,env,coord)
plot(data$x,data$y) #representation of station coordinates
# The spatial configurations of the stations is irregular, so we will compute neighbor relationships with an irregular configuration
# Considering we are studying plant species, we will use Gabriel configuration
flo.gab<-gabrielneigh(as.matrix(coord))# We compute a neighbor relationship using a Gabriel graph method

# Spatial representation
s.label(coord,clabel=0, cpoint=1,neig=nb2neig(graph2nb(flo.gab)))
title("Stations with Gabriel graph method")

# Look into Gabriel configuration
graph2nb(flo.gab)
# Neighbour list object:
#   Number of regions: 97 
# Number of nonzero links: 225 
# Percentage nonzero weights: 2.391327 
# Average number of links: 2.319588 
# 9 regions with no links:
#   5 15 31 38 52 62 73 88 97
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

locm.simpson[locm.simpson[,5]<0.05,]
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
# 7 stations have a significant positive value of local index (Pr(z>0)<0.05):
# 43, 44, 45, 53, 65, 66, 89


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
  geom_point(aes(size=RichnessFlo)) +
  ggtitle("Specific Richness")
rich_map

# Spatial representation of Hill Index
stations_Hill=cbind(coord, H)
Hill_map <- ggplot(stations_Hill, aes(x=x, y=y, label=rownames(stations_Hill))) + 
  geom_point(aes(size=H)) +
  ggtitle("Hill Index")
Hill_map

# Graphic representation of the stations with a significant positive value of local index for Hill Index and Specific Richness
stations <- c('43','44','45','46','53','65')
stations_name <- as.data.frame(stations)
stations_assoc <- cbind(stations_name, coord[c(43,44,45,46,53,65),])
stations_assoc
#map
stations_assoc_map <- ggplot(stations_assoc, aes(x=x, y=y, label=rownames(stations_assoc))) + 
  xlim(c(0,max(coord[,1]))) +
  ylim(c(0,max(coord[,2]))) +
  geom_text() + 
  scale_color_brewer(palette = "Set2") +
  geom_label(size=8) +
  labs(title = "Spatial representation of significative local index stations") +
  theme(plot.title = element_text(size = 30), 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))
stations_assoc_map

############################ canonical analysis : RDA  #####################

# Because PCA is better for species representation, we will preform a RDA. 
# But let's choose between RDA and pRDA !

par(mfrow = c(1,1))
venn <- varpart(flo, env, coord) #variance partitioning of the different explanatory variables. 
# On flo and env untransformed (we have collinearity between some env variables, but this is "required" for this function)
venn
plot(venn, bg=2:3, Xnames = c("env","coord"))
# 20 % of species diversity inertia can be explained by environment or spacialisation.
# 6 % is explained by geographical coordinates, 10% by environmental variables alone
# Only 4 % is explained by the interaction between env and coord. But be careful ! This fraction is written as "testable=FALSE"
# Because the interaction between env and coord is not testable we decide to do a RDA and not a pRDA

rda_VG_vegan <- rda(flo_chi, env_sorted) # rda with occurrences transformed with chi square, and centered and reduced environmental variables (without collinearity)
rda_VG_vegan
anova.cca(rda_VG_vegan)  # test the significance of rda. p-value=0,001 -> it is significant

#test colinearity between parameters
test_rda <- vif.cca(rda_VG_vegan)
test_rda 
# Clay : 63.80, Silt:30.78, Sand:37.71, K2O : 1.64
# Mg++: 1.34,  Na..100g : 2.97,  K+:2.61,  Conductivity : 2.61  
# RetentionCapacity : 1.61, Altitude :  1.32                       

#we have some values greater than ten -> select some variables
rda_VG_vegan <- rda(flo_chi~., env_sorted)
ordistep(rda_VG_vegan, perm.max=500)
#we keep : Mg, K+, Capa_Reten, Altitude, Clay, Silt
selected <- env_sorted[,c("Mg..", "K.", "Reten_Capa", "Altitude","Clay","Silt")] # variables selected with ordistep
rda_VG_vegan <- rda(flo_chi,selected) # rda with selected environment variable
test_rda <- vif.cca(rda_VG_vegan) 
test_rda
# Mg++: 1.238871,  K+: 1.763315, RetentionCapacity: 1.558781, Altitude : 1.193248, 
# Clay: 2.901147, Silt: 2.019730     
#less than 10 so there is no collinearity left        

#now, let's do the RDA with ade4
pca_VG_ade4 <- dudi.pca(flo_chi, scannf = FALSE, nf = 2)
rda_VG_ade4 <- pcaiv(pca_VG_ade4, selected, scannf = FALSE, nf = 2)

#now we have our 2 RDA -> let's analyse !
rda_VG_vegan
# 0.141
rda_VG_ade4
MVA.synt(rda_VG_vegan)

#calculate with ADE4 the percentage of species explained by env:
sum(rda_VG_ade4$eig)/sum(pca_VG_ade4$eig)
# 14.1 % of the flo community is explained by the environnement represented by "Mg..", "K.", "RetentionCapacity", "Altitude","Clay","Silt" (p<0,01)    
# Axis 1 : 1.275/9.043 = 43.85% of constrained environment
# Axis 2 : 0.2250/1.275= 17.64 % of constrained environment
randtest(rda_VG_ade4) # permutation test to see if constrained inertia is significant (i.e, impact of environment on communities is significant)
plot(randtest(rda_VG_ade4)) # graphic representation of the randtest. Result : SIGNIFICANT

# same for vegan : decomposition of the inertia for each axis
inertia.dudi(rda_VG_ade4, col = TRUE, row = FALSE)
inertia_rda = inertia.dudi(rda_VG_ade4, col = TRUE, row = FALSE)
seuil = (1/56)*100
row.names(inertia_rda$col.abs)[which(inertia_rda$col.abs[,1] > seuil)] #  variables which contributes the most to axis 1
row.names(inertia_rda$col.abs)[which(inertia_rda$col.abs[,2] > seuil)] #  variables which contributes the most to axis 2
s.corcircle(rda_VG_ade4$co[c(1,2,3,4,9,10,11,12,13,18,19,20,23,26,27,30,32,34,36,37,38,39,40,41,42,43,44,47,51,54,55),])
plot(rda_VG_ade4) # plot RDA

# plot the RDA with group of association
s.match.class(rda_VG_ade4$ls,rda_VG_ade4$li, as.factor(decoup_asso), col1 = c(2,6,3,4,5,1,8), col2 = c(2,6,3,4,5,1,8), label = c("C1","C2","C3","C4","C5","C6","C7"))


# End of analysis 

#########################################################################
############################### Appendix ################################
#########################################################################

#################### pRDA ###########################################

# in case we considered that the 4 % of species variance explained by the interaction between env and coord is significant
prda_VG_vegan <- rda(flo_chi, env_sorted,coord) # pRDA 
prda_VG_vegan

prda_VG_vegan <- rda(flo_chi~., c(env_sorted,coord))
ordistep(prda_VG_vegan, perm.max=500) 
# selection with ordistep : Sand + Mg++ + K+ + RetentionCapacity + Altitude + x + y

selected <- env[,c("Mg..", "K.", "Reten_Capa", "Altitude","Sand")] # selected environmental variables
prda_VG_vegan <- rda(flo_chi,selected,coord) # pRDA with selected variable

# test
anova.cca(prda_VG_vegan)  # the test is significant

prda_VG_vegan 
# the constrained environment corresponds to 10.45% of the total environment
# Axe 1 : 0.3779/0.94543 = 39.97 % of constrained environment
# Axe 2 : 0.1990/0.94543 = 21.05 % of constrained environment

prda_VG_vegan$CCA$v[,c(1,2)] # contributions of flora to axes 1 and 2
prda_VG_vegan$CCA$biplot[,c(1,2)] # contribution of the environment to the axes

# Visualisation graphique
MVA.synt(prda_VG_vegan)
plot(prda_VG_vegan, scaling=1) # scaling 1

plot(prda_VG_vegan, scaling=2) # scaling 2

# scaling 2 
plot(prda_VG_vegan, type="n")
text(prda_VG_vegan, col="blue",cex = 0.75) # station
text(prda_VG_vegan, dis="cn",col="black",cex = 1.2) # environnement
text(prda_VG_vegan, "species", col="red", cex=0.8) # species

# color the vectors according to associations
colasso = NULL
for (i in 1:97){
  if (asso[i,] == 1){
    colasso[i] = 2
  }else{
    if (asso[i,] == 2){
      colasso[i] = 6
    } else {if (asso[i,] == 3){
      colasso[i] = 3
    } else{if (asso[i,] == 4){
      colasso[i] = 4
    }else{if (asso[i,] == 5){
      colasso[i] = 5
    } else{if (asso[i,] == 6){
      colasso[i] = 1
    }else{if (asso[i,] == 7){
      colasso[i] = 8
    }}}}}}}
  
}
colasso
plot(prda_VG_vegan, type="n")
text(prda_VG_vegan, col=colasso,cex = 0.75) # station


############ Multivariate  spatial autocorrelation - suite #####################

# To continue even further the analysis of the spatial effects on community composition a multispati analysis can be done
# The multispati analysis is a multivariate extension of the univariate method of spatial autocorrelation analysis
# It provides a spatial ordination by maximizing the product of variance by spatial autocorrelation.
multi_flo <- adespatial::multispati(PCA_flo, listw = listwgab, scannf = F) # Creating the multispati object based on the gabriel connexion criterion weigth matrix
summary(multi_flo) # Summary of the analysis
# Output : 
# """ Multivariate Spatial Analysis
# Call: adespatial::multispati(dudi = pca_flo, listw = listwgab, scannf = F)
# 
# Scores from the initial duality diagram:
#   var      cum      ratio     moran
# RS1 4.815240 4.815240 0.08598643 0.7335086
# RS2 4.161075 8.976315 0.16029134 0.3160862
# 
# Multispati eigenvalues decomposition:
#   eig      var     moran
# CS1 3.864330 4.593861 0.8411944
# CS2 2.055961 3.286875 0.6255065
# plot(multi_flo) """

# The RS1 and RS2 tests are the same result as the moran.randtest since RS1 and RS2 correspond to the axes of the PCA
# CS1 and CS2 are the new axes of the multispati ordination resulting from the combination of the PCA and spatial data
# The moran index is also significantly positive in each case, meaning there is a positive autocorrelation of community composition
multispati.randtest(PCA_flo, listw = listwgab) # Multivariate autocorrelation test

plot(multi_flo) # Plot of the multispati analysis

# Plot of the stations in the multispati ordination but clustered with the association factor :
ade4::s.match.class(multi_flo$ls, multi_flo$li, as.factor(decoup_asso), col1 = c(2,6,3,4,5,1,8), col2 = c(2,6,3,4,5,1,8), label = c("","","","","","",""))
title(main = "Multispati.PCA on species")