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
  labs(title = "Spatial representation of asssociations", colour="Vegetal associations") +
  theme(plot.title = element_text(size = 30), 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))
stations_assoc_map

##### Spatial representation of each environemental variables
#create a data frame with stations_name, coord and the deviations to the mean of each observations of each variables (for more lisibility)
mean_mil <- apply(X = mil, c(2), mean)
mil_deviated <- data.frame(matrix(ncol=ncol(mil),nrow=nrow(mil), dimnames=list(NULL, colnames(mil))))
for (j in 1:ncol(mil)) {
  for (i in 1:nrow(mil)) {
    mil_deviated[i,j] <- mil[i,j] - mean_mil[j]}
}
stations_mil_deviated <- cbind(stations_name, coord, mil_deviated)

#map the environmental variables
clay_map <- ggplot(stations_mil_deviated, aes(x=x, y=y, label=rownames(stations_mil_deviated))) + 
  geom_point(shape=21, aes(size=abs(Clay), color = I("black"), fill = as.factor(sign(Clay)))) +
  scale_fill_manual(values = c("white", "black")) +
  ggtitle("Clay")
clay_map

silt_map <- ggplot(stations_mil_deviated, aes(x=x, y=y, label=rownames(stations_mil_deviated))) + 
  geom_point(shape=21, aes(size=abs(Silt), color = I("black"), fill = as.factor(sign(Silt)))) +
  scale_fill_manual(values = c("white", "black")) +
  ggtitle("Silt")
silt_map

sand_map <- ggplot(stations_mil_deviated, aes(x=x, y=y, label=rownames(stations_mil_deviated))) + 
  geom_point(shape=21, aes(size=abs(Sand), color = I("black"), fill = as.factor(sign(Sand)))) +
  scale_fill_manual(values = c("white", "black")) +
  ggtitle("Sand")
sand_map

K2O_map <- ggplot(stations_mil_deviated, aes(x=x, y=y, label=rownames(stations_mil_deviated))) + 
  geom_point(shape=21, aes(size=abs(K2O), color = I("black"), fill = as.factor(sign(K2O)))) +
  scale_fill_manual(values = c("white", "black")) +
  ggtitle("K2O")
K2O_map

Mg_map <- ggplot(stations_mil_deviated, aes(x=x, y=y, label=rownames(stations_mil_deviated))) + 
  geom_point(shape=21, aes(size=abs(Mg.), color = I("black"), fill = as.factor(sign(Mg.)))) +
  scale_fill_manual(values = c("white", "black")) +
  ggtitle("Mg++")
Mg_map

Na100_map <- ggplot(stations_mil_deviated, aes(x=x, y=y, label=rownames(stations_mil_deviated))) + 
  geom_point(shape=21, aes(size=abs(Na..100g), color = I("black"), fill = as.factor(sign(Na..100g)))) +
  scale_fill_manual(values = c("white", "black")) +
  ggtitle("Na+/100g")
Na100_map

K_map <- ggplot(stations_mil_deviated, aes(x=x, y=y, label=rownames(stations_mil_deviated))) + 
  geom_point(shape=21, aes(size=abs(K.), color = I("black"), fill = as.factor(sign(K.)))) +
  scale_fill_manual(values = c("white", "black")) +
  ggtitle("K+")
K_map

Conduc_map <- ggplot(stations_mil_deviated, aes(x=x, y=y, label=rownames(stations_mil_deviated))) + 
  geom_point(shape=21, aes(size=abs(Conduc), color = I("black"), fill = as.factor(sign(Conduc)))) +
  scale_fill_manual(values = c("white", "black")) +
  ggtitle("Conductivity")
Conduc_map

Reten_map <- ggplot(stations_mil_deviated, aes(x=x, y=y, label=rownames(stations_mil_deviated))) + 
  geom_point(shape=21, aes(size=abs(Reten_Capa), color = I("black"), fill = as.factor(sign(Reten_Capa)))) +
  scale_fill_manual(values = c("white", "black")) +
  ggtitle("Retention Capacity")
Reten_map

NaL_map <- ggplot(stations_mil_deviated, aes(x=x, y=y, label=rownames(stations_mil_deviated))) + 
  geom_point(shape=21, aes(size=abs(Na..l), color = I("black"), fill = as.factor(sign(Na..l)))) +
  scale_fill_manual(values = c("white", "black")) +
  ggtitle("Na+/l")
NaL_map

Altitude_map <- ggplot(stations_mil_deviated, aes(x=x, y=y, label=rownames(stations_mil_deviated))) + 
  geom_point(shape=21, aes(size=abs(Altitude), color = I("black"), fill = as.factor(sign(Altitude)))) +
  scale_fill_manual(values = c("white", "black")) +
  ggtitle("Altitude")
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
title(main = "CA with association groups (after Hellinger transformation)") # add title
# We can see which species characterize which associations ! 
# For example : specie E2 caracterize association n°1. 
# It is in adequation with the classification made by the phytosociologists : E2=Sciprus littoralis,   Class 1= facies Sciprus littoralis
# We see that C1 and C2 are close to each other (both are associations of Scirpetum maritimi, the other are Suadetum fructoasae)

######### Second approach : chi square followed by a PCA on the Species occurrence table (with stations, and then associations) ###########

flo <- decostand(flo, "chi.square") # known for a better representation of rare species

PCA_flo <- dudi.pca(flo, scannf = FALSE, nf = 2) # we keep two axis
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
title(main = "PCA on  with association groups (after chi-square transformation)") # add title
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
data_PCA <- env[,-10]
View(data_PCA)

# PCA
PCA_env <- dudi.pca(data_PCA, scannf = FALSE, nf = 2) # we keep 2 axis
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
title(main = "PCA on environmental variables with association group") # add title
# we are gonna group different class according to their possition on the first axis
s.class(PCA_env$li, as.factor(decoup_asso), col=c("blue","blue","blue","green","green","orange","orange"), label = c("C1","C2","C3","C4","C5","C6","C7")) # color = soil type
title(main = "PCA on environmental variables with association group") # add title
# We distinguish 2 ensembles : in blue C1, C2, C3 in humid zone, with a soil of low granulometry, VS in yellow C6, C7 with a sandy soil
# Class C4 and C5 make the transition


#############  coinertia  ####################
coi <- coinertia(PCA_env, PCA_flo, scannf = FALSE, nf=2) # we chose to make it between a PCA on the species table (chi square) and a PCA on environmental variables
coi # we read the RV coefficient :  0.254

randtest(coi,nrepet=1000)
plot(randtest(coi), main = "Monte-Carlo test") # the coinertia analysis is significant

iner=inertia.dudi(coi,col.inertia=T,row.inertia=T)
iner
# The first axis represents 61,72 % of the coinertia, and the second axis 12.99 %

#In particular : 
abscoiE = iner$col.abs #get the absolute contribution of each environmental variable on each axis
abscoiE
selectE = rbind(abscoiE[abscoiE[,1]>(100/10),], abscoiE[abscoiE[,2]>(100/10),]) #variable which contributes the most to the axis (threshold of 10%)
selectE = selectE[-8,] # we get rid of the 8th line because Altitude contributes to axis 1 and 2 and appears two times in the selected varaibles
selectE 
# Four variables contribute to the first axis with a threshold of 10% : Na+/100g   K+   Conduc  on the right part, and Altitude on the left part
# Four variables contribute to this second axis with a threshold of 10% : Sand on the top part, and Altitude, Capa_Reten, Silt on the bottom part

abscoiS=iner$row.abs #get the absolute contribution of each species on each axis
abscoiS
selectS=rbind(abscoiS[abscoiS[,1]>5,], abscoiS[abscoiS[,2]>5,]) #species which contributes the most to the axis (threshold of 5%)
selectS
# Four species contribute to the first axis with a threshold of 5% : E1 E2  on the right part, and E38 E44 on the left part
# Three species contribute to the second axis with a threshold of 5% : E13, E19, E39 on the bottom part

plot(coi)

# clearer representation
#correlation circle of the environmental variables PCA (with the most contributing variables (10% threshold) :
MVA.plot(coi,"corr",space=1) 
s.corcircle(coi$c1[c(2,3,6,7,8,9,10),],box = FALSE, fullcircle = FALSE)
#correlation circle of the species PCA (with the most contributing species (5% threshold) :
MVA.plot(coi,"corr",space=2) 
s.corcircle(coi$l1[c(1,2,38,44,13,19,39),],box=FALSE, fullcircle = FALSE,)

#  coinertia with association #
s.match.class(coi$mX, coi$mY,as.factor(decoup_asso), col1=c(2,6,3,4,5,1,8), col2=c(2,6,3,4,5,1,8), label = c("C1","C2","C3","C4","C5","C6","C7"))
title(main = "Coinertia with association groups", sub = "(PCA on species occurence table and PCA on environmental variables)")
