library(ade4)
library(vegan)
library(ggplot2)
library(tidyr)
library(dplyr)
test
setwd("C:/Users/serva/Google Drive/1-Partage Ordis/M2/ASA/Projet")

assoc <- read.table("association.txt", header = TRUE)
coord <- read.table("coord.txt", header = TRUE)
flo <- read.table("flo.txt",header = TRUE)
mil <- read.table("mil.txt",header = TRUE)

#################### flo INSPECTION ############################################

names(flo)
str(flo)
effectifs=apply(X = flo, 2, sum)

par(mfrow=c(1,1))
hist(effectifs
     ,breaks=20
     , col='red'
     , xlab='Species effectif'
     , ylab='Number of Species')


################# Analyse des données environnementales ----
# Transformation des données pour compatiblité ggplot2
env_data <- mil[-10] %>% 
  pivot_longer(cols = 1:10,  names_to = "Variable", values_to = "Y" ) %>%
  arrange(Variable)


# Affichage ggplot2 boxplots et distribution des données sous forme de scatterplot
ggplot(env_data, aes(x = Variable, y = Y, colour = Variable))+
  geom_point(position = position_jitterdodge(dodge.width = 0.7),size = 2)+
  geom_boxplot(alpha = 0.5)

# ACP : 
# Matrice des corrélations
mcor <- cor(mil)
corrplot(mcor, method = "number", type = "upper", order = "hclust", tl.col = "black")

data_ACP <- mil[,-10]
View(data_ACP)

# Réalisation de l'ACP
acp_env <- dudi.pca(data_ACP, scannf = FALSE, nf = 2)
acp_env$eig # Valeurs propres brutes 
acp_env$eig*100/sum(acp_env$eig) # Valeurs propres en pourcentage
s.label(acp_env$li, boxes = FALSE) # Affichage des stations dans le nouveau plan formé par l'acp
s.corcircle(acp_env$co) ## Affichage des variables environnementales dans le nouveau plan de l'acp

iner_acp_env <- inertia.dudi(acp_env, col.inertia = TRUE, row.inertia = TRUE)
iner_acp_env$col.abs # Récupération des contributions absolues pour les colonnes (ie les variables environnementales)

seuil = 1/length(iner_acp_env$col.abs[,1])*100
seuil

row.names(iner_acp_env$col.abs)[which(iner_acp_env$col.abs[,1] > seuil)] # Récupération des variables qui contribuent à l'axe 1
row.names(iner_acp_env$col.abs)[which(iner_acp_env$col.abs[,2] > seuil)] # Récupération des variables qui contribuent à l'axe 2


# Réalisation de la CAH sur l'ACP
dist_env <- dist.quant(acp_env$li,1)
dendo <- hclust(dist_env,"ward.D2")
plot(dendo$height, type="s")
abline(h = 6.7)

plot(dendo)
abline(h = 6.7)

decoup <- cutree(dendo, 6)
decoup

s.class(acp_env$li,as.factor(decoup), clabel = 0.5) # Remarque : pour le cluster 6, seulement présence de deux stations donc pas bien visible sous l'étiquette. Changer le clabel à 0 pour pouvoir voir les points.

# Affichage des étiquettes qui correspondent 
indice_axe1 <- which(iner_apc_env$col.abs[,1] > seuil)
s.label(acp_env$co[indice_axe1,], boxes = FALSE, add.plot = TRUE)

indice_axe2 <- which(iner_acp_env$col.abs[,2] > seuil)
s.label(acp_env$co[indice_axe2,], boxes=TRUE, add.plot = TRUE, clabel = 0.8)




#################### AFC sur le tableau des occurences d'espèces#################### 

afc_flo <- dudi.coa(flo, scannf = FALSE, nf = 2)
summary(afc_flo)

s.label(afc_flo$li)
s.label(afc_flo$co)

iner_afc_flo <- inertia.dudi(afc_flo, col.inertia = TRUE, row.inertia = TRUE)
a = format(iner_afc_flo$col.abs,scientific = FALSE)
seuil_AFC<-1/min(c(nrow(flo),ncol(AFC_flo)))*100

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

# Affichage des étiquettes qui correspondent 
indice_axe1 <- which(iner_afc_flo$col.abs[,1] > seuil)
s.label(afc_flo$co[indice_axe1,], boxes = FALSE)

indice_axe2 <- which(iner_afc_flo$col.abs[,2] > seuil)
s.label(afc_flo$co[indice_axe2,], boxes=TRUE, add.plot = TRUE)

