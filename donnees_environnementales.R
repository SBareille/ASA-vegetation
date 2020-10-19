library(ade4)
library(vegan)
library(ggplot2)
library(tidyr)
library(dplyr)
library(corrplot)

setwd("~/Etudes/Rennes1/S3/ASA/Project")

assoc <- read.table("association.txt", header = TRUE)
coord <- read.table("coord.txt", header = TRUE)
flo <- read.table("flo.txt",header = TRUE)
mil <- read.table("mil.txt",header = TRUE)


### Analyse des données environnementales ----
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

