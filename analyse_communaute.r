library(ade4)
library(vegan)
library(ggplot2)
library(tidyr)
library(dplyr)
  
setwd("~/Etudes/Rennes1/S3/ASA/Project")

assoc <- read.table("association.txt", header = TRUE)
coord <- read.table("coord.txt", header = TRUE)
flo <- read.table("flo.txt",header = TRUE)
mil <- read.table("mil.txt",header = TRUE)



### AFC sur le tableau des occurences d'espèces ----

afc_flo <- dudi.coa(flo, scannf = FALSE, nf = 2)
summary(afc_flo)

s.label(afc_flo$li)
s.label(afc_flo$co)

iner_afc_flo <- inertia.dudi(afc_flo, col.inertia = TRUE, row.inertia = TRUE)
a = format(iner_afc_flo$col.abs,scientific = FALSE)
seuil <- 1/56
seuil <- seuil*100

which(iner_afc_flo$col.abs[,1] > seuil)

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

 

### Analyse des données environnementales ----
env_data <- mil[-10] %>% 
  pivot_longer(cols = 1:10,  names_to = "Variable", values_to = "Y" ) %>%
  arrange(Variable)


ggplot(env_data, aes(x = Variable, y = Y, colour = Variable))+
  geom_point(position = position_jitterdodge(dodge.width = 0.7),size = 2)+
  geom_boxplot(alpha = 0.5)


