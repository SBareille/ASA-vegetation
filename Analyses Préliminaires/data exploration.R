setwd("C:/Users/serva/Google Drive/1-Partage Ordis/M2/ASA/Projet")

assoc <- read.table("association.txt", header = TRUE)
coord <- read.table("coord.txt", header = TRUE)
flo <- read.table("flo.txt",header = TRUE)
mil <- read.table("mil.txt",header = TRUE)

#################### flo INSPECTION ############################################

names(flo)
str(flo)
effectifs=apply(X = flo, 2, sum)
# get the most represented species independently of associations
effectifs_ordered <- cbind(colnames(flo),as.data.frame(effectifs)) #ATTENTION ! à mettre avant de faire transfo d'hellinger
effectifs_ordered <- effectifs_ordered[order(effectifs_ordered[,2],decreasing=T), ]
effectifs_ordered

par(mfrow=c(1,1))
hist(effectifs
     ,breaks=20
     , col='red'
     , xlab='Species effectif'
     , ylab='Number of Species')


################# Analyse des donnÃ©es environnementales ----
# Transformation des donnÃ©es pour compatiblitÃ© ggplot2
env_data <- mil[-10] %>% 
  pivot_longer(cols = 1:10,  names_to = "Variable", values_to = "Y" ) %>%
  arrange(Variable)


# Affichage ggplot2 boxplots et distribution des donnÃ©es sous forme de scatterplot
ggplot(env_data, aes(x = Variable, y = Y, colour = Variable))+
  geom_point(position = position_jitterdodge(dodge.width = 0.7),size = 2)+
  geom_boxplot(alpha = 0.5)
  
  
################# Assoc inspection ---
summary(assoc)
#get the most represented species in each associations
assoc_number <- assoc
assoc_flo <- cbind(assoc, flo)
assoc_flo_ordered <- assoc_flo[sort(assoc_flo[,1]),]

sort(assoc_flo[,1])

library(dplyr)
assoc_flo %>%
  arrange(association, assoc_flo [,c(2,57)])

assoc_flo %>%
  group_by(association) %>%
  arrange(assoc_flo [,c(2,57)])
  
#Pas fini : faire sommme des effectifs des espèces pour une assoc données (donc tri en fonction de la première colonne)
#Ca a l'air compliqué. A voir si on a le temps
