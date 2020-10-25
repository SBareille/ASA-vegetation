# Import of the packages
library(ade4)
library(vegan)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(spdep)
library(adespatial)

# Setting the working directory 
setwd("~/Etudes/Rennes1/S3/ASA/Projet/ASA-vegetation-main/data")

# Import of the data as data frames
assoc <- read.table("association.txt", header = TRUE) # Associations for each station
coord <- read.table("coord.txt", header = TRUE) # xy coordinates of the stations
flo <- read.table("flo.txt", header = TRUE) # species occurence data for each station
mil <- read.table("mil.txt", header = TRUE) # environmental data for each station


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


# PCA on the species matrix to structure the stations by species (and therefore represents community)
flo_Xhi2 <- decostand(flo, "chi.square") # xhi2 transformation since species that are rare 
pca_flo <- dudi.pca(flo_Xhi2, scannf = FALSE, nf = 2) # pca on the species matrix after chi.square transformation. Note : dudi.coa scales and centers the data.2 axes are kept with the argument nf = 2

# Moran random test on the line coordinates of the pca. 
# What this does is that it checks if stations that are connected spatially also have close coordinates in the pca respectively for each axis, 
# ie that they have close community composition (since the coordinates in the ordination plan depends on the species composition)
moran.randtest(pca_flo$li, listw = listwgab) # Output :
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
g1 <- adegraphics::s.value(coord, pca_flo$li[,1], symbol = "circle", col = c("black", "white"), ppoint.cex = 0.9, plot = FALSE, main = "axis 1") # First axis
g2 <- adegraphics::s.value(coord, pca_flo$li[,2], symbol = "circle", col = c("black", "white"), ppoint.cex = 1, plot = FALSE, main = "axis 2") # Second axis
adegraphics::ADEgS(c(g1,g2),layout = c(2,1)) # General figure of the 2 plots


# To continue even further the analysis of the spatial effects on community composition a multispati analysis can be done
# The multispati analysis is a multivariate extension of the univariate method of spatial autocorrelation analysis
# It provides a spatial ordination by maximizing the product of variance by spatial autocorrelation.
multi_flo <- adespatial::multispati(pca_flo, listw = listwgab, scannf = F) # Creating the multispati object based on the gabriel connexion criterion weigth matrix
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
multispati.randtest(pca_flo, listw = listwgab) # Multivariate autocorrelation test

plot(multi_flo) # Plot of the multispati analysis

# Grouping of the stations by the association factor 
asso <- assoc %>%
  mutate(
    group = dense_rank(association) 
  ) %>%
  select(group)
decoup_asso <- (t(asso))

# Plot of the stations in the multispati ordination but clustered with the association factor :
ade4::s.match.class(multi_flo$ls, multi_flo$li, as.factor(decoup_asso), col1 = c(2,6,3,4,5,1,8), col2 = c(2,6,3,4,5,1,8), label = c("","","","","","",""))
title(main = "Multispati.PCA on species")


