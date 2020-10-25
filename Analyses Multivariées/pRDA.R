################################### RDA ####################################

library(ggplot2)
library(tidyr)
library(dplyr)
library(vegan) #only library for pRDA
library(ade4)
library(RVAideMemoire)

assoc <- read.table("association.txt", header = TRUE)
coord <- read.table("coord.txt", header = TRUE)
flo <- read.table("flo.txt", header = TRUE)
mil <- read.table("mil.txt", header = TRUE)

miltri <- mil[, -10] 

flo <- decostand(flo, "chi.square") # Lots of zero
env= as.data.frame(scale(miltri, center=TRUE, scale=TRUE)) # variables with different units

par(mfrow = c(1,1))
venn <- varpart(flo, mil, coord) #variance partitioning of the different explanatory variables
venn
plot(venn, bg=2:3, Xnames = c("mil","coord"))


rda_VG_vegan <- rda(flo, env) # rda with
rda_VG_vegan
anova.cca(rda_VG_vegan)  # test of rda

#test colinearity between parameters
test_rda <- vif.cca(rda_VG_vegan) 
test_rda # several results above 10 so 
# Clay : 63.798002, Silt:30.779909, Sand:37.710225, K2O : 1.640598
# Mg++: 1.340078,  Na..100g : 2.974436,  K+:2.607641,  Conductivity : 2.614360  
# RetentionCapacity : 1.609666, Altitude :  1.320735                         
          
         
#we have some values greater than ten -> select some variables
rda_VG_vegan <- rda(flo~., env)
ordistep(rda_VG_vegan, perm.max=500)
#we keep : Mg, K+, Capa_Reten, Altitude
selected <- env[,c("Mg..", "K.", "RetentionCapacity", "Altitude","Clay","Silt")] # selection des variables with ordistep
rda_VG_vegan <- rda(flo,selected) # rda with selected environment variable
test_rda <- vif.cca(rda_VG_vegan) 
test_rda
# Mg++: 1.238871,  K+: 1.763315, RetentionCapacity: 1.558781, Altitude : 1.193248, 
# Clay: 2.901147, Silt: 2.019730     
#less than 10 so there is no collinearity              
           

#now, let's do it with ade4
pca_VG_ade4 <- dudi.pca(flo, scannf = FALSE, nf = 2)
rda_VG_ade4 <- pcaiv(pca_VG_ade4, selected, scannf = FALSE, nf = 2)

#now we have our 2 RDA -> let's analyse !
rda_VG_vegan
rda_VG_ade4
MVA.synt(rda_VG_vegan)



#calculate with ADE4 the percentage of species explained by env:
sum(rda_VG_ade4$eig)/sum(pca_VG_ade4$eig)

# 14.1 % of the flo community is explained by the environnement represented by "Mg..", "K.", "RetentionCapacity", "Altitude","Clay","Silt" (p<0,01)    
# Axe 1 : 1.275/9.043 = 43.85% of constrained environment
# Axe 2 : 0.2250/1.275= 17.64 % of constrained environment


randtest(rda_VG_ade4) # permutation test to see if constrained inertia is significant (i.e, impact of environment on communities is significant)
plot(randtest(rda_VG_ade4)) # graphic representation of the randtest. Result : SIGNIFICANT

# same as vegan : decomposition of the inertia for each axis
inertia.dudi(rda_VG_ade4, col = TRUE, row = FALSE)

plot(rda_VG_ade4, boxes = FALSE) # plot RDA

#  RDA with association #
asso <- assoc %>%
  mutate(
    group = dense_rank(association) # crée une colonne si la station "C1" alors donne le rang "1"
  ) %>%
  select(group)
decoup_asso <- (t(asso))

# plot the RDA with group of association
s.match.class(rda_VG_ade4$ls,rda_VG_ade4$li, as.factor(decoup_asso), col1 = c(2,6,3,4,5,1,8), col2 = c(2,6,3,4,5,1,8), label = c("","","","","","",""))


#################### pRDA ###########################################

prda_VG_vegan <- rda(flo, env,coord) # pRDA 
prda_VG_vegan

prda_VG_vegan <- rda(flo~., c(env,coord))
ordistep(prda_VG_vegan, perm.max=500) 
# selection with ordistep : Sand + Mg++ + K+ + RetentionCapacity + Altitude + x + y

selected <- env[,c("Mg..", "K.", "RetentionCapacity", "Altitude","Sand")] # selected environment variable
prda_VG_vegan <- rda(flo,selected,coord) # pRDA with selected variable

# test
anova.cca(prda_VG_vegan)  # the test is significant

prda_VG_vegan # 
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
#interprete the results

