###############################################################################
###############################################################################
#R code for the Drosophila suzukii pesticide resistance tests paper
###############################################################################
###############################################################################

source("droso_data_load.R")


###############################################################################
#What is the effect of exposure time on the evaluation of LD50?
###############################################################################

#we select the data of lambda-cyhalothrin test with the St Foy population
expodata<-dataDroz[dataDroz$expo_comp==1,]

#let's do a model for every repetition
expo_mod<-drm(dead/total~dose,weights=total,
              data=expodata,curveid=repet,
              fct=LN.3u(),
              type="binomial")
plot(expo_mod,type="confidence")
plot(expo_mod,type="obs",add=TRUE)
EDcomp(expo_mod,c(50,50))
exporez<-ED(expo_mod,50,interval="delta",reference="control")


###############################################################################
#END
###############################################################################