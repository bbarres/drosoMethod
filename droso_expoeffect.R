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



#it seems that EDcomp is not working with LN models ? Therefore we can do
#the computation with LL instead (the results are quite similar anyway)
expo_mod<-drm(dead/total~dose,weights=total,
              data=expodata,curveid=repet,
              fct=LL.3u(),
              type="binomial")
plot(expo_mod,type="confidence")
plot(expo_mod,type="obs",add=TRUE)
temp<-EDcomp(expo_mod,c(50,50))
exporez<-ED(expo_mod,50,interval="delta",reference="control")

plot(as.data.frame(exporez)$Estimate,ylim=c(0,0.3))


#performing a logistic regression to analyse both the sex and duration of 
#exposure effects at the same time
logmod<-glm(cbind(expodata$alive,expodata$dead)~dose*sex*exposition,
            family=binomial(link=probit),data=expodata)

###############################################################################
#END
###############################################################################