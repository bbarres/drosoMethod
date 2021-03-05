##############################################################################/
##############################################################################/
#R code for the Drosophila suzukii pesticide resistance tests: Experiment 1
##############################################################################/
##############################################################################/

#this script provides the code for the analyses of the bioassay results
source("droso_data_load.R")


##############################################################################/
#Is there a difference of LD50 between male and female?####
##############################################################################/

#we select the data of phosmet tests with the St Foy population, with
#flies of 24-48 hours age
sexdata<-dataDroz[dataDroz$number_comp==1,]
#aggregate the data
sexdata<-aggregate(cbind(dead,total)~dose+sex,data=sexdata,"sum")
sexdata_f<-sexdata[sexdata$sex=="female",]
sexdata_m<-sexdata[sexdata$sex=="male",]

#let's do a model for every bioassays combined
sex_mod<-drm(dead/total~dose,weights=total,
             data=sexdata,curveid=sex,
             fct=LN.3u(),
             type="binomial")
#testing the goodness-of-fit of the model
modelFit(sex_mod)
#comparing the LD50
compParm(sex_mod,"e")
sexrez<-ED(sex_mod,50,interval="delta",reference="control")

#comparison of LD50 allowing different slope and "natural death"
sex_mode<-drm(dead/total~dose,weights=total,
             data=sexdata,curveid=sex,
             fct=LN.3u(),
             type="binomial",
             pmodels=list(~factor(sex)-1, ~factor(sex)-1, ~1))
plot(sex_mode)
anova(sex_mod,sex_mode)

#because it seems that there is a bug to display the 95CI with models 
#using curveid, we plot the different modality separately
sex_mod_f<-drm(dead/total~dose,weights=total,
               data=sexdata_f,
               fct=LN.3u(),
               type="binomial")
sex_mod_m<-drm(dead/total~dose,weights=total,
               data=sexdata_m,
               fct=LN.3u(),
               type="binomial")


##############################################################################/
#Code for plotting Figure 1####
##############################################################################/

op<-par(mar=c(5,5,4,1))
plot(sex_mod_f,type="confidence",col="black",bty="n",axes=FALSE,ann=FALSE,
     lwd=3)
#we just double it so that the grey appears darker
plot(sex_mod_f,type="confidence",col="black",add=TRUE)
box(lwd=3,lty=1)
axis(1,at=c(1,10,100,500),labels=c("0","10","100","500"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2)
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("0","20","40","60","80","100"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
plot(sex_mod_f,type="obs",add=TRUE,pch=21,cex=2,
     col=rgb(0,0,0,0.5),bg=rgb(0,0,0,0.5))
plot(sex_mod_m,type="confidence",add=TRUE,col="grey40",lty=2,lwd=3)
plot(sex_mod_m,type="obs",add=TRUE,pch=24,cex=2,
     col=rgb(0,0,0,0.3),bg=rgb(0,0,0,0.0))
title(xlab="Dose (mg/l)",ylab="Mortality rate",cex.lab=2,font.lab=2)
par(op)

#export to pdf 10*7 inches


##############################################################################/
#END
##############################################################################/