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



#it seems that EDcomp is not working with LN models (see the help)? Therefore 
#we can do the computation with LL instead (the results are quite similar 
#anyway)
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
#barplot of an example of evolution of the death rate at the dose 0.25mg/l
###############################################################################

temp<-read.table("data/droso_expo.txt",header=TRUE,sep="\t")
tempgraph<-barplot(temp,col=c("black","grey40","grey80"),border=NA,axes=FALSE,
              axisnames=FALSE,space=0.7,xpd=FALSE)
axis(1,at=temp,labels=FALSE,lwd=4,font=2,
     cex.axis=1.1,padj=0.1,xpd=TRUE,las=1)
text(temp,par("usr")[1]-10,labels=names(effectif),srt=25,
     xpd=TRUE,cex=1.2,font=2)
axis(2,lwd=4,font=2,cex.axis=1.2,las=1)
box(bty="l",lwd=4)
title(main=NULL,xlab="HDI class",ylab="% of countries",cex.lab=2,
      line=3.5,font.lab=2)




###############################################################################
#END
###############################################################################