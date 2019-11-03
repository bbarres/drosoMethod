##############################################################################/
##############################################################################/
#R code for the Drosophila suzukii pesticide resistance tests paper
##############################################################################/
##############################################################################/

source("droso_data_load.R")


##############################################################################/
#What is the effect of exposure time on the evaluation of LD50?####
##############################################################################/

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

plot(as.data.frame(exporez)$Estimate,ylim=c(0,0.3))

#it seems that EDcomp is not working with LN models (see the help)? Therefore 
#we can do the computation with LL instead (the results are quite similar 
#anyway). We also split male and female experiments because we shown that 
#there was a significant difference of LD50 between the two sexes

expodata_m<-expodata[expodata$sex=="male",]
expo_m_mod<-drm(dead/total~dose,weights=total,
                data=expodata_m,curveid=repet,
                fct=LL.3u(),
                type="binomial")
plot(expo_m_mod,type="confidence")
plot(expo_m_mod,type="obs",add=TRUE)
temp<-EDcomp(expo_m_mod,c(50,50))
as.data.frame(temp)[,4]<0.05
exporez_m<-ED(expo_m_mod,50,interval="delta",reference="control")
plot(as.data.frame(exporez_m)$Estimate,ylim=c(0,0.3))

expodata_f<-expodata[expodata$sex=="female",]
expo_f_mod<-drm(dead/total~dose,weights=total,
                data=expodata_f,curveid=repet,
                fct=LL.3u(),
                type="binomial")
plot(expo_f_mod,type="confidence")
plot(expo_f_mod,type="obs",add=TRUE)
temp<-EDcomp(expo_f_mod,c(50,50))
as.data.frame(temp)[,4]<0.05
exporez_f<-ED(expo_f_mod,50,interval="delta",reference="control")
plot(as.data.frame(exporez_f)$Estimate,ylim=c(0,0.3))


#performing a logistic regression to analyse both the sex and duration of 
#exposure effects at the same time
logmod<-glm(cbind(expodata$alive,expodata$dead)~dose*sex*exposition,
            family=binomial(link=probit),data=expodata)


##############################################################################/
#Effect of the type of sampling environment on the LD50####
##############################################################################/

#fitting the "null hypothesis model"
SmodB0<-drm(dead/total~dose,exposition,
            weights=total,
            data=expodata_f,
            fct=LN.3u(),
            type="binomial")
summary(SmodB0)

#testing the effect of the environment on LD50 (ie the 'e' parameter)
SmodB1env<-drm(dead/total~dose,exposition,
               weights=total,
               data=expodata_f,
               fct=LN.3u(),
               type="binomial",
               pmodels=list(~1, ~1, ~exposition-1))
summary(SmodB1env)
compParm(SmodB1env,"e")
anova(SmodB1env,SmodB0)

#testing the effect of the population on LD50 (ie the 'e' parameter)
SmodB1e<-drm(dead/total~dose,exposition,
             weights=total,
             data=expodata_f,
             fct=LN.3u(),
             type="binomial",
             pmodels=list(~exposition-1, ~1, ~1))
summary(SmodB1e)
anova(SmodB1e,SmodB0)



##############################################################################/
#barplot of an example of evolution of the death rate at the dose 0.25mg/l####
##############################################################################/

data_expo<-read.table("data/droso_expo.txt",header=TRUE,sep="\t")
data_expo<-t(data_expo[,c(6:3,1)])
colnames(data_expo)<-data_expo[5,]

#the supplementary figure for the exposure duration
op<-par(mfrow=c(2,1))
#plot of the results for the male
expobarplot<-barplot(data_expo[,c(1:11)],
                     col=c("black","grey60","grey85"),border=NA,axes=FALSE,
                     axisnames=FALSE,space=0.2,xpd=FALSE)
axis(1,at=expobarplot,labels=FALSE,lwd=4,font=2,
     cex.axis=1.1,padj=0.1,xpd=TRUE,las=1)
text(expobarplot,par("usr")[1]-1,labels=colnames(data_expo),srt=0,
     xpd=TRUE,cex=1.2,font=2)
axis(2,lwd=4,font=2,cex.axis=1.2,las=1)
box(bty="l",lwd=4)
title(main="Male",xlab=NULL,ylab="Number of flies",cex.lab=1.5,
      line=2,font.lab=2,cex.main=3)
text(expobarplot-0.03,as.numeric(data_expo[1,c(1:11)])/2,
     data_expo[1,c(1:11)],font=2,cex=2,xpd=TRUE,col="white")
text(expobarplot-0.03,as.numeric(data_expo[1,c(1:11)]) + 
       as.numeric(data_expo[2,c(1:11)])/2,
     data_expo[2,c(1:11)],font=2,cex=2,xpd=TRUE,col="black")
text(expobarplot-0.03,as.numeric(data_expo[1,c(1:11)]) + 
       as.numeric(data_expo[2,c(1:11)]) + 
       as.numeric(data_expo[3,c(1:11)])/2,
     data_expo[3,c(1:11)],font=2,cex=2,xpd=TRUE,col="black")

#plot of the results for the female
expobarplot<-barplot(data_expo[,c(12:22)],
                     col=c("black","grey60","grey85"),border=NA,axes=FALSE,
                     axisnames=FALSE,space=0.2,xpd=FALSE)
axis(1,at=expobarplot,labels=FALSE,lwd=4,font=2,
     cex.axis=1.1,padj=0.1,xpd=TRUE,las=1)
text(expobarplot,par("usr")[1]-1,labels=colnames(data_expo),srt=0,
     xpd=TRUE,cex=1.2,font=2)
axis(2,lwd=4,font=2,cex.axis=1.2,las=1)
box(bty="l",lwd=4)
title(main="Female",xlab=NULL,ylab="Number of flies",cex.lab=1.5,
      line=2,font.lab=2,cex.main=3)
text(expobarplot-0.03,as.numeric(data_expo[1,c(12:22)])/2,
     data_expo[1,c(12:22)],font=2,cex=2,xpd=TRUE,col="white")
text(expobarplot-0.03,as.numeric(data_expo[1,c(12:22)]) + 
       as.numeric(data_expo[2,c(12:22)])/2,
     data_expo[2,c(12:22)],font=2,cex=2,xpd=TRUE,col="black")
text(expobarplot-0.03,as.numeric(data_expo[1,c(12:22)]) + 
       as.numeric(data_expo[2,c(12:22)]) + 
       as.numeric(data_expo[3,c(12:22)])/2,
     data_expo[3,c(12:22)],font=2,cex=2,xpd=TRUE,col="black")
par(op)

#export to pdf 7 x 14 inches


##############################################################################/
#END
##############################################################################/