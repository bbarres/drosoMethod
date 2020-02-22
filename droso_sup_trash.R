##############################################################################/
#supplementary code####
##############################################################################/

##############################################################################/
#Effect of the exposure on the LD50 estimation: male####
##############################################################################/

#fitting the "null hypothesis model"
expo_mod0m<-drm(dead/total~dose,exposition,
                weights=total,
                data=expodata_m,
                fct=LN.2(),
                type="binomial")
summary(expo_mod0m)
plot(expo_mod0m,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30))

#testing for equality of slope
expo_mod1em<-drm(dead/total~dose,exposition,
                 weights=total,
                 data=expodata_m,
                 fct=LN.2(),
                 type="binomial",
                 pmodels=list(~1, ~exposition-1))
summary(expo_mod1em)
plot(expo_mod1em,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30))
anova(expo_mod1em,expo_mod0m) #there is a significant effect of the slope

#testing for equality of LD50
expo_mod1bm<-drm(dead/total~dose,exposition,
                 weights=total,
                 data=expodata_m,
                 fct=LN.2(),
                 type="binomial",
                 pmodels=list(~exposition-1, ~1))
summary(expo_mod1bm)
plot(expo_mod1bm,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30))
anova(expo_mod1bm,expo_mod0m) #there is a significant effect of the LD50

#comparing the LD50 between the different time of exposure on the full model
compParm(expo_mod0m,"e")
plot(expo_mod0m,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30))


##############################################################################/
#logistic regression with the three factors: sex, dose and exposition####
##############################################################################/

#performing a logistic regression to analyse both the sex and duration of 
#exposure effects at the same time
logmod0<-glm(cbind(expodata$alive,expodata$dead)~dose*sex*exposition,
             family=binomial(link=probit),data=expodata)

logmodN<-glm(cbind(expodata$alive,expodata$dead)~dose+sex+exposition,
             family=binomial(link=probit),data=expodata)

logmodN1<-glm(cbind(expodata$alive,expodata$dead)~dose+sex,
              family=binomial(link=probit),data=expodata)

anova(logmodN1,logmodN,test="Chisq")


##############################################################################/
#supplementary figure for raw results at dose 0.25 mg/L####
##############################################################################/

op<-par(mfcol=c(2,1))
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
#END SUPPLEMENTARY
##############################################################################/