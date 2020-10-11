##############################################################################/
##############################################################################/
#R code for the Drosophila suzukii pesticide resistance: Experiment 5
##############################################################################/
##############################################################################/

source("droso_data_load.R")


##############################################################################/
#What is the effect of exposure time on the evaluation of LD50?####
##############################################################################/

#we select the data of lambda-cyhalothrin test with the St Foy population
expodata<-dataDroz[dataDroz$expo_comp==1,]

#because there is a strong effect of sex and we are mainly interested in the
#effect on the female, we split the dataset according to sex
expodata_f<-expodata[expodata$sex=="female",]
rep0<-expodata_f[expodata_f$date=="13/12/17",]
rep1<-expodata_f[expodata_f$date=="09/09/20",]
rep2<-expodata_f[expodata_f$date=="24/09/20",]
rep3<-expodata_f[expodata_f$date=="30/09/20",]
rep4<-expodata_f[expodata_f$date=="01/10/20",]
rep5<-expodata_f[expodata_f$date=="07/10/20",]
rep6<-expodata_f[expodata_f$date=="08/10/20",]
expodata_m<-expodata[expodata$sex=="male",]

#loading the data of an example of evolution of the death rate at the dose 
#0.25mg/L
data_expo<-read.table("data/droso_expo.txt",header=TRUE,sep="\t")
data_expo<-t(data_expo[,c(6:3,1)])
colnames(data_expo)<-data_expo[5,]


##############################################################################/
#Effect of the exposure on the LD50 estimation: female####
##############################################################################/

#fitting the "null hypothesis model"
expo_mod0<-drm(dead/total~dose,exposition,
               weights=total,
               data=expodata_f,
               fct=LN.2(),
               type="binomial")
summary(expo_mod0)
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30))

#testing for equality of slope
expo_mod1e<-drm(dead/total~dose,exposition,
                weights=total,
                data=expodata_f,
                fct=LN.2(),
                type="binomial",
                pmodels=list(~1, ~exposition-1))
summary(expo_mod1e)
plot(expo_mod1e,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30))
anova(expo_mod1e,expo_mod0) #there is a significant effect of the slope

#testing for equality of LD50
expo_mod1b<-drm(dead/total~dose,exposition,
                weights=total,
                data=expodata_f,
                fct=LN.2(),
                type="binomial",
                pmodels=list(~exposition-1, ~1))
summary(expo_mod1b)
plot(expo_mod1b,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30))
anova(expo_mod1b,expo_mod0) #there is a significant effect of the LD50

#comparing the LD50 between the different time of exposure on the full model
compParm(expo_mod0,"e")
exporez_f<-ED(expo_mod0,50,interval="delta",reference="control")
op<-par(mar=c(5.1,5.1,4.1,2.1))
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30),lwd=1.5,
     legendPos=c(15,0.7),xlab="dose (mg/L)",cex.axis=1.5,cex.lab=2,cex=2)
arrows(x0=expo_mod0$parmMat[2,1],y0=0.5,
       x1=expo_mod0$parmMat[2,5],y1=0.5,
       length=0.12,angle=25,lwd=3)
par(op)

#export to pdf 7 x 7 inches


##############################################################################/
#Effect of the exposure on the LD50 estimation: female by rep####
##############################################################################/

#fitting the "null hypothesis model" rep0
expo_mod0<-drm(dead/total~dose,exposition,
               weights=total,
               data=rep0,
               fct=LN.2(),
               type="binomial")
summary(expo_mod0)
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30),
     main="rep0 - 13/12/17 - alive")

#fitting the "null hypothesis model" rep1
expo_mod0<-drm(dead/total~dose,exposition,
               weights=total,
               data=rep1,
               fct=LN.2(),
               type="binomial")
summary(expo_mod0)
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30),
     main="rep1 - 09/09/20 - alive")

#fitting the "null hypothesis model" rep2
expo_mod0<-drm(dead/total~dose,exposition,
               weights=total,
               data=rep2,
               fct=LN.2(),
               type="binomial")
#no convergence for several time (21h, 22h and 23h)
#we remove the problematic time
expo_mod0<-drm(dead/total~dose,exposition,
               weights=total,
               data=rep2[rep2$exposition!="21h"
                         & rep2$exposition!="22h"
                         & rep2$exposition!="23h",],
               fct=LN.2(),
               type="binomial")
summary(expo_mod0)
plot(expo_mod0,col=c(1,2,2,1,1,1,1),xlim=c(0,30),
     main="rep2 - 24/09/20 - alive")

#fitting the "null hypothesis model" rep3
expo_mod0<-drm(dead/total~dose,exposition,
               weights=total,
               data=rep3,
               fct=LN.2(),
               type="binomial")
summary(expo_mod0)
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30),
     main="rep3 - 30/09/20 - alive")

#fitting the "null hypothesis model" rep4
expo_mod0<-drm(dead/total~dose,exposition,
               weights=total,
               data=rep4,
               fct=LN.2(),
               type="binomial")
summary(expo_mod0)
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30),
     main="rep4 - 01/10/20 - alive")

#fitting the "null hypothesis model" rep5
expo_mod0<-drm(dead/total~dose,exposition,
               weights=total,
               data=rep5,
               fct=LN.2(),
               type="binomial")
summary(expo_mod0)
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30),
     main="rep5 - 07/10/20 - alive")

#fitting the "null hypothesis model" rep4
expo_mod0<-drm(dead/total~dose,exposition,
               weights=total,
               data=rep6,
               fct=LN.2(),
               type="binomial")
summary(expo_mod0)
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30),
     main="rep6 - 08/10/20 - alive")


##############################################################################/
#Figure 6: final plots exemplifying effect of time of exposure
##############################################################################/

#plot of the raw results for the female at dose 0.25 mg/L
op<-par(mar=c(5.1,5.1,4.1,2.1))
expobarplot<-barplot(data_expo[,c(12:22)],col=c("black","grey60","grey85"),
                     border=NA,axes=FALSE,axisnames=FALSE,space=0.2,
                     xpd=FALSE)
axis(1,at=expobarplot,labels=FALSE,lwd=4,font=2,
     cex.axis=1.5,padj=0.1,xpd=TRUE,las=1)
text(expobarplot,par("usr")[1]-1.3,labels=colnames(data_expo),srt=0,
     xpd=TRUE,cex=1.4,font=2)
axis(2,lwd=4,font=2,cex.axis=1.5,las=1)
box(bty="l",lwd=4)
title(main="dose = 0.25 mg/L",xlab=NULL,ylab="Number of flies",cex.lab=2,
      line=2.5,font.lab=2,cex.main=1.5)
text(expobarplot-0.03,as.numeric(data_expo[1,c(12:22)])/2,
     data_expo[1,c(12:22)],font=2,cex=2,xpd=TRUE,col="white")
text(expobarplot-0.03,as.numeric(data_expo[1,c(12:22)]) + 
        as.numeric(data_expo[2,c(12:22)])/2,
     data_expo[2,c(12:22)],font=2,cex=2,xpd=TRUE,col="black")
text(expobarplot-0.03,as.numeric(data_expo[1,c(12:22)]) + 
        as.numeric(data_expo[2,c(12:22)]) + 
        as.numeric(data_expo[3,c(12:22)])/2,
     data_expo[3,c(12:22)],font=2,cex=2,xpd=TRUE,col="black")
text(-2,27.5,labels=c("A"),cex=4,xpd=TRUE)
par(op)
#export to pdf 7 x 7 inches

#comparison of the regression curves for the different reading time
op<-par(mar=c(5.1,5.1,4.1,2.1))
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30),lwd=1.5,bp=1e-3,
     legendPos=c(15,0.7),xlab="Dose (mg/L)",cex.axis=1.5,cex.lab=2,
     cex=2,axes=FALSE,font.lab=2,font.axis=2,font=2,bty="n")
arrows(x0=expo_mod0$parmMat[2,1],y0=0.5,
       x1=expo_mod0$parmMat[2,4],y1=0.5,
       length=0.12,angle=25,lwd=3)
axis(1,at=c(0.001,0.01,0.1,1,10),labels=c("0","0.01","0.1","1","10"),
     lwd=4,font=2,cex.axis=1.5,padj=0.1,xpd=TRUE,las=1)
axis(2,lwd=4,font=2,cex.axis=1.5,las=1)
box(bty="l",lwd=4)
text(0.17*10^-3,1.148,labels=c("B"),cex=4,xpd=TRUE)
par(op)
#export to pdf 7 x 7 inches

#heatmap displaying the level of significance of the different LD50 at 
#different reading time
temp<-compParm(expo_mod0,"e")
temp<-cbind(matrix(unlist(strsplit(row.names(temp),"/")),45,byrow=TRUE),
            temp[,4])
temp<-rbind(temp,temp[,c(2,1,3)])
temp<-spread(data.frame(temp),2,3,fill="NA",convert=TRUE)
#changing the rownames
row.names(temp)<-as.character(temp[,1])
#removing the first unnecessary column
temp<-temp[,c(-1)]
#reordering the columns and turning the object into a matrix
temp<-as.matrix(temp[c(1,7,8,9,10,2:6),c(1,7,8,9,10,2:6)])
#scaling the p-value so it is easily usable with the LDheatmap function
temp[temp>0.5]<-0.9
temp[temp>0.10 & temp<0.9]<-0.8
temp[temp>0.05 & temp<0.8]<-0.6
temp[temp>0.01 & temp<0.6]<-0.4
temp[temp>0.001 & temp<0.4]<-0.2
temp[temp<0.001]<-0.1
#the actual plotting start here
chaudemap<-LDheatmap(temp,title=NULL,
                     add.map=FALSE,distances=NULL,SNP.name=row.names(temp),
                     color=c(rep(grey(0.8),3),
                             brewer.pal(6,"YlOrRd")[c(2,4,6)]),
                     name="CHR",flip=FALSE,add.key=FALSE)
grid.edit(gPath("CHR","heatMap","heatmap"),gp=gpar(col="white",lwd=1))
grid.edit(gPath("CHR","SNPnames"),
          gp=gpar(col="black",rot="0",cex=0.9,font=2),
          rot=0,hjust=0.6)
grid.lines(x=unit(c(0.1,0.5),"npc"),y=unit(c(0.5,0.5),"npc"),
           gp=gpar(lwd=3))
grid.lines(x=unit(c(0.5,0.5),"npc"),y=unit(c(0.5,0.9),"npc"),
           gp=gpar(lwd=3))
grid.text("C",x=unit(0.048,"npc"),y=unit(0.955,"npc"),gp=gpar(fontsize=50),
          check=TRUE)
#export to pdf 7 x 7 inches


##############################################################################/
#Effect of the exposure on the LD50 estimation: female moribund=alive####
##############################################################################/

#load the dataset
dataDrozbis<-read.table("data/droso_datad_DD.txt",header=T,sep="\t")
# #we remove the two concentrations that were used at the beginning of the 
# #test when we were still adjusting the range of doses for the test
# dataDroz<-dataDroz[dataDroz$dose!=603.70 & dataDroz$dose!=301.85,]
#creation of variable to distinguish between male and female and time 
#of exposure to pesticide
dataDrozbis<-cbind(dataDrozbis,
                   "repet"=paste(dataDrozbis$date,dataDrozbis$sex,
                                 dataDrozbis$exposition))
#we select the data of lambda-cyhalothrin test with the St Foy population
expodatabis<-dataDrozbis[dataDrozbis$expo_comp==1,]

#because there is a strong effect of sex and we are mainly interested in the
#effect on the female, we split the dataset according to sex
expodatabis_f<-expodatabis[expodatabis$sex=="female",]
rep0bis<-expodatabis_f[expodatabis_f$date=="13/12/17",]
rep1bis<-expodatabis_f[expodatabis_f$date=="09/09/20",]
rep2bis<-expodatabis_f[expodatabis_f$date=="24/09/20",]
rep3bis<-expodatabis_f[expodatabis_f$date=="30/09/20",]
rep4bis<-expodatabis_f[expodatabis_f$date=="01/10/20",]
rep5bis<-expodatabis_f[expodatabis_f$date=="07/10/20",]
rep6bis<-expodatabis_f[expodatabis_f$date=="08/10/20",]

#fitting the "null hypothesis model" rep1
expo_mod0<-drm(dead/total~dose,exposition,
               weights=total,
               data=rep1bis,
               fct=LN.2(),
               type="binomial")
summary(expo_mod0)
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,500),
     xt=c(0.01,0.1,1,10,100),legendPos=c(200,0.5),
     main="rep1 - 09/09/20 - dead")

#fitting the "null hypothesis model" rep2
expo_mod0<-drm(dead/total~dose,exposition,
               weights=total,
               data=rep2bis,
               fct=LN.2(),
               type="binomial")
summary(expo_mod0)
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,500),
     xt=c(0.01,0.1,1,10,100),legendPos=c(200,0.5),
     main="rep2 - 24/09/20 - dead")

#fitting the "null hypothesis model" rep3
expo_mod0<-drm(dead/total~dose,exposition,
               weights=total,
               data=rep3bis,
               fct=LN.2(),
               type="binomial")
summary(expo_mod0)
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,500),ylim=c(0,1),
     xt=c(0.001,0.01,0.1,1,10,100),legendPos=c(0.1,1),
     main="rep3 - 30/09/20 - dead")

#fitting the "null hypothesis model" rep4
expo_mod0<-drm(dead/total~dose,exposition,
               weights=total,
               data=rep4bis,
               fct=LN.2(),
               type="binomial")
summary(expo_mod0)
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,500),
     xt=c(0.001,0.01,0.1,1,10,100),legendPos=c(0.1,1),
     main="rep4 - 01/10/20 - dead")

#fitting the "null hypothesis model" rep4
expo_mod0<-drm(dead/total~dose,exposition,
               weights=total,
               data=rep5bis,
               fct=LN.2(),
               type="binomial")
summary(expo_mod0)
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,500),
     xt=c(0.001,0.01,0.1,1,10,100),legendPos=c(0.1,1),
     main="rep5 - 07/10/20 - dead")

#fitting the "null hypothesis model" rep4
expo_mod0<-drm(dead/total~dose,exposition,
               weights=total,
               data=rep6bis,
               fct=LN.2(),
               type="binomial")
summary(expo_mod0)
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,500),
     xt=c(0.001,0.01,0.1,1,10,100),legendPos=c(0.1,1),
     main="rep6 - 08/10/20 - dead")


##############################################################################/
#END
##############################################################################/