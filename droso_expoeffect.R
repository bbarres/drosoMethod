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
              fct=LN.2(),
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
                fct=LN.2(),
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
                fct=LN.2(),
                type="binomial")
plot(expo_f_mod,type="confidence")
plot(expo_f_mod,type="obs",add=TRUE)
temp<-EDcomp(expo_f_mod,c(50,50))
as.data.frame(temp)[,4]<0.05
exporez_f<-ED(expo_f_mod,50,interval="delta",reference="control")
plot(as.data.frame(exporez_f)$Estimate,ylim=c(0,0.3))

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
op<-par(mar=c(5.1,5.1,4.1,2.1))
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30),lwd=1.5,
     legendPos=c(15,0.7),xlab="dose (mg/L)",cex.axis=1.5,cex.lab=2,cex=2)
arrows(x0=expo_mod0$parmMat[2,1],y0=0.5,
       x1=expo_mod0$parmMat[2,4],y1=0.5,
       length=0.12,angle=25,lwd=3)
par(op)

#export to pdf 7 x 7 inches


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
#barplot of an example of evolution of the death rate at the dose 0.25mg/l####
##############################################################################/

data_expo<-read.table("data/droso_expo.txt",header=TRUE,sep="\t")
data_expo<-t(data_expo[,c(6:3,1)])
colnames(data_expo)<-data_expo[5,]

#the supplementary figure for the exposure duration
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
#heatmap to represent the significant differences between LD50####
##############################################################################/

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

chaudemap<-LDheatmap(temp,title=NULL,
                     add.map=FALSE,distances=NULL,SNP.name=row.names(temp),
                     color=c(rep(grey(0.8),3),
                             brewer.pal(6,"YlOrRd")[c(2,4,6)]),
                     name="CHR",flip=FALSE,add.key=FALSE)
grid.edit(gPath("CHR","heatMap","heatmap"),gp=gpar(col="white",lwd=1))
grid.edit(gPath("CHR","SNPnames"),gp=gpar(col="black",rot="0"),
          rot=0,hjust=0.7)
grid.lines(x=unit(c(0.1,0.5),"npc"),y=unit(c(0.5,0.5),"npc"),
           gp=gpar(lwd=3))
grid.lines(x=unit(c(0.5,0.5),"npc"),y=unit(c(0.5,0.9),"npc"),
           gp=gpar(lwd=3))

#export to pdf 7 x 7 inches


##############################################################################/
#END
##############################################################################/



##############################################################################/
#Final plot examplifying effect of time of exposure
##############################################################################/

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
title(main="dose = 0.25 mg/L",xlab=NULL,ylab="Number of flies",cex.lab=2,
      line=2,font.lab=2,cex.main=2)
text(expobarplot-0.03,as.numeric(data_expo[1,c(12:22)])/2,
     data_expo[1,c(12:22)],font=2,cex=2,xpd=TRUE,col="white")
text(expobarplot-0.03,as.numeric(data_expo[1,c(12:22)]) + 
        as.numeric(data_expo[2,c(12:22)])/2,
     data_expo[2,c(12:22)],font=2,cex=2,xpd=TRUE,col="black")
text(expobarplot-0.03,as.numeric(data_expo[1,c(12:22)]) + 
        as.numeric(data_expo[2,c(12:22)]) + 
        as.numeric(data_expo[3,c(12:22)])/2,
     data_expo[3,c(12:22)],font=2,cex=2,xpd=TRUE,col="black")

#export to pdf 7 x 7 inches


#comparison of the regression curves for the different reading time
op<-par(mar=c(5.1,5.1,4.1,2.1))
plot(expo_mod0,col=c(1,1,1,1,1,2,2,2,2,2),xlim=c(0,30),lwd=1.5,
     legendPos=c(15,0.7),xlab="dose (mg/L)",cex.axis=1.5,cex.lab=2,
     cex=2,axes=TRUE)
arrows(x0=expo_mod0$parmMat[2,1],y0=0.5,
       x1=expo_mod0$parmMat[2,4],y1=0.5,
       length=0.12,angle=25,lwd=3)
box(bty="o",lwd=4)
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

chaudemap<-LDheatmap(temp,title=NULL,
                     add.map=FALSE,distances=NULL,SNP.name=row.names(temp),
                     color=c(rep(grey(0.8),3),
                             brewer.pal(6,"YlOrRd")[c(2,4,6)]),
                     name="CHR",flip=FALSE,add.key=FALSE)
grid.edit(gPath("CHR","heatMap","heatmap"),gp=gpar(col="white",lwd=1))
grid.edit(gPath("CHR","SNPnames"),gp=gpar(col="black",rot="0"),
          rot=0,hjust=0.7)
grid.lines(x=unit(c(0.1,0.5),"npc"),y=unit(c(0.5,0.5),"npc"),
           gp=gpar(lwd=3))
grid.lines(x=unit(c(0.5,0.5),"npc"),y=unit(c(0.5,0.9),"npc"),
           gp=gpar(lwd=3))

#export to pdf 7 x 7 inches



