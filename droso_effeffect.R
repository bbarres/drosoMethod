###############################################################################
###############################################################################
#R code for the Drosophila suzukii pesticide resistance tests paper
###############################################################################
###############################################################################


###############################################################################
#What is the effect of number of flies used for a test on LD50 evaluation?
###############################################################################

source("droso_data_load.R")

#we select the data of phosmet test with the St Foy population
numberdata<-dataDroz[dataDroz$number_comp==1,]
#because on the 2016-06-23, the behaviour of the flies was abnormal, 
#this repetition was removed before analysis
numberdata<-numberdata[numberdata$date!="2016-06-23",]

#first we remove unnecessary levels in the data frame
numberdata<-drop.levels(numberdata)
#we then create a dataframe for the results
REZdroz<-data.frame("repet"=as.character(),"ED50"=as.numeric(),
                    "IC_low"=as.numeric(),"IC_up"=as.numeric(),
                    "SE"=as.numeric())
#and here comes the loop
for (i in 1: length(levels(numberdata$repet))) {
  temp.m1<-drm(dead/total~dose,weights=total,
               data=numberdata[numberdata$repet==levels(numberdata$repet)[i],],
               fct=LN.2(),
               type="binomial")
  temp<-ED(temp.m1,50,interval="delta",reference="control")
  tempx<-data.frame("date"=names(table(numberdata$repet))[i],
                    "ED50"=temp[1],"IC_low"=temp[3],"IC_up"=temp[4],
                    "SE"=temp[2])
  REZdroz<-rbind(REZdroz,tempx)
}
REZdroz<-REZdroz[order(as.character(REZdroz$date)),]
results<- merge(checkdat,REZdroz,by.x="repet",by.y="date",all=FALSE)


###############################################################################
#the scatter plot of the LD50 analysis with different number of fly per dose
###############################################################################

legx<-expression(bold("Mean number of ")*bolditalic("D. suzukii ")*
                   bold("per dose"))
op<-par(mar=c(5,5,1,1),mfrow=c(1,2))
plot(results$ED50[results$sex=="female"]~results$total[results$sex=="female"],
     xlab ="Mean number of D. suzukii per dose",ylab="LD50 (mg/l)",
     ylim=c(0,100),xlim=c(50,270),bty="n",ann=FALSE,axes=FALSE,
     cex=2,pch=21,col=rgb(0,0,0,0.0),bg=rgb(0,0,0,0.0))
box(lwd=3,lty=1)
axis(1,at=c(70,105,140,175,210,245),labels=c("10","15","20","25","30","35"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2)
axis(2,at=c(0,20,40,60,80,100),labels=c("0","20","40","60","80","100"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
title(xlab=legx,
      ylab="LD50 (mg/l)",cex.lab=2,font.lab=2)
text(230,y=90,labels='\\VE',vfont=c("sans serif","bold"),cex=5)
plotCI(results$total[results$sex=="female"],
       results$ED50[results$sex=="female"],
       ui=results$ED50[results$sex=="female"]+
         results$SE[results$sex=="female"],
       li=results$ED50[results$sex=="female"]-
         results$SE[results$sex=="female"],
       #ui=results$IC_up[results$sex=="female"],
       #li=results$IC_low[results$sex=="female"],
       add=TRUE,cex=2,pch=21,col=rgb(0,0,0,1),pt.bg=rgb(0,0,0,0.3),
       gap=0.014)
totfem<-drm(dead/total~dose,weights=total,
            data=numberdata[numberdata$sex=="female",],
            fct=LN.2(),
            type="binomial")
totfemREZ<-ED(totfem,50,interval="delta",reference="control")
abline(totfemREZ[1],0,col="red",lwd=2)
abline(totfemREZ[3],0,col="red",lwd=2,lty=2)
abline(totfemREZ[4],0,col="red",lwd=2,lty=2)

par(mar=c(5,2,1,4))
plot(results$ED50[results$sex=="male"]~results$total[results$sex=="male"],
     xlab ="Mean number of D. suzukii per dose",ylab="LD50 (mg/l)",
     ylim=c(0,100),xlim=c(50,270),bty="n",ann=FALSE,axes=FALSE,
     cex=2,pch=24,col=rgb(0,0,0,0.0),bg=rgb(0,0,0,0.0))
box(lwd=3,lty=1)
axis(1,at=c(70,105,140,175,210,245),labels=c("10","15","20","25","30","35"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2)
axis(2,at=c(0,20,40,60,80,100),labels=FALSE,
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
title(xlab=legx,
      ylab="",cex.lab=2,font.lab=2)
text(230,y=90,labels='\\MA',vfont=c("sans serif","bold"),cex=5)
plotCI(results$total[results$sex=="male"],
       results$ED50[results$sex=="male"],
       ui=results$ED50[results$sex=="male"]+
         results$SE[results$sex=="male"],
       li=results$ED50[results$sex=="male"]-
         results$SE[results$sex=="male"],
       #ui=results$IC_up[results$sex=="male"],
       #li=results$IC_low[results$sex=="male"],
       add=TRUE,cex=2,pch=24,col=rgb(0,0,0,1),pt.bg=rgb(0,0,0,0.0),
       gap=0.014)
totmal<-drm(dead/total~dose,weights=total,
            data=numberdata[numberdata$sex=="male",],
            fct=LN.2(),
            type="binomial")
totmalREZ<-ED(totmal,50,interval="delta",reference="control")
abline(totmalREZ[1],0,col="red",lwd=2)
abline(totmalREZ[3],0,col="red",lwd=2,lty=2)
abline(totmalREZ[4],0,col="red",lwd=2,lty=2)
par(op)
#export .pdf 15*7 inches


###############################################################################
#END
###############################################################################