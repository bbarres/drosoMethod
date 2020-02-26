##############################################################################/
##############################################################################/
#R code for the Drosophila suzukii pesticide resistance: Experiment 3
##############################################################################/
##############################################################################/

source("droso_data_load.R")


##############################################################################/
#Figure 3 combining the two radarplot####
##############################################################################/

op<-par(mfrow=c(1,2),mar=c(2,1,2,0))
radarchart(dataHS,
           #axis parameters
           axistype=1,caxislabels=seq(0.00,1.00,0.25),axislabcol="grey50",
           calcex=1.5,
           #grid parameters
           cglty=1,cglwd=2,cglcol="grey70",
           #labels parameters
           vlcex=1.5,
           #polygones parameters
           pcol=colpoly_bord,pfcol=colpoly_area,plwd=5,plty=1)
text(-1.2,1.1,labels=c("A"),cex=4)

radarchart(dataNa,
           #axis parameters
           axistype=1,caxislabels=seq(1,5,1),axislabcol="grey50",calcex=1.5,
           #grid parameters
           cglty=1,cglwd=2,cglcol="grey70",
           #labels parameters
           vlcex=1.5,
           #polygones parameters
           pcol=colpoly_bord,pfcol=colpoly_area,plwd=5,plty=1)
text(-1.2,1.1,labels=c("B"),cex=4)

par(op)
#export to pdf 14*7 inches


##############################################################################/
#Effect of the genetic diversity of the tested population on LD50
##############################################################################/

#we select the data of phosmet test with the St Foy & SF IsoA populations
genDdata<-dataDroz[dataDroz$genediv_comp==1,]
#aggregate the data
genDdata<-aggregate(cbind(dead,total)~dose+sex+population,data=genDdata,"sum")
genDdata_f<-genDdata[genDdata$sex=="female",]
genDdata_m<-genDdata[genDdata$sex=="male",]

#let's model the mortality rate for the females of both populations
genD_mod_f<-drm(dead/total~dose,weights=total,
                data=genDdata_f,curveid=population,
                fct=LN.3u(),
                type="binomial")
#testing the goodness-of-fit of the model
modelFit(genD_mod_f)
#comparing the LD50
compParm(genD_mod_f,"e")
sexrez_f<-ED(genD_mod_f,50,interval="delta",reference="control")

#let's model the mortality rate for the males of both populations
genD_mod_m<-drm(dead/total~dose,weights=total,
                data=genDdata_m,curveid=population,
                fct=LN.3u(),
                type="binomial")
#testing the goodness-of-fit of the model
modelFit(genD_mod_m)
#comparing the LD50
compParm(genD_mod_m,"e")
sexrez_m<-ED(genD_mod_m,50,interval="delta",reference="control")

#because there is a bug to display the 95CI with models using curveid,
#we plot the different modality separately
genD_mod_f_stf<-drm(dead/total~dose,weights=total,
                    data=genDdata_f[genDdata_f$population=="ste-foy",],
                    fct=LN.2(),
                    type="binomial")
genD_mod_f_isa<-drm(dead/total~dose,weights=total,
                    data=genDdata_f[genDdata_f$population=="sf-isoa",],
                    fct=LN.2(),
                    type="binomial")

genD_mod_m_stf<-drm(dead/total~dose,weights=total,
                    data=genDdata_m[genDdata_m$population=="ste-foy",],
                    fct=LN.2(),
                    type="binomial")
genD_mod_m_isa<-drm(dead/total~dose,weights=total,
                    data=genDdata_m[genDdata_m$population=="sf-isoa",],
                    fct=LN.2(),
                    type="binomial")


##############################################################################/
#Code for plotting the Figure 4####
##############################################################################/

#code for the plot comparing the different populations with different
#levels of genetic diversity, for male and female
op<-par(mar=c(0,5,6,1),mfrow=c(2,1))
#female plot for different genetic diversity populations
plot(genD_mod_f_stf,type="confidence",col=rgb(0.4,0.2,0.6,1),
     bty="n",axes=FALSE,ann=FALSE,lwd=3)
plot(genD_mod_f_stf,type="obs",add=TRUE,pch=21,cex=2,
     col=rgb(0.4,0.2,0.6,0.3),bg=rgb(0.4,0.2,0.6,0.3))
box(lwd=3,lty=1)
axis(1,at=c(1,10,50,150),labels=FALSE,
     cex.axis=1.5,font.axis=2,lwd.ticks=2)
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("0","20","40","60","80","100"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
# segments(ED(genD_mod_f_stf,50,interval="delta",reference="control")[1],
#          -0.2,
#          ED(genD_mod_f_stf,50,interval="delta",reference="control")[1],
#          0.5,lwd=3,col=rgb(0.4,0.2,0.6,1),lty=1)
plot(genD_mod_f_isa,type="confidence",add=TRUE,
     col=rgb(0.6,0.2,0.2,1),lwd=3)
plot(genD_mod_f_isa,type="obs",add=TRUE,pch=21,cex=2,
     col=rgb(0.6,0.2,0.2,0.3),bg=rgb(0.6,0.2,0.2,0.3))
# segments(ED(genD_mod_f_isa,50,interval="delta",reference="control")[1],
#          -0.2,
#          ED(genD_mod_f_isa,50,interval="delta",reference="control")[1],
#          0.5,lwd=3,col=rgb(0.6,0.2,0.2,1),lty=1)
# abline(h=0.5,lwd=3,col=grey(0.5))
text(1.5,y=0.85,labels='\\VE',vfont=c("sans serif","bold"),cex=5)
title(ylab="Mortality rate",cex.lab=2,font.lab=2)
#male plot for different genetic diversity populations
par(mar=c(5,5,1,1))
plot(genD_mod_m_stf,type="confidence",col=rgb(0.4,0.2,0.6,1),
     bty="n",axes=FALSE,ann=FALSE,lwd=3,lty=2)
plot(genD_mod_m_stf,type="obs",add=TRUE,pch=24,cex=2,
     col=rgb(0.4,0.2,0.6,1))
box(lwd=3,lty=1)
axis(1,at=c(1,10,50,150),labels=c("0","10","50","150"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2)
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("0","20","40","60","80","100"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
# segments(ED(genD_mod_m_stf,50,interval="delta",reference="control")[1],
#          -0.2,
#          ED(genD_mod_m_stf,50,interval="delta",reference="control")[1],
#          0.5,lwd=3,col=rgb(0.4,0.2,0.6,1),lty=2)
plot(genD_mod_m_isa,type="confidence",add=TRUE,
     col=rgb(0.6,0.2,0.2,1),lwd=3,lty=2)
plot(genD_mod_m_isa,type="obs",add=TRUE,pch=24,cex=2,
     col=rgb(0.6,0.2,0.2,1))
# segments(ED(genD_mod_m_isa,50,interval="delta",reference="control")[1],
#          -0.2,
#          ED(genD_mod_m_isa,50,interval="delta",reference="control")[1],
#          0.5,lwd=3,col=rgb(0.6,0.2,0.2,1),lty=2)
# abline(h=0.5,lwd=3,col=grey(0.5))
title(xlab="Dose (mg/l)",ylab="Mortality rate",cex.lab=2,font.lab=2)
text(1.5,y=0.85,labels='\\MA',vfont=c("sans serif","bold"),cex=5)
par(op)
#export to pdf 10*14 inches


##############################################################################/
#END
##############################################################################/