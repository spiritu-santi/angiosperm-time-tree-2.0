resSTRICT<-read.table("RES_SBphylo_v2/2.Ages_complete.csv",sep=",",header=T)
resRAW<-read.table("RES_SBcomplete_v2/2.Ages_complete.csv",sep=",",header=T)
pdf("1.ComparisonFossils_RC.pdf")

rich<-read.table("SPP.RICHNESS.csv",heade=T,sep=",")
rich_ord<-rich[which(rich$Family==""),]; dim(rich_ord)
resSTRICT$Order_richness<-rich_ord[match(resSTRICT$Order,rich_ord$Order),"Spp_richness"]
rich_fams<-rich[which(rich$Subfamily==""),]; dim(rich_fams)
rich_fams<-rich_fams[which(rich_fams$Family!=""),]; dim(rich_fams)
resSTRICT$Family_richness<-rich_fams[match(resSTRICT$Family,rich_fams$Family),"Spp_richness"]
rich_subfams<-rich[which(rich$Subfamily!=""),]; dim(rich_subfams)
resSTRICT$Subfamily_richness<-rich_subfams[match(resSTRICT$Subfamily,rich_subfams$Subfamily),"Spp_richness"]
non<-which(rownames(resSTRICT) %in% c("Aphloiaceae", "Daphniphyllaceae", 
                                      "Brunelliaceae", "Cynomoriaceae", "Mayacaceae", 
                                      "Mitrastemonaceae", "Oncothecaceae"))
resSTRICT[non,"CG_age_Acal"]<-NA
resSTRICT$Total_Richness <- rep(NA,dim(resSTRICT)[1])
for (i in 1:dim(resSTRICT)[1]){ 
  if(!is.na(resSTRICT$Subfamily_richness[i])){resSTRICT$Total_Richness[i] <- resSTRICT$Subfamily_richness[i];next}
  if(!is.na(resSTRICT$Family_richness[i])){resSTRICT$Total_Richness[i] <- resSTRICT$Family_richness[i];next}
  resSTRICT$Total_Richness[i] <- resSTRICT$Order_richness[i]}
out.group<-which(resSTRICT$Order%in%c("Pinales","Gnetales","Ginkgoales","Cycadales"))
resSTRICT<-resSTRICT[-out.group,]
resRAW$Order_richness<-rich_ord[match(resRAW$Order,rich_ord$Order),"Spp_richness"]
rich_fams<-rich[which(rich$Subfamily==""),]; dim(rich_fams)
rich_fams<-rich_fams[which(rich_fams$Family!=""),]; dim(rich_fams)
resRAW$Family_richness<-rich_fams[match(resRAW$Family,rich_fams$Family),"Spp_richness"]
rich_subfams<-rich[which(rich$Subfamily!=""),]; dim(rich_subfams)
resRAW$Subfamily_richness<-rich_subfams[match(resRAW$Subfamily,rich_subfams$Subfamily),"Spp_richness"]
non<-which(rownames(resRAW) %in% c("Aphloiaceae", "Daphniphyllaceae", 
                                   "Brunelliaceae", "Cynomoriaceae", "Mayacaceae", 
                                   "Mitrastemonaceae", "Oncothecaceae"))
resRAW[non,"CG_age_Acal"]<-NA
resRAW$Total_Richness <- rep(NA,dim(resRAW)[1])
for (i in 1:dim(resRAW)[1]){ 
  if(!is.na(resRAW$Subfamily_richness[i])){resRAW$Total_Richness[i] <- resRAW$Subfamily_richness[i];next}
  if(!is.na(resRAW$Family_richness[i])){resRAW$Total_Richness[i] <- resRAW$Family_richness[i];next}
  resRAW$Total_Richness[i] <- resRAW$Order_richness[i]}
out.group<-which(resRAW$Order%in%c("Pinales","Gnetales","Ginkgoales","Cycadales"))
resRAW<-resRAW[-out.group,]

orders<-1:64
fams<-65:500
subfams<-501:dim(resRAW)[1]
unique(rownames(resRAW)==rownames(resSTRICT))
resRAW<-resRAW[fams,]
resSTRICT<-resSTRICT[fams,]

plot(resRAW$SG_age_Acal,resSTRICT$SG_age_Acal,type="n",ylim=c(0,200),xlim=c(0,200),
     xlab="Stem ages (complete set)",ylab="Stem ages (phylo-only set)");abline(a=0,b=1)
title("Stem age estimates (families)")
for (i in 1:length(resSTRICT$SG_age_Acal)) segments(y0=resRAW$SG_minHPD[i],y1=resRAW$SG_maxHPD[i],
                                                          x0=resSTRICT$SG_age_Acal[i],x1=resSTRICT$SG_age_Acal[i],lwd = 0.5,col="grey40")
for (i in 1:length(resSTRICT$SG_age_Acal)) segments(x0=resSTRICT$SG_minHPD[i],x1=resSTRICT$SG_maxHPD[i],
                                                          y0=resRAW$SG_age_Acal[i],y1=resRAW$SG_age_Acal[i],lwd = 0.5,col="grey40")
points(resSTRICT$SG_age_Acal,resRAW$SG_age_Acal,pch=21,cex=1,bg=rgb(1,0,0,0.8))
plot(resSTRICT$CG_age_Acal,resRAW$CG_age_Acal,type="n",ylim=c(0,200),xlim=c(0,200),
     xlab="Crown ages (complete set)",ylab="Crown ages (phylo-only set)");abline(a=0,b=1)
title("Crown age estimates (families)")
for (i in 1:length(resSTRICT$SG_age_Acal)) segments(y0=resRAW$CG_minHPD[i],y1=resRAW$CG_maxHPD[i],
                                                          x0=resSTRICT$CG_age_Acal[i],x1=resSTRICT$CG_age_Acal[i],lwd = 0.5,col="grey40")
for (i in 1:length(resSTRICT$SG_age_Acal)) segments(x0=resSTRICT$CG_minHPD[i],x1=resSTRICT$CG_maxHPD[i],
                                                          y0=resRAW$CG_age_Acal[i],y1=resRAW$CG_age_Acal[i],lwd = 0.5,col="grey40")
points(resSTRICT$CG_age_Acal,resRAW$CG_age_Acal,pch=21,cex=1,bg=rgb(1,0,0,0.8))
HPD_complete<-resRAW$CG_maxHPD-resRAW$CG_minHPD
HPD_strict<-resSTRICT$CG_maxHPD-resSTRICT$CG_minHPD
HPDsC<-cbind(HPD_complete,HPD_strict)
boxplot(HPDsC,pch=19,cex=0.6,col=c("grey80","grey30"),pch=19,cex=0.5,
        ylab="95% HPD crown ages")
title("HPD crown age estimates")

HPD_complete<-resRAW$SG_maxHPD-resRAW$SG_minHPD
HPD_strict<-resSTRICT$SG_maxHPD-resSTRICT$SG_minHPD
HPDsC<-cbind(HPD_complete,HPD_strict)
boxplot(HPDsC,pch=19,cex=0.6,col=c("grey80","grey30"),pch=19,cex=0.5,
        ylab="95% HPD stem ages")
title("HPD stem age estimates")
dev.off()
