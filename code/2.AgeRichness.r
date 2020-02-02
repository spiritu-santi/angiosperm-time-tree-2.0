# Processes the family age file and a family species richness file to merge them into one table. 
# The species richness estimates per family and subfamily were obtained from the Angiosperm Phylogeny Website. 
# The code plots mean family ages against species richness.
# Continues with R objects generated in previous steps.
# Uses additional files.
richness_file <- "SPP.RICHNESS.csv" #### File in data folder

resB<-read.table(paste(ruta_write,"2.Ages_complete.csv",sep=""),sep=",",header=T)
rich<-read.table(richness_file,heade=T,sep=",")
rich_ord<-rich[which(rich$Family==""),]
resB$Order_richness<-rich_ord[match(resB$Order,rich_ord$Order),"Spp_richness"]
rich_fams<-rich[which(rich$Subfamily==""),]
rich_fams<-rich_fams[which(rich_fams$Family!=""),]
resB$Family_richness<-rich_fams[match(resB$Family,rich_fams$Family),"Spp_richness"]
rich_subfams<-rich[which(rich$Subfamily!=""),]
resB$Subfamily_richness<-rich_subfams[match(resB$Subfamily,rich_subfams$Subfamily),"Spp_richness"]
### Families with suspect crown nodes
non<-which(rownames(resB) %in% c( "Aphloiaceae","Brunelliaceae","Cynomoriaceae", 
                                  "Daphniphyllaceae", "Mayacaceae",
                                  "Mitrastemonaceae","Oncothecaceae"))
resB[non,"CG_age_Acal"]<-NA; resB[non,"CG_minHPD"]<-NA; resB[non,"CG_maxHPD"]<-NA
resB$Total_Richness <- rep(NA,dim(resB)[1])
for (i in 1:dim(resB)[1]){ 
  if(!is.na(resB$Subfamily_richness[i])){resB$Total_Richness[i] <- resB$Subfamily_richness[i];next}
  if(!is.na(resB$Family_richness[i])){resB$Total_Richness[i] <- resB$Family_richness[i];next}
  resB$Total_Richness[i] <- resB$Order_richness[i]}
out.group <- which(resB$Order%in%c("Pinales","Gnetales","Ginkgoales","Cycadales"))
resB <- resB[-out.group,]
orders <- 1:64
fams <- 65:500
subfams <- 501:dim(resB)[1]
resB$Fuse <- resB$SG_age_Acal - resB$CG_age_Acal
resB_fams <- resB[fams,]

sum(resB[fams,"Total_Richness"],na.rm = T)
par(mfrow = c(1, 1), pty = "s")
pdf(paste(ruta_write,"1.Temporal_distirbution_ages.pdf",sep=""))
ages <- resB_fams[which(!is.na(resB_fams[,c("CG_age_Acal")])),c("CG_age_Acal","SG_age_Acal")]
ages <- ages[order((ages$SG_age_Acal-ages$CG_age_Acal)),]
plot(ages$CG_age_Acal,cex=0.3,pch=19,ylim=c(0,150),type="n",xaxt="n",xlab="Families ranked by fuse period",
     ylab="Million years")
for (i in 1:length(ages$CG_age_Acal)){ 
  segments(x0 = i,x1 = i, y0 = 0, y1= ages$SG_age_Acal[i]-ages$CG_age_Acal[i],
           col=adjustcolor("black",alpha.f = 0.7),lwd=2)
}
title("Fuse period (stem to crown age)")
abline(h=mean(ages$SG_age_Acal-ages$CG_age_Acal),col="red",lwd=3,lty=3)
cat(mean(ages$SG_age_Acal-ages$CG_age_Acal),"\n")
hist(ages$SG_age_Acal-ages$CG_age_Acal,xlim=c(0,220),xlab="Fuse period (stem to crown age)",ylab="No. of families",main="",breaks=c(seq(0,220,10)))
ages[which(ages$SG_age_Acal-ages$CG_age_Acal > 120),]
ages[which(ages$SG_age_Acal-ages$CG_age_Acal < 10),]

full_SG <- c()
for (i in 1:nrow(resB_fams)){
  cat(i,"\r")
  if(is.na(resB_fams$SG_age_Acal[i])) next
  full_SG <- c(full_SG,seq(round(resB_fams$SG_maxHPD[i],1),round(resB_fams$SG_minHPD[i],0.5)))
}
cat("DONE ------")
full_CG <- c()
for (i in 1:nrow(resB_fams)){
  cat(i,"\r")
  if(is.na(resB_fams$CG_age_Acal[i])) next
  full_CG <- c(full_CG,seq(round(resB_fams$CG_maxHPD[i],1),round(resB_fams$CG_minHPD[i],0.5)))
}
cat("DONE ------")

hist(na.omit(resB_fams$SG_age_Acal),breaks=seq(0,250,10),freq=F,ylim=c(0,0.018),border="red",col=NULL,main="Family ages (95HPD)",xlab="million years ago")
hist(na.omit(resB_fams$CG_age_Acal),breaks=seq(0,250,10),freq=F,add=T,border="blue",col=NULL)
lines(density(full_SG,bw=5),col="red")
lines(density(full_CG,bw=5),col="blue")
dev.off()

periods <- data.frame("Stem (number)"=rep(NA,4),"Stem (percent)"=rep(NA,4),"Crown (number)"=rep(NA,4),"Crown (percent)"=rep(NA,4),
        row.names = c("Cenozoic","Cretaceous","Jurassic","Prior to first fossil"))
periods[1,3] <- length(which(resB_fams$CG_age_Acal < 66))
periods[1,4] <- round((length(which(resB_fams$CG_age_Acal < 66)) / length(which(!is.na(resB_fams$CG_age_Acal)))) * 100,1)
periods[2,3] <- length(which(resB_fams$CG_age_Acal >= 66 & resB_fams$CG_age_Acal < 145))
periods[2,4] <- round((length(which(resB_fams$CG_age_Acal >= 66 & resB_fams$CG_age_Acal < 145)) / length(which(!is.na(resB_fams$CG_age_Acal)))) * 100,1)
periods[3,3] <- length(which(resB_fams$CG_age_Acal >= 145))
periods[3,4] <- round((length(which(resB_fams$CG_age_Acal >= 145)) / length(which(!is.na(resB_fams$CG_age_Acal)))) * 100,1)
periods[4,3] <- length(which(resB_fams$CG_age_Acal >= 136))
periods[4,4] <- round((length(which(resB_fams$CG_age_Acal >= 136)) / length(which(!is.na(resB_fams$CG_age_Acal)))) * 100,1)
periods[1,1] <- length(which(resB_fams$SG_age_Acal < 66))
periods[1,2] <- round((length(which(resB_fams$SG_age_Acal < 66)) / length(which(!is.na(resB_fams$SG_age_Acal)))) * 100,1)
periods[2,1] <- length(which(resB_fams$SG_age_Acal >= 66 & resB_fams$SG_age_Acal < 145))
periods[2,2] <- round((length(which(resB_fams$SG_age_Acal >= 66 & resB_fams$SG_age_Acal < 145)) / length(which(!is.na(resB_fams$SG_age_Acal)))) * 100,1)
periods[3,1] <- length(which(resB_fams$SG_age_Acal >= 145))
periods[3,2] <- round((length(which(resB_fams$SG_age_Acal >= 145)) / length(which(!is.na(resB_fams$SG_age_Acal)))) * 100,1)
periods[4,1] <- length(which(resB_fams$SG_age_Acal >= 136))
periods[4,2] <- round((length(which(resB_fams$SG_age_Acal >= 136)) / length(which(!is.na(resB_fams$SG_age_Acal)))) * 100,1)
sink(paste(ruta_write,"3.Numbers_by_period.txt",sep=""))
periods
sink()

pdf(paste(ruta_write,"6.ARCs.pdf",sep=""),useDingbats = F)
par(mfrow = c(1, 1), pty = "s")
Data<-data.frame(Order=1:dim(resB_fams)[1],z=log(resB_fams$Total_Richness))
Data<-Data[order(Data$z),]
k=20
breaks <- getJenksBreaks(Data$z,k)
cols_breaks <- (colorRampPalette(c("yellow","gold","tomato2","red4"))(k))
Data$col <- cols_breaks[k]
for (i in k:1){Data$col[which(Data$z <= breaks[i])] <- cols_breaks[i]}
orderedcolors<-Data[order(Data$Order),"col"]
plot(resB_fams$Fuse,log(resB_fams$Total_Richness),pch=21,cex=log(resB_fams$Total_Richness)*0.5,xlab="Fuse period (My)"
     ,ylab="Species Richness (log)",bg=orderedcolors)
abline(lm(log(resB_fams$Total_Richness)~resB_fams$Fuse),lwd=3,lty=3,col="black")
summ_lm <- summary(lm(log(resB_fams$Total_Richness)~resB_fams$Fuse))
legend("bottomright",inset=0.03,legend=paste("b = ",round(summ_lm$coefficients[2,1],2),";","F = ",round(summ_lm$fstatistic[1],2),";",
                                             "p-val = ",round(summ_lm$coefficients[2,4],2),";","R2* = ",round(summ_lm$adj.r.squared,3),sep=""))
dev.off()
