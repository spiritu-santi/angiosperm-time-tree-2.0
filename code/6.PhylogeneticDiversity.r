# Estimates phylogentic species variability (PSV) across a global grid of 2x2 (resolution as in the original estimation of the presence absence matrix). 
# The code also estimates spatial null models for PSV by shuffling the PAM in two ways:
# - NULL 1 randomly shuffles occurrences within sites (i.e., per-row shuffling), keeping the observed species richness.
# - NULL 2 randomly shuffles occurrences within families (i.e., per-column shuffling), keeping the observed family prevalence.
# PSV is estimated on a dated phylogeny that is pruned to includes a sinlge tip per family.
# Continues with R objects generated in previous steps.

tre.pruned <- drop.tip(tre,tip=list[which(duplicated(list_fams))])
tre.pruned<- ladderize(tre.pruned, FALSE)
list <- tre.pruned$tip.label
list_outs <- foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][2]}
outs<-c("Welwitschiaceae","Pinaceae","Ginkgoaceae","Cupressaceae","Gnetaceae")
tre.pruned <- drop.tip(tre.pruned,tip=tre.pruned$tip.label[which(list_outs%in%outs)])
list<-tre.pruned$tip.label
list <- foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][2]}
bb <- resB_fams[match(list,rownames(resB_fams)),]
match(list,rownames(bb))
print(bb[1:5,];tre.pruned$tip.label[1:5])

list <- tre.pruned$tip.label
list_fams <- foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][2]}
tre.pruned$tip.label<-list_fams
is.rooted(tre.pruned); is.ultrametric(tre.pruned)
tre.pruned

iden <- which(rowSums(prPAM_g2[,-c(1:2)]) > 5)
psv <- psv(samp = prPAM_g2[iden,-c(1:2)],tree = tre.pruned)
PSV_random_freq <- matrix(nrow = dim(prPAM_g2[iden,])[1], ncol = 1002)
for (i in 3:ncol(PSV_random_freq)){ 
  cat(i, "\n")
  psv_random <- psv(samp = randomizeMatrix(prPAM_g2[iden,-c(1:2)],"freq",10000000),tree = tre.pruned)
  PSV_random_freq[,i] <- psv_random$PSVs
}
save(PSV_random_freq,file=paste(ruta_write,"PSV_random_freq.Rdata",sep=""))
PSV_random_rich <- matrix(nrow = dim(prPAM_g2[iden,])[1], ncol = 1002)
for (i in 1:ncol(PSV_random_rich)){ 
  cat(i, "\n")
  psv_random <- psv(samp = randomizeMatrix(prPAM_g2[iden,-c(1:2)],"rich",10000000),tree = tre.pruned)
  PSV_random_rich[,i] <- psv_random$PSVs
}
save(PSV_random_rich,file=paste(ruta_write,"PSV_random_rich.Rdata",sep=""))

pdf(paste(ruta_write,"10.PSV_randonm.pdf",sep=""))
ssp <- coordinates(prPAM_g2[iden,1:2])
Data<-data.frame(Order=1:dim(prPAM_g2[iden,])[1],z=psv$PSVs[iden])
Data<-Data[order(Data$z),]
f <- clean_coordinates(as.data.frame(ssp), lon = "LON",lat = "LAT",
                       tests="seas",species=NULL,value="flagged",seas_scale=110)
k=20
breaks <- getJenksBreaks(Data$z,k)
cols_breaks <- rev(viridis(k,begin=0.0,end=1))
Data$col <- cols_breaks[k]
for (i in k:1){ 
Data$col[which(Data$z <= breaks[i])] <- cols_breaks[i]}
orderedcolors<-Data[order(Data$Order),"col"]
maps::map("world",interior=F)
points(ssp[f,],pch=22,lwd=0.1,bg=orderedcolors[f],cex=0.6,col=NULL)
rect(xleft =-180,xright = 190,ybottom = -86,ytop = 85,border = "black")
color.legend(xl=-35,xr=35,yb=-83,yt=-67,align="lt",gradient="x",
             legend=round(range(psv$PSVs,na.rm = T),2),rect.col=cols_breaks)
abline(h=c(23.5,-23.5),lwd=2,lty=3)
title("Phylogenetic Species Variability")

SES_psv <- (psv$PSVs[iden] - apply(PSV_random_rich[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(PSV_random_rich[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
Data<-data.frame(Order=1:dim(prPAM_g2[iden,])[1],z=SES_psv)
Data<-Data[order(Data$z),]
k=20
breaks <- getJenksBreaks(Data$z,k)
cols_breaks <- rev(viridis(k))
Data$col <- cols_breaks[k]
for (i in k:1){ 
  Data$col[which(Data$z <= breaks[i])] <- cols_breaks[i]}
orderedcolors<-Data[order(Data$Order),"col"]
maps::map("world",interior=F)
points(ssp[f,],pch=22,lwd=0.1,bg=orderedcolors[f],cex=0.5,col=NULL)
rect(xleft =-180,xright = 190,ybottom = -86,ytop = 85,border = "black")
color.legend(xl=-35,xr=35,yb=-83,yt=-67,align="lt",gradient="x",
             legend=round(range(SES_psv,na.rm = T),2),rect.col=cols_breaks)
title("SES-Phylogenetic Species Variability (NULL1)")

SES_psv <- (psv$PSVs[iden] - apply(PSV_random_freq[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(PSV_random_freq[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
Data<-data.frame(Order=1:dim(prPAM_g2[iden,])[1],z=SES_psv)
Data<-Data[order(Data$z),]
k=20
breaks <- getJenksBreaks(Data$z,k)
cols_breaks <- rev(viridis(k))
Data$col <- cols_breaks[k]
for (i in k:1){ 
  Data$col[which(Data$z <= breaks[i])] <- cols_breaks[i]}
orderedcolors<-Data[order(Data$Order),"col"]
maps::map("world",interior=F)
points(ssp[f,],pch=22,lwd=0.1,bg=orderedcolors[f],cex=0.5,col=NULL)
rect(xleft =-180,xright = 190,ybottom = -86,ytop = 85,border = "black")
color.legend(xl=-35,xr=35,yb=-83,yt=-67,align="lt",gradient="x",
             legend=round(range(SES_psv,na.rm = T),2),rect.col=cols_breaks)
title("SES-Phylogenetic Species Variability (NULL2)")
dev.off()
