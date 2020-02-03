# Estimates mean age of families and mean phylogenetic fuse of families across a global grid of 2x2 (resolution as in the original estimation of the presence absence matrix). 
# The estimation for every  grid-cell j is done by using the corresponding species richness (i.e., species present in j) for every family i present in j.
# (this procedure can also be performed using presence only (i.e., species richness = 1) of families in j, by transforming the PAM into a binary matrix).
# Continues with R objects generated in previous steps (saves and loads objects for plotting).
                                                      
which(rowSums(prPAM_g2[,-c(1:2)]) > 5) -> iden; prPAM_g2 <- prPAM_g2[iden,]
prPAM_gt <- prPAM_g2
for (i in 3:438){ 
  cat(i,"\n")
  prPAM_gt[which(prPAM_gt[,i]>=1),i] <- resB_fams[rownames(resB_fams)==colnames(prPAM_gt)[i],"CG_age_Acal"]
}
MDT_mat <- prPAM_g2; MDT_mat[,-c(1:2)] <- NA
richness_tot <- apply(prPAM_g2[,-c(1:2)], MARGIN = 1,sum,na.rm=T)
MDT_mat[,-c(1:2)] <- prPAM_gt[,-c(1:2)] * prPAM_g2[,-c(1:2)]
MDT_crown <- apply(MDT_mat[,-c(1:2)],MARGIN = 1, sum,na.rm=T)
MDT_crown <- MDT_crown /richness_tot

prPAM_gt <- prPAM_g2
for (i in 3:438){ 
  cat(i,"\n")
  prPAM_gt[which(prPAM_gt[,i]>=1),i] <- resB_fams[rownames(resB_fams)==colnames(prPAM_gt)[i],"SG_age_Acal"]
}
MDT_mat <- prPAM_g2; MDT_mat[,-c(1:2)] <- NA
richness_tot <- apply(prPAM_g2[,-c(1:2)], MARGIN = 1,sum,na.rm=T)
MDT_mat[,-c(1:2)] <- prPAM_gt[,-c(1:2)] * prPAM_g2[,-c(1:2)]
MDT_stem <- apply(MDT_mat[,-c(1:2)],MARGIN = 1, sum,na.rm=T)
MDT_stem <- MDT_stem /richness_tot

prPAM_gt <- prPAM_g2
for (i in 3:438){ 
  cat(i,"\n")
  prPAM_gt[which(prPAM_gt[,i]>=1),i] <- (resB_fams[rownames(resB_fams)==colnames(prPAM_gt)[i],"SG_age_Acal"] - resB_fams[rownames(resB_fams)==colnames(prPAM_gt)[i],"CG_age_Acal"])
}
MDT_mat <- prPAM_g2; MDT_mat[,-c(1:2)] <- NA
richness_tot <- apply(prPAM_g2[,-c(1:2)], MARGIN = 1,sum,na.rm=T)
MDT_mat[,-c(1:2)] <- prPAM_gt[,-c(1:2)] * prPAM_g2[,-c(1:2)]
MDT_fuse_non <- apply(MDT_mat[,-c(1:2)],MARGIN = 1, sum,na.rm=T)
MDT_fuse_non <- MDT_fuse_non /richness_tot
save(MDT_stem,MDT_crown,MDT_fuse_non,file=paste(ruta_write,"WeightedMAFS.RData",sep=""))
pdf(paste(ruta_write,"5.MAF_MAPS.pdf",sep=""))

par(mfrow = c(1, 1), pty = "s")
iden <- which(rowSums(prPAM_g2[,-c(1:2)]) > 5)
ssp<-coordinates(prPAM_g2[iden,1:2])
f <- clean_coordinates(as.data.frame(ssp), lon = "LON",lat = "LAT",
                       tests="seas",species=NULL,value="flagged",seas_scale = 110)
k <- 20
Data<-data.frame(Order=1:dim(prPAM_g2[iden,])[1],z=MDT_stem)
Data<-Data[order(Data$z),]
breaks <- getJenksBreaks(Data$z[f],k)
cols_breaks <- (colorRampPalette(c("linen","lightpink","red3","darkred"),bias=1.1)(k))
Data$col <- cols_breaks[k]
for (i in k:1){Data$col[which(Data$z <= breaks[i])] <- cols_breaks[i]}
orderedcolors<-Data[order(Data$Order),"col"]
maps::map("world",interior=F); title("Stem Mean Age of Families (MAF)")
points(ssp[f,],pch=22,lwd=0.4,col=orderedcolors[f],bg=orderedcolors[f],cex=0.5)
rect(xleft =-180,xright = 190,ybottom = -86,ytop = 85,border = "black")
abline(h=c(23.5,-23.5),lwd=2,lty=3)
color.legend(xl=-35,xr=35,yb=-83,yt=-67,align="lt",gradient="x",legend=c(round(range(MDT_stem),2)),
             rect.col=cols_breaks)
Data<-data.frame(Order=1:dim(prPAM_g2[iden,])[1],z=MDT_crown)
Data<-Data[order(Data$z),]
breaks <- getJenksBreaks(Data$z[f],k)
cols_breaks <- (colorRampPalette(c("aliceblue","lightblue","dodgerblue3","darkblue"),bias=1.1)(k))
#cols_breaks <- rev(inferno(k))
Data$col <- cols_breaks[k]
for (i in k:1){Data$col[which(Data$z <= breaks[i])] <- cols_breaks[i]}
orderedcolors<-Data[order(Data$Order),"col"]
maps::map("world",interior=F); title("Crown Mean Age of Families (MAF)")
points(ssp[f,],pch=22,lwd=0.3,col=orderedcolors[f],bg=orderedcolors[f],cex=0.5)
rect(xleft =-180,xright = 190,ybottom = -86,ytop = 85,border = "black")
abline(h=c(23.5,-23.5),lwd=2,lty=3)
color.legend(xl=-35,xr=35,yb=-83,yt=-67,align="lt",gradient="x",legend=c(round(range(MDT_crown),2)),
             rect.col=cols_breaks)
Data<-data.frame(Order=1:dim(prPAM_g2[iden,])[1],z=MDT_fuse_non)
Data<-Data[order(Data$z),];ssp <- coordinates(prPAM_g2[iden,1:2])
breaks <- getJenksBreaks(Data$z[f],k)
cols_breaks <- (colorRampPalette(c("honeydew","darkseagreen1","forestgreen","black"),bias=1)(k))
Data$col <- cols_breaks[k]
for (i in k:1){Data$col[which(Data$z <= breaks[i])] <- cols_breaks[i]}
orderedcolors<-Data[order(Data$Order),"col"]
maps::map("world",interior=F);title("Mean Fuse of Families (absolute)")
points(ssp[f,],pch=22,lwd=0.1,bg=orderedcolors[f],cex=0.6,col=orderedcolors[f])
rect(xleft =-180,xright = 190,ybottom = -86,ytop = 85,border = "black")
color.legend(xl=-35,xr=35,yb=-83,yt=-67,align="lt",gradient="x",
             legend=c(round(range(MDT_fuse_non),2)),rect.col=cols_breaks)
dev.off()
