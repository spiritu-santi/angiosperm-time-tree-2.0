# Estimate spatial null models for mean family age and fuse across a global grid of 2x2 (resolution as in the original estimation of the presence absence matrix). 
# The estimation of the null models is performed by shuffling entries in the PAM in two ways:
# - NULL 1 randomly shuffles occurrences within sites (i.e., per-row shuffling), keeping the observed species richness.
# - NULL 2 randomly shuffles occurrences within families (i.e., per-column shuffling), keeping the observed family prevalence.
# Continues with R objects generated in previous steps.

MDT_crown_NULL2 <- matrix(nrow = dim(prPAM_g2)[1], ncol = 1002)
for (j in 3:ncol(MDT_crown_NULL2)){ 
  cat("SIMULATION NULL2",j-2,"\r")
  prPAM_random <- prPAM_g2
  crown <- which(!is.na(resB_fams$CG_age_Acal))
  names_crown <- rownames(resB_fams)[crown]
  mm <- match(names_crown,colnames(prPAM_random))
  prPAM_random[,mm] <- randomizeMatrix(prPAM_random[,mm],null.model="freq")
  prPAM_random_gt <- prPAM_random
  for (i in 3:438){ 
    prPAM_random_gt[which(prPAM_random_gt[,i]>=1),i] <- resB_fams[rownames(resB_fams)==colnames(prPAM_random_gt)[i],"CG_age_Acal"]
  }
  richness_tot <- apply(prPAM_random[,-c(1:2)], MARGIN = 1,sum,na.rm=T)
  MDT_mat <- prPAM_random
  MDT_mat[,-c(1:2)] <- NA
  MDT_mat[,-c(1:2)] <- prPAM_random_gt[,-c(1:2)] * prPAM_random[,-c(1:2)]
  MDT_crown_random <- apply(MDT_mat[,-c(1:2)],MARGIN = 1, sum,na.rm=T) /richness_tot
  MDT_crown_NULL2[,j] <- MDT_crown_random
}
save(MDT_crown_NULL2,file=paste(ruta_write,"MDT_crown_NULL2.Rdata"))

MDT_crown_NULL1 <- matrix(nrow = dim(prPAM_g2)[1], ncol = 1002)
for (j in 3:ncol(MDT_crown_NULL1)){ 
  cat("SIMULATION NULL1:",j-2,"\r")
  prPAM_random <- prPAM_g2
  crown <- which(!is.na(resB_fams$CG_age_Acal))
  names_crown <- rownames(resB_fams)[crown]
  mm <- match(names_crown,colnames(prPAM_random))
  prPAM_random[,mm] <- randomizeMatrix(prPAM_random[,mm],null.model="richness")
  prPAM_random_gt <- prPAM_random
  for (i in 3:438){ 
    prPAM_random_gt[which(prPAM_random_gt[,i]>=1),i] <- resB_fams[rownames(resB_fams)==colnames(prPAM_random_gt)[i],"CG_age_Acal"]
  }
  richness_tot <- apply(prPAM_random[,-c(1:2)], MARGIN = 1,sum,na.rm=T)
  MDT_mat <- prPAM_random
  MDT_mat[,-c(1:2)] <- NA
  MDT_mat[,-c(1:2)] <- prPAM_random_gt[,-c(1:2)] * prPAM_random[,-c(1:2)]
  MDT_crown_random <- apply(MDT_mat[,-c(1:2)],MARGIN = 1, sum,na.rm=T) /richness_tot
  MDT_crown_NULL1[,j] <- MDT_crown_random
}
save(MDT_crown_NULL1,file=paste(ruta_write,"MDT_crown_NULL1.Rdata"))


MDT_stem_NULL2 <- matrix(nrow = dim(prPAM_g2)[1], ncol = 1002)
for (j in 3:ncol(MDT_stem_NULL2)){ 
  cat("SIMULATION NULL2",j-2,"\r")
  prPAM_random <- prPAM_g2
  back_up<-colnames(prPAM_random)[-c(1:2)]
  crown <- which(!is.na(resB_fams$SG_age_Acal))
  names_crown <- rownames(resB_fams)[crown]
  mm <- match(names_crown,colnames(prPAM_random))
  prPAM_random[,mm] <- randomizeMatrix(prPAM_random[,mm],
                                       null.model="freq")
  prPAM_random_gt <- prPAM_random
  for (i in 3:438){ 
    prPAM_random_gt[which(prPAM_random_gt[,i]>=1),i] <- resB_fams[rownames(resB_fams)==colnames(prPAM_random_gt)[i],"SG_age_Acal"]
  }
  back_up->colnames(prPAM_random_gt)[-c(1:2)]
  richness_tot <- apply(prPAM_random[,-c(1:2)], MARGIN = 1,sum,na.rm=T)
  MDT_mat <- prPAM_random
  MDT_mat[,-c(1:2)] <- NA
  MDT_mat[,-c(1:2)] <- prPAM_random_gt[,-c(1:2)] * prPAM_random[,-c(1:2)]
  MDT_stem_random <- apply(MDT_mat[,-c(1:2)],MARGIN = 1, sum,na.rm=T) /richness_tot
  MDT_stem_NULL2[,j] <- MDT_stem_random
}
save(MDT_stem_NULL2,file=paste(ruta_write,"MDT_stem_NULL2.Rdata"))

MDT_stem_NULL1 <- matrix(nrow = dim(prPAM_g2)[1], ncol = 1002)
for (j in 3:ncol(MDT_stem_NULL1)){ 
  cat("SIMULATION NULL1:",j-2,"\r")
  prPAM_random <- prPAM_g2
  back_up<-colnames(prPAM_random)[-c(1:2)]
  crown <- which(!is.na(resB_fams$SG_age_Acal))
  names_crown <- rownames(resB_fams)[crown]
  mm <- match(names_crown,colnames(prPAM_random))
  prPAM_random[,mm] <- randomizeMatrix(prPAM_random[,mm],
                                       null.model="richness")
  prPAM_random_gt <- prPAM_random
  for (i in 3:438){ 
    prPAM_random_gt[which(prPAM_random_gt[,i]>=1),i] <- resB_fams[rownames(resB_fams)==colnames(prPAM_random_gt)[i],"SG_age_Acal"]
  }
  back_up->colnames(prPAM_random_gt)[-c(1:2)]
  richness_tot <- apply(prPAM_random[,-c(1:2)], MARGIN = 1,sum,na.rm=T)
  MDT_mat <- prPAM_random
  MDT_mat[,-c(1:2)] <- NA
  MDT_mat[,-c(1:2)] <- prPAM_random_gt[,-c(1:2)] * prPAM_random[,-c(1:2)]
  MDT_stem_random <- apply(MDT_mat[,-c(1:2)],MARGIN = 1, sum,na.rm=T) /richness_tot
  MDT_stem_NULL1[,j] <- MDT_stem_random
}
save(MDT_stem_NULL1,file=paste(ruta_write,"MDT_stem_NULL1.Rdata"))

MDT_fuse_non_NULL2 <- matrix(nrow = dim(prPAM_g2)[1], ncol = 1002)
for (j in 3:ncol(MDT_fuse_non_NULL2)){ 
  cat("SIMULATION NULL2",j-2,"\r")
  prPAM_random <- prPAM_g2
  back_up<-colnames(prPAM_random)[-c(1:2)]
  crown <- which(!is.na(resB_fams$SG_age_Acal))
  names_crown <- rownames(resB_fams)[crown]
  mm <- match(names_crown,colnames(prPAM_random))
  prPAM_random[,mm] <- randomizeMatrix(prPAM_random[,mm],
                                       null.model="freq")
  prPAM_random_gt <- prPAM_random
  for (i in 3:438){ 
    prPAM_random_gt[which(prPAM_random_gt[,i]>=1),i] <- (resB_fams[rownames(resB_fams)==colnames(prPAM_random_gt)[i],"SG_age_Acal"]-resB_fams[rownames(resB_fams)==colnames(prPAM_random_gt)[i],"CG_age_Acal"])
  }
  back_up->colnames(prPAM_random_gt)[-c(1:2)]
  richness_tot <- apply(prPAM_random[,-c(1:2)], MARGIN = 1,sum,na.rm=T)
  MDT_mat <- prPAM_random
  MDT_mat[,-c(1:2)] <- NA
  MDT_mat[,-c(1:2)] <- prPAM_random_gt[,-c(1:2)] * prPAM_random[,-c(1:2)]
  MDT_fuse_random <- apply(MDT_mat[,-c(1:2)],MARGIN = 1, sum,na.rm=T) /richness_tot
  MDT_fuse_non_NULL2[,j] <- MDT_fuse_random
}
save(MDT_fuse_non_NULL2,file=paste(ruta_write,"MDT_fuse_non_NULL2.Rdata"))
MDT_fuse_non_NULL1 <- matrix(nrow = dim(prPAM_g2)[1], ncol = 1002)
for (j in 3:ncol(MDT_fuse_non_NULL1)){ 
  cat("SIMULATION NULL1:",j-2,"\r")
  prPAM_random <- prPAM_g2
  back_up<-colnames(prPAM_random)[-c(1:2)]
  crown <- which(!is.na(resB_fams$SG_age_Acal))
  names_crown <- rownames(resB_fams)[crown]
  mm <- match(names_crown,colnames(prPAM_random))
  prPAM_random[,mm] <- randomizeMatrix(prPAM_random[,mm],
                                       null.model="richness")
  prPAM_random_gt <- prPAM_random
  for (i in 3:438){ 
    prPAM_random_gt[which(prPAM_random_gt[,i]>=1),i] <- (resB_fams[rownames(resB_fams)==colnames(prPAM_random_gt)[i],"SG_age_Acal"]-resB_fams[rownames(resB_fams)==colnames(prPAM_random_gt)[i],"CG_age_Acal"])
  }
  back_up->colnames(prPAM_random_gt)[-c(1:2)]
  richness_tot <- apply(prPAM_random[,-c(1:2)], MARGIN = 1,sum,na.rm=T)
  MDT_mat <- prPAM_random
  MDT_mat[,-c(1:2)] <- NA
  MDT_mat[,-c(1:2)] <- prPAM_random_gt[,-c(1:2)] * prPAM_random[,-c(1:2)]
  MDT_fuse_random <- apply(MDT_mat[,-c(1:2)],MARGIN = 1, sum,na.rm=T) /richness_tot
  MDT_fuse_non_NULL1[,j] <- MDT_fuse_random
}
save(MDT_fuse_non_NULL1,file=paste(ruta_write,"MDT_fuse_non_NULL1.Rdata"))

pdf(paste(ruta_write,"SES-MAF_NULL2.pdf",sep=""))
par(mfrow = c(1, 1), pty = "s")
null_crown <- MDT_crown_NULL2
null_stem <- MDT_stem_NULL2
null_fuse <- MDT_fuse_non_NULL2
iden <- which(rowSums(prPAM_g2[,-c(1:2)]) > 5)
ssp<-coordinates(prPAM_g2[iden,1:2])
SES_crownMAF <- (MDT_crown - apply(null_crown[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_crown[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
Data<-data.frame(Order=1:dim(prPAM_g2[iden,])[1],z=SES_crownMAF)
Data<-Data[order(Data$z),]
f <- clean_coordinates(as.data.frame(ssp), lon = "LON",lat = "LAT",
                       tests="seas",species=NULL,value="flagged",seas_scale=110)
k=20
breaks <- getJenksBreaks(Data$z,k)
cols_breaks <- (colorRampPalette(c("aliceblue","lightblue","dodgerblue3","darkblue"),bias=1.1)(k))
Data$col <- cols_breaks[k]
for (i in k:1){ 
  Data$col[which(Data$z <= breaks[i])] <- cols_breaks[i]}
orderedcolors<-Data[order(Data$Order),"col"]
maps::map("world",interior=F)
points(ssp[f,],pch=22,lwd=0.1,bg=orderedcolors[f],cex=0.5,col=NULL)
rect(xleft =-180,xright = 190,ybottom = -86,ytop = 85,border = "black")
color.legend(xl=-35,xr=35,yb=-83,yt=-67,align="lt",gradient="x",
             legend=round(range(SES_crownMAF,na.rm = T),2),rect.col=cols_breaks)
title("SES-crown MAF")

SES_stemMAF <- (MDT_stem - apply(null_stem[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_stem[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
hist(SES_stemMAF)
Data<-data.frame(Order=1:dim(prPAM_g2[iden,])[1],z=SES_stemMAF)
Data<-Data[order(Data$z),]
k=20
breaks <- getJenksBreaks(Data$z,k)
cols_breaks <- (colorRampPalette(c("linen","lightpink","red3","darkred"),bias=1.1)(k))
Data$col <- cols_breaks[k]
for (i in k:1){ 
  Data$col[which(Data$z <= breaks[i])] <- cols_breaks[i]}
orderedcolors<-Data[order(Data$Order),"col"]
maps::map("world",interior=F)
points(ssp[f,],pch=22,lwd=0.1,bg=orderedcolors[f],cex=0.5,col=NULL)
rect(xleft =-180,xright = 190,ybottom = -86,ytop = 85,border = "black")
color.legend(xl=-35,xr=35,yb=-83,yt=-67,align="lt",gradient="x",
             legend=round(range(SES_stemMAF,na.rm = T),2),rect.col=cols_breaks)
title("SES-stem MAF")

SES_fuseMAF <- (MDT_fuse_non - apply(null_fuse[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_fuse[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
Data<-data.frame(Order=1:dim(prPAM_g2[iden,])[1],z=SES_fuseMAF)
Data<-Data[order(Data$z),]
k=20
breaks <- getJenksBreaks(Data$z,k)
cols_breaks <- (colorRampPalette(c("honeydew","darkseagreen1","forestgreen","black"),bias=1)(k))
Data$col <- cols_breaks[k]
for (i in k:1){ 
  Data$col[which(Data$z <= breaks[i])] <- cols_breaks[i]}
orderedcolors<-Data[order(Data$Order),"col"]
maps::map("world",interior=F)
points(ssp[f,],pch=22,lwd=0.1,bg=orderedcolors[f],cex=0.5,col=NULL)
rect(xleft =-180,xright = 190,ybottom = -86,ytop = 85,border = "black")
color.legend(xl=-35,xr=35,yb=-83,yt=-67,align="lt",gradient="x",
             legend=round(range(SES_fuseMAF,na.rm = T),2),rect.col=cols_breaks)
title("SES-fuse")
dev.off()

pdf(paste(ruta_write,"SES-MAF_NULL1.pdf",sep=""))
par(mfrow = c(1, 1), pty = "s")
null_crown <- MDT_crown_NULL1
null_stem <- MDT_stem_NULL1
null_fuse <- MDT_fuse_non_NULL1
iden <- which(rowSums(prPAM_g2[,-c(1:2)]) > 5)
ssp<-coordinates(prPAM_g2[iden,1:2])
SES_crownMAF <- (MDT_crown - apply(null_crown[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_crown[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
Data<-data.frame(Order=1:dim(prPAM_g2[iden,])[1],z=SES_crownMAF)
Data<-Data[order(Data$z),]
f <- clean_coordinates(as.data.frame(ssp), lon = "LON",lat = "LAT",
                       tests="seas",species=NULL,value="flagged",seas_scale=110)
k=20
breaks <- getJenksBreaks(Data$z,k)
cols_breaks <- (colorRampPalette(c("aliceblue","lightblue","dodgerblue3","darkblue"),bias=1.1)(k))
Data$col <- cols_breaks[k]
for (i in k:1){ 
  Data$col[which(Data$z <= breaks[i])] <- cols_breaks[i]}
orderedcolors<-Data[order(Data$Order),"col"]
maps::map("world",interior=F)
points(ssp[f,],pch=22,lwd=0.1,bg=orderedcolors[f],cex=0.5,col=NULL)
rect(xleft =-180,xright = 190,ybottom = -86,ytop = 85,border = "black")
color.legend(xl=-35,xr=35,yb=-83,yt=-67,align="lt",gradient="x",
             legend=round(range(SES_crownMAF,na.rm = T),2),rect.col=cols_breaks)
title("SES-crown MAF")

SES_stemMAF <- (MDT_stem - apply(null_stem[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_stem[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
hist(SES_stemMAF)
Data<-data.frame(Order=1:dim(prPAM_g2[iden,])[1],z=SES_stemMAF)
Data<-Data[order(Data$z),]
k=20
breaks <- getJenksBreaks(Data$z,k)
cols_breaks <- (colorRampPalette(c("linen","lightpink","red3","darkred"),bias=1.1)(k))
Data$col <- cols_breaks[k]
for (i in k:1){ 
  Data$col[which(Data$z <= breaks[i])] <- cols_breaks[i]}
orderedcolors<-Data[order(Data$Order),"col"]
maps::map("world",interior=F)
points(ssp[f,],pch=22,lwd=0.1,bg=orderedcolors[f],cex=0.5,col=NULL)
rect(xleft =-180,xright = 190,ybottom = -86,ytop = 85,border = "black")
color.legend(xl=-35,xr=35,yb=-83,yt=-67,align="lt",gradient="x",
             legend=round(range(SES_stemMAF,na.rm = T),2),rect.col=cols_breaks)
title("SES-stem MAF")

SES_fuseMAF <- (MDT_fuse_non - apply(null_fuse[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_fuse[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
Data<-data.frame(Order=1:dim(prPAM_g2[iden,])[1],z=SES_fuseMAF)
Data<-Data[order(Data$z),]
k=20
breaks <- getJenksBreaks(Data$z,k)
cols_breaks <- (colorRampPalette(c("honeydew","darkseagreen1","forestgreen","black"),bias=1)(k))
Data$col <- cols_breaks[k]
for (i in k:1){ 
  Data$col[which(Data$z <= breaks[i])] <- cols_breaks[i]}
orderedcolors<-Data[order(Data$Order),"col"]
maps::map("world",interior=F)
points(ssp[f,],pch=22,lwd=0.1,bg=orderedcolors[f],cex=0.5,col=NULL)
rect(xleft =-180,xright = 190,ybottom = -86,ytop = 85,border = "black")
color.legend(xl=-35,xr=35,yb=-83,yt=-67,align="lt",gradient="x",
             legend=round(range(SES_fuseMAF,na.rm = T),2),rect.col=cols_breaks)
title("SES-fuse")
dev.off()
