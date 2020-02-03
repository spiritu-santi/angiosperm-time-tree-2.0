# The code reads the presence absence matrix (PAM) and matches family names between the PAM and the age file generated in previous steps. 
# Then it estimates and plots species (and family) richness maps across the globe using the grid resolution used for the PAM.
# Continues with R objects generated in previous steps.
# Uses additional files

load("~/prPAM_FAM_grid.r") ### file in data folder

cat("Maximum number of species per grid-cell:",max(rowSums(prPAM_g[,-c(1:2)])),"species","\n")
which(rowSums(prPAM_g[,-c(1:2)]) == 0) -> iden_inf
prPAM_g2 <- prPAM_g[-iden_inf,]
nana <- which(is.na(match(colnames(prPAM_g2),rownames(resB)[fams])));colnames(prPAM_g2)[nana]
which(is.na(match(colnames(prPAM_g2),rownames(resB)[fams])))
prPAM_fam <- prPAM_g2
prPAM_fam[,-c(1:2)] <- decostand(prPAM_fam[,-c(1:2)],method="pa")
ssp <- coordinates(prPAM_fam[,1:2])
which(rowSums(prPAM_g2[,-c(1:2)]) > 5) -> iden
f <- clean_coordinates(as.data.frame(ssp[iden,]), lon = "LON",lat = "LAT",
                       tests="seas",species=NULL,value="flagged",seas_scale = 110)
k <- 20
richness_fam <- rowSums(prPAM_fam[iden,-c(1:2)]);range(richness_fam)
Data<-data.frame(Order=1:dim(prPAM_fam[iden,])[1],z=richness_fam)
Data<-Data[order(Data$z),]; ssp<-coordinates(prPAM_fam[iden,1:2])
breaks <- getJenksBreaks(Data$z[f],k)
cols_breaks <- rev(viridis(k))
Data$col <- cols_breaks[k]
for (i in k:1){Data$col[which(Data$z <= breaks[i])] <- cols_breaks[i]}
orderedcolors<-Data[order(Data$Order),"col"]
pdf(paste(ruta_write,"4.Richness_MAPS.pdf",sep=""))
maps::map("world",interior=F)
title("Family richness per 2x2 grid-cell")
points(ssp[f,],pch=22,lwd=0.4,bg=orderedcolors[f],cex= 0.45,col=orderedcolors[f])
rect(xleft =-180,xright = 190,ybottom = -86,ytop = 85,border = "black")
color.legend(xl=-35,xr=35,yb=-83,yt=-67,align="lt",gradient="x",
             legend=c(range(richness_fam)),rect.col=cols_breaks)
abline(h=c(23.5,-23.5),lwd=2,lty=3)
richness <- rowSums(prPAM_g2[iden,-c(1:2)])
Data <- data.frame(Order=1:dim(prPAM_g2[iden,])[1],z=richness)
Data <- Data[order(Data$z),];ssp<-coordinates(prPAM_g2[iden,1:2])
breaks <- getJenksBreaks(Data$z[f],k)
cols_breaks <- rev(viridis(k))
Data$col <- cols_breaks[k]
for (i in k:1){Data$col[which(Data$z <= breaks[i])] <- cols_breaks[i]}
orderedcolors<-Data[order(Data$Order),"col"]
maps::map("world",interior=F)
title("Species richness per 2x2 grid-cell")
points(ssp[f,],pch=22,lwd=0.4,col=orderedcolors[f],bg=orderedcolors[f],cex= 0.45)
rect(xleft =-180,xright = 190,ybottom = -86,ytop = 85,border = "black")
color.legend(xl=-35,xr=35,yb=-83,yt=-67,align="lt",gradient="x",
             legend=c(range(richness)),rect.col=cols_breaks)
abline(h=c(23.5,-23.5),lwd=2,lty=3)

par(mfrow = c(1, 1), pty = "s")
plot(log(richness),log(richness_fam),cex=0.7,pch=21,bg=rgb(1,0,0,0.5),lwd=0.6,
     ylab="Family Richness (log-trasnformed)",xlab="Species Richness (log-trasnformed)")
abline(lm(log(richness_fam)~log(richness)))
summ_lm <- summary(lm(log(richness_fam)~log(richness)))
legend("bottomright",inset=0.03,legend=paste("b = ",round(summ_lm$coefficients[2,1],2),";","F = ",round(summ_lm$fstatistic[1],2),";",
                                               "p-val < 0.001",";","R2* = ",round(summ_lm$adj.r.squared,2),sep=""))
dev.off()
