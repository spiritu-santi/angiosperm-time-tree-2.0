###### 21. COMPARE PSV ACROSS STRATEGIES (COMPLETE ONLY) ########
#### DO STEPS 1-4 AND 12 FOR EACH TREE, THEN PROCEED
iden <- which(rowSums(prPAM_g2[,-c(1:2)]) > 5)
ssp<-coordinates(prPAM_g2[iden,1:2])
f <- clean_coordinates(as.data.frame(ssp), lon = "LON",lat = "LAT",
                       tests="seas",species=NULL,value="flagged",seas_scale=110)

iden <- which(rowSums(prPAM_g2[,-c(1:2)]) > 5)
psvSA <- psv(samp = prPAM_g2[iden,-c(1:2)],tree = tre.pruned); psvSA <- psvSA$PSVs[iden];psvSA <- psvSA[f]
psvSB <- psv(samp = prPAM_g2[iden,-c(1:2)],tree = tre.pruned); psvSB <- psvSB$PSVs[iden]; psvSB <- psvSB[f]
psvSD <- psv(samp = prPAM_g2[iden,-c(1:2)],tree = tre.pruned); psvSD <- psvSD$PSVs[iden]; psvSD <- psvSD[f]


ruta_write<-"~/Desktop/NATURE_rev/RES_SAcomplete_v3/"
load(paste(ruta_write,"PSV_random_rich.Rdata",sep=""))
PSV_random_rich_SA <- PSV_random_rich
load(paste(ruta_write,"PSV_random_freq.Rdata",sep=""))
PSV_random_freq_SA <- PSV_random_freq
ruta_write<-"~/Desktop/NATURE_rev/RES_SBcomplete_v3/"
load(paste(ruta_write,"PSV_random_rich.Rdata",sep=""))
PSV_random_rich_SB <- PSV_random_rich
load(paste(ruta_write,"PSV_random_freq.Rdata",sep=""))
PSV_random_freq_SB <- PSV_random_freq
ruta_write<-"~/Desktop/NATURE_rev/RES_SDcomplete_v3/"
load(paste(ruta_write,"PSV_random_rich.Rdata",sep=""))
PSV_random_rich_SD <- PSV_random_rich
load(paste(ruta_write,"PSV_random_freq.Rdata",sep=""))
PSV_random_freq_SD <- PSV_random_freq

SESnull1_psvSA <- (psvSA - apply(PSV_random_rich_SA[f,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(PSV_random_rich_SA[f,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
SESnull1_psvSB <- (psvSB - apply(PSV_random_rich_SB[f,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(PSV_random_rich_SB[f,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
SESnull1_psvSD <- (psvSD - apply(PSV_random_rich_SD[f,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(PSV_random_rich_SD[f,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))

SESnull2_psvSA <- (psvSA - apply(PSV_random_freq_SA[f,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(PSV_random_freq_SA[f,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
SESnull2_psvSB <- (psvSB - apply(PSV_random_freq_SB[f,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(PSV_random_freq_SB[f,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
SESnull2_psvSD <- (psvSD - apply(PSV_random_freq_SD[f,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(PSV_random_freq_SD[f,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))

pdf("~/Desktop/NATURE_rev/2.PSV_nulls_comparison_strategies.pdf",useDingbats = F)
plot(psvSB,psvSA,xlim=c(0.5,1),ylim=c(0.5,1),type="n",xlab="Relaxed Calibration (RC)",ylab="Constrained Calibration (CC)")
abline(a=0,b=1)
points(psvSB,psvSA,pch=21,bg="red",cex=0.8)
title("PSV")
plot(psvSB,psvSD,xlim=c(0.5,1),ylim=c(0.5,1),type="n",xlab="Relaxed Calibration (RC)",ylab="Unconstrained Calibration (CC)")
abline(a=0,b=1)
points(psvSB,psvSD,pch=21,bg="red",cex=0.8)
title("PSV")

plot(SESnull1_psvSB,SESnull1_psvSA,xlim=c(-20,20),ylim=c(-20,20),type="n",xlab="Relaxed Calibration (RC)",ylab="Constrained Calibration (CC)")
abline(a=0,b=1)
points(SESnull1_psvSB,SESnull1_psvSA,pch=21,bg="red",cex=0.8)
title("SES PSV (null 1")
plot(SESnull1_psvSB,SESnull1_psvSD,xlim=c(-20,20),ylim=c(-20,20),type="n",xlab="Relaxed Calibration (RC)",ylab="Unconstrained Calibration (CC)")
abline(a=0,b=1)
points(SESnull1_psvSB,SESnull1_psvSD,pch=21,bg="red",cex=0.8)
title("SES PSV (null 1")

plot(SESnull2_psvSB,SESnull2_psvSA,xlim=c(-20,20),ylim=c(-20,20),type="n",xlab="Relaxed Calibration (RC)",ylab="Constrained Calibration (CC)")
abline(a=0,b=1)
points(SESnull2_psvSB,SESnull2_psvSA,pch=21,bg="red",cex=0.8)
title("SES PSV (null 2")
plot(SESnull2_psvSB,SESnull2_psvSD,xlim=c(-20,20),ylim=c(-20,20),type="n",xlab="Relaxed Calibration (RC)",ylab="Unconstrained Calibration (CC)")
abline(a=0,b=1)
points(SESnull2_psvSB,SESnull2_psvSD,pch=21,bg="red",cex=0.8)
title("SES PSV (null 2")
dev.off()
