###### 19. COMPARE MAFS ACROSS STRATEGIES (COMPLETE ONLY) ########
#### DO STEPS 1-4, IRRESPECTIVE OF TREE, THEN PROCEED
iden <- which(rowSums(prPAM_g2[,-c(1:2)]) > 5)
ssp<-coordinates(prPAM_g2[iden,1:2])
f <- clean_coordinates(as.data.frame(ssp), lon = "LON",lat = "LAT",
                       tests="seas",species=NULL,value="flagged",seas_scale = 110)

#### NULL 2 #####
load("~/Desktop/NATURE_rev/RES_SAcomplete_v3/WeightedMAFS.RData")
MAFcrown_SA <- MDT_crown
MAFstem_SA <- MDT_stem
MAFfuse_SA <- MDT_fuse_non
load("~/Desktop/NATURE_rev/RES_SAcomplete_v3/ MDT_crown_NULL2.Rdata")
null_crown <- MDT_crown_NULL2
SES_crownMAF_SA <- (MAFcrown_SA - apply(null_crown[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_crown[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
load("~/Desktop/NATURE_rev/RES_SAcomplete_v3/ MDT_stem_NULL2.Rdata")
null_stem <- MDT_stem_NULL2
SES_stemMAF_SA <- (MAFstem_SA - apply(null_stem[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_stem[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
load("~/Desktop/NATURE_rev/RES_SAcomplete_v3/ MDT_fuse_non_NULL2.Rdata")
null_fuse <- MDT_fuse_non_NULL2
SES_fuseMAF_SA <- (MAFfuse_SA - apply(null_fuse[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_fuse[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))

load("~/Desktop/NATURE_rev/RES_SBcomplete_v3/WeightedMAFS.RData")
MAFcrown_SB <- MDT_crown
MAFstem_SB <- MDT_stem
MAFfuse_SB <- MDT_fuse_non
load("~/Desktop/NATURE_rev/RES_SBcomplete_v3/ MDT_crown_NULL2.Rdata")
null_crown <- MDT_crown_NULL2
SES_crownMAF_SB <- (MAFcrown_SB - apply(null_crown[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_crown[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
load("~/Desktop/NATURE_rev/RES_SBcomplete_v3/ MDT_stem_NULL2.Rdata")
null_stem <- MDT_stem_NULL2
SES_stemMAF_SB <- (MAFstem_SB - apply(null_stem[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_stem[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
load("~/Desktop/NATURE_rev/RES_SBcomplete_v3/ MDT_fuse_non_NULL2.Rdata")
null_fuse <- MDT_fuse_non_NULL2
SES_fuseMAF_SB <- (MAFfuse_SB - apply(null_fuse[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_fuse[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))

load("~/Desktop/NATURE_rev/RES_SDcomplete_v3/WeightedMAFS.RData")
MAFcrown_SD <- MDT_crown
MAFstem_SD <- MDT_stem
MAFfuse_SD <- MDT_fuse_non
load("~/Desktop/NATURE_rev/RES_SDcomplete_v3/ MDT_crown_NULL2.Rdata")
null_crown <- MDT_crown_NULL2
SES_crownMAF_SD <- (MAFcrown_SD - apply(null_crown[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_crown[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
load("~/Desktop/NATURE_rev/RES_SDcomplete_v3/ MDT_stem_NULL2.Rdata")
null_stem <- MDT_stem_NULL2
SES_stemMAF_SD <- (MAFstem_SD - apply(null_stem[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_stem[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
load("~/Desktop/NATURE_rev/RES_SDcomplete_v3/ MDT_fuse_non_NULL2.Rdata")
null_fuse <- MDT_fuse_non_NULL2
SES_fuseMAF_SD <- (MAFfuse_SD - apply(null_fuse[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_fuse[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))

pdf("~/Desktop/NATURE_rev/2.MAF_null2_comparison_strategies.pdf",useDingbats = F)
plot(MAFcrown_SB[f],MAFcrown_SA[f],xlim=c(30,120),ylim=c(30,120),type="n",xlab="Relaxed Calibration (RC)",ylab="Constrained Calibration (CC)")
abline(a=0,b=1)
points(MAFcrown_SB[f],MAFcrown_SA[f],pch=21,bg="blue",cex=0.8)
title("MAF crown")
plot(MAFcrown_SB[f],MAFcrown_SD[f],xlim=c(30,120),ylim=c(30,120),type="n",xlab="Relaxed Calibration (RC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(MAFcrown_SB[f],MAFcrown_SD[f],pch=21,bg="blue",cex=0.8)
title("MAF crown")
plot(MAFcrown_SA[f],MAFcrown_SD[f],xlim=c(30,120),ylim=c(30,120),type="n",xlab="Constrained Calibration (CC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(MAFcrown_SA[f],MAFcrown_SD[f],pch=21,bg="blue",cex=0.8)
title("MAF crown")

plot(MAFstem_SB[f],MAFstem_SA[f],xlim=c(60,180),ylim=c(60,180),type="n",xlab="Relaxed Calibration (RC)",ylab="Constrained Calibration (CC)")
abline(a=0,b=1)
points(MAFstem_SB[f],MAFstem_SA[f],pch=21,bg="red",cex=0.8)
title("MAF stem")
plot(MAFstem_SB[f],MAFstem_SD[f],xlim=c(60,180),ylim=c(60,180),type="n",xlab="Relaxed Calibration (RC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(MAFstem_SB[f],MAFstem_SD[f],pch=21,bg="red",cex=0.8)
title("MAF stem")
plot(MAFstem_SA[f],MAFstem_SD[f],xlim=c(60,180),ylim=c(60,180),type="n",xlab="Constrained Calibration (CC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(MAFstem_SA[f],MAFstem_SD[f],pch=21,bg="red",cex=0.8)
title("MAF stem")

plot(SES_crownMAF_SB[f],SES_crownMAF_SA[f],xlim=c(-5,10),ylim=c(-5,10),type="n",xlab="Relaxed Calibration (RC)",ylab="Constrained Calibration (CC)")
abline(a=0,b=1)
points(SES_crownMAF_SB[f],SES_crownMAF_SA[f],pch=21,bg="blue",cex=0.8)
title("SES MAF crown (NULL 2)")
plot(SES_crownMAF_SB[f],SES_crownMAF_SD[f],xlim=c(-5,10),ylim=c(-5,10),type="n",xlab="Relaxed Calibration (RC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(SES_crownMAF_SB[f],SES_crownMAF_SD[f],pch=21,bg="blue",cex=0.8)
title("SES MAF crown (NULL 2)")
plot(SES_crownMAF_SA[f],SES_crownMAF_SD[f],xlim=c(-5,10),ylim=c(-5,10),type="n",xlab="Constrained Calibration (CC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(SES_crownMAF_SA[f],SES_crownMAF_SD[f],pch=21,bg="blue",cex=0.8)
title("SES MAF crown (NULL 2)")

plot(SES_stemMAF_SB[f],SES_stemMAF_SA[f],xlim=c(-10,5),ylim=c(-10,5),type="n",xlab="Relaxed Calibration (RC)",ylab="Constrained Calibration (CC)")
abline(a=0,b=1)
points(SES_stemMAF_SB[f],SES_stemMAF_SA[f],pch=21,bg="red",cex=0.8)
title("SES MAF stem (NULL 2)")
plot(SES_stemMAF_SB[f],SES_stemMAF_SD[f],xlim=c(-10,5),ylim=c(-10,5),type="n",xlab="Relaxed Calibration (RC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(SES_stemMAF_SB[f],SES_stemMAF_SD[f],pch=21,bg="red",cex=0.8)
title("SES MAF stem (NULL 2)")
plot(SES_stemMAF_SA[f],SES_stemMAF_SD[f],xlim=c(-10,5),ylim=c(-10,5),type="n",xlab="Constrained Calibration (CC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(SES_stemMAF_SA[f],SES_stemMAF_SD[f],pch=21,bg="red",cex=0.8)
title("SES MAF stem (NULL 2)")

plot(SES_fuseMAF_SB[f],SES_fuseMAF_SA[f],xlim=c(-5,10),ylim=c(-5,10),type="n",xlab="Relaxed Calibration (RC)",ylab="Constrained Calibration (CC)")
abline(a=0,b=1)
points(SES_fuseMAF_SB[f],SES_fuseMAF_SA[f],pch=21,bg="green",cex=0.8)
title("SES fuse (NULL 2)")
plot(SES_fuseMAF_SB[f],SES_fuseMAF_SD[f],xlim=c(-5,10),ylim=c(-5,10),type="n",xlab="Relaxed Calibration (RC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(SES_fuseMAF_SB[f],SES_fuseMAF_SD[f],pch=21,bg="green",cex=0.8)
title("SES fuse (NULL 2)")
plot(SES_fuseMAF_SA[f],SES_fuseMAF_SD[f],xlim=c(-5,10),ylim=c(-5,10),type="n",xlab="Constrained Calibration (CC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(SES_fuseMAF_SA[f],SES_fuseMAF_SD[f],pch=21,bg="green",cex=0.8)
title("SES fuse (NULL 2)")

dev.off()
###########

##### NULL 1 ######
load("~/Desktop/NATURE_rev/RES_SAcomplete_v3/WeightedMAFS.RData")
MAFcrown_SA <- MDT_crown
MAFstem_SA <- MDT_stem
MAFfuse_SA <- MDT_fuse_non
load("~/Desktop/NATURE_rev/RES_SAcomplete_v3/ MDT_crown_NULL1.Rdata")
null_crown <- MDT_crown_NULL2
SES_crownMAF_SA <- (MAFcrown_SA - apply(null_crown[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_crown[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
load("~/Desktop/NATURE_rev/RES_SAcomplete_v3/ MDT_stem_NULL1.Rdata")
null_stem <- MDT_stem_NULL2
SES_stemMAF_SA <- (MAFstem_SA - apply(null_stem[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_stem[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
load("~/Desktop/NATURE_rev/RES_SAcomplete_v3/ MDT_fuse_non_NULL1.Rdata")
null_fuse <- MDT_fuse_non_NULL2
SES_fuseMAF_SA <- (MAFfuse_SA - apply(null_fuse[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_fuse[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))

load("~/Desktop/NATURE_rev/RES_SBcomplete_v3/WeightedMAFS.RData")
MAFcrown_SB <- MDT_crown
MAFstem_SB <- MDT_stem
MAFfuse_SB <- MDT_fuse_non
load("~/Desktop/NATURE_rev/RES_SBcomplete_v3/ MDT_crown_NULL1.Rdata")
null_crown <- MDT_crown_NULL2
SES_crownMAF_SB <- (MAFcrown_SB - apply(null_crown[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_crown[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
load("~/Desktop/NATURE_rev/RES_SBcomplete_v3/ MDT_stem_NULL1.Rdata")
null_stem <- MDT_stem_NULL2
SES_stemMAF_SB <- (MAFstem_SB - apply(null_stem[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_stem[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
load("~/Desktop/NATURE_rev/RES_SBcomplete_v3/ MDT_fuse_non_NULL1.Rdata")
null_fuse <- MDT_fuse_non_NULL2
SES_fuseMAF_SB <- (MAFfuse_SB - apply(null_fuse[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_fuse[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))

load("~/Desktop/NATURE_rev/RES_SDcomplete_v3/WeightedMAFS.RData")
MAFcrown_SD <- MDT_crown
MAFstem_SD <- MDT_stem
MAFfuse_SD <- MDT_fuse_non
load("~/Desktop/NATURE_rev/RES_SDcomplete_v3/ MDT_crown_NULL1.Rdata")
null_crown <- MDT_crown_NULL2
SES_crownMAF_SD <- (MAFcrown_SD - apply(null_crown[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_crown[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
load("~/Desktop/NATURE_rev/RES_SDcomplete_v3/ MDT_stem_NULL1.Rdata")
null_stem <- MDT_stem_NULL2
SES_stemMAF_SD <- (MAFstem_SD - apply(null_stem[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_stem[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
load("~/Desktop/NATURE_rev/RES_SDcomplete_v3/ MDT_fuse_non_NULL1.Rdata")
null_fuse <- MDT_fuse_non_NULL2
SES_fuseMAF_SD <- (MAFfuse_SD - apply(null_fuse[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(null_fuse[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))

pdf("~/Desktop/NATURE_rev/2.MAF_null1_comparison_strategies.pdf",useDingbats = F)
plot(MAFcrown_SB[f],MAFcrown_SA[f],xlim=c(30,120),ylim=c(30,120),type="n",xlab="Relaxed Calibration (RC)",ylab="Constrained Calibration (CC)")
abline(a=0,b=1)
points(MAFcrown_SB[f],MAFcrown_SA[f],pch=21,bg="blue",cex=0.8)
title("MAF crown")
plot(MAFcrown_SB[f],MAFcrown_SD[f],xlim=c(30,120),ylim=c(30,120),type="n",xlab="Relaxed Calibration (RC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(MAFcrown_SB[f],MAFcrown_SD[f],pch=21,bg="blue",cex=0.8)
title("MAF crown")
plot(MAFcrown_SA[f],MAFcrown_SD[f],xlim=c(30,120),ylim=c(30,120),type="n",xlab="Constrained Calibration (CC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(MAFcrown_SA[f],MAFcrown_SD[f],pch=21,bg="blue",cex=0.8)
title("MAF crown")

plot(MAFstem_SB[f],MAFstem_SA[f],xlim=c(60,180),ylim=c(60,180),type="n",xlab="Relaxed Calibration (RC)",ylab="Constrained Calibration (CC)")
abline(a=0,b=1)
points(MAFstem_SB[f],MAFstem_SA[f],pch=21,bg="red",cex=0.8)
title("MAF stem")
plot(MAFstem_SB[f],MAFstem_SD[f],xlim=c(60,180),ylim=c(60,180),type="n",xlab="Relaxed Calibration (RC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(MAFstem_SB[f],MAFstem_SD[f],pch=21,bg="red",cex=0.8)
title("MAF stem")
plot(MAFstem_SA[f],MAFstem_SD[f],xlim=c(60,180),ylim=c(60,180),type="n",xlab="Constrained Calibration (CC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(MAFstem_SA[f],MAFstem_SD[f],pch=21,bg="red",cex=0.8)
title("MAF stem")

plot(SES_crownMAF_SB[f],SES_crownMAF_SA[f],xlim=c(-5,10),ylim=c(-5,10),type="n",xlab="Relaxed Calibration (RC)",ylab="Constrained Calibration (CC)")
abline(a=0,b=1)
points(SES_crownMAF_SB[f],SES_crownMAF_SA[f],pch=21,bg="blue",cex=0.8)
title("SES MAF crown (NULL 1)")
plot(SES_crownMAF_SB[f],SES_crownMAF_SD[f],xlim=c(-5,10),ylim=c(-5,10),type="n",xlab="Relaxed Calibration (RC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(SES_crownMAF_SB[f],SES_crownMAF_SD[f],pch=21,bg="blue",cex=0.8)
title("SES MAF crown (NULL 1)")
plot(SES_crownMAF_SA[f],SES_crownMAF_SD[f],xlim=c(-5,10),ylim=c(-5,10),type="n",xlab="Constrained Calibration (CC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(SES_crownMAF_SA[f],SES_crownMAF_SD[f],pch=21,bg="blue",cex=0.8)
title("SES MAF crown (NULL 1)")

plot(SES_stemMAF_SB[f],SES_stemMAF_SA[f],xlim=c(-10,5),ylim=c(-10,5),type="n",xlab="Relaxed Calibration (RC)",ylab="Constrained Calibration (CC)")
abline(a=0,b=1)
points(SES_stemMAF_SB[f],SES_stemMAF_SA[f],pch=21,bg="red",cex=0.8)
title("SES MAF stem (NULL 1)")
plot(SES_stemMAF_SB[f],SES_stemMAF_SD[f],xlim=c(-10,5),ylim=c(-10,5),type="n",xlab="Relaxed Calibration (RC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(SES_stemMAF_SB[f],SES_stemMAF_SD[f],pch=21,bg="red",cex=0.8)
title("SES MAF stem (NULL 1)")
plot(SES_stemMAF_SA[f],SES_stemMAF_SD[f],xlim=c(-10,5),ylim=c(-10,5),type="n",xlab="Constrained Calibration (CC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(SES_stemMAF_SA[f],SES_stemMAF_SD[f],pch=21,bg="red",cex=0.8)
title("SES MAF stem (NULL 1)")

plot(SES_fuseMAF_SB[f],SES_fuseMAF_SA[f],xlim=c(-5,10),ylim=c(-5,10),type="n",xlab="Relaxed Calibration (RC)",ylab="Constrained Calibration (CC)")
abline(a=0,b=1)
points(SES_fuseMAF_SB[f],SES_fuseMAF_SA[f],pch=21,bg="green",cex=0.8)
title("SES fuse (NULL 1)")
plot(SES_fuseMAF_SB[f],SES_fuseMAF_SD[f],xlim=c(-5,10),ylim=c(-5,10),type="n",xlab="Relaxed Calibration (RC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(SES_fuseMAF_SB[f],SES_fuseMAF_SD[f],pch=21,bg="green",cex=0.8)
title("SES fuse (NULL 1)")
plot(SES_fuseMAF_SA[f],SES_fuseMAF_SD[f],xlim=c(-5,10),ylim=c(-5,10),type="n",xlab="Constrained Calibration (CC)",ylab="Unconstrained Calibration (UC)")
abline(a=0,b=1)
points(SES_fuseMAF_SA[f],SES_fuseMAF_SD[f],pch=21,bg="green",cex=0.8)
title("SES fuse (NULL 1)")

dev.off()
