# Estimates mean age of families, mean phylogenetic fuse of families, and phylogenetic species variability (PSV) across three super biome categories: arid, temperate, tropical.
# The super biome categorization follows the grouping of the original biomes of Olson et al. 2001 described in the Supplementary Methods.
# Continues with R objects generated in previous steps (loads objects).
# Uses additional files.

load("~/Documents/BIOMES/super_biomeRaster.Rdata") # file in data folder
load("~/Documents/NATURE/DISTS.july19/prPAM_FAM_grid.r") # file in data folder

load(paste(ruta_write,"WeightedMAFS.RData",sep=""))
cat("Maximum number of species per grid-cell:",max(rowSums(prPAM_g[,-c(1:2)])),"species","\n")
which(rowSums(prPAM_g[,-c(1:2)]) == 0) -> iden_inf
prPAM_g2<- prPAM_g[-iden_inf,]
iden <- which(rowSums(prPAM_g2[,-c(1:2)]) > 5)
pointos <- SpatialPoints(prPAM_g2[iden,1:2])
proj4string(pointos) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
biome_points<- cbind(prPAM_g2[iden,1:2],raster::extract(super_biomes,pointos))
colnames(biome_points)[3]<-"BIOME"
which(is.na(biome_points$BIOME)) -> iden_bim
biome_points <- biome_points[-which(is.na(biome_points$BIOME)),]
psv <- psv(samp = prPAM_g2[iden,-c(1:2)],tree = tre.pruned)
PSV <- psv$PSVs
Stem <- MDT_stem
Crown <- MDT_crown
Fuse <- MDT_fuse_non
load(paste(ruta_write," MDT_stem_NULL2.Rdata",sep=""))
load(paste(ruta_write," MDT_crown_NULL2.Rdata",sep=""))
load(paste(ruta_write," MDT_fuse_non_NULL2.Rdata",sep=""))
load(paste(ruta_write," MDT_stem_NULL1.Rdata",sep=""))
load(paste(ruta_write," MDT_crown_NULL1.Rdata",sep=""))
load(paste(ruta_write," MDT_fuse_non_NULL1.Rdata",sep=""))
SES_stem_N1 <- (MDT_stem - apply(MDT_stem_NULL1[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(MDT_stem_NULL1[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
SES_crown_N1 <- (MDT_crown - apply(MDT_crown_NULL1[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(MDT_crown_NULL1[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
SES_fuse_N1 <- (MDT_fuse_non - apply(MDT_fuse_non_NULL1[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(MDT_fuse_non_NULL1[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
SES_stem_N2 <- (MDT_stem - apply(MDT_stem_NULL2[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(MDT_stem_NULL2[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
SES_crown_N2 <- (MDT_crown - apply(MDT_crown_NULL2[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(MDT_crown_NULL2[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
SES_fuse_N2 <- (MDT_fuse_non - apply(MDT_fuse_non_NULL2[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(MDT_fuse_non_NULL2[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))

load(paste(ruta_write,"PSV_random_rich.Rdata",sep=""))
load(paste(ruta_write,"PSV_random_freq.Rdata",sep=""))
SES_PSV_N1 <- (PSV - apply(PSV_random_rich[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(PSV_random_rich[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
SES_PSV_N2 <- (PSV - apply(PSV_random_freq[,-c(1:2)],MARGIN=1,FUN=mean,na.rm=T))/(apply(PSV_random_freq[,-c(1:2)],MARGIN=1,FUN=sd,na.rm=T))
SR <-  rowSums(prPAM_g2[iden,-c(1:2)]); print(min(SR))
data<-data.frame(LAT=biome_points$LAT,LON=biome_points$LON,Stem=Stem[-iden_bim],Crown=Crown[-iden_bim],
                 FUSE=Fuse[-iden_bim],
                 PSV=PSV[-iden_bim],
                 SES_crown_1 = SES_crown_N1[-iden_bim],
                 SES_stem_1 = SES_stem_N1[-iden_bim],
                 SES_PSV_1 = SES_PSV_N1[-iden_bim],
                 SES_fuse_1 = SES_fuse_N1[-iden_bim],
                 SES_crown_2 = SES_crown_N2[-iden_bim],
                 SES_stem_2 = SES_stem_N2[-iden_bim],
                 SES_PSV_2 = SES_PSV_N2[-iden_bim],
                 SES_fuse_2 = SES_fuse_N2[-iden_bim],
                 Super_biome = as.factor(biome_points$BIOME),
                 SR=SR[-iden_bim])
unique(data$Super_biome)
levels(data$Super_biome) <- c("Temperate","Temperate","Tropical","Arid"); levels(data$Super_biome)
ssp<-coordinates(data[,2:1])
iden <- clean_coordinates(as.data.frame(ssp), lon = "LON",lat = "LAT",
                       tests="seas",species=NULL,value="flagged",seas_scale=110)

stem <- data[iden,] %>%
  mutate(Super_biome = fct_reorder(Super_biome, Stem)) %>%
  ggplot(aes(x=Super_biome, y=Stem,fill=Super_biome)) +
  geom_violin(alpha=0.6) +
  theme(legend.position="none") +
  labs(title = 'Mean family stem age') +
  geom_hline(yintercept = mean(data[iden,"Stem"])) +
  geom_hline(yintercept = mean(data[iden,"Stem"][which(data$Super_biome[iden]=="Arid")]),col="red")+
  geom_hline(yintercept = mean(data[iden,"Stem"][which(data$Super_biome[iden]=="Temperate")]),col="green")+
  geom_hline(yintercept = mean(data[iden,"Stem"][which(data$Super_biome[iden]=="Tropical")]),col="blue")+
  geom_hline(yintercept = mlv(data[iden,"Stem"],method="Venter"),lty=3) +
  geom_hline(yintercept = mlv(data[iden,"Stem"][which(data$Super_biome[iden]=="Arid")],method="Venter"),col="red",lty=3) +
  geom_hline(yintercept = mlv(data[iden,"Stem"][which(data$Super_biome[iden]=="Temperate")],method="Venter"),col="green",lty=3)+
  geom_hline(yintercept = mlv(data[iden,"Stem"][which(data$Super_biome[iden]=="Tropical")],method="Venter"),col="blue",lty=3)
summary(aov(Stem ~ Super_biome,data=data[iden,]))

crown <- data[iden,] %>%
  mutate(Super_biome = fct_reorder(Super_biome, Stem)) %>%
  ggplot(aes(x=Super_biome, y=Crown,fill=Super_biome)) +
  geom_violin(alpha=0.6) +
  theme(legend.position="none") +
  labs(title = 'Mean family Crown age') +
  geom_hline(yintercept = mean(data[iden,"Crown"])) +
  geom_hline(yintercept = mean(data[iden,"Crown"][which(data$Super_biome[iden]=="Arid")]),col="red")+
  geom_hline(yintercept = mean(data[iden,"Crown"][which(data$Super_biome[iden]=="Temperate")]),col="green")+
  geom_hline(yintercept = mean(data[iden,"Crown"][which(data$Super_biome[iden]=="Tropical")]),col="blue")+
  geom_hline(yintercept = mlv(data[iden,"Crown"],method="Venter"),lty=3) +
  geom_hline(yintercept = mlv(data[iden,"Crown"][which(data$Super_biome[iden]=="Arid")],method="Venter"),col="red",lty=3) +
  geom_hline(yintercept = mlv(data[iden,"Crown"][which(data$Super_biome[iden]=="Temperate")],method="Venter"),col="green",lty=3)+
  geom_hline(yintercept = mlv(data[iden,"Crown"][which(data$Super_biome[iden]=="Tropical")],method="Venter"),col="blue",lty=3)
summary(aov(Crown ~ Super_biome,data=data[iden,]))


fuse <- data[iden,] %>%
  mutate(Super_biome = fct_reorder(Super_biome, Stem)) %>%
  ggplot(aes(x=Super_biome, y=FUSE,fill=Super_biome)) +
  geom_violin(alpha=0.6) +
  theme(legend.position="none") +
  labs(title = 'Mean Fuse of Families') +
  geom_hline(yintercept = mean(data[iden,"FUSE"])) +
  geom_hline(yintercept = mean(data[iden,"FUSE"][which(data$Super_biome[iden]=="Arid")]),col="red")+
  geom_hline(yintercept = mean(data[iden,"FUSE"][which(data$Super_biome[iden]=="Temperate")]),col="green")+
  geom_hline(yintercept = mean(data[iden,"FUSE"][which(data$Super_biome[iden]=="Tropical")]),col="blue")+
  geom_hline(yintercept = mlv(data[iden,"FUSE"],method="Venter"),lty=3) +
  geom_hline(yintercept = mlv(data[iden,"FUSE"][which(data$Super_biome[iden]=="Arid")],method="Venter"),col="red",lty=3) +
  geom_hline(yintercept = mlv(data[iden,"FUSE"][which(data$Super_biome[iden]=="Temperate")],method="Venter"),col="green",lty=3)+
  geom_hline(yintercept = mlv(data[iden,"FUSE"][which(data$Super_biome[iden]=="Tropical")],method="Venter"),col="blue",lty=3)
summary(aov(FUSE ~ Super_biome,data=data[iden,]))

psv <- data[iden,] %>%
  mutate(Super_biome = fct_reorder(Super_biome, Stem)) %>%
  ggplot(aes(x=Super_biome, y=PSV,fill=Super_biome)) +
  geom_violin(alpha=0.6) +
  theme(legend.position="none") +
  labs(title = 'PSV') +
  geom_hline(yintercept = mean(data[iden,"PSV"],na.rm=T)) +
  geom_hline(yintercept = mean(data[iden,"PSV"][which(data$Super_biome[iden]=="Arid")],na.rm=T),col="red")+
  geom_hline(yintercept = mean(data[iden,"PSV"][which(data$Super_biome[iden]=="Temperate")],na.rm=T),col="green")+
  geom_hline(yintercept = mean(data[iden,"PSV"][which(data$Super_biome[iden]=="Tropical")],na.rm=T),col="blue")+
  geom_hline(yintercept = mlv(data[iden,"PSV"],method="Venter",na.rm=T),lty=3) +
  geom_hline(yintercept = mlv(data[iden,"PSV"][which(data$Super_biome[iden]=="Arid")],method="Venter",na.rm=T),col="red",lty=3) +
  geom_hline(yintercept = mlv(data[iden,"PSV"][which(data$Super_biome[iden]=="Temperate")],method="Venter",na.rm=T),col="green",lty=3)+
  geom_hline(yintercept = mlv(data[iden,"PSV"][which(data$Super_biome[iden]=="Tropical")],method="Venter",na.rm=T),col="blue",lty=3)
summary(aov(PSV ~ Super_biome,data=data[iden,]))

ses.psv1 <- data[iden,] %>%
  mutate(Super_biome = fct_reorder(Super_biome, Stem)) %>%
  ggplot(aes(x=Super_biome, y=SES_PSV_1,fill=Super_biome)) +
  geom_violin(alpha=0.6) +
  theme(legend.position="none") +
  labs(title = 'SES_PSV_1') +
  geom_hline(yintercept = mean(data[iden,"SES_PSV_1"],na.rm=T)) +
  geom_hline(yintercept = mean(data[iden,"SES_PSV_1"][which(data$Super_biome[iden]=="Arid")],na.rm=T),col="red")+
  geom_hline(yintercept = mean(data[iden,"SES_PSV_1"][which(data$Super_biome[iden]=="Temperate")],na.rm=T),col="green")+
  geom_hline(yintercept = mean(data[iden,"SES_PSV_1"][which(data$Super_biome[iden]=="Tropical")],na.rm=T),col="blue")+
  geom_hline(yintercept = mlv(data[iden,"SES_PSV_1"],method="Venter",na.rm=T),lty=3) +
  geom_hline(yintercept = mlv(data[iden,"SES_PSV_1"][which(data$Super_biome[iden]=="Arid")],method="Venter",na.rm=T),col="red",lty=3) +
  geom_hline(yintercept = mlv(data[iden,"SES_PSV_1"][which(data$Super_biome[iden]=="Temperate")],method="Venter",na.rm=T),col="green",lty=3)+
  geom_hline(yintercept = mlv(data[iden,"SES_PSV_1"][which(data$Super_biome[iden]=="Tropical")],method="Venter",na.rm=T),col="blue",lty=3)
summary(aov(SES_PSV_1 ~ Super_biome,data=data[iden,]))

ses.psv2 <- data[iden,] %>%
  mutate(Super_biome = fct_reorder(Super_biome, Stem)) %>%
  ggplot(aes(x=Super_biome, y=SES_PSV_2,fill=Super_biome)) +
  geom_violin(alpha=0.6) +
  theme(legend.position="none") +
  labs(title = 'SES_PSV_2') +
  geom_hline(yintercept = mean(data[iden,"SES_PSV_2"],na.rm=T)) +
  geom_hline(yintercept = mean(data[iden,"SES_PSV_2"][which(data$Super_biome[iden]=="Arid")],na.rm=T),col="red")+
  geom_hline(yintercept = mean(data[iden,"SES_PSV_2"][which(data$Super_biome[iden]=="Temperate")],na.rm=T),col="green")+
  geom_hline(yintercept = mean(data[iden,"SES_PSV_2"][which(data$Super_biome[iden]=="Tropical")],na.rm=T),col="blue")+
  geom_hline(yintercept = mlv(data[iden,"SES_PSV_2"],method="Venter",na.rm=T),lty=3) +
  geom_hline(yintercept = mlv(data[iden,"SES_PSV_2"][which(data$Super_biome[iden]=="Arid")],method="Venter",na.rm=T),col="red",lty=3) +
  geom_hline(yintercept = mlv(data[iden,"SES_PSV_2"][which(data$Super_biome[iden]=="Temperate")],method="Venter",na.rm=T),col="green",lty=3)+
  geom_hline(yintercept = mlv(data[iden,"SES_PSV_2"][which(data$Super_biome[iden]=="Tropical")],method="Venter",na.rm=T),col="blue",lty=3)
summary(aov(SES_PSV_2 ~ Super_biome,data=data[iden,]))

ses.stem1 <- data[iden,] %>%
  mutate(Super_biome = fct_reorder(Super_biome, Stem)) %>%
  ggplot(aes(x=Super_biome, y=SES_stem_1,fill=Super_biome)) + 
  geom_violin(alpha=0.6) +
  theme(legend.position="none") +
  labs(title = 'SES_stem_1') +
  geom_hline(yintercept = mean(data[iden,"SES_stem_1"])) +
  geom_hline(yintercept = mean(data[iden,"SES_stem_1"][which(data$Super_biome[iden]=="Arid")]),col="red")+
  geom_hline(yintercept = mean(data[iden,"SES_stem_1"][which(data$Super_biome[iden]=="Temperate")]),col="green")+
  geom_hline(yintercept = mean(data[iden,"SES_stem_1"][which(data$Super_biome[iden]=="Tropical")]),col="blue")+
  geom_hline(yintercept = mlv(data[iden,"SES_stem_1"],method="Venter"),lty=3) +
  geom_hline(yintercept = mlv(data[iden,"SES_stem_1"][which(data$Super_biome[iden]=="Arid")],method="Venter"),col="red",lty=3) +
  geom_hline(yintercept = mlv(data[iden,"SES_stem_1"][which(data$Super_biome[iden]=="Temperate")],method="Venter"),col="green",lty=3)+
  geom_hline(yintercept = mlv(data[iden,"SES_stem_1"][which(data$Super_biome[iden]=="Tropical")],method="Venter"),col="blue",lty=3)
summary(aov(SES_stem_1 ~ Super_biome,data=data[iden,]))

ses.stem2 <- data[iden,] %>%
  mutate(Super_biome = fct_reorder(Super_biome, Stem)) %>%
  ggplot(aes(x=Super_biome, y=SES_stem_2,fill=Super_biome)) + 
  geom_violin(alpha=0.6) +
  theme(legend.position="none") +
  labs(title = 'SES_stem_2') +
  geom_hline(yintercept = mean(data[iden,"SES_stem_2"])) +
  geom_hline(yintercept = mean(data[iden,"SES_stem_2"][which(data$Super_biome[iden]=="Arid")]),col="red")+
  geom_hline(yintercept = mean(data[iden,"SES_stem_2"][which(data$Super_biome[iden]=="Temperate")]),col="green")+
  geom_hline(yintercept = mean(data[iden,"SES_stem_2"][which(data$Super_biome[iden]=="Tropical")]),col="blue")+
  geom_hline(yintercept = mlv(data[iden,"SES_stem_2"],method="Venter"),lty=3) +
  geom_hline(yintercept = mlv(data[iden,"SES_stem_2"][which(data$Super_biome[iden]=="Arid")],method="Venter"),col="red",lty=3) +
  geom_hline(yintercept = mlv(data[iden,"SES_stem_2"][which(data$Super_biome[iden]=="Temperate")],method="Venter"),col="green",lty=3)+
  geom_hline(yintercept = mlv(data[iden,"SES_stem_2"][which(data$Super_biome[iden]=="Tropical")],method="Venter"),col="blue",lty=3)
summary(aov(SES_stem_2 ~ Super_biome,data=data[iden,]))

ses.crown1 <- data[iden,] %>%
  mutate(Super_biome = fct_reorder(Super_biome, Stem)) %>%
  ggplot(aes(x=Super_biome, y=SES_crown_1,fill=Super_biome)) + 
  geom_violin(alpha=0.6) +
  theme(legend.position="none") +
  labs(title = 'SES_crown_1') +
  geom_hline(yintercept = mean(data[iden,"SES_crown_1"])) +
  geom_hline(yintercept = mean(data[iden,"SES_crown_1"][which(data$Super_biome[iden]=="Arid")]),col="red")+
  geom_hline(yintercept = mean(data[iden,"SES_crown_1"][which(data$Super_biome[iden]=="Temperate")]),col="green")+
  geom_hline(yintercept = mean(data[iden,"SES_crown_1"][which(data$Super_biome[iden]=="Tropical")]),col="blue")+
  geom_hline(yintercept = mlv(data[iden,"SES_crown_1"],method="Venter"),lty=3) +
  geom_hline(yintercept = mlv(data[iden,"SES_crown_1"][which(data$Super_biome[iden]=="Arid")],method="Venter"),col="red",lty=3) +
  geom_hline(yintercept = mlv(data[iden,"SES_crown_1"][which(data$Super_biome[iden]=="Temperate")],method="Venter"),col="green",lty=3)+
  geom_hline(yintercept = mlv(data[iden,"SES_crown_1"][which(data$Super_biome[iden]=="Tropical")],method="Venter"),col="blue",lty=3)
summary(aov(SES_crown ~ Super_biome,data=data[iden,]))

ses.crown2 <-data[iden,] %>%
  mutate(Super_biome = fct_reorder(Super_biome, Stem)) %>%
  ggplot(aes(x=Super_biome, y=SES_crown_2,fill=Super_biome)) + 
  geom_violin(alpha=0.6) +
  theme(legend.position="none") +
  labs(title = 'SES_crown_2') +
  geom_hline(yintercept = mean(data[iden,"SES_crown_2"])) +
  geom_hline(yintercept = mean(data[iden,"SES_crown_2"][which(data$Super_biome[iden]=="Arid")]),col="red")+
  geom_hline(yintercept = mean(data[iden,"SES_crown_2"][which(data$Super_biome[iden]=="Temperate")]),col="green")+
  geom_hline(yintercept = mean(data[iden,"SES_crown_2"][which(data$Super_biome[iden]=="Tropical")]),col="blue")+
  geom_hline(yintercept = mlv(data[iden,"SES_crown_2"],method="Venter"),lty=3) +
  geom_hline(yintercept = mlv(data[iden,"SES_crown_2"][which(data$Super_biome[iden]=="Arid")],method="Venter"),col="red",lty=3) +
  geom_hline(yintercept = mlv(data[iden,"SES_crown_2"][which(data$Super_biome[iden]=="Temperate")],method="Venter"),col="green",lty=3)+
  geom_hline(yintercept = mlv(data[iden,"SES_crown_2"][which(data$Super_biome[iden]=="Tropical")],method="Venter"),col="blue",lty=3)
summary(aov(SES_crown_2 ~ Super_biome,data=data[iden,]))

ses.fuse1 <- data[iden,] %>%
  mutate(Super_biome = fct_reorder(Super_biome, Stem)) %>%
  ggplot(aes(x=Super_biome, y=SES_fuse_1,fill=Super_biome)) + 
  geom_violin(alpha=0.6) +
  theme(legend.position="none") +
  labs(title = 'SES_fuse_1') +
  geom_hline(yintercept = mean(data[iden,"SES_fuse_1"])) +
  geom_hline(yintercept = mean(data[iden,"SES_fuse_1"][which(data$Super_biome[iden]=="Arid")]),col="red")+
  geom_hline(yintercept = mean(data[iden,"SES_fuse_1"][which(data$Super_biome[iden]=="Temperate")]),col="green")+
  geom_hline(yintercept = mean(data[iden,"SES_fuse_1"][which(data$Super_biome[iden]=="Tropical")]),col="blue")+
  geom_hline(yintercept = mlv(data[iden,"SES_fuse_1"],method="Venter"),lty=3) +
  geom_hline(yintercept = mlv(data[iden,"SES_fuse_1"][which(data$Super_biome[iden]=="Arid")],method="Venter"),col="red",lty=3) +
  geom_hline(yintercept = mlv(data[iden,"SES_fuse_1"][which(data$Super_biome[iden]=="Temperate")],method="Venter"),col="green",lty=3)+
  geom_hline(yintercept = mlv(data[iden,"SES_fuse_1"][which(data$Super_biome[iden]=="Tropical")],method="Venter"),col="blue",lty=3)
summary(aov(SES_fuse_1 ~ Super_biome,data=data[iden,]))

ses.fuse2 <-data[iden,] %>%
  mutate(Super_biome = fct_reorder(Super_biome, Stem)) %>%
  ggplot(aes(x=Super_biome, y=SES_fuse_2,fill=Super_biome)) + 
  geom_violin(alpha=0.6) +
  theme(legend.position="none") +
  labs(title = 'SES_fuse_2') +
  geom_hline(yintercept = mean(data[iden,"SES_fuse_2"])) +
  geom_hline(yintercept = mean(data[iden,"SES_fuse_2"][which(data$Super_biome[iden]=="Arid")]),col="red")+
  geom_hline(yintercept = mean(data[iden,"SES_fuse_2"][which(data$Super_biome[iden]=="Temperate")]),col="green")+
  geom_hline(yintercept = mean(data[iden,"SES_fuse_2"][which(data$Super_biome[iden]=="Tropical")]),col="blue")+
  geom_hline(yintercept = mlv(data[iden,"SES_fuse_2"],method="Venter"),lty=3) +
  geom_hline(yintercept = mlv(data[iden,"SES_fuse_2"][which(data$Super_biome[iden]=="Arid")],method="Venter"),col="red",lty=3) +
  geom_hline(yintercept = mlv(data[iden,"SES_fuse_2"][which(data$Super_biome[iden]=="Temperate")],method="Venter"),col="green",lty=3)+
  geom_hline(yintercept = mlv(data[iden,"SES_fuse_2"][which(data$Super_biome[iden]=="Tropical")],method="Venter"),col="blue",lty=3)
summary(aov(SES_fuse_2 ~ Super_biome,data=data[iden,]))

pdf(paste(ruta_write,"8.BiomeMAF_SBcomplete_20SR.pdf",sep=""))
gridExtra::grid.arrange(stem, crown, fuse, psv, nrow = 2)
gridExtra::grid.arrange(ses.stem1, ses.stem2, ses.crown1, ses.crown2, nrow = 2)
gridExtra::grid.arrange(ses.psv1, ses.psv2, ses.fuse1, ses.fuse2, nrow = 2)
dev.off()

sink(paste(ruta_write,"Superbiome_tests.txt"))
cat("STEM ----------------","\n")
summary(aov(Stem ~ Super_biome,data=data[iden,]))
cat("CROWN ----------------","\n")
summary(aov(Crown ~ Super_biome,data=data[iden,]))
cat("PSV ----------------","\n")
summary(aov(PSV ~ Super_biome,data=data[iden,]))
cat("SES CROWN ----------------","\n")
summary(aov(SES_crown ~ Super_biome,data=data[iden,]))
cat("SES STEM ----------------","\n")
summary(aov(SES_stem ~ Super_biome,data=data[iden,]))
cat("SES PSV ----------------","\n")
summary(aov(SES_PSV ~ Super_biome,data=data[iden,]))
sink()
