######### NEW FIGURES #######
library(ggnewscale)
tree <- tre.pruned
tree <- drop.tip(tree,tip = which(tree$tip.label=="Cycadaceae"))
circ <- ggtree(tree,layout = "circular")
circ

df <- data.frame(stem=resB_fams$SG_age_Acal[match(resB_fams$Family,tree$tip.label)])
rownames(df) <- resB_fams$Family[match(resB_fams$Family,tree$tip.label)]

df2 <- data.frame(crown=resB_fams$CG_age_Acal[match(resB_fams$Family,tree$tip.label)])
rownames(df2) <- resB_fams$Family[match(resB_fams$Family,tree$tip.label)]

df3 <- data.frame(fuse=resB_fams$Fuse[match(resB_fams$Family,tree$tip.label)])
rownames(df3) <- resB_fams$Family[match(resB_fams$Family,tree$tip.label)]

fuses <- resB_fams$Fuse[match(resB_fams$Family,tree$tip.label)]
data.frame(fuses) -> df4
head(df4)
df4$fuses[order(fuses)][1:124] <- "short"
df4$fuses[order(fuses)][125:249] <- "Inter"
df4$fuses[order(fuses)][250:373] <- "Long"
df4$fuses[order(fuses)][c(374:436)]
df4$fuses <- as.factor(df4$fuses)
rownames(df4) <- resB_fams$Family[match(resB_fams$Family,tree$tip.label)]
levels(df4$fuses)
pdf("~/Desktop/Test_phylo.pdf",useDingbats = F)
k=10
p1 <- gheatmap(circ, df, offset=.1, width=.2, colnames=F) +
  scale_fill_gradientn(colours=colorRampPalette(c("linen","lightpink","red3","darkred"),bias=2)(k))
p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, df2, offset=30, width=.2, colnames=F) +
  scale_fill_gradientn(colours=colorRampPalette(c("aliceblue","lightblue","dodgerblue3","darkblue"),bias=2)(k))
p3 <- p2 + new_scale_fill()
p3 <- gheatmap(p3, df3, offset=50, width=.2, colnames=F) +
  scale_fill_gradientn(colours=colorRampPalette(c("honeydew","darkseagreen1","forestgreen","black"),bias=2)(k))
p4 <- p3 + new_scale_fill()
p4 <- gheatmap(circ, df4, offset=80, width=.2, colnames_angle=95, colnames_offset_y = .25) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")
p4

dev.off()

###############
ruta_write <- "~/Desktop/NATURE_rev/RES_SBcomplete_v3/"
############ READ AGE FILE AND MERGE WITH RICHNESS FILE #####
resB<-read.table(paste(ruta_write,"2.Ages_complete.csv",sep=""),sep=",",header=T)
rich<-read.table("SPP.RICHNESS.csv",heade=T,sep=",")
rich_ord<-rich[which(rich$Family==""),]; dim(rich_ord)
resB$Order_richness<-rich_ord[match(resB$Order,rich_ord$Order),"Spp_richness"]
rich_fams<-rich[which(rich$Subfamily==""),]; dim(rich_fams)
rich_fams<-rich_fams[which(rich_fams$Family!=""),]; dim(rich_fams)
resB$Family_richness<-rich_fams[match(resB$Family,rich_fams$Family),"Spp_richness"]
rich_subfams<-rich[which(rich$Subfamily!=""),]; dim(rich_subfams)
resB$Subfamily_richness<-rich_subfams[match(resB$Subfamily,rich_subfams$Subfamily),"Spp_richness"]
############# ELIMINATE CROWN AGES FOR FAMILIES WITH SUSPECT CROWN NODES
non<-which(rownames(resB) %in% c( "Aphloiaceae","Brunelliaceae","Cynomoriaceae", 
                                  "Daphniphyllaceae", "Mayacaceae",
                                  "Mitrastemonaceae","Oncothecaceae"))
resB[non,"CG_age_Acal"]<-NA
resB[non,"CG_minHPD"]<-NA
resB[non,"CG_maxHPD"]<-NA

resB$Total_Richness <- rep(NA,dim(resB)[1])
for (i in 1:dim(resB)[1]){ 
  if(!is.na(resB$Subfamily_richness[i])){resB$Total_Richness[i] <- resB$Subfamily_richness[i];next}
  if(!is.na(resB$Family_richness[i])){resB$Total_Richness[i] <- resB$Family_richness[i];next}
  resB$Total_Richness[i] <- resB$Order_richness[i]}
##### GET RID OF OUTGROUPS
out.group<-which(resB$Order%in%c("Pinales","Gnetales","Ginkgoales","Cycadales"))
resB<-resB[-out.group,]
orders<-1:64
fams<-65:500
resB[490:501,]
subfams<-501:dim(resB)[1]
resB_fams<-resB[fams,]

resB$Fuse <- resB$SG_age_Acal - resB$CG_age_Acal
resB_fams$Fuse <- resB_fams$SG_age_Acal - resB_fams$CG_age_Acal
length(which(resB_fams$Fuse >= 25))
length(which(resB_fams$Fuse >= 25)) / length(resB_fams$Fuse)


##############
bb <- boxplot(resB_fams$Fuse);bb
short.fused <- which(resB_fams$Fuse < 15);length(short.fused)
#resB_fams$Family[short.fused][match(famis.big,resB_fams$Family[short.fused])]
#resB_fams$Fuse[short.fused][match(famis.big,resB_fams$Family[short.fused])]
#resB_fams[which(resB_fams$Family%in%famis.big5),]
sum(resB_fams$Total_Richness[short.fused])/sum(resB_fams$Total_Richness)
sum(resB_fams$Total_Richness[short.fused])

mean(resB_fams$CG_age_Acal[short.fused])
#mean(resB_fams$CG_age_Acal[-short.fused],na.rm=T)
mean(resB_fams$CG_age_Acal,na.rm=T)

mean(resB_fams$SG_age_Acal[short.fused])
#mean(resB_fams$SG_age_Acal[-short.fused],na.rm=T)
mean(resB_fams$SG_age_Acal,na.rm=T)

sum(resB_fams$Total_Richness[which(resB_fams$CG_age_Acal < 66)])/sum(resB_fams$Total_Richness)
sum(resB_fams$Total_Richness[which(resB_fams$CG_age_Acal < 66)])
#sum(resB_fams$Total_Richness[which(resB_fams$CG_age_Acal > 66)])/sum(resB_fams$Total_Richness)
#sum(resB_fams$Total_Richness[which(resB_fams$CG_age_Acal > 66)])

#famis.big
#sum(resB_fams$Total_Richness[which(resB_fams$Family%in%famis.big)])/sum(resB_fams$Total_Richness)
bigi <- which(resB_sub$Family %in% famis.big)
#sum(resB_sub$Total_Richness[bigi][which(resB_sub$CG_age_Acal < 66)])/sum(resB_fams$Total_Richness)
resB_hibrid <- resB_fams
resB_hibrid <- resB_hibrid[-which(resB_hibrid$Family%in%famis.big),]
colnames(resB_hibrid);colnames(resB_sub)
resB_hibrid <- rbind(resB_hibrid,resB_sub[,-c(16:18)])
#hist(resB_hibrid$CG_age_Acal,breaks=seq(0,200,10),col=rgb(0,0,1,0.5),xlim=c(0,200))
#hist(resB_hibrid$SG_age_Acal,breaks=seq(0,200,10),col=rgb(1,0,0,0.5),add=T,xlim=c(0,200))
#hist(resB_hibrid$Fuse,breaks=seq(0,160,10),col=rgb(0,1,0,0.5),freq=F)
#hist(resB_fams$Fuse,breaks=seq(0,160,10),col=rgb(0,0,0,0.5),add=T,freq=F)
#mean(resB_hibrid$Fuse,na.rm=T)
#mean(resB_fams$Fuse,na.rm=T)
sum(resB_hibrid$Total_Richness[which(resB_hibrid$CG_age_Acal < 66)])/sum(resB_hibrid$Total_Richness)
sum(resB_hibrid$Total_Richness[which(resB_hibrid$CG_age_Acal < 66)])
