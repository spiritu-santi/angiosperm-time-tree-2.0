resB<-read.table(paste(ruta_write,"2.Ages_complete.csv",sep=""),sep=",",header=T)
rich<-read.table("SPP.RICHNESS.csv",heade=T,sep=",") #### File in data folder
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
