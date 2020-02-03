#
#
#
#
# Uses additional files

rich<-read.table("SPP.RICHNESS.csv",heade=T,sep=",") # file in data folder
tree = "~/Documents/NATURE/RESULTADOS_V4/TREES/SDstrict_treePL_dated.tre" # file in data folder (the six trees available in a ZIP file).

tre <- read.tree(tree)
##### CORRECT SOME TIP NAMES, ALTHOUGH NOT NECESSARY AT THIS POINT
tre$tip.label[which(tre$tip.label=="Celastrales_Celastraceae_Euonymus_spp_Celastraceae")]<-"Celastrales_Celastraceae_Euonymus_spp_Celastraceae.sstr"
tre$tip.label[which(tre$tip.label=="Laurales_Lauraceae_Hypodaphniszenkeri_Hypodaphnideae")]<-"Laurales_Lauraceae_Hypodaphnis_zenkeri_Hypodaphnideae"
tre$tip.label[which(tre$tip.label=="Brassicales_Batidaceae_Batis_maritima")]<-"Brassicales_Bataceae_Batis_maritima"
tre$tip.label[which(tre$tip.label=="Caryophyllales_Didieraceae_Alluaudia_spp")]<-"Caryophyllales_Didiereaceae_Alluaudia_spp"
tre$tip.label[which(tre$tip.label=="Caryophyllales_Didieraceae_Portulacaria_afra")]<-"Caryophyllales_Didiereaceae_Portulacaria_afra"
tre$tip.label[which(tre$tip.label=="Huerteales_Gerrardiniaceae_Gerrardina_foliosa")]<-"Huerteales_Gerrardinaceae_Gerrardina_foliosa"
tre$tip.label[which(tre$tip.label=="Huerteales_Petenaceae_Petenaea_cordata")]<-"Huerteales_Petenaeaceae_Petenaea_cordata"
tre$tip.label[which(tre$tip.label=="Liliales_Petermaniaceae_Petermannia_cirrosa")]<-"Liliales_Petermanniaceae_Petermannia_cirrosa"
tre$tip.label[which(tre$tip.label=="Liliales_Rhipogonaceae_Rhipogonum_elseyanum")]<-"Liliales_Ripogonaceae_Ripogonum_elseyanum"
tre$tip.label[which(tre$tip.label=="Liliales_Rhipogonaceae_Ripogonum_scandens")]<-"Liliales_Ripogonaceae_Ripogonum_scandens"
tre$tip.label[which(tre$tip.label=="Liliales_Rhipogonaceae_Rhipogonum_elseyanum")]<-"Liliales_Ripogonaceae_Ripogonum_elseyanum"
tre$tip.label[which(tre$tip.label=="Zingiberales_Strelitiziaceae_Ravenala_madagascariensis")]<-"Zingiberales_Strelitziaceae_Ravenala_madagascariensis"
tre$tip.label[which(tre$tip.label=="Zingiberales_Strelitiziaceae_Strelitizia_spp")]<-"Zingiberales_Strelitziaceae_Strelitzia_spp"
tre$tip.label[which(tre$tip.label=="Asterales_Stylidaceae_Donatia_spp")]<-"Asterales_Stylidiaceae_Donatia_spp"
tre$tip.label[which(tre$tip.label=="Asterales_Stylidaceae_Stylidium_spp")]<-"Asterales_Stylidiaceae_Stylidium_spp"
tre <- ladderize(tre,TRUE)
list<-tre$tip.label
list_ords <- foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][1]}
list_fams <- foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][2]}
list_subfams <- foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][5]}
list_ords_un<-sort(unique(list_ords));length(list_ords_un)
list_fams_un<-sort(unique(list_fams));length(list_fams_un)
list_subfams_un<-sort(unique(list_subfams));length(list_subfams_un)
list_un <- (c(list_ords_un,list_fams_un,list_subfams_un))
list_table<-cbind(list_ords,list_fams,list_subfams);list_table<-list_table[!duplicated(list_table),]
res<-matrix(nrow=length(list_fams_un)+length(list_ords_un)+length(list_subfams_un),ncol=6)
rownames(res)<-list_un; colnames(res)<-c("CG_age_Acal","SG_age_Acal","No_tips","Order","Family","Subfamily")
for(i in 1:length(list_un)){
  iden<-grep(list_un[i],x = tre$tip.label)
  res[i,3]<-length(iden)
  cat(i,"----",list_un[i],"\n")
  if(length(iden) == 1) {
    res[i,1]<-NA
    ssNODE <-tre$edge[which(tre$edge[,2]==iden),1]
    subtrea<-extract.clade(tre,ssNODE)
    res[i,2]<-max(branching.times(subtrea))
  }
  else{ 
    subtre<-extract.clade(tre,getMRCA(tre,tre$tip.label[iden]))
    if(length(subtre$tip.label)-length(iden)!=0) {cat("-----------","non-monophyletic","\n"); next}
    res[i,3]<-length(iden)
    res[i,1]<-max(branching.times(subtre))
    NODE <- getMRCA(tre,tre$tip.label[iden])
    subtre<-extract.clade(tre,tre$edge[which(tre$edge[,2]==NODE),1])
    res[i,2]<-max(branching.times(subtre))
  }}
for (i in 1:length(list_un)){
  if(i <= 68){res[i,4] <- as.character(list_un[i])}
  if(i > 68 & i<= 510) {res[i,4]<-as.character(list_table[match(rownames(res),list_table[,2])[i],1]); res[i,5]<-as.character(list_table[match(rownames(res),list_table[,2])[i],2])}
  if(i > 510) {res[i,4]<-as.character(list_table[match(rownames(res),list_table[,3])[i],1]); res[i,5]<-as.character(list_table[match(rownames(res),list_table[,3])[i],2]); res[i,6]<-as.character(list_table[match(rownames(res),list_table[,3])[i],3])
  }
}
res<-as.data.frame(res)
write.table(res,file=paste(ruta_write,2.1.Ages_complete_treePL.csv,sep=""),sep=",",row.names=T,col.names=T,quote = F)

file_beast = "RES_SDcomplete_v2/2.Ages_complete.csv" # generated in previous scripts.
file_treepl = "RES_SDcomplete_v2/2.1.Ages_complete_treePL.csv"

resBEAST<-read.table(file_beast,sep=",",header=T)
resTREEPL<-read.table(file_treepl,sep=",",header=T)
rich_ord<-rich[which(rich$Family==""),]; dim(rich_ord)
resBEAST$Order_richness<-rich_ord[match(resBEAST$Order,rich_ord$Order),"Spp_richness"]
rich_fams<-rich[which(rich$Subfamily==""),]; dim(rich_fams)
rich_fams<-rich_fams[which(rich_fams$Family!=""),]; dim(rich_fams)
resBEAST$Family_richness<-rich_fams[match(resBEAST$Family,rich_fams$Family),"Spp_richness"]
rich_subfams<-rich[which(rich$Subfamily!=""),]; dim(rich_subfams)
resBEAST$Subfamily_richness<-rich_subfams[match(resBEAST$Subfamily,rich_subfams$Subfamily),"Spp_richness"]
non<-which(rownames(resBEAST) %in% c("Aphloiaceae", "Daphniphyllaceae", 
                                      "Brunelliaceae", "Cynomoriaceae", "Mayacaceae", 
                                      "Mitrastemonaceae", "Oncothecaceae"))
resBEAST[non,"CG_age_Acal"]<-NA
resBEAST$Total_Richness <- rep(NA,dim(resBEAST)[1])
for (i in 1:dim(resBEAST)[1]){ 
  if(!is.na(resBEAST$Subfamily_richness[i])){resBEAST$Total_Richness[i] <- resBEAST$Subfamily_richness[i];next}
  if(!is.na(resBEAST$Family_richness[i])){resBEAST$Total_Richness[i] <- resBEAST$Family_richness[i];next}
  resBEAST$Total_Richness[i] <- resBEAST$Order_richness[i]}
out.group<-which(resBEAST$Order%in%c("Pinales","Gnetales","Ginkgoales","Cycadales"))
resBEAST<-resBEAST[-out.group,]
resTREEPL$Order_richness<-rich_ord[match(resTREEPL$Order,rich_ord$Order),"Spp_richness"]
rich_fams<-rich[which(rich$Subfamily==""),]; dim(rich_fams)
rich_fams<-rich_fams[which(rich_fams$Family!=""),]; dim(rich_fams)
resTREEPL$Family_richness<-rich_fams[match(resTREEPL$Family,rich_fams$Family),"Spp_richness"]
rich_subfams<-rich[which(rich$Subfamily!=""),]; dim(rich_subfams)
resTREEPL$Subfamily_richness<-rich_subfams[match(resTREEPL$Subfamily,rich_subfams$Subfamily),"Spp_richness"]
non<-which(rownames(resTREEPL) %in% c("Aphloiaceae", "Daphniphyllaceae", 
                                   "Brunelliaceae", "Cynomoriaceae", "Mayacaceae", 
                                   "Mitrastemonaceae", "Oncothecaceae"))
resTREEPL[non,"CG_age_Acal"]<-NA
resTREEPL$Total_Richness <- rep(NA,dim(resTREEPL)[1])
for (i in 1:dim(resTREEPL)[1]){ 
  if(!is.na(resTREEPL$Subfamily_richness[i])){resTREEPL$Total_Richness[i] <- resTREEPL$Subfamily_richness[i];next}
  if(!is.na(resTREEPL$Family_richness[i])){resTREEPL$Total_Richness[i] <- resTREEPL$Family_richness[i];next}
  resTREEPL$Total_Richness[i] <- resTREEPL$Order_richness[i]}
out.group<-which(resTREEPL$Order %in% c("Pinales","Gnetales","Ginkgoales","Cycadales"))
resTREEPL<-resTREEPL[-out.group,]

orders<-1:64
fams<-65:500
subfams<-501:dim(resTREEPL)[1]
print(unique(rownames(resTREEPL)==rownames(resBEAST)))
resTREEPL<-resTREEPL[fams,]
resBEAST<-resBEAST[fams,]

pdf("1.ComparisonMethods_UC.pdf")
plot(resTREEPL$SG_age_Acal,resBEAST$SG_age_Acal,type="n",ylim=c(0,200),xlim=c(0,200),
     xlab="Stem ages (treePL)",ylab="Stem ages (BEAST)");abline(a=0,b=1)
title("Stem age estimates (families)")
points(resBEAST$SG_age_Acal,resTREEPL$SG_age_Acal,pch=21,cex=1,bg=rgb(1,0,0,0.8))
plot(resBEAST$CG_age_Acal,resTREEPL$CG_age_Acal,type="n",ylim=c(0,200),xlim=c(0,200),
     xlab="Crown ages (complete set)",ylab="Crown ages (phylo-only set)");abline(a=0,b=1)
title("Crown age estimates (families)")
points(resBEAST$CG_age_Acal,resTREEPL$CG_age_Acal,pch=21,cex=1,bg=rgb(1,0,0,0.8))
dev.off()
