tre<-phyloch::read.beast(tree)
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
tre <- ladderize(tre,TRUE); list<-tre$tip.label
list_ords <- foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][1]}
list_fams <- foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][2]}
list_subfams <- foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][5]}
list_ords_un<-sort(unique(list_ords))
list_fams_un<-sort(unique(list_fams))
list_subfams_un<-sort(unique(list_subfams))
list_un <- (c(list_ords_un,list_fams_un,list_subfams_un))
list_table<-cbind(list_ords,list_fams,list_subfams);list_table<-list_table[!duplicated(list_table),]
res<-matrix(nrow=length(list_fams_un)+length(list_ords_un)+length(list_subfams_un),ncol=10)
rownames(res)<-list_un; colnames(res)<-c("CG_age_Acal","SG_age_Acal","No_tips","CG_minHPD","CG_maxHPD","SG_minHPD","SG_maxHPD","Order","Family","Subfamily")
for(i in 1:length(list_un)){
  iden<-grep(list_un[i],x = tre$tip.label)
  res[i,3]<-length(iden)
  cat(i,"----",list_un[i],"\n")
  if(length(iden) == 1) {
    res[i,1]<-NA
    ssNODE <-tre$edge[which(tre$edge[,2]==iden),1]
    subtrea<-extract.clade(tre,ssNODE)
    res[i,2]<-max(branching.times(subtrea))
    stopifnot(round(tre$height_median[ssNODE-1215],digits = 2) == round(res[i,2],digits = 2))
    res[i,6]<-tre$`height_95%_HPD_MIN`[ssNODE-1215];res[i,7]<-tre$`height_95%_HPD_MAX`[ssNODE-1215]
  }
  else{ 
    subtre<-extract.clade(tre,getMRCA(tre,tre$tip.label[iden]))
    if(length(subtre$tip.label)-length(iden)!=0) {cat("-----------","non-monophyletic","\n"); next}
    res[i,3]<-length(iden)
    res[i,1]<-max(branching.times(subtre))
    NODE <- getMRCA(tre,tre$tip.label[iden])
    stopifnot(round(tre$height_median[NODE-1215],digits = 2)==round(res[i,1],digits = 2))
    res[i,4]<-tre$`height_95%_HPD_MIN`[NODE-1215];res[i,5]<-tre$`height_95%_HPD_MAX`[NODE-1215]
    sNODE <-tre$edge[which(tre$edge[,2]==NODE),1]
    subtre<-extract.clade(tre,sNODE)
    res[i,2]<-max(branching.times(subtre))
    stopifnot(round(tre$height_median[sNODE-1215],digits = 2)==round(res[i,2],digits = 2))
    res[i,6]<-tre$`height_95%_HPD_MIN`[sNODE-1215];res[i,7]<-tre$`height_95%_HPD_MAX`[sNODE-1215]
  }}
for (i in 1:length(list_un)){
  if(i <= 68){res[i,8] <- as.character(list_un[i])}
  if(i > 68 & i<= 510) {res[i,8]<-as.character(list_table[match(rownames(res),list_table[,2])[i],1]); res[i,9]<-as.character(list_table[match(rownames(res),list_table[,2])[i],2])}
  if(i > 510) {res[i,8]<-as.character(list_table[match(rownames(res),list_table[,3])[i],1]); res[i,9]<-as.character(list_table[match(rownames(res),list_table[,3])[i],2]); res[i,10]<-as.character(list_table[match(rownames(res),list_table[,3])[i],3])
  }
}
res<-as.data.frame(res)
write.tree(tre,file=paste(ruta_write,"0.MCC_parsed.tree",sep=""))
write.table(res,file=paste(ruta_write,"2.Ages_complete.csv",sep=""),sep=",",row.names=T,col.names=T,quote = F)
