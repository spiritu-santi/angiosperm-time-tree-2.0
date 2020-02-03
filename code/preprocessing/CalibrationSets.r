#############################################################################################
######## calibration.sets()                                                      ############
######## The function takes (1) an input phylogenetic tree (NEWICK), (2) a data  ############
######## frame with a list of species with information on clade (subfamily),     ############
######## family, and order (*.csv), and (3) a list of calibrations as produced   ############                 
######## by PROTEUS (*.csv). Informal clades or higher taxa (above order) are    ############
######## accepted, prompting the user to input two family (or subfamily) names   ############
######## to obtain an MRCA for these calibrations.                               ############
######## The code works directly with the output from PROTEUS.                   ############
######## The function generates a PDF file of the tree with all calibrated       ############
######## nodes depicted with the PROTEUS fossil number, the name of the fossil,  ############
######## and the age used as minimum age.                                        ############
######## Currently there are two ways to generate taxa sets with the             ############
######## 'beast.treepl' option:                                                  ############
######## -- "beast": produce a taxa sets block (NEXUS) with full tip labels for  ############
########     each calibrated clade. This is ready to paste into the NEXUS file   ############
########     for reading in BEAUTI.                                              ############
######## -- "treepl": produces a text file for use in treepl, with alist of two  ############
########     two tip labels (MRCA) for each calibrated clade. This file only     ############
########     needs the final running options implemented in treepl.              ############
#############################################################################################

# 163,114,117,169,177,147,12,183,182,360,336,195,19,251,345,244,346,361,2,215,286,285,288,312,339,124
input= "Data1_RAxML_v.8.0.tre" # available as Supplementary Data in Ramírez-Barahona et al. (2020)
species.table = "Tree_TIPS.csv" # file in data folder
PROTEUS.cal = "PROTEUS_INCLUDED_24-4-18.csv" # available as Supplementary Data in Ramírez-Barahona et al. (2020)
setwd("~/Documents/")
rename.tree.nex <- function (input,species.table,output){
  if(length(grep(".tre",input))==1){
    tre <- read.tree(input)
    tre <- ladderize(tre, TRUE)
    tre$tip.label[which(tre$tip.label=="Lauraceae_Hypodaphniszenkeri")]<-"Lauraceae_Hypodaphnis_zenkeri"
    tre$tip.label[which(tre$tip.label=="Batidaceae_Batis_maritima")]<-"Bataceae_Batis_maritima"
    tre$tip.label[which(tre$tip.label=="Didieraceae_Alluaudia_spp")]<-"Didiereaceae_Alluaudia_spp"
    tre$tip.label[which(tre$tip.label=="Didieraceae_Portulacaria_afra")]<-"Didiereaceae_Portulacaria_afra"
    tre$tip.label[which(tre$tip.label=="Gerrardiniaceae_Gerrardina_foliosa")]<-"Gerrardinaceae_Gerrardina_foliosa"
    tre$tip.label[which(tre$tip.label=="Petenaceae_Petenaea_cordata")]<-"Petenaeaceae_Petenaea_cordata"
    tre$tip.label[which(tre$tip.label=="Petermaniaceae_Petermannia_cirrosa")]<-"Petermanniaceae_Petermannia_cirrosa"
    tre$tip.label[which(tre$tip.label=="Rhipogonaceae_Rhipogonum_elseyanum")]<-"Ripogonaceae_Ripogonum_elseyanum"
    tre$tip.label[which(tre$tip.label=="Rhipogonaceae_Ripogonum_scandens")]<-"Ripogonaceae_Ripogonum_scandens"
    tre$tip.label[which(tre$tip.label=="Rhipogonaceae_Rhipogonum_elseyanum")]<-"Ripogonaceae_Ripogonum_elseyanum"
    tre$tip.label[which(tre$tip.label=="Strelitiziaceae_Ravenala_madagascariensis")]<-"Strelitziaceae_Ravenala_madagascariensis"
    tre$tip.label[which(tre$tip.label=="Strelitiziaceae_Strelitizia_spp")]<-"Strelitziaceae_Strelitzia_spp"
    tre$tip.label[which(tre$tip.label=="Stylidaceae_Donatia_spp")]<-"Stylidiaceae_Donatia_spp"
    tre$tip.label[which(tre$tip.label=="Stylidaceae_Stylidium_spp")]<-"Stylidiaceae_Stylidium_spp"
    list<-tre$tip.label
    list_fams<-foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][1]}
    list_fams_un<-unique(list_fams)
    newlist<-read.table(species.table,sep=",",header=T)
    if(length(which(is.na(match(list,newlist$SPECIES))==TRUE))!=0) {
      cat(red("ERROR: Names in phylogeny do not match names in list !"),"\n")}
    newlist$NEW<-paste(newlist$Order,"_",newlist$SPECIES,"_",newlist$Clade_new,sep="")
    if(length(which((list==newlist$SPECIES[match(list,newlist$SPECIES)])==FALSE))!=0) {
      cat(red("ERROR: Names in phylogeny do not match names in list !"),"\n")}
    list_new<-newlist$NEW[match(list,newlist$SPECIES)]
    tre$tip.label<-as.character(list_new)
    tre$tip.label[grep("_$",tre$tip.label)] <- gsub('.{1}$', '', tre$tip.label[grep("_$",tre$tip.label)])
    write.tree(tre,file=output)}
  if(length(grep(".nex",input))==1){ 
    nexus<-read.nexus.data(input)
    list <- names(nexus)
    if(length(which(is.na(match(list,newlist$SPECIES))==TRUE))!=0) {
      cat(red("ERROR: Names in nexus do not match names in list !"),"\n")}
    newlist$NEW<-paste(newlist$Order,"_",newlist$SPECIES,"_",newlist$Clade_clean,sep="")
    if(length(which((list==newlist$SPECIES[match(list,newlist$SPECIES)])==FALSE))!=0) {
      cat(red("ERROR: Names in phylogeny do not match names in list !"),"\n")}
    list_new<-newlist$NEW[match(list,newlist$SPECIES)]
    names(nexus)<- list_new
    write.nexus.data(nexus,file=output,format="dna",interleaved = F,gap="-",missing="?")
  }
}

calibration.sets<-function(input.tre,species.table,PROTEUS.cal,max_age=154,beast.treepl="beast",output="input_treePL.tre"){
  suffix<-beast.treepl
  cat("Reading phylo tree and changing tip.labels to match species list","\n")
  tre <- read.tree(input.tre)
  tre <- ladderize(tre, TRUE)
  tre$tip.label[which(tre$tip.label=="Lauraceae_Hypodaphniszenkeri")]<-"Lauraceae_Hypodaphnis_zenkeri"
  tre$tip.label[which(tre$tip.label=="Batidaceae_Batis_maritima")]<-"Bataceae_Batis_maritima"
  tre$tip.label[which(tre$tip.label=="Didieraceae_Alluaudia_spp")]<-"Didiereaceae_Alluaudia_spp"
  tre$tip.label[which(tre$tip.label=="Didieraceae_Portulacaria_afra")]<-"Didiereaceae_Portulacaria_afra"
  tre$tip.label[which(tre$tip.label=="Gerrardiniaceae_Gerrardina_foliosa")]<-"Gerrardinaceae_Gerrardina_foliosa"
  tre$tip.label[which(tre$tip.label=="Petenaceae_Petenaea_cordata")]<-"Petenaeaceae_Petenaea_cordata"
  tre$tip.label[which(tre$tip.label=="Petermaniaceae_Petermannia_cirrosa")]<-"Petermanniaceae_Petermannia_cirrosa"
  tre$tip.label[which(tre$tip.label=="Rhipogonaceae_Rhipogonum_elseyanum")]<-"Ripogonaceae_Ripogonum_elseyanum"
  tre$tip.label[which(tre$tip.label=="Rhipogonaceae_Ripogonum_scandens")]<-"Ripogonaceae_Ripogonum_scandens"
  tre$tip.label[which(tre$tip.label=="Rhipogonaceae_Rhipogonum_elseyanum")]<-"Ripogonaceae_Ripogonum_elseyanum"
  tre$tip.label[which(tre$tip.label=="Strelitiziaceae_Ravenala_madagascariensis")]<-"Strelitziaceae_Ravenala_madagascariensis"
  tre$tip.label[which(tre$tip.label=="Strelitiziaceae_Strelitizia_spp")]<-"Strelitziaceae_Strelitzia_spp"
  tre$tip.label[which(tre$tip.label=="Stylidaceae_Donatia_spp")]<-"Stylidiaceae_Donatia_spp"
  tre$tip.label[which(tre$tip.label=="Stylidaceae_Stylidium_spp")]<-"Stylidiaceae_Stylidium_spp"
  list<-tre$tip.label
  list_fams<-foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][1]}
  list_fams_un<-unique(list_fams)
  newlist<-read.table(species.table,sep=",",header=T)
  if(length(which(is.na(match(list,newlist$SPECIES))==TRUE))!=0) {
    cat(red("ERROR: Names in phylogeny do not match names in list !"),"\n")}
  newlist$NEW<-paste(newlist$Order,"_",newlist$SPECIES,"_",newlist$Clade_new,sep="")
  if(length(which((list==newlist$SPECIES[match(list,newlist$SPECIES)])==FALSE))!=0) {
    cat(red("ERROR: Names in phylogeny do not match names in list !"),"\n")}
  list_new<-newlist$NEW[match(list,newlist$SPECIES)]
  tre$tip.label<-as.character(list_new)
  tre$tip.label[grep("_$",tre$tip.label)] <- gsub('.{1}$', '', tre$tip.label[grep("_$",tre$tip.label)])
  write.tree(drop.tip(tre,"_Ophioglossaceae_Ophioglossum_spp",trim.internal = T),file=output)
  proteus<- read.table(PROTEUS.cal,sep=",",header=T,quote="\"");head(proteus);names(proteus)
  proteus$Node.calibrated <- sub("The","",proteus$Node.calibrated)
  proteus$Node.calibrated <- sub("  "," ",proteus$Node.calibrated)
  dima<- dim(proteus)[1]
  if(readline("Exclude calibrations? (y/n):") == "y") { 
    cat("Which Nfos to exclude? (numbers separated by a comma)","\n")
    exc <- readline()
    exc <- unlist(strsplit(exc,","))
    proteus <- proteus[-which(proteus$NFos %in% exc),]
  }
  cat("Excluding:",dima-dim(proteus)[1],"fossils","\n")
  cat("Re-formatting PROTEUS csv output","\n")
  #proteus$Family..of.clade.calibrated. <- rep("",dim(proteus)[1])
  proteus$Clade.calibrated <- rep("",dim(proteus)[1])
  proteus$INFORMAL <- rep("",dim(proteus)[1])
  
  major <- which(proteus$Order..of.clade.calibrated.=="")
  family <- strsplit(as.character(proteus$Node.calibrated),split = " ")
  proteus$Crown.or.stem. <- unlist(lapply(family,"[",1))
  #fams <- unlist(family)[grep("ceae",unlist(family))]
  #proteus$Family..of.clade.calibrated.[grep("ceae",family)] <- fams
  ##fams <- sub("\\(","",fams); fams <- sub(",","",fams)
  subfams <- unlist(lapply(family,"[",2))[grep("ideae",unlist(lapply(family,"[",2)))]
  proteus$Clade.calibrated[match(subfams,unlist(lapply(family,"[",2)))] <- subfams
  duplas<-unique(subfams[duplicated(match(subfams,unlist(lapply(family,"[",2))))])
  for (i in 1:length(duplas)){ 
    proteus$Clade.calibrated[which(unlist(lapply(family,"[",2))==duplas[i])]<-duplas[i]}
  which(grepl(paste(c("core","\\+","neae","polles","clade"), collapse="|"),proteus$Node.calibrated)==TRUE)
  
  iden <- c(major,which(grepl(paste(c("core","\\+","neae","polles","clade"), collapse="|"),proteus$Node.calibrated)==TRUE))
  iden<-unique(iden)
  proteus$INFORMAL[iden] <- "informal"
  proteus$NAME <- paste(unlist(lapply(family,"[",1)),unlist(lapply(family,"[",2)),sep=" ")
  
  iden <- c(major,which(grepl(paste(c("core","\\+","neae","polles","clade","ideae"), collapse="|"),proteus$Node.calibrated)==TRUE))
  subfams <- unlist(lapply(family,"[",2))[-iden]
  proteus$INFORMAL[grep(paste(subfams[-grep(paste(c("ales","ceae"), collapse="|"),subfams)],collapse = "|"),proteus$Node.calibrated)]<-"informal"
  
  
  if (beast.treepl == "treepl") {
    cat("Generating output format for treepl","\n")
    ff<-matrix(ncol=9,nrow=dim(proteus)[1])
    ff[,1] <- "mrca ="; ff[,4] <- "min ="; ff[,7] <- "max =";ff[,9]<-max_age
    colnames(ff)<-c("mrca =","Calibration","tips","Calibration","min =","Min_Age","max =","Calibration","Max_Age")
    head(ff)
    pdf(file=paste(input.tre,"_NODES.pdf",sep=""),width=20,height=180,onefile=F,
        useDingbats=FALSE)
    plot(tre,show.tip.label = T,cex=0.8,no.margin = T,edge.color="grey",
         edge.width = 1,use.edge.length = F,"phylo",node.pos=1,node.depth=2)
    for (i in 1:dim(proteus)[1]) { 
      fos.tax <- str_sub(proteus$Fossil.taxon, 2)
      names <- paste("Nfos.",proteus$NFos,":","\u2020",fos.tax,"(",proteus$Safe.minimum.age,"Mya)")
      names <- sub("  "," ",names);names <- sub("   "," ",names);
      cat(i,"Searching clade for",as.character(proteus$Fossil.taxon[i]),"\n")
      ff[i,2] <- as.character(paste("Nfos",proteus$NFos[i],sep="_"));
      ff[i,5] <- as.character(paste("Nfos",proteus$NFos[i],sep="_"));
      ff[i,8] <- as.character(paste("Nfos",proteus$NFos[i],sep="_"));
      ff[i,6] <- proteus$Safe.minimum.age[i];
      if(proteus$INFORMAL[i]=="informal") {
        cat(green(i,"Specify family (or subfamily) names for MRCA (case sensitive):",as.character(proteus$Fossil.taxon[i]),"\n"))
        cat(green(i,"Calibrated clade",as.character(proteus$Node.calibrated[i]),"\n"))
        names.in<- readline()
        names.in2<- readline()
        names.in<- c(names.in,names.in2)
        cat(green(i,"Specify whether crown (c) or stem (s):",as.character(proteus$Crown.or.stem.[i]),as.character(proteus$Fossil.taxon[i]),"\n"))
        cr.st<- readline()
        tipas<- unique(grep(paste(names.in,collapse="|"),tre$tip.label,value=TRUE));tipas
        #### this one was added on the latest update!!
        if(length(tipas)==1){ancestor <- Ancestors(tre,node=which(tre$tip.label==tipas),"parent")
        subtre<-extract.clade(tre,node=ancestor)
        nodelabels(node=ancestor,frame="none",pch=19,cex=2,col="black")
        nodelabels(node=ancestor,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
        ff[i,3]<-stri_flatten(c(child1,child2)," ") 
        next
        }
        Nnode<- getMRCA(tre,tipas)
        if(cr.st=="c"){
          nodelabels(node=Nnode,frame="none",pch=19,cex=2,col="black")
          nodelabels(node=Nnode,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
          childs<-Children(tre,node=Nnode)
          if(length(Descendants(tre,childs[1])[[1]])==1) {child1<-tre$tip.label[childs[1]]}
          if(length(Descendants(tre,childs[1])[[1]])>1) {child1<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[1])$tip.label[1]]}
          if(length(Descendants(tre,childs[2])[[1]])==1) {child2<-tre$tip.label[childs[2]]}
          if(length(Descendants(tre,childs[2])[[1]])>1) {child2<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[2])$tip.label[1]]}
          ff[i,3]<-stri_flatten(c(child1,child2)," ") 
          next}
        if (cr.st=="s"){  
          ancestor<-Ancestors(tre,node=Nnode,"parent")
          nodelabels(node=ancestor,frame="none",pch=19,cex=2,col="black")
          nodelabels(node=ancestor,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
          childs<-Children(tre,node=ancestor)
          if(length(Descendants(tre,childs[1])[[1]])==1) {child1<-tre$tip.label[childs[1]]}
          if(length(Descendants(tre,childs[1])[[1]])>1) {child1<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[1])$tip.label[1]]}
          if(length(Descendants(tre,childs[2])[[1]])==1) {child2<-tre$tip.label[childs[2]]}
          if(length(Descendants(tre,childs[2])[[1]])>1) {child2<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[2])$tip.label[1]]}
          ff[i,3]<-stri_flatten(c(child1,child2)," ")
          next
        }
      } else { 
        if (proteus$Clade.calibrated[i] != ""){ 
          tipas<- grep(proteus$Clade.calibrated[i],tre$tip.label,fixed=T)
          if(length(tipas)== 0) next
          if (length(tipas)==1) {
            Nnode <- tre$edge[which(tre$edge[,2]==tipas),1] 
            childs<-Children(tre,node=Nnode)
            if(length(Descendants(tre,childs[1])[[1]])==1) {child1<-tre$tip.label[childs[1]]}
            if(length(Descendants(tre,childs[1])[[1]])>1) {child1<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[1])$tip.label[1]]}
            if(length(Descendants(tre,childs[2])[[1]])==1) {child2<-tre$tip.label[childs[2]]}
            if(length(Descendants(tre,childs[2])[[1]])>1) {child2<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[2])$tip.label[1]]}
            ff[i,3]<-stri_flatten(c(child1,child2)," ")
            nodelabels(node=Nnode,frame="none",pch=19,cex=2,col="black")
            nodelabels(node=Nnode,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
            next}
          if(length(tipas) > 1) {Nnode<- getMRCA(tre,tipas)
          if(proteus$Crown.or.stem.[i] == "crown") {
            subtre<- extract.clade(tre,node=Nnode)
            nodelabels(node=Nnode,frame="none",pch=19,cex=2,col="black")
            nodelabels(node=Nnode,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
            childs<-Children(tre,node=Nnode)
            if(length(Descendants(tre,childs[1])[[1]])==1) {child1<-tre$tip.label[childs[1]]}
            if(length(Descendants(tre,childs[1])[[1]])>1) {child1<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[1])$tip.label[1]]}
            if(length(Descendants(tre,childs[2])[[1]])==1) {child2<-tre$tip.label[childs[2]]}
            if(length(Descendants(tre,childs[2])[[1]])>1) {child2<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[2])$tip.label[1]]}
            ff[i,3]<-stri_flatten(c(child1,child2)," ") 
            next}
          if(proteus$Crown.or.stem.[i] == "stem") {
            ancestor<-Ancestors(tre,node=Nnode,"parent")
            nodelabels(node=ancestor,frame="none",pch=19,cex=2,col="black")
            nodelabels(node=ancestor,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
            childs<-Children(tre,node=ancestor)
            if(length(Descendants(tre,childs[1])[[1]])==1) {child1<-tre$tip.label[childs[1]]}
            if(length(Descendants(tre,childs[1])[[1]])>1) {child1<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[1])$tip.label[1]]}
            if(length(Descendants(tre,childs[2])[[1]])==1) {child2<-tre$tip.label[childs[2]]}
            if(length(Descendants(tre,childs[2])[[1]])>1) {child2<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[2])$tip.label[1]]}
            ff[i,3]<-stri_flatten(c(child1,child2)," ")
            next}
          }} 
        else { 
          if (proteus$Family..of.clade.calibrated.[i] != ""){ 
            tipas<- grep(proteus$Family..of.clade.calibrated.[i],tre$tip.label,fixed=T);tipas
            if(length(tipas)== 0) next
            if (length(tipas)==1) {
              Nnode <- tre$edge[which(tre$edge[,2]==tipas),1] 
              childs<-Children(tre,node=Nnode)
              if(length(Descendants(tre,childs[1])[[1]])==1) {child1<-tre$tip.label[childs[1]]}
              if(length(Descendants(tre,childs[1])[[1]])>1) {child1<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[1])$tip.label[1]]}
              if(length(Descendants(tre,childs[2])[[1]])==1) {child2<-tre$tip.label[childs[2]]}
              if(length(Descendants(tre,childs[2])[[1]])>1) {child2<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[2])$tip.label[1]]}
              ff[i,3]<-stri_flatten(c(child1,child2)," ")
              next}
            if(length(tipas)>1) {Nnode<- getMRCA(tre,tipas)
            if(proteus$Crown.or.stem.[i] == "crown") {
              subtre<- extract.clade(tre,node=Nnode)
              nodelabels(node=Nnode,frame="none",pch=19,cex=2,col="black")
              nodelabels(node=Nnode,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
              childs<-Children(tre,node=Nnode)
              if(length(Descendants(tre,childs[1])[[1]])==1) {child1<-tre$tip.label[childs[1]]}
              if(length(Descendants(tre,childs[1])[[1]])>1) {child1<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[1])$tip.label[1]]}
              if(length(Descendants(tre,childs[2])[[1]])==1) {child2<-tre$tip.label[childs[2]]}
              if(length(Descendants(tre,childs[2])[[1]])>1) {child2<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[2])$tip.label[1]]}
              ff[i,3]<-stri_flatten(c(child1,child2)," ")
              next}
            if(proteus$Crown.or.stem.[i] == "stem"){
              ancestor<-Ancestors(tre,node=Nnode,"parent")
              nodelabels(node=ancestor,frame="none",pch=19,cex=2,col="black")
              nodelabels(node=ancestor,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
              childs<-Children(tre,node=ancestor)
              if(length(Descendants(tre,childs[1])[[1]])==1) {child1<-tre$tip.label[childs[1]]}
              if(length(Descendants(tre,childs[1])[[1]])>1) {child1<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[1])$tip.label[1]]}
              if(length(Descendants(tre,childs[2])[[1]])==1) {child2<-tre$tip.label[childs[2]]}
              if(length(Descendants(tre,childs[2])[[1]])>1) {child2<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[2])$tip.label[1]]}
              ff[i,3]<-stri_flatten(c(child1,child2)," ")
              next}
            }}
          else { 
            if (proteus$Order..of.clade.calibrated.[i] != ""){ 
              tipas<- grep(proteus$Order..of.clade.calibrated.[i],tre$tip.label,fixed=T)
              if(length(tipas)== 0) next
              if (length(tipas)==1) {
                Nnode <- tre$edge[which(tre$edge[,2]==tipas),1] 
                childs<-Children(tre,node=Nnode)
                if(length(Descendants(tre,childs[1])[[1]])==1) {child1<-tre$tip.label[childs[1]]}
                if(length(Descendants(tre,childs[1])[[1]])>1) {child1<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[1])$tip.label[1]]}
                if(length(Descendants(tre,childs[2])[[1]])==1) {child2<-tre$tip.label[childs[2]]}
                if(length(Descendants(tre,childs[2])[[1]])>1) {child2<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[2])$tip.label[1]]}
                ff[i,3]<-stri_flatten(c(child1,child2)," ")
                next}
              if(length(tipas)>1) {Nnode<- getMRCA(tre,tipas)
              if(proteus$Crown.or.stem.[i] == "crown") {
                subtre<- extract.clade(tre,node=Nnode)
                nodelabels(node=Nnode,frame="none",pch=19,cex=2,col="black")
                nodelabels(node=Nnode,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
                childs<-Children(tre,node=Nnode)
                if(length(Descendants(tre,childs[1])[[1]])==1) {child1<-tre$tip.label[childs[1]]}
                if(length(Descendants(tre,childs[1])[[1]])>1) {child1<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[1])$tip.label[1]]}
                if(length(Descendants(tre,childs[2])[[1]])==1) {child2<-tre$tip.label[childs[2]]}
                if(length(Descendants(tre,childs[2])[[1]])>1) {child2<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[2])$tip.label[1]]}
                ff[i,3]<-stri_flatten(c(child1,child2)," ")
                next}
              if(proteus$Crown.or.stem.[i] == "stem"){
                ancestor<-Ancestors(tre,node=Nnode,"parent")
                nodelabels(node=ancestor,frame="none",pch=19,cex=2,col="black")
                nodelabels(node=ancestor,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
                childs<-Children(tre,node=ancestor)
                if(length(Descendants(tre,childs[1])[[1]])==1) {child1<-tre$tip.label[childs[1]]}
                if(length(Descendants(tre,childs[1])[[1]])>1) {child1<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[1])$tip.label[1]]}
                if(length(Descendants(tre,childs[2])[[1]])==1) {child2<-tre$tip.label[childs[2]]}
                if(length(Descendants(tre,childs[2])[[1]])>1) {child2<-tre$tip.label[tre$tip.label==extract.clade(tre,node=childs[2])$tip.label[1]]}
                ff[i,3]<-stri_flatten(c(child1,child2)," ")
                next}
              } }
          }
        }}}
    for (i in 1:dim(ff)[1]){
      cat(stri_flatten(c(ff[i,1:3])," "),sep="","\n")
      nodelabels(node=getMRCA(tre,tip = unlist(strsplit(ff[i,3]," ")[1])),frame="none",text=as.character(names[i]),cex=0.5,col="red",adj=1.1)
    }
    sink(paste("Calibration_set_",suffix,".txt",sep=""))
    cat("treefile = input_treePL.tre",sep="","\n")
    cat("smooth = 0.001",sep="","\n")
    cat("numsites = 11395",sep="","\n")
    for (i in 1:dim(ff)[1]){
      nodelabels(node=getMRCA(tre,tip = unlist(strsplit(ff[i,3]," ")[1])),frame="none",text=as.character(names[i]),cex=0.5,col="red",adj=1.1)
      
      cat(stri_flatten(c(ff[i,1:3])," "),sep="","\n")
      cat(stri_flatten(c(ff[i,4:6])," "),sep="","\n")
      cat(stri_flatten(c(ff[i,7:9])," "),sep="","\n")
    }
    cat("outfile = treePL_dated.tre",sep="","\n")
    cat("thorough",sep="","\n")
    cat("opt = 4",sep="","\n")
    cat("moredetail",sep="","\n")
    cat("optad = 4",sep="","\n")
    cat("moredetailad",sep="","\n")
    cat("optcvad = 1",sep="","\n")
    sink()
    dev.off()
    write.table(ff,"Calibrations_raw.csv",sep=",",row.names=F,col.names=T)
  }
  if (beast.treepl == "beast") {
    cat("Generating output format for beast","\n")
    ff<-matrix(ncol=4,nrow=dim(proteus)[1])
    ff[,1] <- "taxset"; ff[,3] <- "="
    colnames(ff)<-c("taxset","Name","=","taxa")
    sub("Data1",input.tre,"Data4") -> out_name
    pdf(file=paste(out_name,"_NODES.pdf",sep=""),width=20,height=180,onefile=F,
        useDingbats=FALSE)
    plot(tre,show.tip.label = T,cex=0.8,no.margin = T,edge.color="grey",
         edge.width = 1,use.edge.length = F,"phylo",node.pos=1,node.depth=2)
    for (i in 1:dim(proteus)[1]) {
      fos.tax <- str_sub(proteus$Fossil.taxon, 2)
      names <- paste("Nfos.",proteus$NFos,":","\u2020",fos.tax,"(",proteus$Safe.minimum.age,"Mya)")
      names <- sub("  "," ",names);names <- sub("   "," ",names)
      cat(i,"Searching clade for",as.character(proteus$Fossil.taxon[i]),"\n")
      ff[i,2] <- as.character(paste("Nfos",proteus$NFos[i],sep="_"));
      if(proteus$INFORMAL[i]=="informal") {
        cat(green(i,"Specify family (or subfamily) names for MRCA (case sensitive):",as.character(proteus$Fossil.taxon[i]),"\n"))
        cat(green(i,"Calibrated clade",as.character(proteus$Node.calibrated[i]),"\n"))
        names.in<- readline()
        names.in2<- readline()
        names.in<- c(names.in,names.in2)
        cat(green(i,"Specify whether crown (c) or stem (s):",as.character(proteus$Crown.or.stem.[i]),as.character(proteus$Fossil.taxon[i]),"\n"))
        cr.st<- readline()
        tipas<- grep(paste(names.in,collapse="|"),tre$tip.label,value=TRUE);tipas
        #### this one was added on the latest update!!
        if(length(tipas)==1){ancestor <- Ancestors(tre,node=which(tre$tip.label==tipas),"parent")
        subtre<-extract.clade(tre,node=ancestor)
        if(proteus$Node.assignment.quality.score.and.method[i] == "5 / phylogenetic analysis") {nodelabels(node=ancestor,frame="none",pch=21,cex=2,bg="white")} else {nodelabels(node=ancestor,frame="none",pch=19,cex=2,col="black")}
        nodelabels(node=ancestor,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
        ff[i,4]<-stri_flatten(subtre$tip.label," ") 
        next
        }
        Nnode<- getMRCA(tre,tipas)
        if(cr.st=="c"){
          subtre<- extract.clade(tre,node=Nnode)
          if(proteus$Node.assignment.quality.score.and.method[i] == "5 / phylogenetic analysis") {nodelabels(node=Nnode,frame="none",pch=21,cex=2,bg="white")} else {nodelabels(node=Nnode,frame="none",pch=19,cex=2,col="black")}
          nodelabels(node=Nnode,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
          ff[i,4]<-stri_flatten(subtre$tip.label," ") 
          next}
        if (cr.st=="s"){  
          ancestor<-Ancestors(tre,node=Nnode,"parent")
          subtre<- extract.clade(tre,node=ancestor)
          if(proteus$Node.assignment.quality.score.and.method[i] == "5 / phylogenetic analysis") {nodelabels(node=ancestor,frame="none",pch=21,cex=2,bg="white")} else {nodelabels(node=ancestor,frame="none",pch=19,cex=2,col="black")}
          nodelabels(node=ancestor,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
          ff[i,4]<-stri_flatten(subtre$tip.label," ") 
          next
        }
      }
      if (proteus$Clade.calibrated[i] != ""){ 
        tipas<- grep(proteus$Clade.calibrated[i],tre$tip.label,fixed=T);tipas
        if(length(tipas) == 0) next
        if (length(tipas) == 1) {Nnode <- tre$edge[which(tre$edge[,2]==tipas),1] 
        subtre<- extract.clade(tre,node=Nnode)
        ff[i,4]<-stri_flatten(subtre$tip.label," ")
        if(proteus$Node.assignment.quality.score.and.method[i] == "5 / phylogenetic analysis") {nodelabels(node=Nnode,frame="none",pch=21,cex=2,bg="white")} else {nodelabels(node=Nnode,frame="none",pch=19,cex=2,col="black")}
        
        nodelabels(node=Nnode,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
        next}
        if(length(tipas) > 1) {Nnode<- getMRCA(tre,tipas)
        if(proteus$Crown.or.stem.[i] == "crown") {
          subtre<- extract.clade(tre,node=Nnode)
          if(proteus$Node.assignment.quality.score.and.method[i] == "5 / phylogenetic analysis") {nodelabels(node=Nnode,frame="none",pch=21,cex=2,bg="white")} else {nodelabels(node=Nnode,frame="none",pch=19,cex=2,col="black")}
          
          nodelabels(node=Nnode,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
          ff[i,4]<-stri_flatten(subtre$tip.label," ") 
          next}
        if(proteus$Crown.or.stem.[i] == "stem") {
          ancestor<-Ancestors(tre,node=Nnode,"parent")
          if(proteus$Node.assignment.quality.score.and.method[i] == "5 / phylogenetic analysis") {nodelabels(node=ancestor,frame="none",pch=21,cex=2,bg="white")} else {nodelabels(node=ancestor,frame="none",pch=19,cex=2,col="black")}
          
          nodelabels(node=ancestor,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
          subtre<- extract.clade(tre,node=ancestor)
          ff[i,4]<-stri_flatten(subtre$tip.label," ") 
          next}}
      } else { 
        if (proteus$Family..of.clade.calibrated.[i] != ""){ 
          tipas<- grep(proteus$Family..of.clade.calibrated.[i],tre$tip.label,fixed=T);tipas
          if(length(tipas)== 0) next
          if (length(tipas)==1) {Nnode <- tre$edge[which(tre$edge[,2]==tipas),1] 
          subtre<- extract.clade(tre,node=Nnode)
          ff[i,4]<-stri_flatten(subtre$tip.label," ")
          if(proteus$Node.assignment.quality.score.and.method[i] == "5 / phylogenetic analysis") {nodelabels(node=Nnode,frame="none",pch=21,cex=2,bg="white")} else {nodelabels(node=Nnode,frame="none",pch=19,cex=2,col="black")}
          nodelabels(node=Nnode,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
          next}
          if(length(tipas)>1) {Nnode<- getMRCA(tre,tipas)
          if(proteus$Crown.or.stem.[i] == "crown") {
            subtre<- extract.clade(tre,node=Nnode)
            if(proteus$Node.assignment.quality.score.and.method[i] == "5 / phylogenetic analysis") {nodelabels(node=Nnode,frame="none",pch=21,cex=2,bg="white")} else {nodelabels(node=Nnode,frame="none",pch=19,cex=2,col="black")}
            nodelabels(node=Nnode,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
            ff[i,4]<-stri_flatten(subtre$tip.label," ") 
            next}
          if(proteus$Crown.or.stem.[i] == "stem"){
            ancestor<-Ancestors(tre,node=Nnode,"parent")
            subtre<- extract.clade(tre,node=ancestor)
            if(proteus$Node.assignment.quality.score.and.method[i] == "5 / phylogenetic analysis") {nodelabels(node=ancestor,frame="none",pch=21,cex=2,bg="white")} else {nodelabels(node=ancestor,frame="none",pch=19,cex=2,col="black")}
            nodelabels(node=ancestor,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
            ff[i,4]<-stri_flatten(subtre$tip.label," ") 
            next}}
        } else { 
          if (proteus$Order..of.clade.calibrated.[i] != ""){ 
            tipas<- grep(proteus$Order..of.clade.calibrated.[i],tre$tip.label,fixed=T);tipas
            if(length(tipas) == 0) next
            if (length(tipas) == 1) {Nnode <- tre$edge[which(tre$edge[,2]==tipas),1] 
            subtre<- extract.clade(tre,node=Nnode)
            ff[i,4]<-stri_flatten(subtre$tip.label," ")
            if(proteus$Node.assignment.quality.score.and.method[i] == "5 / phylogenetic analysis") {nodelabels(node=Nnode,frame="none",pch=21,cex=2,bg="white")} else {nodelabels(node=Nnode,frame="none",pch=19,cex=2,col="black")}
            nodelabels(node=Nnode,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
            next}
            if(length(tipas)>1) {Nnode<- getMRCA(tre,tipas)
            if(proteus$Crown.or.stem.[i] == "crown") {
              subtre<- extract.clade(tre,node=Nnode)
              if(proteus$Node.assignment.quality.score.and.method[i] == "5 / phylogenetic analysis") {nodelabels(node=Nnode,frame="none",pch=21,cex=2,bg="white")} else {nodelabels(node=Nnode,frame="none",pch=19,cex=2,col="black")}
              nodelabels(node=Nnode,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
              ff[i,4]<-stri_flatten(subtre$tip.label," ") 
              next}
            if(proteus$Crown.or.stem.[i] == "stem"){
              ancestor<-Ancestors(tre,node=Nnode,"parent")
              subtre<- extract.clade(tre,node=ancestor)
              if(proteus$Node.assignment.quality.score.and.method[i] == "5 / phylogenetic analysis") {nodelabels(node=ancestor,frame="none",pch=21,cex=2,bg="white")} else {nodelabels(node=ancestor,frame="none",pch=19,cex=2,col="black")}
              nodelabels(node=ancestor,frame="none",text=as.character(names[i]),cex=0.5,col="black",adj=1.1)
              ff[i,4]<-stri_flatten(subtre$tip.label," ")
              next}}
          } 
        }
      }}
    for (i in 1:dim(ff)[1]){
      cat(i,names[i],sep=" ","--")
      nodelabels(node=getMRCA(tre,tip = unlist(strsplit(ff[i,4]," ")[1])),frame="none",text=as.character(names[i]),cex=0.5,col="red",adj=1.1)
      cat("OK","\n")
    }
    sink(paste("Calibration_set_",suffix,".txt",sep=""))
    cat("begin sets;","\n")
    for (i in 1:dim(ff)[1]){ 
      cat(stri_flatten(c(ff[i,1:4])," "),"\n")
    }
    sink()
    dev.off()
  }}
