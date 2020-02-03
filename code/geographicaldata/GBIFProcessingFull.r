######################################################################################################
######################################################################################################
#####         ALL THE FOLLOWING IS NOT AUTOMATED (last update: JUNE 30, 2019)                    #####
##### The following script does all the necessary processing for the geographic data as obtained #####
##### from GBIF. There is a pre-processing step in bash to reduce file size, but row ID can be   #####
##### traces back to the original data as downloaded from GBIF with extra information.           #####
##### The script includes the necessary code to perform all the analyses and figures presented   #####
##### in Ramírez-Barahona et al. (2019) Tracing the evolutionary history of flowering  plants    #####
##### through time and space.                                                                    #####
######################################################################################################
######################################################################################################
# Preprocessing: setting up data files retrieved from GBIF
# Use bash commands to reduce file size by selecting fewer columns
# Original GBIF data download for angiosperms available here: doi.org/10.15468/dl.iftsjs
# New GBIF data download for vascular plants available here: 
# 'DISTASnew.csv' is the name given to the downloaded file
spiritusanti$ awk -F"\t" '{print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$17"\t"$18}’ DISTASnew.csv > distas_COR.txt
spiritusanti$ awk -F"\t" '{print $10"\t"$28}' DISTASnew.csv > distas_COR_dates.txt
spiritusanti$ wc -l DISTASnew.csv
spiritusanti$ wc -l distas_COR.txt

# Prepare file in R to run python script of Edwards et al. (2018) (available here: https://www.nature.com/articles/nature14393?draft=collection) 
# R commands
py <- fread("distas_COR.txt");dim(py);colnames(py)
py <- cbind(py[,7:9],c(1:dim(py)[1]))
py[which(py$species!=""),]->py; dim(py)
write.table(py,"DISTAS_python.csv",col.names = F,row.names = F,quote = F,sep=",")

# First Step: clean the database using the python script of Edwards et al.
spiritusanti$ python cleanGbifCoords.1.0.py DISTAS_python.csv allHerbaria_ADM1_badCoords.txt DISTAS_good.csv DISTAS_bad.csv

# Second step: eliminate unidentified records (species level)
a <- fread("DISTAS_good.csv")
bad <- fread("DISTAS_bad.csv")
aa <- fread("distas_COR.txt")
iden<-as.vector(a$V4); b <- aa[iden,]
dim(b); dim(aa)
aa<-b; aa$species<-sub(" ","_",aa$species)

# Third step: name filtering and homogeneization
list_namesA <- fread("TPL_all_genera.csv") # list of genera in The Plant List
list_names <- list_namesA;dim(list_names)
names(list_names)[11] <- "Taxonomic_status_TPL"
list_names$Complete_NAME <- paste(list_names$Genus,list_names$Species,sep="_")
list_names$New_NAME_ACCEPTED <- list_names$Complete_NAME[match(list_names$`Accepted ID`,list_names$ID,nomatch = NA)]
list_names[is.na(list_names$New_NAME_ACCEPTED)]$New_NAME_ACCEPTED<-list_names[is.na(list_names$New_NAME_ACCEPTED)]$Complete_NAME
match(aa$species,list_names$Complete_NAME)->iden
aa$species_ACCEPTED<-list_names$New_NAME_ACCEPTED[iden]
fams<-read.table("APW_synonyms.txt",sep=",") # list of family names in the Angiosperm Phylogeny Website
colnames(fams)<-c("Name","Synonym","Order")
fams$ACCEPTED<-rep(NA,dim(fams)[1])
for(i in 1:dim(fams)[1]) {if((fams$Synonym[i]=="")==TRUE){fams$ACCEPTED[i]<-as.character(fams$Name[i])}
                          else (fams$ACCEPTED[i]<-as.character(fams$Synonym[i]))}
match(aa$family,fams$Name)->iden
aa$family_ACCEPTED<-fams$ACCEPTED[iden]






###### FOURTH STEP: MINOR EDITS AND SAVE ########
length(unique(aa$family_ACCEPTED))
sort(unique(aa$family_ACCEPTED))
aa[which(aa$species=="Xyris_araracuare"),"species_ACCEPTED"]<-"Xyris_araracuarae"
aa[which(aa$species=="Xyris_araracuare"),"family_ACCEPTED"]<-"Xyridaceae"
aa[which(aa$genus=="Afrothismia"),"family_ACCEPTED"]<-"Afrothismieae"
aa[which(aa$family==""),]
##### REMOVE SOME FOSSILS
aa<-aa[-c(which(aa$genus=="Liliacidites")),]
aa<-aa[-c(which(aa$species=="Lindernia_confusa")),]
aa<-aa[-c(which(aa$genus=="Thalassotaenia")),]
aa[which(aa$family==""),]
cat("Number of families (RAW):",length(unique(aa$family)),"\n")
cat("Number of families (CORRECTED):",length(unique(aa$family_ACCEPTED)),"\n")
cat("Number of species (RAW):",length(unique(aa$species)),"\n")
cat("Number of species (CORRECTED):",length(unique(aa$species_ACCEPTED)),"\n")
unique(aa$class)
save(aa,file="DISTAS_corNAMESangios.r") ##### DATA BASE
write.table(sort(unique(aa$family_ACCEPTED)),file="FamAngios_corNAMES.txt",row.names=F,quote=F)

###### FIFTH STEP: FAMILY RE-CIRCUMSCRIPTION TO MATCH CLASSIFICATION USED IN ANGIOS TIME-TREE 2.0 ########
load("DISTAS_corNAMESangios.r")
tipas<- read.table("TIPAS_v.7.1.csv",sep=",",header=T); fa_ti<-unique(tipas$Family)
fams <- sort(unique(aa$family_ACCEPTED))
fa_ti[which(is.na(match(fa_ti,fams)))] ######## IDENTIFY FAMILIES IN PHYLOGENY THAT ARE NOT IN THE DATA BASE
####### THIS IS A BIT OF A NASTY CODE..... NEED TO HAVE A FILE WITH THE NAMES
length(unique(aa$family_ACCEPTED))
sort(unique(aa$family_ACCEPTED))

######## TACCACEA
aa[grep("Tacca_",aa$species_ACCEPTED),"family_ACCEPTED"]<-"Taccaceae"
length(unique(aa$family_ACCEPTED))

######## CODONACEAE + COLDENIACEAE + HOPLESTIGMATACEAE
aa[grep("Codon_",aa$species_ACCEPTED),"family_ACCEPTED"]<-"Codonaceae"
aa[grep("Coldenia_",aa$species_ACCEPTED),"family_ACCEPTED"]<-"Coldeniaceae"
aa[grep("Hoplestigma_",aa$species_ACCEPTED),"family_ACCEPTED"]<-"Hoplestigmataceae"

####### CORDIACEAE
names.in <- c("Cerdana_","Cienkowskya_","Cordia_","Catonia_","Carpiphea_","Calyptracordia_",
              "Bourgia_","Auxemma_","Ascania_","Acnadena_","Varroniopsis_","Varronia_","Ulmarronia_",
              "Toquera_","Topiaris_","Sebestena_","Salimori_","Saccellium","Quarena_","Rhabdocalyx_",
              "Plethostephia_","Piloisia_","Pilicordia_","Physoclada_","Patagonula_","Paradigmia_",
              "Novella_","Myxa_","Montjolya_","Macria_","Maciela_","Lithocardium_","Hymenesthes_",
              "Hemigymnia_","Gynaion_","Gerascanthus_","Firensia_","Ectemis_","Diacoria_","Cordiopsis_",
              "Cordiada_","Coilanthera_","Collococcus_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED),"family_ACCEPTED"]<-"Cordiaceae"
length(unique(aa$family_ACCEPTED))
sort(unique(aa$family_ACCEPTED))

###### EHRETIACEAE
names.in <- c("Ammobroma_","Antrophora_","Beurreria_","Bourreria_","Carmona_","Corallophyllum_","Cortesia_",
              "Crematomia_","Desmophyla_","Diplostylus_","Eddya_","Ehretia_","Galapagoa_","Gaza_","Halgania_",
              "Hilsenbergia_","Lennoa_","Lepidocordia_","Lithothamnus_","Lutrostylis_","Menais_","Monomesia_",
              "Morelosia_","Ptilocalyx_","Pholisma_","Rhabdia_","Rhabdia_","Rotula_","Stegnocarpus_","Subrisia_",
              "Tetracoccus_","Tiquilia_","Tiquiliopsis_","Traxilum_","Zombiana_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED),"family_ACCEPTED"]<-"Ehretiaceae"
length(unique(aa$family_ACCEPTED))

####### WELLSTEDIACEAE
aa[grep("Wellstedia_",aa$species_ACCEPTED),"family_ACCEPTED"]<-"Wellstediaceae"

###### NAMACEAE
names.in <- c("Andropus_","Conanthus_","Eriodictyon_","Ernstamra_","Lemmonia_","Marilaunidium_","Nama_",
              "Turricula_","Wigandia_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED),"family_ACCEPTED"]<-"Namaceae"
length(unique(aa$family_ACCEPTED))

###### HELIOTROPIACEAE
names.in <- c("Argusia_","Beruniella_","Bourjotia_","Bucanion_","Ceballosia_","Cochranea_",
              "Dialion_","Eliopia_","Euploca_","Heliophytum_","Heliotropium_","Hieranthemum_","Hilgeria_",
              "Ixorhea_","Lithococca_","Mallotonia_","Meladendron_","Messerschmidia_","Messersmidia_","Myriopus_",
              "Nogalia_","Oskampia_","Parabouchetia_","Peristema_","Piptoclaina_","Pittonia_",
              "Schleidenia_","Schobera_","Scorpianthes_","Scorpiurus_","Spilocarpus_","Synzistachium_","Tetrandra_",
              "Tiaridium_","Tournefortia_","Valentina_","Valentiniella_","Verrucaria_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED),"family_ACCEPTED"]<-"Heliotropiaceae"
length(unique(aa$family_ACCEPTED))

######## MAUNDIACEAE
aa[grep("Maundia_",aa$species_ACCEPTED),"family_ACCEPTED"]<-"Maundiaceae"

####### MYODOCARPACEAE
names.in <- c("Delarbrea_","Myodocarpus_","Porospermum_","Pseudosciadium_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED),"family_ACCEPTED"]<-"Myodocarpaceae"
length(unique(aa$family_ACCEPTED))

####### MYSTROPETALACEAE
names.in <- c("Dactylanthus_","Hachettea_","Mystropetalon_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED),"family_ACCEPTED"]<-"Mystropetalaceae"
length(unique(aa$family_ACCEPTED))

###### PELTANTHERACEAE
aa[grep("Peltanthera_",aa$species_ACCEPTED),"family_ACCEPTED"]<-"Peltantheraceae"

###### PETIVERIACEA
names.in <- c("Flueckigera_","Gallesia_","Hilleria_","Ledenbergia_","Mohlana_","Monococcus_",
              "Petiveria_","Rivina_","Schindleria_","Seguieria_","Trichostigma_","Villamilla_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED),"family_ACCEPTED"]<-"Petiveriaceae"
length(unique(aa$family_ACCEPTED))

###### THISMIACEAE + AFROTHISMIA
names.in <- c("Bagnisia_","Cymbocarpa_","Geomitra_","Glaziocharis_","Haplothismia_","Mamorea_",
              "Myostoma_","Ophiomeris_","Oxygyne_","Saionia_","Sarcosiphon_","Scaphiophora_","Thismia_",
              "Tiputinia_","Tribrachys_","Triscyphus_","Triurocodon_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED),"family_ACCEPTED"]<-"Thismiaceae"
aa[grep("Afrothismia_",aa$species_ACCEPTED),"family_ACCEPTED"]<-"Afrothismieae"

######## NYSSACEAE
names.in <- c("Camptotheca_","Davidia_","Diplopanax_","Mastixia_","Nyssa_")
aa[grep(paste(names.in,collapse = "|"),aa$species_ACCEPTED),"family_ACCEPTED"]<-"Nyssaceae"
length(unique(aa$family_ACCEPTED))

####### LAST EDITS
unique(aa[which(is.na(aa$family_ACCEPTED)),"family"])
aa[grep("Stixaceae",aa$family),"family_ACCEPTED"]<-"Resedaceae"
aa[grep("Vivianiaceae",aa$family),"family_ACCEPTED"]<-"Francoaceae"
length(unique(aa$family_ACCEPTED))
sort(unique(aa$family_ACCEPTED))

cat("Number of families (RAW):",length(unique(aa$family)),"\n")
cat("Number of families (CORRECTED):",length(unique(aa$family_ACCEPTED)),"\n")
cat("Number of species (RAW):",length(unique(aa$species)),"\n")
cat("Number of species (CORRECTED):",length(unique(aa$species_ACCEPTED)),"\n")

##### FINAL CHECK AGAINST TREE
tre<-phyloch::read.beast("~/Documents/NATURE/PAPER/v_5/SUPPS/Data9_Rcode/RC_complete_MCCv_2")
##### SPELLING ERRORS IN TREE
{tre$tip.label[which(tre$tip.label=="Celastrales_Celastraceae_Euonymus_spp_Celastraceae")]<-"Celastrales_Celastraceae_Euonymus_spp_Celastraceae.sstr"
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
  tre$tip.label[which(tre$tip.label=="Asterales_Stylidaceae_Stylidium_spp")]<-"Asterales_Stylidiaceae_Stylidium_spp"}
list<-tre$tip.label
list_fams <- foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][2]}
list_fams_un<-sort(unique(list_fams))
match(list_fams_un,sort(unique(aa$family_ACCEPTED)))->matcha
list_fams_un[which(is.na(matcha))]


match(sort(unique(aa$family_ACCEPTED)),list_fams_un)->matcha
sort(unique(aa$family_ACCEPTED))[which(is.na(matcha))]
save(aa,file="DISTAS_angios_Corrected.r") ########## SAVE FILE

###### SIXTH STEP: GEOGRAPHICAL OUTLIER DETECTION BY SPECIES ########
setwd("~/Desktop/DISTS")
load("DISTAS_angios_Corrected.r"); p_clean<-as_tibble(aa);p_clean
names(p_clean)[8:9] <- c("decimallatitude","decimallongitude")
p_clean <- p_clean[-c(which(is.na(p_clean$decimallatitude)),which(is.na(p_clean$decimallongitude))),]
p_clean$IDs <- 1:dim(p_clean)[1]
aa_na <- p_clean[which(is.na(p_clean$species_ACCEPTED)),] ######## SPECIES WITH NO ACCEPTED NAMES IN TPL: 'NAs' IN TABLE
familia <-  sort(unique(p_clean$species_ACCEPTED))
p_clean$Correct <- NA
fami <- (unique(p_clean$family_ACCEPTED))
data("buffland")
#### FIRST THE BASIC TESTS
p_clean
#### THIS CAN BE QUICKER.....
for (i in 1:length(fami)){ 
  fama <- p_clean[which(p_clean$family_ACCEPTED == fami[i]),]
  cat(toupper(fami[i]),"--",i,"--","\n")
  outas_1 <- clean_coordinates(x=fama,lon="decimallongitude",lat="decimallatitude",
                               test=c("duplicates","centroids","capitals","gbif","seas","institutions","equal","zeros"),
                               outliers_method = "distance",inst_rad=0.05,seas_ref=buffland,
                               capitals_rad=0.05,centroids_rad = 0.05,zeros_rad = 0.2,
                               species="species_ACCEPTED",value="flagged",verbose=T)
  p_clean[which(p_clean$family_ACCEPTED == fami[i]),"Correct"] <- outas_1 
}
table(p_clean$Correct)
p_clean_v4 <- p_clean
save(p_clean_v4,file = "DISTAS_angios_Corrected.r")

###### EXTRA OUTLIER CORRECTION ON TWO FAMILIES
load("DISTAS_angios_Corrected.r")
p_clean <- p_clean_v4
table(p_clean$Correct)
p_clean[which(p_clean$family_ACCEPTED=="Tetracarpaeaceae"),"species_ACCEPTED"]<- "Tetracarpaea_tasmannica"
p_clean[which(p_clean$family_ACCEPTED=="Kewaceae"),"species_ACCEPTED"] <- p_clean[which(p_clean$family_ACCEPTED=="Kewaceae"),"species"]
fami <- c("Kewaceae","Tetracarpaeaceae")
for (i in 1:length(fami)){ 
  fama <- p_clean[which(p_clean$family_ACCEPTED == fami[i]),]
  cat(toupper(fami[i]),"--",i,"--","\n")
  outas_1 <- clean_coordinates(x=fama,lon="decimallongitude",lat="decimallatitude",
                                         test=c("duplicates","centroids","capitals","gbif","seas","institutions","equal","zeros"),
                                         outliers_method = "distance",inst_rad=0.05,seas_ref=buffland,
                                         capitals_rad=0.05,centroids_rad = 0.05,zeros_rad = 0.2,
                                         species="species_ACCEPTED",value="flagged",verbose=T)
  p_clean[which(p_clean$family_ACCEPTED == fami[i]),"Correct"] <- outas_1 
}
dim(p_clean);table(p_clean$Correct)

jpeg("Data_raw.jpg",height=900,width=900)
plot(p_clean$decimallongitude,p_clean$decimallatitude,pch=19,cex=0.4,col=rgb(0,0,0,0.3))
maps::map("world",interior=F,add=T)
dev.off()
jpeg("Data_cleaned.jpg",height=900,width=900)
f <- which(p_clean$Correct==TRUE)
plot(p_clean$decimallongitude[f],p_clean$decimallatitude[f],pch=19,cex=0.4,col=rgb(0,0,0,0.3))
maps::map("world",interior=F,add=T)
dev.off()


###### ITERATIVE NAME MATCHING USING JARO-WINKLER DISTANCES
spp_na <- as.data.frame(unique(p_clean[is.na(p_clean$species_ACCEPTED),"species"]))
sort(spp_na[,1])-> spp_na
list_namesA<-fread("FAMILIAS/TPL_all_genera.csv")
list_names<-list_namesA
names(list_names)[11]<-"Taxonomic_status_TPL"
list_names$Complete_NAME<-paste(list_names$Genus,list_names$Species,sep="_")
list_names$New_NAME_ACCEPTED<-list_names$Complete_NAME[match(list_names$`Accepted ID`,list_names$ID,nomatch = NA)]
list_names[is.na(list_names$New_NAME_ACCEPTED)]$New_NAME_ACCEPTED<-list_names[is.na(list_names$New_NAME_ACCEPTED)]$Complete_NAME
spp_full <- as.data.frame(list_names[,c("Complete_NAME","New_NAME_ACCEPTED")])
spp_full <- spp_full[!duplicated(spp_full[,1]),]
genera <- unlist(lapply(strsplit(spp_full[,1],"_"),"[[",1))
spp <- unlist(lapply(strsplit(spp_full[,1],"_"),"[[",2))
spp_full <- cbind(spp_full,genera,spp)
ID_rows_2 <- c();needy<-c(); matched_spp <- c()
for(i in 1:length(spp_na)){
cat(i,"Jaro-Winkler matching on:",as.character(spp_na[i]),"\n")
genera_dists <- jarowinkler(strsplit(spp_na[i],"_")[[1]][1],as.character(spp_full$genera))
if(length(which(genera_dists == 1))==0) {cat(red("No match for genus!"),"\n"); next}
spp_dists <- jarowinkler(strsplit(spp_na[i],"_")[[1]][2],as.character(spp_full$spp)[which(genera_dists == 1)])
if(max(spp_dists) < 0.94) {
  p_clean[which(p_clean$species==spp_na[i]),"species_ACCEPTED"] <- NA
  cat(red("No match"),"\n")
  needy <- c(needy,as.character(spp_na[i]))
  next
}
if(max(spp_dists) >= 0.94) {
  if(length(which(spp_dists==max(spp_dists)))>1) {
    cat(red(bold("Found more than two matches","\n")))
    next}
cat(blue("Found possible match:",paste(as.character(spp_full$genera)[which(genera_dists == 1)][1],as.character(spp_full$spp)[which(genera_dists == 1)][which(spp_dists==max(spp_dists))],sep="_")),"\n")
ID_rows_2 <- c(ID_rows_2, p_clean[which(p_clean$species==spp_na[i]),"IDs"])
p_clean[which(p_clean$species==spp_na[i]),"species_ACCEPTED"] <- paste(as.character(spp_full$genera)[which(genera_dists == 1)][1],as.character(spp_full$spp)[which(genera_dists == 1)][which(spp_dists==max(spp_dists))],sep="_")
matched_spp<-c(matched_spp,as.character(paste(as.character(spp_full$genera)[which(genera_dists == 1)][1],as.character(spp_full$spp)[which(genera_dists == 1)][which(spp_dists==max(spp_dists))],sep="_")))
}
}
spp_nonmatched <- as.data.frame(unique(p_clean[is.na(p_clean$species_ACCEPTED),"species"]))
length(matched_spp);dim(spp_nonmatched);length(spp_na)
p_clean_bu <- p_clean

###### OUTLIER DETECTION ON MATCHED SPECIES
species_ma <-  unique(matched_spp)
fama <- p_clean[which(p_clean$species_ACCEPTED %in% species_ma),]
outas_1 <- clean_coordinates(x=fama,lon="decimallongitude",lat="decimallatitude",
                                       test=c("duplicates","centroids","capitals","gbif","seas","institutions","equal","zeros"),
                                       outliers_method = "distance",inst_rad=0.05,seas_ref=buffland,
                                       capitals_rad=0.05,centroids_rad = 0.05,zeros_rad = 0.2,
                                       species="species_ACCEPTED",value="flagged",verbose=T)
p_clean[which(p_clean$species_ACCEPTED %in% species_ma),"Correct"] <- outas_1
table(p_clean$Correct)

jpeg("Cleaned.jpg",height=900,width=900)
f<-which(p_clean$Correct==TRUE)
plot(p_clean$decimallongitude[f],p_clean$decimallatitude[f],pch=19,cex=0.4,col=rgb(0,0,0,0.4))
maps::map("world",interior=F,add=T)
dev.off()

save(p_clean,file="DISTAS_angios_Corrected_2.r")
######## STATS
p_clean_2 <- p_clean
length(matched_spp);dim(spp_nonmatched);length(spp_na) ##### NO. OF MISPELLED SPECIES: MATCHED, NO MATCH, TOTAL 
length(unique(p_clean_2$species_ACCEPTED)) #### TOTAL ACCEPTED SPECIES
length(unique(p_clean_2$species)) #### TOTAL UNCORRECTED SPECIES
length(unique(p_clean_2$family_ACCEPTED)) #### TOTAL ACCEPTED FAMILIES
length(unique(p_clean_2$family)) #### TOTAL UNCORRECTED FAMILIES
p_clean_2[which(p_clean_2$Correct==TRUE),] #### TOTAL CORRECT RECORDS
p_clean_2[which(p_clean_2$Correct==FALSE),] #### TOTAL OUTLIER RECORDS


#### SECOND THE DISTANCE TEST
load("DISTAS_angios_Corrected_2.r")
p_clean_2 <- p_clean[p_clean$Correct,]
table(p_clean_2$Correct)
dim(p_clean);dim(p_clean_2)
#which(p_clean_2$species_ACCEPTED=="Zea_mays") -> idon;length(idon)
fami <- (unique(p_clean_2$family_ACCEPTED))
p_clean_2
#### THIS CAN BE QUICKER.....
for (i in 1:length(fami)){ 
  fama <- p_clean_2[which(p_clean_2$family_ACCEPTED == fami[i]),]
  cat(toupper(fami[i]),"--",i,"--","\n")
  outas_1 <- clean_coordinates(x=fama,lon="decimallongitude",lat="decimallatitude",
                                         test=c("outliers"),outliers_td=500,
                                         species="species_ACCEPTED",value="flagged",verbose=T)
  p_clean_2[which(p_clean_2$family_ACCEPTED == fami[i]),"Correct"] <- outas_1
}
p_clean_2[which(p_clean_2$Correct==TRUE),] #### TOTAL CORRECT RECORDS
p_clean_2[which(p_clean_2$Correct==FALSE),] #### TOTAL OUTLIER RECORDS
save(p_clean_2,file="DISTAS_angios_Corrected_2.1.r")
###### EXTRA STEP: MINOR CORRECTIONS AND CHECK ########
load("DISTAS_angios_Corrected_2.1.r")
which(is.na(p_clean_2$species_ACCEPTED))
p_clean_2 <- p_clean_2[-which(is.na(p_clean_2$species_ACCEPTED)),];p_clean_2
length(unique(p_clean_2$species_ACCEPTED)) #### TOTAL ACCEPTED SPECIES
length(unique(p_clean_2$family_ACCEPTED)) #### TOTAL ACCEPTED FAMILIES
##### FINAL CHECK AGAINST TREE
  tre<-phyloch::read.beast("~/Documents/NATURE/PAPER/v_5/FINAL_files/DATA_supp/BEAST_mcctrees/RC_complete_MCCv_2.tre")
##### SPELLING ERRORS IN TREE
{tre$tip.label[which(tre$tip.label=="Celastrales_Celastraceae_Euonymus_spp_Celastraceae")]<-"Celastrales_Celastraceae_Euonymus_spp_Celastraceae.sstr"
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
  tre$tip.label[which(tre$tip.label=="Asterales_Stylidaceae_Stylidium_spp")]<-"Asterales_Stylidiaceae_Stylidium_spp"}
list<-tre$tip.label
list_fams <- foreach(i=1:length(list),.combine=c) %dopar% {strsplit(list,"_")[[i]][2]}
list_fams_un<-sort(unique(list_fams))
match(list_fams_un,sort(unique(p_clean_2$family_ACCEPTED)))->matcha
list_fams_un[which(is.na(matcha))]
match(sort(unique(p_clean_2$family_ACCEPTED)),list_fams_un)->matcha
sort(unique(p_clean_2$family_ACCEPTED))[which(is.na(matcha))]

p_clean_2 <- p_clean_2[which(p_clean_2$Correct==TRUE),]; p_clean_2
xy<-as.data.frame(p_clean_2[,c("decimallongitude","decimallatitude")]) ##### EXCLUDING ERRONEOUS RECORDS
spp<-as.data.frame(p_clean_2[,c("species_ACCEPTED","family_ACCEPTED")])
dim(xy);dim(spp)
p_clean_2$genus_ACCEPTED<-unlist(lapply(strsplit(p_clean_2$species_ACCEPTED,"_"),"[",1))
a<-read.table("FINAL_DISTAS/Correction_needed.csv",sep=",");head(a)
for (i in 1:dim(a)[1]){
p_clean_2[which(p_clean_2$genus_ACCEPTED==a[i,2]),"family_ACCEPTED"]<-as.character(a[i,3])}
genera <- unique(p_clean_2$genus_ACCEPTED); length(genera)
genus_fam<- paste(p_clean_2$genus_ACCEPTED,p_clean_2$family_ACCEPTED,sep="_")
p_clean_2$GENUS_FAM<-genus_fam
genera_fam<-(unique(genus_fam))
length(genera);length(genera_fam)
length(unique(p_clean_2$species_ACCEPTED)) #### TOTAL ACCEPTED SPECIES
length(unique(p_clean_2$family_ACCEPTED)) #### TOTAL ACCEPTED FAMILIES
save(p_clean_2,file="DISTAS_angios_Corrected_3.r")

###### EXTRA STEP: FLAG RECORDS USING THE GLONAF DATABASE #####
###### READING THE TAXON LIST AND DEFINING A VECTOR OF ALIEN SPECIES
table1 <- fread("~/Downloads/GLONAF/Taxon_x_List_GloNAF_vanKleunenetal2018Ecology.txt",header=T,sep="\t",fill=T)
alien_taxa <- unique(table1$standardized_name);length(alien_taxa)
alien_taxa<-sub(" ","_",alien_taxa)
alien_taxa <- unlist(lapply(strsplit(alien_taxa," "),"[",1))
###### READING THE REGION ID LIST
table2 <- fread("~/Downloads/GLONAF/Region_GloNAF_vanKleunenetal2018Ecology.txt",header=T,sep="\t",fill=T)
glonaf <- readOGR("~/Downloads/GLONAF/257_2_GloNAF_Shapefile/regions2.shp")
###### LOADING THE GEOGRAPHIC OCCURRENCE DATASET
load("DISTAS_angios_Corrected_3.r")
p_clean_2$flagged_glonaf <- rep(0,dim(p_clean_2)[1])
###### IDENTIFY SUSPECT RECORDS FOR ALIEN SPECIES
flagged_records<-NULL
for (i in 1:length(alien_taxa)){
  cat("Extracting alien taxa",i,":",alien_taxa[i],"\n")
  idon <- unique(which(p_clean_2$species_ACCEPTED == alien_taxa[i]))
  if(length(idon==0)){idon<-unique(which(p_clean_2$species == alien_taxa[i]))}
  test <- p_clean_2[idon,]
  if(dim(test)[1]==0) {next}
  test[,9:8] -> ssp
  pointos<-SpatialPoints(as.data.frame(ssp[,1:2]))
  proj4string(pointos)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  glonaf_points<- sp::over(pointos,glonaf)
  regions_alien <- unique(table1[which(table1$standardized_name == unique(table1$standardized_name)[i]),"region_id"])
  poly_alien <- table2[which(table2$region_id %in% regions_alien$region_id),"OBJIDsic"]
  flagged <- which(glonaf_points$OBJIDsic %in% poly_alien$OBJIDsic)
  if(length(flagged)==0){next}
  flagged_records<-c(flagged_records,test[["IDs"]][flagged])
}
flagged_records;length(flagged_records)
#### FLAG AND 
p_clean_2[match(flagged_records,p_clean_2$IDs),"flagged_glonaf"]<-1
p_clean_glonaf <- p_clean_2
table(p_clean_glonaf$Correct)
table(p_clean_glonaf$flagged_glonaf)
maps::map("world",interior=F,add=F)
antartica <- which(p_clean_glonaf$decimallatitude< -60)
p_clean_glonaf[antartica,"flagged_glonaf"] <- 1
jpeg("Cleaned_last.jpg",height=900,width=900)
f<-which(p_clean_glonaf$flagged_glonaf==0)
plot(p_clean_glonaf$decimallongitude[f],p_clean_glonaf$decimallatitude[f],pch=19,cex=0.4,col=rgb(0,0,0,0.4))
maps::map("world",interior=F,add=T)
dev.off()
save(p_clean_glonaf,file="DISTAS_angios_Corrected_4.r")


###### SEVENTH STEP: DO PROPORTIONAL SPECIES RICHNESS PAMs ########
###### RELOAD FROM SAVE
setwd("~/Desktop/DISTS")
load("DISTAS_angios_Corrected_4.r")
p_clean_2 <- p_clean_glonaf[which(p_clean_glonaf$flagged_glonaf==0),]
cat("Estimating grid-cells prPAMs","\n")
prPAM <- seq(1:16200)
familias <- sort(unique(p_clean_2$family_ACCEPTED))
tot <- length(familias)
p_clean_2$genus_ACCEPTED<-unlist(lapply(strsplit(p_clean_2$species_ACCEPTED,"_"),"[",1))
###### prPAM-grid FOR FAMILIES
length(unique(p_clean_2$species_ACCEPTED)); length(unique(p_clean_2$family_ACCEPTED))
for (i in 1:length(familias)) {
  fama <- p_clean_2[which(p_clean_2$family_ACCEPTED==familias[i]),]
  fama<-as.data.frame(fama)
  cat("PAMs for:",green(toupper(familias[i])),"\n")
  pp <- lets.presab.points(xy=fama[,c(9,8)],species=fama[,10],xmn=-180,xmx = 180,ymn = -90,ymx = 90,resol = 2,remove.sp = F,
      remove.cells = F,show.matrix = F,count=T)
  save(pp,file=paste("PAMs/",familias[i],"_Grid_PAM.r",sep = ""))
  if(dim(pp$Presence_and_Absence_Matrix)[2]==3){prPAM <- cbind(prPAM,pp$Presence_and_Absence_Matrix[,3])}
  if(dim(pp$Presence_and_Absence_Matrix)[2]>3){prPAM <- cbind(prPAM,rowSums(pp$Presence_and_Absence_Matrix[,-c(1:2)]))}
}
head(prPAM);dim(prPAM)
colnames(prPAM)[-1]<-familias
prPAM_g <- as_tibble(prPAM[,-1])
rowSums(prPAM_g)
prPAM_g<-add_column(prPAM_g,LAT=pp$Presence_and_Absence_Matrix[,2],.before = 1)
prPAM_g<-add_column(prPAM_g,LON=pp$Presence_and_Absence_Matrix[,1],.before = 1)
prPAM_g
save(prPAM_g,file="prPAM_FAM_grid.r")
###### prPAM-grid FOR GENERA
p_clean_2$GENUS_FAM<-genus_fam
un_genus_fam<-unique(genus_fam)
prPAM_genera <- seq(1:16200)
count<-length(un_genus_fam)-1;sum<-0
genera_2<-c()
for (i in 1:length(un_genus_fam)) {
  sum<-sum+1; count<-count - 1 
  if(sum==1) {tic("Time taken each 100 genera")}
  fama <- p_clean_2[which(p_clean_2$GENUS_FAM==genus_fam[i]),]
  fama<-as.data.frame(fama)
  cat("PAM for:",green(toupper(un_genus_fam[i])),"\n")
  pp <- lets.presab.points(xy=fama[,c(9,8)],species=fama[,10],xmn=-180,xmx = 180,ymn = -90,ymx = 90,resol = 2,remove.sp = F,
                              remove.cells = F,show.matrix = F,count=F)
  if(dim(pp$Presence_and_Absence_Matrix)[2] == 3){
    prPAM_genera <- cbind(prPAM_genera,pp$Presence_and_Absence_Matrix[,3])
    genera_2<-c(genera_2,un_genus_fam[i])
  }
  if(dim(pp$Presence_and_Absence_Matrix)[2] > 3){
    prPAM_genera <- cbind(prPAM_genera,rowSums(pp$Presence_and_Absence_Matrix[,-c(1:2)]))}
  genera_2<-c(genera_2,un_genus_fam[i])
  if(sum==100) {toc(log=T)
    sum<-0; secs <- tic.log(format = F)[[1]]$toc- tic.log(format = F)[[1]]$tic
    cat(red(count,"species to go --"),((count/100)*secs)/3600,"hrs to go (approx)","\n")
    tic.clearlog()}
}
colnames(prPAM_genera)[-1]<-un_genus_fam
prPAM_genera <- as.tibble(prPAM_genera[,-1])
prPAM_genera<-add_column(prPAM_genera,LAT=pp$Presence_and_Absence_Matrix[,2],.before = 1)
prPAM_genera<-add_column(prPAM_genera,LON=pp$Presence_and_Absence_Matrix[,1],.before = 1)
prPAM_genera
save(prPAM_genera,file="prPAM_GEN_grid.r")

###### prPAM-eco FOR FAMILIES AND GENERA
poly<-readOGR("~/Desktop/official/wwf_terr_ecos.shp")
ecos_area <- cbind(as.character(poly$ECO_NAME),poly$area_km2,as.character(poly$REALM),as.character(poly$BIOME))
ecos_area <- ecos_area[!duplicated(ecos_area),] ###### eco regions with area size
sort(unique(poly$ECO_NAME))
lista_eco <- paste(ecos_area[,1],ecos_area[,3],sep="_")
eco_PAM<- data.frame(Ecoregions=as.character(lista_eco))
cat("Estimating ecoRegions prPAMs","\n")
familias <- sort(unique(p_clean_2$family_ACCEPTED))
tot <- length(familias)
length(unique(p_clean_2$species_ACCEPTED)); length(unique(p_clean_2$family_ACCEPTED))
p_clean_2 <- cbind(p_clean_2,rep("-",dim(p_clean_2)[1]))
p_clean_2 <- cbind(p_clean_2,rep("-",dim(p_clean_2)[1]))
p_clean_2 <- cbind(p_clean_2,rep("-",dim(p_clean_2)[1]))
names(p_clean_2)
names(p_clean_2)[16:18]<- c("eco_NAME","AREA","REALM")
prPAM <- seq(1:dim(eco_PAM)[1])
prPAM_genus<-seq(1:dim(eco_PAM)[1])
genera<-c()
for (i in 1:length(familias)) {
  fama <- p_clean_2[which(p_clean_2$family_ACCEPTED==familias[i]),]
  cat("PAMs for:",green(toupper(familias[i])),"\n")
    ssp <- SpatialPoints(fama[,9:8],proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    pre <- over(ssp,poly,returnList = F)
    fama$eco_NAME <- paste(pre$ECO_NAME,pre$REALM,sep="_")
    fama$AREA <- pre$AREA
    fama$REALM <- pre$REALM
    pp <- table(fama$eco_NAME,fama$species_ACCEPTED)
    save(pp,file=paste(familias[i],"PAM.r",sep = "_"))
    fama <- fama[!duplicated(fama[,c(10,16)]),]
    pg<-table(fama$eco_NAME,fama$GENUS_FAM)
    prPAM_genus<-cbind(prPAM_genus,pg[match(lista_eco,rownames(pg)),])
    genera<-c(genera,colnames(pg))
    p <- table(fama$eco_NAME,fama$family_ACCEPTED)
    prPAM <- cbind(prPAM,p[match(lista_eco,rownames(p))])
    cat(unique(sort(unique(fama$GENUS_FAM))==colnames(pg)),"\n")
}
head(prPAM);dim(prPAM)
prPAM_2 <- as.tibble(prPAM)
prPAM_2[is.na(prPAM_2)]<-0
prPAM_2[,1]<-as.character(lista_eco)
colnames(prPAM_2) <- c("ECO_NAME",familias)
prPAM_2;max(rowSums(prPAM_2[,-1]))
save(prPAM_2,file="prPAM_FAM_ECO.r")
head(prPAM_genus);dim(prPAM_genus)
prPAM_genus <- as.tibble(prPAM_genus)
prPAM_genus[is.na(prPAM_genus)]<-0
prPAM_genus[,1]<-as.character(lista_eco)
colnames(prPAM_genus)<-c("ECO_NAME",genera)
dim(prPAM_genus);max(rowSums(prPAM_genus[,-1]))
save(prPAM_genus,file="prPAM_genus_ECO.r")

###### EIGTH STEP: ANALYSES ON ECO-PAM #########
setwd("~/Desktop/ANGIOSPERMS/DISTS./FINAL_DISTAS")
load("prPAM_FAM_ECO.r")
# prPAM_2<-prPAM_genus
poly<-readOGR("~/Desktop/official/wwf_terr_ecos.shp")
ecos_area <- cbind(as.character(poly$ECO_NAME),poly$area_km2,as.character(poly$REALM),as.character(poly$BIOME))
ecos_area <- ecos_area[!duplicated(ecos_area),] ###### eco regions with area size
unique(prPAM_2[,1]==paste(ecos_area[,1],ecos_area[,3],sep="_"))
max(rowSums(prPAM_2[,-1]))
which(rowSums(prPAM_2[,-1]) <= 100) -> iden_inf
ecos2 <- ecos_area[-iden_inf,]
prPAM_3<- prPAM_2[-iden_inf,]
unique(prPAM_3[,1]==paste(ecos2[,1],ecos2[,3],sep="_"))
which(is.na(ecos2[,3]))
Lake <- which(ecos2[,1]=="Lake")
RockIce<- which(ecos2[,1]=="Rock and Ice")
AN_NA<-c(which(is.na(ecos2[,3])),which(ecos2[,3]=="AN"))
biomes <- (as.vector(ecos2[-c(AN_NA,Lake,RockIce),3]));biomes_COR <- (as.vector(ecos2[-c(AN_NA,Lake,RockIce),3]))
biomes2<-biomes
prPAM_hell <- decostand(prPAM_3[-c(AN_NA,Lake,RockIce),-1],method="hell",MARGIN=1) ####### DOING EUCLIDEAN WITH HELLINGER TRANSFORMED
pcoA <- cmdscale(vegdist(prPAM_hell,"euc"),k = 2);plot(pcoA)
prPAM_4<-add_column(prPAM_3,rep(1,dim(prPAM_3)[1]))
prPAM_hell2 <- decostand(prPAM_3[-c(AN_NA,Lake,RockIce),-1],method="total") ###### DOING BRAY-CURTIS WITH RELATIVE ABUNDANCES
pcoA <- cmdscale(vegdist(prPAM_hell2,"bray"),k = 2);plot(pcoA)
###### REARRANGE REALMS AND DEFINE COLORS
iden<-c(which(biomes2=="AA"))
ecos3<-ecos2[-c(AN_NA,Lake,RockIce),]
pcoA[iden,];sort(as.data.frame(prPAM_3[iden,1])[,1])
ecos3[c(which(biomes2=="AA")),1]
aa_pcoA<-pcoA[iden,]
kmeans(aa_pcoA,2)$cluster->clusts
plot(aa_pcoA,xlab="First PCo",ylab="Second PCo",
     xlim=c(-0.7,0.7),ylim=c(-.5,.5))
points(aa_pcoA,bg=clusts,pch=21)
iden2 <- which(clusts==1)
ecos_AA <- ecos3[iden,1]
ecos_AA[iden2];ecos3[iden,3]
points(aa_pcoA[iden2,],bg="blue",pch=21)
plot(subset(poly,REALM == "AA"),lwd=0.1,xlim=c(60,180))
plot(subset(poly,ECO_NAME %in% ecos_AA[iden2]),add=T,lwd=0.1,col="red")
biomes_COR[which(ecos3[,1] %in% ecos_AA[iden2])] <- "MA"
unique(biomes_COR)
ecos_AA<-ecos_AA[iden2]

iden<-c(which(biomes2=="NT"))
pcoA[iden,];sort(as.data.frame(prPAM_3[iden,1])[,1])
ecos3[c(which(biomes2=="NT")),1]
aa_pcoA<-pcoA[iden,]
kmeans(aa_pcoA,centers=2,nstart = 10)$cluster->clusts
plot(aa_pcoA,xlab="First PCo",ylab="Second PCo",
     xlim=c(-0.7,0.7),ylim=c(-.5,.5))
points(aa_pcoA,bg=clusts,pch=21)
iden2 <- which(clusts==1)
ecos_NT <- ecos3[iden,1]
ecos_NT[iden2]
points(aa_pcoA[iden2,],bg="blue",pch=21)
plot(subset(poly,REALM == "NT"),lwd=0.1,xlim=c(-130,-30))
plot(subset(poly,ECO_NAME %in% ecos_NT[iden2]),add=T,lwd=0.1,col="red")
biomes_COR[which(ecos3[,1] %in% ecos_NT[iden2])] <- "NTZ"
unique(biomes_COR)
ecos_NT<-ecos_NT[iden2]
biomes_COR2<-biomes_COR
unique(biomes_COR);unique(biomes_COR2)
colors <- c("springgreen4","yellow3","lightblue2","blue","tomato1",
            "gold","firebrick3","mediumpurple1",
            "darkorchid2")


for (i in 1:length(unique(biomes_COR))){ 
  biomes_COR[biomes_COR==unique(biomes_COR)[i]]<-colors[i]
}
biomes_COR

o<-cbind(unique(biomes_COR2),colors)


###### PREPARE DATA FOR LEFSE
library(foreach)
prPAM_lefse<-prPAM_2
realms <- foreach(i=1:830,.combine=c) %do% {
  strsplit(as.character(prPAM_lefse[i,1]),"_")[[1]][2]
}
ecos <- foreach(i=1:830,.combine=c) %do% {
  strsplit(as.character(prPAM_lefse[i,1]),"_")[[1]][1]
}
realms[which(ecos%in%ecos_AA)]<-"MA"
realms[which(ecos%in%ecos_NT)]<-"NTS"
table(realms)

max(rowSums(prPAM_lefse[,-1]))
which(rowSums(prPAM_lefse[,-c(1)]) <= 100) -> iden
prPAM_lefse<-prPAM_lefse[-iden,]
which(colSums(prPAM_lefse[,-1])==0)->identa;length(identa)
identa<-identa+1
prPAM_lefse<-prPAM_lefse[,-identa]
prPAM_lefse[,-1] <- decostand(prPAM_lefse[,-1],method="total",MARGIN=1)
dim(ecos_area[-iden,]);dim(prPAM_lefse)
biome <- ecos_area[,4]
#prPAM_lefse<-add_column(biome[-iden],.data=prPAM_lefse,.before = T)
prPAM_lefse<-add_column(realms[-iden],.data=prPAM_lefse,.before = T)
#colnames(prPAM_lefse)<-unlist(lapply(strsplit(colnames(prPAM_lefse),"_"),"[",1))
t_prPAM_lefse<-t(prPAM_lefse)
t_prPAM_lefse[-c(1:2),]<-as.numeric(t_prPAM_lefse[-c(1:2),])
dim(t_prPAM_lefse)
rownames(t_prPAM_lefse)
write.table(t_prPAM_lefse,"~/Desktop/LEFSE/prPAM_FAM_ECO.txt",sep="\t",row.names=T,col.names=F,quote = F)

###### PLOT PCoA BY REALM
re <- match(levels(as.factor(biomes_COR2)),o[,1])
unique(biomes_COR2)
pdf("Ordination_EcoREgions_GENUS.pdf",width=13,height=9)
par(mai=c(1,1,1,2.5));par(xpd=T)
plot(pcoA,type="n",xlab="First PCo",ylab="Second PCo",
     xlim=c(-0.6,0.8),ylim=c(-.6,.6))
points(pcoA,pch=21,cex=0.6,lwd=0.4,bg=biomes_COR)
dataEllipse(x=pcoA,groups=as.factor(biomes_COR2),levels=c(0.9),
            center.pch=F, draw=TRUE, plot.points=F,
robust=T,col=colors[re],lwd=3, fill=FALSE,group.labels = NULL)
title("Eco-regions by propSR (Bray-Curtis)")
colors2 <- c("springgreen4","yellow3","gold","tomato1","firebrick3","mediumpurple1",
            "darkorchid2","lightblue2","blue")
names <- c("Neotropics","Neo-transitional","Afrotropics","Indo-Malayan",
           "Melanesia","Pacific-Oceanic","Australasia","Palearctic","Nearctic")
legend(x=0.9,y=-.095,legend=c(names),pch=21,
       pt.bg=c(colors2),horiz = F)
text(x=1.3,y=-.05,labels=c("BIOGEOGRAPHIC REALMS"),font=2)
dev.off()

colors

rainbow <- colorRampPalette(c("white","gold", "firebrick1"), space ="rgb")(10)
eco.dist <- vegdist(prPAM_hell,"euc")
eco.cust <- hclust(eco.dist,"aver")
heatmap(x=as.matrix(prPAM_hell),Rowv= as.dendrogram(eco.cust),
        Colv=NA,col = rainbow, margins = c(5, 5))

###### DOING ORDINATION ON FAMILIES
#prPAM_5<-as_tibble(t(prPAM_3[,-1]))
#prPAM_5<-add_column(prPAM_5,dummy=rep(1,436))
#prPAM_hell <- decostand(prPAM_5,method="total")
#pcoA <- cmdscale(vegdist(prPAM_hell,"raup"),k = 2);plot(pcoA)
#fams <- names(prPAM_3)[-1]
#ssp_Rich<-read.table("SPP.RICHNESS.csv",sep=",",header=T)
#head(ssp_Rich)
#as.vector(ssp_Rich$Spp_richness)
#ssp_Rich<-ssp_Rich[which(ssp_Rich$Subfamily==""),]
#dim(ssp_Rich)
#f_na<-fams[is.na(match(fams,ssp_Rich$Family))]
#f_na<-which(fams%in%f_na)
#ords <- ssp_Rich[match(fams,ssp_Rich$Family),"Order"]
#fams<-cbind(fams,as.character(ords))
#head(fams)
#fams[f_na,2]<-c("Dioscoreales","Santalales","Brassicales","Caryophyllales","Boraginales","Huerteales","Huerteales","Liliales","Zingiberales","Asterales")
#richness <- ssp_Rich[match(fams[,1],ssp_Rich$Family),"Spp_richness"]
#fams<-cbind(fams,richness)
#fams[is.na(match(fams[,1],ssp_Rich$Family)),1]
#ssp_Rich[is.na(match(ssp_Rich$Family,fams[,1])),"Spp_richness"]
#ssp_Rich[is.na(match(ssp_Rich$Family,fams[,1])),"Family"]
#fams[f_na,3]<-c("16","39","13","16","150","2","1","1","7","245")
#Data<-data.frame(Order=1:length(fams[,1]),z=as.numeric(fams[,3]))
#Data<-Data[order(Data$z),]
#Data$col<-rev(heat.colors(length(fams[,1])))
#orderedcolors<-Data[order(Data$Order),"col"]
#un_ords <- which(table(as.factor(fams[,2]))<=2)
#una <- which(fams[,2] %in% names(table(as.factor(fams[,2])))[un_ords])
#plot(pcoA[-una,],pch=21,bg=orderedcolors,cex=0.5)
#dataEllipse(x=pcoA[-una,],groups=as.factor(fams[-una,2]),levels=c(0.5),
#           center.pch=F, draw=TRUE, plot.points=F,
#            robust=F,col=orderedcolors[-una],lwd=1, fill=FALSE,group.labels = levels(as.factor(fams[-una,2])))

###### EXTRACT BIOS FOR ECO_REGIONS
setwd("~/Desktop/ANGIOSPERMS/DISTS./FINAL_DISTAS")
poly<-readOGR("~/Desktop/official/wwf_terr_ecos.shp")
ecos_area <- cbind(as.character(poly$ECO_NAME),poly$area_km2,as.character(poly$REALM),as.character(poly$BIOME))
ecos_area <- ecos_area[!duplicated(ecos_area),] ###### eco regions with area size
sort(unique(poly$ECO_NAME))
lista_eco <- paste(ecos_area[,1],ecos_area[,3],sep="_")
eco_PAM<- data.frame(Ecoregions=as.character(lista_eco))
cat("Extracting Bios","\n")
listBIO = list.files("~/Documents/MAPOTECA/",pattern=".asc")
for(j in 1:length(listBIO)) { 
  cat("Extracting BIOS",listBIO[j],"\n")
  r<-raster(paste(ruta,listBIO[j],sep="/"))
  x[,(j+3)]<-t(t(extract(r,pp)))
  nombreslistRaster[j]<-gsub(paste(".","asc",sep=""),"",listBIO[j],ignore.case=TRUE)
  x[,2]<-pp[,1]; x[,3]<-pp[,2];print(head(x))}

familias <- sort(unique(p_clean_2$family_ACCEPTED))
prPAM <- seq(1:dim(eco_PAM)[1])
tot <- length(familias)
length(unique(p_clean_2$species_ACCEPTED)); length(unique(p_clean_2$family_ACCEPTED))
p_clean_2 <- cbind(p_clean_2,rep("-",dim(p_clean_2)[1]))
p_clean_2 <- cbind(p_clean_2,rep("-",dim(p_clean_2)[1]))
p_clean_2 <- cbind(p_clean_2,rep("-",dim(p_clean_2)[1]))
names(p_clean_2)[14:16]<- c("eco_NAME","AREA","REALM")

library(raster)
setwd("~/Documents/MAPOTECA/TRABUCCO_etal/Priestley-Taylor Alpha Coefficient/")
lista_arch <- list.dirs()[-c(1,14)];lista_arch
r<-raster(paste(lista_arch[1],"w001001.adf",sep="/"))
r<-aggregate(r,fact=80,fun="mean")
r
for (i in 2:length(lista_arch)){
  cat(lista_arch[i],"\n")
  r1<- raster(paste(lista_arch[i],"w001001.adf",sep="/"))
  r1<- aggregate(r1,fact=40,fun="mean")
r<-stack(r,r1)
}
mean_r <- mean(r)
mean_r; plot(mean_r)
sd_r <- calc(r,fun=sd)
sd_r; plot(sd_r)

###### NINTH STEP: ANALYSES ON grid-PAM #########
setwd("~/Documents/NATURE/RESULTADOS_V4/DISTS./FINAL_DISTAS/PAMs_glonaf/")
load("prPAM_FAM_grid.r")
prPAM_g
max(rowSums(prPAM_g[,-c(1:2)]))
which(rowSums(prPAM_g[,-c(1:2)]) <= 100) -> iden_inf
prPAM_g<- prPAM_g[-iden_inf,]
prPAM_g[,-c(1:2)] <- decostand(prPAM_g[,-c(1:2)],method="hell",MARGIN=1) 
range_size <- colSums(prPAM_g[,-c(1:2)])
rich<-read.table("../SPP.RICHNESS.csv",heade=T,sep=",")
rich <- subset(rich,Subfamily=="");dim(rich);head(rich)
richa <- subset(rich,Family!="");dim(richa);head(richa)
richa$Family<-as.vector(richa$Family)
richa$Family[which(is.na(match(richa$Family,names(range_size))))]<-c("Ehretriaceae","Afrothismia","Rhipogonaceae")
which(is.na(match(richa$Family,names(range_size))))
richa[which(richa$Family=="Aphloiaceae"),5]<-2
richa$Range_size<- range_size[match(richa$Family,names(range_size))]
medias<-c()
varianzas<-c()
nombres<-c()
gridas<-list.files(pattern="_Grid_PAM.r")
for (i in 1:436){
load(gridas[i]);cat(toupper(sub("_Grid_PAM.r","",gridas[i])),"\n")
  nombres<-c(nombres,sub("_Grid_PAM.r","",gridas[i]))
if (dim(pp$Presence_and_Absence_Matrix)[2]<=3){ medias<-c(medias,sum(pp$Presence_and_Absence_Matrix[,3]))
varianzas<-c(varianzas,0)}
  if (dim(pp$Presence_and_Absence_Matrix)[2]>3){
medias<-c(medias,mean(colSums(pp$Presence_and_Absence_Matrix[,-c(1:2)])))
varianzas<-c(varianzas,var(colSums(pp$Presence_and_Absence_Matrix[,-c(1:2)])))
}}
richa$Mean_range_size<-medias[match(richa$Family,nombres)]
richa$Var_range_size<-varianzas[match(richa$Family,nombres)]
richness<-rowSums(prPAM_g[,-c(1:2)])

###### PLOT RANGE V. RICHNESS
par(mfrow=c(1,2))
pdf("~/Desktop/Figura_3.pdf")
plot(log(richa$Range_size),log(richa$Spp_richnes),pch=21,bg="grey50",cex=0.6,xlim=c(0,10),
     xlab="No. de celdas (log) por familia",ylab="Riqueza de especies (log)",las=1)
text(x = 0.2,y=10,labels="(A)",cex=1.5)
plot(log(richa$Mean_range_size),log(richa$Spp_richness),pch=21,bg="grey50",cex=0.6,xlim=c(0,6),
     xlab="Promedio no. de celdas (log) por familia",ylab="",las=1)
text(x = 0.2,y=10,labels="(B)",cex=1.5)
dev.off()

###### PLOT RICHNESS MAP
pdf("~/Desktop/Anexo_2.pdf")
maps::map("world",interior=F)
Data<-data.frame(Order=1:dim(prPAM_g)[1],z=richness)
Data<-Data[order(Data$z),]
range(log(richness))
ssp<-coordinates(prPAM_g[,1:2])
Data$col<-colorRampPalette(colors=c("olivedrab1","olivedrab3","khaki2","yellow2","red","darkred"))(dim(prPAM_g)[1])
orderedcolors<-Data[order(Data$Order),"col"]
points(ssp,pch=22,lwd=0.1,bg=orderedcolors,cex=log(richness)*0.1)
rect(xleft =-180,xright = 190,ybottom = -86,ytop = 85,border = "black")
dev.off()
######## PLOT GRIDS COLORED BY REALM
poly<-readOGR("~/Desktop/official/wwf_terr_ecos.shp")
ssp <- SpatialPoints(prPAM_g[,1:2],proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
pre <- over(ssp,poly,returnList = F)
spp_na <- ssp[which(is.na(pre$REALM))]
plot(spp_na)
#distas <- dist2Line(spp_na[1],poly)
#for (i in 46:47) {cat(i,"\n"); distas <- rbind(distas,dist2Line(spp_na[i],poly))}
#nanas<- which(is.na(pre$REALM))[-c(8,14,28,42,45)]
#prPAM_g3[nanas,1:2]<-distas[,2:3]
prPAM_g<-add_column(prPAM_g,pre$REALM,.before = T)
prPAM_gt<-prPAM_g[,-c(1:3)]
pcoA <- cmdscale(vegdist(prPAM_gt,"bray"),k = 2);plot(pcoA)
realms <- as.data.frame((prPAM_g[,1]))
realms<-as.vector(realms$`pre$REALM`)
unas<- unique(realms)
cols<-c("blue","white","lightblue2","gold","springgreen4","yellow3","mediumpurple1")
for (i in 1:length(unas)) realms[which(realms == unas[i])] <- cols[i]
plot(pcoA,bg=realms,pch=21)

###### TENTH STEP: MAP SPECIES RICHNESS PER FAMILY #######
setwd("~/Desktop/DISTS./FINAL_DISTAS/")
lista_arch<-list.files(pattern="_Grid_PAM.r")
par(mfrow=c(2,2))
pdf("RICHNESS_MAPS.pdf",onefile =T,width=14,height=7)
for (i in 1:length(lista_arch)){
load(lista_arch[i])
name<-lista_arch[i]
cat(lista_arch[i],"\n")
if(dim(pp$Presence_and_Absence_Matrix)[2]>3){ 
rowSums(pp$Presence_and_Absence_Matrix[,-c(1:2)])->SR
pama <- pp$Presence_and_Absence_Matrix[-which(SR==0),]
Data<-data.frame(Order=1:length(pama[,1]),z=rowSums(pama[,-c(1:2)]))
Data<-Data[order(Data$z),]
Data$col<-rev(heat.colors(length(pama[,1])))
orderedcolors<-Data[order(Data$Order),"col"]
maps::map("world",interior=F,col="grey80")
res <- scales::rescale(rowSums(pama[,-c(1:2)]),to = c(0.5,1),from=range(rowSums(pama[,-c(1:2)])))
points(pama[,1:2],cex=res*1.5,pch=21,col=adjustcolor("black",alpha.f = 0.5),lwd=0.5,bg=orderedcolors,xlim=c(-180,180),ylim=c(-90,90))
title(paste(sub("_Grid_PAM.r","",name)))
iden<-c(1+round(length(orderedcolors)/5),1+(round(length(orderedcolors)/5)*2),1+(round(length(orderedcolors)/5)*3),
        1+(round(length(orderedcolors)/5)*4),length(orderedcolors));s<- Data$z[iden];o<-Data$col[iden]
s<-Data$z[iden]
color.legend(xl=-50,xr=50,yb=-75,yt=-65,align="rb",gradient="x",legend=c(s[1],"","","",s[5]),rect.col=(adjustcolor(Data$col,alpha.f = 0.8)))
rect(xleft=-182,xright=190,ybottom =-85,ytop = 85,col=NULL,lwd=2)}
if(dim(pp$Presence_and_Absence_Matrix)[2]==3){ 
  pama <- pp$Presence_and_Absence_Matrix[-which(pp$Presence_and_Absence_Matrix[,3]==0),]
  maps::map("world",interior=F,col="grey80")
  res <- pama[,3]
  points(pama[,1:2],cex=res*1.5,pch=21,col=adjustcolor("black",alpha.f = 0.5),lwd=0.5,bg="red",xlim=c(-180,180),ylim=c(-90,90))
  title(paste(sub("_Grid_PAM.r","",name)))
  color.legend(xl=-50,xr=50,yb=-75,yt=-65,align="rb",gradient="x",legend=c(1),rect.col="red")
  rect(xleft=-182,xright=190,ybottom =-85,ytop = 85,col=NULL,lwd=2)
}
}
dev.off()
