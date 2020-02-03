# Estimates and plots family crown age ranges obtained from the Angiosperm Phylogeny Website. 
# The code also plots the mean family crown ages obtained from the Angiosperm Time Tree 2.0 (ATT 2.0).
# Must run scripts 1, 2 and 3.
# Uses additional files.
APWeb <- "Parsed_APWeb.csv" # file in data folder

####### The following chunk of code generates the "Parsed_APWeb.csv" file.
rawAPWeb <- "CGages_APWeb.csv" # file in data folder
ages<-read.table(rawAPWeb,sep=",",header=T)
colnames(ages)
table(ages$FAMILY)
un_node <- as.vector(unique(ages$NODE))
un_node <- un_node[match(rownames(resB_fams),un_node)]
un_clades <- ages[match(un_node,ages$NODE),"FAMILY"]
edades <- c();maximos<-c();minimos <- c()
for (i in 1: length(un_node)){ 
  cat(i,"\n")
  if(is.na(un_node[i])){ edades <- c(edades, NA);maximos <- c(maximos,NA);minimos <- c(minimos,NA);next}
  if(dim(ages[ages$NODE==un_node[i],])[1]==1){ 
    edades <- c(edades,mean(as.numeric(ages[ages$NODE==un_node[i],9:11]),na.rm = T))
    maximos <- c(maximos,max(as.numeric(ages[ages$NODE==un_node[i],9:11]),na.rm = T))
    minimos <- c(minimos,min(as.numeric(ages[ages$NODE==un_node[i],9:11]),na.rm = T))
    next}
  edades <- c(edades, mean(ages[ages$NODE==un_node[i],"mean"],na.rm=T))
  maximos <- c(maximos, max(ages[ages$NODE==un_node[i],"mean"],na.rm=T))
  minimos <- c(minimos, min(ages[ages$NODE==un_node[i],"mean"],na.rm=T))
}
data <- data.frame(edades=edades,max=maximos,min=minimos,family_APWeb=un_node,clade=un_clades,family_beast=rownames(resB_fams),beast_age=resB_fams$CG_age_Acal)
data$grupos <- rep(NA,dim(data)[1])
head(data)
data[which(data$clade==1),"grupos"] <- "ANA"
data[which(data$clade==2),"grupos"] <- "Magnoliids"
data[which(data$clade==3),"grupos"] <- "Chlor + Cera"
data[which(data$clade==4),"grupos"] <- "Monocots"
data[which(data$clade==5),"grupos"] <- "Monocots"
data[which(data$clade==6),"grupos"] <- "Eudicots"
data[which(data$clade==7),"grupos"] <- "Eudicots"
data[which(data$clade==8),"grupos"] <- "Eudicots"
data[which(data$clade==9),"grupos"] <- "Eudicots"
data[which(data$clade==10),"grupos"] <- "Eudicots"
data[which(data$clade==11),"grupos"] <- "Eudicots"
data[which(data$clade==12),"grupos"] <- "Eudicots"
data[which(is.na(data$clade)),"grupos"] <- "no data"
write.table(data,file="Parsed_APWeb.csv",sep=",",col.names=T, row.names=T)
##########

data<- read.table(APWeb,sep=",",header=T)
data <- data[order(data$grupos,data$clade),]
unique(data$grupos)->un_clados
tre.pruned <- drop.tip(tre,tip=list[which(duplicated(list_fams))])
tre.pruned<- ladderize(tre.pruned, FALSE)
list <- tre.pruned$tip.label
list_outs <- foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][2]}
outs<-c("Welwitschiaceae","Pinaceae","Ginkgoaceae","Cupressaceae","Gnetaceae","Cycadaceae")
tre.pruned <- drop.tip(tre.pruned,tip=tre.pruned$tip.label[which(list_outs%in%outs)])
list<-tre.pruned$tip.label
list <- foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][2]}
bb <- resB_fams[match(list,rownames(resB_fams)),]
match(list,rownames(bb))
print(bb[1:5,];tre.pruned$tip.label[1:5])
list<-tre.pruned$tip.label
list_fams <- foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][2]}
tre.pruned$tip.label<-list_fams
is.rooted(tre.pruned)
is.ultrametric(tre.pruned)
data[match(tre.pruned$tip.label,data$family_beast),"beast_age"]
data[which(data$family_APWeb=="Chloranthaceae"),c(2:3)]<-NA
data<-data[match(tre.pruned$tip.label,data$family_beast),]
HPD_max <- resB_fams[match(data$family_beast,rownames(resB_fams)),"CG_maxHPD"]
HPD_min <- resB_fams[match(data$family_beast,rownames(resB_fams)),"CG_minHPD"]
data$HPD_max <- HPD_max; data$HPD_min <- HPD_min

pdf(paste(ruta_write,"DatedTree_APWED.pdf"),useDingbats = F,width=21,height = 100)
d1 <- data.frame(id=tre.pruned$tip.label, beast=data[match(tre.pruned$tip.label,data$family_beast),"beast_age"],
                 HPDmax=data[match(tre.pruned$tip.label,data$family_beast),"HPD_max"],
                 HPDmin=data[match(tre.pruned$tip.label,data$family_beast),"HPD_min"])
p <- ggtree(tre.pruned) + geom_tiplab() + ggplot2::xlim(0,250)
p2 <- facet_plot(p, panel="Family crown age estimates", data=d1, geom=geom_point, aes(x=beast),pch=21, bg='red',cex=8)
d2 <- data.frame(id=tre.pruned$tip.label, max=data[match(tre.pruned$tip.label,data$family_APWeb),"max"],min=data[match(tre.pruned$tip.label,data$family_APWeb),"min"])
facet_plot(p2, panel='Family crown age estimates', data=d2, geom=geom_segment, aes(x=min, xend=max, y=y, yend=y), size=5, color = rgb(0,0,0,0.6)) + theme_tree2()
dev.off()
