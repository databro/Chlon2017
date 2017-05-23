#Load in image triplets
load(file="../Computed Objects/dvitriplets2.RData")
#Load in lymphocyte markers
load("../rawdata/phenData.RData")
#Load Image Data
load("../rawdata/allData.rdata")
f1 = allData
load("../rawdata/allData_validation.rdata")
f2 = allData
imageData = rbind(f1,f2)
rownames(imageData) = imageData$Sample
imageData = imageData[,c("nLymphozyte","lym")]
plot(density(imageData$nLymphozyte))
plot(density(imageData$lym))

               
#Intersect common general causal triplets between TCGA and METABRIC
dInd = paste(discoveryTriplets[,1],discoveryTriplets[,2],discoveryTriplets$model,sep=":")
vInd = paste(tcgaTriplets[,1],tcgaTriplets[,2],tcgaTriplets$model,sep=":")
int = intersect(dInd,vInd)
iInd = paste(imageTriplets[,1],imageTriplets[,2],imageTriplets$model,sep=":")

tmp = data.frame(Model = unlist(lapply(strsplit(int,":"),function(x){x[[3]]})), 
                 Validated = int %in% iInd,
                 stringsAsFactors = F)
to.plot = data.frame(Model = c("M1","M2"),rbind(table(tmp[tmp$Model == 1,2]), table(tmp[tmp$Model == 2,2])),stringsAsFactors = F)
to.plot = melt(to.plot,id="Model")
colnames(to.plot) = c("Model","Validated_In_Images","value")
ggplot(to.plot, aes(Model, value, fill=Validated_In_Images)) + geom_bar(stat="identity") +
  ggtitle("Validation of METABRIC/TCGA in \n Images")+
  theme(text = element_text(size=18),
        plot.title = element_text(lineheight=.8, face="bold"),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(angle = 90,size=13)) +
  xlab("Model") + ylab("Number of Triplets")


CT = "TC"
immMarker = "lym"
discoveryTriplets = discoveryTriplets[discoveryTriplets$Type == CT,]
tcgaTriplets = tcgaTriplets[tcgaTriplets$Type == CT,]
imageTriplets = imageTriplets[imageTriplets$Type == immMarker,]

#Intersect common general causal triplets between TCGA and METABRIC
dInd = paste(discoveryTriplets[,1],discoveryTriplets[,2],discoveryTriplets$model,sep=":")
vInd = paste(tcgaTriplets[,1],tcgaTriplets[,2],tcgaTriplets$model,sep=":")
int = intersect(dInd,vInd)

#Further intersect this list with the list of causal triplets computed from the image data
iInd = paste(imageTriplets[,1],imageTriplets[,2],imageTriplets$model,sep=":")
vInt = intersect(int,iInd)

#How many image triplets validate in the METABRIC-TCGA validated triplet set?
paste0(round(length(vInt)/length(int) * 100,3), "% of triplets successfully validated in images")


#Isolate CIBERSORT T cell causal triplets for Image Analysis validation
rownames(discoveryTriplets) = dInd
rownames(imageTriplets) = iInd
discoveryTriplets = discoveryTriplets[vInt,]
imageTriplets = imageTriplets[vInt,]

#How well does the variance explained agree with the genomics discovery data?
cTest = cor.test(imageTriplets$VIPER.Phen.Cor,discoveryTriplets$VIPER.Phen.Cor)
#How do the VIPER-IMMUNE correlations and CNV-VIPER correlations look?
to.plot= data.frame(Im = imageTriplets$VIPER.Phen.Cor,Meta = discoveryTriplets$VIPER.Phen.Cor,Model = discoveryTriplets[,"model"])
to.plot = to.plot[to.plot[,3] %in% c("1","2"),]
ggplot(to.plot,aes(x=Im,y=Meta,colour=Model))+
  geom_point(alpha=0.4,size=1.75)+
  xlab("VIPER-Immune Correlation Coefficient \n Microarray-METABRIC")+
  ylab("Lymphocyte Count")+
  ggtitle("VIPER-Immune Correlation Coefficient\n Microarray vs Lymphocyte Count")

#Control for triplets that have the same associative directionality as the discovery cohort
keep = sign(discoveryTriplets$VIPER.Phen.Cor) * sign(imageTriplets$VIPER.Phen.Cor) > 0
paste0(length(which(keep))," out of ", length(keep), " triplets have same associative directionality")
imageTriplets = imageTriplets[keep,]
discoveryTriplets = discoveryTriplets[keep,]

#To what extent do the TFs regulate the influx of immune cells into the tumour?
to.plot = data.frame(P = -log10(discoveryTriplets$VIPER.Phen.Pval),Cor = discoveryTriplets$VIPER.Phen.Cor)
qplot(to.plot[,1],to.plot[,2],label=discoveryTriplets$TF,geom="text",xlab="-log10(p-value)",ylab="Pearson Correlation Coefficient", 
      main="Image T-cell Infiltrate vs TF activity\nModel 1")

#What communities can we identify in the regulons of the positive regulators?
#Load in Overall Regulon for VIPER analysis
load("../rawdata/regulon.RData")
#Get regulons of +ve TFs
discoveryTriplets = discoveryTriplets[which(sign(imageTriplets$VIPER.Phen.Cor) == sign(discoveryTriplets$VIPER.Phen.Cor)),]
pRegTF = as.character(discoveryTriplets$TF[discoveryTriplets$VIPER.Phen.Cor > 0])
pReg = lapply(pRegTF,function(x){
data.frame(Targets = names(regul[[x]]$tfmode),TF = x,
           stringsAsFactors = F)})
pReg = do.call(rbind.data.frame,pReg)

#Construct the graph
g <- graph.data.frame(pReg,directed=F,vertices=NULL)
ig <- igraph::simplify(g,remove.multiple = TRUE)
#How many shared genes are there?
print(length(which(duplicated(pReg[,1]))))
dupGenes = pReg[duplicated(pReg[,1]),1]
#Plot the communities using RedeR
library("RedeR")
library("igraph")
rdp <- RedPort()
calld(rdp, maxlag = 5000)
resetd(rdp)
vnames <- get.vertex.attribute(g, "name")
V(g)$nodeFontSize <- ifelse(vnames %in% unique(pReg$TF), 100,1)
V(g)$nodeSize <- ifelse(vnames %in% unique(pReg$TF), 100,30)
col = rep(NA,length(vnames))
names(col) = vnames

#Examine the TF gene interaction network and highlight duplicated genes
col[c(pReg[pReg[,2] == unique(pReg[,2])[1],1],unique(pReg[,2])[1])] = "deeppink"
col[c(pReg[pReg[,2] == unique(pReg[,2])[2],1],unique(pReg[,2])[2])] = "royalblue"
col[c(pReg[pReg[,2] == unique(pReg[,2])[3],1],unique(pReg[,2])[3])] = "springgreen"
col[dupGenes] = "gray"
V(g)$nodeColor <- col
addGraph(rdp,g,layout = NULL)

#What communities can we identify in the regulons of the negative regulators?
#Get regulons of -ve TFs
nRegTF = as.character(discoveryTriplets$TF[discoveryTriplets$VIPER.Phen.Cor < 0])
nReg = lapply(unique(nRegTF),function(x){
  data.frame(Targets = names(regul[[x]]$tfmode),TF = x,
             stringsAsFactors = F)})
nReg = do.call(rbind.data.frame,nReg)

#Construct the graph
g <- graph.data.frame(nReg,directed=F,vertices=NULL)
ig <- igraph::simplify(g,remove.multiple = TRUE)
#How many shared genes are there?
print(length(which(duplicated(nReg[,1]))))
dupGenes = nReg[duplicated(nReg[,1]),1]

#Plot the communities using RedeR
rdp <- RedPort()
calld(rdp, maxlag = 5000)
resetd(rdp)
vnames <- get.vertex.attribute(g, "name")
V(g)$nodeFontSize <- ifelse(vnames %in% unique(nReg$TF), 100,1)
V(g)$nodeSize <- ifelse(vnames %in% unique(nReg$TF), 100,30)
col = rep(NA,length(vnames))
names(col) = vnames

#Examine the TF gene interaction network and highlight duplicated genes
col[c(nReg[nReg[,2] == unique(nReg[,2])[1],1],unique(nReg[,2])[1])] = "deeppink"
col[c(nReg[nReg[,2] == unique(nReg[,2])[2],1],unique(nReg[,2])[2])] = "royalblue"
col[c(nReg[nReg[,2] == unique(nReg[,2])[3],1],unique(nReg[,2])[3])] = "springgreen"
col[c(nReg[nReg[,2] == unique(nReg[,2])[4],1],unique(nReg[,2])[4])] = "orange"
col[dupGenes] = "gray"
V(g)$nodeColor <- col
addGraph(rdp,g,layout = NULL)



#GO term enrichment from the merged regulons
#Now we examine the functional annotation of all these genes to their respective pathways for the ER positive and ER negative z-scores
xx <- as.list(org.Hs.egALIAS2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx <- xx[!is.na(xx)]
pos_regulators = xx[pReg[,1]]
neg_regulators = xx[nReg[,1]]

geneids = list(unlist(pos_regulators), unlist(neg_regulators))
names(geneids) = c("Positive TF Regulons","Negative TF Regulons")
cep <- compareCluster(geneids, fun="enrichGO", ont="BP", OrgDb='org.Hs.eg.db')
plot(cep)

