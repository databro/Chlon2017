setwd("~/Dropbox/causality/leon/code/")
#In this section we examine the distribution of causal triplets across both the METABRIC and TCGA cohorts
#to infer which models are robust with respect to explaining the variance seen in infiltrating immune cells

#Load in triplet data
load(file="../Computed Objects/dvitriplets2.RData")

#What are the poportion of models in the discovery triplets?
table(discoveryTriplets$model)/sum(table(discoveryTriplets$model))

#What does the distribution of significant models look like between TCGA and METABRIC?
to.plot = data.frame(cohort = c("METABRIC","TCGA"),rbind(table(discoveryTriplets$model),table(tcgaTriplets$model)),stringsAsFactors = F)
colnames(to.plot)[c(2,3)] = c("M1", "M2")
colnames(to.plot) = gsub("[.]"," ",colnames(to.plot))
library("reshape2")
library("ggplot2")

to.plot = melt(to.plot,id="cohort")
ggplot(to.plot, aes(variable, value, fill=cohort)) + geom_bar(stat="identity",position = "dodge") +
  ggtitle("Model Distribution\nMETABRIC & TCGA")+
  theme(text = element_text(size=18),
        plot.title = element_text(lineheight=.8, face="bold"),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(angle = 90,size=13)) +
  xlab("Model") + ylab("Number of Triplets")

#Intersect common general causal triplets
dInd = paste(discoveryTriplets[,1],discoveryTriplets[,2],discoveryTriplets$model,sep=":")
vInd = paste(tcgaTriplets[,1],tcgaTriplets[,2],tcgaTriplets$model,sep=":")
int = intersect(dInd,vInd)

#How many common causal triplets intersect?
paste0(round(length(int)/length(unique(dInd)) * 100,3), "% of general triplets successfully validated")

discoveryTriplets$ind = dInd
tcgaTriplets$ind = vInd

tmp = data.frame(Model = discoveryTriplets$model, 
                     Validated = discoveryTriplets$ind %in% int,
                     stringsAsFactors = F)
to.plot = data.frame(Model = c("M1","M2"),rbind(table(tmp[tmp$Model == 1,2]), table(tmp[tmp$Model == 2,2])),stringsAsFactors = F)
to.plot = melt(to.plot,id="Model")
colnames(to.plot) = c("Model","Validation_Status","value")
ggplot(to.plot, aes(Model, value, fill=Validation_Status)) + geom_bar(stat="identity") +
  ggtitle("Validation of METABRIC Triplets in \n TCGA")+
  theme(text = element_text(size=18),
        plot.title = element_text(lineheight=.8, face="bold"),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(angle = 90,size=13)) +
  xlab("Model") + ylab("Number of Triplets")

#How many models validate on an individual level?
#Model M1
table(tmp[tmp[,1] == 1,2])
table(tmp[tmp[,1] == 1,2])/sum(table(tmp[tmp[,1] == 1,2]))

table(tmp[tmp[,1] == 2,2])
table(tmp[tmp[,1] == 2,2])/sum(table(tmp[tmp[,] == 2,2]))

#Consider validated examples on a phenotype-specific basis
discoveryTriplets$key = discoveryTriplets$ind %in% int
tmp = data.frame(discoveryTriplets$model,discoveryTriplets$key,discoveryTriplets$Type)
colnames(tmp) = c("Model","Validated","Phenotype")
to.plot = data.frame(Model = c("M1","M2","M1","M2"),Phenotype = c("TC","TC","cytoScore","cytoScore"),rbind(rbind(table(tmp[tmp$Model == 1 & tmp$Phenotype == "TC",2]), cytoScore = table(tmp[tmp$Model == 2 & tmp$Phenotype == "TC",2])),rbind(table(tmp[tmp$Model == 1 & tmp$Phenotype == "cytoScore",2]), table(tmp[tmp$Model == 2 & tmp$Phenotype == "cytoScore",2]))),stringsAsFactors = F)
F1 = melt(to.plot[to.plot$Model == "M1",-1],id.vars = "Phenotype")
colnames(F1) = c("Phenotype","Validation_Status","value")
ggplot(F1, aes(Phenotype, value, fill=Validation_Status)) + geom_bar(stat="identity") +
  ggtitle("Validation of M1 METABRIC Triplets in \n TCGA")+
  theme(text = element_text(size=18),
        plot.title = element_text(lineheight=.8, face="bold"),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(angle = 90,size=13)) +
  xlab("Model") + ylab("Number of Triplets")

F2 = melt(to.plot[to.plot$Model == "M2",-1],id.vars = "Phenotype")
colnames(F2) = c("Phenotype","Validation_Status","value")
ggplot(F2, aes(Phenotype, value, fill=Validation_Status)) + geom_bar(stat="identity") +
  ggtitle("Validation of M2 METABRIC Triplets in \n TCGA")+
  theme(text = element_text(size=18),
        plot.title = element_text(lineheight=.8, face="bold"),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(angle = 90,size=13)) +
  xlab("Model") + ylab("Number of Triplets")


#How many triplets are common between the TC phenotype and the cytoScore phenotype?
discoveryTriplets = discoveryTriplets[discoveryTriplets$ind %in% int,]
phenoProp = c(nrow(discoveryTriplets[discoveryTriplets$Type == "TC",]),length(which(duplicated(discoveryTriplets$ind))),nrow(discoveryTriplets[discoveryTriplets$Type == "cytoScore",]))
phenoProp[c(1,3)] = phenoProp[c(1,3)] - phenoProp[2]
to.plot = data.frame(x = c("All Models","All Models","All Models"),c("TC","Intersected","cytoScore"),phenoProp)
colnames(to.plot) = c("X","Phenotype","Frequency")
ggplot(to.plot, aes(x = X, y = Frequency, fill=Phenotype)) + geom_bar(stat="identity") +
  ggtitle("Frequency of Validated Triplets by Phenotype\n METABRIC + TCGA")+
  theme(text = element_text(size=18),
        plot.title = element_text(lineheight=.8, face="bold"),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(angle = 90,size=13)) +
  xlab("Model") + ylab("Number of Triplets")


#Intersect common specific causal triplets
dInd = paste(discoveryTriplets[,1],discoveryTriplets[,2],discoveryTriplets$model,discoveryTriplets$Type,sep=":")
vInd = paste(tcgaTriplets[,1],tcgaTriplets[,2],tcgaTriplets$model,tcgaTriplets$Type,sep=":")
int = intersect(dInd,vInd)

#What proportion of specific lymphocyte trplet METABRIC triplets validates in TCGA?
paste0(round(length(int)/nrow(discoveryTriplets) * 100,3), "% of triplets successfully validated")

#Assign new labels to their own index column
rownames(discoveryTriplets)  = dInd
rownames(tcgaTriplets) = vInd

gint = int
int = gint[grep(":TC",gint)]
#How well do the correlations between VIPER and the immune infiltrate agree between both cohorts?
test = cor.test(tcgaTriplets[int,"VIPER.Phen.Cor"],discoveryTriplets[int,"VIPER.Phen.Cor"])
plot(discoveryTriplets[int,"VIPER.Phen.Cor"],tcgaTriplets[int,"VIPER.Phen.Cor"],main = "VIPER- TCS Correlation Coefficients\n METABRIC vs TCGA",xlab="METABRIC",ylab="TCGA")
abline(lm(tcgaTriplets[int,"VIPER.Phen.Cor"] ~ discoveryTriplets[int,"VIPER.Phen.Cor"]), col='red')
legend("topleft", bty="n", 
       legend=paste("r=",format(test$estimate,digits=4),"\nP=",test$p.value),
       cex=1)
int = gint[grep(":cytoScore",gint)]
#How well do the correlations between VIPER and the immune infiltrate agree between both cohorts?
test = cor.test(tcgaTriplets[int,"VIPER.Phen.Cor"],discoveryTriplets[int,"VIPER.Phen.Cor"])
plot(discoveryTriplets[int,"VIPER.Phen.Cor"],tcgaTriplets[int,"VIPER.Phen.Cor"],main = "VIPER-cytoScore Correlation Coefficients\n METABRIC vs TCGA",xlab="METABRIC",ylab="TCGA")
abline(lm(tcgaTriplets[int,"VIPER.Phen.Cor"] ~ discoveryTriplets[int,"VIPER.Phen.Cor"]), col='red')
legend("topleft", bty="n", 
       legend=paste("r=",format(test$estimate,digits=4),"\nP=",test$p.value),
       cex=1)


#Isolate duplicate triplets and inspect phenotype-specific triplets
toIsolate = discoveryTriplets$ind[duplicated(discoveryTriplets$ind)]
gdiscoveryTriplets = discoveryTriplets
gtcgaTriplets = tcgaTriplets

gint = intersect(rownames(discoveryTriplets),rownames(tcgaTriplets))
int = gint[grep(":TC",gint)]
#Heatmap the correlation between the TF and the Immune Cell Across METABRIC and TCGA
to.heatmap = data.frame(METABRIC= discoveryTriplets[int,"VIPER.Phen.Cor"],TCGA = tcgaTriplets[int,"VIPER.Phen.Cor"])
rownames(to.heatmap) = int
to.heatmap = data.frame(discoveryTriplets[int,1],to.heatmap,stringsAsFactors = F)
to.heatmap = to.heatmap[!duplicated(to.heatmap[,1]),-1]
to.heatmap = to.heatmap[order(to.heatmap[,1],decreasing = TRUE),]
to.heatmap = rbind(to.heatmap[1:10,],to.heatmap[-c(1:(nrow(to.heatmap) - 10)),])
pheatmap(to.heatmap,main="VIPER-TCS Edge Correlation Coefficient\nMETABRIC + TCGA\n")

int = gint[grep(":cytoScore",gint)]
#Heatmap the correlation between the TF and the Immune Cell Across METABRIC and TCGA
to.heatmap = data.frame(METABRIC= discoveryTriplets[int,"VIPER.Phen.Cor"],TCGA = tcgaTriplets[int,"VIPER.Phen.Cor"])
rownames(to.heatmap) = int
to.heatmap = data.frame(discoveryTriplets[int,1],to.heatmap,stringsAsFactors = F)
to.heatmap = to.heatmap[!duplicated(to.heatmap[,1]),-1]
to.heatmap = to.heatmap[order(to.heatmap[,1],decreasing = TRUE),]
to.heatmap = rbind(to.heatmap[1:10,],to.heatmap[-c(1:(nrow(to.heatmap) - 10)),])
pheatmap(to.heatmap,main="VIPER-cytoScore Edge Correlation Coefficient\nMETABRIC + TCGA\n")

#Heatmap the correlation between the TF and the TCS
int = gint[grep(":cytoScore",gint)]
to.heatmap2 = data.frame(CNV.Phen= discoveryTriplets[int,"CNV.Phen.Cor"],TF.Phen= discoveryTriplets[int,"VIPER.Phen.Cor"])
rownames(to.heatmap2) = int
to.heatmap = to.heatmap2[rownames(to.heatmap),]
to.heatmap = rbind(head(to.heatmap),tail(to.heatmap))
pheatmap(to.heatmap,main="CNV/TF Association with CytoScore Trait\n METABRIC")

#Cytoscore analysis 
cytoAnalysis = discoveryTriplets[discoveryTriplets$Type == "cytoScore",]
#Load in Overall Regulon for VIPER analysis
load("../rawdata/regulon.RData")
#Get regulons of +ve TFs
cytoAnalysis = cytoAnalysis[cytoAnalysis$CNV %in% c("NCOR1","EP300"),]

pRegTF = as.character(cytoAnalysis$TF[cytoAnalysis$VIPER.Phen.Cor > 0])
pReg = lapply(pRegTF,function(x){
  data.frame(Targets = names(regul[[x]]$tfmode),TF = x,
             stringsAsFactors = F)})
pReg = do.call(rbind.data.frame,pReg)

nRegTF = as.character(cytoAnalysis$TF[cytoAnalysis$VIPER.Phen.Cor < 0])
nReg = lapply(nRegTF,function(x){
  data.frame(Targets = names(regul[[x]]$tfmode),TF = x,
             stringsAsFactors = F)})
nReg = do.call(rbind.data.frame,nReg)

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



