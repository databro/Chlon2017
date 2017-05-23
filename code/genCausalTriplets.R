#Define package manager
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE, repos="http://cran.rstudio.com/")
    if(!require(x,character.only = TRUE))
    {
      source("http://bioconductor.org/biocLite.R")
      biocLite(x)
      if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }
}
packages = c("Hmisc","MatrixEQTL","DESeq","pheatmap","ppcor","corpcor","preprocessCore","survival","limma","reshape2","gplots","annotate","hgu133a.db","hgu133plus2.db","hgu95av2.db","hgug4110b.db","hgu133b.db", "GEOquery","clusterProfiler","mixtools","doMC","foreach","reshape2","ggplot2","VennDiagram","viper","RTN","gplots","qvalue","mixtools","metap")
sapply(packages,pkgTest)


#Load in METABRIC Transcriptome
load("../rawdata/METABRIC.RData")
#Load in Regulon for VIPER analysis
load("../rawdata/regulon.RData")
#Load in output of CIBERSORT analysis
load("../rawdata/cibersortScoresMETTCG2.RData")
#Load in METABRIC CNV data
load(file = "../rawdata/CNASeg_mean.RData")
#Load in TCGA raw counts data
load("../rawdata/brca-rawcounts.rda")
#Load in TCGA CNV data
rawcnv = read.table("../rawdata/all_data_by_genes.txt",sep="\t",row.names=1,stringsAsFactors = F,header=T)
#Set the format
loci = rawcnv[,c(1,2)]
rawcnv = rawcnv[,-c(1,2)]
colnames(rawcnv)= gsub("[.]","-",colnames(rawcnv))
GAINS = t(GAINS)

#Load Image Data
load("../rawdata/allData.rdata")
f1 = allData
load("../rawdata/allData_validation.rdata")
f2 = allData
imageData = rbind(f1,f2)
rownames(imageData) = imageData$Sample
imageData = imageData[,c("nLymphozyte","lym")]

#Convert the absolute number of lymphocytes in the tissue to a log based data column
imageData[,1] = log10(imageData[,1] + 1)

#To compare first we filter and keep only the primary tumour samples
#Keep only primary tumour samples
rawcnv= rawcnv[,grep("-01A-",colnames(rawcnv))]
#Amend the names
colnames(rawcnv) = substr(colnames(rawcnv),1,12)

#Convert ENTREZ ids to gene symbols using the annotate package for the gene expression matrix
GeneSymbol = getSYMBOL(rownames(rawcounts),data="org.Hs.eg")
rawcounts = data.frame(GeneSymbol = as.character(GeneSymbol),rawcounts)
#Collapse duplicate gene rows according to the row with the highest variance
var = apply(rawcounts[,-1],1,sd)^2
rawcounts = cbind(var=var,rawcounts)
rawcounts<- rawcounts[order(rawcounts[,2], rawcounts[,1], decreasing=TRUE),]
rawcounts <- rawcounts[!duplicated(rawcounts[,2]),]
#Check all duplicated rows have been removed
any(duplicated(rawcounts[,2]))
#Remove all NA genes
rawcounts = rawcounts[!is.na(rawcounts[,2]),]
#Assign them to rownames
rownames(rawcounts) = rawcounts[,2]
rawcounts = rawcounts[,-c(1,2)]

#Keep only primary tumour samples
colnames(rawcounts) = gsub("[.]","-",colnames(rawcounts))
rawcounts = rawcounts[,grep("-01",colnames(rawcounts))]
colnames(rawcounts) = gsub("-01","",colnames(rawcounts))


#Compare TCGA and METABRIC CNV Gene Expression distributions
plot(density(as.matrix(disc)))
lines(density(as.matrix(rawcounts)),col="red")

#Normalise TCGA rna-seq data using VST
vst <- function(countdata){
  condition <- factor(rep("Tumour", ncol(countdata)))
  countdata <- newCountDataSet(countdata,condition )
  countdata <- estimateSizeFactors( countdata )
  cdsBlind <- DESeq::estimateDispersions( countdata, method="blind")
  vstdata <- varianceStabilizingTransformation( cdsBlind )
  return(exprs(vstdata))
}

rawcounts.norm = normalize.quantiles.use.target(vst(as.matrix(rawcounts)),apply(disc,1,median))
plot(density(as.matrix(disc)[1,]))
lines(density(as.matrix(rawcounts.norm)[1,]),col="red")

rownames(rawcounts.norm) = rownames(rawcounts)
colnames(rawcounts.norm) = colnames(rawcounts)
rawcounts = rawcounts.norm

#Load in and manage PPI information
PPI = read.table("../rawdata/STRING_ixn_human_v9.05_ALL.txt",header=T,stringsAsFactors = F)
PPI[,1] = gsub("9606.","",PPI[,1])
PPI[,2] = gsub("9606.","",PPI[,2])
loadPrev = TRUE
if(!loadPrev){
  library(biomaRt)
  # define biomart object
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  
  table = lapply(unique(PPI[,1]),function(x){
    print(x)
    results <- unlist(getBM(attributes = c("hgnc_symbol"),
                            filters = "ensembl_peptide_id", values = x,
                            mart = mart))})
  names(table) = unique(PPI[,1])
  table = unlist(table)
  names(table) = gsub(".hgnc_symbol","",names(table))
}else{
  load(file="../rawdata/Hgnc2Ensembl.RData")}
invtable = names(table)
names(invtable) = table
backupPPI = PPI
metabricVIPER = viper(disc, regul,minsize = 5)
tcgaVIPER = viper(rawcounts, regul,minsize = 5)

geom_mean = function(x){
  tmp = apply(x,1,prod)
  tmp^(1/ncol(x))
}

#Run VIPER and get regulon activities
discoveryTriplets = getTripletMat(metabricVIPER,data.frame(cytoScore = geom_mean(t(disc[c("PRF1","GZMA"),])),TC = geom_mean(t(disc[c("CD3D","CD4","CD8A","CD8B"),]))),GAINS,backupPPI,0.001,TRUE)
#Get TCGA triplets
tcgaTriplets = getTripletMat(tcgaVIPER,data.frame(cytoScore = geom_mean(t(rawcounts[c("PRF1","GZMA"),])),TC = geom_mean(t(rawcounts[c("CD3D","CD4","CD8A","CD8B"),]))),rawcnv,backupPPI,0.05,FALSE)
#Get Image Triplets
imageTriplets = getTripletMat(metabricVIPER,imageData,GAINS,backupPPI,0.05,FALSE)

save(discoveryTriplets, tcgaTriplets, imageTriplets, file="../Computed Objects/dvitriplets2.RData")
