#Generate the image dataset from raw processed file

#Read in image analysis table
imTab <- read.table("../rawdata/cellularityScoreKS1Str1.csv",sep=",",header=T,stringsAsFactors = F)
#Extract relevant information
imTab = imTab[,c(4,5,6,7,8,9,18)]
#Load in clinical to sample mapper
map <- read.table("../rawdata/Metabric1980IntClust.txt",header=T,row.names=2)
imTab$id = map[imTab[,7],1]

#Clean up table
imTab = imTab[!is.na(imTab$id),]

#Remove duplicates
imTab = imTab[!duplicated(imTab[,8]),]

#Remove labels
rownames(imTab) = imTab[,8]
imTab = imTab[,-c(7,8)]

#Calculate lymphocyte proportions
imTab$nlym1 = imTab$nLymphozyte/(imTab$nTumour + imTab$nLymphozyte)
imTab$nlym2 = imTab$nLymphozyte/(imTab$nLymphozyte + imTab$nArtifact + imTab$nTumour)


save(imTab,file="../rawdata/imageData.RData")




#T Cell Metric from METABRIC
load("../rawdata/cibersortScoresMETTCG2.RData")
metaCIB = data.frame(metaCIB)
metaCIB = metaCIB[metaCIB$P.value <= 0.05,]
tCount = rowSums(metaCIB[,c(grep(c("T.cells"),colnames(metaCIB)),grep(c("B.cells"),colnames(metaCIB)))])

#Intersect image file with T metric
int = intersect(names(tCount),rownames(imageData))
cor.test(tCount[int],imageData[int,1])
cor.test(tCount[int],imageData[int,2])

#Examine QTL
GAINS = t(GAINS)
int = intersect(colnames(GAINS),rownames(imageData))
GAINS = GAINS[,int]
imageData = imageData[int,]

outp = apply(GAINS,1,function(x){
  cor.test(x,imageData[,1])$p.value
})

#Examine TF relation
colnames(dout) = gsub("[.]","-",colnames(dout))
int = intersect(colnames(dout),rownames(imageData))
imageData = imageData[int,]
dout = dout[,int]
outq = apply(dout,1,function(x){
  cor.test(x,imageData[,2])$p.value
})

