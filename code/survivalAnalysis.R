#Load in survival data for METABRIC
load("../rawdata/metabricClinical.RData")
#Load in gene expression data for METABRIC
load("../rawdata/METABRIC.RData")
#Load in regulon
load("../rawdata/regulon.RData")
#Load in CNV data
load(file = "../rawdata/CNASeg_mean.RData")

#Housekeeping
colnames(disc) = gsub("[.]","-",colnames(disc))
clinical = clinical[rownames(clinical) %in% colnames(disc),]
disc = disc[,colnames(disc) %in% rownames(clinical)]
disc = disc[,rownames(clinical)]
GAINS = t(GAINS)
GAINS = GAINS[,colnames(disc)]

#Run VIPER
dout = viper(disc, regul,minsize = 5)

#TFs of interest
tfkeep = as.character(unique(discoveryTriplets$TF))
dout = dout[c(tfkeep,"ESR1"),]

#CNVs of interest
cnvKeep = as.character(unique(discoveryTriplets$CNV))
GAINS = GAINS[cnvKeep,]

cnvER = apply(GAINS,1,function(x){
  t.test(x ~ clinical$ER.Expr)$statistic
})

tfER = apply(dout,1,function(x){
  t.test(x ~ clinical$ER.Expr)$statistic
})


to.plot = data.frame(cnvER[as.character(discoveryTriplets$CNV)],
                     tfER[as.character(discoveryTriplets$TF)],
                     paste(discoveryTriplets$CNV,discoveryTriplets$TF,sep="-"))
colnames(to.plot) = c("CNV","TF","Triplet")
to.plot = melt(to.plot,id="Triplet")
to.plot = data.frame(to.plot,erstat = sign(to.plot$value),stringsAsFactors = F)
to.plot$erstat[p.adjust(pt(abs(to.plot$value),df=ncol(dout),lower.tail = F),method="BH") >= 0.05] = 0

to.plot$erstat[to.plot$erstat == 1] = "ER-"
to.plot$erstat[to.plot$erstat == "0"] = "Inconclusive"
to.plot$erstat[to.plot$erstat == "-1"] = "ER+"
to.plot$erstat = as.factor(to.plot$erstat)
colnames(to.plot) = c("Triplet","Variable","T_statistic","Significant_Association")
to.plot= data.frame(to.plot,stringsAsFactors=T)

require(ggplot2)
ggplot(to.plot, aes(y = Triplet,
                  x = Variable)) +        ## global aes      ## to get the rect filled
  geom_point(aes(colour = T_statistic,shape=Significant_Association),size = 5)  +    ## geom_point for circle illusion
  scale_color_gradient(low = "blue",  
                       high = "red")+       ## color of the corresponding aes
  scale_size(range = c(1, 20))+             ## to tune the size of circles
  theme_bw()+
  ggtitle("Triplet Variable Associations between \n ER+/ER- Samples")
