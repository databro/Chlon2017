setwd("~/Dropbox/causality/leon/code/")
#In this section we examine the distribution of causal triplets across both the METABRIC and TCGA cohorts
#to infer which models are robust with respect to explaining the variance seen in infiltrating immune cells

#Load in triplet data
load(file="../Computed Objects/dvitriplets.RData")

#What does the distribution of significant models look like between TCGA and METABRIC?
to.plot = data.frame(cohort = c("METABRIC","TCGA"),rbind(table(discoveryTriplets$model),table(validationTriplets$model)),stringsAsFactors = F)
colnames(to.plot) = gsub("[.]"," ",colnames(to.plot))
library("reshape2")
library("ggplot2")

to.plot = melt(to.plot,id="cohort")
ggplot(to.plot, aes(variable, value, fill=cohort)) + geom_bar(stat="identity",position = "dodge") +
  ggtitle("model Distribution\nMETABRIC & TCGA")+
  theme(text = element_text(size=18),
        plot.title = element_text(lineheight=.8, face="bold"),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(angle = 90,size=13))

#Intersect common causal triplets
rownames(discoveryTriplets) = paste(discoveryTriplets$TF,discoveryTriplets$`CNV`,discoveryTriplets$model,discoveryTriplets$Cell,sep=":")
rownames(validationTriplets) = paste(validationTriplets$TF,validationTriplets$`CNV`,validationTriplets$model,validationTriplets$Cell,sep=":")
int = intersect(rownames(discoveryTriplets),rownames(validationTriplets))

#What proportion of METABRIC validates in TCGA?
paste0(round(length(int)/nrow(discoveryTriplets) * 100,3), "% of triplets successfully validated")

#What does the model distribution look like now?
barplot(table(discoveryTriplets[int,"model"]),main = "Number of METABRIC models validated in TCGA",xlab = "model",ylab="Total")

#How well do the correlations between VIPER and the immune infiltrate agree between both cohorts?
test = cor.test(validationTriplets[int,"VIPER.IMMUNE.Cor"],discoveryTriplets[int,"VIPER.IMMUNE.Cor"])
plot(discoveryTriplets[int,"VIPER.IMMUNE.Cor"],validationTriplets[int,"VIPER.IMMUNE.Cor"],main = "VIPER-CIBERSORT Correlation Coefficients\n METABRIC vs TCGA",xlab="METABRIC",ylab="TCGA")
abline(lm(validationTriplets[int,"VIPER.IMMUNE.Cor"] ~ discoveryTriplets[int,"VIPER.IMMUNE.Cor"]), col='red')
legend("topleft", bty="n", 
       legend=paste("r=",format(test$estimate,digits=4),"\nP=",test$p.value),
       cex=1)

#How do the VIPER-IMMUNE correlations and CNV-VIPER correlations look?
to.plot= data.frame(VI = discoveryTriplets[int,"VIPER.IMMUNE.Cor"],CV = discoveryTriplets[int,"CNV.VIPER.Cor"],model = discoveryTriplets[int,"model"])
to.plot = to.plot[to.plot[,3] %in% c("model 1","model 2"),]
ggplot(to.plot,aes(x=VI,y=CV,colour=model))+
  geom_point(alpha=0.4,size=1.75)+
  xlab("VIPER-Immune Correlation Coefficient")+
  ylab("CNV-VIPER Correlation Coefficient")+
  ggtitle("model 1 & 2 Distribution\nMETABRIC")

to.plot= data.frame(VI = validationTriplets[int,"VIPER.IMMUNE.Cor"],CV = validationTriplets[int,"CNV.VIPER.Cor"],model = validationTriplets[int,"model"])
to.plot = to.plot[to.plot[,3] %in% c("model 1","model 2"),]
ggplot(to.plot,aes(x=VI,y=CV,colour=model))+
  geom_point(alpha=0.4,size=1.75)+
  xlab("VIPER-Immune Correlation Coefficient")+
  ylab("CNV-VIPER Correlation Coefficient")+
  ggtitle("model 1 & 2 Distribution\nTCGA")
plot(discoveryTriplets[int,"CNV.IMMUNE.Cor"],discoveryTriplets[int,"CNV.VIPER.Cor"],xlab="CNV-Immune Correlation",ylab="CNV-VIPER Correlation",main="METABRIC model 3")

#What models have prevalent cell types members?
tmp = data.frame(discoveryTriplets[int,"CellType"],discoveryTriplets[int,"model"])
tmp = rbind(table(tmp[tmp[,2] == "model 1",1]),table(tmp[tmp[,2] == "model 2",1]),table(tmp[tmp[,2] == "model 3",1]))
tmp = tmp/rowSums(tmp)
to.plot = data.frame(model = c("model 1", "model 2", "model 3"), tmp)
colnames(to.plot) = gsub("[.]"," ",colnames(to.plot))
to.plot = melt(to.plot,id="model")
ggplot(to.plot, aes(model, value, fill=variable)) + geom_bar(stat="identity") +
  ggtitle("Proportion of Cells Belonging to each model")+
  theme(text = element_text(size=18),
        plot.title = element_text(lineheight=.8, face="bold"),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(angle = 90,size=13))+
  ylab("Proportion of Cells")

#We now use fisher p value combination to compute a score for the strength of a triplet.
#First take model 1
keep = discoveryTriplets$model == "model 1" & rownames(discoveryTriplets) %in% int
m1 = pchisq(-2 * (log(discoveryTriplets$CNV.VIPER.Pval[keep]) + log(discoveryTriplets$VIPER.IMMUNE.Pval[keep])),df=4,lower.tail = F)
m2 = discoveryTriplets$CNV.Immune.PP[keep]
m3 = pchisq(-2 * (log(m1) + log(1-m2)),df=4,lower.tail = F)

#Plot Causal Triplet Solutions
#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
OverallFisher <- -log10(m3)
qplot(-log10(m1),m2,main="model M1 Edge Correlates",xlab="-log(sumlog(CNV-VIPER.PValue,VIPER-Immune.Pvalue))",ylab="CNV-Immune.Partial.P",colour=OverallFisher)+
  scale_colour_gradient(low="red", high="blue")
discoveryM1 = cbind(discoveryTriplets[keep,],assocFisher =-log10(m1),parFisher = m2, overallFisher = -log10(m3))
q2 = discoveryM1[discoveryM1$assocFisher > 75 & discoveryM1$parFisher > 0.5,]

#How do these match up to those in TCGA?
q2vald = validationTriplets[rownames(q2),]

to.plot = data.frame(Triplet = rep(paste(q2$`CNV Pertubation`,q2$TF,q2$CellType,sep="-"),2),
                     rbind(cbind(q2[,c("CNV.VIPER.Pval","VIPER.IMMUNE.Pval","CNV.Immune.PP")],cohort = "METABRIC"),
                           cbind(q2vald[,c("CNV.VIPER.Pval","VIPER.IMMUNE.Pval","CNV.Immune.PP")],cohort = "TCGA")))
to.plot = melt(to.plot)
to.plot[,4] = -log10(to.plot[,4])
ggplot(to.plot, aes(x=variable, y=value, fill=cohort)) + geom_boxplot() +labs(size= "Nitrogen",x = "Cohort",y = "TF Activity",title = "Q2 Edge -log10(P-value) Comparison")+
  xlab("-log10(P-value)")+
  ylab("Triplet Edge Name")+
  theme(text = element_text(size=18),
        plot.title = element_text(lineheight=.8, face="bold"),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13,angle=90,hjust=1))+
  coord_flip()

#Heatmap the correlation between the TF and the Immune Cell Across METABRIC and TCGA
to.heatmap = data.frame(METABRIC= q2$VIPER.IMMUNE.Cor,TCGA = q2vald$VIPER.IMMUNE.Cor)
rownames(to.heatmap) = paste(q2$`CNV Pertubation`,q2$TF,q2$CellType,sep="-")
pheatmap(to.heatmap,main="VIPER-Immune Edge Correlation Coefficient\nMETABRIC + TCGA\n Q2 model 1")


#Heatmap the correlation between the TF and the Immune Cell
to.heatmap = data.frame(CNV_VIPER= q2$CNV.VIPER.Cor,VIPER_IMMUNE = q2$VIPER.IMMUNE.Cor)
rownames(to.heatmap) = paste(q2$`CNV Pertubation`,q2$TF,q2$CellType,sep="-")
pheatmap(to.heatmap,main="Causal Triplet Edge Correlation Coefficients\nQ2 model 1")

