####################################  Main Function ##########################################
getTripletMat <- function(dout, cibOut, cnvMatrix,backupPPI,thresh,adj){
  PPI = backupPPI
  tfs = na.omit(invtable[rownames(dout)])
  
  PPI = PPI[PPI[,1] %in% tfs,]
  PPI[,1] = table[PPI[,1]]
  PPI[,2] = table[PPI[,2]]
  PPI = PPI[!(is.na(PPI[,1]) | is.na(PPI[,2])),]
  PPI = PPI[PPI[,3] == 1,]
  
  colnames(dout) = gsub("[.]","-",colnames(dout))
  rownames(cibOut) = gsub("[.]","-",rownames(cibOut))
  colnames(cnvMatrix) = gsub("[.]","-",colnames(cnvMatrix))
  
  #Filter CIBERSORT 
  if(ncol(cibOut) > 22){
    cibOut = cibOut[cibOut$P.value <= thresh,1:22]}
  
  #Organise dataframes by matching labels
  int = intersect(colnames(cnvMatrix),colnames(dout))
  int = intersect(int, rownames(cibOut))
  dout = as.matrix(dout[,int])
  cibOut = cibOut [int,]
  cnvMatrix = as.matrix(cnvMatrix[,int])
  
  report = mclapply(1:ncol(cibOut),function(inx){
  #What CNVs correlate with the phenotype?
  cibIn = cibOut[,inx]
  cnvls = unique(PPI[,2])
  cnvls = intersect(cnvls,rownames(cnvMatrix))
  ls = rcorr(t(rbind(cibIn,cnvMatrix[cnvls,])))
  asc = data.frame(Cor = ls$r[1,-1],P = ls$P[1,-1])
  if(adj){
  asc = asc[p.adjust(asc$P,method = "BH") <= thresh,]}else{
    asc = asc[asc$P <= thresh,]
  }
  
  #What CNVs correlate with the TFs?
  ls = rcorr(t(rbind(dout,cnvMatrix[cnvls,])))
  asct = cbind(Cor = melt(ls$r[1:nrow(dout),-c(1:nrow(dout))]),P = melt(ls$P[1:nrow(dout),-c(1:nrow(dout))]))
  asct = asct[,c(1,2,3,6)] 
  asct = asct[!is.na(asct[,4]),]
  if(adj){
    asct = asct[p.adjust(asct$P,method = "BH") <= thresh,]}else{
      asct = asct[asct$P <= thresh,]
    }
  asct = asct[paste0(asct[,1],asct[,2]) %in% paste0(PPI[,1],PPI[,2]),]
  
  #What TFs correlate with the phenotype?
  ls = rcorr(t(rbind(cibIn,dout)))
  ascim = data.frame(Cor = ls$r[1,-1],P = ls$P[1,-1])
  if(adj){
    ascim = ascim[p.adjust(ascim$P,method = "BH") <= thresh,]}else{
      ascim = ascim[ascim$P <= thresh,]
    }
  
  #Filter and intersect
  cnvfinal = intersect(as.character(rownames(asc)),as.character(asct$Cor.Var2))
  tfinal = intersect(as.character(asct$Cor.Var1),as.character(rownames(ascim)))
  asct = asct[as.character(asct$Cor.Var1) %in% tfinal & as.character(asct$Cor.Var2) %in% cnvfinal,]
  
  #Begin building Triplet matrix
  colnames(asc) = c("CNV.Phen.Cor","CNV.Phen.Pval")
  colnames(ascim) = c("VIPER.Phen.Cor","VIPER.Phen.Pval")
  colnames(asct) = c("TF","CNV","CNV.TF.Cor","CNV.TF.Pval")
  asct = data.frame(asct,asc[as.character(asct[,2]),],ascim[as.character(asct[,1]),],stringsAsFactors = F)
  
  #Classify models in matrix line-by-line using loop
  outp = lapply(1:nrow(asct),function(x){
    #Declare general form for marginal bivariate distribution
    inp = c(as.character(asct[x,1]),as.character(asct[x,2]))
    l1 = function(X,Y,mu_1,mu_2,sigma_1,sigma_2,rho){
      sigma  = 2 * sigma_1^2 * (1- rho^2)
      coef = (X - mu_1 - rho * (sigma_1/sigma_2) * (Y - mu_2)^2)^2/(2 * sigma)
      1/sqrt(2 * pi * sigma) * exp(-coef)
    }
    
    #Generate probability function distributions for different CNV genotypes
    #Mutation
    alpha = scale(cnvMatrix[inp[2],])
    #No Mutation
    delta = abs(1/alpha)
    
    #Declare M1 Likelihood Model 
    M1 = function(p){-sum(log(
      (pnorm(alpha) * dnorm(dout[inp[1],],p[1] + p[2] * alpha,p[3])+
         (1-pnorm(alpha)) * dnorm(dout[inp[1],],p[1] - p[2] * alpha,p[3])+
         pnorm(scale(delta)) * dnorm(dout[inp[1],],p[1] + p[4] * delta,p[3]))*
        l1(cibIn,dout[inp[1],],p[5],p[6],p[7],p[3],p[8])))
      
    }
    #Initialise parameters
    params = c(0.5,7,sd(dout[inp[1],]),1,mean(cibIn),mean(dout[inp[1],]),sd(cibIn),0.5)
    
    #Optimise parameters using MLE and compute the AIC for model 1
    AIC1 = optim(params,M1)$value * 2 + length(params) * 2
    
    #Declare M2 Likelihood Model
    M2 = function(p){-sum(log(
      (pnorm(alpha) * dnorm(cibIn,p[1] + p[2] * alpha,p[3])+
         (1-pnorm(alpha)) * dnorm(cibIn,p[1] - p[2] * alpha,p[3])+
         pnorm(scale(delta)) * dnorm(cibIn,p[1] + p[4] * delta,p[3]))*
        l1(dout[inp[1],],cibIn,p[5],p[6],p[7],p[3],p[8])))
      
    }
    #Initialise parameters
    params = c(0.5,7,sd(cibIn),1,mean(dout[inp[1],]),mean(cibIn),sd(dout[inp[1],]),0.5)
    #Optimise parameters using MLE and compute the AIC for model 2
    AIC2 = optim(params,M2)$value * 2 + length(params) * 2
    
    #Declare M3 likelihood model
    M3 = function(p){-sum(log(
      pnorm(alpha) * dnorm(dout[inp[1],],p[1] + p[2] * alpha,p[3]) * l1(cibIn,dout[inp[1],],p[1] + p[5] * alpha,p[6],p[7],p[3],p[8])+
        (1-pnorm(alpha)) * dnorm(dout[inp[1],],p[1] - p[2] * alpha,p[3]) * l1(cibIn,dout[inp[1],],p[1] - p[5] * alpha,p[6],p[7],p[3],p[8])+
        pnorm(scale(delta)) * dnorm(dout[inp[1],],p[1] + p[4] * delta,p[3]) * l1(cibIn,dout[inp[1],],p[1] + p[9] * alpha,p[6],p[7],p[3],p[8])
    ))
      
    }
    #Initialise parameters
    params = c(0.5,7,sd(dout[inp[1],]),1,7,mean(dout[inp[1],]),sd(cibIn),0.5,0)
    #Optimise parameters using MLE and compute the AIC for model 3
    AIC3 = optim(params,M3)$value * 2 + length(params) * 2
    c(AIC1,AIC2,AIC3)})
  
  #Return triplets and their classified states
  to.ret = cbind(asct,model=unlist(lapply(outp,which.min)))
  data.frame(to.ret,Type = colnames(cibOut)[inx],stringsAsFactors = F)
  },mc.cores = 7)
  to.ret = do.call(rbind.data.frame,report)
  }
######################################End Main Function#######################################

#Network RedeR plotting function
rederPlotter = function(inputData){
#Construct the graph
g <- graph.data.frame(inputData,directed=F,vertices=NULL)
ig <- igraph::simplify(g,remove.multiple = TRUE)
#How many shared genes are there?
dupGenes = inputData[duplicated(inputData[,1]),1]
#Plot the communities using RedeR
library("RedeR")
library("igraph")
rdp <- RedPort()
calld(rdp, maxlag = 5000)
resetd(rdp)
vnames <- get.vertex.attribute(g, "name")
V(g)$nodeFontSize <- ifelse(vnames %in% unique(inputData$TF), 100,1)
V(g)$nodeSize <- ifelse(vnames %in% unique(inputData$TF), 100,30)
col = rep("gray",length(vnames))
names(col) = vnames

#Examine the TF gene interaction network and highlight duplicated genes
col[dupGenes] = "springgreen"
col[unique(inputData[,2])] = "royalblue"
V(g)$nodeColor <- col
addGraph(rdp,g,layout = NULL)}
