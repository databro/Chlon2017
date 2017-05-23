map = read.table("../rawdata/hg18.txt")
map = map[,c(3,4,5,6,13)]
map[,1] = gsub("chr","",map[,1])


load("../rawdata/CNA1980.RData")



  samples = unique(CNA[,1])

outp = lapply(samples,function(x){
print(length(samples) - grep(x,samples))
subCNA = CNA[CNA[,1] == x,]
keep = 0
for(i in 1:nrow(subCNA)){
submap = map[map[,1] == subCNA[i,2],]
m8=submap[submap[,3] >= subCNA[i,3] & subCNA[i,4] >= submap[,4],5]
res = rep(subCNA[i,6],length(m8))
names(res) = m8
res = res[!duplicated(names(res))]
keep = c(keep,res)}
keep[-1]})

keep = table(unlist(lapply(outp,names)))
keep = names(keep)[-1]
outp = lapply(outp,function(x){x[keep]})
outp = do.call(rbind.data.frame,outp)
colnames(outp) = keep
rownames(outp) = samples
outp[is.na(outp)] = 0
save(GAINS,file = "../rawdata/CNASeg_mean.RData")
