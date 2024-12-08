rm(list=ls())
memory.limit(size=30000)  #max size in MB

load("//wolftech.ad.ncsu.edu/cos/stat/Redirect/spal4/Documents/Genomics codes and data/Core Gene Data/whole.RData")
snp=unique(whole[,4])
gene=unique(whole[,2])

require(MASS)
library(SetTest)
source("Functions.R")

gene.hc=gene.hc.pval=gene.bj=gene.bj.pval=rep(NA,length(gene))

for (k in 1:21593){
  pval = whole[whole[,2]==gene[k],5] 
  p=length(pval)
  
  if (p>1) {
    bb = test.hc(pval, M=diag(p), k0=1, k1=floor(p/2))
    gene.hc[k] = bb$hcstat
    gene.hc.pval[k] = bb$pvalue
    
    if (gene.hc.pval[k]<0.5){
      cc = test.bj(pval, M=diag(p), k0=1, k1=floor(p/2)) 
      gene.bj[k] = cc$bjstat
      gene.bj.pval[k] = cc$pvalue
    } else {
      gene.bj.pval[k] = 0.99
    }
  } else {
    print(paste(k, "gene associated with only one SNP"))
  }
  print(k)
}

bj1= gene.bj.pval[1:21593]
bj2=bj1
bj2[which(bj2<0)]=0
bj2[which(bj2>0.99)]=0.99
#bj.bh05 = BH_FDR(bj2, 0.05) #sum(bj1<bj.bh05)=928
bj.bh01 = BH_FDR(bj2, 0.01) #sum(bj1<bj.bh01)=582

top.ind = sort(bj1, index.return=T)$ix[1:sum(bj1<bj.bh01)] 
#top.ind = sort(bj1, index.return=T)$ix[1:sum(bj1<bj.bh05)] 
top.gene = gene[top.ind] #length(top.gene)=582 for fdr=0.01. 
top.gene = as.character(top.gene)

# exclude top genes with large negative bj1 
# save top.gene
write(top.gene[-c(1:4)], file="candidate_gene.txt") # 578 genes
#write(top.gene[-c(1:4)], file="candidate_gene_FDR05.txt") # 924 genes

test.bj(pval, M=diag(p), k0=1, k1=floor(p/2))

# HC based procedure select less
hc1 = gene.hc.pval[1:21593]
hc2 = hc1[hc1>0 & hc1<1]
hc.bh05 = BH_FDR(hc2, 0.05)
sum(hc2<hc.bh05) # =566

################################################################################################

## Compare with all-pair FDR of the GWAS SNPs (not significant eQTLs in full scale)
thres.all = BH_FDR(whole[,5], 0.01)
aa = sum(whole[,5]<=thres.all) 
# aa=1175 for fdr=0.05, aa=955 for fdr=0.01, aa=1389 for fdr=0.1

p.whole = sort(whole[whole[,5]<thres.all,5])
ind.whole = match(p.whole, whole[,5])
gene.whole = whole[ind.whole, 2]
gene.whole = unique((gene.whole)) # length(gene.whole)=803
# the top ones have little overlap with top ones from BJ


## Apply FDR to the set of SNPs for each significant gene
neqtl.gene = rep(NA, length(top.gene))
for(j in 1:length(top.gene)){
  pval= whole[whole[,2]==top.gene[j],5] 
  thres.gene = BH_FDR(pval, 0.01)
  neqtl.gene[j] = sum(pval<=thres.gene)
  print(j)
}
