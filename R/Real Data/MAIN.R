rm(list=ls())
memory.limit(size=30000)  #max size in MB

snp=unique(whole[,4])
#> length(snp)
#[1] 2879

gene=unique(whole[,2])
#> length(gene)
#[1] 24371

require(MASS)
library(SetTest)
source("Functions.R")

gene.hc=gene.hc.pval=gene.bj=gene.bj.pval=rep(NA,length(gene))

#gene[8765] generates errors in test.hc and test,bj. 
#gene.hc[8765]=gene.bj[8765]=0;   gene.hc.pval[8765]=gene.bj.pval[8765]=0.99

#gene[21594] and many other genes after 21593 that from cisdata 
#only associated with one snp. Generate errors in test.hc and test,bj. 
#Have checked a few, all have insignificant pval. 
#Tried mannually set bj and hc statistics as follows, not important, not used in follow up analysis
#gene.hc[21594]=gene.bj[21594]=0;   gene.hc.pval[21594]=gene.bj.pval[21594]=0.99
#Question for Karen, the trans-eGene "TPMT", "UNC13B", etc in Civelek(20) are not in the "transdata" datset

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

########
# which(gene=="DGKQ") = 15490 DGKQ is mentioned in Mohlke's paper
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


## check the significant genes and their association p-values with the SNPs
## row entries after 62337646 are cis pairs. 
# top.gene[224] is TLK1. The top eqtl p values are all for distant SNPs. 
# None of the eqtl p-values reach cis or trans significance levels (9.6*10^-6 and 3.4*10^-13 as in Mohlke's paper). This gene should not has been identified in the paper. 
# But from literature: seven genes (MSH6, TLK1, SYCP2L, BRCA1, PGAP3, DIDO1, and DDX17) were previously annotated to be responsible for the association based on distance, biological function, eQTL effect and non-synonymous SNP in high LD   
# top.gene[50] is DGKQ, which is mentioned in the paper
# top.gene[1]-[4] have negative BJ and HC values < (-1). However, they don't have very small/significant p-values. Should exclude them 
# 578 candidate genes are top.gene[5]-[582]. 

j=696
pval = whole[whole[,2]==top.gene[j],5] 
p=length(pval)
#hist(pval,n=500)
#hist(pval[pval<0.01],n=200)

whole.gene=whole[whole[,2]==top.gene[j],]
xx=sort(whole.gene[, 5], index.return=T)$ix
whole.gene[xx[1:15],]

test.bj(pval, M=diag(p), k0=1, k1=floor(p/2))

## This is about detect genes driven by multiple eQTLs signals with both cis and trans effects. 
## Check with Co-Is about the claims of novelty based on their knowledge

# HC based procedure select less
hc1 = gene.hc.pval[1:21593]
hc3 = sort(hc1, index.return=T)$ix
which(gene.hc.pval[hc3] >0)[1] # =108
#hc3[108:...] has p values>0, 
#and the top ones basically agree with top.ind of BJ 

hc2 = hc1[hc1>0 & hc1<1]
hc.bh05 = BH_FDR(hc2, 0.05)
sum(hc2<hc.bh05) # =566


## Recursive testing
#test.bj(pval, M=diag(p), k0=1, k1=floor(p/2))
test.bj(pval, M=diag(p), k0=5, k1=floor(p/2))


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

# sum(neqtl.gene)=1236 (61 more than all-pair fdr) with fdr=0.05
# sum(neqtl.gene)= 827( less than all-pair fdr) with fdr=0.01


#save.image("C:/Dropbox/Grant/NIH/2021/Aim1_FunctionalGene/metsim_21593genes_compare.RData")
#load metsim_21593genes_compare.RData






################ The following not used
#gene.hc[21595]=gene.bj[21595]=0;   gene.hc.pval[21595]=gene.bj.pval[21595]=0.99
#gene.hc[21596]=gene.bj[21596]=0;   gene.hc.pval[21596]=gene.bj.pval[21596]=0.99
#gene.hc[21597]=gene.bj[21597]=0;   gene.hc.pval[21597]=gene.bj.pval[21597]=0.99
#gene.hc[21598]=gene.bj[21598]=0;   gene.hc.pval[21598]=gene.bj.pval[21598]=0.99
#gene.hc[21599]=gene.bj[21599]=0;   gene.hc.pval[21599]=gene.bj.pval[21599]=0.99
#gene.hc[21600]=gene.bj[21600]=0;   gene.hc.pval[21600]=gene.bj.pval[21600]=0.99
#gene.hc[21601]=gene.bj[21601]=0;   gene.hc.pval[21601]=gene.bj.pval[21601]=0.99
#gene.hc[21602]=gene.bj[21602]=0;   gene.hc.pval[21602]=gene.bj.pval[21602]=0.99
#gene.hc[21603]=gene.bj[21603]=0;   gene.hc.pval[21603]=gene.bj.pval[21603]=0.99
#gene.hc[21604]=gene.bj[21604]=0;   gene.hc.pval[21604]=gene.bj.pval[21604]=0.99
#gene.hc[21605]=gene.bj[21605]=0;   gene.hc.pval[21605]=gene.bj.pval[21605]=0.99

#gene.hc[21625]=gene.bj[21600]=0;   gene.hc.pval[21600]=gene.bj.pval[21600]=0.99
#gene.hc[21660]=gene.bj[21600]=0;   gene.hc.pval[21600]=gene.bj.pval[21600]=0.99
#gene.hc[21685]=gene.bj[21600]=0;   gene.hc.pval[21600]=gene.bj.pval[21600]=0.99
#gene.hc[21686]=gene.bj[21600]=0;   gene.hc.pval[21600]=gene.bj.pval[21600]=0.99
#gene.hc[21687]=gene.bj[21600]=0;   gene.hc.pval[21600]=gene.bj.pval[21600]=0.99


#top.gene[1] is RP11-64K12.4, which has only one very small p-value of a cis SNP as follows
#62359525 ENSG00000259368.1 RP11-64K12.4 15:40660515_G/A rs35565646 4.53480e-18  0.75103000
#so are top.gene[2], LINC00202-1, and top.ind[10], RP11-17E13.2  
#Note that entries after 62,337,646 are from cisdata.

#top.gene[3] is gene HSD17B13, which remains BH significant after removing the top 2 cis SNPs
#but not significant after removing the top 3 SNPs
#                   ENSG     Gene             SNP        rsID        pval      beta
#62402725 ENSG00000170509.7 HSD17B13  4:88030261_G/T    rs442177 5.80354e-20  0.589758
#62402724 ENSG00000170509.7 HSD17B13  4:88018991_A/G   rs2035403 1.61395e-08  0.375828
#36900907 ENSG00000170509.7 HSD17B13 11:30760335_T/C   rs3925584 7.38900e-06 -0.215243
#16476395 ENSG00000170509.7 HSD17B13 4:156645513_C/A  rs13139571 1.88785e-04 -0.179942
#47121450 ENSG00000170509.7 HSD17B13 15:63831869_A/T   rs2058913 2.11675e-04  0.178571

#top.gene[4] is gene CLUAP1, which is not significant after removing the top 3 cis SNPs
#62363189  ENSG00000103351.8  CLUAP1  16:3599655_G/A  rs12448257 3.72111e-13  0.670691
#62363188  ENSG00000103351.8  CLUAP1  16:3583173_C/T   rs3751837 4.43968e-12  0.575781
#62363190  ENSG00000103351.8  CLUAP1  16:3627358_C/T    rs758747 3.55344e-09  0.492976
#52133612  ENSG00000103351.8  CLUAP1 17:34942595_G/A rs376057009 6.03512e-04 -0.165524

#top.gene[5] is gene AP3B2, which remains significant (at the boundary) 
#after removing the top 1 cis and more trans SNPs. 
#The bj values, however, remain the same, which is strange
#This gene seems to be associated with eye and nuron related diseases and DNA methylation relatedto diet
#62361859 ENSG00000103723.8         AP3B2 15:83647483_T/A rs199563536 1.90045e-20  0.893307
#51678723 ENSG00000103723.8         AP3B2 17:17420199_G/A   rs4646404 1.75452e-05  0.206409
#6438475  ENSG00000103723.8         AP3B2  2:56603985_T/C  rs13432055 1.51571e-04  0.182544
#8561761  ENSG00000103723.8         AP3B2 2:165567695_A/G   rs6705646 4.30984e-04 -0.169820
#14236032 ENSG00000103723.8         AP3B2 3:186675277_G/A   rs7645517 4.54852e-04  0.169139
#15210786 ENSG00000103723.8         AP3B2  4:65700865_G/A  rs11945861 1.04725e-03  0.158259
#10944210 ENSG00000103723.8         AP3B2  3:49088112_T/C   rs9846123 1.20348e-03 -0.156377
#3969147  ENSG00000103723.8         AP3B2 1:202116238_G/A      rs9077 1.40866e-03 -0.154222

#top.gene[6] is LINC00310, which not significant after excluding top 4 (3 cis + 1 trans) SNPs
#62395762  ENSG00000227456.3 LINC00310  21:35593827_G/A  rs28451064 1.39007e-15 -0.717855
#62395763  ENSG00000227456.3 LINC00310  21:35599128_C/T   rs9982601 3.47246e-12 -0.665816
#62395764  ENSG00000227456.3 LINC00310  21:35710290_G/A   rs2834456 1.35799e-05  0.313684
#58995468  ENSG00000227456.3 LINC00310  19:55698183_T/C  rs11671059 5.40238e-05 -0.194321
#56380752  ENSG00000227456.3 LINC00310  19:33790556_G/A  rs41355649 3.18075e-04 -0.173610

#top.gene[7] is KLF14, which not significant after excluding top 13 SNPs.
#This is an important functional gene with verified trans-eQTLs at  
#one of the largest trans-eQTL hotspots known in the human genome!
#62413730 ENSG00000174595.4 KLF14  7:130433384_C/T  rs4731702 9.63773e-15  0.518518
#62413733 ENSG00000174595.4 KLF14  7:130457931_C/G  rs1562398 2.50211e-11 -0.448859
#62413731 ENSG00000174595.4 KLF14  7:130445981_C/T  rs1364422 1.76504e-09 -0.445866
#62413732 ENSG00000174595.4 KLF14  7:130457914_A/G  rs1562396 8.76346e-07 -0.368251
#15482853 ENSG00000174595.4 KLF14   4:89054667_A/G  rs4148155 7.15970e-05  0.191176
#35063425 ENSG00000174595.4 KLF14 10:113042093_T/G rs10885122 2.06453e-04  0.178871
#16132893 ENSG00000174595.4 KLF14  4:130731284_T/C  rs4864201 3.12114e-04  0.173843
#11540183 ENSG00000174595.4 KLF14   3:57125424_G/A rs17289035 5.17937e-04 -0.167488
#17996231 ENSG00000174595.4 KLF14   5:87730027_A/G  rs7444298 5.47780e-04  0.166771
#26551090 ENSG00000174595.4 KLF14  7:100804430_C/T  rs1048365 8.80861e-04 -0.160573

#top.gene[8] is HLA-B, which is another promising multiple eQTL gene,
#but the bj values with k0>5 are the same, which is strange.
#62406742 ENSG00000234745.5 HLA-B   6:31265490_T/C  rs2247056 1.24924e-18 -0.686785
#62406744 ENSG00000234745.5 HLA-B   6:31524851_C/T  rs2857605 2.22486e-12 -0.526606
#62406740 ENSG00000234745.5 HLA-B   6:31237061_C/T  rs2894204 9.42702e-09 -0.419263
#62406741 ENSG00000234745.5 HLA-B   6:31259579_C/T  rs2524163 1.03818e-07 -0.384894
#20117456 ENSG00000234745.5 HLA-B   6:11042909_G/A  rs9393903 2.75989e-04  0.175354
#39023810 ENSG00000234745.5 HLA-B 11:116581641_G/A  rs2367970 3.06467e-04 -0.174068
#255618   ENSG00000234745.5 HLA-B   1:11866183_C/T rs13306560 4.08370e-04  0.170498
#20550798 ENSG00000234745.5 HLA-B   6:26093141_G/A  rs1800562 4.20708e-04  0.170124

#top.gene[9] is TCP11, which is significant only when excluding the top 1 cis SNP 
#62408424 ENSG00000124678.13        TCP11   6:35033854_G/A   rs820077 9.27823e-21 -1.055030
#25553500 ENSG00000124678.13        TCP11   7:44522617_T/C rs10268050 5.77732e-04  0.166087


