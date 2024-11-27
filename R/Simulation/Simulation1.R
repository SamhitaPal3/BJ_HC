library(doParallel)
library(expm)

# Loading the SNP data
snp_dat       <- load("C:/Users/HP/Downloads/Genomics/Data/SNP_data_HapMap_ch21_n90.RData")

# Selecting 2000 SNPs
X             <- X[,(1:2000)]

num_of_sub    <- nrow(X)
num_of_snp    <- ncol(X)
cor_mat       <- cor(X)
num_of_genes  <- 100L
num_of_core   <- 20L
cis_snp       <- 2

# Numner of active SNPs
tau           <- sample(((cis_snp+1):num_of_snp), 1000)

# Sets of trans SNPs
set1          <- sample(tau, 500)
set2          <- sample(tau, 500)
set3          <- sample(tau, 500)
set4          <- sample(tau, 500)
set5          <- sample(tau, 500)
set6          <- sample(tau, 500)
set7          <- sample(tau, 500)
set8          <- sample(tau, 500)
set9          <- sample(tau, 500)
set10         <- sample(tau, 500)
set11         <- sample(tau, 500)
set12         <- sample(tau, 500)
set13         <- sample(tau, 500)
set14         <- sample(tau, 500)
set15         <- sample(tau, 500)
set16         <- sample(tau, 500)
set17         <- sample(tau, 500)
set18         <- sample(tau, 500)
set19         <- sample(tau, 500)

iters = 100
A_s = 4.25 # Mean of the Normal distribution for Strong Cis effects
A_w = 0.2  # Mean of the Normal distribution for Weak Trans effects
sig_s = sig_w = 1.05 # Signal Variability

cl = makeCluster(4); registerDoParallel(cl)

# Creating test statistics Z for each Gene (in columns)
mu_noise      <- matrix(rep(0,num_of_snp*(100-num_of_core)), nrow = num_of_snp)
list_mu <- foreach(i = 1:iters) %:% foreach(ss=sig_s) %:%
  foreach(sw=sig_w) %dopar% {
    mu_strong       <- c(rnorm(1,A_s,ss), rep(0,(num_of_snp-1)))
    mu_weak1        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak1[set1]  <- rnorm(length(set1),A_w,sw) 
    mu_weak2        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak2[set2]  <- rnorm(length(set2),A_w,sw)
    mu_weak3        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak3[set4]  <- rnorm(length(set1),A_w,sw) 
    mu_weak4        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak4[set5]  <- rnorm(length(set1),A_w,sw) 
    mu_weak5        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak5[set6]  <- rnorm(length(set1),A_w,sw) 
    mu_weak6        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak6[set7]  <- rnorm(length(set1),A_w,sw) 
    mu_weak7        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak7[set8]  <- rnorm(length(set1),A_w,sw)
    mu_weak8        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak8[set9]  <- rnorm(length(set1),A_w,sw)
    mu_weak9        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak9[set10]  <- rnorm(length(set1),A_w,sw)
    mu_weak10        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak10[set11]  <- rnorm(length(set1),A_w,sw)
    mu_weak11        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak11[set12]  <- rnorm(length(set1),A_w,sw)
    mu_weak12        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak12[set13]  <- rnorm(length(set1),A_w,sw)
    mu_weak13        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak13[set14]  <- rnorm(length(set1),A_w,sw)
    mu_weak14        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak14[set15]  <- rnorm(length(set1),A_w,sw)
    mu_weak15        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak15[set16]  <- rnorm(length(set1),A_w,sw)
    mu_weak16        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak16[set17]  <- rnorm(length(set1),A_w,sw)
    mu_weak17        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak17[set18]  <- rnorm(length(set1),A_w,sw)
    mu_weak18        <- c(0, rnorm(1,A_w,sw), rep(0,num_of_snp-2))
    mu_weak18[set19]  <- rnorm(length(set1),A_w,sw)
    mu_mixed        <-  c(0, rnorm(1,A_s,ss), rep(0,num_of_snp-2))
    mu_mixed[set3]  <- rnorm(length(set3),A_w,sw)
    
    mu              <- cbind(mu_strong,mu_mixed,mu_weak1,mu_weak2,
                             mu_weak3,mu_weak4,mu_weak5,mu_weak6,
                             mu_weak7,mu_weak8,mu_weak9,mu_weak10,
                             mu_weak11,mu_weak12,mu_weak13,mu_weak14,
                             mu_weak15,mu_weak16,mu_weak17,mu_weak18,
                             mu_noise) 
  }

mu = lapply(1:iters, function(i){list_mu[[i]][[1]][[1]]})

S = cor_mat %^% (1/2)
Y = foreach(i = 1:iters) %dopar% {matrix(rnorm(2000*100),ncol = 100)}

summary_stat = foreach(i = 1:iters) %dopar% {mu[[i]]+S%*%Y[[i]]}

clusterExport(cl, c("summary_stat", "cor_mat", "mu"))

pval = foreach(i = 1:iters) %dopar% {
  sapply(1:num_of_genes, function(k) {
    (2*pnorm(abs(summary_stat[[i]][,k]), lower.tail = FALSE))})}

install.packages("C:/Users/HP/Downloads/Genomics/Papers/SetTest_0.3.0.tar.gz", 
                 repos = NULL, 
                 type = "source")

library(SetTest)

min.pvalue    <- foreach(i = 1:iters) %dopar% {apply(pval[[i]],2,min)}

bj.stat       <- foreach(i = (1:iters)) %dopar% {
  sapply(1:num_of_genes, function(k) {
    SetTest::stat.bj(p=sort(pval[[i]][,k]),k0 = 1, k1 = 1000)$value })}

hc.stat        <- foreach(i = (1:iters)) %dopar% {
  sapply(1:num_of_genes, function(k) {
    SetTest::stat.hc(p=pval[[i]][,k],k0 = 2, k1 = 1000)$value })}

mean.stat     <- foreach(i = 1:iters) %dopar% {
  apply(abs(summary_stat[[i]]),2,mean)}

bj.rank       <-  foreach(i = 1:iters) %dopar% {
  sort(bj.stat[[i]], decreasing = T, 
       index.return = T) }

hc.rank       <-  foreach(i = 1:iters) %dopar% {
  sort(hc.stat[[i]], decreasing = T, 
       index.return = T) }

mean.rank     <- foreach(i = 1:iters) %dopar% {
  sort(mean.stat[[i]], decreasing = T,
       index.return = T) }

min.rank      <- foreach(i = 1:iters) %dopar% {
  sort(min.pvalue[[i]], decreasing = F,
       index.return = T) }

active = 1:20

Recall.BJ = foreach(i = 1:iters) %:% 
  foreach(k = 1:100) %dopar% {
    length(intersect(bj.rank[[i]]$ix[1:k], 
                     active))/20 }

Recall.HC = foreach(i = 1:iters) %:%
  foreach(k = 1:100) %dopar% {
    length(intersect(hc.rank[[i]]$ix[1:k], 
                     active))/20 }

Recall.mean = foreach(i = 1:iters) %:%
  foreach(k = 1:100) %dopar% {
    length(intersect(mean.rank[[i]]$ix[1:k], 
                     active))/20 }

Recall.min = foreach(i = 1:iters) %:%
  foreach(k = 1:100) %dopar% {
    length(intersect(min.rank[[i]]$ix[1:k], 
                     active))/20 }

BJ_mean_recall = colMeans(do.call(rbind, 
                                  lapply(1:iters, function(k) {
                                    unlist(Recall.BJ[[k]])}))) 


HC_mean_recall = colMeans(do.call(rbind, 
                                  lapply(1:100, function(k) {
                                    unlist(Recall.HC[[k]])}))) 

mean_mean_recall = colMeans(do.call(rbind, 
                                    lapply(1:100, function(k) {
                                      unlist(Recall.mean[[k]])}))) 

min_mean_recall = colMeans(do.call(rbind, 
                                   lapply(1:100, function(k) {
                                     unlist(Recall.min[[k]])}))) 

Prec.BJ   = foreach(i = 1:iters) %:% 
  foreach(k = 1:100) %dopar% {
    length(intersect(bj.rank[[i]]$ix[1:k], active))/k}

Prec.HC   = foreach(i = 1:iters) %:% 
  foreach(k = 1:100) %dopar% {
    length(intersect(hc.rank[[i]]$ix[1:k], active))/k}

Prec.min   = foreach(i = 1:iters) %:% 
  foreach(k = 1:100) %dopar% {
    length(intersect(min.rank[[i]]$ix[1:k], active))/k}

Prec.mean   = foreach(i = 1:iters) %:% 
  foreach(k = 1:100) %dopar% {
    length(intersect(mean.rank[[i]]$ix[1:k], active))/k}


BJ_mean_prec = colMeans(do.call(rbind, 
                                lapply(1:100, function(k) {
                                  unlist(Prec.BJ[[k]])}))) 

HC_mean_prec = colMeans(do.call(rbind, 
                                lapply(1:100, function(k) {
                                  unlist(Prec.HC[[k]])}))) 

min_mean_prec = colMeans(do.call(rbind, 
                                 lapply(1:100, function(k) {
                                   unlist(Prec.min[[k]])}))) 

mean_mean_prec = colMeans(do.call(rbind, 
                                  lapply(1:100, function(k) {
                                    unlist(Prec.mean[[k]])}))) 


plot(BJ_mean_recall,BJ_mean_prec, 
     type = "l", col = "darkslateblue", ylab = "Precision", 
     xlab = "Recall", lwd = 2, main = expression(tau[k]^2*'=1.05, s=500, B=0.2'))
lines(mean_mean_recall,mean_mean_prec, type = "l", 
      col = "darkolivegreen4", lwd = 2, lty = 2)
lines(min_mean_recall,min_mean_prec, type = "l", 
      col = "darkgoldenrod2", lwd = 2, lty = 6)
lines(HC_mean_recall,HC_mean_prec, type = "l", 
      col = "deeppink3", lwd = 2, lty = 3)

legend("bottomleft",
       legend = c("Mean-based","BJ", "HC","MinimumP"), 
       col=c("darkolivegreen4","darkslateblue","deeppink3","darkgoldenrod2"),
       lty = c(2,1,3,6),
       lwd=2, cex=.7, horiz = F, xpd = NA)
