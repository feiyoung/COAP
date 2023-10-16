## setwd("D:\\LearnFiles\\Research paper\\idea\\3. PoisFactor\\Simu\\simu_codes")
source("./simu_util.R")



datList <- gendata_simu_lowrank_select(seed = 1,n=200, p=200, d=50, rank0 = 6, q=5, 
                                       rho=c(3, 6), sigma2_eps = 2)
q_max <- 15
d <- ncol(datList$Z)
datList <- gendata_simu(seed = 1,n=200, p=200, d=50, rank0 = 6, q=5, 
                        rho=c(3, 6), sigma2_eps = 2)
selectParams(X_count=datList$X, Z=datList$Z, verbose=F)
library(COAP)

reslist <- RR_COAP(datList$X, Z = datList$Z, rank_use = floor(d/2), q= q_max,
                   verbose = T,joint_opt_beta=F)

threshold <- 0.1
svalues <- svd(reslist$bbeta)$d
cumsum(svalues)/ sum(svalues)

par(mfrow=c(2,1))
svalues <- svalues[svalues>threshold]
ratio1 <- svalues[-length(svalues)] / svalues[-1]
which.max(ratio1)

dat1 <- data.frame(Ratio=ratio1, r=1:length(ratio1))
library(ggplot2)
p1 <- ggplot(data=dat1, aes(x=r, y=Ratio)) + geom_line(linewidth=0.8) +geom_point(size=1.8)+ 
  scale_x_continuous(breaks=seq(3, length(ratio1), by=3)) +theme_bw(base_size = 20) 
plot(ratio1, type='o', xlab='r', ylab=paste0('ratio of eigenvalue of ',expression(beta)))
abline(v=6, col='red')
which.max(ratio1[-length(ratio1)]) 
svalues <- svd(reslist$B)$d
svalues <- svalues[svalues>1e-2]
ratio_fac <- svalues[-length(svalues)] / svalues[-1]
which.max(ratio_fac)
dat2 <- data.frame(Ratio=ratio_fac, q=1:length(ratio_fac))
p2 <- ggplot(data=dat2, aes(x=q, y=Ratio)) + geom_line(linewidth=0.8) +geom_point(size=1.5) +
  theme_bw(base_size=20)
p12 <- PRECAST::drawFigs(pList=list(p1,p2), layout.dim = c(1,2), legend.position='none')

p12



np_setList <- list(c(150, 200), c(200, 200))
rhoz_set <- c(6) ## (rhoz, rhoB) = (6, 3) is good
rhoB_set <- c(3)  
i_rhoz <- 1; i_rhoB <- 1
sigma2_set <- c(2, 4)
i_np <- 2; 
i_sig <- 2


library(COAP)
N <- 200
rqMat <- matrix(NA,N, 2)
n <- np_setList[[i_np]][1]; p <- np_setList[[i_np]][2]
colnames(rqMat) <- c("r", 'q')
for(i in 1:N){
  # i <- 1
  message("i = ", i)
  datList <- gendata_simu_lowrank_select(seed = i ,n=n, p=p, d=50, rank0 = 6, q=5, 
                                         rho=c(rhoB_set[i_rhoB], rhoz_set[i_rhoz]), 
                                         sigma2_eps = sigma2_set[i_sig])
  q_max <- 15
  d <- ncol(datList$Z)
  
  
  reslist <- RR_COAP(datList$X, Z = datList$Z, rank_use = floor(d/2), q= q_max, verbose = F)
  
  threshold <- 0.1
  beta_svalues <- svd(reslist$bbeta)$d
  beta_svalues <- beta_svalues[beta_svalues>threshold]
  ratio1 <- beta_svalues[-length(beta_svalues)] / beta_svalues[-1]
  hr <- which.max(ratio1[-length(ratio1)])
  
  B_svalues <- svd(reslist$B)$d
  B_svalues <- B_svalues[B_svalues>1e-2]
  ratio_fac <- B_svalues[-length(B_svalues)] / B_svalues[-1]
  hq <- which.max(ratio_fac)
  rqMat[i, ] <- c(hr, hq)
  
}
c(mean(rqMat[,1]==6), mean(rqMat[,2]==5))

save(rqMat, file=paste0("selectqr_simu_rhoz", rhoz_set[i_rhoz], "rhoB", rhoB_set[i_rhoB],"sigma2_", sigma2_set[i_sig], "n", n, "p", p, "_R1.rds"))

