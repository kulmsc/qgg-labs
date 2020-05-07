#lab 10 solution
#Rachel LeCover
library(MASS)
library(ggplot2)
library(ggthemes)


gamma_inv_calc <- function(X_mx, beta_t){
  #initialize gamma
  # K is the part which goes into the exponent
  K <- X_mx %*% beta_t
  gamma_inv <- exp(K)/(1+exp(K))
  return(gamma_inv)
}

W_calc <- function(gamma_inv){
  W <- diag(as.vector(gamma_inv * (1- gamma_inv)))
  return(W)
}

beta_update <- function(X_mx, W, Y, gamma_inv, beta){
  beta_up <- beta + ginv(t(X_mx)%*%W%*%X_mx)%*%t(X_mx)%*%(Y-gamma_inv)
  return(beta_up)
}

dev_calc <- function(Y, gamma_inv){
  deviance <- 2*( sum(Y[Y==1]*log(Y[Y==1]/gamma_inv[Y==1])) + sum((1-Y[Y==0])*log((1-Y[Y==0])/(1-gamma_inv[Y==0]))) )  
  return(deviance)
}


logistic_IRLS_keep_track_betas<- function(Xa,Xd,Y = Y, beta.initial.vec = c(0,0,0), d.stop.th = 1e-6, it.max = 100) {
  #initialize the beta parameter vector at t=0
  #Create the X matrix
  X_mx <- cbind(rep(1,length(Xa)), Xa, Xd)
  beta_t <- beta.initial.vec
  
  # initialize deviance at d[t]
  dt <- 0
  beta_history = as.matrix(t(beta.initial.vec))
  
  #initialize gamma
  # K is the part which goes into the exponent
  gamma_inv <- gamma_inv_calc(X_mx, beta_t)
  
  for(i in 1:it.max) {
    dpt1 <- dt #store previous deviance
    
    # create empty matrix W
    W <- W_calc(gamma_inv)
    
    beta_t <- beta_update(X_mx, W, Y, gamma_inv, beta_t)
    beta_history = rbind(beta_history, t(beta_t))
    
    #update gamma since it's a function of beta
    
    gamma_inv <- gamma_inv_calc(X_mx, beta_t)
    #calculate new deviance
    dt <- dev_calc(Y, gamma_inv)
    
    absD <- abs(dt - dpt1)
    
    if(absD < d.stop.th) {
      cat("Convergence at iteration:", i, "at threshold:", d.stop.th, "\n")
      #logl <- loglik_calc(Y, gamma_inv)
      logl = NA
      return(list(beta_t,logl,beta_history))
    }	
  }
  cat("Convergence not reached after iteration:", i, "at threshold:", d.stop.th, "\n")
  return(list(beta_t= c(NA,NA,NA),logl=NA, beta_history))
}

pheno_data = read.table("phenotypes-lab10.tsv", sep =" ", stringsAsFactors = FALSE, header = TRUE)
geno_data = read.table("genotypes-lab10.tsv", sep = "\t", header = TRUE, row.names = 1)
Y = pheno_data$phenotype
Xa = as.matrix(geno_data)
Xd = ifelse(abs(Xa==1), 1, -1)

#general case
for(j in seq(1, ncol(Xa))){
  #print(j)
  outputs = logistic_IRLS_keep_track_betas(Xa[,j], Xd[,j], Y)
  #print(outputs[[3]])
}

#for genotype #32
idx = 32
outputs = logistic_IRLS_keep_track_betas(Xa[,idx], Xd[,idx], Y)
beta_hisory =data.frame(outputs[[3]], seq(1, nrow(outputs[[3]])))
colnames(beta_hisory) = c("beta_mu", "beta_a", "beta_d", "iter")

p1 = ggplot(data = beta_hisory)
p1 = p1+geom_point(aes(x = iter, y = beta_mu, color = "beta_mu"))
p1 = p1+geom_point(aes(x = iter, y = beta_a, color = "beta_a"))
p1 = p1+geom_point(aes(x = iter, y = beta_d, color = "beta_d"))
p1 = p1+labs(title = "ril34 \n Beta Development", x= "Iteration", y="Beta Estimate")
ggsave("ril34_beta_dev.pdf")
