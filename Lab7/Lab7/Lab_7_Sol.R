library(ggplot2)

pval_calculator_w_covars <- function(pheno_input, xa_input, xd_input, xz_input){
  n_samples <- length(xa_input) #calculate your number of samples
  X_mx <- cbind(rep(1,length(xa_input)),xa_input, xd_input, xz_input) #create your X matrix under H1
  
  MLE_beta <- ginv(t(X_mx) %*% X_mx) %*% t(X_mx) %*% pheno_input #calculate your MLE of the betas
  
  x_h0 =  cbind(rep(1,length(xa_input)), xz_input) #calculate your x under H0
  MLE_h0 = ginv(t(x_h0) %*% x_h0) %*% t(x_h0) %*% pheno_input #calculate your MLE under h0
  y_hat_0 = x_h0 %*% MLE_h0 #calculate y_hat under the null hypothesis
  y_hat_1 = X_mx%*% MLE_beta #calculate y_hat under H1
  
  SSE_theta_0 = sum((pheno_input-y_hat_0)^2) #calculate SSE under null 
  SSE_theta_1 = sum((pheno_input-y_hat_1)^2) #calculate SSE under H1
  
  df_M <- 2
  df_E <- n_samples - 4 
  
  numerator <- (SSE_theta_0-SSE_theta_1) / df_M #calculate your F statistic
  denom <- SSE_theta_1 / df_E
  Fstatistic <-numerator / denom
  
  # to check if it is correct 
  pval <- pf(Fstatistic, df_M, df_E,lower.tail = FALSE) #calculate your p value and return it
  return(pval)
}

covar_data = read.csv("./covar_data.csv", row.names = 1) #tell R we do have row names
geno_import <- read.csv("./genotype_data.csv", 
                        header = TRUE, 
                        stringsAsFactors = FALSE,
                        row.names = 1, colClasses = "character")
codes <- genotype_coder(geno_import, 0)
xa_mat <- codes[[1]]
xd_mat <- codes[[2]]
sim_pheno_mx <- read.csv("./phenotype_data.csv", 
                         header = TRUE, row.names = 1)

pval_mx <- rep(0,ncol(xa_mat))
for(i in 1:ncol(xa_mat)){
  pval_mx[i] <- pval_calculator_w_covars(sim_pheno_mx[,1], xa_mat[,i], xd_mat[,i], covar_data$child)
}

num_sig_snps = length(pval_mx[pval_mx <=.05])

manhattan_df = data.frame(x = seq(1,ncol(xa_mat)), y = -log10(pval_mx))
p1=ggplot(manhattan_df)+geom_point(aes(x, y))
p1 = p1+geom_hline(yintercept = -log10(.05), color = "red")
p1 = p1+ggtitle(paste("ril34 \n", num_sig_snps, "significant SNPs with p=.05"))
p1 = p1+xlab("Position")+ylab("-log p value")
ggsave("Manhattan_plot.png")


