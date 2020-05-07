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
  df_E <- n_samples - 6 
  
  numerator <- (SSE_theta_0-SSE_theta_1) / df_M #calculate your F statistic
  denom <- SSE_theta_1 / df_E
  Fstatistic <-numerator / denom
  
  # to check if it is correct 
  pval <- pf(Fstatistic, df_M, df_E,lower.tail = FALSE) #calculate your p value and return it
  return(pval)
}

pca.result <- prcomp(xa_mat %*% t(xa_mat))

pca.result$sdev
(pca.result$sdev / sum(pca.result$sdev))*100
summary(pca.result)

#$x contains the coordinate values of the data projected to a principal component in each of its columns.
pcaDf <- data.frame(pc1=pca.result$x[,1], pc2=pca.result$x[,2], pc3 = pca.result$x[,3])

pval_mx <- rep(0,ncol(xa_mat))
for(i in 1:ncol(xa_mat)){
  pval_mx[i] <- pval_calculator_w_covars(sim_pheno_mx[,2], xa_mat[,i], xd_mat[,i], cbind(pcaDf$pc1, pcaDf$pc2, pcaDf$pc3))
}

num_tests = ncol(xa_mat)

manhattan_df = data.frame(x = seq(1,ncol(xa_mat)), y = -log10(pval_mx))
p1=ggplot(manhattan_df)+geom_point(aes(x, y))
p1 = p1+geom_hline(yintercept = -log10(.05), color = "red")
p1 = p1+geom_hline(yintercept = -log10(.05/num_tests), color = "blue")
p1 = p1+xlab("Position")+ylab("-log p value")
ggsave("Manhattan_plot.png")

#QQ Plot
qqDf <- data.frame(ps=sort(-log10(pval_mx)), normalQuantiles=sort(seq(0, 1, length.out = num_tests)))
p2 = ggplot(qqDf)+geom_point(aes(normalQuantiles, ps))
p2 = p2+ geom_abline(intercept = 0, slope = 1, color="red")
ggsave("QQ_plot.png")