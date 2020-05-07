#lab 11 solutions
library(MASS)
library(ggplot2)

em_expect <- function(Z, solve_A, sigma_sq_e, sigma_sq_a, X_j, beta){
  S = ginv(t(Z) %*% Z + solve_A * sigma_sq_e/sigma_sq_a) 
  a = S %*% t(Z) %*% (Y - X_j %*% beta)
  V = S * sigma_sq_e
  
  return(list(a, V))
}

em_maximize <- function(X_j, Y, Z, a, V, n, solve_A){
  beta = ginv(t(X_j) %*% X_j) %*% t(X_j) %*% (Y - Z %*% a)
  sigma_sq_a = as.numeric(1/n * (t(a) %*% solve_A %*% a + sum(diag(solve_A %*% V))))
  sigma_sq_e = as.numeric(1/n * (t(Y - X_j %*% beta - Z %*% a) %*% (Y - X_j %*% beta - Z %*% a) +     sum(diag(t(Z) %*% Z %*% V))))
  
  return(list(beta, sigma_sq_a, sigma_sq_e))
}

EM_algorithm = function(Y, X_j, A, max.iter = 100) {
  #Initiate values
  #These values are "shortcuts" for future calculations
  solve_A = ginv(A)
  n = length(Y)
  Z = diag(1, n) 
  log_L = c()
  
  #These are our random guesses of the maximized values
  sigma_sq_a = 70
  sigma_sq_e = 10
  beta = as.vector(rep(0, ncol(X_j)))
  iter = 2
  
  #Calculate initial likelihoodd
  V = A * sigma_sq_a + Z * sigma_sq_e
  #a = ginv(t(Z)%*%Z +solve_A*sigma_sq_e/sigma_sq_a)%*%t(Z)%*%(Y-X_j %*%beta)
  #useful_term = (Y-X_j %*% beta-Z %*%a)
  #log_L[1] = -n/2*log(sigma_sq_e)-n/2*log(sigma_sq_a)-1/(2*sigma_sq_e)%*%(t(useful_term))%*%useful_term-1/(2*sigma_sq_a)%*%t(a)%*%solve_A%*%a
  log_L[1] = -1/2 * determinant(V)$modulus - 1/2 * t(Y - X_j %*% beta) %*% ginv(V) %*% (Y - X_j %*% beta) 
  
  while (iter < max.iter) {
    #Expect
    expect_out <- em_expect(Z, solve_A, sigma_sq_e, sigma_sq_a, X_j, beta)
    a <- expect_out[[1]]
    V <- expect_out[[2]]
    
    #Maximize
    max_out <- em_maximize(X_j, Y, Z, a, V, n, solve_A)
    beta <- max_out[[1]]
    sigma_sq_a <- max_out[[2]]
    sigma_sq_e <- max_out[[3]]
    
    #Recalcuate log-likelihood then compare
    V = A * sigma_sq_a + Z * sigma_sq_e #really means V
    log_L[iter] = -1/2 * determinant(V)$modulus - 1/2 * t(Y - X_j %*% beta) %*% ginv(V) %*% (Y - X_j %*% beta) 
    if (log_L[iter] - log_L[iter - 1] < 1e-05) { break }
    iter = iter + 1
  }
  return(list(beta = beta, sigma_sq_a = sigma_sq_a, sigma_sq_e = sigma_sq_e, log_L = log_L[iter - 1]))
}

pval_calculator_lab7 <- function(pheno_input, xa_input, xd_input){
  n_samples <- length(xa_input)
  
  X_mx <- cbind(1,xa_input,xd_input)
  
  MLE_beta <- ginv(t(X_mx) %*% X_mx) %*% t(X_mx) %*% pheno_input
  y_hat <- X_mx %*% MLE_beta
  
  SSM <- sum((y_hat - mean(pheno_input))^2)
  SSE <- sum((pheno_input - y_hat)^2)
  
  df_M <- 2
  df_E <- n_samples - 3 
  
  MSM <- SSM / df_M
  MSE <- SSE / df_E
  
  Fstatistic <- MSM / MSE
  
  pval <- pf(Fstatistic, df_M, df_E,lower.tail = FALSE)
  
  return(pval)
}

X = as.matrix(read.table("EM_X.txt"))
Y = as.matrix(read.table("EM_Y.txt"))
A = as.matrix(read.table("EM_A.txt"))

Xd = ifelse(abs(X)==1, -1, 1)

# Null model
n_indivs = length(Y)
One = as.matrix(rep(1, n_indivs))
log_L_null = EM_algorithm(Y, One, A)$log_L

# Full model
p_values_EM = c()
times_EM = c()
for (j in 1:ncol(X)) {
  start_time <- Sys.time()
  X_j = cbind(1, X[, j])
  fit = EM_algorithm(Y, X_j, A)
  p_values_EM[j] = pchisq(-2 * (log_L_null - fit$log_L), 1, lower.tail = FALSE)
  end_time =  Sys.time()
  times_EM[j] =end_time-start_time
  cat(".") 
}



ps_linear_reg = c()
times_linear = c()
for(j in 1:ncol(X)){
  start_time= Sys.time()
  ps_linear_reg[j]= pval_calculator_lab7(Y, X[,j], Xd[,j])
  end_time = Sys.time()
  times_linear[j]=end_time-start_time
}

df = data.frame(lm = sort(ps_linear_reg), mixed = sort(p_values_EM), x_plot = seq(1, ncol(X)))
p1 = ggplot(df)+geom_point(aes(x=x_plot, y=lm, color ="linear"))
p1 = p1+geom_point(aes(x=x_plot, y = mixed, color = "mixed model"))
p1 =p1+labs(title = "ril34 \n Comparison of Linear Modeling Vs Mixed Modeling",subtitle = paste("Average time linear model", mean(times_linear), "s \n Average Time Mixed Model", mean(times_EM), "s"), xlab="id", ylab ="pvalue")
ggsave("lab11_solution_ril34.pdf")
