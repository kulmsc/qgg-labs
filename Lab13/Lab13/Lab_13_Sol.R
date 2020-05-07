Xa = as.matrix(read.table("EM_X.txt"))
Y = as.matrix(read.table("EM_Y.txt"))
Xd = ifelse(abs(X)==1, -1, 1)

library(rstan)
idx = 65
simple_data = list(N = dim(Xa)[1], xa = as.vector(Xa[,idx]), xd = as.vector(Xd[,idx]), y = as.vector(Y))
fit = stan(file = 'lab12sol.stan', data =simple_data, iter = 1000, chains = 4)
res = extract(fit)

png('histograms_ril34.png')
par(mfrow = c(2,2))
hist(res$beta_a, main = "beta_a")
hist(res$beta_d, main = "beta_d")
hist(res$sigma, main = "sigma")
hist(res$beta_mu, main = "beta_mu")
dev.off()
