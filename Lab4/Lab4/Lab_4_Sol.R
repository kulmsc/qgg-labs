library(ggplot2)
library(ggthemes)
theme_set(theme_few())
heights = read.csv("heights.csv")
data = heights$x

silly_estimator = function(mu_est, sigma, data){
  num_data_points = length(heights$x)
  prod = 1
  for(j in 1:num_data_points){
    prod = prod*1/sqrt(2*pi*sigma^2)*exp(-(data[j]-mu_est)^2/(2*sigma^2))
  }
  L =prod
  return(L)
}

#now, loop through and look at our liklihoods as a function of mu
mu_range = seq(0, 5, by =.1)
sigma = 2
Ls = c() #an empty vector
for (mu in mu_range){
  L_est = silly_estimator(mu, sigma, data)
  Ls=append(Ls, L_est)
}

#want to find max mu
best_guess_idx = which(Ls ==max(Ls))
best_guess_val = mu_range[best_guess_idx]

mu_mle_smart = mean(data)
df = data.frame(mu_range, Ls)
p1=ggplot(df)+geom_line(aes(mu_range, Ls))+geom_vline(xintercept =best_guess_val, color ="red")
p1 = p1+labs(x="mu", y="Value of Liklihood Function", title = "ril34", subtitle =paste("Best guess of mu using silly estimator", best_guess_val))
print(p1)
