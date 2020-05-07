#lab 9 solution
library(MASS)
library(ggplot2)
pval_calculator_using_lm <- function(pheno_input, xa_input, xd_input, xz_input=NULL){
  #so this will work either with or without covariates
  if(missing(xz_input)){
    regression_df = data.frame(Y = pheno_input, Xa = xa_input, Xd = xd_input)
    model <- lm(Y ~ Xa+Xd, data = regression_df)
    s = summary(model)
    pval <- pf(as.numeric(s$fstatistic[1]),as.numeric(s$fstatistic[2]), as.numeric(s$fstatistic[3]), lower.tail = FALSE)
  } else{
    regression_df = data.frame(Y = pheno_input, Xa = xa_input, Xd = xd_input, Xz= xz_input)
    model <- lm(Y ~ Xa+Xd+Xz, data = regression_df)
    model_h0 = lm(Y~Xz,data = regression_df)
    res=anova(model, model_h0)
    pval = res$`Pr(>F)`[2]
  }
  
  return(pval)
}

hapmap.pheno.mx <- read.table("./HapMap_phenotypes.tsv", sep = "\t")
hapmap.geno.mx <- read.table("./HapMap_genotypes.tsv", sep = "\t")
hapmap.gene.info.df <- read.table("./HapMap_gene_info.tsv", sep = "\t")
hapmap.snp.info.df <- read.table("./HapMap_snp_info.tsv", sep = "\t")
xa_mat=as.matrix(hapmap.geno.mx)
xd_mat = ifelse(abs(xa_mat)==1, -1, 1)

pca.result <- prcomp(xa_mat %*% t(xa_mat))
des_covar = pca.result$x[,1]

pval_mx_w_covar <- rep(0,ncol(xa_mat))
pval_mx_no_covar <- rep(0,ncol(xa_mat))
for(i in 1:ncol(xa_mat)){
  pval_mx_w_covar[i] <- pval_calculator_using_lm(hapmap.pheno.mx[,4], xa_mat[,i], xd_mat[,i], des_covar)
  pval_mx_no_covar[i] = pval_calculator_using_lm(hapmap.pheno.mx[,4], xa_mat[,i], xd_mat[,i])
}


num_tests = ncol(xa_mat)
alpha = .05/num_tests
num_sig_ps_w_covar = length(pval_mx_w_covar[pval_mx_w_covar<alpha])
num_sig_ps_no_covar = length(pval_mx_no_covar[pval_mx_no_covar<alpha])

qqDf <- data.frame(ps=sort(-log10(pval_mx_w_covar)), normalQuantiles=sort(seq(0, 1, length.out = num_tests)))
p1 = ggplot(qqDf)+geom_point(aes(normalQuantiles, ps))
p1 = p1+ geom_abline(intercept = 0, slope = 1, color="red")
p1 = p1+ggtitle(paste("ril34 \n QQ plot with 1st PC as covar \n with alpha =", alpha, num_sig_ps_w_covar, "significant"))
ggsave("QQ_plot_w_covar.png")

qqDf <- data.frame(ps=sort(-log10(pval_mx_no_covar)), normalQuantiles=sort(seq(0, 1, length.out = num_tests)))
p2 = ggplot(qqDf)+geom_point(aes(normalQuantiles, ps))
p2 = p2+ geom_abline(intercept = 0, slope = 1, color="red")
p2 = p2+ggtitle(paste("ril34 \n QQ plot with no covars \n with alpha =", alpha, num_sig_ps_no_covar, "significant"))
ggsave("QQ_plot_no_covar.png")