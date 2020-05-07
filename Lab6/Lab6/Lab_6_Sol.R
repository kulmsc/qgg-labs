

xd_converter_from_xa = function(xa_matrix){
  num_pop = nrow(xa_matrix)
  num_genes = ncol(xa_matrix)
  xd =matrix(nrow = nrow(xa_matrix), ncol = ncol(xa_matrix))
  for(j in 1:num_pop){
    xd[j,]=ifelse(abs(xa_matrix[j,]), -1, 1)
  }
  return(xd)
}

xd_converter <- function(geno_in){
  xd_result <- ifelse(geno_in[,1]==geno_in[,2], -1, 1)
  return(xd_result)
}

filter_alleles=function(geno_dat, threshold=.1){
  to_drop = c()
  for(i in 1:(ncol(geno_import)/2)){
    geno_count <- table(c(geno_dat[,2*i-1],geno_dat[,2*i]))
    minor_allele <- names(geno_count[geno_count == min(geno_count)])[1] #grab first element incase we end up with a tie
    minor_allele_freq = geno_count[minor_allele]/(2*nrow(geno_dat))
    cat("i ", i,"minor allele" , minor_allele, "freq", minor_allele_freq, "\n")
    if(minor_allele_freq<threshold){
      to_drop = append(to_drop, 2*i-1)
      to_drop = append(to_drop, 2*i)
    }
    
  }
  return(geno_dat[,-to_drop])
}

xd_matrix_from_xa = xd_converter_from_xa(xa_matrix)

xd_matrix <- matrix(NA, nrow = nrow(geno_import), ncol = ncol(geno_import)/2)

count = 1
for (i in seq(2, ncol(geno_import), by =2)){
  xd_matrix[,count] <- xd_converter(geno_import[,c(i-1, i)])
  count = count+1
}

filtered_geno_dat = filter_alleles(geno_import, .1)
xd_matrix_filtered = matrix(NA, nrow =nrow(filtered_geno_dat), ncol =ncol(filtered_geno_dat)/2)
count = 1
for (i in seq(2, ncol(filtered_geno_dat), by =2)){
  xd_matrix_filtered[,count] <- xd_converter(filtered_geno_dat[,c(i-1, i)])
  count = count+1
}

p1=ggplot(data.frame(y=xd_matrix_filtered[10,], x = seq(1, ncol(xd_matrix_filtered))))+geom_point(aes(x, y))
p1 = p1+labs(x="genome position", y ="Xd", title = paste("Patient 10 \n ril34 \n ", ncol(xd_matrix_filtered), "genes remaining post filtering"))
