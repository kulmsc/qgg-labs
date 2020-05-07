get_names_from_kegg_list=function(gene_list){
  list_len = length(gene_list)
  gene_names = c()
  for(j in seq(2, list_len, by =2)){
    curr_name = strsplit(gene_list[j], ";")[[1]][1] #get out the usefulpart
    gene_names = append(gene_names, curr_name)
  }
  return(gene_names)
}


circ_info =keggGet("hsa04710")
circ_genes = circ_info[[1]]$GENE #get the gene infomration out of circ_info
circ_gene_names = get_names_from_kegg_list(circ_genes) #now, we only need the gene names (the get_names_from_kegg_list function may be helpful)


#now querry ncbi to get the ids
gene_weights = c()
chromosomes =c()
for (gn in circ_gene_names){
  querry =paste("(", gn,"[Gene Name]) AND homo sapiens[Organism]",sep = "")
  res=entrez_search(db="gene", term=querry)
  #get their geneweight to plot by
  es =  entrez_summary(db="gene", res$id)
  curr_weight = es$geneweight
  curr_chrom = es$chromosome
  #make sure we have something there
  if(length(curr_weight)<1){
    curr_weight = NA
  }
  if(length(curr_chrom)<1){
    curr_chrom = NA
  }
  
  cat("gene ", gn, "weight", curr_weight,"\n")
  gene_weights = append(gene_weights, curr_weight)
  chromosomes = append(chromosomes, curr_chrom)
}

#now plot!
df = data.frame(circ_gene_names, gene_weights, as.numeric(chromosomes))
p1=ggplot(df)+geom_point(aes(chromosomes, gene_weights))
p1 = p1+geom_text(aes(x=chromosomes, y = gene_weights, label = circ_gene_names))
print(p1)
