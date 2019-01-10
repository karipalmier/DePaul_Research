#####################################################################################################################
#
#  GetEnrichedGO.R
#
#  This script compiles all of the GO terms for a list of genes passed in.  It then determines which terms are
#  enriched by performing a hypergeometric test on each term with the number of genes that contain that term,
#  the total number of genes in the sample, 40 as the random number of genes, and 19960 total possible genes.
#
#  Input Arguments:
#     sig_level (optional): signal level proportion to use for significance testing (for example 0.05 = 5%)
#     list of gene symbols: string gene symbol listed in arguments (each gene symbol is one argument)
#
#  Outputs:
#     JSON string of output_data R list (contents below):
#        enriched_genes: list of gene symbols present in enriched terms
#        num_total_terms: number of all terms present (before enrichment testing)
#        num_uniq_terms: number of unique terms present before enrichment testing
#        num_sample_genes: number of genes passed in
#        BP: Biological Process enriched term information - dataframe with following (rows = terms):
#                             id = term ID 
#                             term = term text 
#                             pvalue = pvalue of hypergeometric significance test of term
#                             genes = genes of corresponding term 
#                             num_genes = number of genes of the term
#                             term_count = number of times term is present in list of all terms (for all genes)
#                             type = term type (BP, CC, MF)
#        CC: Cellular Component enriched term information - dataframe with following (rows = terms):
#                             id = term ID 
#                             term = term text 
#                             pvalue = pvalue of hypergeometric significance test of term
#                             genes = genes of corresponding term 
#                             num_genes = number of genes of the term
#                             term_count = number of times term is present in list of all terms (for all genes)
#                             type = term type (BP, CC, MF)
#        MF: Molecular Function enriched term information - dataframe with following (rows = terms):
#                             id = term ID 
#                             term = term text 
#                             pvalue = pvalue of hypergeometric significance test of term
#                             genes = genes of corresponding term 
#                             num_genes = number of genes of the term
#                             term_count = number of times term is present in list of all terms (for all genes)
#                             type = term type (BP, CC, MF)
#        gene_term_list: list of term information for each gene - list with following dataframe per gene:
#                             id = term ID 
#                             term = term text 
#                             pvalue = pvalue of hypergeometric significance test of term
#                             genes = genes of corresponding term 
#                             num_genes = number of genes of the term
#                             term_count = number of times term is present in list of all terms (for all genes)
#                             type = term type (BP, CC, MF)
#        num_enriched_genes
#        num_enriched_terms
#
#
#  Date         Name          Comment
#  4/30/2018    Kari Palmier  File Created
#  5/4/2018     Kari Palmier  Updated to use # genes w term and # genes in sample for test
#  5/14/2018    Kari Palmier  Updated to return number of genes and terms per entry and overall
#  6/8/2018     Kari Palmier  Update to return list of all genes with enriched terms and list of enriched terms 
#                               per gene
#  6/13/2018    Kari Palmier  Updated to allow for first argument to be the significance level for enrichment
#  7/10/2018    Kari Palmier  Added catch to handle if GO BP, CC, or MF only has 1 entry
#  8/29/2018   Kari Palmier   File commented
#
#
#####################################################################################################################


# Only need to load the biodconductor library first run

# load bioconductor gene query library
tryCatch({
  library(mygene)
}, error=function(...){
  tryCatch({
    source("biocLite.R")
  }, error=function(...) {         # no / older BiocInstaller
    source("https://bioconductor.org/biocLite.R") 
  })
  
  biocLite("mygene")
  library(mygene)
})

# Only need to install libraries first run
library(jsonlite)

# Get command line arguements (first is gene symbol and second is the path to same the json output text file)
args = commandArgs(trailingOnly=TRUE)

temp_num = as.numeric(args[1])

if (is.na(temp_num)){
  sig_level = 0.05
  geneList = args
} else {
  sig_level = as.numeric(args[1])
  geneList = args[2:length(args)]
}

#sig_level = 0.05
#geneList = c("POLR2A", "LDHB", "NOSIP", "AQP3", "CD79A", "CD79B", "C15orf48")
# geneList = c("POLR2A", "LDHB", "NOSIP", "AQP3", "CD79A", "CD79B")

num_sample_genes = length(geneList)

gene_go_df = data.frame(id = character(), term = character(), type = character(), gene = character(),
                        stringsAsFactors = FALSE)

# Get gene information from Bioconductor database and compile a dataframe with info for all genes
for (i in 1:length(geneList)){

  ## return the query result
  geneSymbol = geneList[i]
  genequery = query(geneSymbol, size=1)

  # if entrezgene information is present, get it
  if(!is.null(genequery$hits$entrezgene)) {
    temp_gene <- getGene(genequery$hits$entrezgene, fields="all")
    
    # create list of select gene information we want
    temp_go = temp_gene[[1]]$go
    
    # Populate initial dataframe with BP information
    if (length(temp_go$BP) > 0){
      
      if (is.null(temp_go$BP$id)){
        
        for (j in 1:length(temp_go$BP)){
        
          temp_df = data.frame(id = temp_go$BP[[j]]$id, term = temp_go$BP[[j]]$term, type = "BP", 
                               gene = geneSymbol, stringsAsFactors = FALSE)
          gene_go_df = rbind(gene_go_df, temp_df)
        
        }
      } else{
        
        temp_df = data.frame(id = temp_go$BP$id, term = temp_go$BP$term, type = "BP", 
                             gene = geneSymbol, stringsAsFactors = FALSE)
        gene_go_df = rbind(gene_go_df, temp_df)
        
      }
    }
    
    # Populate initial dataframe with CC information
    if (length(temp_go$CC) > 0){
      
      if (is.null(temp_go$CC$id)){
 
        for (x in 1:length(temp_go$CC)){
          
          temp_df = data.frame(id = temp_go$CC[[x]]$id, term = temp_go$CC[[x]]$term, type = "CC", 
                               gene = geneSymbol, stringsAsFactors = FALSE)
          gene_go_df = rbind(gene_go_df, temp_df)
          
        }
      } else{
        
        temp_df = data.frame(id = temp_go$CC$id, term = temp_go$CC$term, type = "CC", 
                             gene = geneSymbol, stringsAsFactors = FALSE)
        gene_go_df = rbind(gene_go_df, temp_df)
        
      }
    }
      
    # Populate initial dataframe with MF information
    if (length(temp_go$MF) > 0){
      
      if (is.null(temp_go$MF$id)){
        
        for (y in 1:length(temp_go$MF)){
          
          temp_df = data.frame(id = temp_go$MF[[y]]$id, term = temp_go$MF[[y]]$term, type = "MF", 
                               gene = geneSymbol, stringsAsFactors = FALSE)
          gene_go_df = rbind(gene_go_df, temp_df)
          
        }
      
      } else{
        
        temp_df = data.frame(id = temp_go$MF$id, term = temp_go$MF$term, type = "MF", 
                             gene = geneSymbol, stringsAsFactors = FALSE)
        gene_go_df = rbind(gene_go_df, temp_df)
        
      }
    }
    
  } 
  
}

# Create variables required for term enrichment process
uniq_terms = unique(gene_go_df$id)
num_total_terms = length(gene_go_df$id)
num_uniq_terms = length(uniq_terms)

# Create output dataframes to be populated
enriched_bp_df = data.frame(id = character(), term = character(), pvalue = numeric(), 
                         genes = character(), num_genes = numeric(), term_count = numeric(),
                         stringsAsFactors = FALSE)
enriched_cc_df = data.frame(id = character(), term = character(), pvalue = numeric(), 
                         genes = character(), num_genes = numeric(), term_count = numeric(),
                         stringsAsFactors = FALSE)
enriched_mf_df = data.frame(id = character(), term = character(), pvalue = numeric(), 
                         genes = character(), num_genes = numeric(), term_count = numeric(),
                         stringsAsFactors = FALSE)
enriched_all_df = data.frame(id = character(), term = character(), pvalue = numeric(), 
                            genes = character(), num_genes = numeric(), term_count = numeric(),
                            type = character(), stringsAsFactors = FALSE)

enriched_genes = c()
num_random = 40
total_pop = 19960

# Test each term to see if it passes the hypergeometric test - if so, add to appropriate dataframes
num_enriched_terms = 0
for (z in 1:num_uniq_terms){
  
  term_id = uniq_terms[z]
  term_ndx = which(gene_go_df$id == term_id)
  term_count = length(term_ndx)
  
  uniq_genes = unique(gene_go_df$gene[term_ndx])
  num_genes = length(uniq_genes)
  
  # Perform hypergeometric test
  term_pval = min(1 - cumsum(dhyper(0:(num_genes - 1), num_random, total_pop, num_sample_genes)))
  
  for (w in 1:length(uniq_genes)){
    if (w == 1){
      gene_str = uniq_genes[w]
    } else {
      gene_str = paste(gene_str, ",", uniq_genes[w], sep = "")
    }
  }
  
  if (term_pval <= sig_level){
    
    enriched_genes = unique(c(enriched_genes, uniq_genes))
    num_enriched_terms = num_enriched_terms + 1

    temp_df = data.frame(id = term_id, term = gene_go_df$term[term_ndx[1]], 
                         pvalue = term_pval, genes = gene_str, num_genes = num_genes, 
                         term_count = term_count, stringsAsFactors = FALSE)
    
    temp_all_df = data.frame(id = term_id, term = gene_go_df$term[term_ndx[1]], 
                         pvalue = term_pval, genes = gene_str, num_genes = num_genes, 
                         term_count = term_count, type = gene_go_df$type[term_ndx[1]],
                         stringsAsFactors = FALSE)

    enriched_all_df = rbind(enriched_all_df, temp_all_df)
    
    if (gene_go_df$type[term_ndx[1]] == "BP"){
      
      enriched_bp_df = rbind(enriched_bp_df, temp_df)
      
    } else if (gene_go_df$type[term_ndx[1]] == "CC"){
      
      enriched_cc_df = rbind(enriched_cc_df, temp_df)
      
    } else if (gene_go_df$type[term_ndx[1]] == "MF"){
      
      enriched_mf_df = rbind(enriched_mf_df, temp_df)
      
    }
  }
  
}

gene_term_list = list()

# Create a list of dataframes, each entry corresponding to an enriched gene (to lookup term info by gene symbol)
for (i in 1:length(enriched_genes)){
  gene_name = enriched_genes[i]
  gene_ndx = which(gene_go_df$gene == gene_name)
  
  temp_gene_df = ""
  for (j in 1:length(gene_ndx)){
    term_ndx = which(enriched_all_df$term == gene_go_df$term[gene_ndx[j]])
    if (length(term_ndx) > 0){
      if (temp_gene_df == ""){
        temp_gene_df = enriched_all_df[term_ndx,]
      } else{
        temp_gene_df = rbind(temp_gene_df, enriched_all_df[term_ndx,])
      }
    }
  }
  
  gene_term_list[[gene_name]] = temp_gene_df
  
}

total_out = list(enriched_genes = enriched_genes, num_total_terms = num_total_terms, num_uniq_terms = num_uniq_terms, 
                 num_sample_genes = num_sample_genes, BP = enriched_bp_df, CC = enriched_cc_df, MF = enriched_mf_df,
                 gene_term_list = gene_term_list, num_enriched_genes = length(enriched_genes), 
                 num_enriched_terms = num_enriched_terms)
print(jsonlite::toJSON(total_out, pretty = TRUE))


