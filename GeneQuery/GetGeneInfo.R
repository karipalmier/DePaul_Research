#####################################################################################################################
#
#  GetGeneInfo.R
#
#  This code gets a gene symbol and output file path from the command line arguments present when the file is called
#  using the Rscript command.  It then gets the Bioconductor gene information if the gene symbol has an
#  entrezgene present.  If the entrezgene is present, this code creates a JSON text file with the genomic_pos, 
#  genomic_pos_hg19, type_of_gene, pathway, generif, other_names, name, summary, wikipedia, and go (genetic
#  oncology) entries from the entrezgene query return.
#
#  Date     Name          Comment
#  1/30/18  Kari Palmier  File Created and integrated
#  2/6/18   Kari Palmier  Updated to print JSON output to STDOUT
#
#
#
#####################################################################################################################


# Only need to load the biodconductor library first run
#tryCatch({
#  #throw(1)
#  source("https://bioconductor.org/biocLite.R") 
#}, error=function(...) {         # no / older BiocInstaller
#   print("reading from local copy")
#   source("biocLite.R")
#})

# Only need to install libraries first run
#biocLite("mygene") 

# load bioconductor gene query library
library(mygene)

# Get command line arguements (first is gene symbol and second is the path to same the json output text file)
args = commandArgs(trailingOnly=TRUE)
geneSymbol = args[1]

## return the query result
genequery = query(geneSymbol, size=1)

# if entrezgene information is present, get it
if(!is.null(genequery$hits$entrezgene)) {
  gene <- getGene(genequery$hits$entrezgene, fields="all")
  
  # create list of select gene information we want
  temp_gene = list(genomic_pos = gene[[1]]$genomic_pos,
                   genomic_pos_hg19 = gene[[1]]$genomic_pos_hg19,
                   type_of_gene = gene[[1]]$type_of_gene,
                   pathway = gene[[1]]$pathway,
                   generif = gene[[1]]$generif,
                   other_names = gene[[1]]$other_names,
                   name = gene[[1]]$name,
                   summary = gene[[1]]$summary,
                   wikipedia = gene[[1]]$wikipedia,
                   go = gene[[1]]$go
                   )
} else {
  temp_gene = list(genomic_pos = {},
                   genomic_pos_hg19 = {},
                   type_of_gene = {},
                   pathway = {},
                   generif = {},
                   other_names = {},
                   name = {},
                   summary = {},
                   wikipedia = {},
                   go = {}
  )
}

print(jsonlite::toJSON(temp_gene, pretty = TRUE))

