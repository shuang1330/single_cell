# This script provides annotation of the function of genes, tested in different datasets

library(biomaRt)  
library(stringr)  
# Preparing Biomart
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
genes <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",c(1:22)), mart = mart)

# Reading expression files

stemi_v3 <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/ontology/GRN_MA/expression_files/STEMI/stemi_v3_lowerres_20210301/t8w/bulk_expression.tsv',sep='\t', header=T)
stemi_v2 <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/ontology/GRN_MA/expression_files/STEMI/stemi_v2_lowerres_20210301/t8w/bulk_expression.tsv',sep='\t', header=T)
ng <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/ontology/GRN_MA/expression_files/NG/bulk_expression.tsv',sep='\t', header=T)
v2_1m <- read.table('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/cfdr_limix/expr_files/v2/UT/bulk_expression.tsv',sep='\t', header=T)
v3_1m <- read.table('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/cfdr_limix/expr_files/v3/UT/bulk_expression.tsv',sep='\t', header=T)

# output for summary
path_with_gene_annotations <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/co-eQTLs/GRN_review/gene_annotation/'

genes_all <- biomaRt::getBM(attributes = c('ensembl_gene_id','transcript_biotype'), mart = mart)
expr_genes_list <- v3_1m$X

getting_percent_protein_coding <- function(expr_genes_list, ds_id){
  genes_ds <- genes_all[genes_all$ensembl_gene_id %in%expr_genes_list, ]
  genes_ds <- genes_ds[genes_ds$transcript_biotype != 'processed_transcript',]
  genes_ds <- genes_ds[genes_ds$transcript_biotype != 'retained_intron',]
  
  genes_ds <- genes_ds[genes_ds$transcript_biotype != 'nonsense_mediated_decay',]
  
  gene_tab <- as.data.frame(table(genes_ds$transcript_biotype))
  gene_tab$perc <- gene_tab$Freq/sum(gene_tab$Freq)
  gene_tab <- gene_tab[order(gene_tab$perc, decreasing = T),]
  head(gene_tab)
  
  write.table(gene_tab, paste0(path_with_gene_annotations, ds_id, '.tsv'),sep='\t',col.names = T, row.names = F)
  
}
getting_percent_protein_coding(expr_genes_list = stemi_v3$X, ds_id= 'stemi_v3')
getting_percent_protein_coding(expr_genes_list= stemi_v2$X, ds_id= 'stemi_v2')
getting_percent_protein_coding(expr_genes_list = rownames(ng), ds_id= 'ng')
getting_percent_protein_coding(expr_genes_list = v2_1m$X, ds_id= 'v2_1m')
getting_percent_protein_coding(expr_genes_list = v3_1m$X, ds_id= 'v3_1m')

