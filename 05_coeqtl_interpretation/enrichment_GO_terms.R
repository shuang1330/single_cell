# ------------------------------------------------------------------------------
# Check for each co-eQTL with at least 5 co-eGenes if there is any enrichment
# ------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(clusterProfiler)

# Using path on gearshift
setwd("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/")

path<-"coeqtl_mapping/output/filtered_results/"
outdir<-"coeqtl_interpretation/go_enrichment/"

#Run the GO enrichment for each cell type
enrichment<-NULL
enrichment_summary<-NULL
for(cell_type in c("CD4T","CD8T","monocyte","NK","B","DC")){

  coeqtls <- fread(paste0(path, "UT_",cell_type, 
                         "/coeqtls_fullresults_fixed.all.tsv.gz"))
  coeqtls$gene1<-gsub(";.*","",coeqtls$Gene)
  coeqtls$gene2<-gsub(".*;","",coeqtls$Gene)
  coeqtls$second_gene<-ifelse(coeqtls$gene1 == coeqtls$eqtlgen, coeqtls$gene2,
                        coeqtls$gene1)
  coeqtls$gene1<-NULL
  coeqtls$gene2<-NULL
  
  # Take all tested genes as background
  background_genes  <- union(coeqtls$eqtlgen,coeqtls$second_gene)
  
  coeqtls_sign<-coeqtls[coeqtls$gene2_isSig,]
  
  print(paste(cell_type,"with",nrow(coeqtls_sign),"co-eQTLs"))
  
  # Identify all eQTLs with at least 5 coeGenes
  coegene_count<-coeqtls_sign%>%
    group_by(snp_eqtlgene)%>%
    summarize(count_coeGenes=n())%>%
    filter(count_coeGenes>4)
  
  enrichment_found<-0
  #Perform GO enrichemt separately for each eQTL
  for(eqtl in coegene_count$snp_eqtlgene){
    
    # Run enrichment analysis with background set
    enrich_out <-enrichGO(gene=coeqtls_sign$second_gene[coeqtls_sign$snp_eqtlgene == eqtl],
                          OrgDb='org.Hs.eg.db',
                          keyType="SYMBOL",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          universe = background_genes,
                          ont="all",
                          minGSSize=5)
    
    if(nrow(enrich_out@result)>0){
      
      # Save if a enrichment was found
      enrichment_found<-enrichment_found+1
      
      # Save result dataframe
      res<-enrich_out@result
      res$cell_type<-cell_type
      res$snp_eGene<-eqtl
      enrichment<-rbind(enrichment,
                        res[,c("cell_type","snp_eGene","ID",
                               "Description","pvalue","p.adjust")])
    }

  }
  
  enrichment_summary<-rbind(enrichment_summary,
                            data.frame(cell_type,
                                       n_eqtls_freq=nrow(coegene_count),
                                       n_enrich=enrichment_found,
                                       freq_enrich=enrichment_found/nrow(coegene_count)))

}


#Format p-values
enrichment$pvalue<-format(enrichment$pvalue,digits=3)
enrichment$p.adjust<-format(enrichment$p.adjust,digits=3)
write.table(enrichment,
            file=paste0(outdir,"GOenrichment_coeGenes_allcelltypes.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
