# ------------------------------------------------------------------------------
# Estimate for the 50% cutoff, how many genes overlap between cell types and
# time points
#
# Requirement: correlation values were calculated with % cutoff separately
#              estimated for each cell type and timepoints
# ------------------------------------------------------------------------------

library(data.table)
library(ggplot2)

theme_set(theme_bw())

setwd("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/")

cell_types<-c("B","CD4T","CD8T","DC","monocyte","NK")
timepoints<-c("UT","3h","24h")

num_genes<-NULL
intersect_union_genes<-NULL
#union_union_genes<-NULL
for(ct in cell_types){
  
  union_expressed_genes<-NULL
  for(tp in timepoints){
    #Load correlation values (filtered)
    corr_ct<-fread(paste0("co-expression_indivs_combined/",ct,"/",ct,"_",tp,
                          "_correlation.csv"))
    expressed_genes<-union(gsub(";.*","",corr_ct$V1),gsub(".*;","",corr_ct$V1))
    
    num_genes<-rbind(num_genes,
                     data.frame(ct,tp,ngenes=length(expressed_genes)))
    
    union_expressed_genes<-union(union_expressed_genes,expressed_genes)
  }
  
  #Save file
  write.table(data.frame(genes=union_expressed_genes),
              file=paste0("co-expression_indivs_combined/coexp_tp_union/genelists_tp_union/",
                          "expressed_gene_",ct,".tsv"),
              quote=FALSE,row.names=FALSE)
  
  #Create the union of all time points
  num_genes<-rbind(num_genes,
                   data.frame(ct,tp="union",ngenes=length(union_expressed_genes)))
  
  if(is.null(intersect_union_genes)){
    intersect_union_genes<-union_expressed_genes
  } else {
    intersect_union_genes<-intersect(intersect_union_genes,union_expressed_genes)
  }
  
  #union_union_genes<-union(union_union_genes,union_expressed_genes)
}

#Save file
write.table(data.frame(genes=intersect_union_genes),
            file=paste0("co-expression_indivs_combined/coexp_tp_union/genelists_tp_union/",
                        "expressed_gene_intersect_cts.tsv"),
            quote=FALSE,row.names=FALSE)

num_genes<-rbind(num_genes,
                 data.frame(ct="intersect",tp="union",
                            ngenes=length(intersect_union_genes)))#,
                 # data.frame(ct="union",tp="union",
                 #            ngenes=length(union_union_genes)))

#Create barplot with number of expressed genes
num_genes$tp<-factor(num_genes$tp,levels=c(timepoints,"union"))
num_genes$ct<-factor(num_genes$ct,levels=c(cell_types,"union","intersect"))

g<-ggplot(num_genes,aes(x=ct,y=ngenes,fill=tp))+
  geom_bar(stat="identity",position="dodge")+
  xlab("Cell type")+ylab("Number genes")+
  scale_fill_discrete("Time point")
ggsave(g,file="co-expression_indivs_combined/plots/expressed_genes.png",width=6,height=3)
ggsave(g,file="co-expression_indivs_combined/plots/expressed_genes.pdf",width=6,height=3)
