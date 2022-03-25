# ------------------------------------------------------------------------------
# Check Pearson correlation between cell types
# Correlation matrices generated with correlation_timepoint_combined_indivs.py
# Checking only UT cells and genes expressed in at least one cell type
#  -----------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(viridis)

theme_set(theme_bw())

setwd("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/")

cell_types<-c("CD4T","CD8T","monocyte","NK","DC","B")

#Full cell type names as reported in the paper
cell_types_corrected<-setNames(c("CD4+ T","CD8+ T","Monocyte","NK","DC","B"),
                               c("CD4T","CD8T","monocyte","NK","DC","B"))

#Check for both v3 and v2 chemistry
path_v3<-"co-expression_indivs_combined/"
path_v2<-"co-expression_indivs_combined/one_million_version2/"

for(version in c("v2","v3")){
  
  if(version=="v2"){
    path<-path_v2
  } else {
    path<-path_v3
  }
  
  corr_comp<-NULL
  for(c1 in 1:c(length(cell_types)-1)){
    
    #Read correlation file one
    cell_type1<-cell_types[c1]
    corr_c1<-fread(paste0(path,cell_type1,"/",cell_type1,"_UT_correlation.csv"))
    
    #Unique genes
    num_genes<-length(union(gsub(";.*","",corr_c1$V1),
                             gsub(".*;","",corr_c1$V1)))
    
    corr_comp<-rbind(corr_comp,
                     data.frame(c1=cell_type1,
                                c2=cell_type1,
                                gene_pairs=nrow(corr_c1),
                                genes_unique=num_genes,
                                corr=1))
    
    for(c2 in (c1+1):length(cell_types)){
      
      #Read correlation file two
      cell_type2<-cell_types[c2]
      corr_c2<-fread(paste0(path,cell_type2,"/",cell_type2,"_UT_correlation.csv"))
      
      corr<-merge(corr_c1,corr_c2,by=c("V1"))
      
      #Unique genes
      num_genes<-length(union(gsub(";.*","",corr$V1),
                              gsub(".*;","",corr$V1)))
      
      corr_comp<-rbind(corr_comp,
                       data.frame(c1=cell_type1,
                                  c2=cell_type2,
                                  gene_pairs=nrow(corr),
                                  genes_unique=num_genes,
                                  corr=cor(corr$UT.x,corr$UT.y,method="pearson")))
    }
  }
  
  #Add last diagonal entry
  cell_type1<-cell_types[length(cell_types)]
  corr_c1<-fread(paste0(path,cell_type1,"/",cell_type1,"_UT_correlation.csv"))
  
  #Unique genes
  num_genes<-length(union(gsub(";.*","",corr_c1$V1),
                          gsub(".*;","",corr_c1$V1)))
  
  corr_comp<-rbind(corr_comp,
                   data.frame(c1=cell_type1,
                              c2=cell_type1,
                              gene_pairs=nrow(corr_c1),
                              genes_unique=num_genes,
                              corr=1))
  
  #Rename cell types to make it coherent with other part of the manuscript
  corr_comp$c1<-cell_types_corrected[corr_comp$c1]
  corr_comp$c2<-cell_types_corrected[corr_comp$c2]
  corr_comp$c1<-factor(corr_comp$c1,levels=cell_types_corrected)
  corr_comp$c2<-factor(corr_comp$c2,levels=cell_types_corrected)
  
  #Create heatmap
  g<-ggplot(corr_comp,aes(x=c1,y=c2,fill=corr))+
    geom_tile()+
    geom_text(aes(label=paste0(round(corr,3),"\n(",genes_unique,")")),size=4)+
    xlab("Cell type")+
    ylab("Cell type")+
    #scale_fill_gradient2("Correlation",limits = c(-1,1),low="darkblue",mid="white",high="darkred")
    scale_fill_viridis("Correlation",limits=c(0,1))+
    theme(axis.title = element_text(size=14),
          axis.text = element_text(size=12),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12))
  ggsave(g,file=paste0("co-expression_indivs_combined/plots/corr_celltypes_",version,".pdf"),
         width=7,height=5)
  
  #Check correlation distribution across cell types
  summary(corr_comp$corr[corr_comp$c1 != corr_comp$c2])
  
}