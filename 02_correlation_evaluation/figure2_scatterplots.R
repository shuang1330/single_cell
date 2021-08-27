################################################################################
# Create scatterplots for Figure 2 (for gearshift)
# * for BIOS
# * for Blueprint
# * the third data set
################################################################################

library(data.table)
library(reticulate) # to read the single cell data (numpy)
library(reshape2)
library(ggplot2)
library(viridis)
library(ggpubr)

use_python("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/tools/miniconda3/envs/r40/bin/python")

np <- import("numpy")

theme_set(theme_bw())

setwd("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing")

# Compare the different single cell data sets
create_scatter_plot_sc<-function(sc_ct, path, xlab_text,ylab_text,plot_path){
  
  # Load single cell data (1 Million cell - v3)
  corr_ct<-fread(paste0("co-expression_indivs_combined/",sc_ct,"/",
                        sc_ct,"_UT_correlation.csv"))
  corr_ct$gene1<-sapply(corr_ct$V1,function(s) strsplit(s,";")[[1]][1])
  corr_ct$gene2<-sapply(corr_ct$V1,function(s) strsplit(s,";")[[1]][2])
  corr_ct$V1<-NULL
  
  #Order so that gene1 is always the one first in alphabet
  corr_ct$swap<-ifelse(corr_ct$gene1 > corr_ct$gene2,corr_ct$gene1,corr_ct$gene2)
  corr_ct$gene1<-ifelse(corr_ct$gene1 > corr_ct$gene2,corr_ct$gene2,corr_ct$gene1)
  corr_ct$gene2<-corr_ct$swap
  corr_ct$swap<-NULL
  
  #Load the STEMI data set
  corr_bios<-fread(path)
  corr_bios$gene1<-sapply(corr_bios$V1,function(s) strsplit(s,";")[[1]][1])
  corr_bios$gene2<-sapply(corr_bios$V1,function(s) strsplit(s,";")[[1]][2])
  corr_bios$V1<-NULL
  
  #Order so that gene1 is always the one first in alphabet
  corr_bios$swap<-ifelse(corr_bios$gene1 > corr_bios$gene2,corr_bios$gene1,corr_bios$gene2)
  corr_bios$gene1<-ifelse(corr_bios$gene1 > corr_bios$gene2,corr_bios$gene2,corr_bios$gene1)
  corr_bios$gene2<-corr_bios$swap
  corr_bios$swap<-NULL
  corr_bios<-corr_bios[corr_bios$gene1!=corr_bios$gene2,]
  colnames(corr_bios)[1]<-"corr"
  
  #Merge both
  corr_bios<-merge(corr_bios,corr_ct,by=c("gene1","gene2"))
  
  print(paste("Overlapping genes:",length(union(corr_bios$gene1,
                                                corr_bios$gene2))))
  
  corr_corr<-cor(corr_bios$UT,corr_bios$corr,
                 method="pearson")
  
  #Plot
  g<-ggplot(corr_bios,aes(UT,corr))+
    geom_bin2d(bins=50)+
    xlab(xlab_text)+
    ylab(ylab_text)+
    xlim(-1,1)+ylim(-1,1)+
    #scale_fill_viridis(trans="log10")+
    scale_fill_continuous("Density",type = "viridis",
                          trans="log10",
                          breaks = c(5, 1500), 
                          labels = c("Low", "High"))+
    annotate(geom="text", x=-0.95, y=0.95,size=4,
             hjust = 0,vjust=1,
             label=paste0("r = ",format(corr_corr,digits=2)))+
    theme(axis.title=element_text(size=12),
          axis.text=element_text(size=11))
  #ggsave(g,file=plot_path)
  return(g)
}

g_1<-create_scatter_plot_sc(sc_ct="monocyte",
                            path=paste0("co-expression_indivs_combined/one_million_version2/monocyte/",
                                       "monocyte_UT_correlation.csv"),
                            xlab_text="Correlation 1 Mio cell V3",
                            ylab_text="Correlation 1 Mio cell V2",
                            plot_path=NULL)

g_2<-create_scatter_plot_sc(sc_ct="monocyte",
                            path=paste0("co-expression_indivs_combined/stemi/version2/monocyte/",
                                        "monocyte_Baseline_correlation.csv"),
                           xlab_text="Correlation 1 Mio cell V3",
                           ylab_text="Correlation STEMI V2",
                           plot_path=NULL)

g_3<-create_scatter_plot_sc(sc_ct="monocyte",
                            path=paste0("co-expression_indivs_combined/stemi/version3/monocyte/",
                                        "monocyte_Baseline_correlation.csv"),
                            xlab_text="Correlation 1 Mio cell V3",
                            ylab_text="Correlation STEMI V3",
                            plot_path=NULL)

g_4<-create_scatter_plot_sc(sc_ct="monocyte",
                            path=paste0("co-expression_indivs_combined/ng_updated_version/monocyte/",
                                        "monocyte_correlation.csv"),
                            xlab_text="Correlation 1 Mio cell V3",
                            ylab_text="Correlation pilot",
                            plot_path=NULL)

g<-ggarrange(g_1,g_2,g_3,g_4,ncol=4,nrow=1,common.legend=TRUE,legend="bottom")
ggsave(g,file="bios/plots/figure2_sc_combined.png",
       width=15,height=4)
ggsave(g,file="bios/plots/figure2_sc_combined.pdf",
       width=15,height=4)

#Iterate over all other cell types
plot_list<-list()
for(ct in c("CD4T","CD8T","B","NK","DC")){
  
  g_1<-create_scatter_plot_sc(sc_ct=ct,
                              path=paste0("co-expression_indivs_combined/one_million_version2/",
                                          ct,"/",ct,"_UT_correlation.csv"),
                              xlab_text=paste("Corr 1 Mio cell V3",ct),
                              ylab_text="Correlation 1 Mio cell V2",
                              plot_path=NULL)
  
  g_2<-create_scatter_plot_sc(sc_ct=ct,
                              path=paste0("co-expression_indivs_combined/stemi/version2/",
                                          ct,"/",ct,"_Baseline_correlation.csv"),
                              xlab_text=paste("Corr 1 Mio cell V3",ct),
                              ylab_text=paste("Corr STEMI V2",ct),
                              plot_path=NULL)
  
  g_3<-create_scatter_plot_sc(sc_ct=ct,
                              path=paste0("co-expression_indivs_combined/stemi/version3/",
                                          ct,"/",ct,"_Baseline_correlation.csv"),
                              xlab_text=paste("Corr 1 Mio cell V3",ct),
                              ylab_text=paste("Corr STEMI V3",ct),
                              plot_path=NULL)
  
  g_4<-create_scatter_plot_sc(sc_ct="monocyte",
                              path=paste0("co-expression_indivs_combined/ng_updated_version/monocyte/",
                                          "monocyte_correlation.csv"),
                              xlab_text=paste("Corr 1 Mio cell V3",ct),
                              ylab_text=paste("Corr pilot",ct),
                              plot_path=NULL)
  
  plot_list<-c(plot_list,list(g_1,g_2,g_3,g_4))
}

g<-ggarrange(plotlist = plot_list,ncol=4,nrow=5,common.legend=TRUE,legend="bottom")
ggsave(g,file="bios/plots/figure2_sc_combined_supplement.png",
       width=15,height=20)
ggsave(g,file="bios/plots/figure2_sc_combined_supplement.pdf",
       width=15,height=20)

################################################################################
# Compare single cell and bulk data (or between different bulk data set)
create_scatter_plot<-function(path,plot_path,ylab_text,
                           sc_ct,
                           rowname_suffix="rows.txt",
                           colname_suffix="cols.txt",
                           scData=TRUE,
                           xlab_text="Correlation single cell",
                           data_type="numpy"){
  
  # Load single cell data
  if(scData){
    corr_ct<-fread(paste0("co-expression_indivs_combined/",sc_ct,"/",sc_ct,"_UT_correlation.csv"))
    corr_ct$gene1<-sapply(corr_ct$V1,function(s) strsplit(s,";")[[1]][1])
    corr_ct$gene2<-sapply(corr_ct$V1,function(s) strsplit(s,";")[[1]][2])
    corr_ct$V1<-NULL
    
    #Order so that gene1 is always the one first in alphabet
    corr_ct$swap<-ifelse(corr_ct$gene1 > corr_ct$gene2,corr_ct$gene1,corr_ct$gene2)
    corr_ct$gene1<-ifelse(corr_ct$gene1 > corr_ct$gene2,corr_ct$gene2,corr_ct$gene1)
    corr_ct$gene2<-corr_ct$swap
    corr_ct$swap<-NULL
    
  #Load ImmunNexUT data
  } else{
    corr_ct<-fread(paste0("imd_paper_rna_data/correlation/",sc_ct,"_correlation.txt"))
    corr_ct$V1<-NULL
    colnames(corr_ct)[3]<-"UT"
  }
  
  if(data_type=="numpy"){
    corr_bios <- np$load(paste0(path,"npy"), mmap_mode="r")
    row_names<-fread(paste0(path,rowname_suffix),header=FALSE)
    rownames(corr_bios)<-row_names$V1
    col_names<-fread(paste0(path,colname_suffix),header=FALSE)
    colnames(corr_bios)<-col_names$V1
    rm(row_names,col_names)
    
    #Filter for single cell data
    sc_genes<-sort(union(corr_ct$gene1,corr_ct$gene2))
    sc_genes<-sc_genes[sc_genes %in% colnames(corr_bios)]
    corr_bios<-corr_bios[sc_genes,sc_genes]
    corr_bios<-reshape2::melt(corr_bios)
    corr_bios$Var1<-as.character(corr_bios$Var1)
    corr_bios$Var2<-as.character(corr_bios$Var2)
    colnames(corr_bios)[1:2]<-c("gene1","gene2")
  
    
    #Order so that gene1 is always the one first in alphabet
    corr_bios$swap<-ifelse(corr_bios$gene1 > corr_bios$gene2,corr_bios$gene1,corr_bios$gene2)
    corr_bios$gene1<-ifelse(corr_bios$gene1 > corr_bios$gene2,corr_bios$gene2,corr_bios$gene1)
    corr_bios$gene2<-corr_bios$swap
    corr_bios$swap<-NULL
    corr_bios<-corr_bios[corr_bios$gene1!=corr_bios$gene2,]
  } else {
    corr_bios<-fread(path)
    corr_bios$V1<-NULL
    colnames(corr_bios)[3]<-"value"
  }
  
  #Merge both
  corr_bios<-merge(corr_bios,corr_ct,by=c("gene1","gene2"))
  
  print(paste("Overlapping genes:",length(union(corr_bios$gene1,
                                                corr_bios$gene2))))
  
  corr_corr<-cor(corr_bios$UT,corr_bios$value,
                 method="pearson")
  
  #Plot
  g<-ggplot(corr_bios,aes(UT,value))+
    geom_bin2d(bins=50)+
    xlab(xlab_text)+
    ylab(ylab_text)+
    xlim(-1,1)+ylim(-1,1)+
    #scale_fill_viridis(trans="log10")+
    scale_fill_continuous("Density",type = "viridis",
                          trans="log10",
                          breaks = c(5, 1500), 
                          labels = c("Low", "High"))+
    annotate(geom="text", x=-0.95, y=0.95,size=4,
             hjust = 0,vjust=1,
             label=paste0("r = ",format(corr_corr,digits=2)))+
    theme(axis.title=element_text(size=12),
          axis.text=element_text(size=11))
  #ggsave(g,file=plot_path)
  return(g)
}

###############################################################################
#Combine plots for Figure 2 - bulk
g_1<-create_scatter_plot(path="blueprint_data/mono_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.1PCAsOverSamplesRemoved.spearmanr.",
                         plot_path=NULL,
                         ylab_text="Correlation Blueprint",
                         sc_ct="monocyte",
                         rowname_suffix = "genes.txt",
                         colname_suffix = "genes.txt")

g_2<-create_scatter_plot(path="imd_paper_rna_data/correlation/CL_Mono_correlation.txt",
                         plot_path=NULL,
                         ylab_text="Correlation ImmuNexUT",
                         sc_ct="monocyte",
                         data_type="txt")

g_3<- create_scatter_plot(path="bios/actual_mono_bios.spearmanR.",
                          plot_path=NULL,
                          ylab_text="Correlation BIOS (measured counts)",
                          sc_ct="monocyte")

g_4<- create_scatter_plot(path="bios/deconvoluted_mono_bios.spearmanR.",
                          plot_path=NULL,
                          ylab_text="Correlation BIOS (estimated counts)",
                          sc_ct="monocyte")

g<-ggarrange(g_1,g_2,g_3,g_4,ncol=4,nrow=1,common.legend=TRUE,legend="bottom")
ggsave(g,file="bios/plots/figure2_bulk_combined.png",
       width=15,height=4)
ggsave(g,file="bios/plots/figure2_bulk_combined.pdf",
       width=15,height=4)

###############################################################################
# Combine plots for associated supplementary figure

#Plot Blueprint data set (T cells)
g<-create_scatter_plot(path="blueprint_data/tcel_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.spearmanR.",
                 plot_path="bios/plots/blueprint_tcells.png",
                 ylab_text="Correlation Blueprint - T cells",
                 sc_ct="CD4T")

plot_list<-list(g)

#Plus all ImmuNexUT cell types (except monocytes)
ct_mapping_woMono<-data.frame(sc_ct=c("CD4T","CD8T","B","NK","DC"),
                       imn_ct=c("Naive_CD4","Naive_CD8","CL_Mono","NK","mDC"))


for(i in 1:nrow(ct_mapping_woMono)){
  g<-create_scatter_plot(path=paste0("imd_paper_rna_data/correlation/",
                                       ct_mapping_woMono$imn_ct[i],
                                       "_correlation.txt"),
                           plot_path=NULL,
                           ylab_text=paste("Correlation ImmuNexUT -",
                                           ct_mapping_woMono$imn_ct[i]),
                           sc_ct=ct_mapping_woMono$sc_ct[i],
                           xlab_text=paste("Correlation single cell -",
                                           ct_mapping_woMono$sc_ct[i]),
                           data_type="txt")

  plot_list<-c(plot_list,list(g))
}

g<-ggarrange(plotlist=plot_list,ncol=3,nrow=2,common.legend=TRUE,legend="bottom")
ggsave(g,file="bios/plots/figure2_bulk_combined_supplement.png",
       width=15,height=8)
ggsave(g,file="bios/plots/figure2_bulk_combined_supplement.pdf",
       width=15,height=8)

##############################################################################
# Comparison bluprint vs ImmuNexUT

#Plot Blueprint vs ImmuNexUT (T cells)
g_1<-create_scatter_plot(path="blueprint_data/tcel_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.spearmanR.",
                 plot_path=NULL,
                 ylab_text="Correlation Blueprint - CD4T cells",
                 sc_ct="Naive_CD4",
                 scData=FALSE,
                 xlab_text="Corr ImmuNexUT - Naive CD4")

#Plot Blueprint vs ImmuNexUT (Monocytes)
g_2<-create_scatter_plot(path="blueprint_data/mono_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.1PCAsOverSamplesRemoved.spearmanr.",
                    plot_path=NULL,
                    ylab_text="Correlation Blueprint - Monocytes",
                    sc_ct="CL_Mono",
                    rowname_suffix = "genes.txt",
                    colname_suffix = "genes.txt",
                    scData=FALSE,
                    xlab_text="Corr ImmuNexUT - CL_Mono")

g_comp<-ggarrange(g_1,g_2,ncol=2,nrow=1,common.legend=TRUE,legend="bottom")
ggsave(g_comp,file="bios/plots/blueprint_immunext.png",width=6,height=4)
ggsave(g_comp,file="bios/plots/blueprint_immunext.pdf",width=6,height=4)

