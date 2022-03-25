################################################################################
# Create scatterplots for Figure 2 and the supplement (for gearshift)
# => show T cell in main figure and other cell types in the supplement
# * for comparing the single cell data set within each other (per cell type)
# * for comparing with bulk (Blueprint, ImmuNexUT, BIOS)
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
create_scatter_plot_sc<-function(sc_ct, path, xlab_text,ylab_text,plot_path,
                                 annot_text_size=4,
                                 annot_text_digits=2){
  
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
    geom_bin2d(bins=50)+geom_abline()+
    xlab(xlab_text)+
    ylab(ylab_text)+
    xlim(-1,1)+ylim(-1,1)+
    #scale_fill_viridis(trans="log10")+
    scale_fill_continuous("Density",type = "viridis",
                          trans="log10",
                          breaks = c(5, 1500), 
                          labels = c("Low", "High"))+
    annotate(geom="text", x=-0.95, y=0.95,size=annot_text_size,
             hjust = 0,vjust=1,
             label=paste0("r = ",format(corr_corr,digits=annot_text_digits)))+
    theme(axis.title=element_text(size=12),
          axis.text=element_text(size=11))
  #ggsave(g,file=plot_path)
  return(g)
}

# Cell type to be shown in main Figure
main_celltype<-"CD4T"

g_1<-create_scatter_plot_sc(sc_ct=main_celltype,
                            path=paste0("co-expression_indivs_combined/one_million_version2/",
                                        main_celltype,"/",main_celltype,
                                       "_UT_correlation.csv"),
                            xlab_text="Correlation Oelen data set (V3)",
                            ylab_text="Correlation Oelen data set (V2)",
                            plot_path=NULL)

g_2<-create_scatter_plot_sc(sc_ct=main_celltype,
                            path=paste0("co-expression_indivs_combined/stemi/version2/",
                                        main_celltype,"/",main_celltype,
                                        "_Baseline_correlation.csv"),
                           xlab_text="Correlation Oelen data set (V3)",
                           ylab_text="Correlation van Blokland data set (V2)",
                           plot_path=NULL)

g_3<-create_scatter_plot_sc(sc_ct=main_celltype,
                            path=paste0("co-expression_indivs_combined/stemi/version3/",
                                        main_celltype,"/",main_celltype,
                                        "_Baseline_correlation.csv"),
                            xlab_text="Correlation Oelen data set (V3)",
                            ylab_text="Correlation van Blokland data set (V3)",
                            plot_path=NULL)

g_4<-create_scatter_plot_sc(sc_ct=main_celltype,
                            path=paste0("co-expression_indivs_combined/ng_updated_version/",
                                        main_celltype,"/",main_celltype,
                                        "_correlation.csv"),
                            xlab_text="Correlation Oelen data set (V3)",
                            ylab_text="Correlation van der Wijst data set",
                            plot_path=NULL)

g<-ggarrange(g_1,g_2,g_3,g_4,ncol=4,nrow=1,common.legend=TRUE,legend="bottom")
ggsave(g,file="bios/plots/figure2_sc_combined.pdf",
       width=15,height=4)

#Iterate over all other cell types
plot_list<-list()
for(ct in c("CD8T","monocyte","NK","B","DC")){

  g_1<-create_scatter_plot_sc(sc_ct=ct,
                              path=paste0("co-expression_indivs_combined/one_million_version2/",
                                          ct,"/",ct,"_UT_correlation.csv"),
                              xlab_text=paste("Corr Oelen data set (V3)",ct),
                              ylab_text=paste("Corr Oelen data set (V2)",ct),
                              plot_path=NULL)

  g_2<-create_scatter_plot_sc(sc_ct=ct,
                              path=paste0("co-expression_indivs_combined/stemi/version2/",
                                          ct,"/",ct,"_Baseline_correlation.csv"),
                              xlab_text=paste("Corr Oelen data set (V3)",ct),
                              ylab_text=paste("Corr van Blokland data set (V2)",ct),
                              plot_path=NULL)

  g_3<-create_scatter_plot_sc(sc_ct=ct,
                              path=paste0("co-expression_indivs_combined/stemi/version3/",
                                          ct,"/",ct,"_Baseline_correlation.csv"),
                              xlab_text=paste("Corr Oelen data set (V3)",ct),
                              ylab_text=paste("Corr van Blokland data set (V3)",ct),
                              plot_path=NULL)

  g_4<-create_scatter_plot_sc(sc_ct=ct,
                              path=paste0("co-expression_indivs_combined/ng_updated_version/monocyte/",
                                          "monocyte_correlation.csv"),
                              xlab_text=paste("Corr Oelen data set (V3)",ct),
                              ylab_text=paste("Corr van der Wijst data set",ct),
                              plot_path=NULL)

  plot_list<-c(plot_list,list(g_1,g_2,g_3,g_4))
}

g<-ggarrange(plotlist = plot_list,ncol=4,nrow=5,common.legend=TRUE,legend="bottom")
ggsave(g,file="bios/plots/figure2_sc_combined_supplement.pdf",
       width=15,height=20)

rm(plot_list)
gc()

################################################################################
# Compare single cell and bulk data (or between different bulk data set)
create_scatter_plot<-function(path,plot_path,ylab_text,
                           sc_ct,
                           rowname_suffix="rows.txt",
                           colname_suffix="cols.txt",
                           scData=TRUE,
                           xlab_text="Correlation Oelen data set (V3)",
                           data_type="numpy",
                           annot_text_size=4,
                           annot_text_digits=2,
                           reverse=FALSE){
  
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
  
  if(! reverse){
    g<-ggplot(corr_bios,aes(UT,value))+
      geom_bin2d(bins=50)
  } else {
    g<-ggplot(corr_bios,aes(value,UT))+
      geom_bin2d(bins=50)
  }
  
  #Plot
  g<-g+geom_abline()+
    xlab(xlab_text)+
    ylab(ylab_text)+
    xlim(-1,1)+ylim(-1,1)+
    #scale_fill_viridis(trans="log10")+
    scale_fill_continuous("Density",type = "viridis",
                          trans="log10",
                          breaks = c(5, 1500), 
                          labels = c("Low", "High"))+
    annotate(geom="text", x=-0.95, y=0.95,size=annot_text_size,
             hjust = 0,vjust=1,
             label=paste0("r = ",format(corr_corr,digits=annot_text_digits)))+
    theme(axis.title=element_text(size=12),
          axis.text=element_text(size=11))
  #ggsave(g,file=plot_path)
  return(g)
}

###############################################################################
#Combine plots for Figure 2 - bulk

main_celltype<-"CD4T"

g_1<-create_scatter_plot(path="blueprint_data/tcel_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.spearmanR.",
                         plot_path=NULL,
                         ylab_text="Correlation Blueprint",
                         sc_ct=main_celltype)

g_2<-create_scatter_plot(path="imd_paper_rna_data/correlation/Naive_CD4_correlation.txt",
                         plot_path=NULL,
                         ylab_text="Correlation ImmuNexUT",
                         sc_ct=main_celltype,
                         data_type="txt")

g_3<- create_scatter_plot(path="bios/bios_correlation_tcellfiltered.tsv",
                          plot_path=NULL,
                          ylab_text="Correlation BIOS",
                          sc_ct=main_celltype,
                          data_type="txt")

g_4<-create_scatter_plot(path="blueprint_data/tcel_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.spearmanR.",
                         plot_path=NULL,
                         ylab_text="Correlation Blueprint",
                         sc_ct="Naive_CD4",
                         scData=FALSE,
                         xlab_text="Correlation ImmuNexUT")

g<-ggarrange(g_1,g_2,g_3,g_4,ncol=4,nrow=1,common.legend=TRUE,legend="bottom")
ggsave(g,file="bios/plots/figure2_bulk_combined.pdf",
       width=15,height=4)

###############################################################################
# Combine plots for associated supplementary figure

#Plot Blueprint data set (T cells)
g<-create_scatter_plot(path="blueprint_data/mono_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.1PCAsOverSamplesRemoved.spearmanr.",
                 plot_path=NULL,
                 xlab_text="Correlation Oelen data set (V3) - monocyte",
                 ylab_text="Correlation Blueprint - monocyte",
                 sc_ct="monocyte",
                 rowname_suffix = "genes.txt",
                 colname_suffix = "genes.txt")

plot_list<-list(g)

#Plus all ImmuNexUT cell types (except monocytes)
ct_mapping_woMono<-data.frame(sc_ct=c("CD8T","monocyte","B","NK","DC"),
                       imn_ct=c("Naive_CD8","CL_Mono","Naive_B","NK","mDC"))


for(i in 1:nrow(ct_mapping_woMono)){
  g<-create_scatter_plot(path=paste0("imd_paper_rna_data/correlation/",
                                       ct_mapping_woMono$imn_ct[i],
                                       "_correlation.txt"),
                           plot_path=NULL,
                           ylab_text=paste("Correlation ImmuNexUT -",
                                           ct_mapping_woMono$imn_ct[i]),
                           sc_ct=ct_mapping_woMono$sc_ct[i],
                           xlab_text=paste("Correlation Oelen data set (V3) -",
                                           ct_mapping_woMono$sc_ct[i]),
                           data_type="txt")

  plot_list<-c(plot_list,list(g))
}

g<-ggarrange(plotlist=plot_list,ncol=3,nrow=2,common.legend=TRUE,legend="bottom")
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
ggsave(g_comp,file="bios/plots/blueprint_immunext.pdf",width=6,height=4)

################################################################################
# Select example plots for main Figure 2 (a,b,c)
################################################################################

create_scatter_plot_vanBlokland_v2<-function(path,plot_path,ylab_text,
                              sc_ct,
                              rowname_suffix="rows.txt",
                              colname_suffix="cols.txt",
                              xlab_text="Correlation Oelen data set (V3)",
                              data_type="numpy",
                              annot_text_size=4,
                              annot_text_digits=2,
                              reverse=FALSE){
  
  # Load single cell data
  corr_ct<-fread(paste0("co-expression_indivs_combined/stemi/version2/",
                        sc_ct,"/",sc_ct,"_t8w_correlation.csv"))
  corr_ct$gene1<-sapply(corr_ct$V1,function(s) strsplit(s,";")[[1]][1])
  corr_ct$gene2<-sapply(corr_ct$V1,function(s) strsplit(s,";")[[1]][2])
  corr_ct$V1<-NULL
  
  #Order so that gene1 is always the one first in alphabet
  corr_ct$swap<-ifelse(corr_ct$gene1 > corr_ct$gene2,corr_ct$gene1,corr_ct$gene2)
  corr_ct$gene1<-ifelse(corr_ct$gene1 > corr_ct$gene2,corr_ct$gene2,corr_ct$gene1)
  corr_ct$gene2<-corr_ct$swap
  corr_ct$swap<-NULL
  
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
  
  corr_corr<-cor(corr_bios$t8w,corr_bios$value,
                 method="pearson")
  
  if(! reverse){
    g<-ggplot(corr_bios,aes(t8w,value))+
      geom_bin2d(bins=50)
  } else {
    g<-ggplot(corr_bios,aes(value,t8w))+
      geom_bin2d(bins=50)
  }
  
  #Plot
  g<-g+geom_abline()+
    xlab(xlab_text)+
    ylab(ylab_text)+
    xlim(-1,1)+ylim(-1,1)+
    #scale_fill_viridis(trans="log10")+
    scale_fill_continuous("Density",type = "viridis",
                          trans="log10",
                          breaks = c(5, 1500), 
                          labels = c("Low", "High"))+
    annotate(geom="text", x=-0.95, y=0.95,size=annot_text_size,
             hjust = 0,vjust=1,
             label=paste0("r = ",format(corr_corr,digits=annot_text_digits)))+
    theme(axis.title=element_text(size=12),
          axis.text=element_text(size=11))
  #ggsave(g,file=plot_path)
  return(g)
}

main_celltype<-"CD4T"

#a: STEMI v2 - Oelen v3
g<-create_scatter_plot_sc(sc_ct=main_celltype,
                          path=paste0("co-expression_indivs_combined/stemi/version2/",
                                        main_celltype,"/",main_celltype,
                                        "_t8w_correlation.csv"),
                          xlab_text="Oelen (v3)",
                          ylab_text="van Blokland (v2)",
                          plot_path=NULL,
                          annot_text_size=9,
                          annot_text_digits=3)

#Save once the legend alone
g<-g+theme(legend.position="bottom")+
  scale_fill_distiller("Density",palette="BuPu",trans="log10",
                       breaks = c(2, 600), 
                       labels = c("Low", "High"))
g_leg<-get_legend(g)
g_leg<-as_ggplot(g_leg)
ggsave(g_leg,file="bios/plots/figure2_legend_inset.pdf",width=3,height=1)

g<-g+
  ggtitle("Pairwise gene correlation")+
  theme(legend.position="none",
        plot.title=element_text(size=25),
        axis.title=element_text(size=25),
        axis.text=element_text(size=20))
ggsave(g,file="bios/plots/figure2a_exampleplot.pdf",width=5,height=5)

#b: ImmuNexUT vs van Blokland (v2)
g<-create_scatter_plot_vanBlokland_v2(path="imd_paper_rna_data/correlation/Naive_CD4_correlation.txt",
                       plot_path=NULL,
                       ylab_text="Correlation ImmuNexUT",
                       sc_ct=main_celltype,
                       data_type="txt",
                       annot_text_size=9,
                       annot_text_digits=3,
                       reverse=TRUE)
g<-g+
  xlab("ImmuNexUT")+
  ylab("van Blokland (v2)")+ 
  ggtitle("Pairwise gene correlation")+
  scale_fill_distiller(palette="BuPu",trans="log10")+
  theme(legend.position="none",
        plot.title=element_text(size=25),
        axis.title=element_text(size=25),
        axis.text=element_text(size=20))
ggsave(g,file="bios/plots/figure2b_exampleplot.pdf",width=5,height=5)

#c: ImmunNexuT - Blueprint
g<-create_scatter_plot(path="blueprint_data/tcel_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.spearmanR.",
                       xlab_text="ImmuNexUT",
                       ylab_text="BLUEPRINT",
                       plot_path=NULL,
                       sc_ct="Naive_CD4",
                       scData=FALSE,
                       annot_text_size=9,
                       annot_text_digits=3)

g<-g+
  ggtitle("Pairwise gene correlation")+
  scale_fill_distiller(palette="BuPu",trans="log10")+
  theme(legend.position="none",
        plot.title=element_text(size=25),
        axis.title=element_text(size=25),
        axis.text=element_text(size=20))
ggsave(g,file="bios/plots/figure2c_exampleplot.pdf",width=5,height=5)