##################################################################################
# Benchmark our correlation results with the perturbation data set of Monique
# * focus on UT cells and compare with CD4T cells (closest cell type)
# * new version: Wilcoxon Test
#
# Overlap between single cell and bulk data not perfect 
# (see script investigate_overlap_problem.R). Main reason is that chrMT and chrX
# is included for the single cell data (but not in the bulk data).
# Decision for now: take all genes for the single cell data set and do NOT remove
# the ones that are not found in the bulk data set.
#################################################################################

library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(gtools)
library(dplyr)

theme_set(theme_bw())

setwd("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing")

#Get colors
cols_brewer <- c(brewer.pal(n = 3, "Set2")[2],"grey78")

#Set MT correction for KO gene identification
MTcorrection<-"FDR" #alternatives: "Bonf","FDR"
print(paste("MT correction:",MTcorrection))

# Load single cell data
ct<-"CD4T" #alternative "CD8T"
cond<-"UT"

corr_sc<-fread(paste0("co-expression_indivs_combined/",ct,"/",ct,"_",cond,
                      "_correlation.csv"))
corr_sc$gene1<-gsub(";.*","",corr_sc$V1)
corr_sc$gene2<-gsub(".*;","",corr_sc$V1)

corr_sc$swap<-ifelse(corr_sc$gene1 > corr_sc$gene2,corr_sc$gene1,corr_sc$gene2)
corr_sc$gene1<-ifelse(corr_sc$gene1 > corr_sc$gene2,corr_sc$gene2,corr_sc$gene1)
corr_sc$gene2<-corr_sc$swap
corr_sc$swap<-NULL

# Load ImmuNexUT data (already preprocessed correctly)
ct<-"Naive_CD4"
corr_imn<-fread(paste0("imd_paper_rna_data/correlation/",
                      ct,"_correlation_extended.txt"))
corr_imn$V1<-NULL
colnames(corr_imn)[3]<-"UT"

# Filter for genes that are expressed in both data sets
expressed_genes_sc<-union(corr_sc$gene1,corr_sc$gene2)
expressed_genes_bulk<-union(corr_imn$gene1,corr_imn$gene2)
expressed_genes<-intersect(expressed_genes_sc,expressed_genes_bulk)

print(paste("Number of genes expressed in both data sets:",length(expressed_genes)))

# Load perturbation data
path<-"perturbation_dataset/perturbation_data/CD4T_GATE2019_MAST_DE/WT_KO/"
path_negControl <- "perturbation_dataset/perturbation_data/CD4T_GATE2019_MAST_DE/WT_NP/"

# Get a list with all DE genes
files<-list.files(path)

# Use the setting without artifical cells
files<-files[!startsWith(files,"artificialCells_")]
genes<-unique(sapply(files,function(fl) strsplit(fl,"\\.")[[1]][1]))
print(paste0("Unique genes:",length(genes)))

genes<-genes[genes %in% expressed_genes]
print(paste0("Unique genes expressed in 50% of cells:",length(genes)))

# Iterate over both data sets and all KO genes to perform Wilcoxon test
p_vals<-NULL
all_comps<-NULL
for(data_type in c("ImmuNexUT","sc")){
  
  if(data_type == "sc"){
    corr_ct<-corr_sc
    expressed_genes<-expressed_genes_sc
   
  } else if (data_type == "ImmuNexUT"){
    corr_ct<-corr_imn
    expressed_genes<-expressed_genes_bulk
  }
  
  #Bonferroni cutoff corrected for the number of expressed genes
  if(MTcorrection == "Bonf"){
    cutoff<-0.05/length(expressed_genes)
  } else {
    cutoff<-0.05
  }
  
  # Go over each gene
  plot_list<-list()
  for(gene in genes){
    
    corr_ct_ko<-corr_ct[corr_ct$gene1==gene,c("gene2","UT")]
    colnames(corr_ct_ko)[1]<-"gene1"
    corr_ct_ko<-rbind(corr_ct_ko,corr_ct[corr_ct$gene2==gene,c("gene1","UT")])
    
    #Use absolute correlation
    corr_ct_ko$UT<-abs(corr_ct_ko$UT)
    
    #Get all knock_out genes
    all_measured_ko_genes<-NULL
    ko_genes_combined<-NULL
    for(fl in files[startsWith(files,gene)]){
      ko_genes<-read.table(paste0(path,fl))
      all_measured_ko_genes<-union(all_measured_ko_genes,rownames(ko_genes))
      
      if(MTcorrection=="FDR"){
        ko_genes<-ko_genes[rownames(ko_genes) %in% expressed_genes,]
        ko_genes$p_val<-p.adjust(ko_genes$p_val,method="BH")
      }
      
      #Filter for expressed genes and significant threshold
      ko_genes<-ko_genes[rownames(ko_genes) %in% expressed_genes &
                           ko_genes$p_val<cutoff,]
      ko_genes_combined<-union(ko_genes_combined,rownames(ko_genes))
    }
    
    print(gene)
    # print(paste("Measured genes:",length(all_measured_ko_genes)))
    # print(paste("Number associated genes:",length(ko_genes_combined)))
    # 
    genes_found_in_both<-intersect(all_measured_ko_genes,expressed_genes)
    print(paste("Number genes measured in both:",length(genes_found_in_both)))
    
    #Remove false positive DE genes
    for(fl in list.files(path_negControl,pattern=gene)){
      fp_genes<-read.table(paste0(path_negControl,fl))
      
      if(MTcorrection=="FDR"){
        fp_genes<-ko_genes[rownames(fp_genes) %in% expressed_genes,]
        fp_genes$p_val<-p.adjust(fp_genes$p_val,method="BH")
      }
      
      fp_genes<-fp_genes[fp_genes$p_val<cutoff,]
      ko_genes_combined<-setdiff(ko_genes_combined,rownames(fp_genes))
    }
    
    corr_ct_ko<-corr_ct_ko[corr_ct_ko$gene1 %in% genes_found_in_both,]
    corr_ct_ko$is_ko<-corr_ct_ko$gene1 %in% ko_genes_combined
    
    corr_ct_ko$ko_gene<-paste0(gene," (",sum(corr_ct_ko$is_ko)," DE genes)")
    corr_ct_ko$data_type<-data_type
    
    all_comps<-rbind(all_comps,corr_ct_ko)
    
    wt<-wilcox.test(corr_ct_ko$UT[corr_ct_ko$is_ko],corr_ct_ko$UT[!corr_ct_ko$is_ko],
                    paired=FALSE,alternative="greater")
    
    p_vals<-rbind(p_vals,data.frame(ko_gene=paste0(gene," (",sum(corr_ct_ko$is_ko)," DE genes)"),
                                    data_type,pval=wt$p.value))
  }
}

all_comps$data_type<-ifelse(all_comps$data_type=="sc","single cell",
                            "ImmuNexUT")
all_comps$data_type<-factor(all_comps$data_type,
                            levels=c("single cell","ImmuNexUT"))
all_comps$is_ko<-ifelse(all_comps$is_ko,"DE genes","Non DE\ngenes")
all_comps$is_ko<-factor(all_comps$is_ko,levels=c("DE genes","Non DE\ngenes"))

p_vals$data_type<-ifelse(p_vals$data_type=="sc","single cell",
                            "ImmuNexUT")
# p_vals$text<-stars.pval(p_vals$pval)
# p_vals$text[p_vals$text=="."]<-"+"

p_vals$pvaltext<-paste0("p = ",format(round(p_vals$pval,3),nsmall=3))
p_vals$pvaltext<-ifelse(p_vals$pvaltext=="p = 0.000","p < 0.001",p_vals$pvaltext)

#Set text always to maximal value
max_corr<-all_comps%>%
  group_by(ko_gene,data_type)%>%
  summarize(max_UT=max(UT))

p_vals<-merge(p_vals,max_corr,by=c("ko_gene","data_type"))
p_vals$max_UT<-p_vals$max_UT*1.1
p_vals$is_ko<-1.5
#p_vals$is_ko<-factor(p_vals$is_ko,levels=c("TRUE","FALSE"))

g_sc<-ggplot()+
  geom_violin(data=all_comps[all_comps$data_type=="single cell",],
               aes(x=is_ko,y=UT,fill=is_ko))+
  geom_boxplot(data=all_comps[all_comps$data_type=="single cell",],
               aes(x=is_ko,y=UT,fill=is_ko),
               width = 0.15, outlier.shape = NA)+
  geom_text(data=p_vals[p_vals$data_type=="single cell",],
            aes(x=is_ko,y=max_UT,label=pvaltext),
            size=6)+
  facet_wrap(~ko_gene,scales ="free",nrow=1)+
  xlab("")+
  ylab("Absolute correlation (single cell)")+
  scale_fill_manual("DE gene\nafter KO",values=cols_brewer)+
  theme(legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        strip.text=element_text(size=15))

g_bulk<-ggplot()+
  geom_violin(data=all_comps[all_comps$data_type=="ImmuNexUT",],
               aes(x=is_ko,y=UT,fill=is_ko))+
  geom_boxplot(data=all_comps[all_comps$data_type=="ImmuNexUT",],
               aes(x=is_ko,y=UT,fill=is_ko),
               width = 0.15, outlier.shape = NA)+
  geom_text(data=p_vals[p_vals$data_type=="ImmuNexUT",],
            aes(x=is_ko,y=max_UT,label=pvaltext),
            size=6)+
  facet_wrap(~ko_gene,scales ="free",nrow=1)+
  xlab("")+
  ylab("Absolute correlation (ImmuNexUT)")+
  scale_fill_manual("DE gene\nafter KO",values=cols_brewer)+
  theme(legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        strip.text=element_text(size=15))

g<-ggarrange(g_sc,g_bulk,ncol=1,align="hv")

ggsave(g,file="perturbation_dataset/plots/wilcoxon_all_combined.pdf",
       width=15,height=8)
                          
# #Save merged plot
# g_comp<-ggarrange(plotlist=plot_list,ncol=5,nrow=1)
# 
# if(data_type == "sc"){
#   ggsave(g_comp,file=paste0("perturbation_dataset/plots/wilcoxon_allgenes_",
#                             MTcorrection,".pdf"),
#          width=15,height=4)
#   ggsave(g_comp,file=paste0("perturbation_dataset/plots/wilcoxon_allgenes_",
#                             MTcorrection,".png"),
#                               width=15,height=4)
# } else if (data_type == "ImmuNexUT"){
#   ggsave(g_comp,file=paste0("perturbation_dataset/plots/wilcoxon_ImmuNexUT_allgenes_",
#                             MTcorrection,".pdf"),
#          width=15,height=4)
#   ggsave(g_comp,file=paste0("perturbation_dataset/plots/wilcoxon_ImmuNexUT_allgenes_",
#                             MTcorrection,".png"),
#          width=15,height=4)
# }

# #Check how similar the knockout of different guide RNAs is
# #gene<-"FOXP1"
# gene<-"PHB2"
# test_files<-files[startsWith(files,gene)]
# res_1<-read.table(paste0(path,test_files[1]))
# res_2<-read.table(paste0(path,test_files[2]))
# res_1$gene<-rownames(res_1)
# res_2$gene<-rownames(res_2)
# res<-merge(res_1,res_2,by="gene")
# g<-ggplot(res,aes(-log10(p_val.x),-log10(p_val.y)))+
#   geom_point()+ggtitle(gene)
# ggsave(g,file=paste0("perturbation_dataset/plots/compare_guideRNA_",gene,"_pval.png"))
# 
# g<-ggplot(res,aes(avg_logFC.x,avg_logFC.y))+
#   geom_point()+ggtitle(paste0(gene," - corr=",round(cor(res$avg_logFC.x,res$avg_logFC.y),2)))
# ggsave(g,file=paste0("perturbation_dataset/plots/compare_guideRNA_",gene,"_fc.png"))
