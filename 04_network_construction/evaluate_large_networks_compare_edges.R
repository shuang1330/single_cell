################################################################################
# Compare number of differing edges
#
# Calculate the sensitivity (how many of the edges from network X are also found
# in network Y)
################################################################################

library(data.table)
library(ggplot2)
library(viridis)

theme_set(theme_bw())

setwd("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/")
#setwd("/groups/umcg-lld/tmp04/projects/1MCellRNAseq/GRN_reconstruction/ongoing/")

celltypes<-c("B","CD4T","CD8T","DC","monocyte","NK")
conditions<-c("UT","3h","24h")


#Load the file with the intersect of genes (to prevent iterating multiple times)
intersect_set<-read.csv(paste0("co-expression_indivs_combined/coexp_tp_union/genelists_tp_union/",
                               "expressed_gene_intersect_cts.tsv"))

#all_tested_pairs<-nrow(intersect_set)*(nrow(intersect_set)-1)/2

#Do all pairwise comparisons
common_edges<-NULL
for(i in 1:17){
  
  #Load correlation values first data set
  ct_1<-celltypes[ceiling(i/3)]
  cond_1<-conditions[((i-1) %% 3)+1]
  corr_ct_1<-fread(paste0("co-expression_indivs_combined/coexp_tp_union/",ct_1,"/",ct_1,"_",cond_1,
                          "_correlation_03filtered.csv"))
  colnames(corr_ct_1)<-c("pair","corr_1")
  
  #Filter for the intersect gene set
  corr_ct_1$gene1<-gsub(";.*","",corr_ct_1$pair)
  corr_ct_1$gene2<-gsub(".*;","",corr_ct_1$pair)
  corr_ct_1<-corr_ct_1[corr_ct_1$gene1 %in% intersect_set$genes &
                         corr_ct_1$gene2 %in% intersect_set$genes,]
  
  common_edges<-rbind(common_edges,
                      data.frame(ct_1,cond_1,ct_2=ct_1,cond_2=cond_1,
                                 total_edges=nrow(corr_ct_1),
                                 edge_sensitivity=1))
  for(j in (i+1):18){
    
    #Load correlation values first data set
    ct_2<-celltypes[ceiling(j/3)]
    cond_2<-conditions[((j-1) %% 3)+1]
    corr_ct_2<-fread(paste0("co-expression_indivs_combined/coexp_tp_union/",ct_2,"/",ct_2,"_",cond_2,
                            "_correlation_03filtered.csv"))
    colnames(corr_ct_2)<-c("pair","corr_2")
    
    #Filter for the intersect gene set
    corr_ct_2$gene1<-gsub(";.*","",corr_ct_2$pair)
    corr_ct_2$gene2<-gsub(".*;","",corr_ct_2$pair)
    corr_ct_2<-corr_ct_2[corr_ct_2$gene1 %in% intersect_set$genes &
                           corr_ct_2$gene2 %in% intersect_set$genes,]
    
    #Merge both
    corr_ct<-merge(corr_ct_1,corr_ct_2,by="pair")
    common_edges<-rbind(common_edges,
                        data.frame(ct_1,cond_1,ct_2,cond_2,
                                   total_edges=nrow(corr_ct),
                                   edge_sensitivity=nrow(corr_ct)/nrow(corr_ct_1)),
                        data.frame(ct_1=ct_2,cond_1=cond_2,ct_2=ct_1,cond_2=cond_1,
                                   total_edges=nrow(corr_ct),
                                   edge_sensitivity=nrow(corr_ct)/nrow(corr_ct_2)))
    
  }
}

#Add for NK cells the diagonal entry
common_edges<-rbind(common_edges,
                    data.frame(ct_1=ct_2,cond_1=cond_2,ct_2,cond_2,
                               total_edges=nrow(corr_ct_2),
                               edge_sensitivity=1))

common_edges$cond_1<-factor(common_edges$cond_1,levels=c("24h","3h","UT"))
common_edges$cond_2<-factor(common_edges$cond_2,levels=c("UT","3h","24h"))
g<-ggplot(common_edges,aes(x=cond_1,y=cond_2,fill=edge_sensitivity))+
  geom_tile()+facet_grid(ct_1~ct_2)+
  xlab("Timepoint data set 1")+ylab("Timepoint data set 2")+
  scale_fill_viridis("Sensitvity",limits=c(0,1))
ggsave(g,file="co-expression_indivs_combined/plots/compare_edges_fraction.png",width=9)
ggsave(g,file="co-expression_indivs_combined/plots/compare_edges_fraction.pdf",width=9)

#Not giving nice pictures
# g<-ggplot(common_edges,aes(x=cond_1,y=cond_2,fill=log(total_edges)))+
#   geom_tile()+facet_grid(ct_1~ct_2)+
#   xlab("Timepoint data set 1")+ylab("Timepoint data set 2")+
#   scale_fill_viridis("# common edges")
# ggsave(g,file="co-expression_indivs_combined/plots/compare_edges_abs.png",width=9)


