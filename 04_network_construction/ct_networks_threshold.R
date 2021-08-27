# ------------------------------------------------------------------------------
# Evaluate for both Blueprint and ImmuNexUT for each cell type independent
# the best matching cutoffs (where MCC is the highest)
#
# Check also how many gene pairs would be selected each time
# ------------------------------------------------------------------------------

library(data.table)
#library(mccr) #to calculate MCC (not working if data set gets too large)
library(yardstick) # to calculate MCC also for larger data sets
library(viridis)
library(ggplot2)
library(ggpubr)

theme_set(theme_bw())

setwd("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing")

estimate_cutoff_mcc<-function(corr,cutoff=seq(0.1,0.9,0.05)){
  
  res<-NULL
  for(c_sc in cutoff){
    
    for(c_bl in cutoff){
      
      res<-rbind(res,
                 data.frame(c_bl,c_sc,
                            #Calculate Matthew's correlation coefficient
                            mcc=mcc_vec(factor(corr$corr_bp>c_bl,
                                               levels=c("FALSE","TRUE")),
                                        factor(corr$corr_sc>c_sc,
                                              levels=c("FALSE","TRUE"))),
                            nedges_sc=sum(corr$corr_sc>c_sc),
                            nedges_bp=sum(corr$corr_bp>c_bl)))
    }
  }
  
  #Set NA values to 0 (happens if there are no single cell levels above the threshold)
  #is.na(res$mcc) == (res$nedges_sc==0)
  res$mcc[is.na(res$mcc)]<-0

  return(res)
}

###############################################################################
# Test for ImmuNexUT cell types
###############################################################################

# Cell type matching
ct_mapping<-data.frame(sc_ct=c("CD4T","CD8T","B","monocyte","NK","DC"),
           imn_ct=c("Naive_CD4","Naive_CD8","Naive_B","CL_Mono","NK","mDC"))

plot_list<-list()
for(i in 1:nrow(ct_mapping)){
  
  print(ct_mapping$imn_ct[i])
  
  #Read single cell correlation
  corr_sc<-fread(paste0("co-expression_indivs_combined/",ct_mapping$sc_ct[i],"/",
                        ct_mapping$sc_ct[i],"_UT_correlation.csv"))
  corr_sc$gene1<-gsub(";.*","",corr_sc$V1)
  corr_sc$gene2<-gsub(".*;","",corr_sc$V1)
  
  #Order so that gene1 is always the one first in alphabet
  corr_sc$swap<-ifelse(corr_sc$gene1 > corr_sc$gene2,corr_sc$gene1,corr_sc$gene2)
  corr_sc$gene1<-ifelse(corr_sc$gene1 > corr_sc$gene2,corr_sc$gene2,corr_sc$gene1)
  corr_sc$gene2<-corr_sc$swap
  corr_sc$swap<-NULL
  corr_sc$V1<-NULL
  
  #Read ImmuNexUT data (already preprocessed correctly)
  corr_bulk<-fread(paste0("imd_paper_rna_data/correlation/",
                          ct_mapping$imn_ct[i],"_correlation.txt"))
  corr_bulk$V1<-NULL
  # all(corr_bulk$gene1<corr_bulk$gene2)
  corr_bulk<-merge(corr_bulk,corr_sc,by=c("gene1","gene2"))
  colnames(corr_bulk)[3:4]<-c("corr_bp","corr_sc")
  
  #Get optimization result
  optim_res<-estimate_cutoff_mcc(corr_bulk)
  print(optim_res[which.max(optim_res$mcc),])
  
  #Plot MCC results
  g<-ggplot(optim_res,aes(x=c_bl,y=c_sc,fill=mcc))+
    geom_tile()+
    xlab(paste("Correlation cutoff ImmuNexUT",
               ct_mapping$imn_ct[i]))+
    ylab(paste("Correlation cutoff single cell",
               ct_mapping$sc_ct[i]))+
    scale_fill_viridis(limits = c(-0.1,0.6))
    #scale_fill_gradient2(limits = c(-1,1),low="darkblue",mid="white",high="darkred")
  # ggsave(g,file=paste0("imd_paper_rna_data/plots/mcc_optim_",
  #                      ct_mapping$imn_ct[i],".png"),
  #        width=6,height=5)
  plot_list<-c(plot_list,list(g))
}

g<-ggarrange(plotlist=plot_list,ncol=3,nrow=2,common.legend=TRUE,legend="bottom")
ggsave(g,file="imd_paper_rna_data/plots/mcc_optim_allcts.pdf",
       width=12,height=10)

