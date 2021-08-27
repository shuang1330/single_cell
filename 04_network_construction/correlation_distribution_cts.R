#-------------------------------------------------------------------------------
# Check overall correlation distribution for each cell type and time point
#-------------------------------------------------------------------------------

library(data.table)
library(ggplot2)

theme_set(theme_bw())

setwd("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/")

cell_types<-c("B","CD4T","CD8T","DC","monocyte","NK")
timepoints<-c("UT","3h","24h")

corr_boxplot<-NULL
all_corrs<-NULL
for(ct in cell_types){
  
  union_expressed_genes<-NULL
  for(tp in timepoints){
    
    #Load correlation values
    corr_ct<-fread(paste0("co-expression_indivs_combined/coexp_tp_union/",ct,"/",
                          ct,"_",tp,"_correlation.csv"))
    colnames(corr_ct)[2]<-"corr"
    
    #Save correlation (with direction) for boxplots
    corr_ct$label<-paste0(ct,"_",tp)
    corr_boxplot<-rbind(corr_boxplot,corr_ct)
    
    #Get absolute correlation
    corr_ct$corr<-abs(corr_ct$corr)
    
    all_corrs<-rbind(all_corrs,
                     data.frame(level=c(">0.05",">0.1",">0.2",">0.3"),
                                values=c(sum(abs(corr_ct$corr)>0.05 &
                                               abs(corr_ct$corr)<0.1),
                                         sum(abs(corr_ct$corr)>0.1 &
                                               abs(corr_ct$corr)<0.2),
                                         sum(abs(corr_ct$corr)>0.2 & 
                                               abs(corr_ct$corr)<0.3),
                                         sum(abs(corr_ct$corr)>0.3)),
                                freq=c(mean(abs(corr_ct$corr)>0.05 &
                                             abs(corr_ct$corr)<0.1),
                                       mean(abs(corr_ct$corr)>0.1 &
                                             abs(corr_ct$corr)<0.2),
                                       mean(abs(corr_ct$corr)>0.2 & 
                                             abs(corr_ct$corr)<0.3),
                                       mean(abs(corr_ct$corr)>0.3)),
                                ct,tp))
    
  }
}

g<-ggplot(corr_boxplot,aes(x=label,y=corr))+
  geom_boxplot()+xlab("Network (cell type - timepoint)")+
  ylab("Correlation")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(g,file="co-expression_indivs_combined/plots/correlation_boxplot.png",width=6,height=4)

all_corrs$level<-factor(all_corrs$level,levels=c(">0.3",">0.2",">0.1",">0.05"))
all_corrs$tp<-factor(all_corrs$tp,levels=c("UT","3h","24h"))
g<-ggplot(all_corrs,aes(x=tp,y=freq,fill=level))+
  geom_bar(stat="identity")+
  xlab("Network (cell type - timepoint)")+ylab("Fraction correlated genes")+
  scale_fill_discrete("Absolute\ncorrelation")+
  facet_wrap(~ct,nrow=1)
ggsave(g,file="co-expression_indivs_combined/plots/correlation_level_freq.png",width=6,height=4)
ggsave(g,file="co-expression_indivs_combined/plots/correlation_level_freq.pdf",width=6,height=4)

