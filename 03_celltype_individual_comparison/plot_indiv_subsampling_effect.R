# ------------------------------------------------------------------------------
# Plot effect of number of cells on correlation between individuals
# Correlation calculation done in python script (correlation_subsampling.py)
# ------------------------------------------------------------------------------

library(ggplot2)

theme_set(theme_bw())

suffix<-"v3"
#suffix<-"v2"

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#Full cell type names as reported in the paper
cell_types_corrected<-setNames(c("CD4+ T","CD8+ T","Monocyte","NK","DC","B"),
                               c("CD4T","CD8T","monocyte","NK","DC","B"))

#Load results
res<-read.csv(paste0("co-expression_indivs_subsampled/",
                     "correlation_individuals_subsampled_1M_",suffix,".csv"))
res$X<-NULL

res$celltype<-cell_types_corrected[res$celltype]

#For v2 B cells are missing
if(length(unique(res$celltype))==5){
  ct_colors<-gg_color_hue(6)[2:6]
} else {
  ct_colors<-gg_color_hue(6)
}

#Filter out some values to make it more visible
#res<-res[res$cell_num<450,]
res<-res[res$cell_num %in% seq(25,500,50),]

res$cell_num<-as.factor(res$cell_num)

g<-ggplot(res,aes(x=cell_num,y=corr,fill=celltype))+
  geom_violin(position = position_dodge(0.9)) +
  # geom_boxplot(width = 0.15, position = position_dodge(0.9),
  #              outlier.shape = NA)+
  xlab("Subsampled number of cells per individual")+
  ylab("Correlation between individuals")+
  ylim(0,1)+
  scale_fill_manual("Cell type",values=ct_colors)+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        legend.title = element_text(size=13),
        legend.text = element_text(size=13),
        legend.position=c(0.9,0.2))+
  guides(fill=guide_legend(nrow=3,byrow=FALSE))
print(g)
ggsave(g,file=paste0("co-expression_indivs_subsampled/plots/subsampling_1M_",
                     suffix,"_filtered.pdf"),
       width=14,height=4)
