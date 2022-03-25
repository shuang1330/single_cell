# ------------------------------------------------------------------------------
# Barplot Blueprint cutoff dependency
# (values read manually from result plots)
# ------------------------------------------------------------------------------

library(ggplot2)
#library(viridis)
library(RColorBrewer)
theme_set(theme_bw())

# vals<-read.table("co-expression_indivs_combined/blueprint_cutoff_eval_mono.txt",
#                  sep=",",header=TRUE)
vals<-read.table("co-expression_indivs_combined/blueprint_cutoff_eval_CD4T.txt",
                 sep=",",header=TRUE)

vals$threshold<-as.factor(vals$threshold)
g<-ggplot(vals,aes(x=threshold,y=corr_pearson,fill=ngenes))+
  geom_bar(stat="identity")+
  geom_text(aes(x = threshold, y = corr_pearson / 2, label = ngenes,
            color=ifelse(ngenes<1000,'white','black')),size=5)+
  scale_color_manual(values=c("black","white"))+
  # geom_text(aes(x = threshold, y = corr_pearson + 0.12, label = 
  #                 paste0("(",ngenes,")")),
  #           size=4,angle = 90)+
  xlab("Expression cutoff")+
  ylab("Correlation between Oelen (v3)\nand BLUEPRINT")+ylim(0,1)+
  scale_fill_distiller("Number of\ngenes",palette="YlOrBr")+
  #scale_fill_viridis("Number selected\ngenes",option="mako")+
  theme(legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        legend.title=element_text(size=13),
        legend.text=element_text(size=12))+
  guides(color=FALSE)
print(g)
# ggsave(g,file="co-expression_indivs_combined/plots/eval_blueprint_cutoff.png",
#        width=5,height=3)
ggsave(g,file="co-expression_indivs_combined/plots/eval_blueprint_cutoff.pdf",
       width=6,height=5)
