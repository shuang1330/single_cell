# ------------------------------------------------------------------------------
# Barplot Blueprint cutoff dependency
# (values read manually from result plots)
# ------------------------------------------------------------------------------

library(ggplot2)
library(viridis)
theme_set(theme_bw())

vals<-read.table("co-expression_indivs_combined/blueprint_cutoff_eval.txt",
                 sep=",",header=TRUE)

vals$threshold<-as.factor(vals$threshold)
g<-ggplot(vals,aes(x=threshold,y=corr_pearson,fill=ngenes))+
  geom_bar(stat="identity")+xlab("Expression cutoff")+
  ylab("Correlation between data sets")+ylim(0,1)+
  scale_fill_viridis("Number selected\ngenes")+
  theme(legend.position = "right")
#print(g)
ggsave(g,file="co-expression_indivs_combined/plots/eval_blueprint_cutoff.png",
       width=5,height=3)
ggsave(g,file="co-expression_indivs_combined/plots/eval_blueprint_cutoff.pdf",
       width=5,height=3)
