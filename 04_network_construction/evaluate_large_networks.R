################################################################################
# Analyse cell type and condition specific networks
# (generated using the python script correlation_timepoint_combined_indivs.py)
#
# updated version: use the same set of genes within each cell type 
# (union of expressed genes for different stimulations)
#
# * check degree distribution (plus betweeness and pagerank)
# * louvain clustering plus GO term enrichment
# * check if there is a connection between correlation and variance
################################################################################

library(data.table)
library(ggplot2)
library(viridis)
library(reshape2)
library(dplyr)
library(igraph) #to plot the graph
library(RColorBrewer) #to color the nodes in the graph
library(leiden)

library(clusterProfiler) #for GO enrichment
library(org.Hs.eg.db) # required for GO enrichment
library(plotrix) #only to get the color ids

theme_set(theme_bw())

setwd("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/")
#setwd("/groups/umcg-lld/tmp04/projects/1MCellRNAseq/GRN_reconstruction/ongoing/")

celltypes<-c("B","CD4T","CD8T","DC","monocyte","NK")
conditions<-c("UT","3h","24h")

convertCorrToMatrix<-function(corr_ct,cond){
  corr_ct$gene1<-gsub(";.*","",corr_ct$V1)
  corr_ct$gene2<-gsub(".*;","",corr_ct$V1)
  
  #Order so that gene1 is always the one first in alphabet
  corr_ct$swap<-ifelse(corr_ct$gene1 > corr_ct$gene2,corr_ct$gene1,corr_ct$gene2)
  corr_ct$gene1<-ifelse(corr_ct$gene1 > corr_ct$gene2,corr_ct$gene2,corr_ct$gene1)
  corr_ct$gene2<-corr_ct$swap
  corr_ct$swap<-NULL
  
  #Reformat to an edge matrix
  corr.mat<-reshape2::dcast(corr_ct,formula=gene1~gene2,value.var=cond)
  rownames(corr.mat)<-corr.mat$gene1
  corr.mat$gene1<-NULL
  
  #Add missing genes
  first.gene<-setdiff(corr_ct$gene1,corr_ct$gene2)
  last.gene<-setdiff(corr_ct$gene2,corr_ct$gene1)
  corr.mat[,first.gene]<-NA
  corr.mat[last.gene,]<-NA
  
  #Order columns in the same way as rows
  corr.mat<-as.matrix(corr.mat)
  corr.mat[corr.mat>0]<-1
  corr.mat[is.na(corr.mat)]<-0
  corr.mat<-corr.mat[sort(row.names(corr.mat)),sort(row.names(corr.mat))]
  
  #Fill missing values
  diag(corr.mat)<-0
  corr.mat[lower.tri(corr.mat)]<-t(corr.mat)[lower.tri(corr.mat)]
  
  return(corr.mat)
}

# Compare degree distribution (between cell types-conditions)
# and if hub genes have a higher variance
degree_dist<-NULL
betweenness_dist<-NULL
pagerank_dist<-NULL
for(ct in celltypes){
  
  print(ct)
  
  #var_degree_dist<-NULL
  for(cond in conditions){
    
    #Load correlation values (filtered)
    corr_ct<-fread(paste0("co-expression_indivs_combined/coexp_tp_union/",ct,"/",ct,"_",cond,
                          "_correlation_03filtered.csv"))
   
    corr.mat<-convertCorrToMatrix(corr_ct,cond)
    
    degree_current<-data.table(gene=rownames(corr.mat),
                               degree=rowSums(corr.mat),
                               timepoint=cond,
                               celltype=ct)
    
    #Get degree distribution 
    degree_dist<-rbind(degree_dist,degree_current)

    #Identify clusters and plot the graph
    graph_object <- graph_from_adjacency_matrix(corr.mat, mode = "undirected")
    partition <- leiden(corr.mat, resolution_parameter=1)

    #Save leiden partition for downstream analysis
    tmp<-data.frame(genes=rownames(corr.mat),
                    cluster=partition)

    write.table(tmp,file=paste0("co-expression_indivs_combined/coexp_tp_union/",
                                "louvain_clustering/clustering_",ct,"_",cond,
                                "_cutoff03.tsv"),
                sep="\t",quote=FALSE,row.names=FALSE)

    #Color graph nodes
    tb<-table(partition)
    #colors<-brewer.pal(sum(tb>1),"Pastel1") (only max 0 colors)
    colors<-brewer.pal(sum(tb>1),"Set3")
    single_genes_partitions<-names(tb[tb==1])
    node_cols<-ifelse(partition %in% single_genes_partitions,"#CDCDCD",
                      colors[partition])

    pdf(paste0("co-expression_indivs_combined/plots/graphs/leiden_",
               ct,"_",cond,"_standard.pdf"))
    plot(graph_object, vertex.label=NA, vertex.size=3, vertex.color = node_cols)
    dev.off()

    # At least for the monocyte cell types perform GO term enrichment
    # (Run it later on all cell types)
    if(ct=="monocyte"){

      #Check if directory exists to save the files in
      dir<-paste0("co-expression_indivs_combined/coexp_tp_union/louvain_clustering/",
                  "GO_enrichment/",ct,"_",cond)
      
      if(! dir.exists(dir)){
        dir.create(dir)
      }
      
      for(part in names(tb[tb>5])){
        print(paste("Partition:",part,"(",color.id(colors[as.numeric(part)])[1],")"))


        go_enrich<-enrichGO(gene=colnames(corr.mat)[partition==part],
                            OrgDb='org.Hs.eg.db',
                            keyType="SYMBOL",
                            pvalueCutoff = 0.05,
                            #pAdjustMethod = "none",
                            universe = colnames(corr.mat),
                            minGSSize=5)
        print(head(go_enrich[,c("ID","Description","p.adjust")]))
        
        write.table(go_enrich,file=paste0(dir,"/cluster_",part,"_color_",
                                          color.id(colors[as.numeric(part)])[1],".tsv"),
                    sep="\t",quote=FALSE,row.names=FALSE)
        
      }
    }
    
    #Get also betweenness and page rank scores
    betweenness_dist<-rbind(betweenness_dist,
                            data.table(betweenness=betweenness(graph_object),
                                       timepoint=cond,
                                       celltype=ct))
    pagerank_dist<-rbind(pagerank_dist,
                         data.table(page_rank=page_rank(graph_object)$vector,
                                    timepoint=cond,
                                    celltype=ct))
    
    # #Get variance across individuals using zscores
    # corr_zscores_indiv<-fread(paste0("co-eQTLs/input/co-expression_zscores/",
    #                                  ct,"/",ct,"_",cond,"_zscores.csv"))
    # 
    # zscores<-as.matrix(corr_zscores_indiv[,2:ncol(corr_zscores_indiv)])
    # rownames(zscores)<-corr_zscores_indiv$V1
    # vars<-apply(zscores,1,var)
    # 
    # var_indiv<-data.table(gene1=gsub(";.*","",names(vars)),
    #                       gene2=gsub(".*;","",names(vars)),
    #                       var=vars)
    # 
    # #Combine variance with degree
    # var_indiv<-merge(var_indiv,degree_current[,c("gene","degree")],
    #                  by.x="gene1",by.y="gene",all.x=TRUE)
    # var_indiv<-merge(var_indiv,degree_current[,c("gene","degree")],
    #                  by.x="gene2",by.y="gene",all.x=TRUE,
    #                  suffixes=c(".gene1",".gene2"))
    # var_indiv$degree.gene1[is.na(var_indiv$degree.gene1)]<-0
    # var_indiv$degree.gene2[is.na(var_indiv$degree.gene2)]<-0
    # var_indiv$degree.min<-pmin(var_indiv$degree.gene1,var_indiv$degree.gene2)
    # var_indiv$degree.max<-pmax(var_indiv$degree.gene1,var_indiv$degree.gene2)
    # var_indiv$degree.mean<-rowMeans(cbind(var_indiv$degree.gene1,var_indiv$degree.gene2))
    # #Get max, min and mean degree for each pair
    # 
    # #Remove NA values
    # var_indiv<-var_indiv[! is.na(var_indiv$var),]
    # 
    # #Split var_indiv into quantiles
    # 
    # quant_separate<-quantile(var_indiv$var)
    # 
    # var_indiv$group<-case_when(var_indiv$var > quant_separate[4] ~ "75%",
    #                            var_indiv$var > quant_separate[3] ~ "50%",
    #                            var_indiv$var > quant_separate[2] ~ "25%",
    #                            TRUE ~ "0%")
    # 
    # var_indiv<-melt(var_indiv[,c("group","degree.min","degree.max","degree.mean")],
    #                 id.vars=c("group"))
    # 
    # var_indiv$timepoint<-cond
    # var_degree_dist<-rbind(var_indiv,var_degree_dist)
  }
  
  # g<-ggplot(var_degree_dist,aes(x=group,y=value,color=variable))+
  #   geom_boxplot(outlier.size=0)+
  #   ylab("Degree")+xlab("Variance quantile")+
  #   facet_wrap(~timepoint)+ggtitle(ct)+
  #   theme(legend.position="bottom")
  # ggsave(g,file=paste0("co-expression_indivs_combined/plots/comp_degree_variance",
  #                      ct,".png"),
  #        width=9,height=5)
  
}

print("Genes with at least one connection per cell type and condition:")
print(table(degree_dist$timepoint,degree_dist$celltype))

#Plot degree distribution
degree_dist$timepoint<-factor(degree_dist$timepoint,levels=c("UT","3h","24h"))
g<-ggplot(degree_dist,aes(x=degree,fill=timepoint))+
  geom_histogram(position="dodge")+facet_wrap(~celltype,ncol=3,scales="free")+
  theme(legend.position = "bottom")+
  xlab("Degree")+ylab("Gene count")
ggsave(g,file="co-expression_indivs_combined/plots/degree_comparison.png",width=9)
ggsave(g,file="co-expression_indivs_combined/plots/degree_comparison.pdf",width=9)

#Plot betweenness distribution
betweenness_dist$timepoint<-factor(betweenness_dist$timepoint,levels=c("UT","3h","24h"))
g<-ggplot(betweenness_dist,aes(x=betweenness,fill=timepoint))+
  geom_histogram(position="dodge")+facet_wrap(~celltype,ncol=3,scales="free")+
  theme(legend.position = "bottom")
ggsave(g,file="co-expression_indivs_combined/plots/betweeness_comparison.png",width=9)

#Plot page rank distribution
pagerank_dist$timepoint<-factor(pagerank_dist$timepoint,levels=c("UT","3h","24h"))
g<-ggplot(pagerank_dist,aes(x=page_rank,fill=timepoint))+
  geom_histogram(position="dodge")+facet_wrap(~celltype,ncol=3,scales="free")+
  theme(legend.position = "bottom")
ggsave(g,file="co-expression_indivs_combined/plots/pagerank_comparison.png",width=9)
