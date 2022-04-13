# ------------------------------------------------------------------------------
# Check enrichment of TFBS among co-eGenes using Remap annotations
# Part 1) filter Remap file for blood related cell lines 
# (due to size of file this takes some time)
# ------------------------------------------------------------------------------

library(rtracklayer)

# Old version of Remap (2015)
# peaks<-import("tfbs_enrichment_remap/filPeaks_all.bed")
# 
# # Filter blood related cell lines
# ann <- t(matrix(unlist(strsplit(values(peaks)[,"name"], ".", fixed=T)), nrow=3))
# colnames(ann) <- c("geo_id", "TF", "condition")
# ann <- as.data.frame(ann,stringsAsFactors=FALSE)
# blood_related_terms<-c("amlpz12_leukemic", "aplpz74_leukemia",
#                        "bcell", "bjab", "bl41",
#                        "blood", "lcl", "erythroid", "gm",
#                        "hbp", "k562", "kasumi",
#                        "lymphoblastoid", "mm1s", "p493",
#                        "plasma", "sem", "thp1", "u937")


#Read peak file from new version (Remap2022) (different format, as all peaks file gets too large)
peaks<-import("tfbs_enrichment_remap/remap2022_nr_macs2_hg19_v1_0.bed.gz")

# Filter blood related cell lines
ann <- t(matrix(unlist(strsplit(values(peaks)[,"name"], ":", fixed=T)), nrow=2))
colnames(ann) <- c("TF", "conditions")
ann <- as.data.frame(ann,stringsAsFactors=FALSE)

conditions<-data.frame(condition=unique(unlist(strsplit(ann$conditions,","))),
                       blood_related=FALSE)
conditions<-conditions[order(conditions$condition),]

blood_related_terms<-c("ALL","AML",
                       "B-cell","BJAB","BL41","blood",
                       "CLL","DC","erythroid",
                       "erythroid-progenitor","GM",
                       "Jurkat","K-562","Kasumi",
                       "LCL","leukemia","lymphoblast","lymphocyte",
                       "macrophage","MM1-S","monocyte",
                       "neutrophil","P493","peripheral-blood",
                       "SEM","T-cell","Th1","Th17","THP-1","U-937",
                       "monocyte")

for(term in blood_related_terms){
  conditions[grep(term, conditions$condition),"blood_related"] <- TRUE
}

write.table(conditions,file="tfbs_enrichment_remap/conditions_remap2022.tsv",
            quote=FALSE,sep="\t")

conditions<-conditions[conditions$blood_related,]

ann$blood_related<-FALSE
for(term in blood_related_terms){
  ann[grep(term, ann$conditions),"blood_related"] <- TRUE
}

values(peaks)<-ann
peaks<-peaks[peaks$blood_related]
peaks$blood_related<-NULL

#Put back into name column as only this column is exported from rtracklayer
peaks$name<-paste0(peaks$TF,":",peaks$conditions)
export(peaks,"tfbs_enrichment_remap/remap2022_nr_macs2_hg19_v1_0_blood_related.bed")
