if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("trackViewer")
library("trackViewer")


setwd("/Users/dilution/Documents/NDAL/prospax/analysis/SACS_analysis/filtration_reports/")

folders <- list.dirs(recursive = TRUE)
for(f in folders[2:length(folders)]){
  sample_name = strsplit(f, "/" )
  sample_name_f = sample_name[[1]]
  print(sample_name_f)
  
  
  table_name = paste(sample_name_f[2], "variants.tsv", sep=".")
  table_path = paste(f, table_name, sep="/")
  df <- read.table(file=table_path, sep = "\t", header = TRUE)
  df2<-df[!(df$protein_pos==-1),]
  if (nrow(df2)==0){
    next
  }
  else{
    SNP <- df2$protein_pos 
    
    sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=df2$hgvsp_short), )
    # color <- c("#48CA50")
    sample.gr$color <- sample.int(1, length(SNP), replace=TRUE)
    
    features <- GRanges("chr1", IRanges(c(1, 143, 145, 167, 250, 306, 508, 544, 797), 
                                        width=c(58, 1, 21, 70, 22, 175, 35, 202, 1),
                                        names=c("MTS", "FtsH", "TM1", "FtsH", "TM2", "AAA+ ATPase",
                                                "AAA+ Lid", "M41 peptidase","")))
    features$fill <- c("#80A5FC", "#DAF7A6", "#F7B520", "#94F379", "#F7B520",
                       "#900C3F",  "#C70039", "#103044", "#111111")
    features$height <- list( unit(24, "points"))
    
    lolliplot(sample.gr, features)
  }
  
}

nrow(df2)
table_lvl1 <- read.table(file="2140X/2140X.variants.tsv", sep = "\t", header = TRUE)
SNP <- table_lvl1$protein_pos
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=table_lvl1$hgvsp_short), )
# color <- c("#48CA50")
sample.gr$color <- sample.int(1, length(SNP), replace=TRUE)

features <- GRanges("chr1", IRanges(c(1, 143, 145, 167, 250, 306, 508, 544, 797), 
                                    width=c(58, 1, 21, 70, 22, 175, 35, 202, 1),
                                    names=c("MTS", "FtsH", "TM1", "FtsH", "TM2", "AAA+ ATPase",
                                            "AAA+ Lid", "M41 peptidase","")))
features$fill <- c("#80A5FC", "#DAF7A6", "#F7B520", "#94F379", "#F7B520",
                   "#900C3F",  "#C70039", "#103044", "#111111")
features$height <- list( unit(24, "points"))

lolliplot(sample.gr, features)
# SNP <- SACS_trackviewer_1$prot_pos
sample.gr$color <- df2$csq_minimal
ifelse(grepl("^rs", names(mutation.frequency)), 
       "lightcyan", "lavender")
legends <- list(labels=c("Missense", "Nonsense", "Synonymous", "Frameshift", "Inframe"),
                fill=c("#bada55", "#aabbcc", "#123661", "#ffa486", "#fcb515"),
                color=c("gray80", "gray80", "gray80", "gray80", "gray80"))


