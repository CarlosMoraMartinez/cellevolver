library(tidyverse)
library(magrittr)#overwrites set_rownames and colnames
library(monocle)
library(pheatmap)

## Followed instructions from 
#https://atlas.gs.washington.edu/worm-rna/docs/
#(Cao et al., 2017)

download.file(
  "http://waterston.gs.washington.edu/sci_RNA_seq_gene_count_data/Cao_et_al_2017_vignette.RData",
  destfile = "Cao_et_al_2017_vignette.RData")
load("Cao_et_al_2017_vignette.RData")

#TFs downloaded from
#https://walhoutlab.umassmed.edu/?page_id=31
# wTF 3.0 (Fuxman et al., 2016)
tfs <- read_csv("Table-S2-wTF-3.0-Fuxman-Bass-Mol-Sys-Biol-2016_WalhoutLab.csv")
tfs2 <- tfs %>% filter(tfs$`Public name` %in% fData(cds)$symbol )
tfs_not_found <- tfs %>% filter(!tfs$`Public name` %in% fData(cds)$symbol ) #size 5

exprnames <- show.expr.info(tfs$`Public name`[1], "neuron type")$facet

#Get expression of TFs
i = 1
for(g in tfs2$`Public name`){
  tryCatch(
    {
      aux <- show.expr.info(g, "neuron type")
      j=1
      for(tis in aux$facet){
        tfs2[i, tis] <- aux$tpm[j]
        j = j + 1
      }
      i = i + 1
      cat(g, round(100*(i - 1)/nrow(tfs2), 2), "%\n")
    },
    error=cat("ERROR: ", g, "\n"))
}

tfs2 <- tfs2 %>% mutate(meanExpr = rowMeans(select(., exprnames)))
write_csv(tfs2, "200919_worm_tf_expr.csv" )


####MAKE HEATMAP WITH SCALED VARIABLES
tfs2 <- read_csv("200919_worm_tf_expr.csv")
exprnames <- names(tfs2)[5:(ncol(tfs2)-1)]
nt <- tfs2
nt <- nt %>% 
  select(exprnames)%>% 
  as.matrix %>% 
  t %>% scale %>% t %>% 
  set_colnames(exprnames) %>% 
  set_rownames(nt$`Public name`) 
nt <- nt[!(apply(nt, MAR=1, FUN=function(x)any(is.na(x)))), ]

#This is just to get clusters easily and see data aggregated by gene cluster
hm1<-pheatmap(nt, show_rownames = FALSE, angle_col = 45, treeheight_row = 50, 
         kmeans_k = ncol(nt),
         cluster_rows = T, cluster_cols = T)
row_order <- names(sort(hm1$kmeans[1][[1]]) )
col_order <- apply(hm1$kmeans$centers, MAR=2, FUN=function(x)which(x==max(x))) %>% sort %>% names
col_order <- c(col_order, colnames(nt)[! colnames(nt) %in% col_order])
nt_sorted <- nt[row_order, col_order]
#Definitive plot
pheatmap(nt_sorted, show_rownames = FALSE, angle_col = 45, treeheight_row = 50, fontsize_col = 13,
         width = 10, height = 10,
         filename = "wormTFs_scaled_ordered_size10x10.pdf",
         cluster_rows = F, cluster_cols = F)
  
