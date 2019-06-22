library(magrittr)
library(ggplot2)
library(gridExtra)
library(scales)# to use personalized scales in ggplot2
library(reshape)

# Script to get motif positions with the most biased gene types


#### script to plot motifs


#map types to numbers and colors
type1define = list('3'='lineage this cell','4'='lineage other cell', '8'='other activators', 
                   '9'='inhibitors', '10'='inhibitor lineage this cell',
                   '11'='inhibitor lineage other cell',
                   '5'='lineage many this', '6'='lineage many other','7'='lineage all', 
                   '12'='terminal specific this', '15'='terminal specific other',
                   '13'='terminal 2 this', '16'='terminal 2 other', 
                   '14'='terminal all', '1'='tf', '2'='non tf', '0'='lin')
col2map <- data.frame(t=unlist(type1define), col = c("limegreen", "purple3", "dodgerblue1", "firebrick3", "red", "red", "chartreuse3" ,"chartreuse3", "orange", "burlywood4", "gray", "cyan2", "darkolivegreen3", "hotpink", "blue", "yellow", "green"), stringsAsFactors = F)
kols=c()
for(i in 1:nrow(col2map))kols[col2map$t[i]]=col2map$col[i]


## relevant filetypes:

threshold <- 100
min_z_score<-2
minprop <- 0.3
maxprop <- 1-minprop

motsum <- "motifSummary"
motprop <- "motifsByType2_prop_motifs"
motcount <- "motifsByType2_count"
propfname <- "proportion_differences_ordered"
#### script to plot motifs


condlist = c('4cell_mce0', '4cell_mce1fix', '4cell_mce2fix', '4cell_mce2inh5', '4cell_mce0fix')
setwd("C:/Users/Carlos/Desktop/cellEvolver_all/1802_cellevolver/all_finaldata/motifs_sep1_sepmot3_selby0_filtTrue_ignoreSelfReg_fixed1/")

files <- list.files(recursive = FALSE, full.names=TRUE)
files <- files[grep(propfname, files)]
files2 <- list.files(recursive = TRUE, full.names=TRUE)
files2 <- files2[grepl(motsum, files2) |grepl(motprop, files2) | grepl(motcount, files2)]

for (c in condlist){
  this_files <- files[grep(paste(c,  propfname, sep="__", collapse="__"), files)]
  this_files2 <- files2[grep(paste('\\/', c, '\\/', sep="", collapse=""), files2)]
  #this_files2 <- this_files2[grep(motsum, this_files2)]
  summ <- preprocessSum(c, this_files2)
  summ <- summ[summ[, c]>= threshold & summ$z_score >= min_z_score, ]

  tt <- read.table(this_files, header=TRUE, stringsAsFactors = F)
  tt$motif <- sapply(tt$variable, FUN=function(x)strsplit(x, '_')[[1]][1])
  tt <- tt[(tt$terminal_all >= maxprop) | (tt$terminal_specific_this >= maxprop), ]
  
  tt <- tt[tt$motif %in% summ$motif_id, ]
  tt <- tt[, c("variable", "inhibitor", "lineage_other_cell", "lineage_this_cell", "other_activator", "terminal_all", "terminal_specific_other", "terminal_specific_this", "z_score")]
  names(tt)[2:ncol(tt)] <- paste(c, names(tt)[2:ncol(tt)], sep=";")
  if(c == condlist[1]){
    bigt <- tt
  }else{
    bigt <- bigt[bigt$variable %in% tt$variable, ]
    bigt <- merge(bigt, tt, by.x="variable", by.y = "variable")
  }
}
#Now calculate mean and variance across conditions
bigt$motif <- sapply(bigt$variable, FUN=function(x)strsplit(x, '_')[[1]][1])

vars <- c("inhibitor", "lineage_other_cell", "lineage_this_cell", "other_activator", "terminal_all", "terminal_specific_other", "terminal_specific_this", "z_score")
for (vv in vars){
  xx <- names(bigt)[grep(vv, names(bigt))]
  m <- as.matrix(bigt[,xx]) 
  means <- apply(m, MAR=1, mean)
  stds <- apply(m, MAR=1, sd)
  bigt[, paste(vv, "condMean", sep="_", collapse="_")] <- means
  bigt[, paste(vv, "condStd", sep="_", collapse="_")] <- stds
}
bigt <- bigt[order(bigt$z_score_condMean, decreasing = TRUE), ]

write.table(bigt, file="180918_mostBiasedMotifsAcrossConditions.csv", sep="\t", row.names = FALSE, quote=FALSE)
