library(magrittr)
library(ggplot2)


makeMotifPlots <- function(zonly){
  
  colInd<-makeHeatmap(zonly)
  library(igraph)
  library(gridExtra)
  ## make plots for each motif
  pdf("DifferentInConditions_motif_graphs.pdf")
  
  for(i in 1:nrow(zonly)){
    name = paste("#",as.character(i), ": ", zonly[i, "motif_id"],", cv = ", as.character(round(zonly$z_coef_var[i],2)), ", mean = ",as.character(round(zonly[i, "zmean"], 2)) , sep="", collapse="")
    plotMotif(zonly$motif_str[i], name)
    plotZscoreByCond(zonly[i, ], colInd)
    
  }
  dev.off()
}


plotMotif <- function(m, motname=2){
  motif = mot2mat(m)
  motmat <- motif$mat
  x = c()
  all_colors <- c("lightblue", "gray", "red", "green", "black")
  colors <- all_colors[motif$type]
  for(i in 1:nrow(motmat)){
    for(j in 1:ncol(motmat)){
      if(motmat[i, j]==1) x <- c(x, j, i)
    }
  }
  g2 <- graph( edges=x, n=nrow(motmat) )
  plot(g2, vertex.color=colors, vertex.size=30, main=motname)
  return()
}

mot2mat <- function(m){
  a <- strsplit(m, "\\|")[[1]]
  #matrix
  motmat <- gsub("[a-zA-Z\\.:]+", "", a[grep("mat", a)])
  motmat <- gsub("\\[\\[","" , motmat)
  motmat <- gsub("\\]\\]","" , motmat)
  motmat <- strsplit(motmat, "\\]\\[")[[1]]
  motmat <- strsplit(split = "", motmat) %>% sapply(as.numeric) %>% t
  #types
  t <- a[grep("types", a)]
  t <- gsub("[a-zA-Z\\.:]+", "", t)
  t <- gsub("\\[","" , t)
  t <- gsub("\\]","" , t)
  t <- strsplit(t, " ")[[1]] %>% as.numeric
  #return:
  return(list(mat=motmat, type=t, type_num=a[3]))
}

makeHeatmap<-function(zonly, min_cv = 10){
  prop <- zonly[, grep("_zscore$", names(zonly))] 
  prop <- prop[zonly$z_coef_var>min_cv, ]
  ## hmap positions
  library(gplots)
  posmat<- prop %>% as.matrix %>% set_colnames(gsub("4cell_|_zscore$", "", names(prop))) %>% set_rownames(zonly$motif_id[zonly$z_coef_var>min_cv])
  #posmat <- posmat[prop ]
  # posmat <- apply(posmat, MAR=1, FUN=function(x){x/mean(x)})
  my_palette <- colorRampPalette(c("red", "blue"))(n = 299)
  png(paste("motDifConditions_ZscoreThreshold", as.character(any_threshold),"_heatmap.png", sep="", collapse=""), width=700, height = 1000)
  hh<-heatmap.2(posmat,
                #keysize=0.5,
                #cellnote = mat_data,  # same data set for cell labels
                main = "Type frequency in motif positions", # heat map title
                #notecol="black",      # change font color of cell labels to black
                density.info="none",  # turns off density plot inside color legend
                trace="none",         # turns off trace lines inside the heat map
                margins =c(12,15),     # widens margins around plot
                col=my_palette,       # use on color palette defined earlier
                #breaks=col_breaks,    # enable color transition at specified limits
                dendrogram="col",     # only draw a row dendrogram
                srtCol = 45,
                cexRow= 2.2,
                cexCol=2.2,
                Colv=TRUE,
                Rowv=TRUE)            # turn off column clustering
  hh
  dev.off()
  return(hh$colInd)
}

plotZscoreByCond <- function(zrow, colInd){
  mname <- zrow$motif_id
  zrow <- zrow[, grep("mce", names(zrow))]
  z2 = data.frame(condition=names(zrow), z_score = zrow %>% as.vector %>% unlist  )
  #names(z2) <- gsub("_", " ", names(z2))
  z2$condition <- gsub("4cell_", "", z2$condition)
  z2$condition <- gsub("_zscore", "", z2$condition)
  z2 <- z2[colInd,]
  z2$mod <- gsub("fix|_mutb|inh5", "", z2$condition)
  g1 <- ggplot(data = z2, aes(x = condition, y=z_score, fill=mod)) + 
    geom_bar( stat = "identity", position="dodge") + coord_flip() +
    
    # scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
    theme(plot.title = element_text(size = rel(1.5), colour = "black", face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(axis.text.y = element_text(size = 15, colour = "black", angle = 0, face = "bold"))+ 
    theme(axis.text.x = element_text(margin = margin(t= 7),size = 15, colour = "black", angle = 0, face = "bold"))+
    # guides(fill=type2) +
    theme(axis.title.x = element_text(size = 15, colour = "black", angle = 0, face = "bold"))+
    scale_fill_brewer(palette="Blues") +
    # scale_color_manual(values=c("black")) +
    theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black", 
                                                                                                            size = 1.5, linetype = "solid"), 
           panel.background = element_blank())+
    #scale_fill_manual(values=colores) +
    labs(title = paste("Z score of ", mname, sep="", collapse=""), 
         y = "Condition", x = "Z score", face = "bold") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title=element_text(size=18,face="bold"))
  
  grid.arrange(g1, nrow = 1)
  
}



setwd("C:/Users/Carlos/Desktop/cellEvolver_all/1802_cellevolver/all_finaldata/motifs_sep1_sepmot3_selby0_filtTrue_ignoreSelfReg_fixed1/180701all_tsne")
f <- "all_Zscores.csv"

any_threshold = 3

zscores <- read.table(f, sep=",", header=TRUE, stringsAsFactors = F, row.names=1)
names(zscores) <- gsub("^X", "", names(zscores))
zscores <- zscores[, -grep("mce0Xss", names(zscores))]
zonly <- zscores[, grep("_zscore", names(zscores))]
zcvar <- apply(zonly %>% as.matrix, MAR=1, FUN=function(x) sd(x)/abs(mean(x)))
zmean <- apply(zonly %>% as.matrix, MAR=1, FUN=mean)
zmean <- apply(zonly %>% as.matrix, MAR=1, FUN=mean)
zany <- apply(zonly %>% as.matrix, MAR=1, FUN=function(x)any(x > any_threshold))



zonly <- cbind(zscores[,1:3], zonly, data.frame(z_coef_var=zcvar, zmean = zmean, z_any_greater_than_threshold = zany))
zonly <- zonly[zany, ]
zonly <- zonly[order(zonly$z_coef_var, decreasing = TRUE), ]

write.table(zonly, file="Zscores_anyGreatherThan3_orderedByVariationBetweenConditions_noXss.csv", sep="\t", row.names=F, quote=F)

makeMotifPlots(zonly)
