library(magrittr)
library(ggplot2)
library(gridExtra)
library(scales)# to use personalized scales in ggplot2

mapcols<-function(types, col2map){#this is not used
  return(sapply(types, FUN=function(x){col2map$col[col2map$t==x]}))
}
makeRegnumberHist <- function(cond, f){
  regnum <- read.table(f, sep=",", dec=".",header=T,stringsAsFactors=F,row.names=1)
  for(t in names(regnum)[grep("type[0-9]", names(regnum))]){
    regnum[, t] <- sapply(regnum[, t], FUN=function(x) type1define[as.character(x)]) %>% unlist
  }
  names(regnum)[names(regnum)== "other_activator"] <- "other_activators"
  names(regnum)[names(regnum)== "inhibitor"] <- "inhibitors"
  tf_cols <-  c("lineage_this_cell", "lineage_other_cell", "other_activators", "inhibitors", "all_activators")
  reg2 <- data.frame()
  for(tt in tf_cols){
    this <- cbind(regnum[, setdiff(names(regnum), tf_cols)], data.frame(tf_type = rep(gsub("_", " ", tt), nrow(regnum)), number_regulators = regnum[, tt], stringsAsFactors = F))
    reg2 <- rbind(reg2, this)
  }
  
  reg2 <- reg2[grep("terminal", reg2$type2 ), ]
  reg2 <- reg2[-grep("other", reg2$type2 ), ]
  #reg2 <- reg2[!(grepl("other", reg2$type2 ) & grepl("lineage", reg2$type2 )) , ]
  mu <- aggregate(reg2$number_regulators, by = list(reg2$type2, reg2$tf_type, reg2$number_regulators), FUN= length)
  names(mu) <- c("type2", "tf_type", "number_regulators", "sum")
  mu$percent <- numeric(nrow(mu))
  for(t in unique(mu$type2)){
    for(tf in unique(mu$tf_type)){
      ind <- which(mu$tf_type == tf & mu$type2==t)
      mu$percent[ind] <- 100*mu$sum[ind]/(sum(mu$sum[ind]))
    }}
  
  mu$number_regulators <- as.numeric(as.character(mu$number_regulators))
  g1 <- ggplot(data = mu, aes(x = number_regulators, y=percent, fill=type2, group=type2, linetype=type2)) + 
    #geom_histogram(position="dodge", binwidth = 1, center=0, bins = 0:10) +
    #geom_bar(aes(y = 100*(..count..)/sum(..count..)), position="dodge") +
    geom_bar(stat = "identity", position="dodge") +
    scale_x_continuous(breaks=seq(0,10,1))+
    # geom_vline(data=mu, aes(xintercept="mean_number_regulators", color=type2),
    #             linetype="dashed")+
    facet_grid( tf_type~ .)+
    scale_fill_manual(values = kols)+
    # scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
    theme(plot.title = element_text(size = rel(1.5), colour = "black", face = "bold")) +
    theme(strip.text.y = element_text(size = 10, colour = "black", angle = 270, face = "bold")) + 
    theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0, face = "bold")) + 
    theme(legend.title = element_text(face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(axis.text.y = element_text(size = 10, colour = "black", angle = 0, face = "bold"))+ 
    theme(axis.text.x = element_text(margin = margin(t= 7),size = 12, colour = "black", angle = 0, face = "bold"))+
    # guides(fill=type2) +
    theme(axis.title.x = element_text(size = 10, colour = "black", angle = 0, face = "bold"))+
    #scale_fill_brewer(palette="Dark2") +
    theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black", 
                                                                                                            size = 1.5, linetype = "solid"), 
           panel.background = element_blank())+
    #scale_fill_manual(values=colores) +
    labs(title = "Regulators by gene in a given cell", 
         y = "Percent", x = "number of active regulators with one or more motifs", fill = "type", face = "bold") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title=element_text(size=14,face="bold"))
  
  mu2 <- mu[mu$tf_type %in% c("all activators", "inhibitors"), ]
  g2 <- ggplot(data = mu2, aes(x = number_regulators, y=percent, fill=type2, group=type2, linetype=type2)) + 
    #geom_histogram(position="dodge", binwidth = 1, center=0, bins = 0:10) +
    #geom_bar(aes(y = 100*(..count..)/sum(..count..)), position="dodge") +
    geom_bar(stat = "identity", position="dodge") +
    scale_fill_manual(values = kols)+
    scale_x_continuous(breaks=seq(0,10,1))+
    # geom_vline(data=mu, aes(xintercept="mean_number_regulators", color=type2),
    #             linetype="dashed")+
    facet_grid( tf_type~ .)+
    # scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
    theme(plot.title = element_text(size = rel(1.5), colour = "black", face = "bold")) +
    theme(strip.text.y = element_text(size = 15, colour = "black", angle = 270, face = "bold")) + 
    theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0, face = "bold")) + 
    theme(legend.title = element_text(face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(axis.text.y = element_text(size = 10, colour = "black", angle = 0, face = "bold"))+ 
    theme(axis.text.x = element_text(margin = margin(t= 7),size = 12, colour = "black", angle = 0, face = "bold"))+
    # guides(fill=type2) +
    theme(axis.title.x = element_text(size = 10, colour = "black", angle = 0, face = "bold"))+
    #scale_fill_brewer(palette="Dark2") +
    theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black", 
                                                                                                            size = 1.5, linetype = "solid"), 
           panel.background = element_blank())+
    #scale_fill_manual(values=colores) +
    labs(title = "Regulators by gene in a given cell", 
         y = "Percent", x = "number of active regulators with one or more motifs", fill = "type", face = "bold") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title=element_text(size=14,face="bold"))
  
  pdf(paste(cond, "_regulatorsByGene_noother.pdf", sep="", collapse=""))
  grid.arrange(g1, nrow = 1)
  grid.arrange(g2, nrow = 1)
  dev.off()
  
}

readCorTab <- function(this_files, tabname, removeTFs=TRUE, removeTerminalOther=TRUE){
  file <- this_files[grep(tabname, this_files)]
  tt <- read.table(file, sep=",", dec=".",header=T,stringsAsFactors=F,row.names=1)
  if(removeTerminalOther){
    tt <- tt[, !grepl("_other", names(tt))]
  }
  tt2 <- data.frame()
  for(n in 4:ncol(tt)){
    aux <- tt[, c(1:3, n)]
    name <- names(aux)[4]
    aux[, "direction"] <- ifelse(grepl("Act", name), "activator set", "inhibitor set")
    names(aux)[4] <- "correlation"
    name <- gsub("Act_|Inh_", "", name)
    if(grepl("selfcor$", name)){
      name <- gsub("_selfcor", "", name)
      aux[, "typeA"] <- gsub("_", " ", name)
      aux[, "typeB"] <- gsub("_", " ", name)
      aux[, "original_var"] <- paste(name, " vs ", name, sep="", collapse="")
    }else{
      name <- strsplit(gsub("_", " ", name),split="\\.")[[1]]
      aux[, "typeA"] <- name[1]
      aux[, "typeB"] <- name[2]
      aux[, "original_var"] <- paste(name[1], " vs ", name[2], sep="", collapse="")
    }
    tt2 <- rbind(tt2, aux)
  }
  return(tt2)
}

plotCorTab<-function(tab, title="Distribution of correlation"){
  g1<-ggplot(data = tab, aes(x = correlation, fill=original_var, alpha=0.4)) + 
    geom_density(size=1, bw=0.1)  + 
    # scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
    facet_grid( direction~ .)+
    theme(strip.text.y = element_text(size = 13, colour = "black", angle = 270, face = "bold")) + 
    theme(strip.text.x = element_text(size = 0.8, colour = "black", angle = 0, face = "bold")) + 
    theme(plot.title = element_text(size = rel(1.5), colour = "black", face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(axis.text.y = element_text(size = 15, colour = "black", angle = 0, face = "bold"))+ 
    theme(axis.text.x = element_text(margin = margin(t= 7),size = 10, colour = "black", angle = 90, face = "bold"))+
    #guides(fill=FALSE) +
    theme(axis.title.x = element_text(size = 13, colour = "black", angle = 0, face = "bold"))+
    #scale_fill_brewer(palette ="Paired", values = c(0, 0.4, 0.6)) +
    theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black", 
                                                                                                            size = 1.5, linetype = "solid"), 
           panel.background = element_blank())+
    #scale_fill_manual(values=colores) +
    labs(title = title, 
         y = "density", x = "Correlation between gene sets in a cell", fill = "type", face = "bold") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title=element_text(size=14,face="bold"))
  
  return(g1)
}


readPhenotypeTab <- function(this_files, tabname, removeTFs=TRUE, removeTerminalOther=TRUE){
  file <- this_files[grep(tabname, this_files)]
  tt <- read.table(file, sep=",", dec=".",header=T,stringsAsFactors=F,row.names=1)
  tt2 <- data.frame()
  for(n in 9:ncol(tt)){
    aux <- tt[, c(1:8, n)]
    type <- names(tt)[n]
    names(aux)[9] <- "phenotype"
    aux[, "tf_mutated_type"] <- type
    tt2 <- rbind(tt2, aux)
  }
  if(removeTFs) tt2 <- tt2[tt2$type3 != 1, ]
  if(removeTerminalOther) tt2 <- tt2[tt2$type2 != 16 & tt2$type2 != 15, ]
  for(i in paste("type", as.character(0:3), sep="")){
    tt2[, i] <- as.character(tt2[, i])
    tt2[, i] <- sapply(tt2[, i], FUN=function(x){type1define[[x]]})
  }
  tt2$tf_mutated_type <- gsub("_", " ", tt2$tf_mutated_type)
  tt2$tf_mutated_type[tt2$tf_mutated_type == "inhibitor"] <- "inhibitors"
  tt2$tf_mutated_type[tt2$tf_mutated_type == "other activator"] <- "other activators"
  return(tt2)
}

plotPhenotypeDifferences<- function(sites, tfs, cond){
  tfs$diff_to_sites <- abs(tfs$phenotype) - abs(sites$phenotype)
  
  g1<-ggplot(data = tfs, aes(x = diff_to_sites, fill=type2, alpha=0.4)) + 
    geom_density(size=1, bw=0.1)  + 
    # scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
    scale_fill_manual(values = kols)+
    facet_grid( tf_mutated_type~.)+
    geom_vline(xintercept=0, size=1,linetype="dotted", color = "red") +
    theme(strip.text.y = element_text(size = 9, colour = "black", angle = 270, face = "bold")) + 
    #theme(strip.text.x = element_text(size = 0.8, colour = "black", angle = 0, face = "bold")) + 
    theme(plot.title = element_text(size = rel(1.5), colour = "black", face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(axis.text.y = element_text(size = 15, colour = "black", angle = 0, face = "bold"))+ 
    theme(axis.text.x = element_text(margin = margin(t= 7),size = 10, colour = "black", angle = 90, face = "bold"))+
    #guides(fill=FALSE) +
    theme(axis.title.x = element_text(size = 13, colour = "black", angle = 0, face = "bold"))+
    #scale_fill_brewer(palette ="Paired", values = c(0, 0.4, 0.6)) +
    theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black", 
                                                                                                            size = 1.5, linetype = "solid"), 
           panel.background = element_blank())+
    #scale_fill_manual(values=colores) +
    labs(title = "Absolute differences between TF and TFBS mutation", 
         y = "density", x = "|TFphen|-|TFBSphen|", fill = "type", face = "bold") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title=element_text(size=14,face="bold")) 
  
  pdf(paste(cond, "mutantPhenotypeDifferenceTF-Sites.pdf"))
  grid.arrange(g1, nrow = 1)
  dev.off()
}
plotPhenotypeDensities <- function(sites, tfs, cond, all=FALSE){
  
  g1<-ggplot(data = sites, aes(x = phenotype, fill=type2, alpha=0.4)) + 
    geom_density(size=1, bw=0.1)  + 
    # scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
    scale_fill_manual(values = kols)+
    facet_grid( tf_mutated_type~.)+
    geom_vline(xintercept=0, size=1,linetype="dotted", color = "red") +
    theme(strip.text.y = element_text(size = 9, colour = "black", angle = 270, face = "bold")) + 
    #theme(strip.text.x = element_text(size = 0.8, colour = "black", angle = 0, face = "bold")) + 
    theme(plot.title = element_text(size = rel(1.5), colour = "black", face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(axis.text.y = element_text(size = 15, colour = "black", angle = 0, face = "bold"))+ 
    theme(axis.text.x = element_text(margin = margin(t= 7),size = 10, colour = "black", angle = 90, face = "bold"))+
    #guides(fill=FALSE) +
    theme(axis.title.x = element_text(size = 13, colour = "black", angle = 0, face = "bold"))+
    #scale_fill_brewer(palette ="Paired", values = c(0, 0.4, 0.6)) +
    theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black", 
                                                                                                            size = 1.5, linetype = "solid"), 
           panel.background = element_blank())+
    #scale_fill_manual(values=colores) +
    labs(title = "Distribution of phenotype - TFBS mutation", 
         y = "density", x = "Difference to wildtype in gene expression", fill = "type", face = "bold") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title=element_text(size=14,face="bold")) 
  
  g2<-ggplot(data = tfs, aes(x = phenotype, fill=type2, alpha=0.4)) + 
    geom_density(size=1, bw=0.1)  + 
    # scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
    scale_fill_manual(values = kols)+
    facet_grid( tf_mutated_type~ .)+
    geom_vline(xintercept=0, size=1,linetype="dotted", color = "red") +
    theme(strip.text.y = element_text(size = 9, colour = "black", angle = 270, face = "bold")) + 
    #theme(strip.text.x = element_text(size = 0.8, colour = "black", angle = 0, face = "bold")) + 
    theme(plot.title = element_text(size = rel(1.5), colour = "black", face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(axis.text.y = element_text(size = 15, colour = "black", angle = 0, face = "bold"))+ 
    theme(axis.text.x = element_text(margin = margin(t= 7),size = 10, colour = "black", angle = 90, face = "bold"))+
    #guides(fill=FALSE) +
    theme(axis.title.x = element_text(size = 13, colour = "black", angle = 0, face = "bold"))+
    #scale_fill_brewer(palette ="Paired", values = c(0, 0.4, 0.6)) +
    theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black", 
                                                                                                            size = 1.5, linetype = "solid"), 
           panel.background = element_blank())+
    #scale_fill_manual(values=colores) +
    labs(title = "Distribution of phenotype - TF mutation", 
         y = "density", x = "Difference to wildtype in gene expression", fill = "type", face = "bold") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title=element_text(size=14,face="bold")) 
  pdfname = ifelse(all,paste(cond, "mutantPhenotypeDistribution_all.pdf") ,paste(cond, "mutantPhenotypeDistribution.pdf"))
  pdf(pdfname)
  grid.arrange(g1, nrow = 1)
  grid.arrange(g2, nrow = 1)
  dev.off()
  
}


getMotifPositionMostBiased<- function(cond, this_files, pairoftypes = c("terminal_specific_this", "terminal_all")){
    library(reshape)
    summ <- preprocessSum(cond, this_files)
    summ$ratio <- (summ[, cond]+1)/(summ[, "random_mean"]+1)
    summ <- summ[summ$z_score >= min_z_score & summ[, cond]>=threshold, ]
    prop <- preprocessProp(cond, this_files, summ$motif_id)
    mots <- summ$motif_id[summ$z_score >= min_z_score]
    prop2 <- prop[, setdiff(names(prop), c("type2"))] %>% melt(id = c("typename")) %>% cast(formula = variable~typename)
    propCp <- prop2
    types <- unique(prop$typename)
    for(i in 1:(length(types)-1)){for(j in (i+1):length(types)){
      prop2[, paste(types[i], types[j], sep=";", collapse=";")] <- abs(prop2[, types[i]] - prop2[, types[j]])
    }}
  col <- sapply(names(prop2), FUN=function(x){grepl(pairoftypes[1], x) & grepl(pairoftypes[2], x)}) %>% which
  prop2 <- prop2[order(prop2[, col], decreasing=T), ]
  prop2$z_score <- sapply(prop2$variable %>% as.character, FUN=function(x, y, z){
    z[y == strsplit(x, "_")[[1]][1]]
  }, y=summ$motif_id, summ$z_score)
  write.table(prop2, file= paste("180709", cond, "_proportion_differences_ordered",sep = "_", collapse="_"))
  pp <- prop2[, names(prop2) %in% c("variable","z_score", pairoftypes)]
  pp$condition <- cond
  pp[, names(prop2)[col]] <- prop2[, col]
  return(pp)
}

makeMotifPlots <- function(cond, this_files, order_by_Zscore=FALSE, max_num=20,vertex.size=50, edge.width=4, edge.arrow.size=2, edge.arrow.width=1.8, label.cex=3){
  summ <- preprocessSum(cond, this_files)
  summ$ratio <- (summ[, cond]+1)/(summ[, "random_mean"]+1)
  summ <- summ[summ$z_score >= min_z_score & summ[, cond]>=threshold, ]
  prop <- preprocessProp(cond, this_files, summ$motif_id)
  makeHeatmap(prop)
  
  #make prop table for join plot; add order from clustering
  df = data.frame()
  for(t in names(prop)[3:ncol(prop)]){
    aux <- cbind(prop[, c(1,2)], data.frame(motif_position_name = t, type_composition = prop[, t], stringsAsFactors=F))
    df <- rbind(df, aux)
  }
  df$motif_name <- unlist(sapply( df$motif_position_name, FUN=function(x)strsplit(x, "_")[[1]][1]))
  
  
  library(igraph)
  library(gridExtra)
  ## make plots for each motif
  pdf(paste(cond, "_motif_graphs.pdf"))
  if(order_by_Zscore){
    summ <- summ[order(summ$z_score, decreasing = TRUE), ]
  }
  for(size in 3:max(summ$motif_size)){
    this_summ <- summ[summ$motif_size==size, ]
    for(i in 1:nrow(this_summ)){
      if(i> max_num) break
      name = paste("#",as.character(i), ": ", this_summ$motif_id[i], ", ",as.character(this_summ[i, cond]) ," matches, Z-score=", as.character(round(this_summ$z_score[i], 2)), sep="", collapse="")
      this_mot <- df[this_summ$motif_id[i]==df$motif_name, ]
      this_mot$motif_position_name <- sapply(this_mot$motif_position_name, FUN=function(x){
        y<-strsplit(x, "_")[[1]]
        return(y[length(y)])
      })
      plotMotif(this_summ$motif_str[i], name, vertex.size, edge.width, edge.arrow.size, edge.arrow.width, label.cex)
      plotMotComposition(this_mot, name)
      
    }}
  dev.off()
}


plotMotif <- function(m, motname=2, vertex.size=50, edge.width=4, edge.arrow.size=2, edge.arrow.width=1.8, label.cex=3){
  motif = mot2mat(m)
  motmat <- motif$mat
  x = c()
  all_colors <- c("Skyblue2", "red", "green", "black", "gray")
  colors <- all_colors[motif$type]
  for(i in 1:nrow(motmat)){
    for(j in 1:ncol(motmat)){
      if(motmat[i, j]==1) x <- c(x, j, i)
    }
  }
  g2 <- graph( edges=x, n=nrow(motmat) )
  plot(g2, vertex.color=colors, vertex.size=vertex.size, edge.width=edge.width, edge.color="black", edge.arrow.size=edge.arrow.size, edge.arrow.width=edge.arrow.width,
       vertex.label.cex = label.cex, main=motname, vertex.label.color="black", vertex.frame.width=10)
  return()
}

plotMotComposition<-function(mottab, name=""){
  kolslocal <- kols
  names(kolslocal) <- gsub("tors$", "tor", names(kolslocal))
  mottab$typename <- gsub("_", " ", mottab$typename)
  g<-ggplot(data = mottab, aes(x = motif_position_name, y = type_composition, fill=typename, group=typename)) + 
    geom_bar( stat = "identity")  + 
    # scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
    scale_fill_manual(values = kolslocal)+
    theme(plot.title = element_text(size = rel(1.5), colour = "black", face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(axis.text.y = element_text(size = 24, colour = "black", angle = 0, face = "bold"))+ 
    theme(axis.text.x = element_text(margin = margin(t= 7),size = 30, colour = "black", angle = 0, face = "bold"))+
    #guides(fill=FALSE) +
    theme(axis.title.x = element_text(size = 10, colour = "black", angle = 0, face = "bold"))+
    #scale_fill_brewer(palette ="Paired", values = c(0, 0.4, 0.6)) +
    theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black", 
                                                                                                            size = 1.5, linetype = "solid"), 
           panel.background = element_blank())+
    #scale_fill_manual(values=colores) +
    labs(title = name, 
         y = "", x = "", fill = "type", face = "bold") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title=element_text(size=14,face="bold")) 
  grid.arrange(g, nrow = 1)
  
}


makeHeatmap<-function(prop){
  ## hmap positions
  library(gplots)
  posmat<-as.matrix(prop[, 3:ncol(prop)])%>% t %>% set_colnames(gsub("_", " ", prop$typename))
  my_palette <- colorRampPalette(c("black", "white"))(n = 299)
  pdf(paste(cond, "_heatmap_propType.pdf", sep="", collapse=""))
  heatmap.2(posmat,
            #keysize=0.5,
            #cellnote = mat_data,  # same data set for cell labels
            main = "Type frequency in motif positions", # heat map title
            #notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(12,9),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier
            #breaks=col_breaks,    # enable color transition at specified limits
            dendrogram="row",     # only draw a row dendrogram
            cexRow= 1,
            cexCol= 1.2,
            Colv="NA")            # turn off column clustering
  dev.off()
}

preprocessProp<-function(cond, this_files, getmots){
  prop <- read.table(this_files[grep(motprop, this_files)], sep=",", dec=".",header=T,stringsAsFactors=F,row.names=1)
  find_mos <- sapply(getmots, FUN=function(x, y)grep(paste(x, "_", sep="", collapse=""),y), y = names(prop)) %>% unlist %>% unique
  propred <- cbind(prop[, grepl("type",names(prop))], prop[, find_mos])
  return(propred)
}

preprocessSum <- function(cond, this_files=NULL, summ=NULL, paste_name = FALSE){
  if(is.null(summ)){
    summ <- read.table(this_files[grep(motsum, this_files)], sep=",", dec=".",header=T,stringsAsFactors=F,row.names=1)
    names(summ) <- gsub("^X","", names(summ))
  }
  to_remove = c('4cell_mce0random_5', '4cell_mce0fixrandom_4', '4cell_mce0fixrandom_3', '4cell_mce0_mutbrandom_2', '4cell_mce0fixbrandom_5', '4cell_mce0fixbrandom_6')
  summ <- summ[, !names(summ) %in% to_remove]
  order_mots <- order(summ[, cond], decreasing=T)
  summ <- summ[order_mots,]
  #retain <- which(summ[, cond] >= threshold)
  rand <- as.matrix(summ[, grepl(cond, names(summ)) & grepl("random_", names(summ))])
  summ$random_mean <- apply(rand,MAR= 1, mean)
  summ$random_sd <-  apply(rand,MAR= 1, sd)
  summ$z_score <- (summ[, cond] - summ$random_mean + 1)/(summ$random_sd + 1)
  summ$z_score[is.na(summ$z_score)] <- 0
  
  
  summ$norm_z_score <- numeric(nrow(summ))
  for(n in unique(summ$motif_size)){
    ind <- which(summ$motif_size == n)
    summ$norm_z_score[ind] <- (summ$z_score[ind]+1)/(sum(sqrt(summ$z_score[ind]^2))+1)
  }
  if(paste_name){
    names(summ)[names(summ) %in% c("random_mean", "random_sd", "z_score", "norm_z_score")] <- paste(c("random_mean", "random_sd", "z_score", "norm_z_score"), cond, sep="_")
  }
  return(summ)
}

plotMotifMatches <-function(cond, this_files){
  summ <- preprocessSum(cond, this_files)
  ## build dataframe to plot
  
  df <- data.frame()
  for(s in 3:max_mot_size){
    this_df <- summ[summ$motif_size==s, ]
    df_a <- data.frame(mot_id = this_df$motif_id, mot_order = 1:nrow(this_df),motif_size = this_df$motif_size, networks = cond, matches = this_df[, cond])
    df_b <- data.frame(mot_id = this_df$motif_id, mot_order = 1:nrow(this_df),motif_size = this_df$motif_size, networks = paste("randomized ", cond, sep="", collapse=""), matches = this_df[, "random_mean"])
    df = rbind(df, df_a, df_b)
    
  }
  
  library(ggplot2)
  glist <- list()
  for(s in 3:max_mot_size){
    dfred <- df[df$motif_size==s, ]
    glist[[as.character(s)]]<-ggplot(data = dfred, aes(x = mot_order, y = matches, color=networks, group=networks)) + 
      geom_line( size=1.5)  + 
      # scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
      theme(plot.title = element_text(size = rel(1.5), colour = "black", face = "bold")) +
      theme(strip.text.y = element_text(size = 15, colour = "black", angle = 270, face = "bold")) + 
      theme(strip.text.x = element_text(size = 18, colour = "black", angle = 0, face = "bold")) + 
      theme(legend.title = element_text(face = "bold")) +
      theme(legend.title = element_text(face = "bold")) +
      theme(axis.text.y = element_text(size = 15, colour = "black", angle = 0, face = "bold"))+ 
      theme(axis.text.x = element_text(margin = margin(t= 7),size = 12, colour = "black", angle = 0, face = "bold"))+
      guides(fill=FALSE) +
      theme(axis.title.x = element_text(size = 10, colour = "black", angle = 0, face = "bold"))+
      #scale_fill_brewer(palette ="Paired", values = c(0, 0.4, 0.6)) +
      theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black", 
                                                                                                              size = 1.5, linetype = "solid"), 
             panel.background = element_blank())+
      #scale_fill_manual(values=colores) +
      labs(title = paste("Number of matches of motifs size ", as.character(s), sep="", collapse=""), 
           y = "num. of matches", x = "motif", fill = "type", face = "bold") +
      theme(plot.title = element_text(hjust = 0.5), 
            axis.title=element_text(size=14,face="bold")) 
  }
  library(gridExtra)
  pdf(paste(cond, '_matches_in_sims.pdf', sep="", collapse=""))
  grid.arrange(glist[[1]], glist[[2]], glist[[3]], nrow = 3)
  dev.off()
  
  #Now plot Z-score
  glist <- list()
  for(s in 3:max_mot_size){
    dfred <- summ[summ$motif_size==s, ]
    dfred <- dfred[order(dfred$z_score, decreasing=TRUE),]
    dfred$mot_order <- order(dfred$z_score, decreasing=TRUE)
    glist[[as.character(s)]]<-ggplot(data = dfred, aes(x = mot_order, y = z_score)) + 
      geom_line( size=1.5, col = "black")  + 
      geom_hline(yintercept=min_z_score, size=1.5,linetype="dashed", color = "red") +
      #geom_rect(xmin=-Inf,xmax=Inf, ymin = -1.5, ymax=1.5, alpha=0.25, fill = "red") +
      #geom_hline(yintercept=2, linetype="dashed", color = "red") +
      #geom_hline(yintercept=5, linetype="dashed", color = "red") +
      # scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
      theme(plot.title = element_text(size = rel(1.5), colour = "black", face = "bold")) +
      theme(strip.text.y = element_text(size = 15, colour = "black", angle = 270, face = "bold")) + 
      theme(strip.text.x = element_text(size = 18, colour = "black", angle = 0, face = "bold")) + 
      theme(legend.title = element_text(face = "bold")) +
      theme(legend.title = element_text(face = "bold")) +
      theme(axis.text.y = element_text(size = 15, colour = "black", angle = 0, face = "bold"))+ 
      theme(axis.text.x = element_text(margin = margin(t= 7),size = 12, colour = "black", angle = 0, face = "bold"))+
      guides(fill=FALSE) +
      theme(axis.title.x = element_text(size = 10, colour = "black", angle = 0, face = "bold"))+
      #scale_fill_brewer(palette ="Paired", values = c(0, 0.4, 0.6)) +
      theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black", 
                                                                                                              size = 1.5, linetype = "solid"), 
             panel.background = element_blank())+
      #scale_fill_manual(values=colores) +
      labs(title = paste("size ", as.character(s), sep="", collapse=""), 
           y = "Z-score", x = "subgraph", fill = "type", face = "bold") +
      theme(plot.title = element_text(hjust = 0.5), 
            axis.title=element_text(size=14,face="bold")) 
  }
  library(gridExtra)
  pdf(paste(cond, '_matches_in_sims_Zscore.pdf', sep="", collapse=""))
  grid.arrange(glist[[1]], glist[[2]], glist[[3]], nrow = 3)
  dev.off()
}


makeZscoreCorrPlot <- function(cond, min_z_score=5){
  to_remove = c('4cell_mce0random_5', '4cell_mce0fixrandom_4', '4cell_mce0fixrandom_3', '4cell_mce0_mutbrandom_2')
  summ <- preprocessSum(cond, this_files, paste_name=TRUE)
  for(c in condlist[condlist != cond]){
    summ <- preprocessSum(c, summ=summ, paste_name=TRUE)
  }
  summ$mean_z_score <- apply(as.matrix(summ[, grep("^z_score", names(summ))]), MAR=1, mean)
  summ$any_greater_z_score <- apply(as.matrix(summ[, grep("^z_score", names(summ))]), MAR=1, FUN=function(x)any(x>min_z_score & ! is.na(x)))
  summ <- summ[summ$any_greater_z_score, ]
  summ <- summ[order(summ$mean_z_score, decreasing = TRUE), ]
  for (i in unique(summ$motif_size)){
    summ$mot_order[summ$motif_size == i] <- 1:nrow(summ[summ$motif_size == i,])
  }

  vtab <- data.frame()
  base <- summ[, c("motif_id", "motif_str", "motif_size", "functional", "mot_order")]
  for (c in condlist){
    regname <- paste("_", c, "$", sep="", collapse="")
    this_df <- cbind(base, summ[, grepl(regname, names(summ)) & grepl("[Zz]_score", names(summ))])
    names(this_df) <- gsub(regname, "", names(this_df))
    this_df$dataset <- c
    vtab <- rbind(vtab, this_df)
  }
  
  glist <- list()
  df <- vtab
  for(s in 3:max_mot_size){
    dfred <- df[df$motif_size==s, ]
    glist[[as.character(s)]]<-ggplot(data = dfred, aes(x = mot_order, y = z_score, color=dataset, group=dataset, linetype=dataset)) + 
      geom_line( size=1.5)  + 
      # scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
      theme(plot.title = element_text(size = rel(1.5), colour = "black", face = "bold")) +
      theme(strip.text.y = element_text(size = 15, colour = "black", angle = 270, face = "bold")) + 
      theme(strip.text.x = element_text(size = 18, colour = "black", angle = 0, face = "bold")) + 
      theme(legend.title = element_text(face = "bold")) +
      theme(legend.title = element_text(face = "bold")) +
      theme(axis.text.y = element_text(size = 15, colour = "black", angle = 0, face = "bold"))+ 
      theme(axis.text.x = element_text(margin = margin(t= 7),size = 12, colour = "black", angle = 0, face = "bold"))+
      guides(fill=FALSE) +
      theme(axis.title.x = element_text(size = 10, colour = "black", angle = 0, face = "bold"))+
      #scale_fill_brewer(palette ="Paired", values = c(0, 0.4, 0.6)) +
      theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black", 
                                                                                                              size = 1.5, linetype = "solid"), 
             panel.background = element_blank())+
      #scale_fill_manual(values=colores) +
      labs(title = paste("Z-score of motifs size ", as.character(s), sep="", collapse=""), 
           y = "Z-score", x = "motif", fill = "type", face = "bold") +
      theme(plot.title = element_text(hjust = 0.5), 
            axis.title=element_text(size=14,face="bold")) +
      theme(legend.text=element_text(size=8, keyheight=4))
  }
  library(gridExtra)
  pdf(paste(cond, '_all_Zscore.pdf', sep="", collapse=""))
  grid.arrange(glist[[1]], glist[[2]], glist[[3]], nrow = 3)
  dev.off()
  
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

motsum <- "motifSummary"
motprop <- "motifsByType2_prop_motifs"
motcount <- "motifsByType2_count"
all_filetypes=list(motsum, motprop, motcount)



condlist = c('4cell_mce0', '4cell_mce0fix_mutb', '4cell_mce1fix', '4cell_mce2fix', '4cell_mce2inh5', '4cell_mce0fix', '4cell_mce0_mutb', '4cell_mce1fix_mutb', '4cell_mce2fix_mutb', '4cell_mce2inh5_mutb', '4cell_mce0Xss')
setwd("C:/Users/Carlos/Desktop/cellEvolver_all/1802_cellevolver/all_finaldata/motifs_sep1_sepmot3_selby0_filtTrue_ignoreSelfReg_fixed1/")

#condlist = c('4cell_mce0fixb')
#setwd("C:/Users/Carlos/Desktop/cellEvolver_all/1802_cellevolver/all_finaldata/180911_rmInhibitors/180914_motifsPartial1")

#condlist = c('4cell_mce0', '4cell_mce0fix_mutb', '4cell_mce1fix', '4cell_mce2fix', '4cell_mce2inh5', '4cell_mce0fix', '4cell_mce0_mutb', '4cell_mce1fix_mutb', '4cell_mce2fix_mutb', '4cell_mce2inh5_mutb')
#setwd("C:/Users/Carlos/Desktop/cellEvolver_all/1802_cellevolver/all_finaldata/motifs_sep1_sepmot3_selby0_filtTrue")

files <- list.files(recursive = TRUE, full.names=TRUE)
files <- files[grepl(motsum, files) |grepl(motprop, files) | grepl(motcount, files)]

threshold <- 10
max_mot_size <- 5
#min_ratio <- 1.5
min_z_score <- 2

## plot motif abundance
for(cond in condlist){
  this_files <- files[grep(paste('\\/', cond, '\\/', sep="", collapse=""), files)]
  plotMotifMatches(cond, this_files)
}

## plot all motifs
for(cond in condlist){ ### cond number 8 was problematic!
  this_files <- files[grep(paste('\\/', cond, '\\/', sep="", collapse=""), files)]
    makeMotifPlots(cond, this_files, TRUE, max_num = 5000, vertex.size=45, edge.width=2.5, edge.arrow.size=1.2, edge.arrow.width=1, label.cex=3)
    cat(cond, "\n")
}

## Get positions that are differentially occupied by types
pairoftypes <- c("terminal_specific_this", "terminal_all")
props <- data.frame()
for (cond in condlist){
  this_files <- files[grep(paste('\\/', cond, '\\/', sep="", collapse=""), files)]
  p <- getMotifPositionMostBiased(cond, this_files,  pairoftypes =pairoftypes)
  props<- rbind(props, p)
}
byzscore <- aggregate(props$z_score, by=list(props$variable), mean) %>% set_names(c("motPosition", "z_score"))
diff <- aggregate(props[, ncol(props)], by=list(props$variable), mean)%>% set_names(c("motPosition", paste(pairoftypes, sep=";", collapse=";")))
finaltab <- merge(byzscore, diff)
finaltab <- finaltab[order(finaltab[, 3], decreasing = T),]
write.table(finaltab, file="180709_motifs_differencesBetweenAllAndSpecific.csv", row.names=F, quote=F, sep="\t")
finfilt <- finaltab[finaltab$z_score>20 & finaltab[, 3] > 0.5, ]
write.table(finaltab, file="180709_motifs_differencesBetweenAllAndSpecific_filtered_z20dif50.csv", row.names=F, quote=F, sep="\t")

## plot Z-score for all conditions
lastcond <- "4cell_mce0Xss"makeZscoreCorrPlot
this_files <- files[grep(paste('\\/', lastcond, '\\/', sep="", collapse=""), files)]
(lastcond)

### plot mean phenotype by gene
library(gridExtra)
sitephen <- "meanPhenotypeByGene_sites.csv"
tfphen <- "meanPhenotypeByGene_tfs.csv"

files <- list.files(recursive = TRUE, full.names=TRUE)
files <- files[grepl(sitephen, files) |grepl(tfphen, files)]

for(cond in condlist){
  this_files <- files[grep(paste('\\/', cond, '\\/', sep="", collapse=""), files)]
  sites <- readPhenotypeTab(this_files, sitephen)
  sitesall <- readPhenotypeTab(this_files, sitephen, removeTerminalOther=FALSE)
  tfs <- readPhenotypeTab(this_files, tfphen)
  tfsall <- readPhenotypeTab(this_files, tfphen, removeTerminalOther=FALSE)
  #plotPhenotypeDensities(sites, tfs, cond)
  plotPhenotypeDensities(sitesall, tfsall, cond, TRUE)
  cat(cond, "\n")
} 
  
#get differences between types
for(cond in condlist){
  this_files <- files[grep(paste('\\/', cond, '\\/', sep="", collapse=""), files)]
  sites <- readPhenotypeTab(this_files, sitephen)
  tfs <- readPhenotypeTab(this_files, tfphen)
  plotPhenotypeDifferences(sites, tfs, condE)
  cat(cond, "\n")
} 
#### Correlations between terminal selector sets


corf <- "regCorrelationsBetweenTypes.csv"
interf <- "regIntersectionBetweenTypes.csv"

files <- list.files(recursive = TRUE, full.names=TRUE)
files <- files[grepl(corf, files) |grepl(interf, files)]

for(cond in condlist){
  this_files <- files[grep(paste('\\/', cond, '\\/', sep="", collapse=""), files)]
  cortab <- readCorTab(this_files, corf)
  interstab <- readCorTab(this_files, interf)
  ga<-plotCorTab(cortab)
  gb<-plotCorTab(interstab, "Intersect(A, B)/(A)")
  pdf(paste(cond, "_correlationBetweenRegulators.pdf", sep="", collapse=""))
    grid.arrange(ga, nrow = 1)
    grid.arrange(gb, nrow = 1)
  dev.off()
  cat(cond, "\n")
  aggregate(interstab$correlation, by=list(interstab$typeA, interstab$typeB), mean, na.rm=T)
} 

### Histogram of number of regulators by gene

fmark <- "regulatorsByGeneByCell.csv"
files <- list.files(recursive = TRUE, full.names=TRUE)
files <- files[grepl(fmark, files)]

for(cond in condlist){
  this_files <- files[grep(paste('\\/', cond, '\\/', sep="", collapse=""), files)]
  makeRegnumberHist(cond, this_files)
  cat(cond, "\n")
} 






