library(magrittr)
library(ggplot2)
library(gridExtra)
library(scales)# to use personalized scales in ggplot2
library(reshape)
library(data.table)
library(pheatmap)


#Plot all phenotypes etc by gene



makeBigPhenotypeTable <- function(allfiles, calc_mode = c("diff", "proportion")){
  bigdf<-data.frame()
  i = 1
  for(sim in unique(allfiles$sim)){
    simfiles<- allfiles[allfiles$sim==sim, ]
    tff <- simfiles[grepl(tfname, simfiles$file), "fullname"]
    sitesf <- simfiles[grepl(sitesname, simfiles$file), "fullname"]
    exprf <- simfiles[grepl(exprname, simfiles$file), "fullname"]
    phh <- processSim(sim, exprf, tff, sitesf, calc_mode)
    bigdf <- rbind(bigdf, phh)
    cat(i, ": ", sim, " successfully attached\n")
    i=i+1
  }
  return(bigdf)
}

processSim <- function(sim, exprf, tff, sitesf, calc_mode = c("diff", "proportion")){
  
  tf <- readFile(tff)
  tf<-tf[order(tf$gene, tf$tf_mutated), ]
  sites <- readFile(sitesf)
  sites<-sites[order(sites$gene, sites$mutated_sites), ]
  expr <- readFile(exprf)
  expr$type <- sapply(1:nrow(expr), FUN=function(i){ #put lineage TF type as a paste of expression in different cells
    aux <- expr[i, grep("init", names(expr))]
    if(any(aux != 0)){
      return(paste(as.character(aux), sep="", collapse=""))
    }else{
      return(expr[i, "type"])
    }
  })
  phh <- getPhenotype(expr, tf, sites, sim,  )
  return(phh)
}


readFile <- function(f2read, var0name="gene"){
  con <- file(f2read,"r")
  first_line <- readLines(con,n=1)
  tf <- read.table(con, skip=0, stringsAsFactors=F, header=F) #first row already read
  close(con)
  
  first_line <- strsplit(first_line, split="\t")[[1]]
  if(first_line[1]=="") first_line[1]<- "gene"
  names(tf) <- first_line
  return(tf)
}

getType2<-function(wt){
  wt$type <- as.character(wt$type)
  res<-sapply(1:nrow(wt), FUN=function(i){
    cell <- as.numeric(gsub("exp", "", wt$cell[i])) 
    if(!grepl("\\[", wt[i, "type"]) & nchar(wt[i, "type"])==length(unique(wt$cell))){
      lin <- which(strsplit(wt$type[i], "")[[1]] != "0") - 1
      if(length(lin)>1){
        return("lineage all")
      }else if(cell == lin){
        return("lineage this")
      }else{
        return("lineage other")
      }
    }else if(!grepl("\\[", wt[i, "type"])){
      if( wt[i, "type"] == "1"){
        return("other activators")
      }else if( wt[i, "type"] == "-1"){
        return("inhibitors")
      }
    }else{
      tt <- gsub("\\[|\\]", "", wt$type[i])
      tt <- strsplit(tt, "")[[1]] %>% as.numeric
      if(length(tt) == length(unique(wt$cell))){
        return("terminal all")
      }else if(length(tt)> 1 & length(tt) < length(unique(wt$cell))){
        this_n_cells<- as.character(length(tt))
        if(cell %in% tt){
          return(paste("terminal ",this_n_cells, " this", sep="", collapse=""))
        }else{
          return(paste("terminal ",this_n_cells, " other", sep="", collapse=""))
        }
      }else{
        if(cell == tt){
          return("terminal specific this")
        }else{
          return("terminal specific other")
        }
      }
    }
  })
  return(res)
}


getPhenotype<- function(expr, tf, sites, simname="", calc_mode = c("diff", "proportion") ){
  wt <- melt(data.table(expr), id.vars=c("gene", "type"), measure.vars = names(expr)[grep("exp", names(expr))], variable_name = "cell")
  names(wt)[3] <- "cell"
  wt$type2 <- getType2(wt)
  
  tf2 <- melt(tf, id.vars=c("gene", "tf_mutated"), measure.vars = names(tf)[grep("exp", names(tf))], variable_name = "cell")
  names(tf2)[3]<-"cell"
  
  s2 <- melt(sites, id.vars=c("gene", "mutated_sites"), measure.vars = names(sites)[grep("exp", names(sites))], variable_name = "cell")
  names(s2)[3]<-"cell"
  
  if(all(tf2$tf_mutated == s2$mutated_sites) & all(tf2$gene == s2$gene)){
    both <- data.frame(gene=tf2$gene, 
                       mutation=tf2$tf_mutated,
                       cell=tf2$cell, 
                       tf_mutated=tf2$value, 
                       mutated_sites=s2$value, 
                       stringsAsFactors=FALSE)
    both2 <- merge(both, wt, by.x = c("gene", "cell"), by.y=c("gene", "cell"))
    if(any(is.na(both2$value))){
      cat("WARNING: NAs in expression data")
    }
    names(both2)[names(both2)=="value"] <- "wt_expression"
    if("diff" %in% calc_mode){
      both2$tf_phenotype <- both2$tf_mutated - both2$wt_expression
      both2$site_phenotype <- both2$mutated_sites - both2$wt_expression
      both2$difference <- abs(both2$tf_phenotype) - abs(both2$site_phenotype)
    }
    if("proportion" %in% calc_mode){
      both2$tf_phenotype_prop <- both2$tf_mutated / both2$wt_expression
      both2$site_phenotype_prop <- both2$mutated_sites / both2$wt_expression
      both2$difference_prop <-  abs(both2$site_phenotype)/ abs(both2$tf_phenotype)
    }
    both2$simname <- as.character(simname)
  }else{
    cat("ERROR: genes or mut TFs different in sites and in TFs ", simname, "\t")
    return(data.frame())
  }
  return(both2)
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


densityplotStruc<-function(df, xname, mainname){
  dfaux <- df[df$tf_mutated_type != "inhibitors", ]
  dfaux$tf_mutated_type <- "all activators"
  df <- rbind(df, dfaux)
  g1<-ggplot(data = df, aes(x = aux_var, fill=type2, alpha=0.4)) + 
    geom_density(size=1, bw=0.1)  + 
    # scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
    scale_fill_manual(values = kols)+
    xlim(c(min(df$aux_var), ifelse(any(df$aux_var < 0),1, 1.25))) +
    facet_grid( tf_mutated_type~.)+
    geom_vline(xintercept=ifelse(any(df$aux_var < 0),0, 1), size=1,linetype="dotted", color = "red") +
    theme(strip.text.y = element_text(size = 8, colour = "black", angle = 270, face = "bold")) + 
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
    labs(title = mainname, y = "density", x = xname, fill = "gene type", face = "bold") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title=element_text(size=14,face="bold")) 
  return(g1)
  
}

BoxplotStruc<-function(df, xname, mainname){
  dfaux <- df[df$tf_mutated_type != "inhibitors", ]
  dfaux$tf_mutated_type <- "all activators"
  df <- rbind(df, dfaux)
  
  g1<-ggplot(data = df, aes(x = type2, y=aux_var, fill=type2)) + 
    #geom_violin(trim=T, alpha=1, bw=0.2)  + 
    stat_boxplot(geom = "errorbar", width = 0.5) +
    geom_boxplot( outlier.size = -1, notch=F)+
    stat_summary(fun.y=mean, geom="point", shape=18, size=3, col="black", na.rm=T)+
    ylim(c(min(df$aux_var), ifelse(any(df$aux_var < 0),1, 1.3))) +
    # scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
    scale_fill_manual(values = kols)+
    facet_grid(.~tf_mutated_type)+
    theme(strip.text.y = element_text(size = 8, colour = "black", angle = 270, face = "bold")) + 
    theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0, face = "bold")) + 
    theme(plot.title = element_text(size = rel(1.5), colour = "black", face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(axis.text.y = element_text(size = 15, colour = "black", angle = 0, face = "bold"))+ 
    theme(axis.text.x = element_text(margin = margin(t= 7),size = 14, colour = "black", angle = 45, face = "bold", hjust=1, vjust=1.08))+
    #guides(fill=FALSE) +
    theme(axis.title.x = element_text(size = 13, colour = "black", angle = 0, face = "bold"))+
    #scale_fill_brewer(palette ="Paired", values = c(0, 0.4, 0.6)) +
    theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black", 
                                                                                                            size = 1.5, linetype = "solid"), 
           panel.background = element_blank())+
    #scale_fill_manual(values=colores) +
    labs(title = mainname, y = "phenotype", x = xname, fill = "gene type", face = "bold") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title=element_text(size=14,face="bold")) 
  return(g1)
  
}

ViolinplotStruc<-function(df, xname, mainname){
  dfaux <- df[df$tf_mutated_type != "inhibitors", ]
  dfaux$tf_mutated_type <- "all activators"
  df <- rbind(df, dfaux)
  
  g1<-ggplot(data = df, aes(x = type2, y=aux_var, fill=type2)) + 
    geom_violin(trim=T, alpha=1, bw=0.2)  + 
    #stat_boxplot(geom = "errorbar", width = 0.5) +
    geom_boxplot(width=0.15, fill="black", outlier.size = -1, notch=F)+
    stat_summary(fun.y=mean, geom="point", shape=18, size=3, col="white", na.rm=T)+
    stat_summary(fun.y=median, geom="point", shape=16, size=3, col="white", na.rm=T)+
    ylim(c(min(df$aux_var), ifelse(any(df$aux_var < 0),1, 1.3))) +
    # scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
    scale_fill_manual(values = kols)+
    facet_grid(.~tf_mutated_type)+
    theme(strip.text.y = element_text(size = 8, colour = "black", angle = 270, face = "bold")) + 
    theme(strip.text.x = element_text(size = 8, colour = "black", angle = 0, face = "bold")) + 
    theme(plot.title = element_text(size = rel(1.5), colour = "black", face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(legend.title = element_text(face = "bold")) +
    theme(axis.text.y = element_text(size = 15, colour = "black", angle = 0, face = "bold"))+ 
    theme(axis.text.x = element_text(margin = margin(t= 7),size = 14, colour = "black", angle = 45, face = "bold", hjust=1, vjust=1.08))+
    #guides(fill=FALSE) +
    theme(axis.title.x = element_text(size = 13, colour = "black", angle = 0, face = "bold"))+
    #scale_fill_brewer(palette ="Paired", values = c(0, 0.4, 0.6)) +
    theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black", 
                                                                                                            size = 1.5, linetype = "solid"), 
           panel.background = element_blank())+
    #scale_fill_manual(values=colores) +
    labs(title = mainname, y = "phenotype", x = xname, fill = "gene type", face = "bold") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title=element_text(size=14,face="bold"))
  return(g1)
  
}

filterPhenotype<-function(vec){
  if(PHENOTYPE_MODE == "proportion"){
    return(vec != 1)
  }else if(PHENOTYPE_MODE == "value"){
    return(vec != 0)
  }
}

makeAllPhenPlots<-function(df, cond){
  
  t0<- df[grep("terminal", df$type2), ]
  
  t1 <- t0[filterPhenotype(t0$tf_phenotype), ]
  t1$aux_var <- t1$tf_phenotype
  g1 <- densityplotStruc(t1, xname = "wt - TF mut", mainname = "Phenotype of TF mutation (with indirect)")
  g1box <- BoxplotStruc(t1, xname = "gene type", mainname = "Phenotype of TF mutation (with indirect)")
  g1vio <- ViolinplotStruc(t1, xname = "gene type", mainname = "Phenotype of TF mutation (with indirect)")
  
  t1prop <- t1
  t1prop$aux_var <- t1$tf_phenotype_prop
  #t1prop<-t1prop[t1prop$tf_mutated_type!= "inhibitors",]
  g1prop <- densityplotStruc(t1prop, xname = "TF mut/wt", mainname = "Phenotype of TF mutation (with indirect)")
  g1boxprop <- BoxplotStruc(t1prop, xname = "gene type", mainname = "Phenotype of TF mutation (with indirect)")
  g1vioprop <- ViolinplotStruc(t1prop, xname = "gene type", mainname = "Phenotype of TF mutation (with indirect)")
  
  
  
  t2 <- t0[filterPhenotype(t0$site_phenotype), ]
  t2$aux_var <- t2$tf_phenotype
  g2 <- densityplotStruc(t2, xname = "wt - TF mut", mainname = "Phenotype of TF mutation (only direct)")
  g2box <- BoxplotStruc(t2, xname = "gene type", mainname = "Phenotype of TF mutation (only direct)")
  g2vio <- ViolinplotStruc(t2, xname = "gene type", mainname = "Phenotype of TF mutation (only direct)")
  
  t2prop <- t2
  t2prop$aux_var <- t2$tf_phenotype_prop
  #t2prop<-t2prop[t2prop$tf_mutated_type!= "inhibitors",]
  g2prop <- densityplotStruc(t2prop, xname = "wt - TF mut", mainname = "Phenotype of TF mutation (only direct)")
  g2boxprop <- BoxplotStruc(t2prop, xname = "gene type", mainname = "Phenotype of TF mutation (only direct)")
  g2vioprop <- ViolinplotStruc(t2prop, xname = "gene type", mainname = "Phenotype of TF mutation (only direct)")
  
  
  
  t2$aux_var <- t2$site_phenotype
  g3 <- densityplotStruc(t2, xname = "wt - TFBS mut", mainname = "Phenotype of TFBS mutation")
  g3box <- BoxplotStruc(t2, xname = "wt - TFBS mut", mainname = "Phenotype of TFBS mutation")
  g3vio <- ViolinplotStruc(t2, xname = "wt - TFBS mut", mainname = "Phenotype of TFBS mutation")
  t2$aux_var <- t2$site_phenotype_prop
  g3prop <- densityplotStruc(t2, xname = "TFBS mut/wt", mainname = "Phenotype of TFBS mutation")
  g3boxprop <- BoxplotStruc(t2, xname = "TFBS mut/wt", mainname = "Phenotype of TFBS mutation")
  g3vioprop <- ViolinplotStruc(t2, xname = "TFBS mut/wt", mainname = "Phenotype of TFBS mutation")
  
  t2$aux_var <- t2$difference
  g7 <- densityplotStruc(t2, xname = "|TF phenotype| - |TFBS phenotype|", mainname = "Difference between TF and TFBS")
  g7box <- BoxplotStruc(t2, xname = "|TF phenotype| - |TFBS phenotype|", mainname = "Difference between TF and TFBS")
  g7vio <- ViolinplotStruc(t2, xname = "|TF phenotype| - |TFBS phenotype|", mainname = "Difference between TF and TFBS")
  
  
  
  t3 <- t2[t2$difference<0, ]
  t3$aux_var <- t3$difference
  tty <- table(t3$type2)
  g9 <- densityplotStruc(t3, xname = "|TF phenotype| - |TFBS phenotype|", mainname = "Positive differences between TF and TFBS")
  g9box <- BoxplotStruc(t3, xname = "|TF phenotype| - |TFBS phenotype|", mainname = "Positive differences between TF and TFBS")
  g9vio <- ViolinplotStruc(t3, xname = "|TF phenotype| - |TFBS phenotype|", mainname = "Positive differences between TF and TFBS")
  
  #### 
  t1<-t1[-grep("other", t1$type2), ]
  g4 <- densityplotStruc(t1, xname = "wt - TF mut", mainname = "Phenotype of TF mutation (with indirect)")  
  g4box <- BoxplotStruc(t1, xname = "gene type", mainname = "Phenotype of TF mutation (with indirect)")
  g4vio <- ViolinplotStruc(t1, xname = "gene type", mainname = "Phenotype of TF mutation (with indirect)")
  
  t1prop <- t1
  t1prop$aux_var <- t1$tf_phenotype_prop
  #t1prop<-t1prop[t1prop$tf_mutated_type!= "inhibitors",]
  g4prop <- densityplotStruc(t1prop, xname = "TF mut/wt", mainname = "Phenotype of TF mutation (with indirect)")
  g4boxprop <- BoxplotStruc(t1prop, xname = "gene type", mainname = "Phenotype of TF mutation (with indirect)")
  g4vioprop <- ViolinplotStruc(t1prop, xname = "gene type", mainname = "Phenotype of TF mutation (with indirect)")
  
  
  t2<-t2[-grep("other", t2$type2), ]
  t2$aux_var <- t2$tf_phenotype
  g5 <- densityplotStruc(t2, xname = "wt - TF mut", mainname = "Phenotype of TF mutation (only direct)")
  g5box <- BoxplotStruc(t2, xname = "gene type", mainname = "Phenotype of TF mutation (only direct)")
  g5vio <- ViolinplotStruc(t2, xname = "gene type", mainname = "Phenotype of TF mutation (only direct)")
  
  t2prop <- t2
  t2prop$aux_var <- t2$tf_phenotype_prop
  #t2prop<-t2prop[t2prop$tf_mutated_type!= "inhibitors",]
  g5prop <- densityplotStruc(t2prop, xname = "wt - TF mut", mainname = "Phenotype of TF mutation (only direct)")
  g5boxprop <- BoxplotStruc(t2prop, xname = "gene type", mainname = "Phenotype of TF mutation (only direct)")
  g5vioprop <- ViolinplotStruc(t2prop, xname = "gene type", mainname = "Phenotype of TF mutation (only direct)")
  
  t2$aux_var <- t2$site_phenotype
  g6 <- densityplotStruc(t2, xname = "wt - TFBS mut", mainname = "Phenotype of TFBS mutation")
  g6box <- BoxplotStruc(t2, xname = "wt - TFBS mut", mainname = "Phenotype of TFBS mutation")
  g6vio <- ViolinplotStruc(t2, xname = "wt - TFBS mut", mainname = "Phenotype of TFBS mutation")
  t2$aux_var <- t2$site_phenotype_prop
  g6prop <- densityplotStruc(t2, xname = "TFBS mut/wt", mainname = "Phenotype of TFBS mutation")
  g6boxprop <- BoxplotStruc(t2, xname = "TFBS mut/wt", mainname = "Phenotype of TFBS mutation")
  g6vioprop <- ViolinplotStruc(t2, xname = "TFBS mut/wt", mainname = "Phenotype of TFBS mutation")
  
  t2$aux_var <- t2$difference
  g8 <- densityplotStruc(t2, xname = "|TF phenotype| - |TFBS phenotype|", mainname = "Difference between TF and TFBS")
  g8box <- BoxplotStruc(t2, xname = "|TF phenotype| - |TFBS phenotype|", mainname = "Difference between TF and TFBS")
  g8vio <- ViolinplotStruc(t2, xname = "|TF phenotype| - |TFBS phenotype|", mainname = "Difference between TF and TFBS")
  
  if(PHENOTYPE_MODE == "proportion"){
    t3 <- t2[t2$difference>1, ]
  }else{
    t3 <- t2[t2$difference<0, ]
  }
  t3$aux_var <- t3$difference
  tty2 <- table(t3$type2)
  g10 <- densityplotStruc(t3, xname = "|TF phenotype| - |TFBS phenotype|", mainname = "Positive differences between TF and TFBS")
  g10box <- BoxplotStruc(t3, xname = "|TF phenotype| - |TFBS phenotype|", mainname = "Positive differences between TF and TFBS")
  g10vio <- ViolinplotStruc(t3, xname = "|TF phenotype| - |TFBS phenotype|", mainname = "Positive differences between TF and TFBS")
  
  gglist <- list(g1, g1box, g1vio, g1prop, g1boxprop, g1vioprop,
                 g2, g2box, g2vio, g2prop, g2boxprop, g2vioprop,
                 g3, g3box, g3vio, g3prop, g3boxprop, g3vioprop,
                 g4, g4box, g4vio, g4prop, g4boxprop, g4vioprop,
                 g5, g5box, g5vio, g5prop, g5boxprop, g5vioprop,
                 g6, g6box, g6vio, g6prop, g6boxprop, g6vioprop,
                 g7, g7box, g7vio,
                 g8, g8box, g8vio,
                 g9, g9box, g9vio,
                 g10, g10box, g10vio
  )
  pdfname<- paste(cond, "mutPhenDistr_all.pdf", sep="_", collapse="_")
  pdf(pdfname)
  for(gg in gglist) grid.arrange(gg,nrow = 1)
  grid.arrange(gg,nrow = 1)
  grid.table(data.frame(tty2))
  dev.off()
  
}



getCountTable<-function(df){
  df$regulator <- ifelse(df$site_phenotype != 0, 1, 0)
  df$case <- paste(as.character(df$gene), df$cell, df$simname, sep="_")
  dfaux <- df[df$tf_mutated_type != "inhibitors", ]
  dfaux$tf_mutated_type <- "all activators"
  df <- rbind(df, dfaux)
  
  tab <- aggregate(df$regulator, by=list(df$gene, df$cell, df$simname, df$tf_mutated_type), FUN=sum)
  
  
  names(tab) <- c("gene", "cell", "simname", "tf_mutated_type", "regulators")
  typeIndex <- unique(df[, c("gene", "type2", "cell")])
  names(typeIndex)[2] <- "type2"
  tab2 <- merge(tab, typeIndex, by.x = c("gene", "cell"), by.y = c("gene", "cell"))
  tab2$aux <- 1
  tab3 <- aggregate(tab2$aux,by=list(tab2$type2, tab2$tf_mutated_type, tab2$regulators), FUN=sum)
  names(tab3) <- c("type2", "tf_mutated_type", "regulators", "count")
  
  tab4 <- tab3
  auxtab <- unique(tab4[, c(1,2)])
  newtab <- data.frame()
  for(i in 1:nrow(auxtab)){ #now proportions
    a <- tab4[tab4$type2== auxtab$type2[i] & tab4$tf_mutated_type== auxtab$tf_mutated_type[i], ]
    a$percent_zero <- 100*a$count/sum(a$count)
    a$percent <- 100*a$count/sum(a$count[a$regulators !=0])
    a$percent[a$regulators==0] <- NA
    newtab <- rbind(newtab, a)
  }
  
  return(newtab)
}


histogramBase<-function(mu){
  g1 <- ggplot(data = mu, aes(x = regulators, y=percent, fill=type2, group=type2, linetype=type2)) + 
    #geom_histogram(position="dodge", binwidth = 1, center=0, bins = 0:10) +
    #geom_bar(aes(y = 100*(..count..)/sum(..count..)), position="dodge") +
    geom_bar(stat = "identity", position=position_dodge2(width = 0.9, preserve = "single" )) + #position="dodge" #, preserve = "single" for bar size inside position_dodge2
    scale_x_continuous(breaks=seq(0,max(mu$regulators),1))+
    # geom_vline(data=mu, aes(xintercept="mean_number_regulators", color=type2),
    #             linetype="dashed")+
    facet_grid( tf_mutated_type~ .)+
    scale_fill_manual(values = kols)+
    # scale_y_continuous(breaks=seq(0,100,by=20), limits=c(0,100)) +
    theme(plot.title = element_text(size = rel(1.5), colour = "black", face = "bold")) +
    theme(strip.text.y = element_text(size = 8, colour = "black", angle = 270, face = "bold")) + 
    theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0, face = "bold")) + 
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
  
  return(g1)
}

makeRegnumHistograms<-function(df, cond, withZeros=FALSE){
  counts<- getCountTable(df)
  if(withZeros)counts$percent <- counts$percent_zero
  g1<- histogramBase(counts[grepl("terminal", counts$type2), ])
  g2<- histogramBase(counts[!grepl("other", counts$type2) & grepl("terminal", counts$type2), ])
  g3<- histogramBase(counts[grepl("other", counts$type2) & grepl("terminal", counts$type2), ])  
  g4<- histogramBase(counts[!grepl("other", counts$type2) & grepl("terminal", counts$type2) & counts$tf_mutated_type %in% c("inhibitors", "all activators"), ])
  g5<- histogramBase(counts[grepl("other", counts$type2) & grepl("terminal", counts$type2) & counts$tf_mutated_type %in% c("inhibitors", "all activators"), ])  
  
  if(withZeros){
    pdfname<- paste(cond, "regByGene_Zeros.pdf", sep="_", collapse="_")
  }else{
    pdfname<- paste(cond, "regByGene.pdf", sep="_", collapse="_")
  }
  pdf(pdfname)
  grid.arrange(g1, nrow = 1)
  grid.arrange(g2, nrow = 1)
  grid.arrange(g3, nrow = 1)
  grid.arrange(g4, nrow = 1)
  grid.arrange(g5, nrow = 1)
  dev.off()
  
  
}
#################################################################################################
## functions to get corregulation
getCorregulationGlobal<-function(df, ncores=6){
  tfs_2use <- c("all activators", "inhibitors")
  df<-df[grep("terminal", df$type2),]
  #df <- df[df$site_phenotype !=0, ]
  df_aux <- df[df$tf_mutated_type != "inhibitors", ]
  df_aux$tf_mutated_type <- "all activators"
  df <- rbind(df, df_aux)
  
  genetypes <- df[, c("gene", "type")]
  genetypes<-genetypes[!duplicated(genetypes),]
  genetypes <- genetypes[order(genetypes$gene), ] %>% data.table
  setkey(genetypes, type)
  
  cnames <- unique(df$type)
  cnames <- cnames[order(sapply(cnames, nchar), cnames)]
  
  cl <- makeCluster(getOption("cl.cores", ncores))
  x<-clusterEvalQ(cl, library("magrittr"))
  x<-clusterEvalQ(cl, library("data.table"))
  
  newdf <- data.table()
  for(t in tfs_2use){
    subtab1<- df[df$tf_mutated_type==t, ]
    simlist<-list()
    for(sim in unique(df$simname)) simlist[[sim]] <- subtab1[subtab1$simname==sim, ] ##separate by simulations
    clusterExport(cl=cl, varlist=c("genetypes", "cnames", "t", "sim", "getCoregMat", "listOfNumeric"), envir= environment())
    
    allmats<-parLapply(cl=cl, X=simlist, fun=function(simsd){
      coregMat <- getCoregMat(simsd, genetypes, cnames)
      coregMat[, tf_type := t]
      coregMat[, sim := sim]
      return(coregMat)
    }) %>% rbindlist
    newdf <- rbind(newdf, allmats)
  }
  newdf2 <- aggregate(as.data.frame(newdf)[, cnames], by=list(newdf$tf_type, newdf$type), mean, na.rm=T)
  names(newdf2)[c(1,2)] <- c("tf_type", "type")
  return(newdf2)
}


getCorregulationByCell<-function(df, ncores=8){
  tfs_2use <- c("all activators", "inhibitors")
  df<-df[grep("terminal", df$type2),]
  #df <- df[df$site_phenotype !=0, ]
  df_aux <- df[df$tf_mutated_type != "inhibitors", ]
  df_aux$tf_mutated_type <- "all activators"
  df <- rbind(df, df_aux)
  
  genetypes <- df[, c("gene", "type")]
  genetypes<-genetypes[!duplicated(genetypes),]
  genetypes <- genetypes[order(genetypes$gene), ] %>% data.table
  setkey(genetypes, type)
  
  cnames <- unique(df$type)
  cnames <- cnames[order(sapply(cnames, nchar), cnames)]
  
  cl <- makeCluster(getOption("cl.cores", ncores))
  x<-clusterEvalQ(cl, library("magrittr"))
  x<-clusterEvalQ(cl, library("data.table"))
  
  newdf <- data.table()
  for(t in tfs_2use){
    subtab1<- df[df$tf_mutated_type==t, ]
    simlist<-list()
    for(sim in unique(df$simname)) simlist[[sim]] <- subtab1[subtab1$simname==sim, ] ##separate by simulations
    clusterExport(cl=cl, varlist=c("genetypes", "cnames", "t", "sim", "getCoregAllcells", "getCoregMatDifCells", "listOfNumeric", "getCoregMatDifCellsIntersect"), envir= environment())
    
    allmats<-parLapply(cl=cl, X=simlist, fun=function(simsd){
      coregMat <- getCoregAllcells(simsd, genetypes, cnames)
      coregMat[, tf_type := t]
      coregMat[, sim := sim]
      return(coregMat)
    }) %>% rbindlist
    newdf <- rbind(newdf, allmats)
  }
  newdf2 <- aggregate(as.data.frame(newdf)[, grep("cell[0-9]", names(newdf))], by=list(newdf$tf_type, newdf$type, newdf$cell), mean, na.rm=T)
  names(newdf2)[1:3] <- c("tf_type", "type", "cell")
  return(newdf2)
}

listOfLists<-function(n=10){
  l <- list()
  for(i in 1:n){
    l[[i]] <- list()
  }
  return(l)
}

listOfNumeric<-function(names){
  l <- list()
  for(n in names){
    l[[n]]<-list()
    for(n2 in names){
      l[[n]][[n2]] <- numeric(0)
    }}
  return(l)
}

getCoregMat<- function(st, genet, cnames){
  st2 <- st[st$site_phenotype != 0, ]
  genet[, reg:=sapply(genet$gene, FUN=function(x, st2){
    res<-st2$mutation[st2$gene==x]
    if(length(res)>0){
      return(unique(res))
    }else{
      return (numeric(0))
    }
  }, st2)]
  mat <- matrix(numeric(nrow(genet)^2), ncol=nrow(genet) ) %>% set_colnames(genet$gene) %>% set_rownames(genet$gene)
  typem <- listOfNumeric(cnames) #stores correlation between each pair of genes
  #first, gene to gene
  for(i in 2:nrow(mat)){
    mat[i, i] <- 1
    for(j in 1:(i-1)){
      common <- intersect(genet[i, reg][[1]], genet[j, reg][[1]])
      allregs <- union(genet[i, reg][[1]], genet[j, reg][[1]])
      mat[i, j] <- length(common)/length(allregs)
      mat[j, i] <- mat[i, j]
      typei <- genet[i, type]
      typej <- genet[j, type]
      typem[[typei]][[typej]] <- c(typem[[typei]][[typej]], mat[i, j])
      typem[[typej]][[typei]] <- c(typem[[typej]][[typei]], mat[i, j])
      
    }
  }
  means <- lapply(typem, FUN=function(x)sapply(x, mean) %>% t %>% as.data.frame) %>% rbindlist
  means[,type:=names(typem)]
  return(means)
}

getCoregMatDifCells<- function(st, st2, genet, cnames){
  
  mat <- matrix(numeric(nrow(genet)^2), ncol=nrow(genet) ) %>% set_colnames(genet$gene) %>% set_rownames(genet$gene)
  typem <- listOfNumeric(cnames) #stores correlation between each pair of genes
  #first, gene to gene
  
  for(i in 2:nrow(mat)){
    aux1 <- st[st$gene==genet[i, gene], ]
    aux1 <- aux1[order(aux1$mutation),]
    for(j in 1:(i-1)){
      aux2 <- st2[st2$gene==genet[j, gene], ]
      aux2 <- aux2[order(aux2$mutation),]
      if(any(aux1$site_phenotype != 0) & any(aux2$site_phenotype != 0) ){
        #nonzero <- (aux1$site_phenotype != 0) | (aux2$site_phenotype != 0)
        #corr <- cor(aux1$site_phenotype[nonzero], aux2$site_phenotype[nonzero])
        corr <- cor(aux1$site_phenotype, aux2$site_phenotype, method="pearson")
      }else{
        corr <- 0
      }
      mat[i, j] <- corr
      mat[j, i] <- corr
      typei <- genet[i, type]
      typej <- genet[j, type]
      typem[[typei]][[typej]] <- c(typem[[typei]][[typej]], corr)
      typem[[typej]][[typei]] <- c(typem[[typej]][[typei]], corr)
    }
  }
  means <- lapply(typem, FUN=function(x)sapply(x, mean) %>% t %>% as.data.frame) %>% rbindlist
  #means[,type:=names(typem)]
  return(means)
}

getCoregMatDifCellsIntersect<- function(st, st2, genet, cnames, threshold=0){
  st <- st[abs(st$site_phenotype) > threshold, ]
  st2 <- st2[abs(st2$site_phenotype) > threshold, ]
  auxfun <- function(x, st){
    res<-st$mutation[st$gene==x]
    if(length(res)>0){
      return(unique(res))
    }else{
      return (numeric(0))
    }}
  genet[, reg1:=sapply(genet$gene, FUN=auxfun, st)]
  genet[, reg2:=sapply(genet$gene, FUN=auxfun, st2)]
  
  mat <- matrix(numeric(nrow(genet)^2), ncol=nrow(genet) ) %>% set_colnames(genet$gene) %>% set_rownames(genet$gene)
  typem <- listOfNumeric(cnames) #stores correlation between each pair of genes
  #first, gene to gene
  
  for(i in 2:nrow(mat)){
    for(j in 1:(i-1)){
      common <- intersect(genet[i, reg1][[1]], genet[j, reg2][[1]])
      allregs <- union(genet[i, reg1][[1]], genet[j, reg2][[1]])
      if(length(allregs) == 0){
        corr <- 0
      }else{
        corr <- length(common)/length(allregs)
      }
      mat[i, j] <- corr
      mat[j, i] <- corr
      typei <- genet[i, type]
      typej <- genet[j, type]
      typem[[typei]][[typej]] <- c(typem[[typei]][[typej]], corr)
      typem[[typej]][[typei]] <- c(typem[[typej]][[typei]], corr)
    }
  }
  means <- lapply(typem, FUN=function(x)sapply(x, mean) %>% t %>% as.data.frame) %>% rbindlist
  #means[,type:=names(typem)]
  return(means)
}

getCoregAllcells<- function(simdf, genet, cnames){
  cells <- unique(simdf$cell)
  alldf<-data.table()
  for(c1 in cells){
    st <- simdf[simdf$cell==c1, ]
    cell1df <- data.table(cell=c1, type=cnames)
    for(c2 in cells){
      st2 <- simdf[simdf$cell==c2, ]
      means<- getCoregMatDifCellsIntersect(st, st2, genet, cnames)
      names(means) <- paste(gsub("exp", "cell", c2), names(means), sep="_")
      cell1df<-cbind(cell1df, means)
    }
    alldf <- rbind(alldf, cell1df) 
  }
  return(alldf)
}


makeCorHeatmaps <- function(corr, cond="", save=T){
  for(t in unique(corr$tf_type)){
    mat <- as.matrix(corr[corr$tf_type==t, 3:ncol(corr)]) %>% set_rownames(corr$type[corr$tf_type==t])
    cnames <- unique(corr$type)
    cnames <- cnames[order(sapply(cnames, nchar), cnames)]
    mat<-mat[cnames, cnames]
    name<- paste(cond, " corregulation - all cells ", t, sep="", collapse="")
    fname <- paste(cond,t, "corrAllCells.pdf", sep="_", collapse="")
    if(save){
      pheatmap(mat, cluster_rows = F, cluster_cols = F, angle_col=45, main=name, fontsize=12, fontsize_col = 15, fontsize_row = 15, legend_labels=seq(0,1, 0.1),breaks=seq(0,1, 0.01), filename = fname)
    }else{
      pheatmap(mat, cluster_rows = F, cluster_cols = F, angle_col=45, main=name, fontsize=12, fontsize_col = 15, fontsize_row = 15, legend_labels=seq(0,1, 0.1),breaks=seq(0,1, 0.01))
    }
    
  }
  fname2 <- paste(cond, "corrAllCells.csv", sep="_", collapse="")
  if(save){
    write.table(corr,file = fname2, sep="\t", row.names=F, dec=".", quote=F)
  }
}
