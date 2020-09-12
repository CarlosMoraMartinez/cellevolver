library(tidyverse)

setwd("/home/carmoma/projects/cellevolver/data/all_tables")
outdir = "/home/carmoma/projects/cellevolver/augustplots/"
if(!dir.exists(outdir)) dir.create(outdir)


datasets <- "4cell_mce0  4cell_mce0fix  4cell_mce0fix_mutb  4cell_mce1fix  4cell_mce2fix  4cell_mce2inh5" %>% 
    strsplit("\\s+") %>% unlist
colorscheme <-c("Activator"="dodgerblue", 
                "Inhibitor"="firebrick3",
                "terminal this"="wheat3",
                "terminal other"="gray",
                "2 cells this" = "cyan2",
                "2 cells other" = "darkolivegreen3",
                "terminal all" = "hotpink",
                "tf" = "blue",
                "non tf" = "yellow",
                "lin" = "green", 
                "lineage this cell" = "limegreen",
                "lineage other cell" = "purple3",
                "other activators" = "dodgeblue1"
                )


tf_pattern <- "tfsPhenotype"
tfbs_pattern <- "sitesPhenotype"
terminalPlotPattern <- "TerminalJoinAndExpr_cell0.csv"
reprPattern <- "-1lin[0]andExpr.csv"

expr_thresh <- c(0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2)


#Store names in variables
#terminal_exprCols <- 0:3 %>% as.character
#terminal_genegroups <- list()
#start=4
#for(i in 1:5){terminal_genegroups[[i]] = start:(start+4) %>% as.character;start = start+5}
newnames <- c("tf", 
              paste("expr", as.character(1:4), sep="_"),  
              paste("A", as.character(1:5), sep="_"),
              paste("B", as.character(1:5), sep="_"),
              paste("C", as.character(1:5), sep="_"),
              paste("D", as.character(1:5), sep="_"),
              paste("E", as.character(1:5), sep="_"))
newnames2 <- c("tf", 
              paste("expr", as.character(1:4), sep="_"),  
              paste("A", as.character(1:5), sep="_"),
              paste("B", as.character(1:5), sep="_"),
              paste("C", as.character(1:5), sep="_"),
              paste("D", as.character(1:5), sep="_"),
              paste("D2", as.character(1:5), sep="_"),
              paste("D3", as.character(1:5), sep="_"),
              paste("E", as.character(1:5), sep="_"))

readHeatmapFile <- function(file){
  raw <- read_csv(file)
  if(ncol(raw) == length(newnames)){names(raw) <- newnames}else{names(raw) <- newnames2}
  #vars2plot <- c()
  expmat <- raw %>% select(matches("expr_[0-9]", perl=T)) %>% as.matrix 
  for (th in expr_thresh){
    vname = paste("thresh", as.character(th), sep="_", collapse="_")
    #vars2plot <- c(vars2plot, vname)
    raw[[vname]] <- apply(expmat > th, MARGIN = 1, sum)
  }
  d <- raw %>% gather(key="gene_type", value="phenotype", matches("[A-E][23]?_[1-5]", perl=T))
  d$phenotype_group <- sapply(d$phenotype, function(x){
    #if(x < -3){
    #  r="< -3"
    #}else 
    if(x < -2){
        r="< -2"
    }else if(x < -1){
      r="[-2, -1["
    }else if(x < -0.5){
      r="[-1, -0.5["
    }else if(x <= 0){
      r="[-0.5, 0]"
    }else{
        r="> 0"
    }
  }) %>% 
    factor(levels=c("< -2", "[-2, -1[", "[-1, -0.5[", "[-0.5, 0]", "> 0")) #"[-3, -2["
  d$expr_group <- ifelse(grepl("E_", d$gene_type), 
                         "terminal all", 
                         ifelse(grepl("A_", d$gene_type), 
                                "terminal this", 
                                ifelse(grepl("D2_", d$gene_type),
                                       "2 cells this", 
                                       ifelse(grepl("D3_", d$gene_type),
                                              "2 cells other",
                                "terminal other"))))
  return(d)
}

getPercentTable <- function(d){
  d2 <- data.frame()
  for(gg in unique(d$active_in_cells)) for(eg in unique(d$expr_group)){
    aux1 <- d[d$active_in_cells == gg & d$expr_group == eg, ]
    aux2 <- aux1 %>% pull( phenotype_group) %>% table  %>% as.data.frame %>% set_names(c("phenotype_group", "n"))
    aux2$percentage <- 100*aux2$n/(sum(aux2$n))
    aux2$active_in_cells <- gg
    aux2$expr_group <- eg
    d2 <- rbind(d2, aux2)
  }
  return(d2)
}
getPercentTable2 <- function(d){
  d2 <- data.frame()
  for(eg in unique(d$expr_group)){
    aux <- d[d$expr_group == eg, ] %>% 
      count(active_in_cells, phenotype_group) %>% 
      mutate(percentage = 100*n/sum(n),
             expr_group = eg)
    d2 <- rbind(d2, aux)
  }
 # d2 <- d2 %>% mutate(phenotype_group = factor(phenotype_group, 
  #              levels=c("< -3", "[-3, -2[", "[-2, -1", "[-1, -0]", "> 0")))
  return(d2)
}

files <- list.files() %>% subset(grepl(terminalPlotPattern,.))
#file <- files[1]

for (file in files){
  d <- readHeatmapFile(file)
  d <- d %>% filter(phenotype != 0)
  vars2plot <- names(d)[grep("thresh_", names(d))]
   for(var in vars2plot){
    d[["active_in_cells"]] <- d[[var]]
    d2 <- getPercentTable(d)
    d2 <- d2[d2$active_in_cells != 0, ]
    d2 <- d2 %>% mutate(expr_group = factor(expr_group, 
                                            levels=c("terminal this", "terminal other", 
                                                     "2 cells this", "2 cells other", "terminal all"))) %>% 
      filter(!grepl("other", expr_group))

    (g1<-ggplot(d2, aes(x=phenotype_group, y=percentage, color=expr_group, fill=expr_group)) + 
      #geom_histogram(aes(y=..count../sum(..count..)*100), breaks=c(-4,-3,-2,-1,0,1,2)) +
      geom_bar(stat = "identity", alpha=0.8)  +
      facet_grid(active_in_cells~expr_group) + #scales = "free_x", space = "free"
      #scale_y_continuous(breaks=seq(0,0.6, 0.1), labels = scales::percent(x=seq(0,0.6, 0.1))) + #labels = scales::percent(10), 
      scale_color_manual(values=colorscheme) + 
      scale_fill_manual(values=colorscheme) + 
      theme(axis.text.x = element_text(color = "#000000", size = 12, angle=90, debug=FALSE, vjust = 0.5, hjust=1 ),
            axis.title.x = element_text(size = 13, face="bold"),
            axis.text.y = element_text(color = "#000000", size = 15),
            axis.title.y = element_text(size = 13, face="bold"),
            plot.title = element_text(hjust = 0.5, size = 17, face="bold"),
            #legend.position = c(0.87, 0.9),
            #legend.key = element_blank(),
            #legend.text = element_text(face="bold", size=15),
            axis.line = element_line(size = 0.5, colour = "black"),
            #strip.background = element_blank(),
            strip.text.x = element_text(size=15, face = "bold"), #element_blank(),
            strip.text.y = element_text(size=15, face = "bold"),
            #panel.spacing.x=unit(0, "lines"),
            panel.background = element_rect(fill = "white")) +
      labs(y="Percentage of genes", x = "Phenotype (mutant - wildtype)") +
      ggtitle(gsub("_", "old = " , var))
    )
    dataset <- gsub("_cell0.csv|.csv", "", file)
    basename <- paste(outdir, dataset, "PhenByTFexpr1_", var, sep="", collapse="")
    ggsave(paste(basename, ".pdf", sep="", collapse=""), g1)
    write_tsv(d2, paste(basename, ".tsv", sep="", collapse=""))
    #Now with dodge
    d2 <- getPercentTable2(d)
    d2 <- d2[d2$active_in_cells != 0, ]
    d2 <- d2 %>% mutate(expr_group = factor(expr_group, levels=c("terminal this", "terminal other",
                                                                 "2 cells this", "2 cells other",  
                                                                 "terminal all")),
                        active_in_cells = factor(active_in_cells, levels = c(1, 2, 3, 4)))%>% 
      filter(!grepl("other", expr_group))
    (g1<-ggplot(d2, aes(x=phenotype_group, y=percentage, color=active_in_cells, fill=active_in_cells)) + 
        #geom_histogram(aes(y=..count../sum(..count..)*100), breaks=c(-4,-3,-2,-1,0,1,2)) +
        geom_bar(stat = "identity", alpha=0.8)  +
        facet_grid(.~expr_group) + #scales = "free_x", space = "free"
        #scale_y_continuous(breaks=seq(0,0.6, 0.1), labels = scales::percent(x=seq(0,0.6, 0.1))) + #labels = scales::percent(10), 
        #scale_color_manual(values=colorscheme) + 
        #scale_fill_manual(values=colorscheme) + 
        theme(axis.text.x = element_text(color = "#000000", size = 12, angle=90, debug=FALSE, vjust = 0.5, hjust=1 ),
              axis.title.x = element_text(size = 13, face="bold"),
              axis.text.y = element_text(color = "#000000", size = 15),
              axis.title.y = element_text(size = 13, face="bold"),
              plot.title = element_text(hjust = 0.5, size = 17, face="bold"),
              #legend.position = c(0.87, 0.9),
              #legend.key = element_blank(),
              #legend.text = element_text(face="bold", size=15),
              axis.line = element_line(size = 0.5, colour = "black"),
              #strip.background = element_blank(),
              strip.text.x = element_text(size=15, face = "bold"), #element_blank(),
              strip.text.y = element_text(size=15, face = "bold"),
              #panel.spacing.x=unit(0, "lines"),
              panel.background = element_rect(fill = "white")) +
        labs(y="Percentage of genes", x = "Phenotype (mutant - wildtype)") +
        ggtitle(gsub("_", "old = " , var))
    )
    dataset <- gsub("_cell0.csv|.csv", "", file)
    basename <- paste(outdir, dataset, "PhenByTFexpr2_", var, sep="", collapse="")
    ggsave(paste(basename, ".pdf", sep="", collapse=""), g1)
    write_tsv(d2, paste(basename, ".tsv", sep="", collapse=""))
  }
}


