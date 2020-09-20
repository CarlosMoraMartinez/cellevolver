library(tidyverse)

setwd("/home/carmoma/projects/cellevolver/data")
outdir = "/home/carmoma/projects/cellevolver/augustplots/"
if(!dir.exists(outdir)) dir.create(outdir)


datasets <- "4cell_mce0  4cell_mce0fix  4cell_mce0fix_mutb  4cell_mce1fix  4cell_mce2fix  4cell_mce2inh5" %>% 
    strsplit("\\s+") %>% unlist
colorscheme <-c("Activator"="dodgerblue", "Inhibitor"="red2")

expr_pattern <- "_finalExpression.csv"
tf_pattern <- "_mutantTFs.csv"
tfbs_pattern <- "_mutantSites.csv"
timexpr_pattern <- "_timeExpression.csv"
tfs_pattern <- "_TFset.csv"
annot_pattern <- "annotation.csv"
  
expr_thresh <- c(0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2)


makeExpressionBarplots<-function(dataset){
  expr_raw <- list.files(recursive = TRUE) %>% subset(grepl(pattern=expr_pattern, .)) %>% 
    map(function(x){
      tryCatch(read_tsv(x ) %>% mutate(seqs=NULL, gene = X1, X1 = NULL), error = function(e)cat("ERROR reading ", x)
      )
      }) %>% bind_rows
  expr <- expr_raw %>% filter(type %in% c("1", "-1")) 
  expmat <- expr %>% select(matches("exp[0-9]", perl=T)) %>% as.matrix 
  vars2plot <- c()
  
  for (th in expr_thresh){
    vname = paste("thresh", as.character(th), sep="_", collapse="_")
    vars2plot <- c(vars2plot, vname)
    expr[[vname]] <- apply(expmat > th, MARGIN = 1, sum)
  }
  
  expr$type <- recode(expr$type, "1" = "Activator", "-1" = "Inhibitor")
  expr$barcolor <- colorscheme[expr$type]
  
  for(var in vars2plot){
    expr[["active_in_cells"]] <- expr[[var]]
    expr_summary <- expr %>% count(type, active_in_cells) %>% group_by(type) %>% 
      mutate(prop = n / sum(n))
    g1<-ggplot(expr_summary, aes(x=active_in_cells, y=prop, color=type, fill=type)) + 
      geom_bar(stat = "identity", alpha=0.6, position="dodge") +
      facet_grid(~type, scales = "free_x", space = "free") +
      scale_y_continuous(breaks=seq(0,0.8, 0.1), labels = scales::percent(x=seq(0,0.8, 0.1))) + #labels = scales::percent(10), 
      scale_color_manual(values=colorscheme) + 
      scale_fill_manual(values=colorscheme) + 
      theme(axis.text.x = element_text(color = "#000000", size = 15, angle=0, debug=FALSE, vjust = 1.0, hjust=1 ),
            axis.title.x = element_text(size = 13, face="bold"),
            axis.text.y = element_text(color = "#000000", size = 15),
            axis.title.y = element_text(size = 13, face="bold"),
            plot.title = element_text(hjust = 0.5, size = 17, face="bold"),
            legend.position = c(0.87, 0.9),
            legend.key = element_blank(),
            legend.text = element_text(face="bold", size=15),
            axis.line = element_line(size = 0.5, colour = "black"),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            panel.spacing.x=unit(0, "lines"),
            panel.background = element_rect(fill = "white")) +
      labs(y="Percentage of TFs", x = "Number of cells in which TFs are expressed") +
      ggtitle(gsub("_", "old = " , var))
    
    basename <- paste(outdir, dataset, "/", "expr_", var, sep="", collapse="")
    ggsave(paste(basename, ".pdf", sep="", collapse=""), g1)
    write_tsv(expr_summary, paste(basename, ".tsv", sep="", collapse=""))
  }
}

mutType=tf_pattern
#plotMutants<-function(mutType=tf_pattern){
  expr_raw <- list.files(recursive = TRUE) %>% subset(grepl(pattern=mutType, .)) %>% .[1:10] %>% 
    map(function(x){
      tryCatch(read_tsv(x ) %>% mutate(seqs=NULL, gene = X1, X1 = NULL, fname=x), error = function(e)cat("ERROR reading ", x)
      )
    }) %>% bind_rows 
  types <- expr_raw %>% select(gene, type) %>% distinct %>% as.data.frame
  expr_raw$tf_type <-  sapply(expr_raw$tf_mutated, function(x) types$type[types$gene==x])
  
#}

for(dataset in datasets){
  if(!dir.exists(paste(outdir, dataset, sep="", collapse=""))) dir.create(paste(outdir, dataset, sep="", collapse=""))
  setwd(dataset)
  makeExpressionBarplots(dataset)
  cat(dataset)
  setwd("../")
}


