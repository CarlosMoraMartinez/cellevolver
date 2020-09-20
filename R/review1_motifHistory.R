

library(tidyverse)
#Load simulation data
setwd("/home/carmoma/projects/cellevolver/cellevolver/history_sims/")

outdir = "/home/carmoma/projects/cellevolver/historyplots/"
if(!dir.exists(outdir)) dir.create(outdir)

annpattern <- "g[0-9]+_annotation.csv"
types <- c(rep("tf", 15), rep("terminal all", 5), rep("terminal 1 cell", 20))
tf_types <- c(rep("activator" ,11), rep("inhibitor", 4))
colorscheme <-c("Activator"="dodgerblue", 
                "Inhibitor"="firebrick3",
                "terminal this"="wheat3",
                "terminal 1 cell"="wheat3",
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

getFullAnnotationHistorySingleInd <- function(cond, sim){
  dir <- paste(cond, sim, sep="/", collapse="")
  files <- list.files(dir) %>% subset(grepl(annpattern, .)) 
  allann <- files %>%
    map(function(x)read_tsv(paste(dir, x, sep="/", collapse="")) %>% 
          mutate(generation=str_extract(x, "g[0-9]+_") %>% gsub("g|_", "", .) %>% as.numeric)
    ) %>%
    bind_rows
  return(allann)
}

getFullAnnotationHistory <- function(cond, sim){
  dir <- paste(cond, sim, sep="/", collapse="")
  files <- list.files(dir) %>% subset(grepl(annpattern, .)) 
  allann <- files %>%
          map(function(x)read_tsv(paste(dir, x, sep="/", collapse="")) %>% 
                mutate(generation=str_extract(x, "g[0-9]+o") %>% gsub("g|o", "", .) %>% as.numeric),
                mutate(organism=str_extract(x, "o[0-9]+_") %>% gsub("o|_", "", .) %>% as.numeric)
              ) %>%
  bind_rows
  return(allann)
}

getAnnotationSummary<-function(sim, cond){
  #allann <- getFullAnnotationHistorySingleInd(cond, sim)
  allann <- getFullAnnotationHistory(cond, sim)
  summbymot <- allann %>% unite(col="motifId", gene, positions, tf_ind, sep="_", remove=FALSE) %>% 
    group_by(motifId) %>% summarise(start=min(generation), 
                                    end=max(generation),
                                    starting_affinity = min(percent[generation==min(generation)]),
                                    final_affinity = max(percent[generation==max(generation)]),
                                    ) %>% 
    mutate(duration = end - start + 1, 
           remained = end == max(end),
           sim_length = max(end)) %>% 
    separate(motifId, into=c("gene", "positions", "tf_ind"), remove=FALSE, convert=TRUE) %>% 
    mutate(gene_type = types[gene + 1], 
         tf_type = tf_types[tf_ind + 1],
         cond = cond,
         sim = str_extract(sim, "s[0-9]+$")
         ) %>% 
    unite("globalMotifId", sim, motifId, remove=FALSE, sep="_")
  return(summbymot)
}

getMeansTable <- function(d, discard_generations = -1){
  d3 <- d %>% filter(start >= discard_generations) %>% 
    group_by(gene_type, tf_type) %>% summarise(num_motifs = n(),
                                             num_motifs_per_gene = n()/length(unique(gene)),
                                             mean_start = mean(start),
                                             median_start = median(start),
                                             relative_mean_start = mean_start/max(end),
                                             relative_median_start = median_start/max(end),
                                             mean_life = mean(duration),
                                             median_life = median(duration),
                                             sd_life = sd(duration),
                                             var_life = var(duration),
                                             num_persisted = length(which(remained)),
                                             prob_persist = length(which(remained))/num_motifs
  )
  return(d3)
}

conditions <- list.files()
cond <- conditions[1]
sims <- list.files(cond) %>% subset(!grepl("\\.", .))
sim <- sims[1]


d <- sims %>% map(getAnnotationSummary, cond) %>% bind_rows
d2 <- d %>% filter(gene_type != "tf", start > 250) 
d3 <- getMeansTable(d)
d3 <- getMeansTable(d, 200)
d3 <- getMeansTable(d, 350)


(g1<-ggplot(d, aes(x=start, color=gene_type, fill=gene_type)) + 
    geom_histogram(aes(y=..count../sum(..count..)*100), position = "dodge") +
    #geom_histogram(position = "dodge") +
    #geom_bar(stat = "identity", alpha=0.8)  +
    facet_grid(gene_type~tf_type) + #scales = "free_x", space = "free"
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


par(mfrow=c(1,2))
aux <- allann %>% filter(gene==15)
plot(aux$generation, aux$positions, col=aux$tf_ind)
aux2 <- allann %>% filter(gene==25)
plot(aux2$generation, aux2$positions, col=aux2$tf_ind)

