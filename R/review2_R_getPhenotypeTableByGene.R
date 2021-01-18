
####################################################################################################

source("~/projects/cellevolver/cellevolver/R/201203_functions.R")

#map types to numbers and colors
type1define = list('3'='lineage this cell','4'='lineage other cell', '8'='other activators', 
                   '9'='inhibitors', '10'='inhibitor lineage this cell',
                   '11'='inhibitor lineage other cell',
                   '5'='lineage 3 this', '6'='lineage 3 other','7'='lineage all', 
                   '12'='terminal specific this', '15'='terminal specific other',
                   '13'='terminal 2 this', '16'='terminal 2 other', 
                   '17'='terminal 3 this', '18'='terminal 3 other', 
                   '14'='terminal all', '1'='tf', '2'='non tf', '0'='lin')
#col2map <- data.frame(t=unlist(type1define), col = c("limegreen", "purple3", "dodgerblue1", "firebrick3", "red", "red", "chartreuse3" ,"chartreuse3", "orange", "burlywood4", "gray", "cyan2", "darkolivegreen3", "hotpink", "blue", "yellow", "green"), stringsAsFactors = F)
col2map <- data.frame(t=unlist(type1define), col = c("limegreen", "purple3", "dodgerblue1", "firebrick3", "red", "red","cyan2", "darkolivegreen3", "orange" ,"wheat3", "gray", "cyan2", "darkolivegreen3","cyan2", "darkolivegreen3", "hotpink", "blue", "yellow", "green"), stringsAsFactors = F)


kols=c()
for(i in 1:nrow(col2map))kols[col2map$t[i]]=col2map$col[i]



#condlist = c('4cell_mce0', '4cell_mce0fix_mutb', '4cell_mce1fix', '4cell_mce2fix', '4cell_mce2inh5', '4cell_mce0fix', '4cell_mce0_mutb', '4cell_mce1fix_mutb', '4cell_mce2fix_mutb', '4cell_mce2inh5_mutb', '4cell_mce0Xss')
#condlist <- c('4cell_mce0Xss', '4cell_mce1inh5Xss', '4cell_mce2inh5Xss')
#setwd("C:/Users/Carlos/Desktop/cellEvolver_all/1802_cellevolver/all_finaldata/motifs_sep1_sepmot3_selby0_filtTrue_ignoreSelfReg_fixed1/")

#condlist = c('4cell_mce0fix')
#setwd("C:/Users/Carlos/Desktop/cellEvolver_all/1802_cellevolver/all_finaldata/180911_rmInhibitors/180914_motifsPartial1")


#setwd("C:/Users/Carlos/Desktop/cellEvolver_all/1802_cellevolver/all_finaldata/181101set_plusAll/sim_tables/")
setwd("/home/carmoma/projects/cellevolver/data/all_tables")
dirs <- list.dirs()
condlist <- dirs[grepl("cell_mce[0-9]", dirs) & !grepl("error", dirs)] %>% gsub(pattern="\\./", replacement="", x=.)
condlist <- condlist[condlist!= "4cell_mce9fix"] #this condition has a mistake

newdir <-"./201223_phenotypes_by_gene"
if(!file.exists(newdir)) dir.create(newdir)
setwd(newdir)

#sourcedir <- "C:/Users/Carlos/Desktop/cellEvolver_all/1802_cellevolver/all_finaldata/181101set_plusAll/sim_tables/"
#sourcedir <- "C:/Users/Carlos/Desktop/cellEvolver_all/1802_cellevolver/all_finaldata/cellevolver_chromatin/simulations/all_tables_chromatin/"
#sourcedir <- "C:/Users/Carlos/Desktop/cellEvolver_all/1802_cellevolver/all_finaldata/180911_rmInhibitors/all_tables_rmInh/"
sourcedir <- "/home/carmoma/projects/cellevolver/data/all_tables/"

sitesname <- "mutantSites.csv"
tfname <- "mutantTFs.csv"
exprname <- "finalExpression.csv"

PHENOTYPE_MODE <- "value"
calc_mode = c("diff", "proportion")
#getType2 fails
for(cond in condlist){
  ### 1) get files in condition, ordered by simulation
  sdir <- paste(sourcedir, cond, sep="", collapse="")
  f <- list.files(sdir)
  f <- f[grep(".csv$", f)]
  if(length(f)==0) next
  sim <- sapply(f, FUN=function(x){
    a<-strsplit(x, "_")[[1]]
    b<- paste(a[1:(length(a)-1)], sep="_", collapse="_")
    b <- gsub("_$", "", b)
    return(b)
    })
  allfiles <- data.frame(file=f, sim=sim)
  allfiles$fullname <- paste(sdir, allfiles$file, sep="/")
  cat(cond, length(unique(sim)), "\n")
  
  df <- makeBigPhenotypeTable(allfiles, calc_mode)
  if(nrow(df)==0) next
  
  typeIndex <- unique(df[, c("gene", "type2", "cell")])
  names(typeIndex)[2] <- "tf_mutated_type"
  df2 <- merge(df, typeIndex, by.x = c("mutation", "cell"), by.y = c("gene", "cell"))
  #names(typeIndex)[2] <- "tf_mutated_type"
  df2$condition <- cond
  write.table(df2, file =paste(cond, "allPhenotypeByGene.csv", sep="_", collapse="_"), sep="\t", row.names=FALSE, quote=FALSE)
  
  makeRegnumHistograms(df2, cond, TRUE) #Last argument: with zeros
  makeAllPhenPlots(df2, cond)
  cat("\n******************************\n", cond, "\n*****************************\n")
}



### If files have been already generated

setwd("C:/Users/Carlos/Desktop/cellEvolver_all/1802_cellevolver/all_finaldata/motifs_sep1_sepmot3_selby0_filtTrue_ignoreSelfReg_fixed1/180829phenotypes_by_gene_good/")
tffiles <- list.files()
tffiles <- tffiles[grep(".csv$", tffiles)]
tffiles <- tffiles[-grep("corrAllC", tffiles)]
tffiles <- tffiles[tffiles != "4cell_mce9fix_allPhenotypeByGene.csv"]

for(bf in tffiles){
 # df <- read.table(bf, header=T, stringsAsFactors = F, sep="\t")
  ## THIS REMOVES TERMINAL OTHER THAT ARE NOT EXPRESSED
  df <- df[df$wt_expression>0.01, ]
  cond <- strsplit(bf, "_")[[1]]
  cond <- paste(cond[1:(length(cond)-1)], sep="_", collapse="_")
  makeAllPhenPlots(df, cond)
  makeRegnumHistograms(df, cond, TRUE)
  cat("\n******************************\n", cond, "\n*****************************\n")
}

## falta por celula con mce7
for(bf in tffiles){
  df <- read.table(bf, header=T, stringsAsFactors = F, sep="\t")
  cond <- strsplit(bf, "_")[[1]]
  cond <- paste(cond[1:(length(cond)-1)], sep="_", collapse="_")
  #corr<-getCorregulationGlobal(df) #without separating cells

  corr2 <- getCorregulationByCell(df)
  fname <- paste(cond, "correlByCell_IntersectThres0.csv", sep="_", collapse="_")
  write.table(corr2,file = fname, sep="\t", dec=".", row.names = F, quote=F)
  #makeCorHeatmaps(corr, cond)
  cat("\n", cond, "\n")
}

################################################
################################################
################################################
################################################
################################################


#Finally, get intersection, correlation etc

terp <- df[df$type2 %in% c("terminal specific this", "terminal all"), ]
terp <- terp[terp$tf_mutated_type != "inhibitors", ]
terp <- terp[terp$site_phenotype !=0, ]


ex <- data.frame()
for (sim in unique(terp$simname)) for(c in unique(terp$cell)){
  tt <- terp[terp$simname==sim & terp$cell==c, ]
  thisg=data.frame()
  thisg[1, "sim"]<-sim
  thisg[1, "cell"]<-c
  thisg[1, "total_regs"] <- length(unique(tt$mutation))
  thisg[1, "terminal_all"] <- length(unique(tt$mutation[tt$type2 == "terminal all"]))
  thisg[1, "terminal_specific_this"] <- length(unique(tt$mutation[tt$type2 == "terminal specific this"]))
  thisg[1, "both"] <- length(intersect(unique(tt$mutation[tt$type2 == "terminal specific this"]),unique(tt$mutation[tt$type2 == "terminal all"])) )
  ex <- rbind(ex, thisg)
}

a<-as.matrix(ex[, c("terminal_all", "terminal_specific_this", "both")]) %>% apply(MAR=2, sum)
b<-a
b[1] <- b[1]-b[3]
b[2] <- b[2]- b[3]
c <- b/sum(b)

#on mce0fix
#          terminal_all terminal_specific_this                   both 
###       0.55753160             0.02485537             0.41761303 
#

