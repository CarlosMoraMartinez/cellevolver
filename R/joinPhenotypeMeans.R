
library(tidyverse)

a <- map(f, function(x){tab <- read_tsv(x);tab$dataset <- x;return(tab)})%>% 
  bind_rows %>% 
  mutate(mutation_type = ifelse(grepl("sitesPhenotype", dataset), "TFBS", "TF"), 
         condition=str_extract(dataset, "4cell_mce[0-9a-zA-Z]+(_mutb)?")) %>% 
  select(-dataset)

 b <- a %>% gather("measure", "value", 2:4) %>%
  unite("mescol", mutation_type, measure, sep="_") %>% 
  spread(mescol, value)

write_tsv(b, "200921_PhenotypeMeans.tsv")
