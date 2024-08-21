library(divo)
library(RColorBrewer)
library(ggraph)
library(igraph)
library(dplyr)
library(tibble)
library(data.table)
library(grid)
library(gridExtra)

compute_mh_mtx <- function(data, cdr3nt_col = 'CDR3_NT', condition = 'celltype', celltype.colSum.filter = 10){
  
  f <- as.formula( paste(cdr3nt_col, condition, sep = " ~ "))
  divo_mtx <- data %>% dplyr::filter(!is.na(!!sym(cdr3nt_col))) %>%
    plyr::count( vars = c(cdr3nt_col, condition)) %>% 
    dcast(formula = f, value.var="freq") %>% 
    replace(is.na(.), 0) %>%
    column_to_rownames(cdr3nt_col) %>% dplyr::select_if(colSums(.) >= celltype.colSum.filter)
  
  if(ncol(divo_mtx) > 1){
    divo_mh <- mh(divo_mtx,  graph = TRUE, csv_output = FALSE, PlugIn = FALSE, size = 1, saveBootstrap = FALSE)
  }else{
    divo_mh <- NULL 
  }
  return(divo_mh)
  
}
morisitae_overlap <- function(data, types){
  mh_list <- list(); mh_g <- list()
  for(type in types){
    print(type)
    cluster_mh <- compute_mh_mtx(data %>% dplyr::filter(Type == type),  cdr3nt_col = 'CDR3_NT', condition = 'celltype')
    if(!is.null(cluster_mh) ){
      mh_list[[type]] <- cluster_mh$Mean
    }
  }
  max_value <- max( unlist(mh_list)[!unlist(mh_list) %in% c(1, Inf, -Inf)] )
  
  return(mh_list)
  
}

# analyse morisita's index for pre and post-vaccination 
merged_tcr <- read.csv(file = 'tcr.table.csv')
merged_tcr <- merged_tcr %>% mutate(Type = Time) 
types <- unique(merged_tcr$Type)

cd3_mh_mtx <- morisitae_overlap(merged_tcr, types)

