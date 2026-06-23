log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


library(tidyverse)
#library(biomaRt)
library(xCell2)
library(parallel)
library(ComplexHeatmap)
library(seriation)

#### deconvolution starts here
### checking for duplicates
data_matrix_2 <- read_tsv(snakemake@input[["gene_counts"]], header = TRUE)

data_input <- data_matrix_2 %>% 
  column_to_rownames(var = "hgnc_symbol") %>% 
  dplyr::select(all_of(contains("rna_")))

get_refs <- function(a){
  readRDS(paste0(path_data, "ref_data_sets/", ref_sets[a])) 
}
ref_sets <- list.files(paste0(path_data,"ref_data_sets/"))
ref_set_lst <- lapply(1:length(ref_sets), get_refs)
names(ref_set_lst) <- ref_sets

get_decode <- function(a){
  tryCatch({
    xcell2_pan_res <-
      xCell2::xCell2Analysis(
        mix = data_input,
        xcell2object = ref_set_lst[[a]],
        minSharedGenes = shared_genes[a]
        )
    },
    error = function(e) {
      # instant if function(a) throws an error
      message("NOT WORKING!")
      message(conditionMessage(e))
      return(NA) # returns NAs to generate an output 
    })
}
shared_genes <- c(0.8, 0.9, 0.9, 0.9, 0.65, 0.9)
tst <- lapply(1:length(ref_set_lst), get_decode)
names(tst) <- ref_sets


### plotting convoluted matrices 
all_heatmaps <- function(a){
  tst_plt <- t(tst[[a]]) %>% 
    scale 
  
  tst_plt <-  tst_plt[, !colSums(is.na(tst_plt))]
  
  col_dis    <- dist(tst_plt, method = "euclidean")
  col_clst   <- hclust(col_dis, method = "ward.D2")
  col_seriat <- reorder(col_clst, col_dis, method = "OLO")
  
  row_dis    <- dist(t(tst_plt), method = "euclidean")
  row_clst   <- hclust(row_dis, method = "ward.D2")
  row_seriat <- reorder(row_clst, row_dis, method = "OLO")
  
  Heatmap(t(tst_plt),
          cluster_columns = as.dendrogram(col_seriat),
          cluster_rows    = as.dendrogram(row_seriat)
  )
}

all_plts <- lapply(1:6, all_heatmaps)

for (i in 1:6) {
  tiff(paste0("decon_heatmap_", paste0(names(ref_set_lst)[i]), ".tif"), res = 150, 
       pointsize = 10, units = "mm", compression = "lzw", width = 150, height = 200)
  print(
    all_plts[[i]]
  )
  dev.off()
  
}

