log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# libs to use
library(tidyverse)
library(xCell2)
library(parallel)
library(ComplexHeatmap)
library(seriation)

# ---- inputs from Snakemake ----
gene_counts_path <- snakemake@input[["gene_counts"]]
transcript_index <- snakemake@input[["transcript_ref"]]
ref_path         <- snakemake@input[["ref"]]
ref_set_name     <- gsub("\\..*", "", basename(ref_path))

# ---- getting data from sleuth model and collapse transcripts to gene symbol and removing non expressed genes ----
data_matrix <- read_tsv(gene_counts_path) %>%
  group_by(gene) |> 
  summarize(
    across(where(is.numeric), sum)) |> 
  ungroup() |> 
  filter(!is.na(gene))

# checking for duplicates ! keep an eye on output in generell there should be no duplicated gene
dupes <- data_matrix[duplicated(data_matrix$gene), ]
if (nrow(dupes) > 0) {
  message("Duplicated gene entries found:")
  message(paste(dupes$gene, collapse = ", "))
}

data_input <- data_matrix %>%
  column_to_rownames(var = "gene") %>%
  mutate_all( ~ log10( . + 3))


# ---- load the .rda reference set (replaces get_refs/readRDS) ----
loaded_name <- load(ref_path)          # load() returns the object name(s) created
ref_obj     <- get(loaded_name[1])     # retrieve generically, whatever it's called

# ---- core decon logic (was: get_decode) ----
gene_ls_ref  <- ref_obj@genes_used
genes_shared <- intersect(gene_ls_ref, data_matrix$gene)
input_shared <- length(genes_shared) / length(gene_ls_ref)
input_shared <- floor(input_shared * 20) / 20

message(paste(ref_set_name, input_shared, sep = ": "))

xcell2_res <- xCell2::xCell2Analysis(
  mix = data_input,
  xcell2object = ref_obj,
  minSharedGenes = input_shared
)

### export of convoluted matrices as table
xcell_xprt <- as_tibble(xcell2_res, rownames = "cell_type")

write_tsv(xcell_xprt, snakemake@output[["decon_scores"]])

### plotting convoluted matrices 
decon_heatmap <- xcell2_res %>% 
t() %>%
scale

decon_clean <- decon_heatmap[, !colSums(is.na(decon_heatmap))]
### calculating sample and gene deistribution by seriate   
col_dis    <- dist(decon_clean, method = "euclidean")
col_clst   <- hclust(col_dis, method = "ward.D2")
col_seriat <- reorder(col_clst, col_dis, method = "OLO")
  
row_dis    <- dist(t(decon_clean), method = "euclidean")
row_clst   <- hclust(row_dis, method = "ward.D2")
row_seriat <- reorder(row_clst, row_dis, method = "OLO")

### plotting and layout of heatmap 
ht <- Heatmap(
      t(decon_clean),
      cluster_columns = as.dendrogram(col_seriat),
      cluster_rows    = as.dendrogram(row_seriat),
      row_names_max_width = unit(200, "mm")
)

### setting dimensions
n_cols <- nrow(decon_clean)
n_rows <- ncol(decon_clean)

### scaling factor 
row_height <- 8
col_width  <- 5

plot_height <- max(10, n_rows * row_height)
plot_width  <- max(6,  n_cols * col_width)


tiff(snakemake@output[["heatmap"]], res = 150, pointsize = 10, units = "mm",
 compression = "lzw", width = plot_width, height = plot_height)
print(
  draw(ht,
  heatmap_legend_side = "left",
  align_heatmap_legend = "heatmap_top")
)
dev.off()

sink(type = "message")
sink(type = "output")