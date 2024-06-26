rm(list=ls())
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(DAseq)
  library(miloR)
  library(tibble)
  library(tidyr)
  library(dplyr)
  library(igraph)
  library(cydar)
  library(pdist)
  library(reshape2)
})


run_milo <- function(sce, condition_col, sample_col, reduced.dim="PCA",
                     k=15, d=30, prop=0.1, returnMilo = TRUE,
                     batch_col=NULL) {
  ## Make design matrix
  design_df <- as.tibble(colData(sce)[c(sample_col, condition_col, batch_col)]) %>%
    distinct() %>%
    column_to_rownames(sample_col)
  if (is.null(batch_col)) {
    design <- formula(paste('~', condition_col, collapse = ' '))  
  } else {
    design <- formula(paste('~', batch_col, "+", condition_col, collapse = ' '))
  }
  
  ## Build graph neighborhoods
  milo <- Milo(sce)
  milo <- buildGraph(milo, k=k, d=d, reduced.dim = reduced.dim)
  milo <- makeNhoods(milo, prop = prop, k=k, d=d, reduced_dims = reduced.dim)
  ## Test DA
  milo <- miloR::countCells(milo, meta.data = data.frame(colData(milo)), sample=sample_col)
  milo <- calcNhoodDistance(milo, d=d, reduced.dim = reduced.dim)
  DA_results <- testNhoods(milo, design = design, design.df = design_df, reduced.dim=reduced.dim)
  if (isTRUE(returnMilo)) {
    return(list(Milo=milo, DAres=DA_results))
  } else {
    DA_results 
  }
}

milo2output <- function(milo, da_res, out_type="continuous", alpha=0.1){
  if (out_type=="continuous") {
    da.cell.mat <- milo@nhoods %*% da_res$logFC
    da.cell <- da.cell.mat[,1]
  } else {
    da.nhoods <- ifelse(da_res$SpatialFDR < alpha, ifelse(da_res$logFC > 0, "PosLFC", 'NegLFC'), "NotDA")
    da.nhoods.mat <- sapply(unique(da.nhoods), function(x) as.numeric(da.nhoods==x))
    da.cell.mat <- milo@nhoods %*% da.nhoods.mat
    da.cell <- apply(da.cell.mat, 1, function(x) colnames(da.cell.mat)[which.max(x)])
  }
  da.cell
}

# main

# read from csv
script_dir <- dirname(sys.frame(1)$ofile)
sce <- readRDS(file.path(script_dir, "../../data/synthetic/linear/linear_data_bm.RDS"))
X_pca <- read.csv(file.path(script_dir, "../../data/synthetic/linear/benchmark_linear_pop_M1_enr0.75_seed44_batchEffect0.75.pca.csv")) %>% column_to_rownames()
coldata <- read.csv(file.path(script_dir, "../../data/synthetic/linear/benchmark_linear_pop_M1_enr0.75_seed44.coldata.csv")) %>% column_to_rownames()
                     
source(paste0(script_dir, "/../../scripts/benchmark_utils.R")) # import calculate_outcome

# X_pca <- read.csv("/fh/fast/setty_m/user/dotto/benchmarkDA/data/real/bcr-csv/X.csv", header = FALSE)
# coldata <- read.csv("/fh/fast/setty_m/user/dotto/benchmarkDA/data/real/bcr-csv/obs.csv") %>% column_to_rownames()
# rownames(X_pca) <- rownames(coldata)
# sce <- SingleCellExperiment(t(X_pca))

## Add reduced dim + coldata to sce
colData(sce) <- DataFrame(coldata)
reducedDim(sce, "pca_batch") <- as.matrix(X_pca)

# run the milo method
condition_col <- "synth_labels"
sample_col <- "synth_samples"

set.seed(42)
milo_res <- run_milo(sce, condition_col = condition_col, sample_col = sample_col,
                     reduced.dim = 'pca_batch', d=50, k=30)
out <- milo2output(milo_res$Milo, milo_res$DAres, out_type = "label", alpha=0.933)

# post-process outputs
bm <- data.frame(bm_out=out)
bm$true_prob <- sce$Condition2_prob
bm$true <- sce$true_labels
long_bm <- pivot_longer(bm, cols = bm_out, names_to='milo', values_to="pred")
long_bm[["method"]] <- 'milo'

# compute metrics
metrics <- calculate_outcome(long_bm)
print(metrics)

