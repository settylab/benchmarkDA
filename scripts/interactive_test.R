args <- list(
  data_RDS = "/fh/fast/setty_m/user/dotto/benchmarkDA/data/real/covid19-pbmc/single-cell-atlas-pbmc-sars-cov2_sce.rds",
  seed = 43,
  population = "RBC",
  pop_enrichment = 0.75,
  pop_col = "cell.type.coarse",
  reduced.dim = "PCA", # Default value, change if needed
  k = 50, # Default value, change if needed
  data_id = "covid19-pbmc",
  make_batch_effect = "no",
  outdir = "/fh/fast/setty_m/user/dotto/benchmarkDA/data/real/covid19-pbmc/"
)

source("make_bm_data.R", echo=TRUE)
