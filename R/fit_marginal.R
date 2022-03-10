rm(list = ls())

source("R/setup.R")

#### Loading the dataset ----
df_eeg <- readr::read_csv("data/mwe_dataset.csv")

#### General Inputs ----
# These outputs are the same for all variables in the analysis
phase_id <- 0 # Phase to which you want to fit the model: 0 -> Pre-seizure / 1 -> Post-seizure
lead_id <- 0 # Lag to which you want to fit the model (0, 10, 25, 50)
phase <- "Phase" # Column name identifying pre- and post-seizure phases
covariate <- "Time" # Column name identifying non-stationary covariate for the model
conditioning <- "T3" # Column name identifying the conditioning variate for the H&T model
leading <- "Lead" # Column name identifying the time lag
associates <- c("F7", "F8") # Column name of the associate variables for the model
blks_n <- 12 # The number of pseudo-stationary blocks to split the covariate
chks_n <- 10 # The number of chunks (folds) to split each pseud-stationary block for CV purpose
lbda_n <- 50 # Quantity of penalization parameters (lambdas) to test at the CV procedure
boot_size <- 25 # Size of the block-bootstrap within each pseudo-stationary block
boot_n <- 500 # The number of bootstrap samples ---> The same for all channels
out_dir <- paste0("out/marginals/lead_", lead_id, "/") # Path where to save outputs
nep = c(.90, .95) # Range to sample the thresholds of the marginal GP fit

#### Individual inputs ----
# Intervals where to search the optimal value of lambda for each variate
lambdas <- rbind(
  c(-11, -1),  # T3  # 1
  c(-10, -.5), # F7  # 2
  c(-12, -1)   # F8  # 3
)


# Adding to the data frame a) time blocks for the PPC model and b) folds for the CV procedure
source("R/fun_block_split.R")
df <- block_dataframe(
  df_eeg, phase, phase_id, leading, lead_id,
  covariate, conditioning, associates, blks_n, chks_n
)

#### Marginal GP fit ----
source("R/fun_bootstrap.R")
source("R/fun_gpd.R")
source("R/fun_laplace_scale.R")
source("R/plots_laplace_qqplot.R")
source("R/plots_marginal.R")
ini <- proc.time()
channels_to_fit <- c(conditioning, associates)
for (i in seq_along(channels_to_fit)) {
  marginal_gpd_fit(
    df, phase_id, lead_id, phase, covariate, conditioning,
    leading, associates, blks_n, chks_n, lbda_n, boot_size, boot_n, out_dir,
    marginal_id = i, lbda_rng = lambdas[i, ], nep
  )
}
fini = proc.time() - ini
round(fini[3] / 60, 4)
