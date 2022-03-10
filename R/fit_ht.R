rm(list = ls())

source("R/setup.R")

#### General Inputs ----
# These outputs are the same for all variables in the analysis
lead_id <- 0 # Lag to which you want to fit the model (0, 10, 25, 50)
phase <- "Phase" # Column name identifying pre- and post-seizure phases
covariate <- "Time" # Column name identifying non-stationary covariate for the model
conditioning <- "T3" # Column name identifying the conditioning variate for the H&T model
leading <- "Lead" # Column name identifying the time lag
associates <- c("F7", "F8") # Column name of the associate variables for the model
lbda_n <- 50 # Quantity of penalization parameters (lambdas) to test at the CV procedure
nep <- c(.88, .92) # Range to sample the thresholds of the H&T fit
boot_n <- 500 # The number of bootstrap samples ---> The same number used for the marginals fit
marginal_fit_dir <- paste0("out/marginals/Lead_", lead_id, "/") # Path where to load the marginals fit outputs
ht_fit_dir <- paste0("out/ht/Lead_", lead_id, "/") # Path where to save outputs

# Ensure the same NEP bootstrap sample for both pre and post-seizure moments
tau <- median(nep)
boot_nep <- c(tau, runif(boot_n, min = nep[1], max = nep[2]))

# Fitting the models for pre and post-seizure moments
source("R/fun_ht.R")
source("R/plots_ht.R")
ini <- proc.time()
ht_lead0_pre <- ht_fit(phase_id = 0, lead_id, phase, covariate, conditioning, leading, associates, lbda_n,
                       lbda_rng = c(2.5, 3.5),
                       nep, tau, boot_nep,  boot_n, marginal_fit_dir, ht_fit_dir)

ht_lead0_post <- ht_fit(phase_id = 1, lead_id, phase, covariate, conditioning, leading, associates, lbda_n,
                        lbda_rng = c(2, 3.25),
                        nep, tau, boot_nep, boot_n, marginal_fit_dir, ht_fit_dir)

fini = proc.time() - ini
round(fini[3] / 60, 4)
