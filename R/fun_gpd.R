##' Maximum-likelihood Function for a scale-blockwise GP model (without penalization parameter)
##'
##' @param ini  numeric: (B + 1) x 1 vector of data to be fitted
##' @param xdat  numeric: N x 1 vector of data to be fitted
##' @param yblock numeric: N x 1 covariate vector, with B pseudo-stationary time blocks labels,
##' for generalized linear modeling of the scale parameter.
##' @param tau scalar: a non-exceedance probability used to choose the threshold within each pseudo-stationary time block.
pccgpd_lik <- function(init, xdat, yblock, tau) {
  # number of shape parameters:
  npsc <- length(unique(yblock))
  # number of points (global):
  n <- length(xdat)
  # data (xdat) splitted by blocks (a list):
  x_blkwise <- split(xdat, yblock)
  # number of points within each block (a list):
  n_blkwise <- Map(f = function(x) length(x),
                   x = x_blkwise, USE.NAMES = FALSE)
  # threshold for each block (a list):
  thr_blkwise <- Map(f = function(x, y) quantile(x, probs = y, names = FALSE),
                     x = x_blkwise, y = tau, USE.NAMES = FALSE)
  # threshold point-wise for the entire time series (a vector of size n):
  u <- unlist(Map(f = function(x, y) rep(x, y),
                  x = thr_blkwise, y = n_blkwise, USE.NAMES = FALSE))
  # exceedances vector:
  xdatu <- xdat[xdat > u]
  # exceedances indexes:
  xind <- (1:n)[xdat > u]
  # threshold value for each exceedance:
  u <- u[xind]
  # stationarity time blocks in factor format (a data frame):
  ydat <- data.frame(Block = factor(yblock, ordered = T))
  # dummy matrix of the stationary time blocks for the non-statinary scale parameter:
  ymat <- model.matrix(~ Block - 1, data = ydat, contrasts.arg = list(Block = diag(nlevels(ydat$Block))))
  # exceedances dummy matrix of the stationary time blocks:
  sigmat <- ymat[xind, ]
  # 1-matrix for the stationary shape parameter:
  shmat <- as.matrix(rep(1, length(xdatu)))

  # PPC GPD log-likelihood:
  gpd.lik <- function(a) {
    sc <- sigmat %*% a[1:npsc]
    xi <- shmat %*% a[npsc + 1]
    y <- (xdatu - u)/sc
    y <- 1 + xi * y
    if (min(sc) <= 0) {
      l <- 10^6
    } else {
      if (min(y) <= 0) {
        l <- 10^6
      } else {
        l <- sum(log(sc)) + sum(log(y) * (1/xi + 1))
      }
    }
    return(l)
  }

  return(gpd.lik(init))

}


##' Maximum-likelihood optimization for the PPC GP
##'
##' @param xdat  numeric: a vector of data to be fitted
##' @param yblock numeric: a covariate vector, with pseudo-stationary time blocks identification,
##' for generalized linear modeling of the scale parameter.
##' The number of rows should be the same as the length of xdat.
##' @param tau scalar: a non-exceedance probability used to choose the threshold within each pseudo-stationary time block.
##' @param lambda scaler: a penalization parameter
##' @method string: the method to be used
##' @maxit scalar: the maximum number of iterations
##' @... further arguments to be passed to optim
ppcgpd_fit <- function(xdat, yblock, tau, lambda = 0,
                       method = "BFGS", maxit = 5000, ...) {

  # number of shape parameters:
  npsc <- length(unique(yblock))
  # number of points (global):
  n <- length(xdat)
  # data (xdat) splitted by blocks (a list):
  x_blkwise <- split(xdat, yblock)
  # number of points within each block (a list):
  n_blkwise <- Map(f = function(x) length(x),
                   x = x_blkwise, USE.NAMES = FALSE)
  # threshold for each block (a list):
  thr_blkwise <- Map(f = function(x, y) quantile(x, probs = y, names = FALSE),
                     x = x_blkwise, y = tau, USE.NAMES = FALSE)
  # threshold point-wise for the entire time series (a vector of size n):
  u = threshold <- unlist(Map(f = function(x, y) rep(x, y),
                              x = thr_blkwise, y = n_blkwise, USE.NAMES = FALSE))
  # exceedances vector:
  xdatu <- xdat[xdat > u]
  # exceedances indexes:
  xind <- (1:n)[xdat > u]
  # threshold value for each exceedance:
  u <- u[xind]
  # exceedances splitted by blocks (a list):
  xdatu_blkwise <- split(xdatu, u)
  # stationarity time blocks in factor format (a data frame):
  ydat <- data.frame(Block = factor(yblock, ordered = T))
  # dummy matrix of the stationary time blocks for the non-statinary scale parameter:
  ymat <- model.matrix(~ Block - 1, data = ydat, contrasts.arg = list(Block = diag(nlevels(ydat$Block))))
  # exceedances dummy matrix of the stationary time blocks:
  sigmat <- ymat[xind, ]
  # 1-matrix for the stationary shape parameter:
  shmat <- as.matrix(rep(1, length(xdatu)))
  # initial values for the scale parameter:
  sig0 <- sqrt(6 * var(xdatu))/pi
  siginit <- rep(sig0, npsc)
  # initial values for the shape parameter:
  shinit <- 0.1
  # initial values for the optim function:
  init <- c(siginit, shinit)
  # PPC GPD log-likelihood:
  gpd.lik <- function(a) {
    sc <- sigmat %*% a[1:npsc]
    xi <- shmat %*% a[npsc + 1]
    y <- (xdatu - u)/sc
    y <- 1 + xi * y
    if (min(sc) <= 0)
      l <- 10^6
    else {
      if (min(y) <= 0)
        l <- 10^6
      else {
        l <- sum(log(sc)) + sum(log(y) * (1/xi + 1)) + lambda * (sum((sc - mean(sc)) ^ 2) / npsc)
      }
    }
    l
  }
  # optimizer:
  x <- optim(init, gpd.lik, hessian = TRUE, method = method,
             control = list(maxit = maxit))
  # ML estimates exceedances-wise:
  sc <- sigmat %*% x$par[1:npsc]
  xi <- shmat %*% x$par[npsc + 1]
  # outputs list:
  z <- list()
  # total number of exceedances
  z$nexc <- length(xdatu)
  # number of exceedances block-wise:
  z$nexc_blkwise <- unlist(lapply(xdatu_blkwise, length), use.names = FALSE)
  # optim-convergence code - 0 indicates successful convergence:
  z$conv <- x$convergence
  # single numeric giving the negative log-likelihood value:
  z$nllh <- x$value
  # numeric vector giving the MLE's for the scale and shape parameters, resp:
  z$mle <- x$par
  # the covariance matrix:
  z$cov <- solve(x$hessian)
  # vector containing the standard errors:
  z$se <- sqrt(diag(z$cov))
  # number of data points (global):
  z$n <- n
  # number of points within each block:
  z$n_blkwise <- unlist(n_blkwise)
  # the proportion of data points that lie above the threshold (global):
  z$rate <- length(xdatu) / n
  # block-wise proportion of data points that lie above the threshold:
  z$rate_blkwise <- z$nexc_blkwise / z$n_blkwise
  # A matrix with three columns containing the ML estimates of the scale and shape parameters,
  # and the threshold, at each data point.
  z$vals <- cbind(sc, xi, u)
  colnames(z$vals) <- c("scale", "shape", "threshold")
  # the data that lie above the threshold. For non-stationary models, the data is standardized.
  z$data <- -log(as.vector((1 + (xi * (xdatu - u)) / sc) ^ (-1 / xi)))
  # vector of thresholds (global):
  z$threshold <- threshold
  # the data that has been fitted.
  z$xdata <- xdat

  return(z)

}


##' cv procedure chunkwise/foldwise
##'
##' @param chk scalar: a given chunk/fold to train and test
##' @param channel_df data frame: a N x 4 a data frame with column names: "Block", "Chunk", "Time", "Channel"
##' @param tau scalar: a NEP
##' @param lambda scalar: the value of the penalization parameter to be used
chuks_fit <- function(chk, channel_df, tau, lambda) {

  chks_trng <- channel_df %>% dplyr::filter(Chunk %in% 1:chk) # trainning set
  chks_test <- channel_df %>% dplyr::filter(Chunk %notin% 1:chk) # test set

  # Optimization
  trng_fit <- ppcgpd_fit(chks_trng$Channel, chks_trng$Block, tau, lambda)
  trng_par <- trng_fit$mle

  # Evaluating (negative) loglikelihood on the test set
  test_fit <- pccgpd_lik(trng_par, chks_test$Channel, chks_test$Block, tau)

  return(test_fit)

}


##' Auxiliar function to the 10-fold walk-forward cross-validation procedure
##'
##' @param lambda scalar: scalar: the value of the penalization parameter to be used
##' @param chks_n scalar: number of chunks in the channel_df data frame
##' @param channel_df data frame: a N x 4 a data frame with column names: "Block", "Chunk", "Time", "Channel"
##' @param tau scalar: a NEP
cv_fit <- function(lambda, chks_n = 10, channel_df, tau) {


  # Removing 1 chunk in turn and computing the loglikelihood
  chks_loglik <- vector("numeric", chks_n - 5)
  for (chk in seq_along(chks_loglik)) {
    chks_loglik[chk] <- chuks_fit(chk + 4, channel_df, tau, lambda)
  }
  cv_output <- c(lambda, chks_loglik)
  names(cv_output) <- c("lambda", paste0("chk_", 1:(chks_n - 5)))

  return(cv_output)

}

##' Bootstrap helper function
##'
##' @param boot_i scalar: a given
##' @param channel_df data frame: a N x 4 a data frame with column names: "Block", "Chunk", "Time", "Channel"
##' @param boot_index data frame: a (boot_n + 1) column data frame with the block-bootstrap samples of the time index
##' @param boot_nep numeric: a (boot_n + 1) numeric vector with a uniform sample of the NEP
##' @param best_lambda scalar: the value of the penalization parameter obtained at the CV procedure
boot_fit <- function(boot_i, channel_df, boot_index, boot_nep, best_lambda) {

  boot_df <- channel_df[boot_index[, boot_i], ]

  fit <- ppcgpd_fit(boot_df$Channel, boot_df$Block, boot_nep[boot_i], best_lambda)

  boot_par <- fit$mle
  names(boot_par) <- c(paste0("nu_", 1:(length(boot_par) - 1)), "xi")
  boot_msg <- fit$conv; names(boot_msg) <- "convergence"
  boot_id <- boot_i; names(boot_id) <- "boot_id"
  return(c(boot_par, boot_msg, boot_id))

}

##' Fit the marginal PPC GP model to a given marginal channel.
##' It applies the CV procedure and run the block bootstrap.
##'
##' @param df data frame: a data frame with column names "Chunk", "Block", "Time", and also the values in c(conditioning, associates)
##' @param phase_id scalar: Phase to which you want to fit the model: 0 -> Pre-seizure / 1 -> Post-seizure
##' @param lead_id scalar: lag to which you want to fit the model (0, 10, 25, 50)
##' @param phase string: column name identifying pre- and post-seizure phases
##' @param covariate string: column name identifying non-stationary covariate for the model
##' @param conditioning string: column name identifying the conditioning variate for the H&T model
##' @param leading string: column name identifying the time lag
##' @param associates string: column name of the associate variables for the model
##' @param blks_n scalar: The number of pseudo-stationary blocks to split the covariate
##' @param chks_n scalar: The number of chunks (folds) to split each pseud-stationary block for CV purpose
##' @param lbda_n scalar: Quantity of penalization parameters (lambdas) to test at the CV procedure
##' @param boot_size scalar: Size of the block-bootstrap within each pseudo-stationary block
##' @param boot_n scalar: The number of bootstrap samples
##' @param out_dir string: Path where to save outputs
##' @param marginal_id scalar: number identifying the marginal channel to be fitted
##' @param lbda_rng numeric: a 1x2 vector where to sample lambda for the marginal channel being fitted
##' @param nep numeric: a 1x2 vector with the range to sample the thresholds of the marginal GP fit
marginal_gpd_fit <- function(df,
                             phase_id,
                             lead_id,
                             phase,
                             covariate,
                             conditioning,
                             leading,
                             associates,
                             blks_n,
                             chks_n,
                             lbda_n,
                             boot_size,
                             boot_n,
                             out_dir,
                             marginal_id,
                             lbda_rng,
                             nep) {


  # Variables names, paths to save outputs, other setups
  phase_name <- ifelse(phase_id == 0, "Pre-seizure", "Post-seizure")
  marginals <- c(conditioning, associates)
  marginals_n <- length(marginals)
  file_root <- paste0(phase_id, "_", lead_id, "_gpd_", marginals[marginal_id], "_")
  file_path <- paste0(out_dir, file_root)
  lambda <- seq(10 ^ lbda_rng[1], 10 ^ lbda_rng[2], length.out = lbda_n)
  tau <- median(nep) # The NEP to be used to the non-bootstrap fit


  # Channelwise approach - Selecting one specific channel to fit the PPC-GPD model
  channel_df <- df[, c("Chunk", "Block", covariate, marginals[marginal_id])]
  names(channel_df)[4] <- "Channel"


  # Cross-validation
  cv <- foreach(l = 1:lbda_n,
                .combine = 'rbind',
                .packages = 'tidyverse',
                .export = c("pccgpd_lik", "ppcgpd_fit", "%notin%", "chuks_fit", "cv_fit")) %dopar% {
                  cv_fit(lambda[l], chks_n, channel_df, tau)
                }


  # CV result
  cv_output <- data.frame(lambda = cv[, 1], loglik = rowSums(cv[, 3:6])/4)
  cv_df <- dplyr::filter(cv_output, loglik < 10 ^ 6)
  if (nrow(cv_df) == 0) {
    warning("CV failed! Penalty parameter forced to be zero.")
    best_lambda <- 0
  } else {
    best_lambda <- cv_df[which.min(cv_df[, 2]), 1]
  }


  # Bootstrap
  if (marginal_id == 1) {
    boot_nep <- c(tau, runif(boot_n, min = nep[1], max = nep[2]))
    boot_index <- bootstrap_df(channel_df, boot_size, boot_n)
  } else {
    ref_out <- read_rds(paste0(out_dir, phase_id, "_", lead_id, "_gpd_", conditioning, "_output.rds"))
    boot_nep <- ref_out$boot_nep
    boot_index <- ref_out$boot_index
  }
  boot <- foreach(
    b = 1:(boot_n + 1),
    .combine = 'rbind',
    .packages = 'tidyverse',
    .export = c("ppcgpd_fit", "boot_fit")) %dopar% {
      boot_fit(b, channel_df, boot_index, boot_nep, best_lambda)
    }


  # ACF and PACF
  plot_p_acf <- acf_plot(channel_df, marginals[marginal_id], phase_name, lag.max = boot_size + 50, ci = .95, boot_size = boot_size,
                         lead_id = lead_id)

  png(file = paste0(file_path, "acf", ".png"),  width = 1200, height = 600)
  plot_acf <- plot_p_acf$acf
  print(plot_acf)
  dev.off()

  png(file = paste0(file_path, "pacf", ".png"),  width = 1200, height = 600)
  plot_pacf <- plot_p_acf$pacf
  print(plot_pacf)
  dev.off()


  # CV Plot
  png(file = paste0(file_path, "cross_validation", ".png"),  width = 600, height = 600)
  plot_cv <- cv_plot(cv_df, best_lambda, marginals[marginal_id], phase_name, lead_id)
  print(plot_cv)
  dev.off()


  # Dispersion Plots - Original and Transformed data
  uniform <- u_scale(par = boot[1, 1:(blks_n + 1)], channel_df, tau)
  laplace <- laplace_scale(uniform)
  fit_df <- data.frame(
    channel_df,
    Uniform = uniform,
    Laplace = laplace
  )

  edges_df <- fit_df %>%
    dplyr::group_by(Block) %>%
    summarise(min = min(Time),
              max = max(Time)) %>%
    ungroup()
  edges_bars <- round(c(edges_df$min, last(edges_df$max)), 1)

  png(file = paste0(file_path, "dispersion", ".png"),  width = 1200, height = 1200/3)
  plot_disp <- disp_plot(fit_df, tau, boot_nep, edges_bars, marginals[marginal_id], phase_name, lead_id)
  print(plot_disp)
  dev.off()


  # GPD Parameters
  nu_lwr <- apply(boot[, 1:blks_n], 2, quantile, probs = .05, names = FALSE)
  nu_upp <- apply(boot[, 1:blks_n], 2, quantile, probs = .95, names = FALSE)
  nu_med <- apply(boot[, 1:blks_n], 2, mean)
  nu_df <- data.frame(
    Min = edges_df$min,
    Max = edges_df$max,
    Par = c(boot[1, 1:blks_n]),
    Med = c(nu_med),
    Lwr = c(nu_lwr),
    Upp = c(nu_upp)
  )

  xi_hat0 <- boot[1, (blks_n + 1)]
  xi_df <- data.frame(Par = boot[-1, (blks_n + 1)])

  png(file = paste0(file_path, "parameters", ".png"),  width = 1200, height = 1200/2)
  plot_par <- par_plot(nu_df, xi_df, xi_hat0, edges_bars, marginals[marginal_id], phase_name, tau, boot_n, lead_id)
  print(plot_par)
  dev.off()


  # Shape Stability Plot
  png(file = paste0(file_path, "stability", ".png"),  width = 600, height = 600)
  plot_stb <- stb_plot(tau, nep, boot_nep, boot, marginals[marginal_id], phase_name, lead_id)
  print(plot_stb)
  dev.off()


  # Q-Q plots

  # Overall
  qqplot_overall <- qqdf(fit_df$Laplace)

  png(file = paste0(file_path, "qqplot_overall", ".png"),  width = 600, height = 600)
  plot_qqall <- allqq_plot(qqplot_overall, marginals[marginal_id], phase_name, lead_id)
  print(plot_qqall )
  dev.off()

  # Blockwise
  qqplot_blkwise <- list()
  for (my_block in 1:blks_n) {
    x_blkwise <- dplyr::filter(fit_df, Block == my_block) %>% pull(Laplace)
    qqplot_blkwise[[my_block]] <- qqdf(x_blkwise)
    qqplot_blkwise[[my_block]]$Block <- my_block
  }
  qqplot_blkwise <- do.call("rbind", qqplot_blkwise)

  n_exc <- unlist(lapply(info_blkwise(channel_df, tau)$Exceedances, length))
  strip_labs <- paste0("Block ", 1:blks_n, " | #Exceedances: ", n_exc)
  names(strip_labs) = 1:blks_n

  png(file = paste0(file_path, "qqplot_blkwise", ".png"),  width = 1200, height = 600)
  plot_qqblk <- blkqq_plot(qqplot_blkwise, strip_labs, marginals[marginal_id], phase_name, lead_id)
  print(plot_qqblk)
  dev.off()


  # Saving outputs
  output <- list(
    inputs = data.frame(
      phase = phase_id,
      lead = lead_id,
      marginals_n = marginals_n,
      marginal_id = marginal_id,
      marginal = marginals[marginal_id],
      blks_n = blks_n,
      chks_n = chks_n,
      lbda_n = lbda_n,
      lbda_lwr = lbda_rng[1],
      lbda_upp = lbda_rng[2],
      nep = tau,
      nep_lwr = nep[1],
      nep_upp = nep[2],
      boot_n = boot_n,
      boot_blk_size = boot_size
    ),
    channel_data = channel_df,
    n_exceedances = n_exc,
    cv_chkwise = cv,
    cv_summary = cv_df,
    cv_best_lambda = best_lambda,
    boot_nep = boot_nep,
    boot_index = boot_index,
    boot_summary = boot,
    nu_hat = nu_df,
    xi_hat = boot[, (blks_n + 1)],
    transformed_data = fit_df
  )

  write_rds(output, paste0(file_path, "output.rds"))

  invisible(output)

}
