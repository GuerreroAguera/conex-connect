laplace_negloglik <- function(y1, yd, lambda, almat, a, b, m, s) {

  BigNumber <- 10^40

  a <- almat %*% a

  mu <- a * y1 + m * y1 ^ b
  sig <- s * y1 ^ b

  res <- sum(0.5 * log(2*pi) + log(sig) + 0.5 * ((yd - mu)/sig)^2) +
    lambda * (sum((a - mean(a)) ^ 2) / length(a))

  if (is.infinite(res)) {
    if (res < 0) {
      res <- -BigNumber
    } else {
      res <- BigNumber
    }
  }

  res
}

laplace_profile_negloglik <- function(y1, yd, lambda, almat, a, b) {

  Z <- (yd - y1 * (almat %*% a)) / (y1 ^ b)

  m <- mean(Z)
  s <- sd(Z)

  res <- laplace_negloglik(y1, yd, lambda, almat, a, b, m = m, s = s)
  res <- list(profLik = res, m = m, s = s)
  res
}

Qpos <- function(ini, y1, yd, lambda, almat) {

  n <- length(ini) - 1

  a <- ini[1:n]
  b <- ini[n + 1]

  res <- laplace_profile_negloglik(y1, yd, lambda, almat, a = a, b = b)
  res$profLik
}


ppcht_test_fit <- function(par_hat, df_wide, conditioning, associate_col, tau) {

  # number of alpha parameters:
  npal <- length(unique(df_wide[, "Block"]))
  # Converting Block to factor.
  df_wide[, "Block"] <- factor(df_wide[, "Block"], ordered = T)
  # conditioning variate vector
  ycond <- df_wide[, conditioning]
  # threshold of the conditioning variate
  thr_cond <- LaplacesDemon::qlaplace(tau)
  # exceedances indexes:
  u_ind <- which(ycond > thr_cond)
  # exceedances of the conditioning variate + associates (data frame)
  df_cond <- df_wide[u_ind, ]
  # stationarity time blocks in factor format (a data frame):
  yblock <- data.frame(Block = factor(
    df_cond[, "Block"],
    levels = levels(df_wide[, "Block"]),
    ordered = T)
  )
  # dummy matrix of the stationary time blocks for the non-statinary scale parameter:
  almat <- model.matrix(~ Block - 1, data = yblock,
                        contrasts.arg = list(Block = diag(nlevels(yblock$Block))))

  # Fit H&T to the associate variate "associate_col":

  fit <- numeric(length(associate_col))
  i_aux <- min(associate_col) - 1
  for (i in associate_col) {
    j <- i - i_aux
    # associate variate exceedances df:
    df_assoc <-  df_cond[ , c(1:4, i)]
    # conditining and associate variates data (vector):
    y1 <- df_assoc[, 4]
    yd <- df_assoc[, 5]
    # computing mle:
    fit[[j]] <- Qpos(par_hat[1:(npal + 1), j], y1, yd, lambda = 0, almat)
  }

  return(mean(fit))

}


ppcht_fit <- function(df_wide, conditioning, associate_col, tau, lambda) {

  # block edges:
  edges_df <- df_wide %>%
    dplyr::group_by(Block) %>%
    summarise(min = min(Time),
              max = max(Time)) %>%
    ungroup()

  # number of alpha parameters:
  npal <- length(unique(df_wide[, "Block"]))
  # number of points (global):
  n <- nrow(df_wide)
  # Converting Block to factor.
  df_wide[, "Block"] <- factor(df_wide[, "Block"], ordered = T)
  # conditioning variate vector
  ycond <- df_wide[, conditioning]
  # threshold of the conditioning variate
  thr_cond <- LaplacesDemon::qlaplace(tau)
  # exceedances indexes:
  u_ind <- which(ycond > thr_cond)
  # exceedances of the conditioning variate + associates (data frame)
  df_cond <- df_wide[u_ind, ]
  # stationarity time blocks in factor format (a data frame):
  yblock <- data.frame(Block = factor(
    df_cond[, "Block"],
    levels = levels(df_wide[, "Block"]),
    ordered = T)
  )
  # dummy matrix of the stationary time blocks for the non-statinary scale parameter:
  almat <- model.matrix(~ Block - 1, data = yblock,
                        contrasts.arg = list(Block = diag(nlevels(yblock$Block))))

  # Fit H&T to the associate variate "associate_col":
  o = res = m = s = q = dt <- vector("list", length = length(associate_col))
  q_semipar <- vector("list", length = length(associate_col))
  alpha_hat <- matrix(NA, nrow = npal, ncol = length(associate_col))
  beta_hat = mu_hat = sigma_hat = convergence = nllh <- vector("numeric", length = length(associate_col))
  i_aux <- min(associate_col) - 1
  for (i in associate_col) {
    j <- i - i_aux
    # associate variate exceedances df:
    df_assoc <- df_cond[ , c(1:4, i)]
    # conditining and associate variates data (vector):
    y1 <- df_assoc[, 4]
    yd <- df_assoc[, 5]

    # initial values for the optim function:
    ini <- c(rep(0.01, npal), 0.5)
    WeeNumber <- 10^(-6)
    low_bound_a <- rep(-1 + WeeNumber , npal)
    upp_bound_a <- rep(1 - WeeNumber, npal)
    o[[j]] <- optim(par = ini, fn = Qpos,
                    y1 = y1, yd = yd, lambda = lambda, almat = almat,
                    method = "L-BFGS-B",
                    lower = c(low_bound_a, 0 + WeeNumber),
                    upper = c(upp_bound_a, 1 - WeeNumber),
                    control = list(maxit = 15000))

    res[[j]] <- (yd - y1 * almat %*% o[[j]]$par[1:npal]) / (y1 ^ o[[j]]$par[npal + 1])

    alpha_hat[, j] <- o[[j]]$par[1:npal]
    beta_hat[j] <- o[[j]]$par[npal + 1]
    mu_hat[j] <- mean(res[[j]])
    sigma_hat[j] <- var(res[[j]])

    high_thr <- LaplacesDemon::qlaplace(.975)
    m[[j]] <- as.vector(high_thr * almat %*% o[[j]]$par[1:npal] + mu_hat[j] * (high_thr ^ o[[j]]$par[npal + 1]))
    s[[j]] <- sqrt(sigma_hat[j]) * (high_thr ^ o[[j]]$par[npal + 1])
    q[[j]] <- qnorm(.975, mean = m[[j]], sd = s[[j]])

    q_aux <- high_thr * almat %*% o[[j]]$par[1:npal] + (high_thr ^ o[[j]]$par[npal + 1]) * res[[j]]
    q_aux2 <- sapply(split(q_aux, df_assoc[, 2]), quantile, probs = .975, USE.NAMES = FALSE)
    names(q_aux2) <- NULL
    q_semipar[[j]] <- rep(q_aux2, as.vector(table(df_assoc[, 2])))

    dt[[j]] <- data.frame(df_assoc[, 1], df_assoc[, 2], y1, yd,
                          res[[j]], m[[j]], s[[j]],
                          q[[j]], q_semipar[[j]])
    names(dt[[j]]) <- c("Time", "Block", conditioning, names(df_wide)[i],
                        "Residuals", "Mean", "Sd",
                        "Q_par", "Q_semipar")

    convergence[j] <- o[[j]]$convergence
    nllh[j] <- o[[j]]$value

  }

  z <- vector("list", length = 10)

  # data frame giving the MLE's for the alpha parameter for all associate variates:
  z$alpha <- alpha_hat
  # numeric vector giving the MLE's for the beta parameter for all associate variates:
  z$beta <- beta_hat
  # numeric vector giving the MLE's for the beta parameter for all associate variates:
  z$mu <- mu_hat
  # numeric vector giving the MLE's for the beta parameter for all associate variates:
  z$sigma <- sigma_hat
  # data frame with all MLE's for all associated variates:
  z$mle <- data.frame(rbind(z$alpha, z$beta, z$mu, z$sigma))
  colnames(z$mle) = colnames(z$alpha) <- names(df_wide)[associate_col]
  z$mle$par <- c(paste0("alpha_", 1:npal), "beta", "mu", "sigma")
  # adding time block edges to the alpha's data frame:
  z$alpha <- cbind(edges_df, z$alpha)
  # optim-convergence code - 0 indicates successful convergence:
  z$conv <- ifelse(sum(convergence) == 0, 0, 1)
  # single numeric giving the negative log-likelihood value:
  z$nllh <- mean(nllh)
  # number of data points (global):
  z$n <- n
  # number of exceedances given that the conditioning variate is high:
  z$nexc_cond <- nrow(df_cond)
  # data frame containing the data and residuals for each associate variate:
  z$data <- dt
  rownames(z$data) <- NULL

  return(z)
}


ht_chuks_fit <- function(chk, df_wide, conditioning, associate_col, tau, lambda) {

  chks_trng <- df_wide %>% dplyr::filter(Chunk %in% 1:chk) # trainning set
  chks_test <- df_wide %>% dplyr::filter(Chunk %notin% 1:chk) # test set

  ## Optimization
  trng_fit <- ppcht_fit(chks_trng, conditioning, associate_col, tau, lambda)
  trng_par <- trng_fit$mle

  ## Evaluating (negative) loglikelihood on the test set
  test_fit <- ppcht_test_fit(trng_par, df_wide, conditioning, associate_col, tau)

  return(test_fit)

}


ht_cv_fit <- function(lambda, df_wide, conditioning, associate_col, tau) {

  chks_n <- 10

  ## Removing 1 chunk in turn and computing the loglikelihood
  chks_loglik <- vector("numeric", chks_n - 5)
  for (chk in seq_along(chks_loglik)) {
    chks_loglik[chk] <- ht_chuks_fit(chk + 4, df_wide, conditioning, associate_col, tau, lambda)
  }
  cv_output <- c(lambda, chks_loglik)
  names(cv_output) <- c("lambda", paste0("chk_", 1:(chks_n - 5)))

  return(cv_output)

}


ht_boot <- function(boot_i, boot_index, df_wide, conditioning, associate_col, boot_nep, best_lambda) {

  boot_df <- df_wide[boot_index[, boot_i], ]

  fit <- ppcht_fit(boot_df, conditioning, associate_col, boot_nep[boot_i], best_lambda)

  par <- data.frame(fit$mle,
                    convergence = fit$conv,
                    boot_id = boot_i)

  for (i in 1:length(fit$data)) {
    fit$data[[i]]$boot_id <- boot_i
  }

  return(list(par = par, res = fit$data))

}


ht_fit <- function(phase_id, lead_id, phase, covariate, conditioning, leading, associates,
                   lbda_n, lbda_rng, nep, tau, boot_nep, boot_n, marginal_fit_dir, ht_fit_dir) {


  ## Variables names
  marginals <- c(conditioning, associates)
  marginals_n <- length(marginals)
  phase_name <- ifelse(phase_id == 0, "Pre-seizure:", "Post-seizure:")
  lambda <- seq(10 ^ lbda_rng[1], 10 ^ lbda_rng[2], length.out = lbda_n)


  ## Reading data from the fit of all marginals
  df <- list()
  for (i in seq_along(marginals)) {
    aux <- read_rds(paste0(marginal_fit_dir, phase_id, "_", lead_id, "_gpd_", marginals[i], "_output.rds"))
    df[[i]] <- aux$transformed_data
    df[[i]]$Marginal <- marginals[i]
  }
  df <- do.call("rbind", df)

  df_wide <- df %>%
    dplyr::select(Time, Block, Chunk, Laplace, Marginal) %>%
    pivot_wider(names_from = Marginal, values_from = Laplace) %>%
    data.frame()
  associate_col <- which(names(df_wide) %in% associates)


  ## Cross-validation ----
  cv <- foreach(l = 1:lbda_n,
                .combine = 'rbind',
                .packages = 'tidyverse',
                .export = c("%notin%", "laplace_negloglik",
                            "laplace_profile_negloglik", "Qpos",
                            "ppcht_test_fit", "ppcht_fit",
                            "ht_chuks_fit", "ht_cv_fit" )) %dopar% {
                              ht_cv_fit(lambda[l], df_wide, conditioning, associate_col, tau)
                            }

  ## CV Result
  cv_output <- data.frame(lambda = cv[, 1], loglik = rowSums(cv[, 3:6])/4)
  cv_df <- dplyr::filter(cv_output, loglik < 10 ^ 6)
  if (nrow(cv_df) == 0) {
    warning("CV failed! Penalty parameter forced to be zero.")
    best_lambda <- 0
  } else {
    best_lambda <- cv_df[which.min(cv_df[, 2]), 1]

    ## CV Plot
    png(file = paste0(ht_fit_dir, phase_id, "_", lead_id, "_ht_cv.png"),  width = 600, height = 600)
    plot_cv <- ht_cv_plot_jointly(cv_df, best_lambda, conditioning, phase_name, lead_id)
    print(plot_cv)
    dev.off()
    plot_cv


    ## Bootstrap ----
    ref_out <- read_rds(paste0(marginal_fit_dir, phase_id, "_", lead_id, "_gpd_", conditioning, "_output.rds"))
    boot_index <- ref_out$boot_index
    boot <- foreach(
      b = 1:(boot_n + 1),
      .combine = 'rbind',
      .packages = 'tidyverse',
      .export = c("laplace_negloglik",
                  "laplace_profile_negloglik", "Qpos",
                  "ppcht_fit", "ht_boot")) %dopar% {
                    ht_boot(b, boot_index, df_wide, conditioning, associate_col, boot_nep, best_lambda)
                  }
    par <- do.call("rbind", boot[, 1])
    rownames(par) <- NULL
    res <- list()
    for (j in 1:length(associates)) {
      aux <- list()
      for (i in 1:length(boot[, 2])) {
        aux[[i]] <- boot[, 2][[i]][[j]]
      }
      res[[j]] <- do.call("rbind", aux)
      rownames(res[[j]]) <- NULL
    }


    ## Saving results
    output <- list(
      data = df_wide,
      cv_chkwise = cv,
      cv_summary = cv_df,
      cv_best_lambda = best_lambda,
      boot_summary = par,
      boot_res = res
    )
    write_rds(output, paste0(ht_fit_dir, phase_id, "_", lead_id, "_output_ht.rds"))


    for (i in seq_along(associates)) {

      png(file = paste0(ht_fit_dir, phase_id, "_", lead_id, "_ht_alpha_stb_", associates[i], ".png"),  width = 250*6 , height = 250*2)
      alpha_plot <- alpha_stb(output$boot_summary[, c(associates[i], "par", "convergence", "boot_id")],
                              boot_nep, nep, tau, associates[i], phase_name, lead_id)
      print(alpha_plot)
      dev.off()

      png(file = paste0(ht_fit_dir, phase_id, "_", lead_id, "_ht_beta_stb_", associates[i], ".png"),  width = 600, height = 600)
      beta_plot <- beta_stb(output$boot_summary[, c(associates[i], "par", "convergence", "boot_id")],
                            boot_nep, nep, tau, associates[i], phase_name, lead_id)
      print(beta_plot)
      dev.off()

      png(file = paste0(ht_fit_dir, phase_id, "_", lead_id, "_ht_got_", associates[i], ".png"),  width = 1200, height = 1200/3)
      res_diag <- diagnostics_plot(dt = output$boot_res[[i]], df_wide, conditioning, associates[i], phase_name, lead_id)
      print(res_diag)
      dev.off()
    }

  }

  invisible(output)

}
