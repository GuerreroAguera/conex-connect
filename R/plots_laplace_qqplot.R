##' Standard Laplace Q-Q Plot
##'
##' @param x 1 x N: data to check against a Laplace(0, 1)
qqdf <- function(x) {

  n <- length(x)

  laplace_sim <- function(l) {
    dt_sim <- LaplacesDemon::rlaplace(n)
    sort(dt_sim)
  }

  envelop_sim <- sapply(1:1000, laplace_sim)
  envelop_ci <- data.frame(
    t(apply(envelop_sim, quantile, MARGIN = 1, probs = c(0.025, 0.975), names = FALSE))
  )
  names(envelop_ci) <- c("Lwr", "Upp")

  qqplot_df <- data.frame(
    Theoretical = LaplacesDemon::qlaplace(ppoints(n)),
    Empirical = sort(x),
    envelop_ci
  )

  return(qqplot_df)

}
