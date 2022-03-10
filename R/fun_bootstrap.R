##' Block-Bootstrap samples, organized by columns, within pseudo-stationary time blocks for a given non-stationary time series
##'
##' @param channel_df  data frame: a N x 4 a data frame with a "Block" column indicating the pseudo stationary time blocks to be bootstrapped
##' @param boot_size scalar: bootstrap block size, must be greater than the time dependence lag in a time series
##' @param boot_n scalar: number of bootstrap samples
bootstrap_df <-  function(channel_df, boot_size, boot_n) {

  ##' Block-Bootstrap as in Lahiri (2003)
  ##'
  ##' @param y  1 x n: numeric vector of size n to be sampled by blocks
  bootstrap_index <- function(y, boot_size, boot_n) {

    gprob = 1 / boot_size
    n = length(y)

    ystar = matrix(nrow = n, ncol = boot_n + 1)
    colnames(ystar) <- paste0("Boot_", 0:boot_n)
    ystar[, 1] <- y
    for (r in 2:(boot_n + 1)) {
      loc = round(runif(1,1,n)) # loc for location
      for (i in 1:n) {
        g1 = runif(1, 0, 1)
        # In probability gprob, we take next observation, otherwise we start a new block
        if (g1 > gprob) {loc = loc + 1} else {loc = round(runif(1, 1, n))}
        if (loc > n) loc = loc - n # wrap the time series as a circle
        ystar[i, r] = y[loc]
      }
    }
    ystar <- data.frame(ystar)

    return(ystar)

  }

  index <- 1:nrow(channel_df); blks <- channel_df$Block
  index_blkwise <- split(index, blks)
  boot_blkwise <- Map(
    function(y, boot_size, boot_n) bootstrap_index(y, boot_size, boot_n),
    y = index_blkwise, boot_size = boot_size, boot_n = boot_n, USE.NAMES = FALSE)
  boot_indexes <- do.call("rbind", boot_blkwise)

  return(boot_indexes)

}
