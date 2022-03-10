##' Transform the EDF-GPD data to uniform scale
##'
##' @param par  numeric: a 1 x (B + 1) vector for the GPD fit - B scale parameters + 1 shape parameter
##' @param channel_df data frame: with a "Block" column indicating blocks to split the "Channel", the variable in analysis
##' @param tau scalar: NEP for selecting the threshold for the GPD model
u_scale <- function(par, channel_df, tau) {

  Channel <- channel_df$Channel
  Block <- channel_df$Block

  blks_n <- length(unique(Block)) # number of blocks
  nu_b <- par[-(blks_n + 1)] # gpd scale per block
  xi <- par[blks_n + 1] # gpd shape

  cha_blkwise <- split(Channel, Block) # data splitted per block

  thr_blkwise <- unlist(
    Map(f = function(x, y) quantile(x, probs = y, names = FALSE),
        x = cha_blkwise, y = tau, USE.NAMES = FALSE)
  )

  u <- list()
  for (b in 1:blks_n) {

    y <- cha_blkwise[[b]]
    u[[b]] <- vector("numeric", length = length(y))

    exceedances_index <- (y > thr_blkwise[b])

    below <- (rank(y, ties.method = "random") / (length(y) + 1))[!exceedances_index]
    u[[b]][!exceedances_index] <- below

    exceedances <- y[exceedances_index] - thr_blkwise[b]
    above <- tau + (1 - tau) * evd::pgpd(exceedances, loc = 0, scale = nu_b[b], shape = xi)
    u[[b]][exceedances_index] <- above
  }
  u <- do.call("c", u)

  return(u)

}

##' Transform data in the uniform scale to data in the standard Laplace scale
##'
##' @param u 1 x n: numeric vector with uniform data
laplace_scale <- Vectorize(
  FUN = function(u) sign(.5 - u) * log(2 * min(1 - u, u)),
  vectorize.args = "u",
  SIMPLIFY = TRUE
)
