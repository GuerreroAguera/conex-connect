##' Split a data frame into blks_n blocks where each block is also splitted into chks_n folds
##'
##' @param df  data frame: to be splited into blks_n blocks
##' @param blks_n scalar: number of blocks to split the data frame in
##' @phase string: df column name identifying pre- and post-seizure phases
##' @phase_id scalar: phase to which you want to fit the model: 0 -> Pre-seizure / 1 -> Post-seizure
##' @leading string: df column name identifying the time lag
##' @lead_id scalar: lag to which you want to fit the model (0, 10, 25, 50)
##' @covariate string: df column name identifying non-stationary covariate for the model
##' @conditioning string: df column name identifying the conditioning variate for the H&T model
##' @associates string: df olumn names of the associate variables for the model
##' @param blks_id character: name to be passed to the columns with the block number identification
block_dataframe <- function(
  df, phase, phase_id, leading, lead_id,
  covariate, conditioning, associates, blks_n, chks_n
  ) {

  # Block split aux function
  block_split <- function(df, blks_n, blks_id) {

    nsplit <- blks_n - 1
    nrows <- nrow(df)
    nperdf <- ceiling(nrows / blks_n)
    start <- seq(1, nsplit * nperdf + 1, by = nperdf)

    mylist <- lapply(start, function(i) df[c(i:min((i + nperdf - 1), nrows)), ])
    final <- dplyr::bind_rows(mylist, .id = blks_id)
    final[, blks_id] <- as.double(pull(final[, blks_id]))

    return(final)

  }

  # Filtered dataset
  df <- df %>% dplyr::filter(
    !!sym(phase) == phase_id, !!sym(leading) == lead_id) %>%
    dplyr::select(all_of(c(covariate, conditioning, associates)))

  # Spliting the data into B = blks_n blocks
  df <- block_split(df, blks_n, blks_id = "Block")

  # Spliting each block into K = chks_n chunks for K-fold cross validation
  chks_df <- list()
  for (b in 1:blks_n) {
    aux <- df %>% dplyr::filter(Block == b)
    chks_df[[b]] <- block_split(aux, chks_n, "Chunk")
  }
  df <- do.call("rbind", chks_df)

  return(df)

}


#' Getting blockwise information: Channel, Threshold, Exceendances and Max (per Block)
##'
##' @param channel_df data frame: with a "Block" column indicating blocks splitting the "Channel" (the variable in analysis)
##' @param tau scalar: NEP for selecting the threshold for the GP model
info_blkwise <- function(channel_df, tau) {

  Channel <- channel_df$Channel; Block <- channel_df$Block
  cha_blkwise <- split(Channel, Block)
  thr_blkwise <- Map(f = function(x, y) quantile(x, probs = y, names = FALSE),
                     x = cha_blkwise, y = tau, USE.NAMES = FALSE)
  exc_blkwise <- Map(f = function(x, y) x[x > y] - y,
                     x = cha_blkwise, y = thr_blkwise, USE.NAMES = FALSE)
  max_blkwise <- unlist(lapply(cha_blkwise, FUN = max))

  return(list(
    Channel = cha_blkwise,
    Thresold = thr_blkwise,
    Exceedances = exc_blkwise,
    Max = max_blkwise)
  )

}
