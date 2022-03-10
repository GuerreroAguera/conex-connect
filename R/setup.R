### Loading packages ###
pcks <- c("tidyverse", "ggfortify", "gridExtra",
          "cowplot", "wesanderson", "png", "grid",
          "lazyeval", "ggsci", "akima", "scales",
          "mgcv", "evd", "texmex", "signal", "Jmisc",
          "LaplacesDemon", "foreach", "doParallel")

if (!require("librarian")) install.packages("librarian")
librarian::shelf(pcks, quiet = TRUE) # ,update_all = TRUE)
# print(.Last.value)

### Registering Clusters ###
cl <- parallel::makeCluster(4, setup_strategy = "sequential")
doParallel::registerDoParallel(cl)

### Basic setup ###
freq <- 100
eeg_start <- 1 / freq
eeg_end <- 50000 / freq
seizure <- 35000 / freq

n_points <- 15000 # Maximum of 15000 points

### Helper function to pass to the filter() function in the dyplr package.
`%notin%` <- Negate(`%in%`)

### Helper function to control axis ticks
number_ticks <- function(n) {function(limits) pretty(limits, n)}

### Helper function to get the legend of a ggplot2 object
get_legend <- function(myggplot) {

  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]

  return(legend)

}
