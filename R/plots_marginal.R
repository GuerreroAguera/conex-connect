# General theme for the plots
my_theme <- theme_linedraw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold")
  )


#### ACF and PACF Plots ----

acf_plot <- function(channel_df, marginals, phase_name,
                     lag.max = 100, ci = .95, boot_size = 75,
                     lead_id) {

  df_acf <- channel_df %>%
    group_by(Block) %>%
    summarise(list_acf = list(acf(Channel, lag.max = lag.max, plot = FALSE))) %>%
    mutate(acf_vals = purrr::map(list_acf, ~as.numeric(.x$acf))) %>%
    dplyr::select(-list_acf) %>%
    unnest(cols = acf_vals) %>%
    group_by(Block) %>%
    mutate(lag = row_number() - 1)

  df_pacf <- channel_df %>%
    group_by(Block) %>%
    summarise(list_acf = list(pacf(Channel, lag.max = lag.max, plot = FALSE))) %>%
    mutate(acf_vals = purrr::map(list_acf, ~as.numeric(.x$acf))) %>%
    dplyr::select(-list_acf) %>%
    unnest(cols = acf_vals) %>%
    group_by(Block) %>%
    mutate(lag = row_number() - 1)

  df_ci <- channel_df %>%
    group_by(Block) %>%
    summarise(ci = qnorm((1 + ci) / 2) / sqrt(n()))

  strips <- paste0("Block ", 1:blks_n)
  names(strips) = 1:blks_n

  acf_plot <- ggplot(df_acf, aes(x = lag, y = acf_vals)) +
    geom_bar(stat = "identity", width = .55) +
    geom_hline(yintercept = 0) +
    geom_hline(data = df_ci, aes(yintercept = -ci), color = "blue", linetype = "dotted", size = 1) +
    geom_hline(data = df_ci, aes(yintercept = ci), color = "blue", linetype = "dotted", size = 1) +
    geom_vline(xintercept = boot_size, col = "red", linetype = "dashed", size = .75) +
    ggtitle(paste0("Channel ", marginals, " ", phase_name, ": ACF"),
            subtitle = paste0("Block Bootstrap Size = ", boot_size, " / Lag = ", lead_id)) +
    labs(x = "Lag", y = "ACF") +
    facet_wrap(~Block, ncol = 4, labeller = labeller(Block = strips)) +
    scale_x_continuous(breaks = number_ticks(6), expand = c(0, 2)) +
    scale_y_continuous(breaks = number_ticks(6),expand = c(0, 0)) +
    my_theme +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 12, hjust = 1)
    )

  pacf_plot <- ggplot(df_pacf, aes(x = lag, y = acf_vals)) +
    geom_bar(stat = "identity", width = .75) +
    geom_hline(yintercept = 0) +
    geom_hline(data = df_ci, aes(yintercept = -ci), color = "blue", linetype = "dotted", size = 1) +
    geom_hline(data = df_ci, aes(yintercept = ci), color = "blue", linetype = "dotted", size = 1) +
    geom_vline(xintercept = boot_size, col = "red", linetype = "dashed", size = .75) +
    ggtitle(paste0("Channel ", marginals, " ", phase_name, ": PACF"),
            subtitle = paste0("Block Bootstrap Size = ", boot_size, " / Lag = ", lead_id)) +
    labs(x = "Lag", y = "ACF") +
    facet_wrap(~Block, ncol = 4, labeller = labeller(Block = strips)) +
    scale_x_continuous(breaks = number_ticks(6), expand = c(0, 2)) +
    scale_y_continuous(breaks = number_ticks(6),expand = c(0, 0)) +
    my_theme +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 12, hjust = 1)
    )

  return(list(acf = acf_plot, pacf = pacf_plot))

}


#### Cross-validation Plot ----

cv_plot <- function(cv_df, best_lambda, marginals, phase_name, lead_id) {

  cv_plot <- cv_df %>%
    ggplot(aes(x = lambda, y = loglik)) +
    geom_line(size = 1) + geom_point(size = 2) +
    geom_vline(xintercept = best_lambda,
               size = 1.2, col = "red", linetype = "dashed") +
    geom_label(aes(x = best_lambda, y = max(loglik),
                   label = signif(best_lambda, 2)),
               colour = "white", fill = "red", fontface = "bold" ) +
    scale_y_log10(breaks = number_ticks(8)) +
    scale_x_continuous(breaks = number_ticks(10)) +
    ggtitle(paste0("Channel ", marginals, " ", phase_name, ": GPD Cross-Validation"),
            subtitle = paste0("Roughness Penalty Selection / Lag = ", lead_id)) +
    labs(x = expression(lambda), y = "Negative Log-likelihood \n (log10 Scale)") +
    my_theme

  return(cv_plot)

}


#### Original and Transfomed Data Plots ----

disp_plot <- function(fit_df, tau, boot_nep, edges_bar, marginals, phase_name, lead_id) {

  scaleFUN <- function(x) sprintf("%.1f", x)
  edges_line <- geom_vline(xintercept = edges_bar, col = "red", size = 0.75, linetype = "dashed")
  edges_scale <- scale_x_continuous(breaks = edges_bar, labels = scaleFUN, expand = c(0.025, 0.025))


  thr_interval <- fit_df %>%
    dplyr::group_by(Block) %>%
    mutate(
      Lwr = quantile(Channel, min(boot_nep)),
      Thr = quantile(Channel, boot_nep[1]),
      Upp = quantile(Channel, max(boot_nep))
    ) %>% ungroup()

  plot_raw_data <- thr_interval %>%
    ggplot(aes(x = Time, y = Channel, color = Channel > Thr)) + geom_point() +
    geom_step(data = thr_interval, aes(x = Time, y = Thr), col = "navyblue", size = 1.2) +
    geom_step(data = thr_interval, aes(x = Time, y = Lwr),
              col = "blue", size = 1, linetype = "dotdash") +
    geom_step(data = thr_interval, aes(x = Time, y = Upp),
              col = "blue", size = 1, linetype = "dotdash") +
    scale_color_manual(values = c("gray", "black")) +
    ggtitle(paste0("Channel ", marginals, " ", phase_name, ": Raw data"),
            subtitle = paste0("NEP = ", tau, " / Lag = ", lead_id)) +
    labs(y = "EEG Amplitude (Absolute Values)") +
    edges_line + edges_scale +
    scale_y_continuous(breaks = number_ticks(8), expand = c(0, 0)) +
    my_theme + theme(legend.position = "none",
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 12, angle = 45, hjust = 1))

  plot_u_scale <- fit_df %>%
    ggplot(aes(x = Time, y = Uniform)) + geom_point() +
    ggtitle("On Uniform Margins", subtitle = paste0("NEP = ", tau, " / Lag = ", lead_id)) +
    labs(y = paste0(marginals, " - Uniform Scale")) +
    edges_line + edges_scale +
    scale_y_continuous(breaks = number_ticks(8), expand = c(0, 0)) +
    my_theme +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1))

  plot_laplace_scale <- fit_df %>%
    ggplot(aes(x = Time, y = Laplace)) + geom_point() +
    ggtitle("On Laplace Margins", subtitle = paste0("NEP = ", tau, " / Lag = ", lead_id)) +
    labs(y = paste0(marginals, " - Laplace Scale")) +
    edges_line + edges_scale +
    scale_y_continuous(breaks = number_ticks(8), expand = c(0, 0.)) +
    my_theme +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1))

  return(
    cowplot::plot_grid(plot_raw_data, plot_u_scale, plot_laplace_scale, ncol = 3)
  )

}


#### GP Parameters ----

# Scale (steps plot) and Shape (histogram)
par_plot <- function(nu_df, xi_df, xi_hat0, edges_bar, marginals, phase_name, tau, boot_n, lead_id) {

  aux_df <- nu_df[, -2]
  names(aux_df) <- c("Max", "Par", "Med", "Lwr", "Upp")
  aux_df <- rbind(aux_df, nu_df[nrow(nu_df), -1])

  edges_line <- geom_vline(xintercept = edges_bar, col = "red", size = 0.75, linetype = "dashed")
  edges_scale <- scale_x_continuous(breaks = edges_bar, labels = edges_bar, expand = c(0.025, 0.025))

  cols <- c("no_boot" = "blue", "mean" = "black", "ci" = "black")
  cols_labs <- c("no_boot" = "Estimate", "mean" = "Bootstrap Mean", "ci" = "95% C.I.")
  lines_labs <- c("no_boot" = "solid", "mean" = "solid", "ci" = "dotdash")

  nu_plot <- aux_df %>% ggplot() +
    geom_step(aes(x = Max, y = Med, color = "mean", linetype = "mean"), size = 1.25) +
    geom_step(aes(x = Max, y = Lwr, color = "ci", linetype = "ci"), size = 1.35) +
    geom_step(aes(x = Max, y = Upp, color = "ci", linetype = "ci"), size = 1.35) +
    geom_segment(data = nu_df,
                 aes(x = Min, y = Par, xend = Max, yend = Par,
                     color = "no_boot", linetype = "no_boot"),
                 size = 1.15) +
    scale_color_manual(name = "my_legend", values = cols, labels = cols_labs) +
    scale_linetype_manual(name = "my_legend", values = lines_labs, labels = cols_labs) +
    ggtitle(paste0("Channel ", marginals, " ", phase_name, ": GPD Scale"),
            subtitle = paste0("NEP = ", tau, " / ", boot_n, " Bootstrap Samples / Lag = ", lead_id)) +
    labs(color = NULL, y = expression(hat(nu)), x = "Time") +
    edges_line + edges_scale +
    scale_y_continuous(breaks = number_ticks(8)) +
    my_theme +
    theme(
      legend.position = c(0.860, 0.925),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      panel.grid.minor = element_blank()
    )

  xi_mean <- round(mean(c(xi_hat0, xi_df$Par)), 4)
  xi_sd <- round(sd(c(xi_hat0, xi_df$Par)), 4)
  xi_plot <- xi_df %>% ggplot(aes(x = Par)) +
    geom_histogram(position = "identity", fill = "black", alpha = 0.925, bins = 15) +
    geom_vline(xintercept = xi_hat0, col = "red", size = 1.25) +
    scale_x_continuous(breaks = number_ticks(8), expand = c(0, 0)) +
    scale_y_continuous(breaks = number_ticks(8), expand = c(0, 0)) +
    ggtitle(paste0("Channel ", marginals, " ", phase_name, ": GPD Shape"),
            subtitle = paste0("NEP = ", tau, " / ", boot_n, " Bootstrap Samples / Lag = ", lead_id)) +
    labs(y = "", x = bquote(hat(xi))) +
    my_theme + theme(panel.grid.minor = element_blank())

  yhist_max <- max(ggplot_build(xi_plot)$data[[1]]$y)
  yhist_min <- min(ggplot_build(xi_plot)$data[[1]]$y[ggplot_build(xi_plot)$data[[1]]$y != 0])
  xhist_max <- max(ggplot_build(xi_plot)$data[[1]]$x)
  xi_plot <- xi_plot +
    geom_label(aes(x = xi_hat0, y = 3 * yhist_min, label = round(xi_hat0, 4)),
               colour = "white", fill = "red", fontface = "bold") +
    annotate("label", x = xi_hat0, y = 0.975 * yhist_max,
             label = paste0("mean = ", xi_mean, " / sd = ", xi_sd),
             colour = "white", fill = "black", fontface = "bold")

  return(
    cowplot::plot_grid(nu_plot, xi_plot, ncol = 2, align = "h")
  )

}


# Shape stability
stb_plot <- function(tau, nep, boot_nep, boot, marginals, phase_name, lead_id) {

  stb_brks <- c(seq(nep[1], tau, length.out = 6),
                seq(tau, nep[2], length.out = 6)[-1])
  stb_df <- data.frame(NEP = boot_nep, Par = boot[, (blks_n + 1)])

  stb_plot <- stb_df %>%
    ggplot(aes(x = NEP, y = Par)) +
    geom_point() +
    geom_smooth(method = "loess", formula = "y ~ x", col = "red") +
    scale_x_continuous(breaks = stb_brks, labels = stb_brks, expand = c(0, 0)) +
    scale_y_continuous(breaks = number_ticks(8), expand = c(0, 0.05)) +
    coord_cartesian(xlim = nep) +
    ggtitle(paste0("Channel ", marginals, " ", phase_name, ": GPD Shape Stability by Threshold"),
            subtitle = paste0(boot_n, " Bootstrap Samples / Lag = ", lead_id)) +
    labs(y = expression(hat(xi)), x = "NEP") +
    my_theme +
    theme(
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(5.5, 12.5, 5.5, 5.5), "pt")
    )


  return(stb_plot)

}


#### Q-Q Plots ----

# Overall
allqq_plot <- function(qqplot_overall, marginals, phase_name, lead_id) {

  qqplot <- qqplot_overall %>%
    ggplot(aes(Theoretical, Empirical)) +
    geom_point(col = "darkred", fill = "red", shape = 21) +
    geom_abline(slope = 1, intercept = 0, size = 0.85) +
    geom_line(aes(x = Theoretical, y = Lwr), size = 0.75, linetype = "dashed") +
    geom_line(aes(x = Theoretical, y = Upp), size = 0.75, linetype = "dashed") +
    scale_x_continuous(breaks = number_ticks(8), expand = c(0, 0)) +
    scale_y_continuous(breaks = number_ticks(8), expand = c(0, 0)) +
    ggtitle(paste0("Channel ", marginals, " ", phase_name, ": Laplace Scale Q-Q Plot"),
            subtitle = paste0("Overall Goodness of Fit / Lag = ", lead_id)) +
    labs(y = "Empirical Quantile", x = "Theoretical Quantile") +
    my_theme + theme(panel.grid.minor = element_blank())

  return(qqplot)

}

# Blockwise
blkqq_plot <- function(qqplot_blkwise, strip_labs, marginals, phase_name, lead_id) {

  qqplot <- qqplot_blkwise %>%
    ggplot(aes(Theoretical, Empirical)) +
    geom_point(col = "darkred", fill = "red", shape = 21) +
    geom_abline(slope = 1, intercept = 0, size = 0.85) +
    geom_line(aes(x = Theoretical, y = Lwr), size = 0.75, linetype = "dashed") +
    geom_line(aes(x = Theoretical, y = Upp), size = 0.75, linetype = "dashed") +
    facet_wrap(~Block, ncol = 4, scales = "free",
               labeller = labeller(Block = strip_labs)) +
    scale_x_continuous(breaks = number_ticks(8), expand = c(0, 0)) +
    scale_y_continuous(breaks = number_ticks(8), expand = c(0, 0)) +
    ggtitle(paste0("Channel ", marginals, " ", phase_name, ": Laplace Scale Q-Q Plot"),
            subtitle = paste0("Blockwise Goodness of Fit / Lag = ", lead_id)) +
    labs(y = "Empirical Quantile", x = "Theoretical Quantile") +
    my_theme + theme(panel.grid.minor = element_blank())

  return(qqplot)

}
